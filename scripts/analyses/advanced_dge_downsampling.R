## advanced dge downsampling based on umi_observations.txt files. allows to downsample on an average
## or fixed number of reads, that align to all or only a set of target genes.

# attach required packages
library(tidyverse)
library(here)

# define functions ---------------------------------------------------------------------------------

# calculate transcripts-per-transcript for umi observations to filter for chimeric reads
calc_tpt <- function(umi_obs) {
  
  # calculate total number of reads for every cbc - umi combination
  total_reads <- umi_obs %>%
    group_by(Cell.Barcode, Molecular_Barcode) %>%
    summarize(total_reads = sum(Num_Obs))
  
  # add total reads to umi_obs and calculate TPT
  umi_obs %>%
    left_join(total_reads, by = c("Cell.Barcode", "Molecular_Barcode")) %>%
    mutate(tpt = Num_Obs / total_reads) %>%
    select(-total_reads)
  
}

# filter according to a specified tpt threshold
filter_tpt <- function(umi_obs_tpt, tpt_threshold) {
  umi_obs_tpt %>%
    filter(tpt >= tpt_threshold) %>%
    select(-tpt)
}

# downsample umi_obs for a specifc number of genic reads
downsample_reads <- function(umi_obs, reads, level = c("cell", "sample")) {
  
  level <- match.arg(level)
  
  # "uncount" umi_observations so that each read is a row by repeating each line number-of-reads
  # times (i.e. num_obs)
  reps <- rep(seq_len(nrow(umi_obs)), times = umi_obs$Num_Obs)
  umi_obs_reads <- umi_obs[reps, -4]
  
  # randomly sample reads per cell or the whole sample, depending on level parameter
  if (level == "cell") {
    sampled_reads <- umi_obs_reads %>%
      group_by(Cell.Barcode) %>%
      sample_n(size = reads, replace = FALSE)
  }else{
    sampled_reads <- sample_n(umi_obs_reads, size = reads, replace = FALSE)
  }
  
  # count observations of each sample-cell-umi combination
  sampled_umi_obs <- sampled_reads %>%
    group_by(Cell.Barcode, Gene, Molecular_Barcode) %>%
    summarize(Num_Obs = n()) %>%
    ungroup()
  
}

# extract dge data for each sample by calculating number of umis per gene and  cell
extract_dge <- function(umi_obs){
  umi_obs %>%
    count(Cell.Barcode, Gene) %>%
    rename(txs = n)
}

# convert dge from long format to wide matrix
create_dge_matrix <- function(dge_long) {
  dge_long %>%
    spread(key = Cell.Barcode, value = txs, fill = 0) %>%
    as.data.frame(stringsAsFactors == FALSE) %>%
    rename(GENE = Gene)
}


# load input data ----------------------------------------------------------------------------------

# get parameters
sample_reads <- as.numeric(snakemake@wildcards$reads)
sampling <- match.arg(snakemake@wildcards$sampling, choices = c("fixed", "avg"))
panel <- match.arg(snakemake@wildcards$panel, choices = c("onGenes", "onTarget"))

# set sed for reproducible sampling
set.seed(snakemake@params$seed)

# load input data
umi_obs <- read.table(here(snakemake@input$umi_obs), header = TRUE, stringsAsFactors = FALSE,
                      sep = "\t")
target_genes <- read.csv(here(snakemake@input$target_genes), stringsAsFactors = FALSE)
cbc_whitelist <- readLines(here(snakemake@input$whitelist))

# filter for cell barcodes on CBC whitelist
umi_obs_filt <- filter(umi_obs, Cell.Barcode %in% cbc_whitelist)

# remove any vector transcripts if vector_pattern is provided
if (!is.null(snakemake@params$vector_pattern)) {
  
  umi_obs_filt <- filter(umi_obs_filt, !grepl(Gene, pattern = snakemake@params$vector_pattern))
  
}

# only retain target genes if panel is set to onTarget
if (panel == "onTarget") {
  
  umi_obs_filt <- filter(umi_obs_filt, Gene %in% target_genes$gene)
  
}


# downsample genic reads ---------------------------------------------------------------------------

if (sampling == "fixed") {
  
  # calculate number of reads per cell
  reads_per_cell <- umi_obs_filt %>%
    group_by(Cell.Barcode) %>%
    summarize(reads = sum(Num_Obs))
  
  # filter for cells that have at least the desired sampling size of reads
  cells_filt <- reads_per_cell %>%
    filter(reads >= sample_reads) %>%
    pull(Cell.Barcode)
  
  # retain umi observations ony for these cells
  umi_obs_filt <- filter(umi_obs_filt, Cell.Barcode %in% cells_filt)
  
  # downsample umi observations
  ds_umi_obs <- downsample_reads(umi_obs_filt, reads = sample_reads, level = "cell")
  
}else{
  
  # calculate number of reads to sample from umi_observation
  ncells <- n_distinct(umi_obs_filt$Cell.Barcode)
  sample_reads <- ncells * sample_reads
  
  # abort if scaled sampling size is larger than the number of available reads
  if (sample_reads > sum(umi_obs_filt$Num_Obs)) {
    stop("Not enough input reads for desired sampling size per sample!")
  }
  
  # downsample umi observations
  ds_umi_obs <- downsample_reads(umi_obs_filt, reads = sample_reads, level = "sample")

}

# extract dge --------------------------------------------------------------------------------------

# calculate tpt and filter for chimeric reads using the specified tpt threshold
umi_obs_tpt <- calc_tpt(ds_umi_obs)
umi_obs_tpt <- filter_tpt(umi_obs_tpt, tpt_threshold = snakemake@params$tpt_threshold)

# extract dge and create dge matrix
dge <- extract_dge(umi_obs_tpt)
dge_mat <- create_dge_matrix(dge)

# output file
if (tools::file_ext(snakemake@output[[1]]) == "gz") {
  outfile <- gzfile(here(snakemake@output[[1]]))
}else{
  outfile <- here(snakemake@output[[1]])
}

# save dge matrix to output file
write.table(dge_mat, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
