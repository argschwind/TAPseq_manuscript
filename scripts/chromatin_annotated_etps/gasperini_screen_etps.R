## Create chromatin annotated ETPs for the enhancer screen in Gasperini et al, 2019

# attach required packages
library(GenomicAlignments)
library(GenomicFeatures)
library(here)
library(HiTC)
library(Matrix)
library(rtracklayer)
library(tidyverse)
library(BiocParallel)

# backend for parallel computing
register(MulticoreParam(workers = snakemake@threads))

# Expressed genes ==================================================================================

# Whole-Tx CROP-seq data is used to find expressed genes in K562 cells.

# import gencode hg38 annotations
annot_url <- "ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.chr.gtf.gz"
annot <- import(annot_url, format = "gtf")

# extract exon annotations for protein-coding genes on autosomes and chromosome X
genes <- annot[annot$type == "exon" & annot$gene_biotype == "protein_coding" &
                 annot$transcript_biotype == "protein_coding" &
                 !seqnames(annot) %in% c("Y", "MT")]

# genome-wide K562 expression data
dge <- read.table(here(snakemake@input$wholeTx_dge), header = TRUE, row.names = "GENE",
                  stringsAsFactors = FALSE)

# only retain dge data on protein-coding genes
gene_expr <- dge[rownames(dge) %in% genes$gene_name, ]

# compute average gene expression
avg_expr <- data.frame(gene = rownames(gene_expr), avg_txs = rowSums(gene_expr),
                       stringsAsFactors = FALSE, row.names = NULL)

# extract genes with 100 or more transcripts on average
expr_gene_names <- avg_expr %>% 
  filter(avg_txs >= 100) %>% 
  pull(gene)

# get annotations for these genes and split by gene name
expr_genes <- genes[genes$gene_name %in% expr_gene_names]
expr_genes <- split(expr_genes, f = expr_genes$gene_name)

# Select enhancers - gene pairs ====================================================================

# Tested enhancers from the enhancer screen in Gasperini et al. 2019 are used to create ETPs with
# all expressed genes within 300kb.

# import enhancer sites
enh_coords <- import(snakemake@input$enhancers)

# remove duplicates
enh_coords <- unique(enh_coords)

# create enhancer id
enh_coords$perturbation <- paste0(seqnames(enh_coords), ":", start(enh_coords), "-",
                                  end(enh_coords))

# download chain file for hg19 to hg38 liftover
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_file <- tempfile("hg19ToHg38", fileext = ".gz")
download.file(chain_url, chain_file)
system(paste("gunzip", chain_file))

# import chain file
hg19_to_hg38_chain <- import.chain(tools::file_path_sans_ext(chain_file))

# liftover enhancers from hg38 to hg19
enh_coords <- enh_coords %>%
  liftOver(chain = hg19_to_hg38_chain)

# filter out enhancers that did not lift over correctly (split into 2 loci)
enh_loci <- sapply(enh_coords, length)
enh_loci_good <- which(enh_loci == 1)
enh_coords <- enh_coords[enh_loci_good]

# convert to GRanges object
enh_coords <- unlist(enh_coords)

# function to get a genes TSS
get_tss <- function(x) {
  if (all(strand(x) == "+")) {
    tss <- min(start(x))
  }else if (all(strand(x) == "-")) {
    tss <- max(end(x))
  }else{
    stop("Inconsistent strand information!")
  }
  GRanges(seqnames = unique(seqnames(x)), ranges = IRanges(start = tss, end = tss),
          strand = unique(strand(x)), gene_name = unique(x$gene_name))
}

# get tss for all genes
gene_tss <- endoapply(expr_genes, FUN = get_tss)
gene_tss <- unlist(gene_tss)

# adjust chromosome names of enhancer coordinates
seqlevels(enh_coords) <- sub("chr", "", seqlevels(enh_coords))

# find enhancers that overlap a 600kb window centered on each genes tss
tss_windows <- promoters(gene_tss, upstream = 300000, downstream = 300000)
overlaps <- findOverlaps(query = enh_coords, subject = tss_windows)

# get overlapping gene and gene tss
ol_enh <- enh_coords[from(overlaps)]
ol_tss <- gene_tss[to(overlaps)]

# create data.frame containing etps
enh_df <- data.frame(ol_enh)
tss_df <- data.frame(ol_tss)
etps <- data.frame(perturbation = enh_df$perturbation, gene = tss_df$gene_name,
                   enh_chr = as.character(enh_df$seqnames), enh_start = enh_df$start,
                   enh_end = enh_df$end, gene_chr = as.character(tss_df$seqnames),
                   gene_tss = tss_df$start, gene_strand = as.character(tss_df$strand),
                   stringsAsFactors = FALSE)

# calculate distance to tss for every enhancer - gene pair
etps <- etps %>%
  mutate(enh_center = round((enh_start + enh_end) / 2)) %>% 
  mutate(dist_to_tss = if_else(gene_strand == "+", true  = enh_center - gene_tss,
                               false = gene_tss - enh_center)) %>% 
  select(-enh_center)

# reformat chromosome ids
etps <- mutate(etps, enh_chr = paste0("chr", enh_chr), gene_chr = paste0("chr", gene_chr))

# Chromatin activity ===============================================================================

# Each ETP is annotated for chromatin activity at the involved enhancer. The RPKM is calculated to
# normalize for both the sequencing depth and enhancer size.

# chromatin assays and input files
chrom_infiles <- here(snakemake@input$encode_bam)
names(chrom_infiles) <- sub("_.+", "", basename(chrom_infiles))

# load all bam files
chrom_reads <- bplapply(chrom_infiles, FUN = readGAlignments)

# extract hg38 coordinates of all enhancers
enh_coords <- etps %>%
  select(enh_chr, enh_start, enh_end, perturbation) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# count the number of reads within each enhancer
enh_chrom_reads <- chrom_reads %>%
  bplapply(FUN = countOverlaps, query = enh_coords) %>%
  bind_cols()

# combine with enhancer coordinates and convert to long format
enh_chrom <- enh_coords %>%
  data.frame(stringsAsFactors = FALSE) %>%
  select(-strand) %>%
  rename(enh_chr = seqnames, enh_start = start, enh_end = end) %>%
  bind_cols(enh_chrom_reads) %>%
  gather(key = assay, value = reads, -c(1:5))

# normalize for sequencing depth and enhancer size by calculating rpkm
enh_chrom_norm <- enh_chrom %>%
  group_by(assay) %>%
  mutate(rpkm = reads / (sum(reads) / 1e6) / (width / 1000))

# convert to wide format add to enhancer gene pair results
etps <- enh_chrom_norm %>%
  select(-c(enh_chr, enh_start, enh_end, width, reads)) %>%
  spread(key = "assay", value = "rpkm") %>%
  left_join(x = etps, y = ., by = "perturbation")

# remove data to free up memory
rm(enh_chrom_reads, chrom_reads, enh_chrom, enh_chrom_norm)
rm(annot, dge, expr_genes, gene_expr, genes)
invisible(gc())

# Annotate Hi-C ====================================================================================

# Hi-C contact matrices from Rao et al., 2014 are imported and the contact frequency for each
# enhancer - gene pair is calculated.

# Liftover enhancer and gene TSS coordinates -------------------------------------------------------

# Hi-C contact matrices are in hg19 annotation, so enhancer and gene TSS coordinates need to be
# lifted over from hg38 to hg19.

# download chain file for hg38 to hg19 liftover
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
chain_file <- tempfile("hg38ToHg19", fileext = ".gz")
download.file(chain_url, chain_file)
system(paste("gunzip", chain_file))

# import chain file
hg38_to_hg19_chain <- import.chain(tools::file_path_sans_ext(chain_file))
closeAllConnections()  # close any open file connection

# extract hg38 coordinates of all enhancers
enh_coords_hg38 <- etps %>%
  select(enh_chr, enh_start, enh_end, perturbation) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# liftover enhancers from hg38 to hg19 and convert to data.frame
enh_coords_hg19 <- enh_coords_hg38 %>%
  liftOver(chain = hg38_to_hg19_chain)

# filter out enhancers that did not lift over correctly (split into 2 loci)
enh_loci <- sapply(enh_coords_hg19, length)
enh_loci_good <- which(enh_loci == 1)
enh_coords_hg19 <- enh_coords_hg19[enh_loci_good]

# reformat
enh_coords_hg19 <- enh_coords_hg19 %>%
  unlist() %>%
  as.data.frame() %>%
  select(seqnames, start, end, perturbation) %>%
  rename(enh_chr = seqnames, enh_start = start, enh_end = end) %>%
  mutate(enh_chr = as.character(enh_chr))

# replace tested enhancer - gene pairs hg38 enhancer coordinates with hg19 coordinates
etps <- etps %>%
  select(-c(enh_chr, enh_start, enh_end)) %>%
  left_join(enh_coords_hg19, by = "perturbation")

# extract gene TSS coordinates
gene_tss_coords <- etps %>%
  select(gene, gene_chr, gene_tss, gene_strand) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "gene_chr",
                           start.field = "gene_tss", end.field = "gene_tss",
                           strand.field = "gene_strand")

# liftover tss coordinates to hg19
gene_tss_coords_hg19 <- gene_tss_coords %>%
  liftOver(chain = hg38_to_hg19_chain) %>%
  unlist() %>%
  as.data.frame() %>%
  select(seqnames, start, strand, gene) %>%
  rename(gene_chr = seqnames, gene_tss = start, gene_strand = strand) %>%
  mutate(gene_chr = as.character(gene_chr), gene_strand = as.character(gene_strand))

# replace gene tss hg38 coordinates with hg19 coordinates
etps <- etps %>%
  select(-c(gene_chr, gene_tss, gene_strand)) %>%
  left_join(gene_tss_coords_hg19, by = "gene")

# recalculate distance to tss for hg19 coordinates
etps <- etps %>%
  mutate(enh_center = round((enh_start + enh_end) / 2)) %>% 
  mutate(dist_to_tss = if_else(gene_strand == "+", true  = enh_center - gene_tss,
                               false = gene_tss - enh_center)) %>% 
  select(-enh_center)

# remove enhancers that were not lifted over or switched chromosomes...
etps <- filter(etps, enh_chr == gene_chr)

# Prepare Hi-C data --------------------------------------------------------------------------------

# Contact matrices for K562 cells from Rao et al. 2014 are imported (5kb resolution). The provided
# normalization vectors are used to normalize the observed contacts (read counts). Data from the two
# chromosonal regions on chromosome 8 and 11 is extracted for any further analyses.

# function to import HiC data from Rao et al for one chromosome and create a HTCexp object
import_hic <- function(sparse_cm_file, chromosome, resolution, bins) {
  
  # load sparse contact matrix file (only observed contacts)
  obs_contacts <- read.table(sparse_cm_file, col.names = c("i", "j", "M_ij"), sep = "\t")
  
  # get starting coordinates of assessed genomic bins at 5kb resolution
  max_bin <- (bins - 1) * resolution
  bin_starts <- seq(from = 0, to = max_bin, by = resolution)
  
  # create GRanges object containing all assessed bins for that chromosome
  bin_coords <- GRanges(seqnames = chromosome,
                        ranges = IRanges(start = bin_starts, end = bin_starts + resolution - 1,
                                         names = paste0("bin_", seq_len(length(bin_starts))))
  )
  
  # convert starting coordinates of bins in sparse matrix input to bin ids by dividing by the
  # resolution (and add 1 to get correct index)
  obs_contacts_bins <- data.frame(i = round(obs_contacts$i / resolution + 1),
                                  j = round(obs_contacts$j / resolution + 1), 
                                  M_ij = obs_contacts$M_ij)
  
  # create sparse contact matrix from observed contacts
  sparse_cm <- sparseMatrix(i = obs_contacts_bins$i, j = obs_contacts_bins$j,
                            x = obs_contacts_bins$M_ij, symmetric = TRUE, dims = c(bins, bins))
  
  # create HTCexp object containing data for the given chromosome
  HTCexp(intdata = sparse_cm, xgi = bin_coords, ygi = bin_coords)
  
}

# k562 intrachromosomal sparse matrix files
scm_files <- here(snakemake@input$hic_raw)
names(scm_files) <- sub("_5kb.RAWobserved", "", basename(scm_files))

# normalization vector files
krnorm_files <- here(snakemake@input$hic_norm)
names(krnorm_files) <- sub("_5kb.KRnorm", "", basename(krnorm_files))

# import normalization vectors
krnorm <- lapply(krnorm_files, FUN = function(x) as.numeric(readLines(x)) )

# infer number of bins per chromosome based on the normalization vectors
hic_bins <- lapply(krnorm, FUN = length)

# import hi-c data for these chromosomes
hic_data <- bpmapply(FUN = import_hic,
                     sparse_cm_file = scm_files,
                     chromosome = names(scm_files),
                     bins = hic_bins,
                     MoreArgs = list(resolution = 5000))

## Enhancer - gene pairs ---------------------------------------------------------------------------

# The interaction frequency of all tested enhancer - gene pair computed, defined as the interaction
# frequency of the genomic bins overlapping the enhancer and the target genes TSS.

# function to compute the interaction frequency for enhancer - gene pairs, by finding the hi-c 
# genomic bins with overlap with the enhancer and the target gene tss. the interaction frequency of
# the pair is then defined as the interaction frequency of these bins
compute_int_freq <- function(pairs, htc_object, norm_vector) {
  
  # add pair identifier
  pairs$pair <- seq_len(nrow(pairs))
  
  # get coordinates of enhancer centers as GRanges object
  enh_coords <- pairs %>%
    mutate(enh_center = round((enh_start + enh_end) / 2)) %>%
    select(enh_chr, enh_center) %>%
    makeGRangesFromDataFrame(., seqnames.field = "enh_chr", start.field = "enh_center",
                             end.field = "enh_center")
  
  # get gene tss coordinates as GRanges object
  tss_coords <- pairs %>%
    select(gene_chr, gene_tss) %>%
    makeGRangesFromDataFrame(., seqnames.field = "gene_chr", start.field = "gene_tss",
                             end.field = "gene_tss")
  
  # get bins overlapping for all enhancers and gene tss
  hic_bins <- x_intervals(htc_object)
  enh_bins <- findOverlaps(query = enh_coords, subject = hic_bins)
  tss_bins <- findOverlaps(query = tss_coords, subject = hic_bins)
  
  # combine into one data.frame
  enh_bins <- data.frame(pair = from(enh_bins), enh_bin = to(enh_bins))
  tss_bins <- data.frame(pair = from(tss_bins), tss_bin = to(tss_bins))
  bins <- full_join(enh_bins, tss_bins, by = "pair")
  
  # extract distance matrix between bins from htc object
  dists <- intervalsDist(htc_object)
  
  # calculate distances between bins of enhancer gene pairs
  dist_pairs <- dists[as.matrix(bins[, 2:3])]
  dist_pairs <- data.frame(pair = bins$pair, dist_bins = dist_pairs)
  
  # extract hi-c interaction matrix from htc object
  intdata <- intdata(htc_object)
  
  # get interaction frequencies and normalization factors for all bins
  int_freq_pairs <- intdata[as.matrix(bins[, 2:3])]
  
  # compute normalization factors
  norm_vector_values <- matrix(norm_vector[as.matrix(bins[, 2:3])], ncol = 2)
  norm_factors <- norm_vector_values[, 1] * norm_vector_values[, 2]
  
  # normalize data and add and add pair id
  int_freq_pairs <- int_freq_pairs / norm_factors
  int_freq_pairs <- data.frame(pair = bins$pair, int_freq = int_freq_pairs)
  
  # add interaction frequencies and bin distances to pairs to create output
  pairs %>%
    left_join(dist_pairs, by = "pair") %>%
    left_join(int_freq_pairs, by = "pair") %>%
    select(-c(pair))
  
}

# compute interaction frequencies for all tested enhancer - gene pairs
etps_chrs <- split(etps, f = etps$enh_chr)
hic_data <- hic_data[names(etps_chrs)]
krnorm <- krnorm[names(etps_chrs)]
int_freqs <- bpmapply(FUN = compute_int_freq, pairs = etps_chrs, htc_object = hic_data,
                      norm_vector = krnorm, SIMPLIFY = FALSE)

# combine into one data.frame
chrom_annot_etps <- bind_rows(int_freqs)

# Activity by contact ==============================================================================

# The activity by contact score is computed as proposed by Fulco et al. The geometric mean between
# H3K27ac and DNAse-seq is then calculatated for every site to compute an enhancer activity score
# proposed by Fulco et al. The ABC score is then computed by multiplying the activity score with the
# Hi-C interaction frequency.

# function to calculate geometric mean
gm_mean <- function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# compute enhancer activity score
chrom_annot_etps <- chrom_annot_etps %>%
  rowwise() %>%
  mutate(activity_score = gm_mean(c(`Dnase-seq`, H3K27ac))) %>%
  ungroup()

# compute abc score
chrom_annot_etps <- chrom_annot_etps %>%
  group_by(gene) %>%
  mutate(abc_score = activity_score * (int_freq + 1),
         abc_score_norm = abc_score / sum(abc_score))

# Finalize output ==================================================================================

# The output is finalized by adding the average gene expression, transforming the distance to TSS to
# absolute values and pruning pairs based on the distance to TSS.

# add average gene expression
chrom_annot_etps <- left_join(chrom_annot_etps, avg_expr, by = "gene")

# transform distance to TSS
chrom_annot_etps <- mutate(chrom_annot_etps, dist_to_tss = abs(dist_to_tss))

# only retain pairs with distance to TSS between 1 and 300 kb
chrom_annot_etps <- filter(chrom_annot_etps, dist_to_tss >= 1000, dist_to_tss <= 300000)

# rearrange columns
columns <- c("perturbation", "gene", "avg_txs", "enh_chr", "enh_start", "enh_end", "gene_chr",
             "gene_tss", "gene_strand", "dist_to_tss", "Dnase-seq", "H3K27ac", "H3K27me3",
             "H3K4me1", "H3K4me3", "POLR2A", "dist_bins", "int_freq", "activity_score",
             "abc_score", "abc_score_norm")
chrom_annot_etps <- chrom_annot_etps[, columns]

# write to output file
write.csv(chrom_annot_etps, file = here(snakemake@output[[1]]), quote = FALSE, row.names = FALSE)
