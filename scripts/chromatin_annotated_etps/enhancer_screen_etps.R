## Create chromatin annotated ETPs for the enhancer screen

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

# load processed differential expression output (includes confidence level and distance to TSS)
results <- here(snakemake@input$processed_results) %>%
  read.csv(stringsAsFactors = FALSE) %>%
  mutate(sample = paste0("chr", sub("iScreen1", "", sample)))

# extract cis-interactions on chromosome 8 and 11 and add column identifying significant hits
cis_enh_perts <- results %>%
  filter(enh_type == "cis", abs(dist_to_tss) >= 1000) %>%
  mutate(enh_chr = paste0("chr", enh_chr), gene_chr = paste0("chr", gene_chr),
         significant = if_else(pval_adj_allTests < 0.05, true = "sig", false = "non_sig"))

# remove validation controls on other chromosomes than the samples' target region and significant
# hits that increase gene expression
cis_enh_perts <- cis_enh_perts %>%
  filter(sample == enh_chr) %>%
  filter(!(significant == "sig" & manual_lfc > 0))

# Annotate chromatin ===============================================================================

# Every assessed enhancer is annotated for the number of chromatin assay ChIP-seq and DNAse-seq
# reads. The RPKM is calculated to normalize for both the sequencing depth and enhancer size.

# extract hg38 enhancer coordinates
enh_coords <- cis_enh_perts %>%
  select(enh_chr, enh_start, enh_end, perturbation) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# get coordinates of genomic regions and extend by 10 kb
region_coords <- range(enh_coords)
start(region_coords) <- start(region_coords) - 10000
end(region_coords) <- end(region_coords) + 10000

# bam files containing aligned reads
bam_files <- here(snakemake@input$encode_bam)
names(bam_files) <- sub("_encode_rep1_alignment.bam", "", basename(bam_files))

# load aligned reads within the regions from bam files
chrom_reads <- bplapply(bam_files, FUN = readGAlignments,
                        param = ScanBamParam(which = region_coords))

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
cis_enh_perts <- enh_chrom_norm %>%
  select(-c(enh_chr, enh_start, enh_end, width, reads)) %>%
  spread(key = "assay", value = "rpkm") %>%
  left_join(x = cis_enh_perts, y = ., by = "perturbation")

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

# extract hg38 enhancer coordinates
enh_coords_hg38 <- cis_enh_perts %>%
  select(perturbation, enh_chr, enh_start, enh_end) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "enh_chr",
                           start.field = "enh_start", end.field = "enh_end")

# liftover enhancers from hg38 to hg19 and convert to data.frame
enh_coords_hg19 <- enh_coords_hg38 %>%
  liftOver(chain = hg38_to_hg19_chain) %>%
  unlist() %>%
  as.data.frame() %>%
  select(seqnames, start, end, perturbation) %>%
  rename(enh_chr = seqnames, enh_start = start, enh_end = end)

# replace tested enhancer - gene pairs hg38 enhancer coordinates with hg19 coordinates
cis_enh_perts <- cis_enh_perts %>%
  select(-c(enh_chr, enh_start, enh_end)) %>%
  left_join(enh_coords_hg19, by = "perturbation")

# extract gene TSS coordinates
gene_tss_hg38 <- cis_enh_perts %>%
  select(gene, gene_chr, gene_tss, gene_strand) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "gene_chr",
                           start.field = "gene_tss", end.field = "gene_tss",
                           strand.field = "gene_strand")

# liftover tss coordinates to hg19
gene_tss_hg19 <- gene_tss_hg38 %>%
  liftOver(chain = hg38_to_hg19_chain) %>%
  unlist() %>%
  as.data.frame() %>%
  select(seqnames, start, strand, gene) %>%
  rename(gene_chr = seqnames, gene_tss = start, gene_strand = strand)

# replace gene tss hg38 coordinates with hg19 coordinates
cis_enh_perts <- cis_enh_perts %>%
  select(-c(gene_chr, gene_tss, gene_strand)) %>%
  left_join(gene_tss_hg19, by = "gene")

# recalculate distance to tss for hg19 coordinates
cis_enh_perts <- cis_enh_perts %>%
  mutate(enh_center = round((enh_start + enh_end) / 2)) %>%
  mutate(dist_to_tss = if_else(gene_strand == "+", true  = enh_center - gene_tss,
                               false = gene_tss - enh_center)) %>%
  select(-enh_center)

# Prepare Hi-C data --------------------------------------------------------------------------------

# Contact matrices for K562 cells are imported (5kb resolution). The provided normalization vectors
# are used to normalize the observed contacts (read counts). Data from the two chromosonal regions
# on chromosome 8 and 11 is extracted for any further analyses.

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

# k562 intrachromosomal sparse matrix and krnorm vector files for chromosomes 8 and 11
scm_files <- here(snakemake@input$hic_raw)
names(scm_files) <- sub("_5kb.RAWobserved", "", basename(scm_files))

krnorm_files <- here(snakemake@input$hic_norm)
names(krnorm_files) <- sub("_5kb.KRnorm", "", basename(krnorm_files))

# import normalization vectors
chr8_KRnorm  <- as.numeric(readLines(krnorm_files["chr8"]))
chr11_KRnorm <- as.numeric(readLines(krnorm_files["chr11"]))

# infer number of bins per chromosome based on the normalization vectors
chr8_bins  <- length(chr8_KRnorm)
chr11_bins <- length(chr11_KRnorm)

# import hi-c data for these chromosomes
chr8_hic  <- import_hic(scm_files["chr8"],  chromosome = "chr8", resolution = 5000,
                        bins = chr8_bins)
chr11_hic <- import_hic(scm_files["chr11"], chromosome = "chr11", resolution = 5000,
                        bins = chr11_bins)

# function to normalize Hi-C data based on provided normalization vectors
normalize_hic <- function(htc_obj, norm_vector) {
  
  # extract raw observed interaction matrix
  raw_obs <- intdata(htc_obj)
  
  # create normalization matrix by computing the outer product of the normalization vector
  norm_mat <- outer(norm_vector, norm_vector)
  
  # multiply observed interactions by normalization matrix and add back to HTC object
  intdata(htc_obj) <- raw_obs / norm_mat
  return(htc_obj)
  
}

# normalize HiC data
chr8_hic_norm  <- normalize_hic(chr8_hic,  norm_vector = chr8_KRnorm)
chr11_hic_norm <- normalize_hic(chr11_hic, norm_vector = chr11_KRnorm)

# infer chromosomal region range
region_coords <- cis_enh_perts %>%
  select(sample, enh_start, enh_end, gene_tss) %>%
  gather(key = "key", value = "coord", -sample) %>%
  group_by(sample) %>%
  summarize(start = min(coord), end = max(coord))

# calculate bin coordinates that contain the entire regions
region_bins <- region_coords %>%
  mutate(start = floor(start / 5000) * 5000,
         end = ceiling(end / 5000) * 5000)

# extract data for assessed regions
chr8_region_hic <- extractRegion(chr8_hic_norm, MARGIN = c(1, 2), chr = "chr8",
                                 from = pull(filter(region_bins, sample == "chr8"), start),
                                 to = pull(filter(region_bins, sample == "chr8"), end))

chr11_region_hic <- extractRegion(chr11_hic_norm, MARGIN = c(1, 2), chr = "chr11",
                                  from = pull(filter(region_bins, sample == "chr11"), start),
                                  to = pull(filter(region_bins, sample == "chr11"), end))

# Enhancer - gene pairs ----------------------------------------------------------------------------

# The interaction frequency of all tested enhancer - gene pair computed, defined as the interaction
# frequency of the genomic bins overlapping the enhancer and the target genes TSS.

# function to compute the interaction frequency for enhancer - gene pairs, by finding the hi-c 
# genomic bins with overlap with the enhancer and the target gene tss. the interaction frequency of
# the pair is then defined as the interaction frequency of these bins
compute_int_freq <- function(pairs, htc_object) {
  
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
  
  # get interaction frequencies for all bins and add pair id
  int_freq_pairs <- intdata[as.matrix(bins[, 2:3])]
  int_freq_pairs <- data.frame(pair = bins$pair, int_freq = int_freq_pairs)
  
  # add interaction frequencies and bin distances to pairs to create output
  pairs %>%
    left_join(dist_pairs, by = "pair") %>%
    left_join(int_freq_pairs, by = "pair") %>%
    select(-c(pair))
  
}

# compute interaction frequencies for all tested enhancer - gene pairs
chr8_pairs <- cis_enh_perts %>%
  filter(enh_chr == "chr8", enh_chr == "chr8") %>%
  compute_int_freq(., htc_object = chr8_region_hic)

chr11_pairs <- cis_enh_perts %>%
  filter(enh_chr == "chr11", enh_chr == "chr11") %>%
  compute_int_freq(., htc_object = chr11_region_hic)

# combine into one data.frame
pairs <- rbind(chr8_pairs, chr11_pairs)

# Activity by contact ==============================================================================

# The activity by contact score is computed as proposed by Fulco et al., 2019

# function to calculate geometric mean
gm_mean <- function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# compute enhancer activity score
pairs <- pairs %>%
  rowwise() %>%
  mutate(activity_score = gm_mean(c(`Dnase-seq`, H3K27ac)))

# compute abc score
pairs <- pairs %>%
  ungroup() %>% 
  group_by(sample, gene) %>%
  mutate(abc_score = activity_score * (int_freq + 1),
         abc_score_norm = abc_score / sum(abc_score))

# convert distance to TSS to absolute value and change label of significant hits
pairs <- pairs %>%
  mutate(dist_to_tss = abs(dist_to_tss),
         significant = if_else(significant == "sig", true = 1, false = 0))

# write to output file
write.csv(pairs, file = here(snakemake@output[[1]]), quote = FALSE, row.names = FALSE)
