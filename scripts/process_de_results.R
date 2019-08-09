## process differential expression results. calculate FDR across all experiments (samples), compute
# confidence score based on significant single gRNA hits, add manually computed log fold change and
# infer distance to TSS for all cis enhancer - gene pairs

library(here)
library(rtracklayer)
library(tidyverse)

# prepare data -------------------------------------------------------------------------------------

# files containing differential expression results
result_files <- here(snakemake@input$results)

# set name for each file
samples <- basename(dirname(dirname(result_files)))
names(result_files) <- basename(result_files) %>%
  sub("output_MAST_", "", .) %>%
  sub("_.*", "", .) %>%
  paste(samples, ., sep = "_")

# read files and convert into one data.frame
results <- result_files %>%
  lapply(FUN = read.csv, stringsAsFactors = FALSE) %>%
  bind_rows(.id = "id") %>%
  separate(id, into = c("sample", "strategy"), sep = "_")

# files containing manually calculated log fold change
lfc_files <- here(snakemake@input$lfc)
names(lfc_files) <- basename(dirname(dirname(lfc_files)))

# load lfc files and convert into one data.frame
lfc <- lfc_files %>%
  lapply(FUN = read.csv, stringsAsFactors = FALSE) %>%
  bind_rows(.id = "sample")

# replace infinite values (pert or control is 0) by NA
lfc <- mutate(lfc, lfc = replace(lfc, is.infinite(lfc), NA))

# load target gene annotations
annot <- lapply(here(snakemake@input$annot), FUN = import, format = "gtf")

# merge into one GRanges object
annot <- annot %>%
  do.call("c", .) %>%
  unique()

# filter for exons of protein-coding genes (no processed transcripts etc)
annot <- annot[annot$type == "exon" &
                 annot$gene_biotype == "protein_coding" &
                 annot$transcript_biotype == "protein_coding"
               ]

# split by gene name into a GRangesList
genes <- split(annot, f = annot$gene_name)

# load enhancer coordinates files
enh_coords <- here(snakemake@input$enh_coords) %>%
  lapply( FUN = read.table, header = FALSE, stringsAsFactors = FALSE) %>%
  bind_rows() %>%
  distinct()

colnames(enh_coords) <- c("chr", "start", "end", "name", "score", "strand")


# FDR and confidence levels ------------------------------------------------------------------------

# correct for multiple testing across samples using FDR
results <- results %>%
  group_by(strategy) %>%
  rename(pval_adj_sample = pval_adj_allTests) %>%
  mutate(pval_adj_allTests = p.adjust(pvalue, method = "fdr"))

# extract results for 'per enhancer' tests
enh_perts <- filter(results, strategy == "perEnh")

# function to collapse gRNA ids per targeted enhancer
collapse_grnas <- function(grnas) {
  grnas <- sub("_.+", "", grnas)
  sub("\\.[A|B|C|D]", "", grnas)
}

# count number of gRNA hits per enhancer-gene pair
gRNA_hits <- results %>%
  filter(strategy == "perGRNA") %>%
  mutate(perturbation = collapse_grnas(perturbation)) %>%
  group_by(sample, perturbation, gene) %>%
  summarize(grna_hits = sum(pval_adj_allTests < snakemake@params$confidence_fdr),
            prop_grna_hits = mean(pval_adj_allTests < snakemake@params$confidence_fdr))

# add single gRNA hits to enhancer hits
enh_perts <- left_join(enh_perts, gRNA_hits, by = c("sample", "perturbation", "gene"))

# log fold change ----------------------------------------------------------------------------------

# add manually calculate log fold change to de results
enh_perts <- lfc %>%
  select(sample, perturbation, gene, manual_lfc = lfc) %>%
  left_join(enh_perts, ., by = c("sample", "perturbation", "gene"))

# add enhancer and gene TSS coordinates ------------------------------------------------------------

# compute center of every enhancer
enh_centers <- enh_coords %>%
  mutate(enh_center = round((start + end) / 2)) %>%
  select(name, chr, start, end, enh_center) %>%
  mutate(name = gsub(":|-", ".", name),
         chr = sub("chr", "", chr)) %>%
  rename(enh_chr = chr, enh_start = start, enh_end = end)

# add enhancer centers to de results
enh_perts <- enh_perts %>%
  left_join(enh_centers, by = c("perturbation" = "name"))

# get chromosome and strand for each gene
gene_chr <- unlist(unique(seqnames(genes)))
gene_strand <- unlist(unique(strand(genes)))

# function to get a genes TSS
get_tss <- function(x) {
  if (all(strand(x) == "+")) {
    min(start(x))
  }else if (all(strand(x) == "-")) {
    max(end(x))
  }else{
    stop("Inconsistent strand information!")
  }
}

# get TSS for each gene
gene_tss <- sapply(genes, FUN = get_tss)

# create data.frame with strand and tss coordinates for every gene
tss_coords <- data.frame(gene_chr, gene_tss, gene_strand, check.rows = TRUE) %>%
  rownames_to_column(var = "gene") %>%
  mutate_if(is.factor, as.character)

# add to discovery perturbation results
enh_perts <- left_join(enh_perts, tss_coords, by = "gene")


# add distance to TSS for cis associations ---------------------------------------------------------

# add type for cis or trans enhancer interactions
enh_perts <- enh_perts %>%
  mutate(enh_type = if_else(enh_chr == gene_chr, true = "cis", false = "trans"))

# calculate distance to tss for every cis enhancer - gene pair
enh_perts <- enh_perts %>%
  mutate(dist_to_tss = if_else(enh_type == "cis",
                               true = if_else(gene_strand == "+",
                                              true  = enh_center - gene_tss,
                                              false = gene_tss - enh_center),
                               false = as.numeric(NA)))

# reformat for output
enh_perts <- select(enh_perts, -enh_center)

# save processed output to file
write.csv(enh_perts, file = here(snakemake@output), row.names = FALSE)
