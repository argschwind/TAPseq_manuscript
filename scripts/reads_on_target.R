# simple script to count the number of total reads and reads on target for chr8 and chr11 crop-seq
# validation samples. this script contains a lot of hard-coded things like panel and screen
# identifiers and sample naming patterns, so be careful when changing upstream data!

library(tidyverse)
library(here)

# get tap-seq target genes -------------------------------------------------------------------------

# load tap-seq target genes
target_genes <- read.csv(here(snakemake@input$target_genes), stringsAsFactors = FALSE)

# known enhancer target genes
known_enh <- target_genes %>%
  filter(screen == "validation") %>%
  pull(gene)

# chromosome 11 inhibition panel
chr11_genes <- target_genes %>%
  filter(screen == "inhibition", panel == "chr11_hs2") %>%
  pull(gene) %>%
  c(known_enh)

# chromosome 8 inhibition panel
chr8_genes <- target_genes %>%
  filter(screen == "inhibition") %>%
  filter(panel %in% c("chr8_zfpm2", "chr8_myc")) %>%
  pull(gene) %>%
  c(known_enh)

# string patterns for all target genes (incl CROP-seq vectors)
vector_annot <- snakemake@params$vector_annot
targets_chr11 <- c(chr11_genes, vector_annot)
targets_chr8 <- c(chr8_genes, vector_annot)

# total number of reads ----------------------------------------------------------------------------

# fastq files containing genome reads of each sample
fastq_files <- snakemake@input$fastq

# infer sample names from fastq files and set as names
bn <- sub(".*lane1", "", basename(fastq_files))
names(fastq_files) <- sub("_2_sequence.txt.gz", "", bn)

# count total reads in fastq files
raw_reads <- sapply(X = fastq_files, FUN = function(x) {
  command <- paste("echo $(zcat", x, "| wc -l) / 4 | bc")
  as.numeric(system(command, intern = TRUE))
})

# convert to data.frame
raw_reads <- data.frame(sample = names(raw_reads), total_reads = raw_reads,
                        stringsAsFactors = FALSE)

# reads on target ----------------------------------------------------------------------------------

# function to count reads on target genes
reads_on_target <- function(bam_files, targets) {
  
  # write target genes into text file for grep
  write(targets, file = "./targets.txt")
  
  # count reads on targets
  target_reads <- sapply(bam_files, FUN = function(x) {
    command <- paste("samtools view", x, "| fgrep -f", "./targets.txt", "| wc -l")
    as.numeric(system(command, intern = TRUE))
  })
  
  # delete target genes tempfile
  unlink("./targets.txt")
  
  # create final output
  data.frame(sample = names(bam_files), target_reads = target_reads, stringsAsFactors = FALSE)
  
}

# get bam files
bam_files <- here(snakemake@input$bam)
names(bam_files) <- basename(dirname(bam_files))

# bam files for chr11 and chr8 targets
chr11_bam <- c(bam_files[grep(names(bam_files), pattern = "^11")], bam_files["Sample10X"])
chr8_bam <- c(bam_files[grep(names(bam_files), pattern = "^8")], bam_files["Sample10X"])

# count reads mapping to target genes
chr11_target_reads <- reads_on_target(chr11_bam, targets = targets_chr11)
chr8_target_reads <- reads_on_target(chr8_bam, targets = targets_chr8)

# finalize output and save to file -----------------------------------------------------------------

# combine read counts into data frame and compute percentage of reads on target genes
read_counts <- bind_rows(chr11 = chr11_target_reads, chr8 = chr8_target_reads, .id = "panel") %>%
  left_join(raw_reads, by = "sample") %>%
  mutate(perc = target_reads / total_reads)

# save to output file
write.csv(read_counts, file = here(snakemake@output), row.names = FALSE)
