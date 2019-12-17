#!/usr/bin/env Rscript

# script to create custom tap-seq annotation data for read alignment

# parse command line arguments =====================================================================

library(optparse)

# create arguments list
option_list = list(
  make_option(c("-t", "--txs_file"), type = "character", default = NULL,
              help = "Path to transcripts annotation file (.gtf)",
              metavar = "character"),
  make_option(c("-b", "--BSgenome_id"), type = "character", default = NULL,
              help = "Identifier for used BSgenome", metavar = "character"),
  make_option(c("-f", "--fasta_outfile"), type = "character", default = NULL,
              help = "Path to .fasta output file", metavar = "character"),
  make_option(c("-o", "--gtf_outfile"), type = "character", default = NULL,
              help = "Path to .gtf output file", metavar = "character")
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
  }
}

# check that all required parameters are provided
required_args <- c("txs_file", "BSgenome_id", "fasta_outfile", "gtf_outfile")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

# define functions =================================================================================

# attach required packages
library(rtracklayer)
library(BSgenome)
library(GenomicRanges)
library(GenomicFeatures)

# main function to create gene loci tap-seq alignment reference based on input file paths
create_tapseq_ref <- function(txs_file, BSgenome_id, fasta_outfile, gtf_outfile) {
  
  message("\nLoading data...")
  
  # load transcripts
  txs <- import(txs_file, format ="gtf")
  
  # get BSgenome object containing the genome sequence
  genome <- getBSgenome(BSgenome_id)
  seqlevelsStyle(genome) <- seqlevelsStyle(txs)[1]  # set seqlevels to same style as txs input
  
  # create tap-seq alignment references
  message("Creating tap-seq alignment reference...")
  ref <- tapseq_alignment_ref(txs, genome = genome)
  
  # save alignment reference data
  writeXStringSet(ref$contig_seqs, filepath = fasta_outfile, format = "fasta")
  export_gtf(ref$annot, filepath = gtf_outfile)

  message("Done!")
  
}

# helper functions to create tap-seq alignment reference data ======================================

# create alignment reference data for gene loci
tapseq_alignment_ref <- function(txs, genome, gene_id = "gene_name") {
  
  # convert txs to GRanges if it's a GRangesList
  if (is(txs, "GRangesList")) {
    txs <- unlist(txs)
  }
  
  # abort if genome is not a BSgenome or DNAStringSet object
  if (!any(is(genome, "BSgenome") | is(genome, "DNAStringSet"))) {
    stop("genome must be of class BSgenome or DNAStringSet!", call. = FALSE)
  }
  
  # sort txs by coordinates
  txs <- sort(txs)
  
  # group transcripts by gene names
  genes <- split(txs, f = mcols(txs)[, gene_id])
  
  # create contigs for each locus ------------------------------------------------------------------
  
  # obtain ranges of each gene locus (start of gene to end of gene)
  loci <- range(genes)
  
  # merge loci to create non-overlapping contigs
  contigs <- reduce(unlist(loci), ignore.strand = TRUE)
  strand(contigs) <- "+"
  
  # split contigs into GRangesList
  contigs <- split(contigs, f = 1:length(contigs))
  
  # find genes overlapping with each contig
  overlaps <- as.data.frame(findOverlaps(contigs, genes, ignore.strand = TRUE))
  
  # split by query hit to have all overlapping genes per contig in one table
  overlaps <- split(overlaps, f = overlaps$queryHits)
  
  # create annotations for genes overlapping each contig
  annot <- lapply(overlaps, FUN = create_contig_annot, contigs = contigs, genes = genes)
  annot <- unlist(GRangesList(annot))
  names(annot) <- NULL
  
  # get contig names based on names of genes overlapping each contig
  names(contigs) <- vapply(overlaps, FUN = contig_names, genes = genes, FUN.VALUE = character(1))
  
  # get sequences for contigs
  contig_seqs <- extractTranscriptSeqs(genome, transcripts = contigs)
  
  # create and return output
  list(contig_seqs = contig_seqs, annot = annot)
  
}

# create annotations for one contig
create_contig_annot <- function(overlap, contigs, genes) {
  
  # get contig and overlapping genes
  contig <- contigs[[unique(overlap$queryHits)]]
  ol_genes <- unlist(genes[overlap$subjectHits])
  
  # convert genome coordinates of overlapping genes into coordinates relative
  # to the contig
  start(ol_genes) <- start(ol_genes) - start(contig) + 1
  end(ol_genes) <- end(ol_genes) - start(contig) + 1
  
  # create new contig name to use as seqnames of output annotations
  gene_names <- names(genes[overlap$subjectHits])
  contig_name <- paste(gene_names, collapse = "_")
  
  # create new ouptut with contig name as seqnames
  annot <- GRanges(seqnames = contig_name, ranges = ranges(ol_genes), strand = strand(ol_genes))
  
  # add metadata and return output
  mcols(annot) <- mcols(ol_genes)
  annot
  
}

# create contig names by exracting the names of all overlapping genes per contig
contig_names <- function(overlap, genes) {
  
  # get gene names of overlapping genes
  gene_names <- names(genes[overlap$subjectHits])
  
  # create contig name
  paste(gene_names, collapse = "_")
  
}

# export GRanges object to gtf file. rtracklayers 'export()' function causes an error with STAR 
export_gtf <- function(x, filepath) {
  
  # transform x into data.frame
  x <- data.frame(x)
  
  # transform all columns to character
  x[] <- lapply(x, FUN = as.character)
  
  # remove the width column (from converting GRanges to data.frame)
  x <- x[,-which(colnames(x) == "width")]
  
  # replace potential NA for score and phase by '.'
  x$score[is.na(x$score)] <- "."
  x$phase[is.na(x$phase)] <- "."
  
  # replace potential '*' for strand by '.'
  x$strand[x$strand == "*"] <- "."
  
  # get first 8 fields for gtf output
  gtf_fields <- x[,c("seqnames", "source", "type", "start",
                     "end", "score", "strand", "phase")]
  
  ## the remaining fields are processed into the attribute field
  
  # get remaining fields
  attr_cols <- x[,colnames(x)[!colnames(x) %in% colnames(gtf_fields)]]
  
  # get attribute names
  attr_names <- colnames(attr_cols)
  
  # get mandatory attributes "gene_id" and "transcript_id"
  mand <- c(which(attr_names == "gene_id"),
            which(attr_names == "transcript_id"))
  
  # abort if mandatory fields not found
  if (length(mand) != 2) {
    
    stop("Mandatory attributes \"gene_id\" and \"transcript_id\" not found!",
         call. = FALSE)
    
  }
  
  # order attributes so that the mandatory attributes come first
  attr_names <- c(attr_names[mand], attr_names[-mand])
  attr_cols <- attr_cols[attr_names]
  
  # create matrix with names
  attr_names <- matrix(rep(attr_names, nrow(attr_cols)),
                       ncol=length(attr_names), byrow = TRUE)
  
  # paste attribute names to values
  attr <- matrix(paste(attr_names, as.matrix(attr_cols)),
                 nrow = nrow(attr_names))
  
  # collapse all attributes into one string to create the attribute field
  # and add to gtf_fields
  gtf_fields$attributes <- do.call(paste, c(as.data.frame(attr), sep="; "))
  gtf_fields$attributes <- paste0(gtf_fields$attributes, ';')
  
  # write output to file
  write.table(gtf_fields, file = filepath, row.names = FALSE,
              col.names = FALSE, quote = FALSE, sep = "\t") 
  
}

# create tap-seq alignment references ==============================================================

create_tapseq_ref(txs_file = opt$txs_file, BSgenome_id = opt$BSgenome_id,
                  fasta_outfile = opt$fasta_outfile, gtf_outfile = opt$gtf_outfile)
