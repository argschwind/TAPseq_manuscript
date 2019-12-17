#!/usr/bin/env python

# construct fasta sequence files and gtf annotations for vector transcripts that will be added to
# the alignment reference data

import argparse
from Bio import SeqIO

# define functions to generate reference data ------------------------------------------------------

# function to generate alignment reference data from a fasta file containing vector sequences
def cropseq_alignment_reference(input_fasta, output_bn, fasta_ext = ".fasta", prefix = ""):

  # output files
  fasta_outfile = output_bn + fasta_ext
  gtf_outfile = output_bn + ".gtf"

  # open output files
  output_fasta = open(fasta_outfile, "w")
  output_gtf  = open(gtf_outfile,  "w")

  # process each input sequence
  for record in SeqIO.parse(input_fasta, "fasta"):
    
    # add prefix to sequence id
    record.id = prefix + record.id
    
    # remove name and description to remove them from header in fasta output
    record.name = ""
    record.description = ""
    
    # create gtf entry
    gtf = gtf_entry(record)
    
    # write gtf entry and modified fasta record to respective output files
    output_gtf.write("%s\n" % gtf)
    SeqIO.write(record, output_fasta, "fasta")
 
  # close open files
  output_fasta.close()
  output_gtf.close()

# function to create gtf entry from a crop-seq vector fasta record
def gtf_entry(fasta_record):
  
  # get sequence name and length
  seq_name = fasta_record.id
  seq_len = len(fasta_record.seq)

  # create gtf attribute field
  attr_names = ["gene_id", "transcript_id", "gene_name", "transcript_name"]
  attr = [s + " " + seq_name for s in attr_names]
  attr = "; ".join(attr) + ";"

  # create gtf entry
  gtf_fields = [seq_name, "VECTOR", "exon", "1", str(seq_len), ".", "+", ".", attr]
  gtf_line = "\t".join(gtf_fields)

  return gtf_line

# create reference files ---------------------------------------------------------------------------

# create crop-seq vector references based on command line arguments if the  is called as main
# program
if __name__ == "__main__":

  # parse command line arguments
  parser = argparse.ArgumentParser(description = ("Create CROP-seq vector "
                                    "alignment references"))
  parser.add_argument("-i", "--input_fasta", type = str, required = True,
                      help = ("Input fasta file containing sequences of CROP-seq vectors."))
  parser.add_argument("-o", "--output_bn", type = str, required = True,
                      help = "Basename for output files.")
  parser.add_argument("--fasta_ext", type = str, default = ".fasta",
                      help = "Filename extension of fasta files "
                      "(default: .fasta).")
  parser.add_argument("--prefix", type = str, default = "",
                      help = "Optional prefix to be added to sequence name.")
  args = parser.parse_args()

  # create crop-seq vector reference
  cropseq_alignment_reference(input_fasta = args.input_fasta,
                              output_bn = args.output_bn,
                              fasta_ext = args.fasta_ext,
                              prefix = args.prefix)
