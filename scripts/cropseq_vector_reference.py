#!/usr/bin/env python

# Construct fasta sequence files and gtf annotations for vector transcripts
# that will be added to the alignment reference data

import argparse, glob

# define functions to generate reference data ----------------------------------

# function to generate alignment reference data for a list of fasta files
# containing crop-seq vector sequences
def cropseq_alignment_reference(input_dir, output_bn, fasta_ext = ".fasta",
                                fasta_line = 60, prefix = ""):

  # get all fasta files in input directory
  infiles = glob.glob(input_dir + "/*" + fasta_ext)

  # output files
  fasta_outfile = output_bn + fasta_ext
  gtf_outfile = output_bn + ".gtf"

  # open output files
  output_fasta = open(fasta_outfile, "w")
  output_gtf  = open(gtf_outfile,  "w")

  # process each input sequence file
  for file in infiles:

    with open(file, "r") as infile:

  	  fasta_entry = infile.readlines()  # read file

  	  # create gtf entry
  	  gtf = gtf_entry(fasta_entry, prefix = prefix)

  	  # reformat fasta entry
  	  fasta = format_fasta(fasta_entry, line = fasta_line, prefix = prefix)

  	  # add gtf and fasta to respective output files
  	  output_gtf.write("%s\n" % gtf)

  	  for i in fasta:
  	    output_fasta.write("%s\n" % i)
  	    
  # close open files
  output_fasta.close()
  output_gtf.close()

# function to create gtf entry from fasta sequence
def gtf_entry(fasta_entry, prefix = ""):

  # create sequence name
  seq_name = prefix + fasta_entry[0][1:].strip()

  # get sequence
  seq = list(map(lambda x: x.rstrip(), fasta_entry[1:]))
  seq = "".join(seq)

  # create gtf attribute field
  attr_names = ["gene_id", "transcript_id", "gene_name", "transcript_name"]
  attr = [s + " " + seq_name for s in attr_names]
  attr = "; ".join(attr) + ";"

  # create gtf entry
  gtf = [seq_name, "VECTOR", "exon", "1", str(len(seq)), ".", "+", ".", attr]
  gtf = "\t".join(gtf)

  return gtf

# function to format fasta sequence entry
def format_fasta(fasta_entry, line = 60, prefix = ""):

  # create sequence name
  seq_id = ">" + prefix + fasta_entry[0][1:].strip()

  # get sequence
  seq = list(map(lambda x: x.rstrip(), fasta_entry[1:]))
  seq = "".join(seq)

  # split seq into array of strings with correct line length
  chunks = list(range(0, len(seq), line))
  seq = [seq[i: i + line] for i in chunks]

  # create new fasta entry
  fasta = [seq_id] + seq

  return fasta

# create reference files -------------------------------------------------------

# create crop-seq vector references based on command line arguments if the
# is called as main program
if __name__ == "__main__":

  # parse command line arguments
  parser = argparse.ArgumentParser(description = ("Create CROP-seq vector "
                                    "alignment references"))
  parser.add_argument("-i", "--input_dir", type = str, required = True,
                      help = ("input directory containing sequences of "
                        "CROP-seq vectors in separate fasta files."))
  parser.add_argument("-o", "--output_bn", type = str, required = True,
                      help = "basename for output files.")
  parser.add_argument("--fasta_ext", type = str, default = ".fasta",
                      help = "filename extension of fasta files "
                      "(default: .fasta)")
  parser.add_argument("--fasta_line", type = int, default = 60,
                      help = "line length of output fasta file (default : 60)")
  args = parser.parse_args()

  # create crop-seq vector reference
  cropseq_alignment_reference(input_dir = args.input_dir,
                              output_bn = args.output_bn,
                              fasta_ext = args.fasta_ext,
                              fasta_line = args.fasta_line)
