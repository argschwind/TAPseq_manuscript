## run all workflows for the project

# config file containing samples and parameters
configfile: "config.yml"

# import workflows
include: "rules/create_alignment_refs.smk"
include: "rules/align_reads.smk"
include: "rules/extract_dge.smk"
include: "rules/analyses.smk"

# define required python functions -----------------------------------------------------------------

# create all dge_report output files for a sample id based on the number cells specified for that
# sample under 'dge_ncells' in config file
def dge_report_outfiles(sample):
  ncells = config["dge_ncells"][sample]  # get number of cells for sample
  files = expand("reports/dge/{sample}_{ncells}_cells_dge_report.html",
             sample = sample, ncells = ncells)
  return files
 
# ALL FUNCTION -------------------------------------------------------------------------------------

# run whole workflow to align reads of all samples and extract dge for specified number of cells
# (samples and cells are read from config file)
rule all:
  input:
   align = expand("reports/alignment/{sample}_align_report.html", sample = config["samples"]),
   dge = [dge_report_outfiles(sample) for sample in config["dge_ncells"]]

# functions to run workflow only partially ---------------------------------------------------------

# create whole-genome and tap-seq alignment references
rule alignment_reference:
  input:
    "data/genome_reference/cropseq_ref.dict",
    "data/genome_reference/cropseq_ref.refFlat",
    "data/genome_reference/genomeDir",
    "data/tapseq_reference/cropseq_ref.dict",
    "data/tapseq_reference/cropseq_ref.refFlat",
    "data/tapseq_reference/genomeDir"

# run workflow until read alignment. this allows to determine the number of cells before extracting
# dge data and finishing the workflow
rule align:
  input:
    expand("reports/alignment/{sample}_align_report.html", sample = config["samples"])
