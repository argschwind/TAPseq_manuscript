## run all workflows for the project

# config file containing samples and parameters
configfile: "config.yml"

# import workflows
include: "rules/create_alignment_refs.smk"
include: "rules/align_reads.smk"
include: "rules/extract_dge.smk"
include: "rules/analyses.smk"

# ALL FUNCTION -------------------------------------------------------------------------------------

# run whole workflow to align reads of all samples and extract dge for specified number of cells
# (samples and cells are read from config file)
rule all:
  input:
   align = expand("results/alignment/{sample}_align_report.html", sample = config["samples"]),
   dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"])

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
    expand("results/alignment/{sample}_align_report.html", sample = config["samples"])
