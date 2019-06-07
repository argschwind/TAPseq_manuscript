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
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"]),
    analyses = ["results/tapseq_vs_cropseq.html",
                "results/downsampled_dge.html",
                "results/downsampled_target_reads.html",
                "results/8iScreen1_qc.html",
                "results/11iScreen1_qc.html",
                "results/compare_covariates.html",
                "results/map_enhancers.html",
                "results/chromatin_analyses_screen.html"]

# functions to run workflow only partially ---------------------------------------------------------

# create whole-genome and tap-seq alignment references
rule alignment_reference:
  input:
    "data/genome_reference/cropseq_ref.dict",
    "data/genome_reference/cropseq_ref.refFlat",
    "data/genome_reference/genomeDir",
    "data/tapseq_ref_validation/cropseq_ref.dict",
    "data/tapseq_ref_validation/cropseq_ref.refFlat",
    "data/tapseq_ref_validation/genomeDir",
    "data/tapseq_ref_chr8_screen/cropseq_ref.dict",
    "data/tapseq_ref_chr8_screen/cropseq_ref.refFlat",
    "data/tapseq_ref_chr8_screen/genomeDir",
    "data/tapseq_ref_chr11_screen/cropseq_ref.dict",
    "data/tapseq_ref_chr11_screen/cropseq_ref.refFlat",
    "data/tapseq_ref_chr11_screen/genomeDir",
    "data/genome_reference_v2/cropseq_ref.dict",
    "data/genome_reference_v2/cropseq_ref.refFlat",
    "data/genome_reference_v2/genomeDir",
    "data/tapseq_ref_validation_v2/cropseq_ref.dict",
    "data/tapseq_ref_validation_v2/cropseq_ref.refFlat",
    "data/tapseq_ref_validation_v2/genomeDir"

# run workflow until read alignment. this allows to determine the number of cells before extracting
# dge data and finishing the workflow
rule align:
  input:
    expand("results/alignment/{sample}_align_report.html", sample = config["samples"])
    
# run workflow until digital gene expression (DGE) extraction. this concludes the data processing
# part of the project
rule dge:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["samples"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"])
