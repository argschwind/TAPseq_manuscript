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
                "results/screen_data_qc.html",
                "results/compare_covariates.html",
                "results/map_enhancers.html",
                "results/chromatin_analyses_screen.html",
                "results/hic_analysis.html",
                "results/abc_analysis.html"]

# functions to run workflow only partially ---------------------------------------------------------

# create whole-genome and tap-seq alignment references
rule alignment_references:
  input:
    expand("data/alignment_references/{align_ref}/{ref_file}",
      align_ref = ["hg38_genome_ref", "hg38_tapseq_ref_validation", "hg38_tapseq_ref_chr8_screen",
        "hg38_tapseq_ref_chr11_screen", "hg38_genome_ref_v2", "hg38_tapseq_ref_validation_v2",
        "hg38_tapseq_ref_validation_v3", "hg38_genome_ref_dropseq", "hg38_tapseq_ref_dropseq",
        "mm10_genome_ref", "hg38_genome_ref_rev", "hg38_tapseq_ref_rev", "mm10_tapseq_ref_mix",
        "hg38_tapseq_ref_l1000"],
      ref_file = ["cropseq_ref.dict", "cropseq_ref.refFlat", "genomeDir"])

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
