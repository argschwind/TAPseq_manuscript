## run all workflows for the project

# config file containing samples and parameters
configfile: "config.yml"

# import workflows
include: "rules/download_fastq.smk"
include: "rules/create_alignment_refs.smk"
include: "rules/align_reads.smk"
include: "rules/extract_dge.smk"
include: "rules/analyses.smk"
include: "rules/chromatin_annotated_etps.smk"

# ALL RULE -----------------------------------------------------------------------------------------

# run whole workflow to align reads of all samples and extract dge for specified number of cells
# (samples and cells are read from config file)
rule all:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["samples"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"]),
    analyses = ["results/analyses/tapseq_vs_cropseq.html", "results/analyses/screen_data_qc.html",
      "results/analyses/compare_covariates.html", "results/analyses/map_enhancers.html",
      "results/analyses/chromatin_analyses_screen.html", "results/analyses/hic_analysis.html",
      "results/analyses/abc_analysis.html"],
    ds_dge = [expand("data/Sample10X/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["Sample10X"]),
      expand("data/11iv210ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["11iv210ng"]),
      expand("data/11iv22ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["11iv22ng"]),
      expand("data/8iv210ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["8iv210ng"]),
      expand("data/8iv22ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["8iv22ng"]),
      expand("data/wtxmmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["wtxmmix"]),
      expand("data/wtxlung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["wtxlung"]),
      expand("data/tapmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["tapmix"]),
      expand("data/taplung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["taplung"]),
      expand("data/taphumanmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["taphumanmix"]),
      expand("data/perturbchr8v2/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbchr8v2"]),
      expand("data/perturbchr8alt1/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbchr8alt1"]),
      expand("data/perturbchr8alt2/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbchr8alt2"]),
      expand("data/perturbL1000/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbL1000"]),
      expand("data/tapk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["tapk562deep"]),
      expand("data/wtxk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["wtxk562deep"]),
      expand("data/W4ea/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["W4ea"]),
      expand("data/T4ea/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["T4ea"]),
      expand("data/TAP1/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["TAP1"]),
      expand("data/TAP2/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["TAP2"]),
      expand("data/WholeTx/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["WholeTx"]),
      expand("data/TAPtotalBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["TAPtotalBM"]),
      expand("data/TAPkitBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["TAPkitBM"]),
      expand("data/WholeTotalBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["WholeTotalBM"]),
      expand("data/WholeKitBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["WholeKitBM"]),
      expand("data/taplung/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsampled_150_genes"]["reads_per_cell"]["taplung"]),
      expand("data/tapmix/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsampled_150_genes"]["reads_per_cell"]["tapmix"])]

# RULES TO RUN WORKFLOW PARTIALLY ------------------------------------------------------------------

# all input for figure 1 vignettes
rule Figure1:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["figure1"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["figure1"]),
    ds_dge = [expand("data/tapk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["tapk562deep"]),
              expand("data/wtxk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["wtxk562deep"])]

# all input for figure 2 vignettes              
rule Figure2:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["figure2"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["figure2"]),
    ds_dge = [expand("data/Sample10X/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["Sample10X"]),
              expand("data/TAP1/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
                rpc = config["downsample"]["reads_per_cell"]["TAP1"]),
              expand("data/TAP2/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
                rpc = config["downsample"]["reads_per_cell"]["TAP2"]),
              expand("data/WholeTx/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
                rpc = config["downsample"]["reads_per_cell"]["WholeTx"])]
                
# all enhancer analyses
rule enhancer_screen:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["screen"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["screen"]),
    analyses = ["results/analyses/screen_data_qc.html",
                "results/analyses/compare_covariates.html",
                "results/analyses/map_enhancers.html",
                "results/analyses/chromatin_analyses_screen.html",
                "results/analyses/hic_analysis.html",
                "results/analyses/abc_analysis.html"]

# create all data for mouse bone marrow analyses
rule bone_marrow_cell_types:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["bone_marrow"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["bone_marrow"]),
    ds_dge = [expand("data/TAPtotalBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["TAPtotalBM"]),
      expand("data/TAPkitBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["TAPkitBM"]),
      expand("data/WholeTotalBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["WholeTotalBM"]),
      expand("data/WholeKitBM/downsampled/dge_{rpc}_avg_reads_per_cell_onGenes.txt.gz",
        rpc = config["downsample"]["reads_per_cell"]["WholeKitBM"])]
      
# additional validation data and input
rule tapseq_validation:
  input:
    align = expand("results/alignment/{sample}_align_report.html",
      sample = config["tapseq_validation"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["tapseq_validation"]),
    analyses = "results/analyses/tapseq_vs_cropseq.html",
    ds_dge = [expand("data/Sample10X/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["Sample10X"]),
      expand("data/11iv210ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["11iv210ng"]),
      expand("data/11iv22ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["11iv22ng"]),
      expand("data/8iv210ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["8iv210ng"]),
      expand("data/8iv22ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["8iv22ng"]),
      expand("data/wtxmmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["wtxmmix"]),
      expand("data/wtxlung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["wtxlung"]),
      expand("data/tapmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["tapmix"]),
      expand("data/taplung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["taplung"]),
      expand("data/taphumanmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["taphumanmix"]),
      expand("data/perturbchr8v2/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbchr8v2"]),
      expand("data/perturbchr8alt1/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbchr8alt1"]),
      expand("data/perturbchr8alt2/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbchr8alt2"]),
      expand("data/perturbL1000/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["perturbL1000"]),
      expand("data/tapk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["tapk562deep"]),
      expand("data/wtxk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["wtxk562deep"]),
      expand("data/W4ea/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["W4ea"]),
      expand("data/T4ea/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsample"]["reads_per_cell"]["T4ea"]),
      expand("data/taplung/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsampled_150_genes"]["reads_per_cell"]["taplung"]),
      expand("data/tapmix/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
        rpc = config["downsampled_150_genes"]["reads_per_cell"]["tapmix"])]
   
# create whole-genome and tap-seq alignment references
rule alignment_references:
  input:
    expand("data/alignment_references/{align_ref}/{ref_file}",
      align_ref = ["hg38_genome_ref", "hg38_tapseq_ref_validation", "hg38_tapseq_ref_chr8_screen",
        "hg38_tapseq_ref_chr11_screen", "hg38_genome_ref_v2", "hg38_tapseq_ref_validation_v2",
        "hg38_tapseq_ref_validation_v3", "hg38_genome_ref_dropseq", "hg38_tapseq_ref_dropseq",
        "mm10_genome_ref", "hg38_genome_ref_rev", "hg38_tapseq_ref_rev", "mm10_tapseq_ref_mix",
        "hg38_genome_ref_mix", "hg38_tapseq_ref_mix", "hg38_tapseq_ref_l1000"],
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
