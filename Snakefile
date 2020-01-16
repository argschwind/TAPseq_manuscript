## run all workflows for the project

# config file containing samples and parameters
configfile: "config.yml"

# import workflows
include: "rules/create_alignment_refs.smk"
include: "rules/align_reads.smk"
include: "rules/extract_dge.smk"
include: "rules/analyses.smk"
include: "rules/chromatin_annotated_etps.smk"
include: "rules/revision_analyses.smk"
include: "rules/process_rnaseq.smk"

# ALL RULE -----------------------------------------------------------------------------------------

# run whole workflow to align reads of all samples and extract dge for specified number of cells
# (samples and cells are read from config file)
rule all:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["samples"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"]),
    analyses = ["results/analyses/tapseq_vs_cropseq.html",
                "results/analyses/screen_data_qc.html",
                "results/analyses/compare_covariates.html",
                "results/analyses/map_enhancers.html",
                "results/analyses/chromatin_analyses_screen.html",
                "results/analyses/hic_analysis.html",
                "results/analyses/abc_analysis.html"]

# RULES TO RUN WORKFLOW PARTIALLY ------------------------------------------------------------------

# all validation analyses
rule tapseq_validation:
  input:
    analyses ="results/analyses/tapseq_vs_cropseq.html",
    ds_dge = [expand("data/re-wholetx/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["re-wholetx"]),
              expand("data/re11iv210ng/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["re11iv210ng"]),
              expand("data/W4ea/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["W4ea"]),
              expand("data/T4ea/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["T4ea"]),
              expand("data/wtxmmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["wtxmmix"]),
              expand("data/wtxlung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["wtxlung"]),
              expand("data/tapmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["tapmix"]),
              expand("data/taplung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["taplung"]),
              expand("data/19s005246/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["19s005246"]),
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
              expand("data/perturbL1000a/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["perturbL1000a"]),
              expand("data/perturbL1000b/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["perturbL1000b"]),
              expand("data/perturb10xWTX/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["perturb10xWTX"]),
              expand("data/tapk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["tapk562deep"]),
              expand("data/wtxk562deep/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = config["downsample"]["reads_per_cell"]["wtxk562deep"]),
              expand("data/taplung/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = [1000, 1500, 3500, 5000, 10000, 15000, 20000, 40000]),
              expand("data/tapmix/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
                rpc = [1000, 1500, 3500, 5000, 10000, 15000, 20000, 40000, 80000, 150000])]
                
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
    ds_dge = [expand("data/{sample}/downsampled/dge_{reads}_avg_reads_per_cell_onGenes.txt",
        sample = ["TAPtotalBM", "TAPkitBM"],
        reads = [100, 500, 1000, 1500, 2000, 2500, 5000, 5500]),
      expand("data/{sample}/downsampled/dge_{reads}_avg_reads_per_cell_onGenes.txt",
        sample = ["WholeTotalBM", "WholeKitBM"],
        reads = [100, 500, 1000, 1500, 2000, 2500, 5000, 5500, 10000, 20000, 30000]),
      "data/WholeTotalBM/downsampled/dge_50000_avg_reads_per_cell_onGenes.txt"]
    
    
    
    
    

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
    
# process all bulk RNA-seq samples
rule rnaseq:
  input:
    expand("data/bulk_rnaseq/{sample}/abundance.{format}", sample = config["rnaseq_samples"],
      format = ["h5", "tsv"])
