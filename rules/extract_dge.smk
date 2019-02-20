## rules to extract DGE data from aligned CROP-seq reads

### input, output and shell paths are all relative to the project directory ###

# workflow rules -----------------------------------------------------------------------------------

# count UMI observations per cell barcode and gene tag
rule umi_observations:
  input:
    "data/{sample}/filt_gene_tagged_aligned.bam"
  output:
    "data/{sample}/umi_observations.txt"
  log:
    "data/{sample}/logs/umi_observations.log"
  params:
    ncells = lambda wildcards: config["dge_ncells"][wildcards.sample],
    edit_distance = config["umi_observations"]["edit_distance"],
    read_mq = config["umi_observations"]["read_mq"],
    min_umi_reads = config["umi_observations"]["min_umi_reads"],
    rare_umi_filter = config["umi_observations"]["rare_umi_filter_threshold"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "GatherMolecularBarcodeDistributionByGene "
    "INPUT={input} "
    "OUTPUT={output} "
    "NUM_CORE_BARCODES={params.ncells} "
    "EDIT_DISTANCE={params.edit_distance} "
    "READ_MQ={params.read_mq} "
    "MIN_BC_READ_THRESHOLD={params.min_umi_reads} "
    "RARE_UMI_FILTER_THRESHOLD={params.rare_umi_filter} "
    "2> {log}"

# extract dge with and filter for chimeric reads
rule extract_dge:
  input:
    "data/{sample}/umi_observations.txt"
  output:
    dge = "data/{sample}/dge.txt",
    dge_stats = "data/{sample}/dge_summary.txt",
    tpt_hist = "data/{sample}/dge_tpt_histogram.txt"
  log:
    "data/{sample}/logs/extract_dge.log"
  params:
    tpt_threshold = config["extract_dge"]["tpt_threshold"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "python scripts/extract_dge.py -i {input} -o {output.dge} "
    "--tpt_threshold {params.tpt_threshold} 2> {log}"
    
# infer perturbation status of each cell
rule perturbation_status:
  input:
    "data/{sample}/dge.txt"
  output:
    "data/{sample}/perturb_status.txt"
  log:
    "data/{sample}/logs/perturbation_status.log"
  params:
    vector_prefix = config["create_vector_ref"]["vector_prefix"],
    min_txs = lambda wildcards: config["perturbation_status"]["min_txs"][wildcards.sample]
  conda:
    "../envs/r_dropseq_tools.yml"
  shell:
    "Rscript scripts/perturbation_status.R --infile {input} --outfile {output} "
    "--vector_patter {params.vector_prefix} --min_txs {params.min_txs} --trim 2> {log}"

# compile dge report
rule dge_report:
  input:
    dge = "data/{sample}/dge.txt",
    tpt_hist = "data/{sample}/dge_tpt_histogram.txt",
    dge_stats = "data/{sample}/dge_summary.txt",
    perturb_stats = "data/{sample}/perturb_status.txt"
  output:
    "results/dge/{sample}_dge_report.html"
  params:
    vector_prefix = config["create_vector_ref"]["vector_prefix"]
  conda:
    "../envs/r_dropseq_tools.yml"
  script:
    "../scripts/dge_report.Rmd"
