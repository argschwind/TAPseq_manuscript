## rules to extract DGE data from aligned CROP-seq reads

### input, output and shell paths are all relative to the project directory ###

# workflow rules -----------------------------------------------------------------------------------

# count UMI observations per cell barcode and gene tag
rule umi_observations:
  input:
    "data/{sample}/filt_gene_tagged_aligned.bam"
  output:
    "data/{sample}/umi_observations_top_{ncells}_cells.txt"
  log:
    "data/{sample}/logs/umi_observations_{ncells}_cells.log"
  params:
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
    "NUM_CORE_BARCODES={wildcards.ncells} "
    "EDIT_DISTANCE={params.edit_distance} "
    "READ_MQ={params.read_mq} "
    "MIN_BC_READ_THRESHOLD={params.min_umi_reads} "
    "RARE_UMI_FILTER_THRESHOLD={params.rare_umi_filter} "
    "2> {log}"

# extract dge with and filter for chimeric reads
rule extract_dge:
  input:
    "data/{sample}/umi_observations_top_{ncells}_cells.txt"
  output:
    dge = "data/{sample}/dge_top_{ncells}_cells.txt",
    dge_stats = "data/{sample}/dge_top_{ncells}_cells_summary.txt",
    tpt_hist = "data/{sample}/dge_top_{ncells}_cells_tpt_histogram.txt"
  log:
    "data/{sample}/logs/extract_dge_{ncells}_cells.log"
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
    "data/{sample}/dge_top_{ncells}_cells.txt"
  output:
    "data/{sample}/perturb_status_top_{ncells}_cells.txt"
  log:
    "data/{sample}/logs/perturbation_status_{ncells}_cells.log"
  params:
    vector_annot = config["perturbation_status"]["vector_annot"],
    min_txs = lambda wildcards: config["perturbation_status"]["min_txs"][wildcards.sample]
  conda:
    "../envs/r_dropseq_tools.yml"
  shell:
    "Rscript scripts/perturbation_status.R -d {input} -o {output} -v {params.vector_annot} "
    "-t {params.min_txs} 2> {log}"

# compile dge report
rule dge_report:
  input:
    dge = "data/{sample}/dge_top_{ncells}_cells.txt",
    tpt_hist = "data/{sample}/dge_top_{ncells}_cells_tpt_histogram.txt",
    dge_stats = "data/{sample}/dge_top_{ncells}_cells_summary.txt",
    perturb_stats = "data/{sample}/perturb_status_top_{ncells}_cells.txt"
  output:
    "results/dge/{sample}_{ncells}_cells_dge_report.html"
  params:
    vector_annot = config["perturbation_status"]["vector_annot"]
  conda:
    "../envs/r_dropseq_tools.yml"
  script:
    "../scripts/dge_report.Rmd"
