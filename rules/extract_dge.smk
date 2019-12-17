## rules to extract DGE data from aligned CROP-seq reads

### input, output and shell paths are all relative to the project directory ###

# some samples do not contain any perturbations and therefore the dge report is different and the
# rule "perturbation_status" does not need to be executed. the wildcards for the rule 
# "perturbation_status" are resticted to samples with perturbations (see rule) so that it can only
# be executed for these samples. snakemake's ruleorder feature is used to compile the regular
# dge report unless "perturbation_status" can't be executed and therefore a dge report without
# perturbations is created.
ruleorder: dge_report > dge_report_no_perts

# python function(s) to infer more complex input files ---------------------------------------------

# get cell barcode whitelist file for a given sample
def get_cbc_whitelist(wildcards):
  sample = wildcards.sample
  if sample in config["screen"]:
    whitelist = "meta_data/10x_cell_barcode_whitelists/screen_10x_bc_whitelist_737k_201608.txt.gz"
  elif sample == "WholeTx":
    whitelist = "meta_data/10x_cell_barcode_whitelists/wholeTx_10x_bc_whitelist_737k_201608.txt.gz"
  elif sample == "TAP":
    whitelist = "meta_data/10x_cell_barcode_whitelists/TAP_10x_bc_whitelist_737k_201608.txt.gz"
  elif sample in config["drop-seq"]:
    whitelist = []  # no cell barcode whitelist available for drop-seq samples
  else:
    whitelist = "meta_data/10x_cell_barcode_whitelists/10x_bc_whitelist_737k_201608.txt"
  return whitelist

# create whitelist argument for extract_dge.py. required because whitelist can be an empty list
# (see above), which needs to be translated to an empty ('') string argument
def get_whitelist_arg(whitelist):
  if isinstance(whitelist, str):
    whitelist_arg = "-w " + whitelist + " "
  else:
    whitelist_arg = ""
  return whitelist_arg

# workflow rules -----------------------------------------------------------------------------------

# count UMI observations per cell barcode and gene tag
rule umi_observations:
  input:
    "data/{sample}/gene_tagged_aligned.bam"
  output:
    "data/{sample}/umi_observations.txt"
  log:
    "data/{sample}/logs/umi_observations.log"
  params:
    ncells = lambda wildcards: config["cell_numbers"][wildcards.sample],
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
    umi_obs = "data/{sample}/umi_observations.txt",
    whitelist = get_cbc_whitelist
  output:
    dge = "data/{sample}/dge.txt",
    dge_stats = "data/{sample}/dge_summary.txt",
    tpt_hist = "data/{sample}/dge_tpt_histogram.txt"
  log:
    "data/{sample}/logs/extract_dge.log"
  params:
    tpt_threshold = config["extract_dge"]["tpt_threshold"],
    whitelist_arg = lambda wildcards, input: get_whitelist_arg(input.whitelist)
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "python scripts/processing/extract_dge.py -i {input.umi_obs} -o {output.dge} "
    "{params.whitelist_arg} --tpt_threshold {params.tpt_threshold} 2> {log}"
    
# infer perturbation status of each cell
rule perturbation_status:
  input:
    "data/{sample}/dge.txt"
  output:
    "data/{sample}/perturb_status.txt"
  log:
    "data/{sample}/logs/perturbation_status.log"
  wildcard_constraints:
     sample = "(" + "|".join(list(config["perturbation_status"]["min_txs"].keys())) + ")"
  params:
    vector_prefix = config["create_vector_ref"]["vector_prefix"],
    min_txs = lambda wildcards: config["perturbation_status"]["min_txs"][wildcards.sample]
  conda:
    "../envs/r_dropseq_tools.yml"
  shell:
    "Rscript scripts/processing/perturbation_status.R --infile {input} --outfile {output} "
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
    vector_prefix = config["create_vector_ref"]["vector_prefix"],
    min_txs = lambda wildcards: config["perturbation_status"]["min_txs"][wildcards.sample]
  conda:
    "../envs/r_dropseq_tools.yml"
  script:
    "../scripts/processing/dge_report.Rmd"
    
# compile dge report for samples without perturbations
rule dge_report_no_perts:
  input:
    dge = "data/{sample}/dge.txt",
    tpt_hist = "data/{sample}/dge_tpt_histogram.txt",
    dge_stats = "data/{sample}/dge_summary.txt"
  output:
    "results/dge/{sample}_dge_report.html"
  conda:
    "../envs/r_dropseq_tools.yml"
  script:
    "../scripts/processing/dge_report_no_perts.Rmd"
