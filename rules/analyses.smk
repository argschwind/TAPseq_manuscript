## perform all downstream analyses

### input, output and shell paths are all relative to the project directory ###

# python function(s) to infer more complex input files

# fastq read 2 (genome read)
def fastq_read2(sample):
  fastq_dir = config["samples"][sample]
  f2 = glob.glob(fastq_dir + "/*" + sample  + "_2_sequence.txt.gz")
  return f2

# get cell barcode whitelist file for a given sample
def get_cbc_whitelist(wildcards):
  sample = wildcards.sample
  if sample in config["screen"]:
    whitelist = "meta_data/10x_cell_barcode_whitelists/screen_10x_bc_whitelist_737k_201608.txt.gz"
  elif sample == "WholeTx":
    whitelist = "meta_data/10x_cell_barcode_whitelists/wholeTx_10x_bc_whitelist_737k_201608.txt.gz"
  elif sample == "TAP":
    whitelist = "meta_data/10x_cell_barcode_whitelists/TAP_10x_bc_whitelist_737k_201608.txt.gz"
  else:
    whitelist = "meta_data/10x_cell_barcode_whitelists/10x_bc_whitelist_737k_201608.txt"
  return whitelist

# get correct perturbation status input file for DE strategies
def perturb_status_file(wildcards):
  sample = wildcards.sample
  strategy = wildcards.strategy
  if strategy == "perGRNA":
    pert_file = "data/" + sample + "/perturb_status.txt"
  elif strategy == "perEnh":
    pert_file = "data/" + sample + "/perturb_status_collapsed.txt"
  else:
    raise ValueError("strategy mus be 'perGRNA' or 'perEnh'!")
  return(pert_file)

# validation experiments ---------------------------------------------------------------------------

# rule to downsample genic reads to the same sequencing depth per sample and then extract dge
rule downsample:
  input:
    umi_obs = "data/{sample}/umi_observations.txt",
    whitelist = get_cbc_whitelist
  output:
    dge = "data/{sample}/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
    dge_stats = "data/{sample}/downsampled/dge_{rpc}_avg_reads_per_cell_summary.txt",
    tpt_hist = "data/{sample}/downsampled/dge_{rpc}_avg_reads_per_cell_tpt_histogram.txt"
  params:
    reads = lambda wildcards: config["dge_ncells"][wildcards.sample] * int(wildcards.rpc),
    tpt_threshold = config["extract_dge"]["tpt_threshold"],
    seed = 20190607
  log: "data/{sample}/logs/downsample_{rpc}_avg_reads_per_cell.log"
  conda: "../envs/dropseq_tools.yml"
  shell:
    "python scripts/extract_dge.py -i {input.umi_obs} -o {output.dge} -w {input.whitelist} "
    "--sample {params.reads} --seed {params.seed} --tpt_threshold {params.tpt_threshold} 2> {log}"
  
# calculate reads on target genes enrichment
rule reads_on_target:
  input:
    bam = ["data/" + sample + "/gene_tagged_aligned.bam" for sample in config["validation"]],
    fastq = [fastq_read2(sample) for sample in config["validation"]],
    target_genes = "meta_data/target_genes_validation.csv"
  output:
    "data/reads_on_target.csv"
  params:
    vector_prefix = config["create_vector_ref"]["vector_prefix"]
  conda: "../envs/r_analyses.yml"
  script:
    "../scripts/reads_on_target.R"
  
# TAP-seq vs CROP-seq using downsampled dge to same sequencing depth per sample
rule tapseq_vs_cropseq:
  input:
    reads_on_target = "data/reads_on_target.csv",
    target_genes = "meta_data/target_genes_validation.csv",
    dge = [expand("data/{sample}/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             sample = "Sample10X", rpc = config["downsample"]["reads_per_cell"]["Sample10X"]),
           expand("data/{sample}/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             sample = ["11iv210ng", "11iv22ng", "8iv210ng", "8iv22ng"],
             rpc = config["downsample"]["reads_per_cell"]["tap-seq"])]
  output:
    "results/tapseq_vs_cropseq.html"
  conda: "../envs/r_analyses.yml"
  script:
    "../scripts/tapseq_vs_cropseq.Rmd"
  
# analysis of downsampled dge to same sequencing depth per sample for power comparison between
# TAP-seq and CROP-seq
rule downsampled_dge:
  input:
    dge = expand("data/{sample}/downsampled/dge_3500_avg_reads_per_cell.txt",
      sample = config["validation"]),
    pert_status = expand("data/{sample}/perturb_status.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
  output:
    "results/downsampled_dge.html"
  conda: "../envs/r_analyses.yml"
  script:
    "../scripts/downsampled_dge.Rmd"
  
# downsample TAP-seq and CROP-seq dge to a fixed number of reads on target genes per cell
rule downsampled_target_reads:
  input:
    umi_obs = expand("data/{sample}/umi_observations.txt", sample = config["validation"]),
    pert_status = expand("data/{sample}/perturb_status.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
    whitelist = "meta_data/10x_cell_barcode_whitelists/10x_bc_whitelist_737k_201608.txt"
  output:
    "results/downsampled_target_reads.html"
  conda: "../envs/r_analyses.yml"
  script:
    "../scripts/downsampled_target_reads.Rmd"
  
# screen experiments -------------------------------------------------------------------------------

# QC of screen data (dge data and detected perturbations)
rule screen_data_qc:
  input:
    dge = expand("data/{sample}/dge.txt", sample = config["screen"]),
    perturb_status = expand("data/{sample}/perturb_status.txt", sample = config["screen"]),
    dge_stats = expand("data/{sample}/dge_summary.txt", sample = config["screen"]),
    valid_dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
    exp_data = "meta_data/screen_experimental_info.csv",
    vctr_seqs = expand("meta_data/cropseq_vectors_{chr}_screen.fasta", chr = ["chr8", "chr11"])
  output:
    "results/screen_data_qc.html"
  params:
    vector_pattern = config["create_vector_ref"]["vector_prefix"]
  conda: "../envs/r_analyses.yml"
  script:
    "../scripts/screen_data_qc.Rmd"
    
# collapse perturbation status by gRNA vector targets
rule collapse_perturbations:
  input:
    pert_status = "data/{sample}/perturb_status.txt",
    grna_targets = lambda wildcards: config["collapse_perturbations"]["targets"][wildcards.sample]
  output:
    "data/{sample}/perturb_status_collapsed.txt"
  log: "data/{sample}/logs/collapse_perturbations.log"
  conda: "../envs/dropseq_tools.yml"
  shell:
    "python scripts/collapse_perturbations.py -i {input.pert_status} -t {input.grna_targets} "
    "-o {output} 2> {log}"
  
# map enhancer-gene pairs using MAST for differential gene expression testing. strategy can be
# either "perEnh" or "perGRNA", specifying whether DE tests should be performed using perturbations
# collapsed per enhancer or per gRNA. covars specifies the cell-level covariates that should be used
# for differential expression tests. Can be one of: 'noCovar', 'nGenesCovar', or 'XpcCovar', where
# X is the number of principal components (e.g. 2 = PCs 1+2).
rule diff_expr:
  input:
    dge = "data/{sample}/dge.txt",
    perturb_status = perturb_status_file
  output:
    results = "data/{sample}/diff_expr/output_{strategy}_{covars}.csv",
    ncells =  "data/{sample}/diff_expr/ncells_{strategy}_{covars}.csv"
  params:
    vector_pattern = config["create_vector_ref"]["vector_prefix"],
    exclude_lanes = lambda wildcards: config["map_enhancers"]["remove_lanes"][wildcards.sample],
    min_cells = lambda wildcards: config["map_enhancers"]["min_cells"][wildcards.strategy],
    seed = 20190324
  log: "data/{sample}/logs/diff_expr_{strategy}_{covars}.log"
  threads: config["map_enhancers"]["threads"]
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/differential_expression.R"
  
# compare different covariates for differential expression tests
rule compare_covariates:
  input:
    results = expand("data/{sample}/diff_expr/output_{strategy}_{covars}.csv",
      sample = config["screen"], strategy = ["perEnh", "perGRNA"],
      covars = ["noCovar", "1pcCovar", "2pcCovar", "nGenesCovar"]),
    dge = expand("data/{sample}/dge.txt", sample = config["screen"]),
    perturb_status = expand("data/{sample}/perturb_status_collapsed.txt",
      sample = config["screen"]),
    annot = ["meta_data/target_genes_chr8_screen.gtf", "meta_data/target_genes_chr11_screen.gtf"]
  output:
    "results/compare_covariates.html"
  params:
    vector_pattern = config["create_vector_ref"]["vector_prefix"]
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/compare_covariates.Rmd"
  
# process de results and calculate useful stats such as confidence levels and distance to TSS.
# generates input files for detailed enhancer analyses
rule process_de_results:
  input:
    results = expand("data/{sample}/diff_expr/output_{strategy}_{{covars}}.csv",
      sample = config["screen"], strategy = ["perEnh", "perGRNA"]),
    annot = ["meta_data/target_genes_chr8_screen.gtf", "meta_data/target_genes_chr11_screen.gtf"],
    enh_coords = ["meta_data/enhancers_chr8_screen.bed", "meta_data/enhancers_chr11_screen.bed"]
  output:
    "data/diff_expr_screen_{covars}.csv"
  params:
    confidence_fdr = 0.05  # fdr threshold for single gRNA hits used to assign confidence levels
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/process_de_results.R"
  
# perform basic analyses of differential expression results and enhancer - target gene pairs.
rule map_enhancers:
  input:
    de_output = expand("data/{sample}/diff_expr/output_{strategy}_nGenesCovar.csv",
      sample = config["screen"], strategy = ["perEnh", "perGRNA"]),
    ncells = expand("data/{sample}/diff_expr/ncells_{strategy}_nGenesCovar.csv",
      sample = config["screen"], strategy = ["perEnh", "perGRNA"]),
    processed_results = "data/diff_expr_screen_nGenesCovar.csv",
    perturb_status = expand("data/{sample}/perturb_status_collapsed.txt",
      sample = config["screen"]),
    vector_targets = expand("meta_data/vector_targets_{chr}_screen.csv", chr = ["chr8", "chr11"])
  output:
    "results/map_enhancers.html"
  params:
    vector_pattern = config["create_vector_ref"]["vector_prefix"]
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/map_enhancers.Rmd"
  
# download chromatin data
rule download_chromatin_data:
  output:
    "data/k562_chromatin_data/{assay}_encode_chipseq_fcOverCtrl.bigWig"
  params:
    url = lambda wildcards: config["chromatin_analyses"]["encode_chip"][wildcards.assay]
  conda: "../envs/dropseq_tools.yml"
  shell:
    "wget -O {output} {params.url}"
    
# download Hi-C data
rule download_hic_data:
  output:
    expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_{file}",
      chr = ["chr8", "chr11"], file = ["5kb.RAWobserved", "5kb.KRnorm"])
  params:
    url = config["chromatin_analyses"]["rao_hic"]
  conda: "../envs/dropseq_tools.yml"
  shell:
    "wget -c {params.url} -O - | tar -xz "
    "K562/5kb_resolution_intrachromosomal/chr8/MAPQG0/ "
    "K562/5kb_resolution_intrachromosomal/chr11/MAPQG0/ ; "
    "mv K562/5kb_resolution_intrachromosomal/chr8/MAPQG0/* "
    "data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/chr8/MAPQG0 ; "
    "mv K562/5kb_resolution_intrachromosomal/chr11/MAPQG0/* "
    "data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/chr11/MAPQG0 ; "
    "rm -r K562"
    
# perform in depth chromatin analyses of enhancer mapping results
rule chromatin_analyses:
  input:
    processed_results = "data/diff_expr_screen_nGenesCovar.csv",
    encode_chip = expand("data/k562_chromatin_data/{assay}_encode_chipseq_fcOverCtrl.bigWig",
      assay = config["chromatin_analyses"]["encode_chip"])
  output:
    "results/chromatin_analyses_screen.html"
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/chromatin_analyses_screen.Rmd"
    
# perform hi-c analysis
rule hic_analysis:
  input:
    processed_results = "data/diff_expr_screen_nGenesCovar.csv",
    hic_raw = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.RAWobserved",
      chr = ["chr8", "chr11"]),
    hic_norm = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.KRnorm",
      chr = ["chr8", "chr11"])
  output:
    "results/hic_analysis.html"
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/hic_analysis.Rmd"
