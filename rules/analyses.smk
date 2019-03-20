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
    whitelist = "meta_data/screen_10x_bc_whitelist_737k_201608.txt.gz"
  else:
    whitelist = "meta_data/10x_bc_whitelist_737k_201608.txt"
  return whitelist
  
# validation experiments ---------------------------------------------------------------------------

# rule to downsample genic reads to the same sequencing depth per sample and then extract dge
rule downsample:
  input:
    umi_obs = "data/{sample}/umi_observations.txt",
    whitelist = get_cbc_whitelist
  output:
    dge = "data/{sample}/downsampled_dge.txt",
    dge_stats = "data/{sample}/downsampled_dge_summary.txt",
    tpt_hist = "data/{sample}/downsampled_dge_tpt_histogram.txt"
  log:
    "data/{sample}/logs/downsample.log"
  params:
    reads = lambda wildcards: config["dge_ncells"][wildcards.sample]
              * config["downsample"]["reads_per_cell"],
    tpt_threshold = config["extract_dge"]["tpt_threshold"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "python scripts/extract_dge.py -i {input.umi_obs} -o {output.dge} -w {input.whitelist} "
    "--sample {params.reads} --seed 20190204 --tpt_threshold {params.tpt_threshold} 2> {log}"

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
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/reads_on_target.R"

# TAP-seq vs CROP-seq using downsampled dge to same sequencing depth per sample
rule tapseq_vs_cropseq:
  input:
    reads_on_target = "data/reads_on_target.csv",
    target_genes = "meta_data/target_genes_validation.csv",
    dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"])
  output:
    "results/tapseq_vs_cropseq.html"
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/tapseq_vs_cropseq.Rmd"

# analysis of downsampled dge to same sequencing depth per sample for power comparison between
# TAP-seq and CROP-seq
rule downsampled_dge:
  input:
    dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"]),
    pert_status = expand("data/{sample}/perturb_status.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv"
  output:
    "results/downsampled_dge.html"
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/downsampled_dge.Rmd"
    
# downsample TAP-seq and CROP-seq dge to a fixed number of reads on target genes per cell
rule downsampled_target_reads:
  input:
    umi_obs = expand("data/{sample}/umi_observations.txt", sample = config["validation"]),
    pert_status = expand("data/{sample}/perturb_status.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
    whitelist = "meta_data/10x_bc_whitelist_737k_201608.txt"
  output:
    "results/downsampled_target_reads.html"
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/downsampled_target_reads.Rmd"

# screen experiments -------------------------------------------------------------------------------

# chromosome 8 screen QC
rule qc_8iScreen1:
  input:
    dge = "data/8iScreen1/dge.txt",
    perturb_status = "data/8iScreen1/perturb_status.txt",
    valid_dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
    exp_data = "meta_data/screen_experimental_info.csv",
    vctr_seqs = "meta_data/cropseq_vectors_chr8_screen.fasta"
  output:
    "results/8iScreen1_qc.html"
  params:
    vector_pattern = config["create_vector_ref"]["vector_prefix"]
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/8iScreen1_qc.Rmd"
    
# chromosome 11 screen QC
rule qc_11iScreen1:
  input:
    dge = "data/11iScreen1/dge.txt",
    perturb_status = "data/11iScreen1/perturb_status.txt",
    valid_dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
    exp_data = "meta_data/screen_experimental_info.csv",
    vctr_seqs = "meta_data/cropseq_vectors_chr11_screen.fasta"
  output:
    "results/11iScreen1_qc.html"
  params:
    vector_pattern = config["create_vector_ref"]["vector_prefix"]
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/11iScreen1_qc.Rmd"
