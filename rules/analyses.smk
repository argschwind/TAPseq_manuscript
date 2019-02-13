## perform all downstream analyses

### input, output and shell paths are all relative to the project directory ###

# python functions to infer more complex input files -----------------------------------------------

# fastq read 2 (genome read)
def fastq_read2(sample):
  fastq_dir = config["samples"][sample]
  f2 = glob.glob(fastq_dir + "/*" + sample  + "_2_sequence.txt.gz")
  return f2
  
# rules --------------------------------------------------------------------------------------------

# rule to downsample data sets to a specified number of genic reads and extract dge
rule downsample:
  input:
    "data/{sample}/umi_observations.txt"
  output:
    dge = "data/{sample}/downsampled_dge.txt",
    dge_stats = "data/{sample}/downsampled_dge_summary.txt",
    tpt_hist = "data/{sample}/downsampled_dge_tpt_histogram.txt"
  log:
    "data/{sample}/logs/downsample.txt"
  params:
    reads = lambda wildcards: config["dge_ncells"][wildcards.sample]
              * config["downsample"]["reads_per_cell"],
    tpt_threshold = config["extract_dge"]["tpt_threshold"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "python scripts/extract_dge.py -i {input} -o {output.dge} --sample {params.reads} "
    "--seed 20190204 --tpt_threshold {params.tpt_threshold} 2> {log}"

# calculate reads on target genes enrichment
rule reads_on_target:
  input:
    bam = ["data/" + sample + "/filt_gene_tagged_aligned.bam" for sample in config["validation"]],
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

# tap-seq vs crop-seq
rule tapseq_vs_dropseq:
  input:
    reads_on_target = "data/reads_on_target.csv",
    target_genes = "meta_data/target_genes_validation.csv",
    dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"]),
    whitelist = "meta_data/10x_bc_whitelist_737k_201608.txt"
  output:
    "results/tapseq_vs_cropseq.html"
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/tapseq_vs_cropseq.Rmd"

# downsample dge for power comparison between tap-seq and crop-seq
rule downsampled_dge:
  input:
    dge = expand("data/{sample}/downsampled_dge.txt", sample = config["validation"]),
    pert_status = expand("data/{sample}/perturb_status.txt", sample = config["validation"]),
    target_genes = "meta_data/target_genes_validation.csv",
    whitelist = "meta_data/10x_bc_whitelist_737k_201608.txt"
  output:
    "results/downsampled_dge.html"
  conda:
    "../envs/r_analyses.yml"
  script:
    "../scripts/downsampled_dge.Rmd"
