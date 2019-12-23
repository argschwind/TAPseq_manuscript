## perform nature methods revision analyses

### input, output and shell paths are all relative to the project directory ###

# downsample new data for revisions    
rule downsample_revision_data:
  input:
    dge = [expand("data/wtxmmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["wtxmmix"]),
           expand("data/wtxlung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["wtxlung"]),
           expand("data/tapmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["tapmix"]),
           expand("data/taplung/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
             rpc = config["downsample"]["reads_per_cell"]["taplung"]),
#           expand("data/19s005246/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
#             rpc = config["downsample"]["reads_per_cell"]["19s005246"]),
#           expand("data/taphumanmix/downsampled/dge_{rpc}_avg_reads_per_cell.txt",
#             rpc = config["downsample"]["reads_per_cell"]["taphumanmix"]),
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
             rpc = config["downsample"]["reads_per_cell"]["perturb10xWTX"])]

# filter mouse cell type mix data  (umi obs) for the 150 genes that are expressed in all cells
rule filter_mmix:
  input:
    umi_obs = "data/{sample}/umi_observations.txt",
    target_genes = "meta_data/target_gene_panels/target_genes_mouse_mix.csv"
  output:
    "data/{sample}/umi_obs_150_common_genes.txt"
  conda: "../envs/r_analyses.yml"
  script:
    "../scripts/analyses/filter_mmix_umi_obs.R"

# rule to downsample mouse cell type mix genic reads on the 150 common genes to the same sequencing
# depth per sample and then extract dge data
rule downsample_150_genes:
  input:
    umi_obs = "data/{sample}/umi_obs_150_common_genes.txt",
    whitelist = get_cbc_whitelist
  output:
    dge = "data/{sample}/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
    dge_stats = "data/{sample}/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell_summary.txt",
    tpt_hist = "data/{sample}/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell_tpt_histogram.txt"
  params:
    reads = lambda wildcards: config["cell_numbers"][wildcards.sample] * int(wildcards.rpc),
    tpt_threshold = config["extract_dge"]["tpt_threshold"],
    whitelist_arg = lambda wildcards, input: get_whitelist_arg(input.whitelist),
    seed = 20191220
  log: "data/{sample}/logs/downsample_150_genes_{rpc}_avg_reads_per_cell.log"
  conda: "../envs/dropseq_tools.yml"
  shell:
    "python scripts/processing/extract_dge.py -i {input.umi_obs} -o {output.dge} "
    "{params.whitelist_arg} --sample {params.reads} --seed {params.seed} --tpt_threshold "
    "{params.tpt_threshold} 2> {log}"
    
# downsample mmix 150 common genes dge data
rule downsample_mmix_150_genes:
  input:
    expand("data/taplung/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
      rpc = [1000, 1500, 3500, 5000, 10000, 15000, 20000, 40000]),
    expand("data/tapmix/downsampled_150_genes/dge_{rpc}_avg_reads_per_cell.txt",
      rpc = [1000, 1500, 3500, 5000, 10000, 15000, 20000, 40000, 80000, 150000])
