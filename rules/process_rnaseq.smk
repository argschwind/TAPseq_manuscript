## rules to process bulk RNA-seq data

### input, output and shell paths are all relative to the project directory ###

# download kallisto transcript index
rule download_kallisto_index:
  output:
    index = "data/bulk_rnaseq/kallisto_txs_indices/{species}/transcriptome.idx",
    txs_to_genes = "data/bulk_rnaseq/kallisto_txs_indices/{species}/transcripts_to_genes.txt"
  params:
    url = lambda wildcards: config["download_kallisto_index"][wildcards.species]
  conda: "../envs/kallisto.yml"
  shell:
    "wget -nv -c {params.url} -O - | tar -xvz; "
    "mv {wildcards.species}/* data/bulk_rnaseq/kallisto_txs_indices/{wildcards.species}; "
    "rmdir {wildcards.species}"

# quantify transcripts using kallisto
rule kallisto_quant:
  input:
    fastq = lambda wildcards: config["rnaseq_samples"][wildcards.sample],
    index = lambda wildcards: config["kallisto_index"][wildcards.sample]
  output:
    abundance_h5 =  "data/bulk_rnaseq/{sample}/abundance.h5",
    abundance_tsv = "data/bulk_rnaseq/{sample}/abundance.tsv",
    run_info = "data/bulk_rnaseq/{sample}/run_info.json"
  log: "data/bulk_rnaseq/{sample}/logs/kallisto_quant.log"
  params:
    output_dir = "data/bulk_rnaseq/{sample}",
    bootstrap = config["kallisto_quant"]["bootstrap"],
    frag_len = lambda wildcards: config["kallisto_quant"]["frag_len"][wildcards.sample],
    frag_sd  = lambda wildcards: config["kallisto_quant"]["frag_sd"][wildcards.sample]
  threads: config["kallisto_quant"]["threads"]
  conda: "../envs/kallisto.yml"
  shell:
    "kallisto quant -t {threads} -i {input.index} -o {params.output_dir} -b {params.bootstrap} "
    "--single -l {params.frag_len} -s {params.frag_sd} {input.fastq} 2> {log}"
