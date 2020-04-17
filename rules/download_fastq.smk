## rules to download fastq files from SRA

### input, output and shell paths are all relative to the project directory ###

# rule orders for different processing rules: some samples require that multiple fastq files are
# merged into single files. rules to merge files use the respective samples as wildcard constraints
# and have priority over the simple renaming rule
ruleorder: merge_fastq > rename_fastq
ruleorder: merge_fastq_i7 > rename_fastq
ruleorder: merge_fastq_i7_trim > rename_fastq

# python function(s) to infer complex parameters  --------------------------------------------------

# get all temporary files for a sample used for download from SRA (2 filea per run)
def get_run_files(wildcards):
  sample=wildcards.sample
  files = expand("raw_data/{sample}/{run}_{read}.fastq.gz", sample = sample,
    run = config["samples"][sample], read = [1, 2])
  return files

# workflow rules -----------------------------------------------------------------------------------

# helper rule to download fastq files for one paired-end SRA run
rule download_sra_run:
  output:
    temp(expand("raw_data/{{sample}}/{{run}}_{read}.fastq.gz", read = [1, 2]))
  params:
    outdir = "raw_data/{sample}"
  log: "raw_data/{sample}/logs/download_{run}.log"
  conda: "../envs/fastq_download.yml"
  shell:
    "cmd='fastq-dump --split-files --gzip --outdir {params.outdir} {wildcards.run}' ; "
    "echo $cmd > {log} ; $cmd >> {log}"

# simple renaming of downloaded fastq files for a sample with only one associated run
rule rename_fastq:
  input:
    get_run_files
  output:
    fwd = "raw_data/{sample}/{sample}_1_sequence.txt.gz",
    rev = "raw_data/{sample}/{sample}_2_sequence.txt.gz"
  params:
    access = lambda wildcards: config["samples"][wildcards.sample],
    outdir = "raw_data/{sample}"
  log: "raw_data/{sample}/logs/rename_fastq.log"
  conda: "../envs/fastq_download.yml"
  shell:
    "mv -v {params.outdir}/{params.access}_1.fastq.gz {output.fwd} > {log} & "
    "mv -v {params.outdir}/{params.access}_2.fastq.gz {output.rev} >> {log}"
    
# simple merge of multiple fastq files for a sample into one 'fwd' and one 'rev' read file
rule merge_fastq:
  input:
    get_run_files
  output:
    fwd = "raw_data/{sample}/{sample}_1_sequence.txt.gz",
    rev = "raw_data/{sample}/{sample}_2_sequence.txt.gz"
  params:
    access = lambda wildcards: config["samples"][wildcards.sample],
    outdir = "raw_data/{sample}"
  log: "raw_data/{sample}/logs/merge_fastq.log"
  wildcard_constraints:
    sample = "WholeTotalBM|WholeKitBM"
  conda: "../envs/fastq_download.yml"
  shell:
    "bash scripts/download_fastq/merge_fastq.sh {params.outdir} {output.fwd} {output.rev} "
    "> {log}"

# merge multiple fastq files for a sample with a lane specific i7 barcode added to the cell barcodes
# to avoid cell barcode collisions
rule merge_fastq_i7:
  input:
    get_run_files
  output:
    fwd = "raw_data/{sample}/{sample}_1_sequence.txt.gz",
    rev = "raw_data/{sample}/{sample}_2_sequence.txt.gz"
  params:
    access = lambda wildcards: config["samples"][wildcards.sample],
    outdir = "raw_data/{sample}"
  log: "raw_data/{sample}/logs/merge_fastq.log"
  wildcard_constraints:
    sample = "11iScreen1|WholeTx|TAPtotalBM|TAPkitBM"
  conda: "../envs/fastq_download.yml"
  shell:
    "bash scripts/download_fastq/merge_fastq_i7.sh {params.outdir} {output.fwd} {output.rev} "
    "> {log}"

# merge multiple fastq files with lane specific i7 barcode added to cell barcodes and trim all
# genome (read 2) reads to 58 bp
rule merge_fastq_i7_trim:
  input:
    get_run_files
  output:
    fwd = "raw_data/{sample}/{sample}_1_sequence.txt.gz",
    rev = "raw_data/{sample}/{sample}_2_sequence.txt.gz"
  params:
    access = lambda wildcards: config["samples"][wildcards.sample],
    outdir = "raw_data/{sample}"
  log: "raw_data/{sample}/logs/merge_fastq.log"
  wildcard_constraints:
    sample = "8iScreen1"
  conda: "../envs/fastq_download.yml"
  shell:
    "bash scripts/download_fastq/merge_fastq_i7.sh {params.outdir} {output.fwd} {output.rev} "
    "> {log} ; trimmomatic SE {output.rev} {params.outdir}/tmp.txt.gz CROP:58 2>> {log} ; "
    "mv -v {params.outdir}/tmp.txt.gz {output.rev} >> {log}"
