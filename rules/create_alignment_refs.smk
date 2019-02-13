## rules to create whole-genome and tap-seq alignment references

### input, output and shell paths are all relative to the project directory ###

# workflow rules -----------------------------------------------------------------------------------
rule download_genome_annot:
  output:
    fasta = temp("data/genome_reference/genome.fasta"),
    gtf = temp("data/genome_reference/genome.gtf")
  params:
    fasta_url = config["download_genome_annot"]["fasta_url"],
    gtf_url = config["download_genome_annot"]["gtf_url"]
  log:
    "data/genome_reference/logs/download_genome_annot.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "wget -O {output.gtf}.gz {params.gtf_url} 2> {log} && gzip -d {output.gtf}.gz ; "
    "wget -O {output.fasta}.gz {params.fasta_url} 2>> {log} && gzip -d {output.fasta}.gz"
    
rule create_tapseq_annot:
  input:
    lambda wildcards: config["create_tapseq_ref"]["target_genes"][wildcards.align_ref]
  output:
    fasta = temp("data/{align_ref}/genome.fasta"),
    gtf = temp("data/{align_ref}/genome.gtf")
  log:
    "data/{align_ref}/logs/create_tapseq_ref.log"
  params:
    BSgenome = config["create_tapseq_ref"]["BSgenome_id"]
  conda:
    "../envs/r_alignment_ref.yml"
  shell:
    "Rscript scripts/create_tapseq_annot.R -t {input} -b {params.BSgenome} -f {output.fasta} "
    "-o {output.gtf} 2> {log}"

rule create_vector_ref:
  input:
    lambda wildcards: config["create_vector_ref"]["vector_fasta"][wildcards.align_ref]
  output:
    fasta = temp("data/{align_ref}/cropseq_vectors.fasta"),
    gtf = temp("data/{align_ref}/cropseq_vectors.gtf")
  params:
    output_bn = "data/{align_ref}/cropseq_vectors",
    vector_prefix = config["create_vector_ref"]["vector_prefix"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "python scripts/cropseq_vector_reference.py -i {input} -o {params.output_bn} "
    "--prefix {params.vector_prefix}"

rule create_cropseq_ref:
  input:
    genome_fasta = "data/{align_ref}/genome.fasta",
    genome_gtf = "data/{align_ref}/genome.gtf",
    vectors_fasta= "data/{align_ref}/cropseq_vectors.fasta",
    vectors_gtf = "data/{align_ref}/cropseq_vectors.gtf"
  output:
    cropseq_ref_fasta = "data/{align_ref}/cropseq_ref.fasta",
    cropseq_ref_gtf = "data/{align_ref}/cropseq_ref.gtf"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "cat {input.genome_gtf} > {output.cropseq_ref_gtf} && "
    "cat {input.vectors_gtf} >> {output.cropseq_ref_gtf} ; "
    "cat {input.genome_fasta} > {output.cropseq_ref_fasta} && "
    "cat {input.vectors_fasta} >> {output.cropseq_ref_fasta} ; "

rule create_dict:
  input:
    "data/{align_ref}/cropseq_ref.fasta"
  output:
    "data/{align_ref}/cropseq_ref.dict"
  log:
    "data/{align_ref}/logs/create_dict.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard CreateSequenceDictionary "
    "REFERENCE={input} "
    "OUTPUT={output} "
    "2> {log}"

rule create_refflat:
  input:
    annot = "data/{align_ref}/cropseq_ref.gtf",
    seq_dict = "data/{align_ref}/cropseq_ref.dict"
  output:
    "data/{align_ref}/cropseq_ref.refFlat"
  log:
    "data/{align_ref}/logs/create_refflat.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "ConvertToRefFlat "
    "ANNOTATIONS_FILE={input.annot} "
    "SEQUENCE_DICTIONARY={input.seq_dict} "
    "OUTPUT={output} "
    "2> {log}"

rule create_genomedir:
  input:
    fasta = "data/{align_ref}/cropseq_ref.fasta",
    gtf = "data/{align_ref}/cropseq_ref.gtf"
  output:
    directory("data/{align_ref}/genomeDir")
  log:
    "data/{align_ref}/logs/create_genomedir.log"
  params:
    outprefix = "data/{align_ref}/star.",
    sjdb_overhang = config["create_genomedir"]["sjdb_overhang"]
  threads: config["create_genomedir"]["star_threads"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "mkdir -p {output} && "  # create genomeDir directory
    "STAR --runThreadN {threads} "
    "--runMode genomeGenerate "
    "--genomeDir {output} "
    "--genomeFastaFiles {input.fasta} "
    "--sjdbGTFfile {input.gtf} "
    "--sjdbOverhang {params.sjdb_overhang} "
    "--outFileNamePrefix {params.outprefix} ; "
    "mv data/{wildcards.align_ref}/star.Log.out {log}"
