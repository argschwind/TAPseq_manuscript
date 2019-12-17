## rules to create whole-genome and tap-seq alignment references

### input, output and shell paths are all relative to the project directory ###

# some alignment references might not contain any CROP-seq vectors. therefore empty vector reference
# files are generated via the rule 'create_empty_vector_ref'. this rule is only executed if no
# vectors are provided for a given alignment reference and the rule 'create_vector_ref' would fail.
# this behavior is achieved by setting execution priorities for these rules.
ruleorder: create_vector_ref > create_empty_vector_ref

# workflow rules -----------------------------------------------------------------------------------
rule download_genome_annot:
  output:
    fasta = temp("data/alignment_references/{species}_genome_{ref_version}/genome.fasta"),
    gtf = temp("data/alignment_references/{species}_genome_{ref_version}/genome.gtf")
  params:
    fasta_url = lambda wildcards: config["download_genome_annot"][wildcards.species]["fasta_url"],
    gtf_url = lambda wildcards: config["download_genome_annot"][wildcards.species]["gtf_url"]
  log: "data/alignment_references/{species}_genome_{ref_version}/logs/download_genome_annot.log"
  conda: "../envs/dropseq_tools.yml"
  group: "align_ref"
  shell:
    "wget -O {output.gtf}.gz {params.gtf_url} 2> {log} && gzip -d {output.gtf}.gz ; "
    "wget -O {output.fasta}.gz {params.fasta_url} 2>> {log} && gzip -d {output.fasta}.gz"
    
rule create_tapseq_annot:
  input:
    lambda wildcards: config["create_tapseq_ref"]["target_genes"]
      [wildcards.species + "_tapseq_" + wildcards.ref_version]
  output:
    fasta = temp("data/alignment_references/{species}_tapseq_{ref_version}/genome.fasta"),
    gtf = temp("data/alignment_references/{species}_tapseq_{ref_version}/genome.gtf")
  params:
    BSgenome = lambda wildcards: config["create_tapseq_ref"]["BSgenome_id"][wildcards.species]
  log: "data/alignment_references/{species}_tapseq_{ref_version}/logs/create_tapseq_ref.log"
  conda: "../envs/r_alignment_ref.yml"
  group: "align_ref"
  shell:
    "Rscript scripts/processing/create_tapseq_annot.R -t {input} -b {params.BSgenome} "
    "-f {output.fasta} -o {output.gtf} 2> {log}"

rule create_vector_ref:
  input:
    lambda wildcards: config["create_vector_ref"]["vector_fasta"][wildcards.align_ref]
  output:
    fasta = temp("data/alignment_references/{align_ref}/cropseq_vectors.fasta"),
    gtf = temp("data/alignment_references/{align_ref}/cropseq_vectors.gtf")
  params:
    output_bn = "data/alignment_references/{align_ref}/cropseq_vectors",
    vector_prefix = config["create_vector_ref"]["vector_prefix"]
  conda: "../envs/dropseq_tools.yml"
  group: "align_ref"
  shell:
    "python scripts/processing/cropseq_vector_reference.py -i {input} -o {params.output_bn} "
    "--prefix {params.vector_prefix}"

rule create_empty_vector_ref:
  output:
    fasta = temp("data/alignment_references/{align_ref}/cropseq_vectors.fasta"),
    gtf = temp("data/alignment_references/{align_ref}/cropseq_vectors.gtf")
  group: "align_ref"
  shell:
    "touch {output.fasta} {output.gtf}"

rule create_cropseq_ref:
  input:
    genome_fasta = "data/alignment_references/{align_ref}/genome.fasta",
    genome_gtf = "data/alignment_references/{align_ref}/genome.gtf",
    vectors_fasta= "data/alignment_references/{align_ref}/cropseq_vectors.fasta",
    vectors_gtf = "data/alignment_references/{align_ref}/cropseq_vectors.gtf"
  output:
    cropseq_ref_fasta = "data/alignment_references/{align_ref}/cropseq_ref.fasta",
    cropseq_ref_gtf = "data/alignment_references/{align_ref}/cropseq_ref.gtf"
  conda: "../envs/dropseq_tools.yml"
  group: "align_ref"
  shell:
    "cat {input.genome_gtf} > {output.cropseq_ref_gtf} && "
    "cat {input.vectors_gtf} >> {output.cropseq_ref_gtf} ; "
    "cat {input.genome_fasta} > {output.cropseq_ref_fasta} && "
    "cat {input.vectors_fasta} >> {output.cropseq_ref_fasta} ; "

rule create_dict:
  input:
    "data/alignment_references/{align_ref}/cropseq_ref.fasta"
  output:
    "data/alignment_references/{align_ref}/cropseq_ref.dict"
  log: "data/alignment_references/{align_ref}/logs/create_dict.log"
  conda: "../envs/dropseq_tools.yml"
  group: "align_ref"
  shell:
    "picard CreateSequenceDictionary "
    "REFERENCE={input} "
    "OUTPUT={output} "
    "2> {log}"

rule create_refflat:
  input:
    annot = "data/alignment_references/{align_ref}/cropseq_ref.gtf",
    seq_dict = "data/alignment_references/{align_ref}/cropseq_ref.dict"
  output:
    "data/alignment_references/{align_ref}/cropseq_ref.refFlat"
  log: "data/alignment_references/{align_ref}/logs/create_refflat.log"
  conda: "../envs/dropseq_tools.yml"
  group: "align_ref"
  shell:
    "ConvertToRefFlat "
    "ANNOTATIONS_FILE={input.annot} "
    "SEQUENCE_DICTIONARY={input.seq_dict} "
    "OUTPUT={output} "
    "2> {log}"

rule create_genomedir:
  input:
    fasta = "data/alignment_references/{align_ref}/cropseq_ref.fasta",
    gtf = "data/alignment_references/{align_ref}/cropseq_ref.gtf"
  output:
    directory("data/alignment_references/{align_ref}/genomeDir")
  params:
    outprefix = "data/alignment_references/{align_ref}/star.",
    sjdb_overhang = config["create_genomedir"]["sjdb_overhang"]
  log: "data/alignment_references/{align_ref}/logs/create_genomedir.log"
  threads: config["create_genomedir"]["threads"]
  conda: "../envs/dropseq_tools.yml"
  shell:
    "mkdir -p {output} && "  # create genomeDir directory
    "STAR --runThreadN {threads} "
    "--runMode genomeGenerate "
    "--genomeDir {output} "
    "--genomeFastaFiles {input.fasta} "
    "--sjdbGTFfile {input.gtf} "
    "--sjdbOverhang {params.sjdb_overhang} "
    "--outFileNamePrefix {params.outprefix} ; "
    "mv data/alignment_references/{wildcards.align_ref}/star.Log.out {log}"
