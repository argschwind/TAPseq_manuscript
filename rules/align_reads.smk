## rules to process and align CROP-seq reads

### input, output and shell paths are all relative to the project directory ###

# python function(s) to infer more complex input files ---------------------------------------------

import glob

# infer input fastq files from dirs in config file and sample wildcard
def get_fastq_files(wildcards):
  indir = config["samples"][wildcards.sample]
  f1 = glob.glob(indir + "/*" + wildcards.sample  + "_1_sequence.txt.gz")
  f2 = glob.glob(indir + "/*" + wildcards.sample  + "_2_sequence.txt.gz")
  return {"fastq1" : f1, "fastq2" : f2}

# workflow rules -----------------------------------------------------------------------------------

# convert fastq input files into one unmapped bam file
rule fastq_to_bam:
  input:
    unpack(get_fastq_files)
  output:
    temp("data/{sample}/unmapped.bam")
  log:
    "data/{sample}/logs/fastq_to_bam.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard FastqToSam "
    "FASTQ={input.fastq1} "
    "FASTQ2={input.fastq2} "
    "OUTPUT={output} "
    "SAMPLE_NAME={wildcards.sample} "
    "2> {log}"

# make sure that input is sorted according to query name
rule sort_unmapped:
  input:
    "data/{sample}/unmapped.bam"
  output:
    temp("data/{sample}/sorted_unmapped.bam")
  log:
    "data/{sample}/logs/sort_unmapped.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard SortSam "
    "INPUT={input} "
    "OUTPUT={output} "
    "SORT_ORDER=queryname "
    "2> {log}"

# tag genome reads with CELL barcodes
rule tag_cell_barcodes:
  input:
    "data/{sample}/sorted_unmapped.bam"
  output:
    bam = temp("data/{sample}/cell_tagged_unmapped.bam"),
    summary = "data/{sample}/cell_tags_summary.txt"
  log:
    "data/{sample}/logs/tag_cell_barcodes.log"
  params:
    base_range = lambda wildcards: config["bc_structure"][wildcards.sample][0],
    base_qual = config["tag_cell_barcodes"]["base_quality"],
    bases_below_qual = config["tag_cell_barcodes"]["num_bases_below_quality"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "TagBamWithReadSequenceExtended "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "SUMMARY={output.summary} "
    "BASE_QUALITY={params.base_qual} "
    "NUM_BASES_BELOW_QUALITY={params.bases_below_qual} "
    "BASE_RANGE={params.base_range} "
    "BARCODED_READ=1 "
    "DISCARD_READ=false "
    "TAG_NAME=XC "
    "2> {log}"

# tag genome reads with MOLECULE barcodes
rule tag_molecule_barcodes:
  input:
    "data/{sample}/cell_tagged_unmapped.bam"
  output:
    bam = temp("data/{sample}/mol_tagged_unmapped.bam"),
    summary = "data/{sample}/mol_tags_summary.txt"
  log:
    "data/{sample}/logs/tag_molecule_barcodes.log"
  params:
    base_range = lambda wildcards: config["bc_structure"][wildcards.sample][1],
    base_qual = config["tag_cell_barcodes"]["base_quality"],
    bases_below_qual = config["tag_cell_barcodes"]["num_bases_below_quality"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "TagBamWithReadSequenceExtended "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "SUMMARY={output.summary} "
    "BASE_QUALITY={params.base_qual} "
    "NUM_BASES_BELOW_QUALITY={params.bases_below_qual} "
    "BASE_RANGE={params.base_range} "
    "BARCODED_READ=1 "
    "DISCARD_READ=true "
    "TAG_NAME=XM "
    "2> {log}"

# filter reads marked as 'rejected' by TagBamWithReadSequenceExtended
rule filter_bam:
  input:
    "data/{sample}/mol_tagged_unmapped.bam"
  output:
    temp("data/{sample}/filt_unmapped.bam")
  log:
    "data/{sample}/logs/filter_bam.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "FilterBAM "
    "INPUT={input} "
    "OUTPUT={output} "
    "TAG_REJECT=XQ "
    "2> {log}"

# trim SMART adapter sequences from 5'
rule trim_starting_sequence:
  input:
    "data/{sample}/filt_unmapped.bam"
  output:
    bam = temp("data/{sample}/adapter_trimmed_unmapped.bam"),
    summary = "data/{sample}/adapter_trimming_report.txt"
  log:
    "data/{sample}/logs/trim_starting_sequence.log"
  params:
    adapter_sequence = config["trim_starting_sequence"]["adapter_sequence"],
    mismatches = config["trim_starting_sequence"]["mismatches"],
    num_bases = config["trim_starting_sequence"]["num_bases"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "TrimStartingSequence "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "OUTPUT_SUMMARY={output.summary} "
    "SEQUENCE={params.adapter_sequence} "
    "MISMATCHES={params.mismatches} "
    "NUM_BASES={params.num_bases} "
    "2> {log}"

# trim polyA sequences from 3'
rule trim_polyA:
  input:
    "data/{sample}/adapter_trimmed_unmapped.bam"
  output:
    bam = temp("data/{sample}/polyA_trimmed_unmapped.bam"),
    summary = "data/{sample}/polyA_trimming_report.txt"
  log:
    "data/{sample}/logs/trim_polyA.log"
  params:
    mismatches = config["trim_polyA"]["mismatches"],
    num_bases = config["trim_polyA"]["num_bases"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "PolyATrimmer "
    "INPUT={input} "
    "OUTPUT={output.bam} "
    "OUTPUT_SUMMARY={output.summary} "
    "MISMATCHES={params.mismatches} "
    "NUM_BASES={params.num_bases} "
    "2> {log}"

# convert to fastq for STAR read aligner
rule sam_to_fastq:
  input:
    "data/{sample}/polyA_trimmed_unmapped.bam"
  output:
    temp("data/{sample}/polyA_trimmed_unmapped.fastq")
  log:
    "data/{sample}/logs/sam_to_fastq.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard SamToFastq "
    "INPUT={input} "
    "FASTQ={output} "
    "2> {log}"

# align reads using STAR
rule star_align:
  input:
    fastq = "data/{sample}/polyA_trimmed_unmapped.fastq",
    genomedir = lambda wildcards: config["align_ref"][wildcards.sample] + "/genomeDir"
  output:
    temp("data/{sample}/star.Aligned.out.bam"),
    "data/{sample}/star.Log.final.out"
  params:
    outprefix = "data/{sample}/star."
  threads: config["star_align"]["threads"]
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "STAR --runThreadN {threads} "
    "--genomeDir {input.genomedir} "
    "--readFilesIn {input.fastq} "
    "--outFileNamePrefix {params.outprefix} "
    "--outSAMtype BAM Unsorted ; "
    # move STAR "progress" logs into log directory
    "mv data/{wildcards.sample}/star.Log.progress.out "
    "data/{wildcards.sample}/logs ; "
    "mv data/{wildcards.sample}/star.Log.out "
    "data/{wildcards.sample}/logs"

# sort aligned reads
rule sort_aligned:
  input:
    "data/{sample}/star.Aligned.out.bam"
  output:
    temp("data/{sample}/star.Aligned.sorted.bam")
  log:
    "data/{sample}/logs/sort_aligned.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard SortSam "
    "INPUT={input} "
    "OUTPUT={output} "
    "SORT_ORDER=queryname "
    "2> {log}"

# merge aligned and unaligned reads to add tags to aligned reads. this also removes secondary
# alignments!
rule merge_bam:
  input:
    aligned = "data/{sample}/star.Aligned.sorted.bam",
    unaligned = "data/{sample}/polyA_trimmed_unmapped.bam",
    reference = lambda wildcards: config["align_ref"][wildcards.sample] + "/cropseq_ref.fasta",
    dict = lambda wildcards: config["align_ref"][wildcards.sample] + "/cropseq_ref.dict"
  output:
    temp("data/{sample}/merged_aligned.bam")
  log:
    "data/{sample}/logs/merge_bam.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard MergeBamAlignment "
    "ALIGNED_BAM={input.aligned} "
    "UNMAPPED_BAM={input.unaligned} "
    "REFERENCE_SEQUENCE={input.reference} "
    "OUTPUT={output} "
    "INCLUDE_SECONDARY_ALIGNMENTS=false "
    "PAIRED_RUN=false "
    "2> {log}"

# tag reads with gene exons
rule tag_with_gene_exon:
  input:
    mapped_bam = "data/{sample}/merged_aligned.bam",
    annot = lambda wildcards: config["align_ref"][wildcards.sample] + "/cropseq_ref.refFlat"
  output:
    protected("data/{sample}/gene_tagged_aligned.bam")
  log:
    "data/{sample}/logs/tag_with_gene_exon.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "TagReadWithGeneExon "
    "INPUT={input.mapped_bam} "
    "OUTPUT={output} "
    "ANNOTATIONS_FILE={input.annot} "
    "TAG=GE "
    "2> {log}"
    
# create index file for final bam file
rule create_bam_index:
  input:
    "data/{sample}/gene_tagged_aligned.bam"
  output:
    "data/{sample}/gene_tagged_aligned.bai"
  log:
    "data/{sample}/logs/create_bam_index.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "picard BuildBamIndex "
    "INPUT={input} "
    "OUTPUT={output} "
    "2> {log}"

# calculate reads per cell barcode
rule reads_per_cell:
  input:
    "data/{sample}/gene_tagged_aligned.bam"
  output:
    "data/{sample}/reads_per_cell_barcode.txt"
  params:
    read_quality = config["reads_per_cell"]["read_quality"]
  log:
    "data/{sample}/logs/reads_per_cell.log"
  conda:
    "../envs/dropseq_tools.yml"
  shell:
    "BAMTagHistogram "
    "INPUT={input} "
    "OUTPUT={output} "
    "TAG=XC "
    "READ_QUALITY={params.read_quality} "
    "2> {log}"

# compile R Markdown report in html format
rule align_report:
  input:
    cell_bcs = "data/{sample}/cell_tags_summary.txt",
    mol_bcs = "data/{sample}/mol_tags_summary.txt",
    star_smry = "data/{sample}/star.Log.final.out",
    adapt_trim = "data/{sample}/adapter_trimming_report.txt",
    polyA_trim = "data/{sample}/polyA_trimming_report.txt",
    reads_per_cell = "data/{sample}/reads_per_cell_barcode.txt"
  output:
    "results/alignment/{sample}_align_report.html"
  params:
    expect_cells = lambda wildcards: config["expect_cell_numbers"][wildcards.sample],
    bc_structure = lambda wildcards: config["bc_structure"][wildcards.sample]
  conda:
    "../envs/r_dropseq_tools.yml"
  script:
    "../scripts/align_report.Rmd"
