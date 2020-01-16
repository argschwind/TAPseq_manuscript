## create chromatin annotated enhancer - target pairs for predictive models 

### input, output and shell paths are all relative to the project directory ###

rule enh_screen_etps:
  input:
    processed_results = "data/diff_expr_screen_nGenesCovar.csv",
    encode_bam = expand("data/k562_chromatin_data/{assay}_encode_rep1_alignment.bam",
      assay = config["chromatin_analyses"]["encode_bam"]),
    hic_raw = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.RAWobserved",
      chr = ["chr8", "chr11"]),
    hic_norm = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.KRnorm",
      chr = ["chr8", "chr11"])
  output:
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_screen.csv"
  threads: 5
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/chromatin_annotated_etps/enhancer_screen_etps.R"
    
rule genome_wide_etps:
  input:
    wholeTx_dge = "data/Sample10X/dge.txt",
    enhancers = "meta_data/allEnhancers.hg38.txt",
    encode_bam = expand("data/k562_chromatin_data/{assay}_encode_rep1_alignment.bam",
      assay = config["chromatin_analyses"]["encode_bam"]),
    hic_raw = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.RAWobserved",
      chr = ["chr" + str(i)  for i in range(1, 23)]),
    hic_norm = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.KRnorm",
      chr = ["chr" + str(i)  for i in range(1, 23)])
  output:
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_GW.csv"
  threads: 5
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/chromatin_annotated_etps/genome_wide_etps.R"
  
rule gasperini_enh_screen_etps:
  input:
    wholeTx_dge = "data/Sample10X/dge.txt",
    enhancers = "meta_data/gasperini_enhancers_screen.bed",
    encode_bam = expand("data/k562_chromatin_data/{assay}_encode_rep1_alignment.bam",
      assay = config["chromatin_analyses"]["encode_bam"]),
    hic_raw = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.RAWobserved",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]]),
    hic_norm = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.KRnorm",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]])
  output:
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_gasperini_screen.csv"
  threads: 5
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/chromatin_annotated_etps/gasperini_screen_etps.R"
    
rule gasperini_pilot_etps:
  input:
    wholeTx_dge = "data/Sample10X/dge.txt",
    enhancers = "meta_data/gasperini_enhancers_pilot.bed",
    encode_bam = expand("data/k562_chromatin_data/{assay}_encode_rep1_alignment.bam",
      assay = config["chromatin_analyses"]["encode_bam"]),
    hic_raw = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.RAWobserved",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]]),
    hic_norm = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.KRnorm",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]])
  output:
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_gasperini_pilot.csv"
  threads: 5
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/chromatin_annotated_etps/gasperini_pilot_etps.R"
    
rule fulco_tiling_etps:
  input:
    wholeTx_dge = "data/Sample10X/dge.txt",
    enhancers = "meta_data/fulco_tiling_enhancers.bed",
    encode_bam = expand("data/k562_chromatin_data/{assay}_encode_rep1_alignment.bam",
      assay = config["chromatin_analyses"]["encode_bam"]),
    hic_raw = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.RAWobserved",
      chr = ["chr8", "chrX"]),
    hic_norm = expand("data/k562_chromatin_data/HiC/5kb_resolution_intrachromosomal/{chr}/MAPQG0/{chr}_5kb.KRnorm",
      chr = ["chr8", "chrX"])
  output:
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_fulco_tiling.csv"
  threads: 5
  conda: "../envs/r_map_enhancers.yml"
  script:
    "../scripts/chromatin_annotated_etps/fulco_tiling_etps.R"

# create all chromatin annotated etps
rule chromatin_annotated_etps:
  input:
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_screen.csv",
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_GW.csv",
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_gasperini_screen.csv",
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_gasperini_pilot.csv",
    "data/chromatin_annotated_etps/chromatin_annotated_pairs_fulco_tiling.csv"
