# Analysis code for "Targeted Perturb-seq enables genome-scale genetic screens in single cells"

This repository contains code to reproduce the results presented in the article. It builds upon a
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html#) workflow that handles all data
processing from fastq files to transcript counts etc. It also performs some high-level analyses, in 
particular analyses underlying figure 3 of the manuscript and quality control steps. 

High level analyses underlying figures 1, 2, 3b,h,i and 4 are described in detail in the *vignettes*
subfolder. Input data required for running these vignettes can either be produced by running the
relevant parts of the snakemake pipeline, or they can be downloaded as described in the vignettes.  

The workflow can be downloaded by simply cloning the repository into a location of choice:
```
git clone https://github.com/argschwind/tapseq_manuscript.git
```

## Download data
**This is currently not implemented, as GEO links are not active yet!**
Download of raw data from GEO is handled by the pipeline and necessary fastq files are downloaded
automatically when executing the workflow. All fastq files are saved under raw_data.

## Data processing
All data processing and analyses steps can be executed through snakemake and conda. This only
requires that conda and snakemake are installed. All other required dependencies will be installed
through conda.

Data processing steps (and part of the analyses) can be executed by running parts of the workflow.
This will create the subdirectories data and results. Data contains all processed data in one
subdirectory per sample. All executed analyses will be saved as .html output in results. Alignment
reports for all samples will for instance be saved under
results/alignment/'sample'_align_report.html.

However, some of these steps might take a very long time and use substantial free storage!

### Reproducing analyses 
Input data for the following thematic analyses can be procuced by executing following snakemake
rules.

#### TAP-seq quality control (Figure 1)

Key plots of figure 1 can be reproduced by running the Vignette for Figure1e,f. This requires the
following data processing step:
```
# align all TAP-seq data for figure 1 and create dge downsamplings (--jobs = number of threads to
# use in parallel, please adjust; -n = dryrun, remove it to execute)
snakemake --use-conda Figure1 --jobs 4 -n
```

To obtain all data used for quality control, i.e. also the experiments detailed in supplementary
figure 2-4 and the supplementary note, the following data processing step is required:
```
# create all input data for supplementary figures 2-4
snakemake --use-conda tapseq_validation --jobs 4 -n
```

#### Evaluation of differential expression performance (Figure 2)

The analyses underlyng figure 2 can be reproduced using the corresponding vignette. This requires
the following data processing step; furthermore, it requires the execution of a relatively
compute-expensive sampling and DE testing workflow, documented in the vignette.
```
# create all input data for figure 2
snakemake --use-conda Figure2 --jobs 4 -n
```


#### Enhancer screen (Figure 3)
Data and most analyses for the chromosome 8 and 11 enhancer screen can be reproduced with this
command.
```
snakemake --use-conda enhancer_screen --jobs 4 -n
```

Code to reproduce the enhancer prediction analyses can be found in the Vignette for Figure 3b,h,i.
This requires generating chromatin annotated enhancer - target pairs. This can be achieved by
running the following snakemake rule:
```
snakemake --use-conda chromatin_annotated_etps --jobs 4 -n
```

#### Bone marrow cell type identification (Figure 4)
Input data for mouse bone-marrow cell type identification analyses can be produced by this command.
Downstream analyses to repoduce plots shown in the article are found in the Vignette for Figure 4.
```
snakemake --use-conda bone_marrow_cell_types --jobs 4 -n
```

### Executing different steps of data processing
Following commands can be used to execute data processing steps for **all samples**.
```
# generate all required alignment references
snakemake --use-conda alignment_reference -n

# align reads for all samples
snakemake --use-conda align -n

# create digital gene expression data for all samples
snakemake --use-conda dge -n
```

### Executing full workflow
If you feel brave you can also execute the entire workflow at once. Godspeed, this will take a long
time!
```
# execute entire project (-k: don't abort independent jobs if a job fails)
snakemake --use-conda --jobs 4 -k -n
```
