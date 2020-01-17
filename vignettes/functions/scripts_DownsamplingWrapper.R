#big run: one test (MAST?)
#many different numbers of cells
#many different numbers of reads average

input <- commandArgs(trailingOnly = T)
TEMPDIR <- input[1]
module <- input[2]


ncell_range <- c(10,50,100,150)
reads.TAP <-   as.integer(gsub("ds_dge_(\\d+).+","\\1",list.files("../data/TAP_DIFFEX_sample1/downsampled", pattern = "ds_dge_\\d+_reads_average_perSample.txt.gz")))
reads.Whole <-   as.integer(gsub("ds_dge_(\\d+).+","\\1",list.files("../data/WTX_DIFFEX/downsampled", pattern = "ds_dge_\\d+_reads_average_perSample.txt.gz")))
reads.Whole.Genome <- c(250,2500,10000)


runSlurm <- function(cells, reads, assay, test, onPanel) {
  
  params <- c(sprintf("#SBATCH -t %s","36:00:00"), 
              sprintf("#SBATCH -n %d",1), 
              sprintf("#SBATCH -J %s_%s_%d_%d",test,assay,reads,cells),
              sprintf("#SBATCH -e runs/logs/%s_%s_%d_reads_%d_cells.err",test,assay,reads,cells), 
              sprintf("#SBATCH -o runs/logs/%s_%s_%d_reads_%d_cells.err",test,assay,reads,cells))
  if (assay == "Whole") params <- c(params, "#SBATCH --mem 20000")
  
  command <- sprintf("Rscript script_runDownsamplingTest.R %d %d %s %s %s", cells, reads, assay, test, as.character(onPanel))
  
  tmpfile.s <- tempfile(tmpdir = TEMPDIR)
  
  writeLines(c("#!/bin/bash",params,module,command), con=tmpfile.s)
  
  system(sprintf("sbatch  %s",tmpfile.s))
}




for (nc in ncell_range) {
  for (nr in reads.TAP) runSlurm(nc, nr, "TAP", "MAST.cov", TRUE)
}

for (nc in ncell_range) {
  for (nr in reads.Whole) runSlurm(nc, nr, "Whole", "MAST.cov", TRUE) 
}

for (nc in ncell_range) {
  for (nr in reads.Whole.Genome) runSlurm(nc, nr, "Whole", "MAST.cov", FALSE) 
}
