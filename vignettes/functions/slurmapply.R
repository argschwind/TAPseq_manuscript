getAttachedPackages <- function() {
  toattach <- search()
  toattach <- toattach[grepl("package:",toattach)]
  toattach <- gsub("package:", "", toattach)
  toattach
}


slurmlapply <- function(X, FUN, ..., s.walltime = "1:00:00", s.memory = NULL, s.ncores = 1, s.tmpdir = "/scratch/velten/slurmapply", s.collect = TRUE, s.pckgs = getAttachedPackages(), s.module= "") {
  args <- list(...)
  paths <- .libPaths()
  tmpfile <- tempfile(tmpdir = s.tmpdir)
  tmpfile.o <- tempfile(tmpdir = s.tmpdir)
  tmpfile.e <- tempfile(tmpdir = s.tmpdir)
  
  cat("Writing environment to ", tmpfile,"\n")
  cat("Writing output to ", tmpfile.o,"\n")
  cat("Writing error to ", tmpfile.e,"\n")
  save(args, s.pckgs, FUN, X, paths, file = tmpfile)
  ran <- c()
  for (i in 1:length(X)) {
    params <- c(sprintf("#SBATCH -t %s",s.walltime), sprintf("#SBATCH -n %d",s.ncores), sprintf("#SBATCH -e %s.%d", tmpfile.e,i), sprintf("#SBATCH -o %s.%d",tmpfile.o, i))
    if (!is.null(s.memory))params <- c(params, sprintf("#SBATCH --mem %d", s.memory))
    
    command <- sprintf("Rscript /g/steinmetz/project/singcellTxn/CRISPRdrop/LV/Vignettes/functions/slurm.R %s %d", tmpfile, i)
    tmpfile.s <- tempfile(tmpdir = s.tmpdir)
    writeLines(c("#!/bin/bash",params,s.module,command), con=tmpfile.s)
    ran <- c(ran,system(sprintf("sbatch  %s",tmpfile.s),intern=T))
  }
  
  jobs <- gsub("Submitted batch job ","", ran)
  
  if (!s.collect) return(list(jobs = jobs, tmpfile = tmpfile, tmpfile.e = tmpfile.e) ) else return(slurmcollect(jobs, tmpfile, tmpfile.e))
}

slurmcollect <- function(jobs, tmpfile, tmpfile.e) {
  
  status <- rep("RUNNING",length(jobs))
  result <- list()
  while(any(status == "RUNNING") | any(status == "PENDING")) {
    for (j in 1:length(jobs)) {
      if(file.exists(sprintf("%s-%d",tmpfile,j))) {
        status[j] <- "COMPLETED"
        load(sprintf("%s-%d",tmpfile,j))
        #cat("Collected ", j, "\n")
        result[[j]] <- out
      } else {
      cn <- pipe(sprintf("sacct -j %s", jobs[j]))
      est <- read.fwf(cn, width=c(12,11,11,11,11,11,10), skip=2)
      if (!any(grepl("RUNNING",est[,6])) & !any(grepl("PEND",est[,6]))) {
          status[j] <- "FAILURE"
          warning("Error in core ", j)
          errorMsg <- readLines(sprintf("%s.%d", tmpfile.e,j))
          warning(tail(errorMsg))
          result[[j]] <-NULL
      }
    }    
    
    }
    if (any(status=="RUNNING")) {
      Sys.sleep(15)
      cat("Cores ", which(status=="RUNNING"), " are running.\n")
    } 
  }
  return(result)
}
