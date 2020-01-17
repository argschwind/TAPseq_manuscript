#big run: one test (MAST?)
#many different numbers of cells
#many different numbers of reads average
.libPaths("/g/steinmetz/velten/Software/RLibs-seurat3/")
require(Seurat)
require(ggplot2)
require(plyr)
require(Matrix)
require(SingleCellExperiment)
require(gamlss)
require(VGAM)
require(maxLik)
require(MAST)


#script takes as input: the number of cells and the number of reads desired, as well as whether whole tx or TAP
input <- commandArgs(trailingOnly = T)
ncells <- as.integer(input[1])
nreads <- as.integer(input[2])
assay <- input[3]
test <- input[4]
OnPanel <- as.logical(input[5])


#other options are hardcoded
norm <- function(x) sum(x[x < quantile(x, probs = 0.9)]) + 1
normcounts <- T
nrepeats <- 5 #was 100
maxcells <- 150
FDR <- T

TAP.nods <- readRDS(url("http://steinmetzlab.embl.de/TAPdata//TAP.nods.RDS"))
TAP.nods <- TAP.nods[!grepl("CROPseq",rownames(TAP.nods)),]

knowns <- c("HBG1","HBG2","HBE1","HBD","HBB","GATA1","ZFPM2","MYC")
#panels.in <- read.csv("/g/steinmetz/project/singcellTxn/CRISPRdrop/k562_enhancer_cropseq/tapseq_manuscript/meta_data/target_genes_validation.csv", stringsAsFactors = F)
panels.in <- read.csv(url("http://steinmetzlab.embl.de/TAPdata/target_genes_validation.csv"), stringsAsFactors = F)

panel <- panels.in$gene[grepl("chr8",panels.in$panel) | panels.in$gene %in% knowns]


#data in/out
if (assay == "TAP") {
  tap1.con <- gzfile(sprintf("../data/TAP_DIFFEX_sample1/downsampled/ds_dge_%d_reads_average_perSample.txt.gz",nreads))
  TAP1 <- read.table(tap1.con, header=T, row.names = 1)
  tap2.con <- gzfile(sprintf("../data/TAP_DIFFEX_sample2/downsampled/ds_dge_%d_reads_average_perSample.txt.gz",nreads))
  TAP2 <- read.table(tap2.con, header=T, row.names = 1)
  colnames(TAP1) <- paste0("TAP1-",colnames(TAP1))
  colnames(TAP2) <- paste0("TAP2-",colnames(TAP2))
  TAP1 <- TAP1[intersect(rownames(TAP1),rownames(TAP2)),]
  TAP2 <- TAP2[intersect(rownames(TAP1),rownames(TAP2)),]
  
  DGE <- cbind(TAP1,TAP2)
  PER <- t(readRDS(url("http://steinmetzlab.embl.de/TAPdata//TAP.per.RDS")))
  rm(TAP1,TAP2)
} else {
  PER <- t(readRDS(url("http://steinmetzlab.embl.de/TAPdata//Whole.per.RDS")))
  

    #DGE <- readRDS(sprintf("data/ds_dge_%d_reads_average_perSample.rds",nreads))
    whole.con <- gzfile(sprintf("/g/steinmetz/project/singcellTxn/CRISPRdrop/k562_enhancer_cropseq/tapseq_manuscript/data/WholeTx/downsampled/ds_dge_%d_reads_average_perSample.txt.gz",nreads))
    DGE <- read.table(whole.con, header=T, row.names=1)
    PER <- PER[colnames(DGE),]
    if (onPanel) DGE <- DGE[panel,]
}





scramble.cols <- grep("non-targeting",colnames(PER))
is.scramble <- apply(PER[,scramble.cols],1,function(x) any(x > 0)) & apply(PER[,-scramble.cols],1,function(x) all(x == 0)) # <<- CHECK!
source("functions_runSeuratTest.R")

DGE <- DGE[!grepl("CROPseq",rownames(DGE)),]
cell.per.per <- apply(PER,2,sum)
useper <- colnames(PER)[cell.per.per > max(c(maxcells, ncells)) & !grepl("non-target", colnames(PER))]
#run on guide level first

if (test == "MAST.cov") covariate <- apply(DGE>0,2,sum) else covariate <- NULL

result <- lapply(useper, runSeuratTest,test=test, DGE = DGE, pert= PER, scrcols = scramble.cols,ncell = ncells, normfun =  norm,repeats =nrepeats,normcounts = normcounts, covariate = covariate,returnRaw = FDR, allgenes = rownames(TAP.nods))


if (FDR) result <- do.call(rbind,lapply(result,"[[","raws")) else result <- do.call(rbind, result)
result$reads <- nreads

saveRDS(result, file = sprintf("runs/%s.%s/%s_%d_reads_%d_cells.RDS",test,ifelse(OnPanel,"targeted","wholegenome"),assay,nreads,ncells))
