getDS.wrapper <- function(experiment = "TAP", panel = "perSample", sampling = "fixed", reference = UseCells, reads = NULL) {
  
  if (experiment =="TAP")  folders <- c("TAP_CT_kitBM/downsampled/", "TAP_CT_totalBM/downsampled/") else folders <- c("WTX_CT_kitBM/downsampled/","WTX_CT_totalBM/downsampled/")
  if (experiment =="TAP")  compareTo <- TAPseq else compareTo <- UseCells
  
  
  folders <- file.path(dir, folders)
  if (is.null(reads)) reads <- as.integer(unique(gsub("ds_dge_(\\d+)_reads_.+","\\1",list.files(folders, sprintf("%s_%s",sampling,panel)))))
  kappas <- lapply(reads, function(r) {
    if (r == 50000) return(NULL)
    #if (experiment == "TAP" & !(sampling =="fixed"))nc <- 4000 else nc <- NULL
    projected.id <- getDS(folders, panel, sampling, r, reference = reference, minFeature = ifelse(experiment == "TAP",10, min(c(500, round(r / 20 )))),ncells = NULL)
    
    compVector <- Idents(compareTo) 
    projected.id <- projected.id[names(projected.id) %in% names(compVector)]
    compVector <- compVector[names(projected.id)]
    alln <- names(projected.id)
    
    projected.id <- paste0("A",projected.id); names(projected.id) <- alln
    compVector <- paste0("B",compVector); names(compVector) <- alln
    
    precision_by_class <- sapply(unique(compVector), function(x) {
      cells <- compVector == x
      table(factor(projected.id[cells], level = unique(projected.id)))
    })
    precision_by_class.pl <- t(t(precision_by_class) / apply(precision_by_class,2,sum))
 
    k <- adjustedRandIndex(projected.id, compVector)
    data.frame(experiment = experiment, panel = panel, sampling = sampling, reads = r, kappa = k , ncells = length(projected.id))
  })
  do.call(rbind,kappas)
  
}

getDS <- function(folders, panel, sampling, r, reference = UseCells, output = "prediction", minFeature = 10, ncells = NULL) {
  cat(panel, sampling, r, "\n")
  
  DGE <-
    lapply(
      file.path(folders, sprintf("ds_dge_%d_reads_%s_%s.txt",r, sampling,panel)),
      read.table,
      header = T,
      row.names = 1
    )
  names(DGE) <- sapply(strsplit(folders,"/"), function(x) x[length(x)-1])
  DGE <-
    lapply(names(DGE), function(n) {
      x <- DGE[[n]]
      colnames(x) <- paste(n, colnames(x), sep = "-")
      x
    })
  
  common.genes <- rownames(DGE[[1]])
  if (length(DGE) > 1)
    for (i in 2:length(DGE))
      common.genes <- intersect(common.genes, rownames(DGE[[i]]))
  
  DGE <- lapply(DGE, function(x)
    x[common.genes,])
  DGE <- do.call(cbind, DGE)
  
  
  DGE.seurat <- CreateSeuratObject(DGE, min.features = minFeature)
  if (!is.null(ncells)) DGE.seurat <- subset(DGE.seurat, cells = sample(colnames(DGE.seurat),min(c(ncells,ncol((DGE.seurat)))), replace = F))
  DGE.seurat <- NormalizeData(object = DGE.seurat)
  DGE.seurat <- FindVariableFeatures(object = DGE.seurat)
  DGE.seurat <- ScaleData(object = DGE.seurat)
  DGE.seurat <- RunPCA(object = DGE.seurat)
  DGE.seurat <- FindNeighbors(DGE.seurat)
  DGE.seurat <- FindClusters(DGE.seurat)
  
  if (output == "prediction") Idents(DGE.seurat) else DGE.seurat
  
}
