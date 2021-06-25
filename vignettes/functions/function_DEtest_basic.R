

runSeuratTest <- function(g = "DS42_eGATA1_D", DGE, pert, scrcols, scalefun= log1p, normfun = function(x) 1,normcounts=T, covariate) {
  
  is.scramble <- apply(pert[,scrcols],1,function(x) any(x > 0)) & apply(pert[,-scrcols],1,function(x) all(x == 0)) # <<- CHECK!
  is.genes <- !grepl("CROPseq", rownames(DGE))
  
  object <- CreateSeuratObject(DGE[is.genes,]) #normalization causes problems, only log-transform
  sf <-  apply(as.matrix(DGE[is.genes,]) , 2, normfun)
  object <- SetAssayData(object, new.data = scalefun(t(t(as.matrix(DGE[is.genes,])) /sf )))
  if (normcounts) object <- SetAssayData(object, new.data = round(t(t(as.matrix(DGE[is.genes,])) /sf ) * mean(sf)), slot = "counts")  else object <- SetAssayData(object, new.data = as.matrix(DGE[is.genes,]), slot = "counts") 
  
  #is.scramble <- apply(pert[,1:4],1,function(x) any(x > 0)) & apply(pert[,5:ncol(pert)],1,function(x) all(x == 0)) # <<- CHECK! and remove
  Idents(object) <- ifelse(pert[,g] == 1, g, ifelse(is.scramble, "scramble", "other"))
  object$covariate <- covariate
  a <- colnames(object)[Idents(object) == g]
  b <- colnames(object)[Idents(object) == "scramble"]
  
  use.a <- a
  use.b <- b
  testobject <- subset(object, cells = c(use.a, use.b))
  
  a.obj <- subset(object, cells = use.a)
  b.obj <- subset(object, cells = use.b)
  
  rowmeans.a <- apply(GetAssayData(a.obj, slot="counts"),1,mean)
  rowmeans.b <- apply(GetAssayData(b.obj, slot="counts"),1,mean)
  abs.diff <- rowmeans.a - rowmeans.b
  
  source("functions/MAST.R", local = TRUE)
  
  sce <- SingleCellExperiment(list(counts =GetAssayData(testobject, slot = "counts")  ) , colData = data.frame(row.names = colnames(testobject), ident = Idents(testobject), cov = testobject$covariate))
  MAST.result <- de_MAST(sce)
  m <- data.frame(row.names =MAST.result$gene, p_val = MAST.result$pvalue, p_val_adj = p.adjust( MAST.result$pvalue, method = "fdr"), avg_logFC = MAST.result$logFC)
  
  #m <- FindMarkers(testobject, ident.1 = g, ident.2 = "scramble", test.use = "MAST", logfc.threshold =0 , min.pcr = 0)  
  m$guide <- g
  m$gene <- rownames(m)
  m$ncells <- length(a)
  m$absdiff <- abs.diff[as.character(m$gene)]
  m
}