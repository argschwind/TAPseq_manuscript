runSeuratTest <- function(g = "DS42_eGATA1_D", test = "wilcox", DGE, pert, scrcols, scalefun= log1p, normfun = function(x) 1, repeats = 100, ncell = 50, ncell.scramble = 500, fdr.thresh = 0.1, normcounts=T, covariate = NULL, returnRaw = F, allgenes = NULL) {
  
  is.scramble <- apply(pert[,scrcols],1,function(x) any(x > 0)) & apply(pert[,-scrcols],1,function(x) all(x == 0)) # <<- CHECK!
  is.genes <- !grepl("CROPseq", rownames(DGE))
  
  object <- CreateSeuratObject(DGE[is.genes,]) #normalization causes problems, only log-transform
  sf <-  apply(as.matrix(DGE[is.genes,]) , 2, normfun)
  object <- SetAssayData(object, new.data = scalefun(t(t(as.matrix(DGE[is.genes,])) /sf )))
  if (normcounts) object <- SetAssayData(object, new.data = round(t(t(as.matrix(DGE[is.genes,])) /sf ) * mean(sf)), slot = "counts")  else object <- SetAssayData(object, new.data = as.matrix(DGE[is.genes,]), slot = "counts") 
  
  #is.scramble <- apply(pert[,1:4],1,function(x) any(x > 0)) & apply(pert[,5:ncol(pert)],1,function(x) all(x == 0)) # <<- CHECK! and remove
  Idents(object) <- ifelse(pert[,g] == 1, g, ifelse(is.scramble, "scramble", "other"))
  if (!is.null(covariate)) object$covariate <- covariate
  a <- colnames(object)[Idents(object) == g]
  b <- colnames(object)[Idents(object) == "scramble"]
  if (ncell > length(a)) return(NULL)
  # TP.set <-  subset(summary.total, guide == g & p <= 0.05)$gene
  # TN.set <-  subset(summary.total, guide == g &( p > 0.05 | is.na(p)))$gene
  #parfun <- ifelse(test == "DEsingle", lapply, mclapply)
  tests <- lapply(1:repeats, function(i) {
    use.a <- sample(a, size =  ncell, replace=F)
    use.b <- sample(b, size =  ncell.scramble, replace=F)
    testobject <- subset(object, cells = c(use.a, use.b))
    
    
    
    if (test == "DEsingle") {
      desingle.result <- DEsingle(counts = GetAssayData(testobject, slot = "counts"), group = Idents(testobject))
      m <- data.frame(row.names = rownames(desingle.result), p_val = desingle.result$pvalue, p_val_adj = desingle.result$pvalue.adj.FDR, avg_logFC = -log2(desingle.result$foldChange))
      
    } else if (test == "scDD") {
      
      sce <- SingleCellExperiment(list(normcounts =GetAssayData(testobject)  ) , colData = data.frame(row.names = colnames(testobject), ident = as.integer(Idents(testobject))))
      prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
      scdd.result <- scDD(sce, prior_param=prior_param, testZeroes=T, param=BiocParallel::MulticoreParam(workers=1),categorize=F, condition = "ident")
      scdd.result <- metadata(scdd.result)$Genes
      m1 <- FindMarkers(testobject, ident.1 = g, ident.2 = "scramble", test.use = "MAST", logfc.threshold =0 , min.pcr = 0)  
      m <- na.omit(data.frame(row.names = rownames(scdd.result), p_val = scdd.result$combined.pvalue, p_val_adj = scdd.result$combined.pvalue.adj, avg_logFC = m1[rownames(scdd.result),"avg_logFC"]))
      
    } else if (test == "MAST.cov") {
      source("functions/MAST.R", local = TRUE)
      if (is.null(covariate)){
      featurevar <- apply(GetAssayData(testobject),1,var)
      VariableFeatures(testobject) <- names(featurevar)[featurevar >0 ]
      testobject <- ScaleData(testobject)
      testobject <- RunPCA(testobject, pc.genes = rownames(testobject), npcs = 2 )
      covariate <- Embedding(testobject, "pca")[,1]
      } else covariate <- testobject$covariate
      sce <- SingleCellExperiment(list(counts =GetAssayData(testobject, slot = "counts")  ) , colData = data.frame(row.names = colnames(testobject), ident = Idents(testobject), cov = covariate))
      MAST.result <- de_MAST(sce)
      m <- data.frame(row.names =MAST.result$gene, p_val = MAST.result$pvalue, p_val_adj = p.adjust( MAST.result$pvalue, method = "fdr"), avg_logFC = MAST.result$logFC)
      
    } else {
      m <- FindMarkers(testobject, ident.1 = g, ident.2 = "scramble", test.use = test, logfc.threshold =0 , min.pcr = 0)  
    }
    
    #compute ROC
    # t <- gsub("DS.+_e([^_]+)_.","\\1",g)
    # if (t =="HS2") t <- c("HBG2","HBG1","HBE1")
    # useforroc <- m; useforroc$gene <- rownames(useforroc)
    # useforroc <- subset(useforroc, !gene %in% as.character(TP.set) | gene %in% t)
    # predob <- prediction(1-useforroc$p_val, labels = useforroc$gene %in% t)
    # auc <- performance(predob, measure = "auc")
    hits <- rownames(m)[m$p_val_adj <= fdr.thresh & m$avg_logFC < 0 ]
    nonhits <- rownames(testobject)[! rownames(testobject) %in% hits]
    #precision.tpr = sum(hits %in% TP.set) / length(hits),
    #recall = sum(hits %in% TP.set) / length(TP.set),
    #specificity.tnr = sum(nonhits %in% TN.set) / length(TN.set),
    #fpr = sum(hits %in% TN.set) / length(hits),
    #fnr = sum(nonhits %in% TP.set) / length(nonhits),

    notused <- allgenes[!allgenes %in% rownames(m)]
    raw <- data.frame(guide =g, run =i , test = test, ncells= ncell,norm = substr(paste(deparse(normfun),collapse= " "),1,500),
                      scale = substr(paste(deparse(scalefun),collapse= " "),1,500), gene = c(rownames(m),notused), p_val = c(m$p_val, rep(1,length(notused))), p_val_adj = c(m$p_val_adj, rep(1,length(notused))))
    
    list(summarized = data.frame(guide =g,
               run = i,
               test = test,
               ncells = ncell,
               norm = substr(paste(deparse(normfun),collapse= " "),1,500),
               scale = substr(paste(deparse(scalefun),collapse= " "),1,500),
               hits = paste(hits, collapse = ";")),
         raw = raw)
  })
  tests2 <- do.call(rbind,lapply(tests,"[[","summarized"))
  raws <- do.call(rbind,lapply(tests,"[[","raw"))
  if (returnRaw) list(summarized = tests2, raws = raws) else tests2
}