#---------------------------------------------------------------------
# Title:         IRIS - Functions
# Author:        Brandon Monier
# Created:       2018-01-26 at 11:31:20
# Last Modified: 2018-01-26 at 11:31:34
#---------------------------------------------------------------------

# HOUSE KEEPING FUNCTIONS

## CRAN Install
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(
      new.pkg, 
      dependencies = TRUE,
      repos="http://mirror.las.iastate.edu/CRAN/"
    )
  sapply(pkg, require, character.only = TRUE)
}

## Size factor extraction (for normalization)
getSizeFact <- function(rc.data) {
	geomMean <- function(x) { 
		prod(x)^(1 / length(x)) 
	}
	gm.mean <- apply(rc.data, 1, geomMean)
	gm.mean[gm.mean == 0] <- NA
	rc.data <- sweep(rc.data, 1, gm.mean, FUN = "/") 
	sf <- apply(rc.data, 2, median, na.rm = TRUE)
	return(sf)
}

## Get normalized counts from different object classes
getNormCounts <- function(rc.data) {
  if (class(rc.data) == "DESeqDataSet") {
    nc.data <- BiocGenerics::counts(rc.data, normalize = TRUE)
    return(nc.data)
  } else {
    sf <- getSizeFact(rc.data)
    nc.data <- sweep(rc.data, 2, sf, FUN = "/")
    return(nc.data)
  }
}

## Get normalized counts for selected gene (heatmap interactivity)
getGenes <- function(rc.data, id, coldata) {
  nc.data <- getNormCounts(rc.data)
  nc.data <- as.data.frame(nc.data[rownames(nc.data) == id, ])
  names(nc.data) <- "counts"
  dat.l <- list(coldata, nc.data)
  nc.data <- Reduce(
    merge, lapply(dat.l, function(x) data.frame(x, sample = row.names(x)))
  )
  return(nc.data)
}

## Get contrast tables from different object classes
getContTable <- function(
  de.genes, coef, cts, expset, design, fact, fact5, fact6) {
  if (class(de.genes) == "MArrayLM") {
    de.genes2 <- topTable(
      fit = de.genes,
      coef = coef,
      number = nrow(cts)
    )
    de.genes2 <- de.genes2[order(rownames(de.genes2)), ]
    de.genes2 <- as.data.frame(de.genes2)
    de.genes2$baseMean <- rowMeans(cts)
    names(de.genes2) <- c(
      "log2FoldChange",
      "avgexpr",
      "t",
      "pvalue",
      "padj",
      "B",
      "baseMean"
    )
    de.genes2 <- subset(de.genes2, baseMean != 0)
    de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
    names(de.genes2)[1] <- "id"    
  } else if (class(de.genes) == "DGEGLM") {
    if (expset == "exp1" | expset == "exp2") {
      de.genes2 <- glmLRT(
        glmfit = de.genes,
        contrast = design[, coef]
      )
      de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "logCPM",
        "LR",
        "pvalue",
        "padj",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"
    } else if (expset == "exp3" | expset == "exp4") {
      de.genes2 <- glmLRT(
        glmfit = de.genes,
        coef = coef
      )
      de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "logCPM",
        "LR",
        "pvalue",
        "padj",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"      
    } else if (expset == "exp5" | expset == "exp6") {
      de.genes2 <- glmLRT(
        glmfit = de.genes,
        contrast = design[, coef]
      )
      de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "logCPM",
        "LR",
        "pvalue",
        "padj",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"
    } else if (expset == "exp7") {
      de.genes2 <- glmLRT(
        glmfit = de.genes,
        coef = coef
      )
      de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "logCPM",
        "LR",
        "pvalue",
        "padj",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()      
    }
  } else if (class(de.genes) == "DESeqDataSet") {
    if (expset == "exp1") {
      coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
      coef2 <- unlist(coef2)
      de.genes2 <- results(
        object = de.genes,
        contrast = c(fact, coef2[1], coef2[2])
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"      
    } else if (expset == "exp2") {
      coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
      coef2 <- unlist(coef2)
      de.genes2 <- results(
        object = de.genes,
        contrast = c("group", coef2[1], coef2[2])
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"         
    } else if (expset == "exp3" | expset == "exp4") {
      names(mcols(de.genes))[grep(
        "log2 fold change", 
        mcols(mcols(de.genes))$description
      )] <- colnames(design)
      de.genes2 <- results(
        object = de.genes,
        contrast = list(coef)
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"         
    } else if (expset == "exp5") {
      coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
      coef2 <- unlist(coef2)
      de.genes2 <- results(
        object = de.genes,
        contrast = c(fact5, coef2[1], coef2[2])
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
    } else if (expset == "exp6") {
      coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
      coef2 <- unlist(coef2)
      de.genes2 <- results(
        object = de.genes,
        contrast = c(fact6, coef2[1], coef2[2])
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
    } else if (expset == "exp7") {
      de.genes2 <- results(
        object = de.genes,
        name = coef
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
    }
  }
  return(de.genes2)
}



# LIMMA RETURN FIT functions

## LIMMA - EXP 1 - two group comparisons
limma.exp1 <- function(fact, coldata, cts, perm.h) {
  design <- model.matrix(~ 0 + coldata[, fact])
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  v <- voom(cts, design)
  fit <- lmFit(v)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  return(list(fit.cont, cont))  
}

## LIMMA - EXP 2 - multiple factor comparisons
limma.exp2 <- function(fact1, fact2, coldata, cts, perm.h) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  design <- model.matrix(~ 0 + group.c)
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(group.c)
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  v <- voom(cts, design)
  fit <- lmFit(v)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  return(list(fit.cont, cont)) 
}

## LIMMA - EXP 3 - classical interactions
limma.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[-1] <- paste0(colnames(tmp1)[-1], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)
  v <- voom(cts, design)
  fit <- lmFit(cts, design)
  fit.cont <- eBayes(fit)
  return(list(fit.cont, design))
}

## LIMMA - EXP 4 - added effects blocking and paired
limma.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[-1:-f1n] <- paste0(
    colnames(design)[-1:-f1n],
    "_VS_", 
    fact2.rlvl
  )
  v <- voom(cts, design)
  fit <- lmFit(cts, design)
  fit.cont <- eBayes(fit)
  return(list(fit.cont, design))
}

## LIMMA - EXP 7 - user inputs
limma.exp7 <- function(cts, mod.matrix) {
  design <- mod.matrix
  fit <- lmFit(cts, design)
  fit <- eBayes(fit)
  return(list(fit, design))
}



# EDGER RETURN FIT functions

## EDGER - EXP1 - two group comparisons
edger.exp1 <- function(fact, coldata, cts, perm.h, norm) {
  design <- model.matrix(~ 0 + coldata[, fact])
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont)) 
}

## EDGER - EXP2 - multiple factor comparisons
edger.exp2 <- function(fact1, fact2, coldata, cts, perm.h, norm) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  design <- model.matrix(~ 0 + group.c)
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(group.c)
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont)) 
}

## EDGER - EXP3 - classical interactions
edger.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl, 
                       norm) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[-1] <- paste0(colnames(tmp1)[-1], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, design))
}

## EDGER - EXP4 - added effects blocking and paired
edger.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl,
                       norm) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[-1:-f1n] <- paste0(
    colnames(design)[-1:-f1n],
    "_VS_", 
    fact2.rlvl
  )
  # design <- cbind(tmp1, tmp2, tmp3)
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, design))
}

## EDGER - EXP5 - main effects only
edger.exp5 <- function(fact, fact.levl, cts, coldata, perm.h, norm) {

  coldata[, fact] <- relevel(x = coldata[, fact], ref = fact.levl)

  design <- model.matrix(~ 0 + coldata[, fact])
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h

  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  # dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)

  return(list(fit.edger, cont))
}

## EDGER - EXP6 - main effects + grouping factor
edger.exp6 <- function(
  me.fact, me.levl, gp.fact, gp.levl, cts, coldata, perm.h, norm) {
  
  cts <- cts[, which(coldata[, gp.fact] == gp.levl)]
  coldata <- coldata[which(coldata[, gp.fact] == gp.levl), ]
  
  coldata[] <- lapply(coldata, function(x) if(is.factor(x)) factor(x) else x)
  coldata[, me.fact] <- relevel(x = coldata[, me.fact], ref = me.levl)
  
  design <- model.matrix(~ 0 + coldata[, me.fact])
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(coldata[, me.fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h

  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  # dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)

  return(list(fit.edger, cont))
}

## EDGER - EXP7 - user input
edger.exp7 <- function(cts, mod.matrix) {
  design <- mod.matrix
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  return(list(fit, design))
}




# DESeq2

## DESeq2 - EXP1 - Two group comparisons
deseq.exp1 <- function(fact, coldata, cts, perm.h) {
  design0 <- model.matrix(~ 0 + coldata[, fact])
  rownames(design0) <- rownames(coldata)
  colnames(design0) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h

  design <- paste("~", fact)
  design <- as.formula(design)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design
  )
  dds <- DESeq(dds)

  return(list(dds, cont0))
}

## DESeq2 - EXP2 - multiple factor comparisons
deseq.exp2 <- function(fact1, fact2, coldata, cts, perm.h) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  design <- model.matrix(~ 0 + group.c)
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(group.c)
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ 1
  )
  dds$group <- group.c
  design(dds) <- ~ group
  dds <- DESeq(dds)
  return(list(dds, cont))
}

## DESeq2 - EXP3 - classical interactions
deseq.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[-1] <- paste0(colnames(tmp1)[-1], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)

  design.dds <- paste("~", fact1, "*", fact2)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds)
  return(list(dds, design))
}

## DESeq2 - EXP4 - added effects blocking and paired
deseq.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[-1:-f1n] <- paste0(
    colnames(design)[-1:-f1n],
    "_VS_", 
    fact2.rlvl
  )
  # design <- cbind(tmp1, tmp2, tmp3)
  design.dds <- paste("~", fact1, "+", fact2)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds)
  return(list(dds, design))
}

## DESeq2 - EXP5 - main effects only
deseq.exp5 <- function(fact, fact.levl, cts, coldata, perm.h) {
  coldata[, fact] <- relevel(x = coldata[, fact], ref = fact.levl)
  design0 <- model.matrix(~ 0 + coldata[, fact])
  rownames(design0) <- rownames(coldata)
  colnames(design0) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h
  design.dds <- paste("~", fact)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  return(list(dds, cont0))
}

## DESeq2 - EXP6 - main effects only + grouping factor
deseq.exp6 <- function(
  me.fact, me.levl, gp.fact, gp.levl, cts, coldata, perm.h) {
  cts <- cts[, which(coldata[, gp.fact] == gp.levl)]
  coldata <- coldata[which(coldata[, gp.fact] == gp.levl), ]
  coldata[] <- lapply(coldata, function(x) if(is.factor(x)) factor(x) else x)
  coldata[, me.fact] <- relevel(x = coldata[, me.fact], ref = me.levl)
  design0 <- model.matrix(~ 0 + coldata[, me.fact])
  rownames(design0) <- rownames(coldata)
  colnames(design0) <- levels(coldata[, me.fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h
  design.dds <- paste("~", me.fact)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  return(list(dds, cont0))
}

## DESeq2 - EXP7 - user input
deseq.exp7 <- function(cts, coldata, mod.matrix) {
  design <- mod.matrix
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ 1
  )
  dds <- DESeq(dds, full = design, modelMatrixType = "standard")
  return(list(dds, design))
}





# PLOT Functions

## Bicluster plots
bicPlot <- function(n, res, cts.var) {
  n <- as.numeric(n)
  par(mar = c(10, 4, 3, 5) + 0.1, yaxt = "n")
  quheatmap(
    x = cts.var,
    bicResult = res,
    number = n, 
    showlabel = TRUE
  )  
}

## Sample Distance Matrix
sampdistPlot <- function(cts) {
  cts <- assay(cts)
  sampledists <- dist(t(cts))
  sdm <- as.matrix(sampledists)
  pheatmap(sdm)
}


## QC PLOTS

### Label conversion
qcLabConvert <- function(lab, tran) {
  if (tran == "log") {
    lab <- expression(paste("log"[2], " (counts + 1)"))
    return(lab)
  } else {
    return(lab)
  }
}


### Boxplot
qcBoxPlot <- function(tmp, lab, tran) {
  tmp <- assay(tmp)
  box <- as.data.frame(tmp)
  box <- tidyr::gather(box)
  
  lab <- qcLabConvert(lab, tran)
  
  p <- ggplot(box, aes(key, value)) +
    geom_boxplot() +
    xlab("") +
    ylab(lab) +
    ggtitle("Count data distributions - box and whisker") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
}


### Histogram
qcHist <- function(tmp, lab, tran) {
  tmp <- assay(tmp)
  hist <- as.data.frame(tmp)
  hist <- tidyr::gather(hist)
  
  lab <- qcLabConvert(lab, tran)
  
  p <- ggplot(hist, aes(value, color = key)) +
    geom_freqpoly() +
    xlab(lab) +
    ylab("Frequency") +
    ggtitle("Count data distributions - histogram") +
    theme_light() +
    scale_color_discrete(name = "Sample")
  print(p)
}


### Barplot
qcBarplot <- function(tmp) {
  tmp <- assay(tmp)
  bar <- as.data.frame(tmp)
  bar <- colSums(bar)
  bar <- as.data.frame(t(bar))
  bar <- gather(bar)
  
  p <- ggplot(bar, aes(key, value)) +
    geom_bar(stat = "identity", color = "#3d3d3d", fill = "#dddddd") +
    xlab("") +
    ylab("Counts") +
    ggtitle("Total reads") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
}


### PCA plot
qcPCAPlot <- function(tmp, pcafact) {
  pca.lab <- plotPCA(
    tmp,
    intgroup = pcafact
  )
  pca <- plotPCA(
    tmp,
    intgroup = pcafact,
    returnData = TRUE
  )
  p <- ggplot(pca, aes(PC1, PC2)) +
    geom_point(aes(shape = group), size = 2) +
    xlab(pca.lab$labels$x) +
    ylab(pca.lab$labels$y) +
    ggtitle("Principle Component Analysis") +
    theme_light() +
    scale_shape_discrete(name = pcafact)
  print(p)
}


### MDS plot
qcMDSPlot <- function(tmp, mdsfact) {
  mds <- dist(t(assay(tmp)))
  mds <- as.matrix(mds)
  mds <- as.data.frame(colData(tmp)) %>%
    cbind(cmdscale(mds))
  
  p <- ggplot(mds, aes(mds[, "1"], mds["2"])) +
    geom_point(aes(shape = mds[, mdsfact]), size = 2) +
    xlab("MDS coordinate 1") +
    ylab("MDS coordinate 2") +
    ggtitle("MDS Analysis") +
    theme_light() +
    scale_shape_discrete(name = mdsfact)
  print(p)
}



## HEAT

### Heatmap (1) array
qcHeatMap <- function(heat, n) {
  n <- as.numeric(n)
  pheatmap(
    mat = heat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste0("Heatmap of the ", n, " most expressed IDs")
  )
}


### Heatmap (2) counts
qcHeatCount <- function(s, rc.data, coldata, heatfactor) {
  heat.c <- getGenes(
    rc.data = rc.data, 
    id = s[["y"]],
    coldata = coldata
  )
  p <- ggplot(heat.c, aes(heat.c[, heatfactor], counts)) +
    geom_point(
      position = position_jitter(height = 0, width = 0.05),
      size = 2,
      color = "#3d3d3d"
    ) +
    stat_summary(
      fun.data = mean_cl_boot, 
      geom = "errorbar",
      width = 0.03, 
      colour = "#ff5454", 
      alpha = 0.7
    ) +
    stat_summary(
      fun.y = mean, 
      geom = "point", 
      fill = "#ff5454", 
      pch = 21, 
      size = 3
    ) +
    xlab("") +
    ylab("Normalized counts") +
    ggtitle(paste(s[["y"]], "Counts")) +
    theme_light()
  print(p)
}



## COR

### Correlation matrix
corMatPlot <- function(cor.mat) {
  cor.mat <- as.matrix(cor.mat)
  pheatmap(
    mat = cor.mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Correlation analysis"
  )
}

### Label converter
corLabConvert <- function(lab, tran, x, y) {
  if (tran == "log") {
    xlab <- bquote(.(x)*": log"[2]*"(counts + 1)")
    ylab <- bquote(.(y)*": log"[2]*"(counts + 1)")
    return(list(xlab, ylab))
  } else if (tran == "rlog") {
    xlab <- paste0(x, ": rlog(counts)")
    ylab <- paste0(y, ": rlog(counts)")
    return(list(xlab, ylab))
  } else if (tran == "vst") {
    xlab <- paste0(x, ": vst(counts)")
    ylab <- paste0(y, ": vst(counts")
    return(list(xlab, ylab))
  } else if (tran == "raw") {
    xlab <- paste0(x, ": raw counts")
    ylab <- paste0(y, ": raw counts")
    return(list(xlab, ylab))
  }
}


### Scatterplot
corScatter <- function(s.cor, cts.tran, tran, lab) {
  cts.tran <- as.data.frame(cts.tran)
  x <- s.cor[["x"]]
  y <- s.cor[["y"]]
  
  lab <- corLabConvert(lab, tran, x, y)
  xlab <- lab[[1]]
  ylab <- lab[[2]]
  
  p <- ggplot(cts.tran, aes(cts.tran[, x], cts.tran[, y])) +
    geom_point(size = 0.8) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(paste(y, "vs.", x)) +
    theme_light()
  print(p)
}



## DGE

### Overview plot
dgeOverPlot <- function(comp) {
  p <- ggplot(comp, aes(contrast, value)) +
    geom_bar(
      stat = "identity", 
      aes(fill = variable), 
      color = "#3d3d3d",
      position = position_dodge()
    ) +
    xlab("Comparison") +
    ylab("Number of IDs") +
    ggtitle("DGE Comparisons Overview") +
    theme_light() +
    scale_fill_manual(
      name = "Expression",
      values = c("down" = "#3884ff", "up" = "#ff4949")
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
}


### Get DGE Overview table 
dgeOverTbl <- function(cont.ls, lf, p) {
  p <- as.numeric(p)
  lf <- as.numeric(lf)
  for (i in 1:length(cont.ls)) {
    cont.ls[[i]] <- cont.ls[[i]][abs(cont.ls[[i]]$log2FoldChange) >= lf, ]
    cont.ls[[i]] <- cont.ls[[i]][cont.ls[[i]]$padj <= p, ]
  }
  cont.regup <- lapply(cont.ls, subset, log2FoldChange > 0)
  cont.regdn <- lapply(cont.ls, subset, log2FoldChange < 0)
  cont.regup <- lapply(cont.regup, nrow)
  cont.regup <- unlist(cont.regup)
  cont.regdn <- lapply(cont.regdn, nrow)
  cont.regdn <- unlist(cont.regdn)
  up.df <- data.frame(cont.regup)
  dn.df <- data.frame(cont.regdn)
  al.df <- merge(up.df, dn.df, by = 0, all = TRUE)
  colnames(al.df) <- c("contrast", "up", "down")
  al.df <- melt(al.df)
  al.df$contrast <- factor(al.df$contrast)
  al.df$variable <- factor(al.df$variable)
  return(al.df)
}


### Get MA plot
dgeMAPlot <- function(dgeout2, p, l, cont) {
  dge <- dgeout2
  dge$isPADJ <- ifelse(dge$padj <= p, TRUE, FALSE)
  dge$isPADJ[is.na(dge$isPADJ)] <- FALSE
  dge$isLFC <- ifelse(abs(dge$log2FoldChange) >= l, TRUE, FALSE)
  dge$isDGE <- "No differential expression"
  dge$isDGE[which(dge$isLFC & dge$isPADJ) ] <- paste(
    "padj <=", p, "& LFC >=", l
  )

  p <- ggplot(dge, aes(log10(baseMean), log2FoldChange)) +
    geom_point(size = 0.7, aes(color = isDGE)) +
    xlab(bquote("log"[10]*"(base mean)")) +
    ylab(bquote("log"[2]*"(fold change)")) +
    ggtitle(paste("MA Plot:", cont)) +
    theme_light() +
    scale_color_manual(
      name = "",
      values = c("#999999", "#3884ff")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    geom_hline(
      yintercept = 0, 
      color = "#3d3d3d", 
      linetype = "dashed"
    )
  print(p)
}


### Volcano plot
dgeVolPlot <- function(dgeout2, p, l, cont) {
  dge <- dgeout2
  dge$isPADJ <- ifelse(dge$padj <= p, TRUE, FALSE)
  dge$isPADJ[is.na(dge$isPADJ)] <- FALSE
  dge$isLFC <- ifelse(abs(dge$log2FoldChange) >= l, TRUE, FALSE)
  dge$isDGE <- "No differential expression"
  dge$isDGE[which(dge$isLFC & dge$isPADJ) ] <- paste(
    "padj <=", p, "& LFC >=", l
  )
  
  p <- ggplot(dge, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(size = 0.6, aes(color = isDGE)) +
    xlab(bquote("log"[2]*"(fold change)")) +
    ylab(bquote("-log"[10]*"(p-value)")) +
    ggtitle(paste("Volcano Plot:", cont)) +
    theme_light() +
    scale_color_manual(
      name = "",
      values = c("#999999", "#3884ff")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    geom_vline(
      xintercept = 0, 
      color = "#3d3d3d", 
      linetype = "dashed"
    )
  print(p)
}



## CTS Overview
trubble <- function(cts) {
  tmp <- as.data.frame(
    rbind(
      cts[1:5, ],
      ... = rep("...", length(cts[1, ])),
      cts[(nrow(cts) - 4):(nrow(cts)), ]
    )
  )
  
  if (ncol(tmp) > 10) {
    tmp2 <- tmp[, 1:10]
  } else {
    tmp2 <- tmp
  }
  
  nr <- nrow(cts)
  nc <- ncol(cts)
  
  if (ncol(tmp) > 10) {
    output <- paste(
      "Your pre-processed data contains", nr, "IDs and", nc, "samples. Showing the first 10 samples:\n"
    )
  } else {
    output <- paste(
      "Your pre-processed data contains", nr, "IDs and", nc, "samples.\n"
    )
  }
  
  test <- as.matrix(tmp2)
  test <- rbind(colnames(tmp2), test)
  # test <- cbind(c("", rownames(tmp2)), test)
  y <- sprintf(paste0("%",max(nchar(test)),"s"), test)
  y <- matrix(y, nrow = 12)
  
  gen <- c("", rownames(tmp2))
  gen <- gsub("\\s", " ", format(gen, width = max(nchar(gen))))
  
  if (ncol(tmp) > 10) {
    output2 <- paste("\nSamples not shown:\n\n")
    output3 <- paste(colnames(tmp[11:ncol(tmp)]))
  } else {
    output2 <- NULL
    output3 <- NULL
  }
  
  
  cat(output, "\n")
  for(i in 1:nrow(y)) {
    cat(gen[i], y[i, ], "\n")
  }
  
  cat(output2)
  for (i in 1:length(output3)) {
    cat("  ", output3[i], "\n")
  }
  
}
