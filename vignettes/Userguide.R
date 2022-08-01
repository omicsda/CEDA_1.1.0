## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----data---------------------------------------------------------------------
library(CEDA)
data("mda231")
class(mda231)
length(mda231)
names(mda231)

## ----sgRNA--------------------------------------------------------------------
dim(mda231$sgRNA)
length(mda231$neGene$Gene)
head(mda231$sgRNA)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  UN <- apply((mda231$sgRNA[, 6:8]), 1, mean)
#  TR <- apply(log2(mda231$sgRNA[, 3:8]), 1, mean)
#  foo <- TR/UN

## ----neGene-------------------------------------------------------------------
dim(mda231$neGene)
head(mda231$neGene)

## ----normalization------------------------------------------------------------
mda231.ne <- mda231$sgRNA[mda231$sgRNA$Gene %in% mda231$neGene$Gene,]
cols <- c(3:8)
mda231.norm <- medianNormalization(mda231$sgRNA[,cols], mda231.ne[,cols])[[2]]

## ----design-------------------------------------------------------------------
group <- gl(2, 3, labels=c("Control","Baseline"))
design <- model.matrix(~  0 + group)
colnames(design) <- sapply(colnames(design), function(x) substr(x, 6, nchar(x)))
contrast.matrix <- makeContrasts("Control-Baseline", levels=design)

## ----limfit-------------------------------------------------------------------
limma.fit <- runLimma(log2(mda231.norm+1),design,contrast.matrix)

## ----merge--------------------------------------------------------------------
mda231.limma <- data.frame(mda231$sgRNA, limma.fit)
head(mda231.limma)

## ----betanull-----------------------------------------------------------------
betanull <- permuteLimma(log2(mda231.norm + 1), design, contrast.matrix, 20)
theta0 <- sd(betanull)
theta0

## ----mm, results='hide'-------------------------------------------------------
nmm.fit <- normalMM(mda231.limma, theta0)

## ----fig1, fig.cap = "Log fold ratios of sgRNAs vs. gene expression level"----
scatterPlot(nmm.fit$data,fdr=0.05,xlim=c(-0.5,12),ylim=c(-8,5))

## ----pval---------------------------------------------------------------------
mda231.nmm <- nmm.fit[[1]]
p.gene <- calculateGenePval(exp(mda231.nmm$log_p), mda231.nmm$Gene, 0.05)
fdr.gene <- stats::p.adjust(p.gene$pvalue, method = "fdr")
lfc.gene <- calculateGeneLFC(mda231.nmm$lfc, mda231.nmm$Gene)

