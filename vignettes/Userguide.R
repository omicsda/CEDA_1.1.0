## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----data---------------------------------------------------------------------
library(CEDA)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggprism)
set.seed(1)
data("mda231")
class(mda231)
length(mda231)
names(mda231)

## ----sgRNA--------------------------------------------------------------------
dim(mda231$sgRNA)
length(mda231$neGene$Gene)
head(mda231$sgRNA)

## ----neGene-------------------------------------------------------------------
dim(mda231$neGene)
head(mda231$neGene)

## ----filter-------------------------------------------------------------------
count_df <- mda231$sgRNA
mx.count = apply(count_df[,c(3:8)],1,function(x) sum(x>=1))
table(mx.count)
# keep sgRNA with none zero count in at least one sample
count_df2 = count_df[mx.count>=1,]

## ----normalization------------------------------------------------------------
mda231.ne <- count_df2[count_df2$Gene %in% mda231$neGene$Gene,]
cols <- c(3:8)
mda231.norm <- medianNormalization(count_df2[,cols], mda231.ne[,cols])[[2]]

## ----design-------------------------------------------------------------------
group <- gl(2, 3, labels=c("Control","Baseline"))
design <- model.matrix(~  0 + group)
colnames(design) <- sapply(colnames(design), function(x) substr(x, 6, nchar(x)))
contrast.matrix <- makeContrasts("Control-Baseline", levels=design)

## ----limfit-------------------------------------------------------------------
limma.fit <- runLimma(log2(mda231.norm+1),design,contrast.matrix)

## ----merge--------------------------------------------------------------------
mda231.limma <- data.frame(count_df2, limma.fit)
head(mda231.limma)

## ----betanull-----------------------------------------------------------------
betanull <- permuteLimma(log2(mda231.norm + 1), design, contrast.matrix, 20)
theta0 <- sd(betanull)
theta0

## ----mm, results='hide'-------------------------------------------------------
nmm.fit <- normalMM(mda231.limma, theta0, n.b=5)

## ----fig1, fig.cap = "Log fold ratios of sgRNAs vs. gene expression level"----
scatterPlot(nmm.fit$data,fdr=0.05,xlim=c(-0.5,12),ylim=c(-8,5))

## ----pval---------------------------------------------------------------------
mda231.nmm <- nmm.fit[[1]]
p.gene <- calculateGenePval(exp(mda231.nmm$log_p), mda231.nmm$Gene, 0.05)
gene_fdr <- stats::p.adjust(p.gene$pvalue, method = "fdr")
gene_lfc <- calculateGeneLFC(mda231.nmm$lfc, mda231.nmm$Gene)

gene_summary <- data.frame(gene_pval=unlist(p.gene$pvalue), gene_fdr=as.matrix(gene_fdr), gene_lfc = as.matrix(gene_lfc))
gene_summary$gene <- rownames(gene_summary)
gene_summary <- gene_summary[,c(4,1:3)]

## ----summary------------------------------------------------------------------
#extract gene expression data
gene.express <- mda231.nmm %>% group_by(Gene) %>% summarise_at(vars(exp.level.log2), max)
#merge gene summary with gene expression
gdata <- left_join(gene_summary, gene.express, by = c("gene" = "Gene"))
gdata <- gdata %>% filter(is.na(exp.level.log2)==FALSE)
# density plot and ridge plot
gdata$gene.fdr <- gdata$gene_fdr
data <- preparePlotData(gdata, gdata$gene.fdr)

## ----fig2, fig.cap = "2D density plot of gene log fold ratios vs. gene expression level for different FDR groups"----
densityPlot(data)

## ----fig3, fig.cap = "Ridge plot of gene expression for different FDR groups"----
ridgePlot(data)

