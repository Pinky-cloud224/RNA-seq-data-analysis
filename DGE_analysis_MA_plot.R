library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library(ggplot2)

samples <- read.table("~/Desktop/tif3_urea_final/tif3_urea.csv",sep=",",header=TRUE)
samples # To read the selected flat csv file, sep is the separator character,
# header is a logical value indicating whether the file contains the names of the variables its first line

# specification of the path to the files using the appropriate columns of samples (here, SRR_Number) and we read in a table that links the 
# and to read the quantification of abundances of transcripts
files <-file.path("~/Desktop/tif3_urea_final",samples$SRR_Number,"abundance.tsv")
files

txi <- tximport (files,type="kallisto",txOut=TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,colData=samples,design=~Treatment)
dds <- DESeq(ddsTxi)
# txOut is for whether the function should just output transcript level.

######### Differential expression analysis

dds <- DESeq(dds)
res <- results(dds,contrast=c("Treatment","wild_none","wild_3%-urea")) # the contrast is between wild type without any treatment and the mutant (tif3[delta]ntd) treated with 3% urea
res
head(res)

# The log2fold change is the log-ratio of a gene's expression values in two different conditions
# the positive values of log2 Fold change means the gene is more expressed in Treatment and the negative values means gene is more expressed in Control



###### p-values and adjusted p-values
# ordering the results table by the samllest p value
resOrdered <- res[order(res$padj),]
summary(res) #summarizing the basic tallies

### How many adjusted p-values were less than 0.05?
sum(res$padj <0.05,na.rm=TRUE)


#### Exploring and exporting results
####### generation of MA plot
library(ggpubr)

ggmaplot(res, main = expression("scd6-wild_none" %->% "scd6-wild_3%-urea"),
         fdr = 0.05, fc=1,size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res$name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         xlab = "Normalized mean expression",
         ylab = "Log2 fold change",
         ggtheme = ggplot2::theme_minimal())

#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1
# The MA plot visualizes the differences between measurements taken in two samples, by transforming data onto M (log ratio) and A (mean average) scales
# genes with simiar expression values in both normal and treated samples will cluster around M=0 values, points above and below from M=0 is upregulated and downregulated respectively


# Saving the log2Fold change data as a CSV file format
write.csv(res,file="~/Desktop/scd6/scd6_urea/wild_none_wild_3%-urea.csv")
