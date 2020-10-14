#Aneeta Uppal
#DESeq2Code #adapted from DESeq2 tutorial

#load libraries
library(dplyr)
library(tools)
library(DESeq2)
library("RColorBrewer")
library("gplots")
library("ggpubr")
library("data.table")


#list path to where files of the count data are based on their naming format 
files=list.files(path="path",pattern="\\d*_Count2rnaseqD1.txt$",full.names=T)

#where the metadata/trait data is stored
metadata <- read.csv("path")

dl=lapply(files, read.delim2, header=FALSE)

# give names to each data frame within the list
names(dl)=tools::file_path_sans_ext(basename(files))

# Generate a data frame by mering all the data frame within the list
dl.df=do.call(cbind, dl)
final.count=cbind(dl.df[,1],select(dl.df,contains("V2")))

# Provide column (sample) names for each column
colnames(final.count)=c("GeneSymbol",basename(files))
final.count

# Save the table to hard disk
write.table(final.count, "pathfile", sep="\t", row.names=FALSE)


# Extract file names  without the extension
basename(files)
fnames=file_path_sans_ext(basename(files))
filenames=(basename(files))

flnames=sub("htseq.","\\1", fnames)
flnames


# Extract conditions
conditions = metadata[,2]

# Prepare sample details as data frame
stable <- data.frame(sampleName = flnames, fileName = basename(files), condition =conditions)
stable

#Prepared Deseq object
htdex=DESeqDataSetFromHTSeqCount(
  sampleTable = stable,
  directory = "/path/", 
  design= ~ condition)

#relevel with the reference-  positive control 
htdex$condition <- relevel(htdex$condition, ref = "LPSalone")

# Calculate differential expression
dhtdex=DESeq(htdex)

# Extract and store the results
res.dhtdex=results(dhtdex)
res<-results(dhtdex)

resOrdered <- res[order(res$pvalue),]
summary(res)

resultsNames(dhtdex)


#save results for each comparison analysis
results2 <- results(dhtdex, contrast = c("condition", "Boswelliasystemic", "LPSalone"  ) )
results3 <- results(dhtdex, contrast = c("condition", "Boswelliatopical", "LPSalone"  ) )
results4 <- results(dhtdex, contrast = c("condition", "Coconutoilsystemic", "LPSalone"  ) )
results5 <- results(dhtdex, contrast = c("condition", "Coconutoiltopical", "LPSalone"  ) )
results6 <- results(dhtdex, contrast = c("condition", "Unstimulated", "LPSalone"  ) )

#save results to file
write.table(as.data.frame(results2), file="bossys.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results3), file="bostop.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results4), file="cosys.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results5), file="cotop.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results6), file="unstim.txt", col.names = TRUE, sep='\t')

#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)


#FDR stands for padjust values FC = fold change NOT log2foldchange, but log base 2 of whatever is in the table
#graph MA plots using ggmaplot - i preferred these graphs of DEGs
franksys <- ggmaplot(results2, main = "MA-plot: Frankincense systemic application",
                     fdr = 0.05, fc = 2.25, size = 1.5,
                     palette = c("#B31B21", "#1465AC", "darkgray"),
                     genenames = as.vector(row.names(results2)),
                     legend = "top", top = 0,
                     font.label = c("bold", 18),
                     font.legend = c("bold", 18),
                     font.main = c("bold", 20),
                     ylim = c(-10,10),
                     ggtheme = ggplot2::theme_grey())

franksys +  font("xy", size = 20)+ font("xy.text", size = 20)

franktop <- ggmaplot(results3, main = "MA-plot: Frankincense topical application",
                     fdr = 0.05, fc = 2.25, size = 1.5,
                     palette = c("#B31B21", "#1465AC", "darkgray"),
                     genenames = as.vector(row.names(results3)),
                     legend = "top", top = 0,
                     font.label = c("bold", 18),
                     font.legend = c("bold", 18),
                     font.main = c("bold", 20),
                     ylim = c(-10,10),
                     ggtheme = ggplot2::theme_grey())

franktop + font("xy", size = 20)+ font("xy.text", size = 20)

FCOtop <- ggmaplot(results5, main = "MA-plot: FCO topical application",
                   fdr = 0.05, fc = 2.25, size = 1.5,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(row.names(results5)),
                   legend = "top", top = 0,
                   font.label = c("bold", 18),
                   font.legend = c("bold", 18),
                   font.main = c("bold", 20),
                   ylim = c(-10,10),
                   ggtheme = ggplot2::theme_grey())

FCOtop + font("xy", size = 20)+ font("xy.text", size = 20)


Unstim <- ggmaplot(results6, main = "MA-plot: Untreated",
                   fdr = 0.05, fc = 2.25, size = 1.5,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(row.names(results6)),
                   legend = "top", top = 0,
                   font.label = c("bold", 18),
                   font.legend = c("bold", 18),
                   font.main = c("bold", 20),
                   ylim = c(-10,10),
                   ggtheme = ggplot2::theme_grey())

Unstim + font("xy", size = 20)+ font("xy.text", size = 20)

FCOsys <- ggmaplot(results4, main = "MA-plot: FCO Systemic",
                   fdr = 0.05, fc = 2.25, size = 1.3,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(row.names(results4)),
                   legend = "top", top = 0,
                   font.label = c("bold", 12),
                   font.legend = c("bold", 18),
                   font.main = c("bold", 20),
                   ylim = c(-10,10),
                   ggtheme = ggplot2::theme_grey())
FCOsys + font("xy", size = 20)+ font("xy.text", size = 20)


idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

