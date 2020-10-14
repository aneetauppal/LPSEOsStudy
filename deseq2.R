#Aneeta Uppal
#Deseq2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages("gplots")
install.packages("dplyr")
install.packages("broom")
install.packages("ggpubr")

#load libraries
library(dplyr)
library(tools)
library(DESeq2)
library("RColorBrewer")
library("gplots")
library("ggpubr")
library("data.table")
#library(broom)


files=list.files(path="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/counted2",pattern="\\d*_Count2rnaseqD1.txt$",full.names=T)
files
metadata <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/datainfo.csv")
metadata2 <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/datainfounordered2.csv")

dl=lapply(files, read.delim2, header=FALSE)
# give names to each data frame within the list

names(dl)=tools::file_path_sans_ext(basename(files))
# Generate a data frame by mering all the data frame within the list
dl.df=do.call(cbind, dl)
dl

final.count=cbind(dl.df[,1],select(dl.df,contains("V2")))
# Provide column (sample) names for each column
colnames(final.count)=c("GeneSymbol",basename(files))
final.count

# Save the table to hard disk
write.table(final.count, "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/htseq.r.results.tsv", sep="\t", row.names=FALSE)


# Extract file names  without the extension
basename(files)
fnames=file_path_sans_ext(basename(files))
filenames=(basename(files))


flnames=sub("htseq.","\\1", fnames)
flnames


# Extract conditions
cond=do.call(rbind,strsplit(flnames, '_'))[,1]
cond
conditions = metadata2[,2]


#????
#filnames = file_path_sans_ext(filenames)

# Prepare sample details as data frame
stable <- data.frame(sampleName = flnames, fileName = basename(files), condition =conditions)
stable

#example code
#dds$ALL <- factor(paste0(dds$time, dds$strain)) design(dds) <- ~ ALL
#dds <- DESeq(dds) resultsNames(dds) results(dds, contrast=c("ALL","12", "11"))


# Prepare a DEXseq object -played with this in dexseq
#htdex=DESeqDataSetFromHTSeqCount(
#  sampleTable = stable,
#  directory = "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted", 
#  design= ~ sample + exon + condition:exon)

#Prepared Deseq object
htdex=DESeqDataSetFromHTSeqCount(
  sampleTable = stable,
  directory = "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/counted2", 
  design= ~ condition)




###htdex$condition <- factor(htdex$condition, levels =c("LPSalone", "Boswelliasystemic", "Unstimulated", "Boswelliatopical", "Coconutoilsystemic","Coconutoiltpoical"))

htdex$condition <- relevel(htdex$condition, ref = "LPSalone")
htdex
colData(dxd)


# Calculate differential expression
dhtdex=DESeq(htdex)



# Extract and store the results
res.dhtdex=results(dhtdex)

res<-results(dhtdex)

resOrdered <- res[order(res$pvalue),]
summary(res)

resultsNames(dhtdex)

#pull out rows that only code for transcript mrna - aka leave ncrna's and other bs
newres <- res[rownames(res) %like% "rna-NM",]


results2_2 <- results(dhtdex, contrast = c("condition", "Boswelliasystemic", "LPSalone"  ) )

#results2_3 <- results(dhtdex, contrast = c("condition", "Boswelliasystemic", "Unstimulated"  ))

results3_2 <- results(dhtdex, contrast = c("condition", "Boswelliatopical", "LPSalone"  ) )
results4_2 <- results(dhtdex, contrast = c("condition", "Coconutoilsystemic", "LPSalone"  ) )
results5_2 <- results(dhtdex, contrast = c("condition", "Coconutoiltopical", "LPSalone"  ) )
results6_2 <- results(dhtdex, contrast = c("condition", "Unstimulated", "LPSalone"  ) )

#only potting transcript mrna's - aka rna-NM___ 
results2_cut <- results2_2[rownames(results2_2) %like% "rna-NM",]
results3_cut <- results3_2[rownames(results3_2) %like% "rna-NM",]
results4_cut <- results4_2[rownames(results4_2) %like% "rna-NM",]
results5_cut <- results5_2[rownames(results5_2) %like% "rna-NM",]
results6_cut <- results6_2[rownames(results6_2) %like% "rna-NM",]




#write.table(as.data.frame(resOrdered),file="/Users/aneetauppal/counted/treated_untreated_BosSys1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(resOrdered),file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/treated_untreated_unstim1.txt", col.names = TRUE, sep='\t')


write.table(as.data.frame(results2_2), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/bossys3_2.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results3_2), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/bostop_2.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results4_2), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/cosys_2.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results5_2), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/cotop_2.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results6_2), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/unstim_2.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results7), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/lpsalone.txt", col.names = TRUE, sep='\t')


write.table(as.data.frame(results2_3), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/bosysUnstim.txt", col.names = TRUE, sep='\t')


#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dhtdex, alpha=0.05)
res02 <- results(results2_3, lfcThreshold = 2, alpha=0.05)

summary(res05)

resultNames(res)

#MA-plot from base means and log fold changes ofdifferentially expressed genes in the frankincense systemic application"

plotMA(results2, ylim=c(-10,10), main = "MA-plot: Frankincense Systemic Application", alpha = 0.05)

#FYI FDR stands for padjust values FC = fold change NOT log2foldchange, but log base 2 of whatever is in the table

franksys <- ggmaplot(results2_cut, main = "MA-plot: Frankincense systemic application",
        fdr = 0.05, fc = 2.25, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(row.names(results2_cut)),
         legend = "top", top = 0,
         font.label = c("bold", 18),
         font.legend = c("bold", 18),
         font.main = c("bold", 20),
         ylim = c(-10,10),
         ggtheme = ggplot2::theme_grey())

franksys +  font("xy", size = 20)+ font("xy.text", size = 20)

franktop <- ggmaplot(results3_cut, main = "MA-plot: Frankincense topical application",
         fdr = 0.05, fc = 2.25, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(row.names(results3_cut)),
         legend = "top", top = 0,
         font.label = c("bold", 18),
         font.legend = c("bold", 18),
         font.main = c("bold", 20),
         ylim = c(-10,10),
         ggtheme = ggplot2::theme_grey())

franktop + font("xy", size = 20)+ font("xy.text", size = 20)

FCOtop <- ggmaplot(results5_cut, main = "MA-plot: FCO topical application",
         fdr = 0.05, fc = 2.25, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(row.names(results5_cut)),
         legend = "top", top = 0,
         font.label = c("bold", 18),
         font.legend = c("bold", 18),
         font.main = c("bold", 20),
         ylim = c(-10,10),
         ggtheme = ggplot2::theme_grey())

FCOtop + font("xy", size = 20)+ font("xy.text", size = 20)


Unstim <- ggmaplot(results6_cut, main = "MA-plot: Untreated",
         fdr = 0.05, fc = 2.25, size = 1.5,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(row.names(results6_cut)),
         legend = "top", top = 0,
         font.label = c("bold", 18),
         font.legend = c("bold", 18),
         font.main = c("bold", 20),
         ylim = c(-10,10),
         ggtheme = ggplot2::theme_grey())

Unstim + font("xy", size = 20)+ font("xy.text", size = 20)

FCOsys <- ggmaplot(results4_cut, main = "MA-plot: FCO Systemic",
                   fdr = 0.05, fc = 2.25, size = 1.3,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(row.names(results4_cut)),
                   legend = "top", top = 0,
                   font.label = c("bold", 12),
                   font.legend = c("bold", 18),
                   font.main = c("bold", 20),
                   ylim = c(-10,10),
                   ggtheme = ggplot2::theme_grey())
FCOsys + font("xy", size = 20)+ font("xy.text", size = 20)


idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# MA plot
png("ma.png")
plotMA(res.dhtdex, ylim=c(-2,2))
dev.off()
# Comparison between two samples for a given gene (i.e gene that shows highest fold change difference)
png("plotcounts.png")
plotCounts(dhtdex, gene=which.max(res.dhtdex$log2FoldChange), intgroup="condition")
dev.off()
# Log transform the data
rld.dhtdex=rlog(dhtdex)
# PCA analysis
png("pca.png")
plotPCA(rld.dhtdex)
dev.off()
# Prepare data for histogram
select <- order(rowMeans(counts(dhtdex,normalized=TRUE)),decreasing=TRUE)[1:30]
# Select the color from preexisting colors within R
hmcol <- colorRampPalette(brewer.pal(3, "BrBG"))(100)
# Draw heatmap for selected 30 genes
png("counts.png")
heatmap.2(counts(dhtdex,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.off()
png("rld.png")
heatmap.2(assay(rld.dhtdex)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
# Calculate the distance betwen the samples/groups
distsRL <- dist(t(assay(rld.dhtdex)))
# Store the distance as matrix
mat <- as.matrix(distsRL)
# supply row names and column names for this matrix
rownames(mat) <- colnames(mat)<-row.names(colData(dhtdex))
# Calculate clustering distances
hc <- hclust(distsRL)
# Plot the distances and dendrogram together
png("distances.png")
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
dev.off()
# Plot Dispersion estimates
png("disp.png")
plotDispEsts(dhtdex)
dev.off()
# Save the image for further processing
save.image("./htseqresults/htseq.results.Rdata")






#####################3
file <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Knowledgebase/publishedstudiesofEOs.csv")
file
par(mar=c(7,8,5,2)+.1)
plot(file[,1],file[,2],      
     col="black",
     pch=16,
     xlab="Year", ylab="Number of Published Studies", cex.lab=2.5, cex.main=2.5, cex.axis=2.0,
     main="The Increasing Number of Published Studies on EOs \nby Year (1833-2019)")


###################################
#Deseq rerun of the code after rerunning alignment/count files due to cluster losing all of the data
#within no back up april 16th 2020

files=list.files(path="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/counted2",pattern="\\d*_Count2rnaseqD1.txt$",full.names=T)
files
metadata <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/datainfo.csv")
metadata2 <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/datainfounordered.csv")

dl=lapply(files, read.delim2, header=FALSE)
# give names to each data frame within the list

names(dl)=tools::file_path_sans_ext(basename(files))
# Generate a data frame by mering all the data frame within the list
dl.df=do.call(cbind, dl)
dl

final.count=cbind(dl.df[,1],select(dl.df,contains("V2")))
# Provide column (sample) names for each column
colnames(final.count)=c("GeneSymbol",basename(files))
final.count

# Save the table to hard disk
write.table(final.count, "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/counted2/htseq_2.results.tsv", sep="\t", row.names=FALSE)


# Extract file names  without the extension
basename(files)
fnames=file_path_sans_ext(basename(files))
filenames=(basename(files))


flnames=sub("htseq.","\\1", fnames)
flnames
fileName

# Extract conditions
cond=do.call(rbind,strsplit(flnames, '_'))[,1]
cond


conditions = metadata2[,2]


#????
#filnames = file_path_sans_ext(filenames)

# Prepare sample details as data frame
stable <- data.frame(sampleName = flnames, fileName = basename(files),condition =conditions)
stable

#example code
#dds$ALL <- factor(paste0(dds$time, dds$strain)) design(dds) <- ~ ALL
#dds <- DESeq(dds) resultsNames(dds) results(dds, contrast=c("ALL","12", "11"))


# Prepare a DEXseq object -played with this in dexseq
#htdex=DESeqDataSetFromHTSeqCount(
#  sampleTable = stable,
#  directory = "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted", 
#  design= ~ sample + exon + condition:exon)

#Prepared Deseq object
htdex=DESeqDataSetFromHTSeqCount(
  sampleTable = stable,
  directory = "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/counted2", 
  design= ~ condition)


#Get variance stabilizing transformation for WGCNA application july 8
Wgcnaobj=varianceStabilizingTransformation(htdex, blind = FALSE, fitType="mean")


###htdex$condition <- factor(htdex$condition, levels =c("LPSalone", "Boswelliasystemic", "Unstimulated", "Boswelliatopical", "Coconutoilsystemic","Coconutoiltpoical"))

htdex$condition <- relevel(htdex$condition, ref = "LPSalone")
htdex
colData(dxd)


# Calculate differential expression
dhtdex=DESeq(htdex)


# Extract and store the results
res.dhtdex=results(dhtdex)

res<-results(dhtdex)

resOrdered <- res[order(res$pvalue),]
summary(res)

resultsNames(dhtdex)



results2 <- results(dhtdex, contrast = c("condition", "Boswelliasystemic", "LPSalone"  ) )

results2_3 <- results(dhtdex, contrast = c("condition", "Boswelliasystemic", "Unstimulated"  ))

results3 <- results(dhtdex, contrast = c("condition", "Boswelliatopical", "LPSalone"  ) )
results4 <- results(dhtdex, contrast = c("condition", "Coconutoilsystemic", "LPSalone"  ) )
results5 <- results(dhtdex, contrast = c("condition", "Coconutoiltopical", "LPSalone"  ) )
results6 <- results(dhtdex, contrast = c("condition", "Unstimulated", "LPSalone"  ) )
results7 <- results(dhtdex, contrast = c("condition", "LPSalone", "Unstimulated") )


#write.table(as.data.frame(resOrdered),file="/Users/aneetauppal/counted/treated_untreated_BosSys1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(resOrdered),file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/treated_untreated_unstim1.txt", col.names = TRUE, sep='\t')


write.table(as.data.frame(results2), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/bossys3_1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results3), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/bostop_1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results4), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/cosys_1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results5), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/cotop_1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results6), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/unstim_1.txt", col.names = TRUE, sep='\t')
write.table(as.data.frame(results7), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/lpsalone.txt", col.names = TRUE, sep='\t')


write.table(as.data.frame(results2_3), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/bosysUnstim.txt", col.names = TRUE, sep='\t')


#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dhtdex, alpha=0.05)
summary(res05)

resultNames(res)

#MA-plot from base means and log fold changes ofdifferentially expressed genes in the frankincense systemic application"

plotMA(results2, ylim=c(-10,10), main = "MA-plot: Frankincense Systemic Application", alpha = 0.05)

#FYI FDR stands for padjust values FC = fold change NOT log2foldchange, but log base 2 of whatever is in the table
##so log base 2 of 4 = 2 

franksys <- ggmaplot(results2, main = "MA-plot: Frankincense systemic application",
                     fdr = 0.05, fc = 4, size = 1.3,
                     palette = c("#B31B21", "#1465AC", "darkgray"),
                     genenames = as.vector(row.names(results2)),
                     legend = "top", top = 0,
                     font.label = c("bold", 12),
                     font.legend = c("bold", 16),
                     font.main = c("bold", 16),
                     ggtheme = ggplot2::theme_grey())

franksys +  font("xy", size = 18)+ font("xy.text", size = 16)

franktop <- ggmaplot(results3, main = "MA-plot: Frankincense topical application",
                     fdr = 0.05, fc = 3, size = 1.3,
                     palette = c("#B31B21", "#1465AC", "darkgray"),
                     genenames = as.vector(row.names(results3)),
                     legend = "top", top = 0,
                     font.label = c("bold", 12),
                     font.legend = c("bold", 16),
                     font.main = c("bold", 16),
                     ggtheme = ggplot2::theme_grey())

franktop + font("xy", size = 18)+ font("xy.text", size = 16)

FCOtop <- ggmaplot(results5, main = "MA-plot: FCO topical application",
                   fdr = 0.05, fc = 3, size = 1.3,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(row.names(results5)),
                   legend = "top", top = 0,
                   font.label = c("bold", 12),
                   font.legend = c("bold", 16),
                   font.main = c("bold", 16),
                   ggtheme = ggplot2::theme_grey())

FCOtop + font("xy", size = 18)+ font("xy.text", size = 16)

Unstim <- ggmaplot(results6, main = "MA-plot: Unstimulated",
                   fdr = 0.05, fc = 4, size = 1.3,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(row.names(results6)),
                   legend = "top", top = 0,
                   font.label = c("bold", 12),
                   font.legend = c("bold", 16),
                   font.main = c("bold", 16),
                   ggtheme = ggplot2::theme_grey())

Unstim + font("xy", size = 18)+ font("xy.text", size = 16)





idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]












#### playing with data for wgcna

#Prepared Deseq object
htdex=DESeqDataSetFromHTSeqCount(
  sampleTable = stable,
  directory = "/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/counted/counted2", 
  design= ~ condition)

htdex$condition <- relevel(htdex$condition, ref = "LPSalone")

write.table(as.data.frame(htdex), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/wgcna/wcgnaprep.txt", col.names = TRUE, sep='\t')

#Get variance stabilizing transformation for WGCNA application july 8
Wgcnaobj=varianceStabilizingTransformation(htdex, blind = FALSE, fitType="mean")



#for gsea
dds <- dhtdex
normcount <- vst(dds, blind=FALSE)
head(assay(normcount),3)
write.table(assay(normcount), file="/Users/aneetauppal/Graduate_PhD/DissertationResearch/RNAseq/normalized_counts.txt", col.names = TRUE, sep='\t')
