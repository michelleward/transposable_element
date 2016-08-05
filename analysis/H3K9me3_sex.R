library(gplots)
library(RColorBrewer)
library(DESeq2)
library(GenomicRanges)

# read table, filter out chrX, chrY, random etc only keep autosomes.

autosome<-paste("chr",c(seq(1:22), "2A","2B"),sep="")

countdata.raw<-read.table("data.txt", header=F)
header<-read.table("header.txt",header=F)

metapeaks.hg19.raw<-read.table("comm.hg19.bed.ortho.hg19",header=F)
metapeaks.panTro3.raw<-read.table("comm.hg19.bed.ortho.panTro3",header=F)

metapeaks.hg19<-metapeaks.hg19.raw[metapeaks.hg19.raw[,1] %in% autosome & metapeaks.panTro3.raw[,1] %in% autosome,]
metapeaks.panTro3<-metapeaks.panTro3.raw[metapeaks.hg19.raw[,1] %in% autosome & metapeaks.panTro3.raw[,1] %in% autosome,]

countdata<-countdata.raw[metapeaks.hg19.raw[,1] %in% autosome & metapeaks.panTro3.raw[,1] %in% autosome,]

## DESeq2

# interspecies differences

# design input data
colnames(countdata)<-as.matrix((header[1,]))
countData<-countdata[,header[2,]=="H3K9me3"]

colData<-data.frame(t(header[3,header[2,]=="H3K9me3"]))
colnames(colData)<-c("species")
rownames(colData)<-colnames(countData)

# create DDS object

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ species)

# filter out rows with 0

dds100 <- dds [apply(counts(dds)>0,1,sum) > 10, ]

# run DEseq2
dds100 <- DESeq(dds100)
res <- results(dds100)

dim(res)
sum(res$padj < 0.01, na.rm=TRUE)

# volcano plot

plot(res$log2FoldChange[res$padj > .01], -log10(res$pvalue[res$padj > .01]), xlab= "log fold change", ylab="-log10 p-val", main="HvC \n 19632/166830 0.12  DE (FDR=0.01)", xlim=c(-10,10), ylim=c(0,50))
points(res$log2FoldChange[res$padj < .01], -log10(res$pvalue[res$padj < .01]), col = "red")
legend(6,50,  legend=c("adj p-val > 0.01", "adj p-val < 0.01"), col=c("black", "red"), cex=0.75, bty = "n", pch =16)

# human sex differences

colData<-data.frame(t(header[3,header[2,]=="H3K9me3"]))
colnames(colData)<-c("species")
rownames(colData)<-colnames(countData)

sex <- c("F","F","M","M","M","F","F","M","F","M","F","M","M","F","F","F","M")

colData2<-data.frame(t(header[3,header[2,]=="H3K9me3"]))
colData2$sex <- sex
colnames(colData2)<-c("species", "sex")
rownames(colData2)<-colnames(countData)

colData_human <- colData2[c(1,3,6,7,8,9,10,12,16,17), 2]
countData_human <-countData[,c(1,3,6,7,8,9,10,12,16,17)]
colData_human <- as.matrix(colData_human)
rownames(colData_human)<-colnames(countData_human)
colnames(colData_human) <- c("sex")

dds_human <- DESeqDataSetFromMatrix(countData = countData_human, colData = colData_human, design = ~ sex)
dds_human100 <- dds_human[apply(counts(dds_human)>0,1,sum) > 5, ]
dds_human100$sex <- relevel(dds_human100$sex, ref="F")
dds_human100 <- DESeq(dds_human100)
res_human <- results(dds_human100)
sum(res_human$padj < 0.01, na.rm=TRUE)

plot(res_human$log2FoldChange[res_human$padj > .01], -log10(res_human$pvalue[res_human$padj > .01]), xlab= "log fold change", ylab="-log10 p-val", main="Human M/F (4v4) \n 10382/168363 0.07  DE (FDR=0.01)", xlim=c(-10,10), ylim=c(0,50))
points(res_human$log2FoldChange[res_human$padj < .01], -log10(res_human$pvalue[res_human$padj < .01]), col = "red")
legend(6,50,  legend=c("adj p-val > 0.01", "adj p-val < 0.01"), col=c("black", "red"), cex=0.75, bty = "n", pch =16)

# chimp sex differences

colData_chimp <- colData2[c(2,4,5,11,13,14,15), 2]
countData_chimp <-countData[,c(2,4,5,11,13,14,15)]
colData_chimp <- as.matrix(colData_chimp)
rownames(colData_chimp)<-colnames(countData_chimp)
colnames(colData_chimp) <- c("sex")

dds_chimp <- DESeqDataSetFromMatrix(countData = countData_chimp, colData = colData_chimp, design = ~ sex)
dds_chimp100 <- dds_chimp[apply(counts(dds_chimp)>0,1,sum) > 4, ]
dds_chimp100$sex <- relevel(dds_chimp100$sex, ref="F")
dds_chimp100 <- DESeq(dds_chimp100)
res_chimp <- results(dds_chimp100)
sum(res_chimp$padj < 0.01, na.rm=TRUE)

par(mfrow=c(1,1))
plot(res_chimp$log2FoldChange[res_chimp$padj > .01], -log10(res_chimp$pvalue[res_chimp$padj > .01]), xlab= "log fold change", ylab="-log10 p-val", main="Chimp M/F (3v4) \n 6204/153555 0.04  DE (FDR=0.01)", xlim=c(-10,10), ylim=c(0,50))
points(res_chimp$log2FoldChange[res_chimp$padj < .01], -log10(res_chimp$pvalue[res_chimp$padj < .01]), col = "red")
legend(6,50,  legend=c("adj p-val > 0.01", "adj p-val < 0.01"), col=c("black", "red"), cex=0.75, bty = "n", pch =16)

# overlap sig MF human with sig MF chimp

res_human <- res_human[rowSums(is.na(res_human)) ==0, ]
sig_res_human <- res_human[res_human$padj < 0.01, ]

res_chimp <- res_chimp[rowSums(is.na(res_chimp)) ==0, ]
sig_res_chimp <- res_chimp[res_chimp$padj < 0.01, ]

sig_res_chimpdf <- as.data.frame(sig_res_chimp)
sig_res_humandf <- as.data.frame(sig_res_human)
sig_res_merge_all <-merge(sig_res_chimpdf, sig_res_humandf, by=0, all=TRUE)
sig_res_merge <-merge(sig_res_chimpdf, sig_res_humandf, by=0, all=FALSE)

1- phyper(5514, 6024, 153555, 12262, lower.tail=TRUE)

# plot M/F fold changes in human vs chimp

lm(sig_res_merge$log2FoldChange.x ~ sig_res_merge$log2FoldChange.y)
plot(sig_res_merge[, c(3,9)], main= " Sig. MF Human vs Sig. MF Chimp", xlab="Chimp M/F log2FC ",ylab="Human M/F log2FC")
legend(3.5,-6, legend="R2=0.51", cex=0.75, bty = "n" )
abline(0,1, col="blue")

# Venn diagram

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1= 6204, area2 = 12262, cross.area = 5514,  category   = c("Chimp", "Human"), fill=c("dark gray", "cornflower blue")  

# print out HvC data

colnames(metapeaks.hg19)<-c("hchr","hstart","hend")
colnames(metapeaks.panTro3)<-c("cchr","cstart","cend")

metapeakname<-rownames(res)

metapeaks.hg19.deseq2 <- metapeaks.hg19[metapeakname,]
metapeaks.panTro3.deseq2 <- metapeaks.panTro3[metapeakname,]

count.chimp100<-counts(dds100,normalized=T)[,c(2,4,5,11,13,14,15)]
count.human100<-counts(dds100,normalized=T)[,c(1,3,6,7,8,9,10,12,16,17)]

deseq2all_HvC_0.01<-cbind(metapeaks.hg19.deseq2,metapeaks.panTro3.deseq2,count.human100,count.chimp100,res)

padjcutoff_0.01<-0.01
deseq2all_HvC_0.01$tag<-ifelse(deseq2all_HvC_0.01$log2FoldChange > 0 & deseq2all_HvC_0.01$padj < padjcutoff_0.01 ,"human", "shared")
deseq2all_HvC_0.01$tag[which(deseq2all_HvC_0.01$log2FoldChange < 0 & deseq2all_HvC_0.01$padj < padjcutoff_0.01)] <- "chimp"

write.table(deseq2all_HvC_0.01,file='DEseq2_alldata_HvC_0.01.txt'
            , quote=F,row.names=F, col.names=T, sep="\t")
write.table(deseq2all_HvC_0.01[,c(1,2,3,7:30)],file='DEseq2_all_HvC_0.01.hg19'
            , quote=F,row.names=F, col.names=T, sep="\t")
write.table(deseq2all_HvC_0.01[,c(4,5,6,7:30)],file='DEseq2_all_HvC_0.01.panTro3'
            , quote=F,row.names=F, col.names=T, sep="\t")


# HvC pie chart

sum(grepl("shared",deseq2all_HvC_0.01$tag))
sum(grepl("human",deseq2all_HvC_0.01$tag))
sum(grepl("chimp",deseq2all_HvC_0.01$tag))
pievalues <- c(147198, 8972, 10660)
pienamevalues <- c("shared (88%)", "human-enriched (5%)", "chimp-enriched (6%)")
piecolours <- c("gray", "cornflower blue", "black")
pie(pievalues, labels = pienamevalues, col = piecolours, main="H3K9me3 divergence \n 166830 (FDR=0.01)")

# print out all M vs F data  

# all human MvF
head(res_human)

# all peaks in human
metapeakname_human <-rownames(res_human)
metapeaks.hg19.deseq2_human <- metapeaks.hg19[metapeakname_human,]
metapeaks.panTro3.deseq2_human <- metapeaks.panTro3[metapeakname_human,]

# MvF sig DE in both human & chimp
head(sig_res_merge)

# combine peak regions with all human MvF 
count.human100_human<-counts(dds_human100,normalized=T)
# find counts only in res_human regions? count.human100_human<-count.human100_human[res_human[,1] %in% autosome,]

deseq2all_human <-cbind(metapeaks.hg19.deseq2_human, metapeaks.panTro3.deseq2_human, res_human)

# all MvF differences in human that are conserved in chimp
# if sig_res_merge not in deseq2all_human tag with "shared"
# if sig_res_merge in deseq2all_human & deseq2all_human$log2FoldChange < 0, tag with "F"
# if sig_res_merge in deseq2all_human & deseq2all_human$log2FoldChange > 0, tag with "M"




                                                  
