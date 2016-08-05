## RNAseq data analysis

library(edgeR)
library(gplots)
library(Mfuzz)
library(sva)
library(RColorBrewer)
library(limma)
library(statmod)

# Read in data and plot heatmap

RNAcountdata.raw<-read.table("RNA_counts.txt", header=F)
RNAheader<-read.table("RNA_header.txt",header=F)
RNAgenes<-read.table("RNA_rownames.txt", header=F)
RNAautogenes_H<-read.table("autoorthogenes_H.txt", header=F)
RNAautogenes_C<-read.table("autoorthogenes_C.txt", header=F)
RNAgenes_H<-read.table("orthogenes_H.txt", header=F)
RNAgenes_C<-read.table("orthogenes_C.txt", header=F)

# Remove sex chr

autosome<-paste("chr",c(seq(1:22), "2A","2B"),sep="")

RNAautogenes<-RNAautogenes_H[RNAautogenes_H[,2] %in% autosome & RNAautogenes_C[,2] %in% autosome,]
RNAcounts <- RNAcountdata.raw[RNAautogenes_H[,2] %in% autosome & RNAautogenes_C[,2] %in% autosome,]

# format data

colnames(RNAcounts)<-as.matrix((RNAheader[1,]))
rownames(RNAcounts)<-as.matrix((RNAautogenes[,1]))

colnames(RNAcounts)<-as.matrix((RNAheader[1,]))
rownames(RNAcounts)<-as.matrix((RNAautogenes[,1]))

RNAcounts <- RNAcounts[rowSums(is.na(RNAcounts)) == 0,]

# correlation heatmap & dendrogram

RNAcortable<-cor(as.matrix(RNAcounts),method="spearman")

par(mfrow=c(1,1))
heatmap.2(as.matrix(RNAcortable), key=T, trace="none", col=brewer.pal(9, "Greens"), main="Counts at auto ortho genes", margins = c(11, 11))

plot(hclust(dist(t(RNAcortable))), main = "Read counts at orthogenes")

# Filter data 10 individ>0 RNAcounts100 or 6H+4C>0 RNAcounts640

RNAcounts100 = RNAcounts[apply(RNAcounts>0,1,sum) > 10, ]

RNAcounts640 <-RNAcounts[rowSums(RNAcounts_human > 0) >= 6 | rowSums(RNAcounts_chimp > 0) >= 4,]

RNAcounts_human <- RNAcounts[, c(1,2,3,4,12,13,14,15,16,17)]
RNAcounts_chimp <-RNAcounts[, c(5,6,7,8,9,10,11)]

# Convert to CPM

cpms=cpm(RNAcounts, log =T)
cpms100=cpm(RNAcounts100, log = T)
cpms640=cpm(RNAcounts640, log = T)

cpmsordered <-cpms[,c(5,6,7,8,9,10,11,1,2,3,4,12,13,14,15,16,17)]

# plot cpm values

Leg=c(5,5,5,5,1,1,1,1,1,1,1,4,4,4,4,4,4)
boxplot(cpms[,order(Leg)], las =2, main = "log2 CPM")
boxplot(cpms100[,order(Leg)], las =2, main = "log2 CPM filter 10>0")
boxplot(cpms640[,order(Leg)], las =2, main = "log2 CPM filter 6H>0 + 4C>0")

# Normalise data

cpms.quantile = normalizeQuantiles(cpm(RNAcounts.more.100, log =T), ties=TRUE)
Leg = c("Light Sky Blue", "Light Sky Blue", "Light Sky Blue", "Light Sky Blue", 1,1,1,1,1,1,1,"Royal BLue",
        "Royal Blue", "Royal Blue",  "Royal Blue", "Royal Blue",  "Royal Blue", "Royal Blue")
boxplot(cpms.quantile[,order(Leg)], las =2, main = "Norm CPMs") 

# plot cpm correlations

RNAcpmcortable<-cor(as.matrix(cpms),method="spearman")
RNAcpm100cortable<-cor(as.matrix(cpms100),method="spearman")
RNAcpm640cortable<-cor(as.matrix(cpms640),method="spearman")

par(cex.main=1)
heatmap.2(as.matrix(RNAcpmcortable), key=T, trace = "none", col=brewer.pal(9, "Greens"), main=" Spearman correlation log2cpm \n ortho genes", margins = c(11, 11))
heatmap.2(as.matrix(RNAcpm100cortable), key=T, trace = "none", col=brewer.pal(9, "Greens"), main=" Spearman correlation log2cpm \n filtered 10>0 ortho genes", margins = c(11, 11))

plot(hclust(dist(t(cpms))), main = "log2cpm read counts \n orthogenes")
plot(hclust(dist(t(cpms100))), main = "log2cpm filtered read counts \n 10>0 orthogenes")

# PCA plot of filtered CPM

par(mfrow=c(1,1))
x <- cpms100
plot(x.pca.cpms_filter100)
x.pca <- prcomp(t(x), scale.=TRUE, center=TRUE)
plot(x.pca)
summary(x.pca)

Leg = c("Light Sky Blue", "Light Sky Blue", "Light Sky Blue", "Light Sky Blue", 1,1,1,1,1,1,1,"Royal BLue",
        "Royal Blue", "Royal Blue",  "Royal Blue", "Royal Blue",  "Royal Blue", "Royal Blue")

ind = c(16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16)

plot(x.pca$x[, 1], x.pca$x[, 2], col = Leg, pch = ind, xlab ="PC1 (40.23% of variance)",   ylab = "PC2 (9.32% of variance)", main = "Gene expression (log2 CPM) \n filter 10>0 (n=17354)"); legend(50,70, pch = c(16,16,16), bty="n", col = c("Light Sky Blue", "Royal Blue", 1), c("Human-YRI", "Human-Caucasian", "Chimp" ))

# construct log2cpm table

RNAlog2cpm <- cpms100

# Pull out gene name and plot counts

RNAlog2cpm[rownames(RNAlog2cpm) == "ENSG00000181449", ]

cpms[rownames(cpms) == "ENSG00000179059",  ]

cpmsordered <-cpms[,c(5,6,7,8,9,10,11,1,2,3,4,12,13,14,15,16,17)]
cpmsordered_human <-cpms[,c(1,2,3,4,12,13,14,15,16,17)]
cpmsordered_chimp <-cpms[,c(5,6,7,8,9,10,11)]
colours <- c(1,1,1,1,1,1,1,"cornflower blue", "cornflower blue", "cornflower blue", "cornflower blue",
             "cornflower blue", "cornflower blue","cornflower blue","cornflower blue","cornflower blue",
             "cornflower blue")

barplot(cpmsordered[rownames(cpmsordered) == "ENSG00000001460",  ], las=2, ylab = "log2 cpm", main = "STPG1 expression",
        col=colours)

cpmsordered_human_REX1 <- cpmsordered_human[rownames(cpmsordered_human) == "ENSG00000179059",  ]
cpmsordered_chimp_REX1 <- cpmsordered_chimp[rownames(cpmsordered_chimp) == "ENSG00000179059",  ]
boxplot(cpmsordered_human_REX1, cpmsordered_chimp_REX1, ylab= "log2 cpm", names=c("human", "chimp"), main= "REX1")

# Get log2cpm for each gene with hg19 coordinates
# need genome coordinates of exons in count table from SZ

OrthoExon_hg19Exons<-read.table(".....", header=F)

OrthoExon_hg19Exons_cpms <-merge(RNAlog2cpm, OrthoExon_hg19Exons, by=0, all=FALSE)









