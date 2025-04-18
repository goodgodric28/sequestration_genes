####################################
# Ho_transcriptome_DESeq2_CNIvsSUB.R
# Written by: Jessica A. Goodheart
# Last Updated: 17 April 2025
# Purpose: To analyze differential expression data from Hermissenda and compare with Berghia
####################################

####################################
# Initial setup
####################################
# Set the working directory
directory <- "/Users/jessicagoodheart/Library/CloudStorage/OneDrive-AMNH/Goodheart/1_Research/1_Projects_and_Papers/1_Berghia/berghia_cnidosac_genes/Hermissenda/DE_analysis"
setwd(directory)

outputPrefix <- "Ho_tissues_DESeq2_CNIvsSUB"

# Call packages
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsubread)
library(Rsamtools)
library(rtracklayer)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library("genefilter")
library(dplyr)
library(edgeR)
library(stringr)
library(tidyr)
library(cowplot)
library(ggrepel)
library(splitstackshape)
library(GO.db)
library(statmod)
library(ggpolypath)
library(venn)
library("grDevices")

## Load in metadata
# Conditions
conditions<-as.data.frame(read.csv("inputs/sample_info_JAG.csv",sep="\t"))
conditions

# Counts
counts<-read.csv("inputs/hermissenda_counts.csv",row.names = 1)

####################################
# DESeq2 Analyses
# Modified from: https://benchtobioinformatics.wordpress.com/category/dexseq/
####################################
## Prep conditions and counts tables
# Subset conditions to remove foot and cer samples
c<-subset(conditions,Tissue !='foot')
c<-subset(c,Tissue !='cer')

# subset using counts table
subcounts<-counts[,c(grep("cni|sub",colnames(counts)))]

# remove single cni sample C_Ho061_cni, no matching sub sample
subcounts<-subset(subcounts, select = c(-C_Ho061_cni_He_op_SRR1950939_t_c95_filt_full.sorted.bam))
c<-c[!grepl("C_Ho061_cni",c$SampleName),]
rownames(subcounts) <- sapply(strsplit(rownames(subcounts), "~~"), "[[", 2)

## Create DESeq2 data set from counts
# DESeq2 matrix and analysis
ddsMat <- DESeqDataSetFromMatrix(subcounts,colData=c,design = ~Treatment+Tissue)
dds <- DESeq(ddsMat)
res <- results(dds)

# Save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =FALSE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-full-normalized-counts.csv"), quote=FALSE, row.names=FALSE)

resdata.2 <- resdata[order(resdata$padj),]
write.csv(resdata.2, file = paste0(outputPrefix, "-full-normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

resdata.3 <- subset(resdata.2,resdata.2$log2FoldChange < -2)
write.csv(resdata.3, file = paste0(outputPrefix, "-DCupreg-normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

resdata.4 <- subset(resdata.2,resdata.2$log2FoldChange > 2)
write.csv(resdata.4, file = paste0(outputPrefix, "-PCupreg-normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

resdata.5 <- subset(resdata.3,resdata.3$padj < 0.05)
write.csv(resdata.5, file = paste0(outputPrefix, "-DCupreg-padj_normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

# order results by padj value (most significant to least)
res <- res[order(res$padj), ]

# Produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"), quote=FALSE)

# Replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. Recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates. 
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)

write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"), quote=FALSE)

resClean <- results(ddsClean)
resClean = subset(resClean, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"), quote=FALSE)

# List and blast results of upregulated genes in the distal ceras
resClean.DC <- subset(resClean, log2FoldChange < -2)
resClean.PC <- subset(resClean, log2FoldChange > 2)

####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
####################################################################################

# Transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# Plot to show effect of transformation
# Axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()

# Save normalized values
write.table(cbind(as.data.frame(resdata$gene), as.data.frame(assay(rld))), 
            file = paste0(outputPrefix, "-rlog-transformed-counts-filtered.txt"), sep = '\t', quote=FALSE)
write.table(cbind(as.data.frame(resdata$gene), as.data.frame(assay(vsd))), 
            file = paste0(outputPrefix, "-vst-transformed-counts-filtered.txt"), sep = '\t', quote=FALSE)

# Clustering analysis
# Excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")

select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:500]
distsRL <- dist(t(assay(rld)[select,]))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),rld$SampleName)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(paste0(outputPrefix, "-clustering.png"), width=1100, height=1000)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13), cexRow = 1.7, cexCol = 1.7, dendrogram="row", revC=TRUE)
dev.off()

tiff(paste0(outputPrefix, "-clustering.tif"), width=1100, height=1000)
heatmap.2(mat, trace = "none", col = hmcol, margin = c(13,13), cexRow = 1.7, cexCol = 1.7, revC=TRUE)
dev.off()

# Principal components plot shows additional but rough clustering of samples
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))

sampleCondition2 <- c$Tissue
sampleCondition2 <- sub("cni", "Distal Ceras", sampleCondition2)
sampleCondition2 <- sub("sub", "Proximal Ceras", sampleCondition2)
scores <- data.frame(pc$x, sampleCondition2)

pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = sampleCondition2)) +
  geom_point(size = 5) +
  scale_color_manual(values=c("darkgoldenrod1","blue3")) +
  labs(color = "Tissue") +
  theme(
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA),
    text = element_text(size = 20)
  )

ggsave(pcaplot,file=paste0(outputPrefix, "-pcaplot.png"),device=png)
dev.off()

## Write results for a subset of genes
sum(res$padj < 0.05, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
resdataSig <- subset(resdata, padj < 0.1)
head(resSig[ order( resSig$log2FoldChange ), ])
write.csv(resdataSig, file=paste0(outputPrefix, "-diffexpr_results_padj0.1.csv"))

rld.sub <- rld[c(rownames(resSig)),]
sampleDists.sub <- dist( t( assay(rld.sub) ) )
sampleDists.sub

sampleDistMatrix.sub <- as.matrix( sampleDists.sub )
rownames(sampleDistMatrix.sub) <- paste( rld.sub$condition )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc.sub <- hclust(sampleDists.sub)

pdf(paste0(outputPrefix, '-DE_heatmap.pdf'))
heatmap.2( sampleDistMatrix.sub, Rowv=as.dendrogram(hc.sub),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labRow=c$SampleName, labCol="")
dev.off()

tiff(paste0(outputPrefix, '-DE_heatmap2.pdf'), width=1200, height=1000)
heatmap.2( sampleDistMatrix.sub, Rowv=as.dendrogram(hc.sub),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labRow=c$SampleName, labCol="")
dev.off()

# Top gene
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup="Tissue")


data <- plotCounts(dds, gene=topGene, intgroup=c("Treatment","Tissue"), returnData=TRUE)
ggplot(data, aes(x=Tissue, y=count, fill=Tissue)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")
ggplot(data, aes(x=Tissue, y=count, color=Tissue, group=Tissue)) +
  scale_y_log10() + 
  geom_point() + geom_line()



################################################
### Venn Diagram from presence/absence data ####
################################################
# Create presence/absence data frames for venn diagram construction
dds_df <- as.data.frame(counts(dds,normalized=T)) 
colnames(dds_df) <- c$SampleName
venn_df <- data.frame(matrix(ncol = 0, nrow = nrow(dds_df)))
rownames(venn_df) <- rownames(dds_df)

# Calculate average expression for each tissue type in new data frame
venn_df$cni <- (dds_df$C_Ho057_cni + dds_df$C_R0901_cni + dds_df$C_R0911_cni + 
                 dds_df$N_Ho046_cni + dds_df$N_Ho047_cni + dds_df$N_R0900_cni +
                 dds_df$N_R0909_cni)/7
venn_df$sub <- (dds_df$C_Ho057_sub + dds_df$C_R0901_sub + dds_df$C_R0911_sub + 
                 dds_df$N_Ho046_sub + dds_df$N_Ho047_sub + dds_df$N_R0900_sub +
                 dds_df$N_R0909_sub)/7

# Determine presence absence using average expression: Normalized counts > 0.5 considered present
venn_df[venn_df > 0.5] <- as.integer(1) 
venn_df[venn_df <= 0.5] <- as.integer(0) 
colnames(venn_df)<-c("Distal Ceras",  "Proximal Ceras")

# Plot venn diagram
cols = c("darkgoldenrod1", "blue3")

png(file=paste0(outputPrefix, "-venn-diagram.png"), width=1000, height=1000)
venn::venn(venn_df, zcolor=cols, ilcs = 3, sncs = 3, plotsize=100, ilabels="counts")
dev.off()

#######################################
## Connect Hermissenda to Berghia
#######################################
# Hermissenda data
diff_genes <- read.csv("Ho_tissues_DESeq2_CNIvsSUB-DCupreg-padj_normalized-counts-sorted.csv", row.names = NULL)
colnames(diff_genes)[1] <- "protein_id"

# Pull in full OrthoFinder results with Berghia and modify df for use
Bs.diffexp.data<- read.csv("inputs/Bs_tissues_DCvPC_DESeq2-DCupreg-normalized-counts-sorted.csv")
OF.data <- read.delim("inputs/Orthogroups.tsv", sep="\t")
OF.data.sub.allData <- OF.data[,c(1,10,29)]
OF.data.sub.bothSpOnly <- OF.data.sub.allData[!OF.data.sub.allData$Berghia_stephanieae_brakerRMplusISO==""|OF.data.sub.allData$Hermissenda_crassicornis=="", ] 
OF.data.sub.bothSpOnly.split <- separate_longer_delim(OF.data.sub.bothSpOnly, 
                                                      "Berghia_stephanieae_brakerRMplusISO",
                                                      ", ")
OF.data.sub.bothSpOnly.split <- separate_longer_delim(OF.data.sub.bothSpOnly.split, 
                                                      "Hermissenda_crassicornis",
                                                      ", ")
OF.data.sub.bothSpOnly.split$Berghia_stephanieae_brakerRMplusISO <- gsub("\\..*","",OF.data.sub.bothSpOnly.split$Berghia_stephanieae_brakerRMplusISO)
                                                        
# Subset df for only upregulated genes in one taxon - Hermissenda
OF.data.sub.diffgenes.Ho <- OF.data.sub.bothSpOnly.split[OF.data.sub.bothSpOnly.split$Hermissenda_crassicornis %in% diff_genes$protein_id,]
write.csv(OF.data.sub.diffgenes.Ho, file=paste0(outputPrefix, "-orthofinder_subset_Ho_diffexp.csv"))

diff_exp_genes.subOF.Ho <- diff_genes[diff_genes$protein_id %in% OF.data.sub.diffgenes.Ho$Hermissenda_crassicornis,]
write.csv(diff_exp_genes.subOF.Ho, file=paste0(outputPrefix, "-diff_exp_genes_subOF_Ho.csv"))

# Subset df for only upregulated genes in one taxon - Berghia
OF.data.sub.diffgenes.Bs <- OF.data.sub.bothSpOnly.split[OF.data.sub.bothSpOnly.split$Berghia_stephanieae_brakerRMplusISO %in% Bs.diffexp.data$gene,]
write.csv(OF.data.sub.diffgenes.Bs, file=paste0(outputPrefix, "-orthofinder_subset_Bs_diffexp.csv"))

diff_exp_genes.subOF.Bs <- Bs.diffexp.data[Bs.diffexp.data$gene %in% OF.data.sub.diffgenes.Bs$Berghia_stephanieae_brakerRMplusISO,]
write.csv(diff_exp_genes.subOF.Bs, file=paste0(outputPrefix, "-diff_exp_genes_subOF_Bs.csv"))

# Subset df for only upregulated genes in both taxa
OF.data.sub.diffgenes.Ho.Bs <- OF.data.sub.diffgenes.Ho[OF.data.sub.diffgenes.Ho$Berghia_stephanieae_brakerRMplusISO %in% Bs.diffexp.data$gene,]
write.csv(OF.data.sub.diffgenes.Ho.Bs, file=paste0(outputPrefix, "-orthofinder_subset_HoBs_diffexp.csv"))

#### FOR BOTH SPECIES ###
# List of orthogroups with their respective annotations based on Berghia stephanieae genome
Bs.annot.data <- read.delim("inputs/Bs_protein-blasthit-ALL-info.txt", "\t", header=TRUE)
Bs.annot.data$prot_id <- gsub("\\..*","",Bs.annot.data$prot_id)

Bs.annot.data.sub <- Bs.annot.data[Bs.annot.data$prot_id %in% OF.data.sub.diffgenes.Ho.Bs$Berghia_stephanieae_brakerRMplusISO,c(1:3,5:8)]
Bs.annot.data.sub <- Bs.annot.data.sub[!duplicated(Bs.annot.data.sub$prot_id), ]
write.csv(Bs.annot.data.sub, file=paste0(outputPrefix, "-Bs_annot_subset_diffexpHoBs.csv"), row.names = F)

# Pull upregulated genes from both taxa together with Berghia annotations
colnames(Bs.annot.data.sub)[1] <- "Berghia_stephanieae_brakerRMplusISO"
OF.data.HoBs.annot <- left_join(OF.data.sub.diffgenes.Ho.Bs,Bs.annot.data.sub,
                                by = "Berghia_stephanieae_brakerRMplusISO")
OF.data.HoBs.annot.sorted <- OF.data.HoBs.annot[order(match(OF.data.HoBs.annot$Berghia_stephanieae_brakerRMplusISO,
                                                            diff_exp_genes.subOF.Bs$gene)),]
logdata <- c()
for (x in OF.data.HoBs.annot.sorted$Berghia_stephanieae_brakerRMplusISO) {
  l <- diff_exp_genes.subOF.Bs$log2FoldChange[grep(paste0("^",x,"$"), diff_exp_genes.subOF.Bs$gene)]
  logdata <- c(logdata,l)
}

OF.data.HoBs.annot.sorted$log2FoldChange <- logdata
OF.data.HoBs.annot.sorted <- OF.data.HoBs.annot.sorted[,c(1:3,10,4:9)]

write.csv(OF.data.HoBs.annot.sorted, file=paste0(outputPrefix, "orthofinder_subset_HoBs_diffexp_with_annot_sorted.csv"), row.names = F)
write.table(unique(OF.data.HoBs.annot$Orthogroup), file=paste0(outputPrefix, "-OF_clusters_HoBs.txt"), row.names=F, col.names=F, quote=F)

OF.data.HoBs.annot.sorted.clean <- OF.data.HoBs.annot.sorted[,c(1,4:5,8)]
OF.data.HoBs.annot.sorted.clean <- OF.data.HoBs.annot.sorted %>% 
  group_by(Orthogroup) %>% 
  summarise(log2FoldChange=max(log2FoldChange),
            db_ids = paste(unique(db_id), collapse = ","),
            Protein.names = paste(unique(Protein.names), collapse = ","))

OF.data.HoBs.annot.sorted.clean <-OF.data.HoBs.annot.sorted.clean[order(OF.data.HoBs.annot.sorted.clean$log2FoldChange),]
write.csv(OF.data.HoBs.annot.sorted.clean, file=paste0(outputPrefix, "orthofinder_subset_HoBs_diffexp_with_annot_sorted_clean.csv"), row.names = F)


#### FOR HO ONLY ###
# List of orthogroups with their respective annotations based on Berghia stephanieae genome
Bs.annot.data.sub.Ho <- Bs.annot.data[Bs.annot.data$prot_id %in% OF.data.sub.diffgenes.Ho$Berghia_stephanieae_brakerRMplusISO,c(1:3,5:8)]
Bs.annot.data.sub.Ho <- Bs.annot.data.sub.Ho[!duplicated(Bs.annot.data.sub.Ho$prot_id), ]
write.csv(Bs.annot.data.sub.Ho, file=paste0(outputPrefix, "-Bs_annot_subset_diffexpHo.csv"), row.names = F)

# Pull upregulated genes from both taxa together with Berghia annotations
colnames(Bs.annot.data.sub.Ho)[1] <- "Berghia_stephanieae_brakerRMplusISO"
OF.data.Ho.annot <- left_join(OF.data.sub.diffgenes.Ho,Bs.annot.data.sub.Ho,
                                by = "Berghia_stephanieae_brakerRMplusISO")
OF.data.Ho.annot.sorted <- OF.data.Ho.annot[order(match(OF.data.Ho.annot$Berghia_stephanieae_brakerRMplusISO,
                                                            diff_exp_genes.subOF.Bs$gene)),]

write.csv(OF.data.Ho.annot.sorted, file=paste0(outputPrefix, "-orthofinder_subset_HoONLY_diffexp_with_annot_sorted.csv"), row.names = F)
write.table(unique(OF.data.Ho.annot$Orthogroup), file=paste0(outputPrefix, "-OF_clusters_HoONLY.txt"), row.names=F, col.names=F, quote=F)

#######################################
## Venn Diagram for Both Species
#######################################
# Venn diagram of upregulated orthologs overlap (singleton clusters excluded)
OF.data.sub.allData.venn <- OF.data.sub.allData
OF.data.sub.allData.venn[OF.data.sub.allData.venn == ""] <- 0
OF.data.sub.allData.venn$Berghia_stephanieae_brakerRMplusISO[OF.data.sub.allData.venn$Berghia_stephanieae_brakerRMplusISO != 0] <- 1
OF.data.sub.allData.venn$Hermissenda_crassicornis[OF.data.sub.allData.venn$Hermissenda_crassicornis != 0] <- 1

# Plot all genes
cols = c("#66C2A5", "#FC8D62")

tiff("figures/BsHo-venn-diagram-all.tiff", width=1000, height=1000)
venn::venn(OF.data.sub.allData.venn[,2:3], zcolor=cols, ilcs = 5.5, sncs = 2, plotsize=100, ilabels="counts")
dev.off()

OF.data.sub.upreg.venn <- as.data.frame(OF.data.sub.allData.venn[,1])
colnames(OF.data.sub.upreg.venn) <- "Orthogroups"
OF.data.sub.upreg.venn$Berghia <- as.integer(OF.data.sub.upreg.venn[,1] %in% OF.data.sub.diffgenes.Bs[,1])
OF.data.sub.upreg.venn$Hermissenda <- as.integer(OF.data.sub.upreg.venn[,1] %in% OF.data.sub.diffgenes.Ho[,1])
OF.data.sub.upreg.venn <- OF.data.sub.upreg.venn[rowSums(OF.data.sub.upreg.venn[,2:3])>0,]

tiff("figures/BsHo-venn-diagram-upreg.tiff", width=1000, height=1000)
venn::venn(OF.data.sub.upreg.venn[,2:3], zcolor=cols, ilcs = 5.5, sncs = 2, plotsize=100, ilabels="counts")
dev.off()

#######################################
## Gene expression plots for orthologs
#######################################

for (x in 1:nrow(OF.data.HoBs.annot.sorted)) {
  og <- OF.data.HoBs.annot.sorted$Orthogroup[x]
  Bs <- OF.data.HoBs.annot.sorted$Berghia_stephanieae_brakerRMplusISO[x]
  Ho <- OF.data.HoBs.annot.sorted$Hermissenda_crassicornis[x]
  
  # Berghia
  r.Bs <- grep(Bs, Bs.diffexp.data$gene)
  c.Bs <- as.vector(t(Bs.diffexp.data[r.Bs,c(8:13)]))
  
  # Hermissenda
  r.Ho <- grep(Ho, diff_genes$protein_id)
  c.Ho <- as.vector(t(diff_genes[r.Ho,8:21]))
  counts.data <- append(c.Bs, c.Ho)
  
  # All Counts
  counts.data <- data.frame(counts.data)
  counts.data$Species <- c(rep("Bs",6),rep("Ho",14))
  counts.data$Tissue <- c(rep("Distal",3),rep("Proximal",3),
                          "Distal","Proximal","Distal","Proximal","Distal","Proximal","Distal","Proximal",
                          "Distal","Proximal","Distal","Proximal","Distal","Proximal")
  
  # Summary
  counts.summary <- counts.data %>% group_by(Species, Tissue) %>% summarize(mean_counts=mean(counts.data), sd_counts=sd(counts.data))
  counts.summary <- as.data.frame(counts.summary)
  
  # Figure
  tiff(file=paste0("figures/genes/gene", x, "-", og, ".tiff"), width=800, height=1000)
  print(ggplot(counts.data, aes(x=Tissue, y=counts.data, color=Species, group=Tissue)) +
    scale_y_log10(limits=c(1,1e5)) + 
    geom_point(size=10)  + geom_line() +
    ylab("Normalized Counts (by species)") +
    theme_bw() +
    theme(text = element_text(size = 40)))
  dev.off()
  
  tiff(file=paste0("figures/genes/gene", x, "-", og, "-mean.tiff"), width=800, height=1000)
  print(ggplot(counts.summary, aes(x=Tissue, y=mean_counts, color=Species, group=Species)) +
          scale_y_log10(limits=c(1,1e5)) + 
          geom_point(size=10) + 
          geom_line() +
          ylab("Normalized Counts (by species)") +
          geom_errorbar(aes(ymin=mean_counts-sd_counts, ymax=mean_counts+sd_counts), width=.2) +
          theme_bw() +
          theme(text = element_text(size = 40)))
  dev.off()
}

