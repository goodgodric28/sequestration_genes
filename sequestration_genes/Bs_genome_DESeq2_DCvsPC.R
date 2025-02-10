####################################
# Bs_genome_DESeq2_DCvsPC.R
# Written by: Jessica A. Goodheart
# Last Updated: 9 February 2025
# Purpose: To analyze differential expression data from Berghia distal vs proximal ceras
####################################

####################################
# Initial setup
####################################

# Set the working directory
directory <- "[DIRECTORY]"
read_dir <- "[DIRECTORY]/read_files/"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Bs_tissues_DCvPC_DESeq2"

####################################
# DESeq2 Analyses
# Modified from: https://benchtobioinformatics.wordpress.com/category/dexseq/
####################################
# Load DESeq2 library for differential expression
library("DESeq2")
library("dplyr")

# Read in counts files for each tissue
sampleFiles<- c("D1-BsV1_genome.counts",
                "D2-BsV1_genome.counts",
                "D4-BsV1_genome.counts",
                "P1-BsV1_genome.counts",
                "P2-BsV1_genome.counts",
                "P3-BsV1_genome.counts")

# Set sample names for counts
sampleNames <- c("Distal Ceras 1", "Distal Ceras 2", "Distal Ceras 3",
                 "Proximal Ceras 1", "Proximal Ceras 2", "Proximal Ceras 3")

# Set sample conditions for counts
sampleCondition <- c("distCeras","distCeras","distCeras",
                     "proxCeras","proxCeras","proxCeras")

# Create DESeq2 sample table
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

# Set treatment levels for DESeq2 analysis
treatments = c("distCeras","proxCeras")

# Create DESeq2 data set from counts
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = read_dir,
                                       design = ~1 + condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# Calculate differential expression based on the negative binomial (i.e., Gamma-Poisson) distribution
dds <- DESeq(ddsHTSeq)
res <- results(dds)

# Save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =FALSE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-full-normalized-counts.csv"), quote=FALSE, row.names=FALSE)

resdata.2 <- resdata[order(resdata$padj),]
write.csv(resdata.2, file = paste0(outputPrefix, "-full-normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

resdata.3 <- subset(resdata.2,resdata.2$log2FoldChange < 0)
write.csv(resdata.3, file = paste0(outputPrefix, "-DCupreg-normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

resdata.4 <- subset(resdata.2,resdata.2$log2FoldChange > 0)
write.csv(resdata.4, file = paste0(outputPrefix, "-PCupreg-normalized-counts-sorted.csv"), quote=FALSE, row.names=FALSE)

# order results by padj value (most significant to least)
res <- res[order(res$padj),]

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

# Pull in blast annotations on the filtered proteome dataset from Goodheart et al 2024
blastp = read.csv("inputs/Bs_protein-blasthit-ALL-info.txt", sep="\t", 
                  header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)
blastp$prot_id = gsub(".t[0-9]*", "", as.character(blastp$prot_id))
blastp.sub = dplyr::distinct(blastp, prot_id, .keep_all=TRUE)
blastp.sub = blastp.sub[,c(1:4,6,10:11)]

blastp.annot = data.frame(matrix(ncol = ncol(blastp.sub), nrow = 0))
colnames(blastp.annot) = colnames(blastp.sub)
for (i in rownames(resClean.DC)) {
  ie = paste("\\b", i, "\\b", sep="")
  s = subset(blastp.sub, grepl(ie, blastp.sub$prot_id))
  if(nrow(s)==0) 
    s = data.frame("prot_id"=i, "db_id"=NA, "evalue"=NA, "Entry"=NA,
                   "Protein.names"=NA, "Refseq"=NA, "GO.terms"=NA,stringsAsFactors = FALSE)
  blastp.annot = rbind(blastp.annot, s)
}

resClean.DC <- data.frame(resClean.DC, blastp.annot)
write.csv(as.data.frame(resClean.DC),file = paste0(outputPrefix, "-replaceoutliers-DCupreg-blasthits.csv"), quote=FALSE)

stats.DC <- c(nrow(resClean.DC), nrow(resClean.DC)-sum(is.na(resClean.DC$Protein.names)), 100*((nrow(resClean.DC)-sum(is.na(resClean.DC$Protein.names)))/nrow(resClean.DC)))
names(stats.DC) <- c("# of upreg genes in DC", "# DC upreg genes annotated", "% DC upreg genes annotated")
write.table(stats.DC, file = paste0(outputPrefix, "-DCupreg-stats.txt"), sep = '\t', quote=FALSE, col.names=FALSE)

blastp.annot.PC = data.frame(matrix(ncol = ncol(blastp.sub), nrow = 0))
colnames(blastp.annot.PC) = colnames(blastp.sub)
for (i in rownames(resClean.PC)) {
  ie = paste("\\b", i, "\\b", sep="")
  s = subset(blastp.sub, grepl(ie, blastp.sub$prot_id))
  if(nrow(s)==0) 
    s = data.frame("prot_id"=i, "db_id"=NA, "evalue"=NA, "Entry"=NA,
                   "Protein.names"=NA, "Refseq"=NA, "GO.terms"=NA,stringsAsFactors = FALSE)
  blastp.annot.PC = rbind(blastp.annot.PC, s)
}

resClean.PC <- data.frame(resClean.PC, blastp.annot.PC)
write.csv(as.data.frame(resClean.PC),file = paste0(outputPrefix, "-replaceoutliers-PCupreg-blasthits.csv"), quote=FALSE)

stats.PC <- c(nrow(resClean.PC), nrow(resClean.PC)-sum(is.na(resClean.PC$Protein.names)), 100*((nrow(resClean.PC)-sum(is.na(resClean.PC$Protein.names)))/nrow(resClean.PC)))
names(stats.PC) <- c("# of upreg genes in PC", "# PC upreg genes annotated", "% PC upreg genes annotated")
write.table(stats.PC, file = paste0(outputPrefix, "-PCupreg-stats.txt"), sep = '\t', quote=FALSE, col.names=FALSE)

####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) rlog stabilization and variance stabiliazation
# 2) variance stabilization plot
# 3) heatmap of clustering analysis
# 4) PCA plot
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
rownames(mat) <- colnames(mat) <- with(colData(dds),sampleNames)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(paste0(outputPrefix, "-clustering.png"), width=1100, height=1000)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13), cexRow = 1.7, cexCol = 1.7, dendrogram="row", revC=TRUE)
dev.off()

tiff(paste0(outputPrefix, "-clustering.tif"), width=1100, height=1000)
heatmap.2(mat, trace = "none", col = hmcol, margin = c(13,13), cexRow = 1.7, cexCol = 1.7, revC=TRUE)
dev.off()

# Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))

condition <- treatments
sampleCondition2 <- c("Distal Ceras","Distal Ceras","Distal Ceras",
                      "Proximal Ceras","Proximal Ceras","Proximal Ceras")
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
           margins=c(2,10), labRow=sampleNames, labCol="")
dev.off()

tiff(paste0(outputPrefix, '-DE_heatmap2.pdf'), width=1200, height=1000)
heatmap.2( sampleDistMatrix.sub, Rowv=as.dendrogram(hc.sub),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labRow=sampleNames, labCol="")
dev.off()

################################################
### Venn Diagram from presence/absence data ####
################################################
library(ggpolypath)
library(RColorBrewer)
library(venn)

# Create presence/absence data frames for venn diagram construction
dds_df <- as.data.frame(counts(dds,normalized=T)) 
venn_df <- data.frame(matrix(ncol = 0, nrow = nrow(dds_df)))
rownames(venn_df) <- rownames(dds_df)

# Calculate average expression for each tissue type in new data frame
venn_df$DC <- (dds_df$`Distal Ceras 1` + dds_df$`Distal Ceras 2` + dds_df$`Distal Ceras 3`)/3
venn_df$PC <- (dds_df$`Proximal Ceras 1` + dds_df$`Proximal Ceras 2` + dds_df$`Proximal Ceras 3`)/3

# Determine presence absence using average expression: Normalized counts > 0.5 considered present
venn_df[venn_df > 0.5] <- as.integer(1) 
venn_df[venn_df <= 0.5] <- as.integer(0) 
colnames(venn_df)<-c("Distal Ceras",  "Proximal Ceras")

# Plot venn diagram
cols = c("darkgoldenrod1", "blue3")

png(file=paste0(outputPrefix, "-venn-diagram.png"), width=1000, height=1000)
venn::venn(venn_df, zcolor=cols, ilcs = 3, sncs = 3, plotsize=100, ilabels="counts")
dev.off()

################################################
##### Gene list from presence/absence data #####
################################################
library(R.utils)
library(ProtDomSeq)

# Change presence/absence data set column names
colnames(venn_df) <- c("DC",  "PC")

# Pull out gene names for genes only expressed in certain tissues and save
dc_only <- rownames(filter(venn_df, DC=="1" & PC=="0"))
pc_only <- rownames(filter(venn_df, DC=="0" & PC=="1"))

# Write out tables with genes only expressed in (a) certain tissue(s)
write.table(dc_only, paste0("gene_annotations/", outputPrefix, "-DC-only-genelist.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(pc_only, paste0("gene_annotations/", outputPrefix, "-PC-only-genelist.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE) 

# Pull in interproscan annotation data 
annot <- read.csv("inputs/berghia_RM_iso_anysupport_ipscan_2021_11.tsv", sep="\t")
annot_edit <- annot[-c(2:3,7:11)]
colnames(annot_edit) <- c("prot_id", "source", "type", "label", "ipr_id", "go_term")

# Pull in interproscan database accessions
annot_db <- read.table("inputs/functional_annotation.txt", header=TRUE, sep="\t")

# Pull in blastp hits separately
blastp <- read.table("inputs/blastp_annotation.gff", sep="\t")
blastp <- blastp[-c(4:5,7:8)]
colnames(blastp) <- c("chromosome","source","match","length","hit")

# Pull out annotation data for upregulated genes
dc_up <- rownames(resClean.DC)
pc_up <- rownames(resClean.PC)

### Distal Ceras Upregulated ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
dc_up_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in dc_up) {
  colnames(dc_up_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  dc_up_annot <- rbind(dc_up_annot, s)
}

# Sort annotations by number of times present
dc_label_summary_up <- sort(table(dc_up_annot$label), decreasing=TRUE)
dc_label_summary_up <- dc_label_summary_up[!dc_label_summary_up %in% 0]

# Write out list of annotations to files
write.table(dc_up_annot, paste0("gene_annotations/", outputPrefix, "-DC-up-annotations.txt"), quote=FALSE, row.names = FALSE)
write.table(dc_label_summary_up, paste0("gene_annotations/", outputPrefix, "-DC-up-label-summary.txt"), quote=FALSE, col.names = FALSE)

# Sort go terms by number of times present
dc_goterm_summary <- sort(table(dc_up_annot$go_term), decreasing=TRUE)
dc_goterm_summary <- dc_goterm_summary[!dc_goterm_summary %in% 0]

# Write out list of go terms to file
write.table(dc_goterm_summary, paste0("gene_annotations/", outputPrefix, "-DC-up-goterms-summary.txt"), quote=FALSE, col.names = FALSE)

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
dc_up_acc <- data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                          "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in dc_up) {
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  dc_up_acc <- rbind(dc_up_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(dc_up_acc, paste0("gene_annotations/", outputPrefix, "-DC-up-accessions-summary.txt"), quote=FALSE, row.names = FALSE)

### Proximal Ceras Upregulated ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
pc_up_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in pc_up) {
  colnames(pc_up_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  pc_up_annot <- rbind(pc_up_annot, s)
}

# Sort annotations by number of times present
pc_label_summary_up <- sort(table(pc_up_annot$label), decreasing=TRUE)
pc_label_summary_up <- pc_label_summary_up[!pc_label_summary_up %in% 0]

# Write out list of annotations to files
write.table(pc_up_annot, paste0("gene_annotations/", outputPrefix, "-PC-up-annotations.txt"), quote=FALSE, row.names = FALSE)
write.table(pc_label_summary_up, paste0("gene_annotations/", outputPrefix, "-PC-up-label-summary.txt"), quote=FALSE, col.names = FALSE)

# Sort go terms by number of times present
pc_goterm_summary <- sort(table(pc_up_annot$go_term), decreasing=TRUE)
pc_goterm_summary <- pc_goterm_summary[!pc_goterm_summary %in% 0]

# Write out list of go terms to file
write.table(pc_goterm_summary, paste0("gene_annotations/", outputPrefix, "-PC-up-goterms-summary.txt"), quote=FALSE, col.names = FALSE)

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
pc_up_acc <- data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                          "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in pc_up) {
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  pc_up_acc <- rbind(pc_up_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(pc_up_acc, paste0("gene_annotations/", outputPrefix, "-PC-up-accessions-summary.txt"), quote=FALSE, row.names = FALSE)

################################################
######### TopGO Analysis - Upregulated #########
################################################
library(topGO)
library(patchwork)
library(stringr)

# Gene to GO term mappings
go.mappings <- readMappings(file = "inputs/go_terms_full_edited.txt") 
all_genes <- names(go.mappings)

## DC go term analysis ##
# Genes to go.mappings
dc_geneList <- factor(as.integer(all_genes %in% dc_up))
names(dc_geneList) <- all_genes

# TopGO - MF
GOdataMF.dc <- new("topGOdata", ontology = "MF", allGenes = dc_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherMF.dc <- runTest(GOdataMF.dc, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabMF.dc <- GenTable(GOdataMF.dc, raw.p.value = resultFisherMF.dc, topNodes = length(resultFisherMF.dc@score),
                   numChar = 120)
write.table(tabMF.dc, paste0("gene_annotations/", outputPrefix, "-DC-up-topgoMF.txt"), sep = '\t', quote=FALSE)

pvalFisherMF.dc <- score(resultFisherMF.dc)
hist(pvalFisherMF.dc, 50, xlab = "p-values")
allResMF.dc <- GenTable(GOdataMF.dc, parentchild = resultFisherMF.dc, topNodes = 20)


# TopGO - BP
GOdataBP.dc <- new("topGOdata", ontology = "BP", allGenes = dc_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherBP.dc <- runTest(GOdataBP.dc, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabBP.dc <- GenTable(GOdataBP.dc, raw.p.value = resultFisherBP.dc, topNodes = length(resultFisherBP.dc@score),
                  numChar = 120)
write.table(tabBP.dc, paste0("gene_annotations/", outputPrefix, "-DC-up-topgoBP.txt"), sep = '\t', quote=FALSE)

pvalFisherBP.dc <- score(resultFisherBP.dc)
hist(pvalFisherBP.dc, 50, xlab = "p-values")
allResBP.dc <- GenTable(GOdataBP.dc, parentchild = resultFisherBP.dc, topNodes = 20)

allResBP.dc$parentchild <- as.numeric(allResBP.dc$parentchild)
allResBP.dc <- allResBP.dc[allResBP.dc$parentchild < 0.05,] # filter terms for KS p<0.05
allResBP.dc <- allResBP.dc[,c("GO.ID","Term","parentchild")]
allResBP.dc

allResBP.dc$Term <- factor(allResBP.dc$Term, levels = rev(allResBP.dc$Term)) # fixes order
allResBP.dc$GO.ID <- factor(allResBP.dc$GO.ID, levels = rev(allResBP.dc$GO.ID)) # fixes order

DC.BP.plot <- ggplot(allResBP.dc,
       aes(x = Term, y = -log10(as.numeric(parentchild)), 
           size = -log10(as.numeric(parentchild)), fill = -log10(as.numeric(parentchild)))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(limits = c(0.5,14)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4', limits = c(0.5,14)) +
  
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'Distal Ceras - Biological Processes',
    subtitle = 'Top 20 terms ordered by ParentChild p-value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1, 2)) +
  
  scale_y_continuous(limits = c(0, 14)) +
  scale_x_discrete(labels=paste(allResBP.dc$GO.ID,"-",allResBP.dc$Term)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()

png(file=paste0(outputPrefix, "-DC-BP-GOterm-plot.png"), width=1000, height=500)
DC.BP.plot
dev.off()

 # TopGO - CC
GOdataCC.dc <- new("topGOdata", ontology = "CC", allGenes = dc_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherCC.dc <- runTest(GOdataCC.dc, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabCC.dc <- GenTable(GOdataCC.dc, raw.p.value = resultFisherCC.dc, topNodes = length(resultFisherCC.dc@score),
                  numChar = 120)
write.table(tabCC.dc, paste0("gene_annotations/", outputPrefix, "-DC-up-topgoCC.txt"), sep = '\t', quote=FALSE)


pvalFisherCC.dc <- score(resultFisherCC.dc)
hist(pvalFisherCC.dc, 50, xlab = "p-values")
allResCC.dc <- GenTable(GOdataCC.dc, parentchild = resultFisherCC.dc, topNodes = 20)

## PC go term analysis ##
# Genes to go.mappings
pc_geneList <- factor(as.integer(all_genes %in% pc_up))
names(pc_geneList) <- all_genes

# TopGO - MF
GOdataMF.pc <- new("topGOdata", ontology = "MF", allGenes = pc_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherMF.pc <- runTest(GOdataMF.pc, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabMF.pc <- GenTable(GOdataMF.pc, raw.p.value = resultFisherMF.pc, topNodes = length(resultFisherMF.pc@score),
                  numChar = 120)
write.table(tabMF.pc, paste0("gene_annotations/", outputPrefix, "-PC-up-topgoMF.txt"), sep = '\t', quote=FALSE)

pvalFisherMF.pc <- score(resultFisherMF.pc)
hist(pvalFisherMF.pc, 50, xlab = "p-values")
allResMF.pc <- GenTable(GOdataMF.pc, parentchild = resultFisherMF.pc, topNodes = 20)

# TopGO - BP
GOdataBP.pc <- new("topGOdata", ontology = "BP", allGenes = pc_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherBP.pc <- runTest(GOdataBP.pc, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabBP.pc <- GenTable(GOdataBP.pc, raw.p.value = resultFisherBP.pc, topNodes = length(resultFisherBP.pc@score),
                  numChar = 120)
write.table(tabBP.pc, paste0("gene_annotations/", outputPrefix, "-PC-up-topgoBP.txt"), sep = '\t', quote=FALSE)

pvalFisherBP.pc <- score(resultFisherBP.pc)
hist(pvalFisherBP.pc, 50, xlab = "p-values")
allResBP.pc <- GenTable(GOdataBP.pc, parentchild = resultFisherBP.pc, topNodes = 20)

allResBP.pc$parentchild <- as.numeric(allResBP.pc$parentchild)
allResBP.pc <- allResBP.pc[allResBP.pc$parentchild < 0.05,] # filter terms for KS p<0.05
allResBP.pc <- allResBP.pc[,c("GO.ID","Term","parentchild")]
allResBP.pc

allResBP.pc$Term <- factor(allResBP.pc$Term, levels = rev(allResBP.pc$Term)) # fixes order


PC.BP.plot <- ggplot(allResBP.pc,
       aes(x = Term, y = -log10(as.numeric(parentchild)), 
           size = -log10(as.numeric(parentchild)), fill = -log10(as.numeric(parentchild)))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(limits = c(0.5,14)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4', limits = c(0.5,14)) +
  
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'Proximal Ceras - Biological Processes',
    subtitle = 'Top 20 terms ordered by ParentChild p-value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1, 2)) +
  
  scale_y_continuous(limits = c(0, 14)) +
  scale_x_discrete(labels=paste(allResBP.pc$GO.ID,"-",allResBP.pc$Term)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()

png(file=paste0(outputPrefix, "-PC-BP-GOterm-plot.png"), width=1000, height=500)
PC.BP.plot
dev.off()

# TopGO - CC
GOdataCC.pc <- new("topGOdata", ontology = "CC", allGenes = pc_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherCC.pc <- runTest(GOdataCC.pc, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabCC.pc <- GenTable(GOdataCC.pc, raw.p.value = resultFisherCC.pc, topNodes = length(resultFisherCC.pc@score),
                  numChar = 120)
write.table(tabCC.pc, paste0("gene_annotations/", outputPrefix, "-PC-up-topgoCC.txt"), sep = '\t', quote=FALSE)

pvalFisherCC.pc <- score(resultFisherCC.pc)
hist(pvalFisherCC.pc, 50, xlab = "p-values")
allResCC.pc <- GenTable(GOdataCC.pc, parentchild = resultFisherCC.pc, topNodes = 20)


# GO Terms figure - both
BP.both <- (DC.BP.plot/PC.BP.plot) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides="collect") &
  theme(plot.tag = element_text(size = 40), 
        axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        legend.position = "right") 
BP.both

pdf(file=paste0(outputPrefix, "-BP-GOterm-plot-both.pdf"), height = 10, width = 14)
BP.both
dev.off()

################################################
######### Annotation Analysis - Only ##########
################################################

### Distal Ceras Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
dc_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in dc_only) {
  colnames(dc_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  dc_only_annot <- rbind(dc_only_annot, s)
}

# Sort annotations by number of times present
dc_label_summary <- sort(table(dc_only_annot$label), decreasing=TRUE)
dc_label_summary <- dc_label_summary[!dc_label_summary %in% 0]

# Write out list of annotations to files
write.table(dc_only_annot, paste0("gene_annotations/", outputPrefix, "-DC-only-annotations.txt"), quote=FALSE, row.names = FALSE)
write.table(dc_label_summary, paste0("gene_annotations/", outputPrefix, "-DC-only-label-summary.txt"), quote=FALSE, col.names = FALSE)

# Sort go terms by number of times present
dc_goterm_summary <- sort(table(dc_only_annot$go_term), decreasing=TRUE)
dc_goterm_summary <- dc_goterm_summary[!dc_goterm_summary %in% 0]

# Write out list of go terms to file
write.table(dc_goterm_summary, paste0("gene_annotations/", outputPrefix, "-DC-only-goterms-summary.txt"), quote=FALSE, col.names = FALSE)

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
dc_only_acc <- data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                          "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in dc_only) {
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  dc_only_acc <- rbind(dc_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(dc_only_acc, paste0("gene_annotations/", outputPrefix, "-DC-only-accessions-summary.txt"), quote=FALSE, row.names = FALSE)

### Proximal Ceras Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
pc_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in pc_only) {
  colnames(pc_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  pc_only_annot <- rbind(pc_only_annot, s)
}

# Sort annotations by number of times present
pc_label_summary <- sort(table(pc_only_annot$label), decreasing=TRUE)
pc_label_summary <- pc_label_summary[!pc_label_summary %in% 0]

# Write out list of annotations to files
write.table(pc_only_annot, paste0("gene_annotations/", outputPrefix, "-PC-only-annotations.txt"), quote=FALSE, row.names = FALSE)
write.table(pc_label_summary, paste0("gene_annotations/", outputPrefix, "-PC-only-label-summary.txt"), quote=FALSE, col.names = FALSE)

# Sort go terms by number of times present
pc_goterm_summary <- sort(table(pc_only_annot$go_term), decreasing=TRUE)
pc_goterm_summary <- pc_goterm_summary[!pc_goterm_summary %in% 0]

# Write out list of go terms to file
write.table(pc_goterm_summary, paste0("gene_annotations/", outputPrefix, "-PC-only-goterms-summary.txt"), quote=FALSE, col.names = FALSE)

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
pc_only_acc <- data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                          "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in pc_only) {
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  pc_only_acc <- rbind(pc_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(pc_only_acc, paste0("gene_annotations/", outputPrefix, "-PC-only-accessions-summary.txt"), quote=FALSE, row.names = FALSE)

################################################
############ TopGO Analysis - Only #############
################################################

## DC go term analysis ##
# Genes to go.mappings
dc_only_geneList <- factor(as.integer(all_genes %in% dc_only))
names(dc_only_geneList) <- all_genes

# TopGO - MF
GOdataMF.dc_only <- new("topGOdata", ontology = "MF", allGenes = dc_only_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherMF.dc_only <- runTest(GOdataMF.dc_only, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabMF.dc_only <- GenTable(GOdataMF.dc_only, raw.p.value = resultFisherMF.dc_only, topNodes = length(resultFisherMF.dc_only@score),
                     numChar = 120)
write.table(tabMF.dc_only, paste0("gene_annotations/", outputPrefix, "-DC-only-topgoMF.txt"), sep='\t', quote=FALSE)

pvalFisherMF.dc_only <- score(resultFisherMF.dc_only)
hist(pvalFisherMF.dc_only, 50, xlab = "p-values")
allResMF.dc_only <- GenTable(GOdataMF.dc_only, parentchild = resultFisherMF.dc_only, topNodes = 20)


# TopGO - BP
GOdataBP.dc_only <- new("topGOdata", ontology = "BP", allGenes = dc_only_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherBP.dc_only <- runTest(GOdataBP.dc_only, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabBP.dc_only <- GenTable(GOdataBP.dc_only, raw.p.value = resultFisherBP.dc_only, topNodes = length(resultFisherBP.dc_only@score),
                     numChar = 120)
write.table(tabBP.dc_only, paste0("gene_annotations/", outputPrefix, "-DC-only-topgoBP.txt"), sep='\t', quote=FALSE)

pvalFisherBP.dc_only <- score(resultFisherBP.dc_only)
hist(pvalFisherBP.dc_only, 50, xlab = "p-values")
allResBP.dc_only <- GenTable(GOdataBP.dc_only, parentchild = resultFisherBP.dc_only, topNodes = 20)

# TopGO - CC
GOdataCC.dc_only <- new("topGOdata", ontology = "CC", allGenes = dc_only_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherCC.dc_only <- runTest(GOdataCC.dc_only, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabCC.dc_only <- GenTable(GOdataCC.dc_only, raw.p.value = resultFisherCC.dc_only, topNodes = length(resultFisherCC.dc_only@score),
                     numChar = 120)
write.table(tabCC.dc_only, paste0("gene_annotations/", outputPrefix, "-DC-only-topgoCC.txt"), sep='\t', quote=FALSE)

pvalFisherCC.dc_only <- score(resultFisherCC.dc_only)
hist(pvalFisherCC.dc_only, 50, xlab = "p-values")
allResCC.dc_only <- GenTable(GOdataCC.dc_only, parentchild = resultFisherCC.dc_only, topNodes = 20)

## PC go term analysis ##
# Genes to go.mappings
pc_only_geneList <- factor(as.integer(all_genes %in% pc_only))
names(pc_only_geneList) <- all_genes

# TopGO - MF
GOdataMF.pc_only <- new("topGOdata", ontology = "MF", allGenes = pc_only_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherMF.pc_only <- runTest(GOdataMF.pc_only, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabMF.pc_only <- GenTable(GOdataMF.pc_only, raw.p.value = resultFisherMF.pc_only, topNodes = length(resultFisherMF.pc_only@score),
                     numChar = 120)
write.table(tabMF.pc_only, paste0("gene_annotations/", outputPrefix, "-PC-only-topgoMF.txt"), sep='\t', quote=FALSE)

pvalFisherMF.pc_only <- score(resultFisherMF.pc_only)
hist(pvalFisherMF.pc_only, 50, xlab = "p-values")
allResMF.pc_only <- GenTable(GOdataMF.pc_only, parentchild = resultFisherMF.pc_only, topNodes = 20)

# TopGO - BP
GOdataBP.pc_only <- new("topGOdata", ontology = "BP", allGenes = pc_only_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherBP.pc_only <- runTest(GOdataBP.pc_only, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabBP.pc_only <- GenTable(GOdataBP.pc_only, raw.p.value = resultFisherBP.pc_only, topNodes = length(resultFisherBP.pc_only@score),
                     numChar = 120)
write.table(tabBP.pc_only, paste0("gene_annotations/", outputPrefix, "-PC-only-topgoBP.txt"), sep='\t', quote=FALSE)

pvalFisherBP.pc_only <- score(resultFisherBP.pc_only)
hist(pvalFisherBP.pc_only, 50, xlab = "p-values")
allResBP.pc_only <- GenTable(GOdataBP.pc_only, parentchild = resultFisherBP.pc_only, topNodes = 20)

# TopGO - CC
GOdataCC.pc_only <- new("topGOdata", ontology = "CC", allGenes = pc_only_geneList, annot = annFUN.gene2GO, gene2GO = go.mappings)
resultFisherCC.pc_only <- runTest(GOdataCC.pc_only, algorithm="parentchild", statistic="fisher", scoreOrder="increasing")
tabCC.pc_only <- GenTable(GOdataCC.pc_only, raw.p.value = resultFisherCC.pc_only, topNodes = length(resultFisherCC.pc_only@score),
                     numChar = 120)
write.table(tabCC.pc_only, paste0("gene_annotations/", outputPrefix, "-PC-only-topgoCC.txt"), sep='\t', quote=FALSE)

pvalFisherCC.pc_only <- score(resultFisherCC.pc_only)
hist(pvalFisherCC.pc_only, 50, xlab = "p-values")
allResCC.pc_only <- GenTable(GOdataCC.pc_only, parentchild = resultFisherCC.pc_only, topNodes = 20)

################################################
############ Comparisons - GO Terms ############
################################################

### Make tidyverse object with counts data for comparisons ###
library(tidyr)
library(ggplot2)
library(GOfuncR)
library(rstatix)
library(ggpubr)
library(patchwork)

go.table <- read.table("inputs/go_terms_full_edited.txt", sep="\t")
colnames(go.table) <- c("genes","go_terms")

resdata.counts <- cbind(rownames(assay(vsd)),as.data.frame(assay(vsd)))
colnames(resdata.counts) <- c("genes","DC_1","DC_2","DC_3","PC_1","PC_2","PC_3")

counts.tidy <- as.data.frame(resdata.counts %>% 
  pivot_longer(
    cols = !genes,
    names_to = c("tissue", "sample"),
    names_pattern = "(.*)_(.)",
    values_drop_na = TRUE))

de.genes <- c(rownames(resClean.DC),rownames(resClean.PC))

### Analysis of GO terms - DC ###
# List of GO terms
gos.of.interest <- c("GO:0006897","GO:0006955","GO:0007586")
names(gos.of.interest) <- c("endocytosis","immune response","digestion")

# Get child lists for each broader term
i<-1
for (go in gos.of.interest) {
  n<-names(gos.of.interest)[i]
  assign(paste0(substr(n, start = 1, stop = 4)),get_child_nodes(go)[,2])
  i<-i+1
}

# Go terms analysis 
go.genes.table <- data.frame(GO.Term=character(),GO.ID=character(),genes.all=integer(),genes.de=integer())
my_comparisons <- list(c("DC","PC"))
i<-1
for (go in gos.of.interest) {
  n<-names(gos.of.interest)[i]
  symb <- substr(n, start = 1, stop = 4)
  print(c(symb, go))
  
  go.genes <- vector()
  for (g in eval(parse(text=paste0(symb)))) {
    go.genes <- c(go.genes, filter(go.table, grepl(g, go_terms))[,1])
  }
  go.genes <- unique(go.genes)
  
  counts.go <- subset(counts.tidy, counts.tidy$genes %in% go.genes)
  stats.go <- counts.go %>% group_by(genes, tissue) %>%
    get_summary_stats(value, type="mean_sd")
  
  go <- gsub("GO:", "", go)
  assign(paste0(symb, ".genes"),counts.go)
  
  plot.go <- ggplot(counts.go, aes(x=factor(tissue), y=value, fill=tissue)) + 
    #scale_y_log10() + 
    geom_violin()  +
    scale_x_discrete(labels=c("Distal Ceras", "Proximal Ceras")) +
    scale_fill_manual(values=c("darkgoldenrod1", "blue3")) +
    scale_y_continuous(limits = c(0, 21)) +
    labs(x ="Tissue", y = "VST transformed counts") +
    stat_compare_means(comparisons = my_comparisons, paired = TRUE, label="p.signif", method = "t.test")
  
  assign(paste0(symb, ".all"),plot.go)
  
  # Supplementary figures
  plot.go.line <- ggplot(counts.go, aes(x=factor(tissue), y=value, fill=tissue)) + 
    #scale_y_log10() + 
    geom_violin()  +
    geom_jitter(shape=16, position=position_jitter(0.1)) +
    geom_line(aes(group=genes)) +
    scale_x_discrete(labels=c("Distal Ceras", "Proximal Ceras")) +
    scale_fill_manual(values=c("darkgoldenrod1", "blue3")) +
    scale_y_continuous(limits = c(0, 21)) +
    labs(x ="Tissue", y = "VST transformed counts") +
    stat_compare_means(comparisons = my_comparisons, paired = TRUE, label="p.signif", method = "t.test")
  
  assign(paste0(symb, ".all.line"),plot.go.line)
  
  ggsave(file=paste0(outputPrefix, symb, "-GO_", go, "-plot-all-NL.png"), plot.go, width=10, height=10)
  ggsave(file=paste0(outputPrefix, symb, "-GO_", go, "-plot-all.png"), plot.go.line, width=10, height=10)
  
  # Only upregulated in either tissue
  counts.go.de <- subset(counts.go, counts.go$genes %in% de.genes)
  stats.go.de <- counts.go.de %>% group_by(genes, tissue) %>%
    get_summary_stats(value, type="mean_sd")

  plot.go.de <- ggplot(counts.go.de, aes(x=factor(tissue), y=value, fill=tissue)) + 
    #scale_y_log10() + 
    geom_violin()  +
    scale_x_discrete(labels=c("Distal Ceras", "Proximal Ceras")) +
    scale_fill_manual(values=c("darkgoldenrod1", "blue3")) +
    scale_y_continuous(limits = c(0, 21)) +
    labs(x ="Tissue", y = "VST transformed counts") +
    stat_compare_means(comparisons = my_comparisons, paired = TRUE, label="p.signif", method = "t.test")
  
  assign(paste0(symb, ".de"),plot.go.de)
  
  # Supplementary figures
  plot.go.de.line <- ggplot(counts.go.de, aes(x=factor(tissue), y=value, fill=tissue)) + 
    #scale_y_log10() + 
    geom_violin()  +
    geom_jitter(shape=16, position=position_jitter(0.1)) +
    geom_line(aes(group=genes)) +
    scale_x_discrete(labels=c("Distal Ceras", "Proximal Ceras")) +
    scale_fill_manual(values=c("darkgoldenrod1", "blue3")) +
    scale_y_continuous(limits = c(0, 21)) +
    labs(x ="Tissue", y = "VST transformed counts") +
    stat_compare_means(comparisons = my_comparisons, paired = TRUE, label="p.signif", method = "t.test")
  
  assign(paste0(symb, ".de.line"),plot.go.de.line)
  
  ggsave(file=paste0(outputPrefix, symb, "-GO_", go, "-plot-DE-NL.png"), plot.go.de, width=10, height=10)
  ggsave(file=paste0(outputPrefix, symb, "-GO_", go, "-plot-DE.png"), plot.go.de.line, width=10, height=10)
  
  go.genes.table <- rbind(go.genes.table, list(n,go,length(unique(counts.go$genes)),length(unique(counts.go.de$genes))))

  i<-i+1
}

# For GO terms that are significant with with all genes and DE genes
nested <- ((endo.all|dige.all|immu.all)/(endo.de|dige.de|immu.de)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides="collect") &
  theme(legend.position="none", plot.tag = element_text(size = 30), 
        axis.text = element_text(size = 13), axis.title = element_text(size = 15)) &
  guides(fill=guide_legend(nrow=2))
nested 
    
# For GO terms that are significant with with all genes and DE genes
nested2 <- ((endo.all.line|dige.all.line|immu.all.line)/(endo.de.line|dige.de.line|immu.de.line)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides="collect") &
  theme(legend.position="none", plot.tag = element_text(size = 30), 
        axis.text = element_text(size = 13), axis.title = element_text(size = 15)) &
  guides(fill=guide_legend(nrow=2))
nested2

ggsave(file=paste0(outputPrefix, "-plot-GOterm-expression.tif"), nested, device="tiff")
ggsave(file=paste0(outputPrefix, "-plot-GOterm-expression.png"), nested)

ggsave(file=paste0(outputPrefix, "-plot-GOterm-expression-line.tif"), nested2, device="tiff")
ggsave(file=paste0(outputPrefix, "-plot-GOterm-expression-line.png"), nested2)
