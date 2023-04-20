# RNA-Seq-Pipieline
I share the code used in a Linux server and then RStudio to perform a standard RNA-seq experiment.

--- Linux Server ---
# Step 1: FastQC to perform quality control

fastqc *gz -q --outdir = ~/project/fastqc_results

# Step 2: MultiQC to summarise FastQC output

multiqc . -n MultiQC_FastQC -o ~/project/multiqc_results

# Step 3: STAR alignment to reference genome

reference_genome_directory = "/media/newdrive/data/Reference_genomes/Human/UCSC/STAR_UCSC_refseq/"

for i in *1.fastq.gz
do
    STAR --genomeDir $reference_genome_directory 
    --runThreadN 12 
    --readFilesIn $i ${i%_1.fastq.gz}_2.fastq.gz 
    --readFilesCommand zcat 
    --outFileNamePrefix STAR_${i%_1.fastq.gz} 
    --outSAMtype BAM SortedByCoordinate 
    --outSAMunmapped Within 
    --outSAMattributes Standard
done

# Step 4: Post-allignment QC with bam_stat.py, SAMtools Flagstat, genebody_coverage.py, read_distribution.py

for i in *bam 
do
    bam_stat.py -i $i > ~/project/bam_stat_results/bam_stat_${i%Aligned.sortedByCoord.out.bam}.txt 
    samtools flagstat $i > ~/project/flagstat_results/flagstat_${i%Aligned.sortedByCoord.out.bam}.txt
    geneBody_coverage.py -r /media/newdrive/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12 
    -i $i
    -o ${i%Aligned.sortedByCoord.out.bam}
    read_distribution.py -r /media/newdrive/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12
    -i $i > ~/project/star_results/${i%Aligned.sortedByCoord.out.bam}.read_list.txt
done

# Step 5: Create index files and sort with samtools

for i in *bam
do 
    samtools index ${i} 
    samtools sort $i -i sorted_${i%Aligned_SortedByCoord.out.bam} -T ~/project/star_results
done

# Step 6: Use MultiQC to visualise QC step 4

multiqc ~/project/qc_results -n MultQC_PostQC

# Step 7 (Optional): Determine strandedness if you don't already know it with infer_experiment.py

for i in *bam
do
    infer_experiment.py -i $i -r /media/newdrive/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12 >> ~/project/read_quantification/strandedness.txt
done

# Step 8: Quantify number of times every gene was expressed with Feature Counts

featureCounts -T 12 \
-s 1 \
-p \
-a  /media/newdrive/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.gtf \
-o ~/project/read_quantification/featureCountsStranded.txt ~/project/star_results/sorted*

# Step 9: Visualise results of Feature Counts with MultiQC

multiqc ~/project/read_quantification/feature* -n MultiQc_FeatureCounts -o ~/project/multiqc_results

# If you do not wish to use alignment with STAR, use pseudo-alignment with Salmon

for i in *_1.fastq.gz
do
    salmon quant -i ~/project/salmon_test/index -l A \
    -1 $i \
    -2 ${i%_1.fastq.gz}_2.fastq.gz \
    -p 8 \
    --validateMappings \
    -o ~/project/salmon_test/${i%_1.fastq.gz}
done

# Whichever method was selected, a quantification file will be the final output of either

--- RStudio ---
# DGE with DESEQ2 

setwd("C:/Users/lukar/OneDrive/Documents/Project/R")
install.packages("BiocManager")
install.packages("NMF")
install.packages("dplyr")
BiocManager::install("DESeq2")
library (DESeq2)
library (magrittr)
library(NMF)
library (grDevices)

# Read in feature counts results
featureCountsData <- read.table ("featureCounts.txt", header = T)

# Set the Gene ID column as the row names
row.names (featureCountsData) <- featureCountsData$Geneid

# Remove rows with no reads
featureCountsData <- featureCountsData[ , -c (1:6) ]

# For loop to create names of the samples
sampleNames <- c()
nameAb <- c()
nameIgg <- c()

for (i in 1:4){
    nameAb <- paste ("WT_Ab_", i, sep = "")
    sampleNames <- c(sampleNames, nameAb)
    nameIgg <- paste ("WT_IgG_", i, sep = "")
    sampleNames <- c(sampleNames, nameIgg)
}
sampleNames <- sort (sampleNames)
names (featureCountsData) <- sampleNames

# Create a new data frame with sample information
sampleConditions <- data.frame ( condition = gsub( "_[0 -9]+","" , 
                                 names(featureCountsData)), 
                                 row.names = names(featureCountsData))


deseqRawMatrix <- DESeqDataSetFromMatrix (countData = featureCountsData,
                                            colData = sampleConditions,
                                            design = ~ condition)
# Check data
colData (deseqRawMatrix) %>% head
assay (deseqRawMatrix, "counts") %>% head
rowData (deseqRawMatrix) %>% head
counts (deseqRawMatrix) %>% str

# Remove rows with no data
deseqRawMatrix <- deseqRawMatrix[ rowSums(counts(deseqRawMatrix)) > 0, ]

# These two should be the same
colSums (counts(deseqRawMatrix))
colSums (featureCountsData)

# Normalization and transformation
normalizedMatrix <- estimateSizeFactors (deseqRawMatrix)
sizeFactors (deseqRawMatrix)
colData (deseqRawMatrix) %>% head # now includes size factors

normalizedCounts <- counts (normalizedMatrix, normalized = T)
logNormalizedCounts <- log2 (normalizedCounts + 1)
plot(logNormalizedCounts[,1:2], cex = 1, main = "Normalized and Log2 Transformed")
plot(normalizedCountsRlog[,1:2], cex = 1, main = "Regularised Log-Transformed (Rlog)")

# Homoskedaiety
deseqRlog <- rlog(normalizedMatrix, blind = T)
normalizedCountsRlog <- assay(deseqRlog)

# Hierarchical clustering
distanceRlog <- as.dist(1- cor(normalizedCountsRlog, method = "spearman"))
plot(hclust(distanceRlog), labels = colnames(normalizedCountsRlog))

# PCA
pcaPlot <- plotPCA(deseqRlog)
pcaPlot

# 6.3 Running DGE tools 
colData(deseqRawMatrix)$condition <- relevel(colData(deseqRawMatrix)$condition, "WT_IgG")

deseqRawMatrix <- DESeq(deseqRawMatrix)
deseqDgeResults <- results(deseqRawMatrix, independentFiltering = T, alpha = 0.05)

# DGE checks
hist(deseqDgeResults$pvalue)
plotMA(deseqDgeResults[deseqDgeResults$baseMean > 15,], alpha = 0.05, ylim = c(-4,4))
summary(deseqDgeResults)
head (deseqDgeResults)
table (deseqDgeResults$padj < 0.05)
rownames (subset (deseqDgeResults, padj < 0.05))

# DE genes

sortedDeseqDgeResults <- deseqDgeResults[order(deseqDgeResults$padj), ]
degenesPadjDeseq <- subset(sortedDeseqDgeResults, padj < 0.05)
degenesPadjLogDeseq <- subset(sortedDeseqDgeResults,
                                       padj < 0.05 &
                                       abs(log2FoldChange) >= 1)
heatmapGenes <- logNormalizedCounts[rownames(degenesPadjLogDeseq), ]
top20GenesDeseq <- head(heatmapGenes,20)

# Heatmaps
aheatmap(heatmapGenes, Rowv = NA, Colv = NA)
aheatmap(top20GenesDeseq, Rowv = NA, Colv = NA)
aheatmap(heatmapGenes, Rowv = T, Colv = T, distfun = "euclidean", hclustfun = "average")
aheatmap(top20GenesDeseq, Rowv = T, Colv = T, distfun = "euclidean", hclustfun = "average")
aheatmap(heatmapGenes, Rowv = T, Colv = T, distfun = "euclidean", hclustfun = "average", scale = "row")
aheatmap(top20GenesDeseq, Rowv = T, Colv = T, distfun = "euclidean", hclustfun = "average", scale = "row")

# Volcano Plot
results_tibble <- deseqDgeResults %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble()
results_tibble <- results_tibble %>% mutate(threshold = padj <= 0.05)
results_tibble <- results_tibble %>% arrange(padj) %>% mutate(genelabels = "")
results_tibble$genelabels[1:10] <- results_tibble$gene[1:10]

ggplot(results_tibble, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = threshold)) +
    geom_text_repel(aes(label = genelabels), max.overlaps = 20) +
    ggtitle("Volcano Plot of DE Genes")
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) 

# Compare gene expression

plotCounts(dds = deseqRawMatrix,
           gene = "SIGLEC11",
           normalized = T, transform = F)

DeseqCommonGenes <- c()
for (i in 1:length(commonDeGenes)){
    if (commonDeGenes[i] %in% rownames(heatmapGenes)){
        DeseqCommonGenes <- c(DeseqCommonGenes, commonDeGenes[i])
    }
}

DEGdf <- heatmapGenes [ rownames(heatmapGenes) %in% DeseqCommonGenes, ] 
aheatmap(head(DEGdf, 20), Rowv = T, Colv = T, distfun = "euclidean", hclustfun = "average", scale = "row")

# DGE with EdgeR

library(edgeR)
setwd("C:/Users/lukar/OneDrive/Documents/Project/R")

sampleInfo <- factor (c( rep("WT_IgG",4), rep("WT_Ab", 4)))
sampleInfo <- relevel (sampleInfo, ref = "WT_IgG")

listEdgeR <- DGEList (counts = featureCountsData,
                      group = sampleInfo)
head (listEdgeR$counts)
listEdgeR$samples

cpm <- log2 (rowSums(cpm(listEdgeR)))
hist (cpm)
summary (cpm)

keep <- rowSums (cpm(listEdgeR) >= 1) >= 1 
listEdgeR <- listEdgeR[keep,]
listEdgeR$samples$lib.size <- colSums (listEdgeR$counts)
head (listEdgeR$samples)

listEdgeR <- calcNormFactors (listEdgeR, method = "TMM")

designEdgeR <- model.matrix (~sampleConditions$condition)
listEdgeR <- estimateDisp (listEdgeR, designEdgeR)
fitEdgeR <- glmFit (listEdgeR, designEdgeR)
lrtEdgeR <- glmLRT (fitEdgeR)
resultsEdgeR <- topTags (lrtEdgeR, n = Inf, sort.by = "PValue", adjust.method = "BH")

#DGE with Limma-voom

rownames(designEdgeR) <- colnames(listEdgeR)
voomTransformed <- voom(listEdgeR, designEdgeR, plot = F)
voomFitted <- lmFit(voomTransformed, design = designEdgeR)
voomFitted <- eBayes(voomFitted)
resultsLimma <- topTable(voomFitted, coef="sampleConditions$conditionWT_IgG", 
                              number = Inf, 
                              adjust.method = "BH",
                              sort.by = "logFC")
resultsLimmaSorted <- resultsLimma[order(resultsLimma$adj.P.Val),]

#DGE of Salmon output with DESEQ2

setwd("C:/Users/lukar/OneDrive/Documents/Project/Salmon")
library (tximport)
library (DESeq2)

txgene <- read.csv("tx2gene.csv", sep = "\t", header = F)
colnames(txgene) <- c("Transcript_id", "Gene_id")

filepaths <- c()
for (i in sampleNames){
    filepaths <- c(filepaths,file.path(getwd(), i, paste(i,"quant.sf", sep = "_")))
}

txi <- tximport(files = filepaths,
                type = "salmon",
                tx2gene = txgene,
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = T)

ddsSalmon <- DESeqDataSetFromTximport(txi, sampleConditions, design = ~condition)
ddsSalmon <- ddsSalmon[rowSums(counts(ddsSalmon))>1,]
ddsSalmon <- DESeq(ddsSalmon)
resultsSalmon <- results (ddsSalmon, independentFiltering = T, alpha = 0.05)
degenesSalmon <- resultsSalmon[
    which(resultsSalmon$padj < 0.05 & abs(resultsSalmon$log2FoldChange) < 1), ]

# DE genes
library(NMF)
dgeSalmonSorted <- resultsSalmon[order(resultsSalmon$padj),]
degenesPadjSalmon <- rownames (subset(dgeSalmonSorted, padj < 0.05))
degenesPadjLogSalmon <- rownames (subset(dgeSalmonSorted,
                                          padj < 0.05 &
                                              abs(log2FoldChange) >= 1))

normalizedCountsSalmon <- estimateSizeFactors(ddsSalmon)
salmonNormalized <- counts (normalizedCountsSalmon, normalized = T)
salmonLogNormalized <- log2 (salmonNormalized + 1)

library(org.Hs.eg.db)
data <- as.vector (degenesPadjSalmon)
annots <- select (org.Hs.eg.db, keys=data, columns="SYMBOL", keytype="ENSEMBL")
annots <- na.omit(annots)
sigSalmonResults <- subset (resultsSalmon, padj <= 0.05)

for (i in annots$ENSEMBL){
    for (j in 1:nrow(sigSalmonResults)){
        if (i == rownames(sigSalmonResults)[j]){
            rownames(sigSalmonResults)[j] <- annots$SYMBOL[annots$ENSEMBL == i]
        }
    }
}

# Heatmaps
heatmapGenesSalmon <- salmonLogNormalized[degenesPadjLogSalmon, ]
top20GenesHeatmap <- head (heatmapGenesSalmon, 20)
aheatmap(heatmapGenesSalmon, Rowv = NA, Colv = NA)
aheatmap(top20GenesHeatmap, Rowv = NA, Colv = NA)

-----------------------------------------------
# Compare the results of all DGE tools

library(gplots)
DE_list <- list(edger = rownames(subset(resultsEdgeR$table, FDR<=0.05)), 
                deseq2 = rownames(subset(deseqDgeResults, padj <= 0.05)),
                limma = rownames(subset(resultsLimma, adj.P.Val <= 0.05)),
                salmon = rownames(sigSalmonResults))
gplots::venn(DE_list)

commonDeGenes <- Reduce(intersect, DE_list)
length(commonDeGenes)
head(commonDeGenes, 100)

library(UpSetR)
DEGenes <- fromList(DE_list)
upset(DEGenes, order.by = "freq")

# Gene Set Enrichment Analysis

setwd("C:/Users/lukar/OneDrive/Documents/Project/R")

testGSEA <- deseqDgeResults[order(-deseqDgeResults$log2FoldChange),]
testGSEA <- data.frame(rownames(testGSEA), testGSEA[,"stat"])
testGSEA <- sort (testGSEA$stat, decreasing = T)
write.table(testGSEA, "log2FoldChangeRanking.rnk", sep = "\t",
            col.names = FALSE, row.names = FALSE)
