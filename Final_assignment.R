rm(list = ls())

# External packages loading (these will be the only ones used along the whole exam)
library("ggplot2")
library("reshape2")

## EXERCISE 1 ##

# Data reading
setwd("~/HT/Exam/")
peaks <- read.table("merged_sorted_txnChIP.bed")

# Re-name the columns
names(peaks) <- c("CHR", "START", "END", "PEAKS")

# Plot of the ChIP peaks distribution in each merged region
ggplot(data = peaks, aes(x = peaks$PEAKS)) +
  geom_histogram(binwidth = 1, boundary = -0.5, col = "black",
                 fill = "red", alpha = I(0.8)) +
  labs(title = "Distribution of merged regions per number of ChIP peaks") + 
  labs(x = "Number of ChIP peaks in log2 scale", y = "Number of merged regions") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(breaks = seq(0, 350000, 25000))

# Sort merged regions by number of ChIP peaks in descending order
sorted_peaks <- peaks[order(peaks$PEAKS, decreasing = T), ]

# Extract the top 10 merged regions formed by the highest number of peaks
peaks_top_10 <- head(sorted_peaks, 10)

# Add a 4th column of zeros (it should be the definition of each .bed line, but it is not included in 
# this experiment, its only purpose is to adjust the score column -number of peaks in this case- to its 
# correct position at 5th column to keep the proper .bed format of the file)
peaks_top_10$EMPTY <- rep(0, nrow(peaks_top_10))
peaks_top_10 <- peaks_top_10[c("CHR", "START", "END", "EMPTY", "PEAKS")]

# Export top 10 regions as .bed file
write.table(peaks_top_10, file = "peaks_top_10.bed", row.names = F, col.names = F, quote = F)


## EXERCISE 3 ##

# Data loading
setwd("~/HT/Exam/data/part3/")
data3 <- load("part3.Rdata")

# Number of different isoforms in countDF
length(unique(rownames(countDF)))

## Calculation of RPKM values ##

# Vectorized analysis to calculate RPKM values: countDF = number of reads mapping each isoform
# annotationDF$length = length of each isoform; sum(countDF) = total library size
RPKM_vals <- scale(countDF, center = F, scale = colSums(countDF) / 1e6) / (annotationDF$length / 1e3)

# Calculate the RPKM value of each sample
colMeans(RPKM_vals)

# Melting of the DF (association of isoforms with a Sample and a RPKM_value)
melted_logRpkm <- melt(data = logRpkmDF, variable.name = "Sample", value.name = "RPKM_value")

# Binary classification of the tool used: Kallisto or Salmon
melted_logRpkm$Tool <- ifelse(melted_logRpkm$Sample == "K_WT1_RPKM" | melted_logRpkm$Sample == "K_WT2_RPKM" |
                                 melted_logRpkm$Sample == "K_WT3_RPKM", "Kallisto", "Salmon")

# Classification of the replicate number
index_replicate <- c("K_WT1_RPKM", "K_WT2_RPKM", "K_WT3_RPKM", "S_WT1_RPKM", "S_WT2_RPKM", "S_WT3_RPKM")
values_replicate <- c("WT1", "WT2", "WT3", "WT1", "WT2", "WT3")
melted_logRpkm$Replicate <- values_replicate[match(melted_logRpkm$Sample, index_replicate)]

# Histograms of the logRPKM values, sepparating tool by colour and biological replicates by facets
ggplot() + geom_histogram(data = melted_logRpkm, alpha = I(0.6), bins = 25,
                          aes(x = RPKM_value, col = Tool, fill = Tool)) +
  facet_grid(facets = ~ Replicate) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(breaks = seq(0, 20000, 2500)) +
  labs(title = "Distribution of log10 (RPKM values) obtained in the three replicates",
       x = "log10 (RPKM values)", y = "Count") 

# Calculate all pairwise replicate Pearson correlation coefficients of the logRPKM expression values
# obtained with each tool
cor(logRpkmDF$K_WT1_RPKM, logRpkmDF$K_WT2_RPKM)
cor(logRpkmDF$K_WT1_RPKM, logRpkmDF$K_WT3_RPKM)
cor(logRpkmDF$K_WT2_RPKM, logRpkmDF$K_WT3_RPKM)

cor(logRpkmDF$S_WT1_RPKM, logRpkmDF$S_WT2_RPKM)
cor(logRpkmDF$S_WT1_RPKM, logRpkmDF$S_WT3_RPKM)
cor(logRpkmDF$S_WT2_RPKM, logRpkmDF$S_WT3_RPKM)


# Generic function to calculate the CV
CV_calc <- function(x) {
  cv <- var(as.numeric(x)) / mean(as.numeric(x))
  return(cv)
}

# Example: isoform TCONS_00034685 logRPKM values (for the three replicates analyzed with Kallisto) 
# CV calculation
CV_calc(logRpkmDF[1, c(1, 2, 3)])

# Use of apply() and the customized function to calculate the CVs based on logRpkm values for the 3 
# replicates with Kallisto and Salmon separately
CV_kallisto <- apply(logRpkmDF[, c(1, 2, 3)], 1, function(x) CV_calc(x))
CV_salmon <- apply(logRpkmDF[, c(4, 5, 6)], 1, function(x) CV_calc(x))

# Merge both CV sets
merged_CV <- data.frame(CV_kallisto, CV_salmon)

# Re-name columns and melt merged DF
colnames(merged_CV) <- c("Kallisto", "Salmon")
melted_CV <- melt(merged_CV, variable.name = "Tool", value.name = "CV_value")

ggplot(melted_CV, aes(x = CV_value, fill = Tool)) + geom_density(alpha = 0.3) + 
  labs(title = "Distribution of CV values in log10 scale", x = "CV values") +
  scale_x_continuous(trans = "log10")


## EXERCISE 4 ##

set.seed(2017)

# Data reading
setwd("~/HT/Exam/data/part4/")
data4 <- read.table("ExpressionMatrix.tab", h = T)
data4_design <- read.table("StudyDesign.tab", h = T)

# Number of samples (columns) and genes (rows)
length(colnames(data4))
length(rownames(data4))

# PCA on the scaled samples and amount of variance contained in the first 5 PCs. 
# To have consistent dimensions between the expression matrix and the information gathered in the study 
# design, a transposition of the "data4" matrix is required
pca_data4_scaled <- prcomp(t(data4), scale = T)
summary(pca_data4_scaled)$importance[2, 1:5]

# Proportion of variance extraction
prop_var <- summary(pca_data4_scaled)$importance[2, ]

# Save data and design information on sample groups in a single DF
plot_data4_scaled <- data.frame(pca_data4_scaled$x, data4_design)

# Dimensionality reduction plot differentiating by KD efficiency and KD target
qplot(data = plot_data4_scaled, x = PC1, y = PC2, color = KnockDownEfficiency, shape = KnockDownTarget) +
  labs(x = paste("PC1", round(prop_var[1] * 100, digits = 2), "%", sep = " "),
       y = paste("PC2", round(prop_var[2] * 100, digits = 2), "%", sep = " "))

# K-means clustering of the data with 3 clusters and 10 random starting points
groups <- kmeans(t(data4), centers = 3, nstart = 10)

# Merged DF with PCA output and cluster information, changing the name of clusters column
data4_clustered <- cbind(plot_data4_scaled, as.character(groups$cluster))
colnames(data4_clustered)[80] <- "cluster"

# Plot of PC1 vs PC2, differentiating clusters by color
ggplot(data = data4_clustered, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(x = paste("PC1", round(prop_var[1] * 100, digits = 2), "%", sep = " "),
       y = paste("PC2", round(prop_var[2] * 100, digits = 2), "%", sep = " "))

# Plot of PC1 vs PC2, differentiating by user
ggplot(data = data4_clustered, aes(x = PC1, y = PC2, color = PostDoc)) +
  geom_point() +
  labs(x = paste("PC1", round(prop_var[1] * 100, digits = 2), "%", sep = " "),
       y = paste("PC2", round(prop_var[2] * 100, digits = 2), "%", sep = " "))


## EXERCISE 5 ##

# Data reading and first lines
setwd("~/HT/Exam/data/part5/")
data5 <- read.table("DifferentialExpression.tab", h = T)
head(data5, 5)

# DF massage to perform the plot

# TF_A subset merging
data5_A <- data.frame(data5[, c(3, 4, 5)], rep("Transcription Factor A KD", 10000))
colnames(data5_A) <- c("logFC", "pvalue", "FDR", "TF")

# TF_B subset merging
data5_B <- data.frame(data5[, c(6, 7, 8)], rep("Transcription Factor B KD", 10000))
colnames(data5_B) <- c("logFC", "pvalue", "FDR", "TF")

# Duplicating the mean expression data
data5_exp <- rbind(data5[, c(1, 2)], data5[, c(1, 2)])

# Final DF merging (for plotting purposes)
merged_DF <- cbind(data5_exp, rbind(data5_A, data5_B))

# MA plots for each KD
ggplot(data = merged_DF, aes(x = meanExpression, y = logFC)) + 
  geom_point(alpha = I(0.4), size = 1.2) +
  facet_grid(facets = ~ TF) + 
  labs(title = "MA plots for TFs A and B KDs", x = "A (Mean expression values)", y = "M (log(FC))") +
  geom_smooth(color = "red", size = 0.5) +
  geom_hline(aes(yintercept = mean(logFC, na.rm = T)),
             color = "blue", size = 0.5)

# Number of genes significantly DE (after correcting for multiple testing correction)
sum(data5$A.FDR < 0.05)
sum(data5$B.FDR < 0.05)

# Total number of genes significantly DE 
sum(data5$A.pvalue < 0.05)
sum(data5$B.pvalue < 0.05)

# Classify genes in a new column by their significant DE in none, one or both KD
data5$type <- ifelse(data5$A.FDR < 0.05 & data5$B.FDR < 0.05, "KD_A & KD_B",
                ifelse(data5$A.FDR < 0.05 & data5$B.FDR >= 0.05, "KD_A",
                  ifelse(data5$A.FDR >= 0.05 & data5$B.FDR < 0.05, "KD_B", "NONE")))

# Plot of logFCs A VS logFCs B, differentiating points significantly DE after multiple testing 
# correction in KDA, KDB, both, or none of them
ggplot(data = data5, aes(x = B.logFC, y = A.logFC, color = type)) + 
  geom_point(alpha = I(0.5), size = 1.4) +
  labs(title = "logFCs A VS logFCs B", x = "logFCs B", y = "logFCs A")
  
# Pearson correlation coefficient between the logFCs of the two KDs
cor(data5$A.logFC, data5$B.logFC, method = "pearson")

# Check normality of the logFCs distributions
ggplot() + geom_histogram(data = merged_DF, alpha = I(0.6), bins = 40, col = "black",
                          aes(x = logFC, fill = TF)) +
  labs(title = "logFCs A and B KD subsets distributions") + labs(x = "logFC", y = "Count") +
  scale_x_continuous(limits = c(-4, 4))

# Test whether this correlation is significant
t.test(data5$A.logFC, data5$B.logFC, alternative = "two.sided")

# Number of genes DE and non-DE in both subsets
DE_A <- sum(data5$A.FDR < 0.05)
NDE_A <- sum(data5$A.FDR >= 0.05)
DE_B <- sum(data5$B.FDR < 0.05)
NDE_B <- sum(data5$B.FDR >= 0.05)

# Build contingency table and perform Fisher exact test
contingency_table <- matrix(c(DE_A, NDE_A, DE_B, NDE_B), nrow = 2, byrow = T)
fisher.test(contingency_table)


