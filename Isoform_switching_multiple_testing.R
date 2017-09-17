####  1.2  #### 

# Data reading
setwd("~/HT/HW/HW3/data_for_hw3/data_for_part_1/")

df <- read.table("normalized_data.txt", h = F)
names(df) <- c("D1", "D2", "D3", "D4", "D5", "C1", "C2", "C3", "C4","C5")

# Fancier way of naming the columns
names(df) <- paste(rep(c("HIV", "Control"), each = 5), 1:5, sep = "-")

# We need to test the difference of means for two distributions of expression values for every gene
# in the dataset.
# As there is no evidence showing that the values do not follow a normal distribution, we do not 
# reject it, and thus we will assume that the mean of each sample is following a student-T distribu-
# tion. Therefore, we will use a student-T test to find whereas there are significant differences
# in gene expression.

my_pvals <- apply(df, 1, function (x) t.test(x[1:5], x[6:10])$p.value)


####  1.3  #### 

# Genes with a p-value less than 0.05
diff_genes <- sum(my_pvals < 0.05)                # 1911

## CONCEPTUAL MISTAKE ABOUT FP HERE ## 

# Expected number of FP with a threshold of 0.05 -> 0.05% of the t.tests are wrong
FP <- round(0.05 * length(my_pvals), 0)           # 1114

# Expected FP = threshold% used in the t.test -> 0.05% of the t. tests that reject the H0 are wrong 
FP1 <- round(0.05 * diff_genes, 0)                # 96


####  1.4  #### 

# Bonferroni multiple testing correction
my_bonferroni_pvals <- p.adjust(my_pvals, method = "bonferroni")

# Number of genes with p-value < 0.2 with Bonferroni correction
sum(my_bonferroni_pvals < 0.2) # 0

# BH multiple testing correction
my_BH_pvals <- p.adjust(my_pvals, method = "BH")

# Number of genes with p-value < 0.2 with BH correction
sum(my_BH_pvals < 0.2) # 12

# Expected FP = threshold% used in the t.test -> 0.2% of the t. tests are expected to be FP
FP2 <- round(0.2 * sum(my_BH_pvals < 0.2), 0) # 2

# The number of significant bonferroni adjusted genes was 0, while the FDR method yielded 12 
# significant genes. The expected number of false positive genes under FDR was 2, or 20% of the 
# significant FDR adjusted genes (p-value <0.2). 


####  1.5  #### 

# Addition of means
df$D_mean <- rowMeans(df[, c(1:5)])
df$C_mean <- rowMeans(df[, c(6:10)])

# Calculate the log2 foldchange for each gene using this formula:

            # foldchange = log2(mean(hiv)) - log2(mean(control)) # 

my_foldchanges <- log2(df[, 11]) - log2(df[, 12])


####  1.6  #### 

# Report the fold changes for the genes with a FDR < 0.2 (using the BH method)
# Are there most up or down-regulated genes in the HIV subset?
genes_FDR_lower_0.2 <- log2(df[which(my_BH_pvals < 0.2), 11]) - 
                       log2(df[which(my_BH_pvals < 0.2), 12])

genes_FDR_lower_0.2

# It is interesting to point out that the 12 genes with a significative fold change are upregulated 
# in the HIV patients, where the gene with the highest fold change has a 0.6 fold change, meaning it
# has a 50% higher expression level than the same gene in the control condition (because 2 ** 0.6 =
# = 1.51). 
# However, some of the genes with the lowest fold changes might not be biologically significant


library("dplyr")

####  3.1  #### 

# Data reading
setwd("~/HT/HW/HW3/data_for_hw3/data_for_part_3/")

gene <- read.table("cuffdiff_gene_differential_expression.txt", h = T)
transcript <- read.table("cuffdiff_transcript_differential_expression.txt", h = T)
splicing <- read.table("cuffdiff_splicing_differential_expression.txt", h = T)

head(gene); head(transcript); head(splicing)

names(transcript)[1] <- "transcript_id"


####  3.2  #### 

# Make two new dataframes from gene and transcript containing all rows with expression in at least 
# one of the conditions, that is, with an expression value < 0 in column 8 or 9 (value 1 or value2)
gene_expressed <- gene[gene$value_1 | gene$value_2 != 0, ]                                
transcript_expressed <- transcript[which(transcript$value_1 | transcript$value_2 != 0), ] 


####  3.3  #### 

# Number of genes and transcripts that were expressed and how many were significantly differentially
# expressed between conditions (column 14 or significant = "yes")
nrow(gene_expressed)        # 26
nrow(transcript_expressed)  # 98

# dplyr incorporation
gene_expressed %>% filter(significant == "yes") %>% nrow                        # 12
nrow(transcript_expressed[which(transcript_expressed$significant == "yes"), ])  # 11


####  3.4  #### 

# Make two reduced DF from gene_expressed and transcript_expressed and merge them based on the g. IDs
names(gene_expressed); names(transcript_expressed)

gene_reduced <- gene_expressed[, c(2, 3, 8, 9)]                     
transcript_reduced <- transcript_expressed[, c(1, 2, 8, 9)]

genes_and_transcripts <- merge(gene_reduced, transcript_reduced, by = "gene_id", 
                              suffixes = c(".GENE", ".TRANSCRIPT"))

dim(genes_and_transcripts)    # 98 rows, 7 columns


####  3.5  #### 

# Calculate the IF for each transcript in both conditions (IF = isoform_exp / gene_exp) and the dIF
# (IF2 - IF1)

genes_and_transcripts$IF1 <- genes_and_transcripts$value_1.TRANSCRIPT / 
                             genes_and_transcripts$value_1.GENE

genes_and_transcripts$IF2 <- genes_and_transcripts$value_2.TRANSCRIPT / 
  genes_and_transcripts$value_2.GENE

genes_and_transcripts$dIF <- genes_and_transcripts$IF2 - genes_and_transcripts$IF1

# Sanity check  
head(genes_and_transcripts[, c(3,4,6,7,8,9,10)])

# If the expression value of a certain gene and transcript in a specific condition is 0, when we try to
# obtain the IF dividing one by the other, we obtain a non-sense mathematical result due to the fact
# that we are dividing 0 by 0, and hence, R traduces this non-sense result into a NaN. As we are 
# considering only the isoform switching information, these rows can be safely removed, as there is no 
# transition to the second isoform (because this one does not appear at all in the experiment) and no 
# comparison of isoform fractions across the conditions can be performed


# Overwrite the IF1 and IF2 columns of the dataframe, eliminating rows with missing values
genes_and_transcripts <- genes_and_transcripts[complete.cases(genes_and_transcripts),]

# Recalculate the dIF with the corrected IF1 and IF2 values and overwrite
genes_and_transcripts$dIF <- genes_and_transcripts$IF2 - genes_and_transcripts$IF1

# Sanity check
head(genes_and_transcripts[, c(3,4,6,7,8,9,10)])


####  3.6  #### 

# Calculate the mean and median dIF value and discuss
mean(genes_and_transcripts$dIF)
median(genes_and_transcripts$dIF) 

# By definition the average dif is 0 and deviations from this number are due to numerical errors. 
# The median should also be equal to 0 if all genes had an even number of isoforms, which is hardly 
# the case. However, this will generate small deviations from 0.

# Histogram of the distribution of dIF values
hist(genes_and_transcripts$dIF, breaks = 40, main = "Distribution of dIF values", xlab = "dIF")

# We obtain a mean and a median very close to 0. This agrees with an expected mean dIF of 0 since 
# adding up all the dIF values for the different isoforms of the same gene is 0 by definition
# (although the median does not have to be 0). However, we are not obtaining a 0 mean due to numerical
# errors. In conclusion, the mean and the median do not provide any summary information of the data.


####  3.7  #### 

# Subset the merged dataframe to only contain genes with dIF > +0.25 and dIF < -0.25 
# Then, add the p-value from the splicing data frame using the match() function
genes_and_transcripts_potential_switch <- genes_and_transcripts[which(genes_and_transcripts$dIF 
                                          > 0.25 | genes_and_transcripts$dIF < -0.25), ] 
# p-values from splicing dataframe which same gene_id as the ones in genes_and_transcripts_potential_
# _switch dataframe
genes_and_transcripts_potential_switch$pvals <- splicing[match(genes_and_transcripts_potential_switch$
                                                                 gene_id, splicing$gene_id), 12]
# Sanity check
head(genes_and_transcripts_potential_switch, 3)


####  3.8  #### 

# TRANSCRIPT ID of the gene with the lowest p-value in the genes_and_transcripts_potential_switch DF
genes_and_transcripts_potential_switch[which.min(genes_and_transcripts_potential_switch$pvals), 5]
 # TCONS_00000021 #

# GENE NAME of the gene with the lowest p-value in the genes_and_transcripts_potential_switch DF
genes_and_transcripts_potential_switch[which.min(genes_and_transcripts_potential_switch$pvals), 2]
 # uc001aeo.3 #

# dIF VALUE of the gene with the lowest p-value in the genes_and_transcripts_potential_switch DF
genes_and_transcripts_potential_switch[which.min(genes_and_transcripts_potential_switch$pvals), 10]
 # -0.3426094 #

# P-VALUE of the gene with the lowest p-value in the genes_and_transcripts_potential_switch DF
genes_and_transcripts_potential_switch[which.min(genes_and_transcripts_potential_switch$pvals), 11]
 # 0.00075 #


