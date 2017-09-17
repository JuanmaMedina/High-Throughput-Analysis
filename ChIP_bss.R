# Data reading
setwd("~/HT/HW/HW2")
ERA_coverage <- read.table(file = "ERa.txt")
ERB_coverage <- read.table(file = "ERb.txt")

# Add a column with the site name <---
ERA_coverage$site <- "ERa" 	
ERB_coverage$site <- "ERb" 	

# Combine both DF
df <- rbind(ERA_coverage, ERB_coverage)

# Name the columns
colnames(df) <- c("chr", "depth", "coverage_size", "chr_size", "fraction", "site")

# Keep rows with coverage in the genome (depth = 1) 
df <- df[which(df$depth == 1), ]

# Remove rows with genomic information
df <- df[which(df$chr != "genome"), ]

library(ggplot2)

ggplot(df) + 
  geom_bar(aes(x = chr, y = fraction, fill = site), stat = "identity", position = "dodge")
       
ggplot(df, aes(x = chr, y = fraction)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = site))

# Reading the ERa and ERb overlapping-with-genome files
AtoB <- read.table(file = "AtoBoverlap.bed")
BtoA <- read.table(file = "BtoAoverlap.bed")

totalA <- nrow(AtoB)
totalB <- nrow(BtoA)

# Number of overlapping sites between ERa and ERb
overlapping <- sum(AtoB$V4)

library(VennDiagram)

# Venn diagram
draw.pairwise.venn(totalA, totalB, overlapping, category = c("ERa", "ERb"), fill = c("red", "blue")) 

