setwd("~/HT/HW/HW1/")
d <- read.table("gene_lengths_v2.txt", header = T)
head(d)

## EXERCISE 1 ##

# Vector of breakpoints so that each possible exon number will fall into a bin
# The number of bins is always one less than the number of breakpoints
my_hist <- hist(d$exon_count, breaks = c(min(d$exon_count - 1):max(d$exon_count)), 
     xaxt = "n", col = "red", xlab = "Number of exons", ylab = "Number of genes",
     main = "Number of exons per gene")

# More straight-forward option
my_hist <- hist(d$exon_count, breaks = c(0:max(d$exon_count)), 
                xaxt = "n", col = "red", xlab = "Number of exons", ylab = "Number of genes",
                main = "Number of exons per gene")

# x-axis
axis(1, at = seq(0, 150, 5), labels=seq(0, 150, 5))

# The most common number of exons seems to be between 3 and 5

# Abline showing the threshold (optional)
abline(h = max(my_hist$counts))

# Reduced datasets: only genes with 20 or less exons, to clearly distinguish the number of exons MODE
d_red <- d[which(d$exon_count <= 20), ]
head(d_red)

# Reduced histogram
my_red_hist <- hist(d_red$exon_count, breaks = c(0:20), xaxt = "n", col = "blue", 
                    xlab = "Number of exons", ylab = "Number of genes", 
                    main = "Number of exons per gene (from 0 to 20)")

# x-axis (reduced)
axis(1, at = seq(0, 20), labels = seq(0, 20))

# Fancier plots
library(ggplot2)

# Simple one
ggplot() + geom_histogram(data = d, aes(x = exon_count), breaks = seq(0, 150, by = 1), col = "black",
                          fill = "red", alpha = I(0.8)) +
  labs(title = "Number of exons per gene") + labs(x = "Number of exons", y = "Number of genes")

# With gradient 
ggplot(data = d, aes(d$exon_count)) + 
  geom_histogram(breaks = seq(0, 150, by = 1), aes(fill = ..count..)) + 
  scale_fill_gradient("Count", low = "green", high = "red") + 
  labs(title = "Number of exons per gene") + labs(x = "Number of exons", y = "Number of genes") + 
  scale_x_continuous(breaks = seq(0, 150, 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Reduced with gradient, expanded x and y labels
ggplot(data = d_red, aes(d_red$exon_count)) + 
  geom_histogram(breaks = seq(0, 20), aes(fill = ..count..), col = "black") + 
  scale_fill_gradient("Count", low = "green", high = "red") + 
  labs(title = "Number of exons per gene") + labs(x = "Number of exons", y = "Number of genes") + 
  scale_x_continuous(breaks = seq(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1500, 100)) +
  
  # Optional abline showing the threshold
  geom_hline(yintercept = max(table(d_red$exon_count)), linetype = "dashed") +
  geom_text(mapping = aes(x = which.max(table(d$exon_count)), y = max(table(d$exon_count))),
            label = max(table(d$exon_count)), vjust = -.3)


## EXERCISE 2 ##

d$intron_length <- d$genome_length - d$mrna_length


## EXERCISE 3 ##

par(mfrow = c(2,2))

# Set maximum value of intron and exon lengths as max number of bins
my_max <- max(c(d$intron_length, d$mrna_length))
my_break_points <- seq(0, my_max, length.out = 51)

hist(d$mrna_length, col = "blue", breaks = my_break_points, xlim = c(0, my_max),
     ylim = c(0, 12000), xlab = "Length of exons", ylab = "Numbe of exons", 
     main = "Length of exons distribution")
hist(d$intron_length, col = "red", breaks = my_break_points, xlim = c(0, my_max),
     ylim = c(0, 12000), xlab = "Length of introns", ylab = "Numbe of introns", 
     main = "Length of introns distribution")

boxplot(d$mrna_length, col = "blue", ylim = c(0, 60000))
boxplot(d$intron_length, col = "red", ylim = c(0, 60000))

par(mfrow = c(1,1))

library("reshape2")

# DF with mRNA and intron lengths (using the MELT function)
df <- subset(x = melt(data.frame(mrna_length = d$mrna_length, intron_length = d$intron_length)),
              subset = value != 0)

# Lengths given in megabases
ggplot(df, aes(x = value / 10**6, fill = variable)) +
  geom_histogram(bins = 50, binwidth = 0.2, position = "dodge") +
  labs(title = "Length of exons and introns (in Mb)") + labs(x = "Length (Mb)", y = "Count") +
  scale_fill_discrete(name = "",
                      breaks = c("mrna_length", "intron_length"),
                      labels = c("Exon", "Intron")) +
  scale_y_continuous(breaks = seq (0, max(df$value), 2000))  

# Using logarithmic lengths to show the differences in length between introns and exons
ggplot(df, aes(x = log10(value), fill = variable)) +
  geom_histogram(bins = 50, binwidth = 0.2, position = "dodge") +
  labs(title = "log10 (Length of exons and introns)") + labs(x = "log10 (Length)", y = "Count") +
  scale_fill_discrete(name = "",
                      breaks = c("mrna_length", "intron_length"),
                      labels = c("Exon", "Intron")) +
  scale_y_continuous(breaks = seq (0, max(df$value), 500))

# Boxplots with distances in Mbases
ggplot(df, aes(y = value / 10**6, x = factor(1), fill = variable)) +
  geom_boxplot(lwd = 0.25, outlier.size = 0.5) +
  labs(title = "Length of exons and introns (in Mb)") + labs(x = "", y = "Length (Mb)") +
  theme(strip.text.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_discrete(name = "",
                      breaks = c("mrna_length", "intron_length"),
                      labels = c("Exon", "Intron")) +
  coord_flip()

# Using logarithmic lengths again
ggplot(df, aes(y = log10(value), x = factor(1), fill = variable)) +
  geom_boxplot(lwd = 0.25, outlier.size = 0.5) +
  labs(title = "log10(Length of exons and introns)") + labs(x = "", y = "log10(Length)") +
  theme(strip.text.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_discrete(name = "",
                      breaks = c("mrna_length", "intron_length"),
                      labels = c("Exon", "Intron")) +
  coord_flip()


## EXERCISE 4 ##

wilcox.test(d$intron_length, d$mrna_length, alternative = "two.sided")
median(d$intron_length); median(d$mrna_length)
mean(d$intron_length); mean(d$mrna_length)


## EXERCISE 5 ##

# Correlation indices
r1 <- round(cor(d$mrna_length, d$intron_length, method = "pearson"), digits = 2)
r2 <- round(cor(d$mrna_length, d$exon_count, method = "pearson"), digits = 2)

par(mfrow = c(1,2))

# Linear model 1
model_1 <- lm(d$intron_length / 10**6 ~ d$mrna_length)
plot(d$mrna_length, d$intron_length / 10**6, pch = ".", xlab = "mRNA length", 
     ylab = "intron length (Mb)")

# Trendline 1 
abline(model_1, col = "blue")

# Correlation coefficient r1
text(x = 30000, y = 1.5, labels = paste("r2 = ", r1, sep = ""))

# Linear model 2
model_2 <- lm(d$exon_count ~ d$mrna_length)
plot(d$mrna_length, d$exon_count, pch = ".", xlab = "mRNA length", ylab = "exon count")

# Trendline 2
abline(model_2, col = "blue")

# Correlation coefficient r2
text(x = 30000, y = 150, labels = paste("r2 = ", r2, sep = ""))


## EXERCISE 6 ##

print(d[which.max(d$mrna_length), c(1,2,4)], row.names = F)
print(d[which.max(d$mrna_length), c("name", "mrna_length", "exon_count")], row.names = F)


## EXERCISE 7 ##

count_genes <- function(my_vector = d$mrna_length, 
                        x1 = 0, 
                        x2 = max(my_vector))
  {
  my_count <- sum(my_vector > x1 & my_vector <= x2)
  return(my_count/length(my_vector))
}

count_genes()
count_genes(x1 = 10000)
count_genes(x1 = 1000, x2 = 10000)
count_genes(x1 = 100, x2 = 1000)
count_genes(x2 = 100)
  


