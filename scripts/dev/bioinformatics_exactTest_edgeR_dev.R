##
# Setup
##

# set the working directory
#setwd("/YOUR/FILE/PATH/")
setwd("/Users/bamflappy/Repos/wpPlayground/")

# install libraries, if necessary
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("ggplot2")
#install.packages("ghibli")
#install.packages("ggVennDiagram")

# import libraries
library(edgeR)
library(ggplot2)
library(ghibli)
library(ggVennDiagram)

# import gene count data
tribolium_counts <- read.csv("data/TriboliumCounts.csv", row.names="X")

#Add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("treat_4h",3), rep("cntrl_24h",3), rep("treat_24h",3)))

#Create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)

# change the graphical parameters
png("plots/dev/ghibliPalettes.png")
par(mfrow=c(9,3))

# view all available ghibli palettes
for(i in names(ghibli_palettes)) print(ghibli_palette(i))

# close the plot and return the display to the default graphical parameters
dev.off()

# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")

# view the selected color palette
png("plots/dev/ghibliPalette_ponyoMedium.png")
ghibli_colors
dev.off()

##
# Normalization
##

#Plot the library sizes before normalization and write to a png file
png("plots/dev/exactTest_tribolium_librarySizes.png")
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")
dev.off() 

#There is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Calculate normalized factors
list <- calcNormFactors(list)
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Write the normalized counts to a file
write.table(normList, file="data/tribolium_normalizedCounts.csv", sep=",", row.names=TRUE)

#View normalization factors
list$samples

# vector of shape numbers for the MDS plot
points <- c(0,1,15,16)

# vector of colors for the MDS plot
colors <- rep(c(ghibli_colors[3], ghibli_colors[6]), 2)

# MDS plot with distances approximating log2 fold changes
png("plots/dev/exactTest_tribolium_MDS.png")
# add extra space to right of plot area and change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(list, col=colors[group], pch=points[group])
# place the legend outside the right side of the plot
legend("topright", inset=c(-0.4,0), legend=levels(group), pch=points, col=colors)
dev.off()

#Calculate the log CPM of the gene count data
logcpm <- cpm(list, log=TRUE)

#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
png("plots/dev/exactTest_tribolium_hclust.png")
heatmap(logcpm)
dev.off()

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#View dispersion estimates and biological coefficient of variation
png("plots/dev/exactTest_tribolium_BCV.png")
plotBCV(list)
dev.off()


##
# Pairwise Tests
##

# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

## treat_4h vs cntrl_4h
#Perform an exact test for treat_4h vs cntrl_4h
tested_4h <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_4h))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
png("plots/dev/exactTest_tribolium_4h_DE.png")
plotMD(tested_4h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Create results table of DE genes
resultsTbl_4h <- topTags(tested_4h, n=nrow(tested_4h$table))$table

#Identify significantly DE genes
resultsTbl_4h$topDE <- "NA"
resultsTbl_4h$topDE[resultsTbl_4h$logFC > 1 & resultsTbl_4h$FDR < 0.05] <- "UP"
resultsTbl_4h$topDE[resultsTbl_4h$logFC < -1 & resultsTbl_4h$FDR < 0.05] <- "DOWN"

#Create volcano plot
png("plots/dev/exactTest_tribolium_4h_volcano.png")
ggplot(data=resultsTbl_4h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))
dev.off()

#Create a table of DE genes filtered by FDR
resultsTbl_4h.keep <- resultsTbl_4h$FDR <= 0.05
resultsTbl_4h_filtered <- resultsTbl_4h[resultsTbl_4h.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_4h_filtered, file="data/exactTest_tribolium_4h_filtered.csv", sep=",", row.names=TRUE)

## treat_24h vs cntrl_24h
#Perform an exact test for treat_24h vs cntrl_24h
tested_24h <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_24h))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
png("plots/dev/exactTest_tribolium_24h_DE.png")
plotMD(tested_24h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Create a table of DE genes filtered by FDR
resultsTbl_24h <- topTags(tested_24h, n=nrow(tested_24h$table))$table

#Identify significantly DE genes
resultsTbl_24h$topDE <- "NA"
resultsTbl_24h$topDE[resultsTbl_24h$logFC > 1 & resultsTbl_24h$FDR < 0.05] <- "UP"
resultsTbl_24h$topDE[resultsTbl_24h$logFC < -1 & resultsTbl_24h$FDR < 0.05] <- "DOWN"

#Create volcano plot
png("plots/dev/exactTest_tribolium_24h_volcano.png")
ggplot(data=resultsTbl_24h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))
dev.off()

#Create filtered results table of DE genes
resultsTbl_24h.keep <- resultsTbl_24h$FDR <= 0.05
resultsTbl_24h_filtered <- resultsTbl_24h[resultsTbl_24h.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_24h_filtered, file="data/exactTest_tribolium_24h_filtered.csv", sep=",", row.names=TRUE)

## treat_4h vs treat_24h
#Perform an exact test for treat_4h vs treat_24h
tested_treat <- exactTest(list, pair=c("treat_24h", "treat_4h"))

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_treat))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
png("plots/dev/exactTest_tribolium_treat_DE.png")
plotMD(tested_treat)
abline(h=c(-1, 1), col="blue")
dev.off()

#Create a table of DE genes filtered by FDR
resultsTbl_treat <- topTags(tested_treat, n=nrow(tested_treat$table))$table

#Identify significantly DE genes
resultsTbl_treat$topDE <- "NA"
resultsTbl_treat$topDE[resultsTbl_treat$logFC > 1 & resultsTbl_treat$FDR < 0.05] <- "UP"
resultsTbl_treat$topDE[resultsTbl_treat$logFC < -1 & resultsTbl_treat$FDR < 0.05] <- "DOWN"

#Create volcano plot
png("plots/dev/exactTest_tribolium_treat_volcano.png")
ggplot(data=resultsTbl_treat, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))
dev.off()

#Create filtered results table of DE genes
resultsTbl_treat.keep <- resultsTbl_treat$FDR <= 0.05
resultsTbl_treat_filtered <- resultsTbl_treat[resultsTbl_treat.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_treat_filtered, file="data/exactTest_tribolium_treat_filtered.csv", sep=",", row.names=TRUE)

## cntrl_4h vs cntrl_24h
#Perform an exact test for cntrl_4h vs cntrl_24h
tested_cntrl <- exactTest(list, pair=c("cntrl_24h", "cntrl_4h"))

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_cntrl))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
png("plots/dev/exactTest_tribolium_cntrl_DE.png")
plotMD(tested_cntrl)
abline(h=c(-1, 1), col="blue")
dev.off()

#Create a table of DE genes filtered by FDR
resultsTbl_nctrl <- topTags(tested_cntrl, n=nrow(tested_cntrl$table))$table

#Identify significantly DE genes
resultsTbl_nctrl$topDE <- "NA"
resultsTbl_nctrl$topDE[resultsTbl_nctrl$logFC > 1 & resultsTbl_nctrl$FDR < 0.05] <- "UP"
resultsTbl_nctrl$topDE[resultsTbl_nctrl$logFC < -1 & resultsTbl_nctrl$FDR < 0.05] <- "DOWN"

#Create volcano plot
png("plots/dev/exactTest_tribolium_cntrl_volcano.png")
ggplot(data=resultsTbl_nctrl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))
dev.off()

#Create filtered results table of DE genes
resultsTbl_ctrl.keep <- resultsTbl_nctrl$FDR <= 0.05
resultsTbl_cntrl_filtered <- resultsTbl_nctrl[resultsTbl_ctrl.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_cntrl_filtered, file="data/exactTest_tribolium_cntrl_filtered.csv", sep=",", row.names=TRUE)


##
# Results Exploration
##
# retrieve set of DE gene names for 24h contrast
geneSet_24h <- rownames(resultsTbl_24h_filtered)

# retrieve set of DE gene names for treat contrast
geneSet_treat <- rownames(resultsTbl_treat_filtered)

# retrieve set of DE gene names for cntrl contrast
geneSet_cntrl <- rownames(resultsTbl_cntrl_filtered)

# create combined list of DE gene names
list_venn <- list(h24 = geneSet_24h, 
                  treat = geneSet_treat, 
                  cntrl = geneSet_cntrl)

# create venn diagram
png("plots/dev/exactTest_tribolium_venn.png")
ggVennDiagram(list_venn, label_alpha=0.25, category.names = c("24h","treat","cntrl")) +
  scale_color_brewer(palette = "Paired")
dev.off()
