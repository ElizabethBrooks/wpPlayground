##
# Setup
##

# set the working directory
setwd("/YOUR/FILE/PATH/")

# install BiocManager, if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# install edgeR using BiocManager, if not already installed
BiocManager::install("edgeR")

# install ggplot2 and ghibli, if not already installed
install.packages("ggplot2")
install.packages("ghibli")

# import installed libraries
library(edgeR)
library(ggplot2)
library(ghibli)

# import gene count data
tribolium_counts <- read.csv("TriboliumCounts.csv", row.names="X")

# add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("treat_4h",3), rep("cntrl_24h",3), rep("treat_24h",3)))

# verify grouping factor ordering
group
colnames(tribolium_counts)

# create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)

# change the graphical parameters
par(mfrow=c(9,3))

# view all available ghibli palettes
for(i in names(ghibli_palettes)) print(ghibli_palette(i))

# close the plot and return the display to the default graphical parameters
dev.off()

# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")

# view the selected color palette
ghibli_colors


##
# Normalization
##

# plot the library sizes before normalization and write to a jpg file
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")

# filter the list of gene counts based on expression levels
keep <- filterByExpr(list)

# view the number of filtered genes
table(keep)

# remove genes that are not expressed in either experimental condition
list <- list[keep, , keep.lib.sizes=FALSE]

# calculate scaling factors
list <- calcNormFactors(list)

# view normalization factors
list$samples

# compute counts per million (CPM) using normalized library sizes
normList <- cpm(list, normalized.lib.sizes=TRUE)

# write the normalized counts to a csv file
write.table(normList, file="tribolium_normalizedCounts.csv", sep=",", row.names=TRUE)


##
# Data Exploration
##

# vector of shape numbers for the MDS plot
points <- c(0,1,15,16)

# vector of colors for the MDS plot
colors <- rep(c(ghibli_colors[3], ghibli_colors[6]), 2)

# add extra space to right of plot area and change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

# MDS plot with distances approximating log2 fold changes
plotMDS(list, col=colors[group], pch=points[group])

# place the legend outside the right side of the plot
legend("topright", inset=c(-0.4,0), legend=levels(group), pch=points, col=colors)

# calculate the log CPM of the gene count data
logcpm <- cpm(list, log=TRUE)

# draw a heatmap of individual RNA-seq samples using moderated log CPM
heatmap(logcpm)

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
list <- estimateDisp(list)

# plot dispersion estimates and biological coefficient of variation
plotBCV(list)


##
# Pairwise Test Contrasts
##

# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

## treat_4h vs ctrl_4h
# perform an exact test for treat_4h vs ctrl_4h
tested_4h <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_4h))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_4h)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_4h <- topTags(tested_4h, n=nrow(tested_4h$table))$table

# add column for identifying direction of DE gene expression
resultsTbl_4h$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_4h$topDE[resultsTbl_4h$logFC > 1 & resultsTbl_4h$FDR < 0.05] <- "UP"

# identify significantly down DE genes
resultsTbl_4h$topDE[resultsTbl_4h$logFC < -1 & resultsTbl_4h$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=resultsTbl_4h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))

# identify significantly DE genes by FDR
resultsTbl_4h.keep <- resultsTbl_4h$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_4h_filtered <- resultsTbl_4h[resultsTbl_4h.keep,]

# write the filtered results of the exact tests to a csv file
write.table(resultsTbl_4h_filtered, file="exactTest_4h_filtered.csv", sep=",", row.names=TRUE)

## treat_24h vs ctrl_24h
# perform an exact test for treat_24h vs ctrl_24h
tested_24h <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_24h))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_24h)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_24h <- topTags(tested_24h, n=nrow(tested_24h$table))$table

# add column for identifying direction of DE gene expression
resultsTbl_24h$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_24h$topDE[resultsTbl_24h$logFC > 1 & resultsTbl_24h$FDR < 0.05] <- "UP"

# identify significantly down DE genes
resultsTbl_24h$topDE[resultsTbl_24h$logFC < -1 & resultsTbl_24h$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=resultsTbl_24h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))

# identify significantly DE genes by FDR
resultsTbl_24h.keep <- resultsTbl_24h$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_24h_filtered <- resultsTbl_24h[resultsTbl_24h.keep,]

# write the filtered results of the exact tests to a csv file
write.table(resultsTbl_24h_filtered, file="exactTest_24h_filtered.csv", sep=",", row.names=TRUE)

## treat_4h vs treat_24h
# perform an exact test for treat_4h vs treat_24h
tested_treat <- exactTest(list, pair=c("treat_24h", "treat_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_treat))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_treat)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_treat <- topTags(tested_treat, n=nrow(tested_treat$table))$table

# add column for identifying direction of DE gene expression
resultsTbl_treat$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_treat$topDE[resultsTbl_treat$logFC > 1 & resultsTbl_treat$FDR < 0.05] <- "UP"

# identify significantly down DE genes
resultsTbl_treat$topDE[resultsTbl_treat$logFC < -1 & resultsTbl_treat$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=resultsTbl_treat, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))

# identify significantly DE genes by FDR
resultsTbl_treat.keep <- resultsTbl_treat$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_treat_filtered <- resultsTbl_treat[resultsTbl_treat.keep,]

# write the filtered results of the exact tests to a csv file
write.table(resultsTbl_treat_filtered, file="exactTest_treat_filtered.csv", sep=",", row.names=TRUE)

## cntrl_4h vs ctrl_24h
# perform an exact test for cntrl_4h vs ctrl_24h
tested_cntrl <- exactTest(list, pair=c("cntrl_24h", "cntrl_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_cntrl))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_cntrl)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_nctrl <- topTags(tested_cntrl, n=nrow(tested_cntrl$table))$table

# add column for identifying direction of DE gene expression
resultsTbl_nctrl$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_nctrl$topDE[resultsTbl_nctrl$logFC > 1 & resultsTbl_nctrl$FDR < 0.05] <- "UP"

# identify significantly down DE genes
resultsTbl_nctrl$topDE[resultsTbl_nctrl$logFC < -1 & resultsTbl_nctrl$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=resultsTbl_nctrl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("UP", "DOWN"))

# identify significantly DE genes by FDR
resultsTbl_ctrl.keep <- resultsTbl_nctrl$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_cntrl_filtered <- resultsTbl_nctrl[resultsTbl_ctrl.keep,]

# write the filtered results of the exact tests to a csv file
write.table(resultsTbl_cntrl_filtered, file="exactTest_cntrl_filtered.csv", sep=",", row.names=TRUE)


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
ggVennDiagram(list_venn, label_alpha=0.25, category.names = c("24h","treat","cntrl")) +
  scale_color_brewer(palette = "Paired")
