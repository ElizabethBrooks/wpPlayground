##
# General Setup
##

# set the working directory
setwd("/YOUR/FILE/PATH/")


##
# Packages
##

# install packages, if necessary
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("ggplot2")
install.packages("ghibli")
install.packages("ggVennDiagram")

# import libraries
library(edgeR)
library(ggplot2)
library(ghibli)
library(ggVennDiagram)


##
# Plotting Palettes
##

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

# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])


##
# Data
##

# import gene count data
tribolium_counts <- read.csv("TriboliumCounts.csv", row.names="X")


##
# Pairwise Setup
##

# add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("treat_4h",3), rep("cntrl_24h",3), rep("treat_24h",3)))

# create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)


##
# Pairwise Normalization
##

# plot the library sizes before normalization
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")

# filter the list of gene counts based on expression levels
keep <- filterByExpr(list)

# view the number of filtered genes
table(keep)

# remove genes that are not expressed in either experimental condition
list <- list[keep, , keep.lib.sizes=FALSE]

# calculate scaling factors
list <- calcNormFactors(list)

# compute counts per million (CPM) using normalized library sizes
normList <- cpm(list, normalized.lib.sizes=TRUE)


##
# Pairwise Data Exploration
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

# close the plot
dev.off()

# calculate the log CPM of the gene count data
logcpm <- cpm(list, log=TRUE)

# draw a heatmap of individual RNA-seq samples using moderated log CPM
heatmap(logcpm)


##
# Pairwise Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
list <- estimateDisp(list)

# plot dispersion estimates and biological coefficient of variation
plotBCV(list)


##
# Pairwise Contrasts
##

### 
## treat_4h vs cntrl_4h
###

# perform an exact test for treat_4h vs cntrl_4h
tested_4h <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_4h))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_4h)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_4h <- topTags(tested_4h, n=nrow(tested_4h$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_4h$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_4h$topDE[resultsTbl_4h$logFC > 1 & resultsTbl_4h$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_4h$topDE[resultsTbl_4h$logFC < -1 & resultsTbl_4h$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_4h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

###
## treat_24h vs cntrl_24h
###

# perform an exact test for treat_24h vs cntrl_24h
tested_24h <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_24h))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_24h)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_24h <- topTags(tested_24h, n=nrow(tested_24h$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_24h$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_24h$topDE[resultsTbl_24h$logFC > 1 & resultsTbl_24h$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_24h$topDE[resultsTbl_24h$logFC < -1 & resultsTbl_24h$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_24h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
resultsTbl_24h.keep <- resultsTbl_24h$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_24h_filtered <- resultsTbl_24h[resultsTbl_24h.keep,]

###
## treat_4h vs treat_24h
###

# perform an exact test for treat_4h vs treat_24h
tested_treat <- exactTest(list, pair=c("treat_24h", "treat_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_treat))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_treat)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_treat <- topTags(tested_treat, n=nrow(tested_treat$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_treat$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_treat$topDE[resultsTbl_treat$logFC > 1 & resultsTbl_treat$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_treat$topDE[resultsTbl_treat$logFC < -1 & resultsTbl_treat$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_treat, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
resultsTbl_treat.keep <- resultsTbl_treat$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_treat_filtered <- resultsTbl_treat[resultsTbl_treat.keep,]

###
## cntrl_4h vs cntrl_24h
###

# perform an exact test for cntrl_4h vs cntrl_24h
tested_cntrl <- exactTest(list, pair=c("cntrl_24h", "cntrl_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_cntrl))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_cntrl)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_ncntrl <- topTags(tested_cntrl, n=nrow(tested_cntrl$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_ncntrl$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_ncntrl$topDE[resultsTbl_ncntrl$logFC > 1 & resultsTbl_ncntrl$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_ncntrl$topDE[resultsTbl_ncntrl$logFC < -1 & resultsTbl_ncntrl$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_ncntrl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
resultsTbl_cntrl.keep <- resultsTbl_ncntrl$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_cntrl_filtered <- resultsTbl_ncntrl[resultsTbl_cntrl.keep,]


##
# Pairwise Results Exploration
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


##
# GLM Design
##

# import grouping factor
glm_targets <- read.csv(file="groupingFactors_tribolium.csv", row.names="sample")

# setup a design matrix
glm_group <- factor(paste(glm_targets$treatment, glm_targets$hours, sep="."))

# create DGE list object
glm_list <- DGEList(counts=tribolium_counts, group=glm_group)

# add the sample names
colnames(glm_list) <- rownames(glm_targets)

# parametrize the experimental design with a one-way layout 
glm_design <- model.matrix(~ 0 + glm_group)

# add group names
colnames(glm_design) <- levels(glm_group)


##
# GLM Normalization
##

# filter the list of gene counts based on expression levels
glm_keep <- filterByExpr(glm_list)

# view the number of filtered genes
table(glm_keep)

# remove genes that are not expressed in either experimental condition
glm_list <- glm_list[glm_keep, , keep.lib.sizes=FALSE]

# calculate scaling factors
glm_list <- calcNormFactors(glm_list)

# compute counts per million (CPM) using normalized library sizes
norm_glm_list <- cpm(glm_list, normalized.lib.sizes=TRUE)


##
# GLM Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
glm_list <- estimateDisp(glm_list, glm_design, robust=TRUE)

# estimate the QL dispersions
glm_fit <- glmQLFit(glm_list, glm_design, robust=TRUE)

# plot the QL dispersions
plotQLDisp(glm_fit)



##
# GLM Contrasts
##

###
## treatment
###

# examine the overall effect of treatment
con.treatment <- makeContrasts(set.treatment = 
                               (treat.4h + treat.24h)/2
                               - (cntrl.4h + cntrl.24h)/2,
                               levels=glm_design)

# conduct gene wise statistical tests
anov.treatment <- glmTreat(glm_fit, contrast=con.treatment)

# view summary of DE genes
summary(decideTests(anov.treatment))

# create MD plot of DE genes
plotMD(anov.treatment)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_treatment <- topTags(anov.treatment, n=nrow(anov.treatment$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_treatment$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_treatment$topDE[tagsTbl_treatment$logFC > 1 & tagsTbl_treatment$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_treatment$topDE[tagsTbl_treatment$logFC < -1 & tagsTbl_treatment$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_treatment, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_treatment.glm_keep <- tagsTbl_treatment$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_treatment_filtered <- tagsTbl_treatment[tagsTbl_treatment.glm_keep,]

###
## hours
###

# examine the overall effect of hours
con.hours <- makeContrasts(set.hours = 
                           (cntrl.24h + treat.24h)/2
                           - (cntrl.4h + treat.4h)/2,
                           levels=glm_design)

# conduct gene wise statistical tests
anov.hours <- glmTreat(glm_fit, contrast=con.hours)

# view summary of DE genes
summary(decideTests(anov.hours))

# create MD plot of DE genes
plotMD(anov.hours)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_hours <- topTags(anov.hours, n=nrow(anov.hours$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_hours$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_hours$topDE[tagsTbl_hours$logFC > 1 & tagsTbl_hours$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_hours$topDE[tagsTbl_hours$logFC < -1 & tagsTbl_hours$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_hours, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_hours.glm_keep <- tagsTbl_hours$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_hours_filtered <- tagsTbl_hours[tagsTbl_hours.glm_keep,]

###
## interaction
###

# examine any interaction effect
con.interaction <- makeContrasts(set.interaction = 
                                 ((treat.4h + treat.24h)/2
                                 - (cntrl.4h + cntrl.24h)/2)
                                 - ((cntrl.24h + treat.24h)/2
                                 - (cntrl.4h + treat.4h)/2),
                                 levels=glm_design)

# conduct gene wise statistical tests
anov.interaction <- glmTreat(glm_fit, contrast=con.interaction)

# view summary of DE genes
summary(decideTests(anov.interaction))

# create MD plot of DE genes
plotMD(anov.interaction)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_inter <- topTags(anov.interaction, n=nrow(anov.interaction$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_inter$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_inter$topDE[tagsTbl_inter$logFC > 1 & tagsTbl_inter$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_inter$topDE[tagsTbl_inter$logFC < -1 & tagsTbl_inter$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_inter, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_inter.glm_keep <- tagsTbl_inter$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_inter_filtered <- tagsTbl_inter[tagsTbl_inter.glm_keep,]


##
# GLM Results Exploration
##

# retrieve set of DE gene names for hours contrast
geneSet_hours <- rownames(tagsTbl_hours_filtered)

# retrieve set of DE gene names for interaction contrast
geneSet_interaction <- rownames(tagsTbl_inter_filtered)

# create combined glm_list of DE gene names
glm_list_venn <- list(hours = geneSet_hours, 
                          interaction = geneSet_interaction)

# create venn diagram
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("hours","interaction")) +
  scale_color_brewer(palette = "Paired")


## 
# Saving Plots
##
# save the glm results venn diagram plot to a file
png("glm_tribolium_venn.png")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("hours","interaction")) +
  scale_color_brewer(palette = "Paired")
dev.off()


##
# Saving Tables
##
# write the table of normalized counts to a file
write.table(norm_glm_list, file="Tribolium_normalizedCounts.csv", sep=",", row.names=TRUE)
