# set the working directory
#setwd("/YOUR/FILE/PATH/")
setwd("/Users/bamflappy/Repos/rPlayground/")

# install libraries, if necessary
#https://github.com/ewenme/ghibli
#BiocManager::install("edgeR")
#install.packages("statmod")
#install.packages("ghibli")

# import libraries
library(edgeR)
library(ggplot2)
library(ghibli)

# import gene count data
tribolium_counts <- read.csv("data/TriboliumCounts.csv", row.names="X")

# trim the data table
#countsTable <- head(inputTable, - 5)

# import grouping factor
targets <- read.csv(file="data/groupingFactors_tribolium.csv")

# view available palettes
par(mfrow=c(9,3))
for(i in names(ghibli_palettes)) print(ghibli_palette(i))
dev.off()

# retrieve the vector of colors associated with Zissou1
#Error in check_for_XQuartz() : 
#  X11 library is missing: install XQuartz from www.xquartz.org
#brew install --cask xquartz
(ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete"))
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])


##
# Prep Stage
##

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$hours,sep="."))
#cbind(targets,Group=group)

#Create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)
colnames(list) <- rownames(targets)
head(list)

#Plot the library sizes before normalization
jpeg("plots/dev/glm_tribolium_librarySizes.jpg")
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")
dev.off()

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)

#Retrieve normalized counts
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Write the normalized counts to a file
write.table(normList, file="data/tribolium_normalizedCounts.csv", sep=",", row.names=TRUE)

#View normalization factors
list$samples
dim(list)

#Verify TMM normalization using a MD plot
# and write the plot to jpg file
jpeg("plots/dev/glm_tribolium_MD_afterNorm.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4], ghibli_colors[1]), 2)

#Write plot without legend to file
jpeg("plots/dev/glm_tribolium_MDS_withoutLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
dev.off()

#Write plot with legend to file
jpeg("plots/dev/glm_tribolium_MDS_withLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion

#Visualize the dispersion estimates with a BCV plot
#Write plot to file
jpeg("glm_tribolium_BCV.jpg")
plotBCV(list)
dev.off()

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)

#Write plot to file
jpeg("glm_tribolium_QLDisp.jpg")
plotQLDisp(fit)
dev.off()


##
# DE Stage
##

##
#Test whether the average across all cntrl groups is equal to the average across
#all treat groups, to examine the overall effect of treatment
con.treatment <- makeContrasts(set.treatment = (treat.4h + treat.24h)/2
                                 - (cntrl.4h + cntrl.24h)/2,
                                 levels=design)

#Look at genes with significant expression across all UV groups
anov.treatment <- glmTreat(fit, contrast=con.treatment, lfc=log2(1.2))
summary(decideTests(anov.treatment))

#Write MD plot to file
jpeg("plots/dev/glm_tribolium_treatment_MD.jpg")
plotMD(anov.treatment)
abline(h=c(-1, 1), col="blue")
dev.off()

#Generate table of DE genes
tagsTbl_treatment <- topTags(anov.treatment, n=nrow(anov.treatment$table), adjust.method="fdr")$table
write.table(tagsTbl_treatment, file="data/glm_tribolium_treatment.csv", sep=",", row.names=TRUE)

#Identify significantly DE genes
tagsTbl_treatment$topDE <- "NA"
tagsTbl_treatment$topDE[tagsTbl_treatment$logFC > 1 & tagsTbl_treatment$FDR < 0.05] <- "UP"
tagsTbl_treatment$topDE[tagsTbl_treatment$logFC < -1 & tagsTbl_treatment$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/glm_tribolium_treatment_volcano.jpg")
ggplot(data=tagsTbl_treatment, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()

##
#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.hours <- makeContrasts(set.hours = (cntrl.24h + treat.24h)/2
                            - (cntrl.4h + treat.4h)/2,
                            levels=design)

#Look at genes with significant expression across all UV groups
anov.hours <- glmTreat(fit, contrast=con.hours, lfc=log2(1.2))
summary(decideTests(anov.hours))

#Write plot to file
jpeg("plots/dev/glm_tribolium_hours_MD.jpg")
plotMD(anov.hours)
abline(h=c(-1, 1), col="blue")
dev.off()

#Generate table of DE genes
tagsTbl_hours <- topTags(anov.hours, n=nrow(anov.hours$table), adjust.method="fdr")$table
write.table(tagsTbl_hours, file="data/glm_tribolium_hours.csv", sep=",", row.names=TRUE)

#Identify significantly DE genes
tagsTbl_hours$topDE <- "NA"
tagsTbl_hours$topDE[tagsTbl_hours$logFC > 1 & tagsTbl_hours$FDR < 0.05] <- "UP"
tagsTbl_hours$topDE[tagsTbl_hours$logFC < -1 & tagsTbl_hours$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/glm_tribolium_hours_volcano.jpg")
ggplot(data=tagsTbl_hours, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()

##
#Test whether there is an interaction effect
con.interaction <- makeContrasts(set.interaction = ((treat.4h + treat.24h)/2
                                                    - (cntrl.4h + cntrl.24h)/2)
                                 - ((cntrl.24h + treat.24h)/2
                                    - (cntrl.4h + treat.4h)/2),
                                 levels=design)

#Look at genes with significant expression
anov.interaction <- glmTreat(fit, contrast=con.interaction, lfc=log2(1.2))
summary(decideTests(anov.interaction))

#Write plot to file
jpeg("plots/dev/glm_tribolium_interaction_MD.jpg")
plotMD(anov.interaction)
abline(h=c(-1, 1), col="blue")
dev.off()

#Generate table of DE genes
tagsTbl_inter <- topTags(anov.interaction, n=nrow(anov.interaction$table), adjust.method="fdr")$table
write.table(tagsTbl_inter, file="data/glm_tribolium_interaction.csv", sep=",", row.names=TRUE)

#Identify significantly DE genes
tagsTbl_inter$topDE <- "NA"
tagsTbl_inter$topDE[tagsTbl_inter$logFC > 1 & tagsTbl_inter$FDR < 0.05] <- "UP"
tagsTbl_inter$topDE[tagsTbl_inter$logFC < -1 & tagsTbl_inter$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/glm_tribolium_interaction_volcano.jpg")
ggplot(data=tagsTbl_inter, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()


##
# Pairwise Results Exploration
##
# retrieve set of DE gene names for 24h contrast
geneSet_treatment <- rownames(tagsTbl_treatment)

# retrieve set of DE gene names for treat contrast
geneSet_hours <- rownames(tagsTbl_hours)

# retrieve set of DE gene names for cntrl contrast
geneSet_interaction <- rownames(tagsTbl_inter)

# create combined list of DE gene names
list_venn <- list(treatment = geneSet_treatment, 
                  hours = geneSet_hours, 
                  interaction = geneSet_interaction)

# create venn diagram
png("plots/dev/glm_tribolium_venn.png")
ggVennDiagram(list_venn, label_alpha=0.25, category.names = c("treatment","hours","interaction")) +
  scale_color_brewer(palette = "Paired")
dev.off()

