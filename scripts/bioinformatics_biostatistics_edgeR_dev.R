##
# Setup
##

# set the working directory
#setwd("/YOUR/FILE/PATH/")
setwd("/Users/bamflappy/Repos/wpPlayground/")

# install libraries, if necessary
#https://github.com/ewenme/ghibli
#BiocManager::install("edgeR")
#install.packages("ggplot2")
#install.packages("ghibli")

# import libraries
library(edgeR)
library(ggplot2)
library(ghibli)

# import gene count data
tribolium_counts <- read.csv("data/TriboliumCounts.csv", row.names="X")

#Add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("treat_4h",3), rep("cntrl_24h",3), rep("treat_24h",3)))

# view available palettes
par(mfrow=c(9,3))
for(i in names(ghibli_palettes)) print(ghibli_palette(i))
dev.off()

# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

##
# Normalization
##

#Create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)

#Plot the library sizes before normalization and write to a jpg file
jpeg("plots/dev/exactTest_librarySizes.jpg")
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

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4], ghibli_colors[1]), 2)

# Add extra space to right of plot area; change clipping to figure
jpeg("plots/dev/exactTest_MDS.jpg")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotMDS(list, col=colors[group], pch=points[group])
legend("topright", inset=c(-0.28,0), legend=levels(group), pch=points, col=colors)
dev.off()

#Calculate the log CPM of the gene count data
logcpm <- cpm(list, log=TRUE)

#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
jpeg("plots/dev/exactTest_hclust.jpg")
heatmap(logcpm)
dev.off()

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#View dispersion estimates and biological coefficient of variation
jpeg("plots/dev/exactTest_BCV.jpg")
plotBCV(list)
dev.off()


##
# Pairwise Tests
##

##
#Perform an exact test for treat_4h vs ctrl_4h
tested_4h <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

#Create results table of DE genes
resultsTbl_4h <- topTags(tested_4h, n=nrow(tested_4h$table))$table

#Create a table of DE genes filtered by FDR
resultsTbl_4h.keep <- resultsTbl_4h$FDR <= 0.05
resultsTbl_4h_filtered <- resultsTbl_4h[resultsTbl_4h.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_4h_filtered, file="data/exactTest_4h_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_4h$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_4h))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("plots/dev/exactTest_4h_DE.jpg")
plotMD(tested_4h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("plots/dev/exactTest_4h_smear.jpg")
plotSmear(tested_4h)
dev.off()

#Identify significantly DE genes
resultsTbl_4h$topDE <- "NA"
resultsTbl_4h$topDE[resultsTbl_4h$logFC > 1 & resultsTbl_4h$FDR < 0.05] <- "UP"
resultsTbl_4h$topDE[resultsTbl_4h$logFC < -1 & resultsTbl_4h$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/exactTest_4h_volcano.jpg")
ggplot(data=resultsTbl_4h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()

##
#Perform an exact test for treat_24h vs ctrl_24h
tested_24h <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

#Create a table of DE genes filtered by FDR
resultsTbl_24h <- topTags(tested_24h, n=nrow(tested_24h$table))$table

#Create filtered results table of DE genes
resultsTbl_24h.keep <- resultsTbl_24h$FDR <= 0.05
resultsTbl_24h_filtered <- resultsTbl_24h[resultsTbl_24h.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_24h_filtered, file="data/exactTest_24h_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_24h$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_24h))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("plots/dev/exactTest_24h_DE.jpg")
plotMD(tested_24h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("plots/dev/exactTest_24h_smear.jpg")
plotSmear(tested_24h)
dev.off()

#Identify significantly DE genes
resultsTbl_24h$topDE <- "NA"
resultsTbl_24h$topDE[resultsTbl_24h$logFC > 1 & resultsTbl_24h$FDR < 0.05] <- "UP"
resultsTbl_24h$topDE[resultsTbl_24h$logFC < -1 & resultsTbl_24h$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/exactTest_24h_volcano.jpg")
ggplot(data=resultsTbl_24h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()

##
#Perform an exact test for treat_4h vs treat_24h
tested_treat <- exactTest(list, pair=c("treat_24h", "treat_4h"))

#Create a table of DE genes filtered by FDR
resultsTbl_treat <- topTags(tested_treat, n=nrow(tested_treat$table))$table

#Create filtered results table of DE genes
resultsTbl_treat.keep <- resultsTbl_treat$FDR <= 0.05
resultsTbl_treat_filtered <- resultsTbl_treat[resultsTbl_treat.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_treat_filtered, file="data/exactTest_treat_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_treat$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_treat))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("plots/dev/exactTest_treat_DE.jpg")
plotMD(tested_treat)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("plots/dev/exactTest_treat_smear.jpg")
plotSmear(tested_treat)
dev.off()

#Identify significantly DE genes
resultsTbl_treat$topDE <- "NA"
resultsTbl_treat$topDE[resultsTbl_treat$logFC > 1 & resultsTbl_treat$FDR < 0.05] <- "UP"
resultsTbl_treat$topDE[resultsTbl_treat$logFC < -1 & resultsTbl_treat$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/exactTest_treat_volcano.jpg")
ggplot(data=resultsTbl_treat, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()

##############
#Perform an exact test for cntrl_4h vs ctrl_24h
tested_cntrl <- exactTest(list, pair=c("cntrl_24h", "cntrl_4h"))

#Create a table of DE genes filtered by FDR
resultsTbl_nctrl <- topTags(tested_cntrl, n=nrow(tested_cntrl$table))$table

#Create filtered results table of DE genes
resultsTbl_ctrl.keep <- resultsTbl_nctrl$FDR <= 0.05
resultsTbl_cntrl_filtered <- resultsTbl_nctrl[resultsTbl_ctrl.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_cntrl_filtered, file="data/exactTest_cntrl_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_cntrl$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_cntrl))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("plots/dev/exactTest_cntrl_DE.jpg")
plotMD(tested_cntrl)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("plots/dev/exactTest_cntrl_smear.jpg")
plotSmear(tested_cntrl)
dev.off()

#Identify significantly DE genes
resultsTbl_nctrl$topDE <- "NA"
resultsTbl_nctrl$topDE[resultsTbl_nctrl$logFC > 1 & resultsTbl_nctrl$FDR < 0.05] <- "UP"
resultsTbl_nctrl$topDE[resultsTbl_nctrl$logFC < -1 & resultsTbl_nctrl$FDR < 0.05] <- "DOWN"

#Create volcano plot
jpeg("plots/dev/exactTest_cntrl_volcano.jpg")
ggplot(data=resultsTbl_nctrl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset)
dev.off()
