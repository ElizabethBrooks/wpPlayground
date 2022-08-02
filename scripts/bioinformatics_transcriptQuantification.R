# R script to quantify transcriptomic data for the red flour beetle

# set the working directory
setwd("/YOUR/DIRECTORY/PATH")

# install the BiocManager followed by the Rsubread package, if necessary
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

# load the Rsubread library
library("Rsubread")

# checkout information for the featureCounts command
?featureCounts

# quantify read fragments that align (map) to general features
# contrl samples at 4h
cntrl1_fc_4h <- featureCounts(files="SRR8288561_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl2_fc_4h <- featureCounts(files="SRR8288562_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl3_fc_4h <- featureCounts(files="SRR8288563_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
# contrl samples at 24h
cntrl1_fc_24h <- featureCounts(files="SRR8288558_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl2_fc_24h <- featureCounts(files="SRR8288567_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl3_fc_24h <- featureCounts(files="SRR8288568_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
# treat samples at 4h
treat1_fc_4h <- featureCounts(files="SRR8288564_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat2_fc_4h <- featureCounts(files="SRR8288557_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat3_fc_4h <- featureCounts(files="SRR8288560_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
# treat samples at 24h
treat1_fc_24h <- featureCounts(files="SRR8288559_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat2_fc_24h <- featureCounts(files="SRR8288565_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat3_fc_24h <- featureCounts(files="SRR8288566_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)

# merged all of the sample counts into one data frame
tribolium_counts <- data.frame(
  cntrl1_4h = unname(cntrl1_fc_4h$counts),
  cntrl2_4h = unname(cntrl2_fc_4h$counts),
  cntrl3_4h = unname(cntrl3_fc_4h$counts),
  treat1_4h = unname(cntrl1_fc_24h$counts),
  treat2_4h = unname(cntrl2_fc_24h$counts),
  treat3_4h = unname(cntrl3_fc_24h$counts),
  cntrl1_24h = unname(treat1_fc_4h$counts),
  cntrl2_24h = unname(treat2_fc_4h$counts),
  cntrl3_24h = unname(treat3_fc_4h$counts),
  treat1_24h = unname(treat1_fc_24h$counts),
  treat2_24h = unname(treat2_fc_24h$counts),
  treat3_24h = unname(treat3_fc_24h$counts)
)

# look at the first few lines of the data frame
head(tribolium_counts)

# update the row names for the merged counts data frame
rownames(tribolium_counts) <- rownames(cntrl1_fc_4h$counts)

# look at the first few lines of the updated data frame
head(tribolium_counts)

# export the merged and formatted gene count data
write.csv(tribolium_counts, "TriboliumCounts.csv")
