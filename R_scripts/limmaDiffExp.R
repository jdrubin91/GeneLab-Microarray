#!/usr/bin/env Rscript

#source("http://bioconductor.org/biocLite.R")
#biocLite("Risa")
#biocLite("limma")

require(Risa)
require(limma)

#Read in ISA tab file and extract assay file
rsaIN = readISAtab("~/Google\ Drive/NASA/home/metadata/GLDS-4/metadata/")
assayFactors = rsaIN@assay.files[[1]]

#From assay file, extract column containing 'Factor Value' (the rest of the script assumes there's only one of these columns...)
factorValues <- assayFactors[,grepl("Factor Value",colnames(assayFactors))]
#Also extract sample names from assay file, these should be in the same order as factor values
sampNames <- assayFactors[,grepl("Sample Name",colnames(assayFactors))]

#Read in an expression value txt file
eset <- read.delim("~/Google\ Drive/NASA/home/batch_out/GLDS-4/microarray/exprsValues.txt")
#Modify matrix so that the first column becomes the row names
rownames(eset) <- eset[,1]
eset <- eset[,-1]

#From the eset matrix, determine which columns correspond to which factor values
esetSamples <- colnames(eset)
esetNewColNames <- c(1:length(esetSamples))
j <- 1
for(i in sampNames){
  factor <- factorValues[grepl(i,esetSamples)]
  esetNewColNames[j] <- factor
  j <- j + 1
}

#Create a design matrix based on the ordering of the columns within eset
design <- model.matrix(~0+esetNewColNames)

#Here we assume we're dealing with flight and ground samples... this part may not be necessary..
colnames(design) <- c("flight", "ground")

#This part of the script is straight from Limma documentation
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(flight-ground,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Here we write the results to a tab delimited text file that is ordered by adjusted p-value
#coef refers to which column is of interest (1 is log2FC), adjust refers to multiple hypothesis testing method ("BH" = Benjamini & Hochberg)
table <- data.frame(topTable(fit2, coef=1, n=Inf, adjust="BH"))
results <- decideTests(fit2)
write.table(table,file="~/Google\ Drive/NASA/home/batch_out/GLDS-4/microarray/diffExpression.txt",sep="\t")