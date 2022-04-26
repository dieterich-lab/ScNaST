#! /usr/bin/env Rscript

## module load R/4.1.1_deb10

# for data import
library(zellkonverter)
# for analysis
library('SPARK')
library(tidyverse)
library(SummarizedExperiment)


parent <- '/prj/Florian_Leuschner_spatial/analysis/Illumina'
raw_file <- 'concatenated_raw_bbknn.h5ad'

outname <- 'raw_bbkn'

## Data pre-processing and wrangling

sce <- readH5AD(file.path(parent, 'data', raw_file))
rawcount <- as.matrix(sce@assays@data$X)
metaData <- data.frame(sce@colData)


## Test

# run each slice separately

for (sample in levels(metaData$library_id)) {
    
    rawcount.subset <- rawcount[,grep(paste(sample, "$", sep=""), colnames(rawcount))]
    metaData.subset <- metaData[metaData$library_id==sample,][,c('array_row', 'array_col')]
    info <- cbind.data.frame(x=as.numeric(metaData.subset$array_row),
                             y=as.numeric(metaData.subset$array_col),
                             total_counts=apply(rawcount.subset,2,sum))
                             
    stopifnot(all.equal(rownames(info), colnames(rawcount.subset)))
    
    # use default
    spark <- CreateSPARKObject(counts=rawcount.subset, 
                               location=info[,1:2],
                               percentage = 0.1, 
                               min_total_counts = 10)

    # total counts for each cell/spot
    spark@lib_size <- apply(spark@counts, 2, sum)
    
    # fit statistical model under the null hypothesis
    spark <- spark.vc(spark, 
                      covariates = NULL, 
                      lib_size = spark@lib_size, 
                      num_core = 40,
                      verbose = F)
    
    # test genes that vary spatially 
    # by default kernel matrices are computed by coordinates

    # calculate pval
    spark <- spark.test(spark, 
                        check_positive = T, 
                        verbose = F)
                        
    saveRDS(spark, file=file.path(parent, 'spark', 'results', paste("spark-", sample, ".RDS", sep="")))
                        
}
