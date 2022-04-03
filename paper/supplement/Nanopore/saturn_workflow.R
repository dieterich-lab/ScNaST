#! /usr/bin/env Rscript

## module load R/4.1.1_deb10

# for data import
library(AnnotationHub)
library(ensembldb)
library(zellkonverter)
library(openxlsx)
# for analysis
library(satuRn)
library(edgeR)
library(tidyverse)
library(SummarizedExperiment)


parent <- '/prj/Florian_Leuschner_spatial/analysis/Nanopore'
raw_file <- 'concatenated_raw_scnast.h5ad'
clusters_file <- 'anatomical_regions_and_cell_props_scnast.h5ad'

outname <- 'scnast'

# annotations
config <- yaml::yaml.load_file(file.path(parent, 'ScNaST', 'workflow', 'config.yaml', fsep=.Platform$file.sep))
gtf <- file.path(config$strg_base, 'gffcompare', 'gffcmp.multi_exons.annotated.gtf')
db_file <- gsub(pattern = "\\gtf$", "sqlite", basename(gtf))

## Data pre-processing and wrangling

# start by reading the AnnData to a SingleCellExperiment object
sce <- readH5AD(file.path(parent, 'Nanopore', 'data', clusters_file))

# get metadata
metaData <- data.frame(sce@colData)
metaData <- metaData[,c("anatomical_region"), drop=FALSE]

# here we need the raw counts
sce <- readH5AD(file.path(parent, 'Nanopore', 'data', raw_file))

# dgCMatrix sce@assays@data$counts
counts <- as.matrix(sce@assays@data$X)

stopifnot(all.equal(rownames(metaData), colnames(counts)))

# load the db
txdb <- loadDb(file.path(parent, 'Nanopore', 'dtu', db_file))

# transcript mapping to their corresponding genes
# this is based on the assignment made from gffcompare output, and stored in the txdb 

txs <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))

txInfo <- as.data.frame(matrix(data = NA, nrow = length(txs), ncol = 2)) 
# these names are hard coded in satuRn
colnames(txInfo) <- c("isoform_id", "gene_id")
txInfo$isoform_id <- txs$tx_name
txInfo$gene_id <- unlist(txs$gene_id) # is a CharacterList
rownames(txInfo) <- txInfo$isoform_id

# used features
# remove transcripts that are the only isoform expressed of a certain gene
txInfo <- txInfo[txInfo$isoform_id %in% rownames(counts),]
txInfo <- subset(txInfo, duplicated(gene_id) | duplicated(gene_id, fromLast=TRUE))
counts <- counts[which(rownames(counts) %in% txInfo$isoform_id),]
# satuRn requires the transcripts in the rowData and the transcripts in the count matrix to be in the same order.
txInfo <- txInfo[match(rownames(counts), rownames(txInfo)),]

stopifnot(all.equal(rownames(counts), rownames(txInfo)))


## Functions

# feature-level filtering
filter_func <- function(counts, txInfo, metaData, group=NULL, design=NULL, lib.size=NULL,
                        min.count=2.5, min.total.count=5, large.n=0, min.prop=0.7){
                        
    if(is.null(group)){
        if(is.null(design)) {
            stop("filter_func: group or design must be given.\n", call.=FALSE) 
        }
        group <- group
    } else {
        # overrides silently design
        group <- metaData[[group]]
        design <- NULL
    }

    filter_edgeR <- filterByExpr(counts,
                                 design=design,
                                 group=group,
                                 lib.size=lib.size,
                                 min.count=min.count,
                                 min.total.count=min.total.count,
                                 large.n=large.n,
                                 min.prop=min.prop) 

    counts <- counts[filter_edgeR, ]
    txInfo <- txInfo[which(txInfo$isoform_id %in% rownames(counts)),]
    txInfo <- subset(txInfo, duplicated(gene_id) | duplicated(gene_id, fromLast=TRUE))
    counts <- counts[which(rownames(counts) %in% txInfo$isoform_id),]
    # satuRn requires the transcripts in the rowData and the transcripts in the count matrix to be in the same order.
    txInfo <- txInfo[match(rownames(counts), rownames(txInfo)),]
    list(counts=counts, txInfo=txInfo)

}

# custom functions test for DTU 
test_region <- function(metaData){
    
    anatomical_region <- as.factor(metaData$anatomical_region)
    design <- model.matrix(~ 0 + anatomical_region)
    colnames(design) <- levels(anatomical_region)
    L <- matrix(0, ncol=6, nrow=ncol(design)) # initialize contrast matrix
    rownames(L) <- colnames(design)
    colnames(L) <- c("IvsRZ1", "IvsRZ2", "IvsBZ",
                     "BZvsRZ1", "BZvsRZ2", 
                     "RZ1vsRZ2")
    L[,1] <- c(-1,0,1,0)
    L[,2] <- c(0,-1,1,0)
    L[,3] <- c(0,0,1,-1)
    L[,4] <- c(-1,0,0,1)
    L[,5] <- c(0,-1,0,1)
    L[,6] <- c(1,-1,0,0)
    L
}

# Two-stage testing procedure with stageR                         
two_stage_func <- function(sumExp, L, alpha=0.05){

    # for each contrast, add to SummarizedExperiment (rowData)

    # prepare stageR input
    tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    colnames(tx2gene) <- c("transcript", "gene")
    # this is a trick... it seems that stageR strips after colon (used internally already)
    tx2gene$transcript <- gsub(":", "@", tx2gene$transcript)
    rownames(tx2gene) <- gsub(":", "@", rownames(tx2gene))

    
    for (c in colnames(L)) {
        tmp <- unlist(str_split(c, "vs"))
        num <- tmp[1]
        denom <- tmp[2]
        
        name <- paste("fitDTUResult", c, sep="_")
        pvals <- rowData(sumExp)[[name]]$empirical_pval

        # compute gene level q-values
        geneID <- factor(rowData(sumExp)$gene_id)
        geneSplit <- split(seq(along = geneID), geneID)
        pGene <- sapply(geneSplit, function(i) min(pvals[i]))
        pGene[is.na(pGene)] <- 1
        theta <- unique(sort(pGene))
        
        # gene-level significance testing
        q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
        qScreen <- rep(NA_real_, length(pGene))
        qScreen <- q[match(pGene, theta)]
        qScreen <- pmin(1, qScreen)
        names(qScreen) <- names(geneSplit)
        
        pConfirmation <- matrix(matrix(pvals),
                                ncol = 1,
                                dimnames = list(rownames(tx2gene), "transcript"))
                                
        # create a stageRTx object
        stageRObj <- stageR::stageRTx(pScreen = qScreen,
                                      pConfirmation = pConfirmation,
                                      pScreenAdjusted = TRUE,
                                      tx2gene = tx2gene)
        
        # perform the two-stage testing procedure
        stageRObj <- stageR::stageWiseAdjustment(object = stageRObj,
                                                 method = "dtu",
                                                 alpha = alpha,
                                                 allowNA = TRUE)
                                                
        # retrieves the adjusted p-values from the stageRTx object
        padj <- stageR::getAdjustedPValues(stageRObj,
                                           order = FALSE,
                                           onlySignificantGenes = FALSE)
        
        # fix
        padj$txID <- gsub("@", ":", padj$txID)
        # and make sure the order is the same
        padj <- padj[match(rownames(rowData(sumExp)[[name]]), padj$txID),]
        padj <- padj[,3:4]
        colnames(padj) <- c('padj_gene', 'padj_transcript')
        # add to SummarizedExperiment
        sumExp@elementMetadata[[name]] <- cbind(sumExp@elementMetadata[[name]], padj)
        
    }
    sumExp
}


## Test

out <- filter_func(counts, txInfo, metaData, group='anatomical_region',
                   min.count=1, min.total.count=2, min.prop=0.1)
stopifnot(all.equal(rownames(out$counts), rownames(out$txInfo)))

sumExp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=out$counts),
                                                     colData=metaData,
                                                     rowData=out$txInfo)
metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$anatomical_region)

# fit quasi-binomial generalized linear models models
sumExp <- satuRn::fitDTU(object=sumExp,
                         formula=~ 0 + anatomical_region,
                         parallel=TRUE,
                         BPPARAM=BiocParallel::bpparam(),
                         verbose=TRUE)
                         
L <- test_region(metaData)

sumExp <- satuRn::testDTU(object=sumExp,
                          contrasts=L,
                          diagplot1=FALSE,
                          diagplot2=FALSE,
                          sort=FALSE)

sumExp <- two_stage_func(sumExp, L)

saveRDS(sumExp, file=file.path(parent, 'Nanopore', 'dtu', 'results', paste("satuRn", outname, "fact_regs.RDS", sep="_")))


# also test region vs all (similar to rank markers)

for (reg in levels(metaData$anatomical_region)){
    
    metaData_loc <- metaData
    metaData_loc$region <- "rest"
    metaData_loc$region[metaData_loc$anatomical_region==reg] <- reg
    metaData_loc$region <- as.factor(metaData_loc$region)
    
    sumExp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=out$counts),
                                                         colData=metaData_loc,
                                                         rowData=out$txInfo)
    metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$region)
    sumExp <- satuRn::fitDTU(object=sumExp,
                             formula=~ 0 + region,
                             parallel=TRUE,
                             BPPARAM=BiocParallel::bpparam(),
                             verbose=TRUE)
                             
                             
    L <- matrix(c(1, -1), ncol=1, nrow=2) 
    rownames(L) <- c(reg, "rest")
    colnames(L) <- paste(reg, "_vs_rest", sep="")
    
    sumExp <- satuRn::testDTU(object=sumExp,
                              contrasts=L,
                              diagplot1=FALSE,
                              diagplot2=FALSE,
                              sort=FALSE)

    sumExp <- two_stage_func(sumExp, L)

    saveRDS(sumExp, file=file.path(parent, 'Nanopore', 'dtu', 'results', paste("satuRn", outname, gsub(' ', '', reg), "vs_rest.RDS", sep="_")))
}

