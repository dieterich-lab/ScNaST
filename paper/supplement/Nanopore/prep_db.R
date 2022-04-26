#! /usr/bin/env Rscript

## module load R/4.1.1_deb10

# for data import
library(AnnotationHub)
library(ensembldb)
library(zellkonverter)
library(tidyverse)
library(SummarizedExperiment)


parent <- '/prj/Florian_Leuschner_spatial/analysis/Nanopore'
raw_file <- 'concatenated_raw_scnast.h5ad'

# annotations
config <- yaml::yaml.load_file(file.path(parent, 'ScNaST', 'workflow', 'config.yaml', fsep=.Platform$file.sep))
gtf <- file.path(config$strg_base, 'gffcompare', 'gffcmp.multi_exons.annotated.gtf')
db_file <- gsub(pattern = "\\gtf$", "sqlite", basename(gtf))

## Data pre-processing and wrangling

# start by reading the AnnData to a SingleCellExperiment object
sce <- readH5AD(file.path(parent, 'Nanopore', 'data', raw_file))

# first generate custom database from our transcriptome annotation 
# concatenated var contains duplicated ref_id and ref_gene_id...
var <- data.frame(rowData(sce)) %>% rownames_to_column() %>% select(c(1, 4:6, 8:9)) %>% rename(ref_id=ref_id.A, ref_gene_id=ref_gene_id.A)
var[] <- lapply(var, as.character)
var$ref_id <- ifelse(is.na(var$ref_id), var$qry_id, var$ref_id)
var$ref_gene_id <- ifelse(is.na(var$ref_gene_id), var$qry_gene_id, var$ref_gene_id)
var$extra <- paste(var$ref_id, var$ref_gene_id, sep=":")


txdb <- makeTxDbFromGFF(gtf) 

# dump everything and work component-wise
txdb_dump <- as.list(txdb)

# we cannot add new columns...
# so we just replace tx_name, and use tx_type to add ref_id:ref_gene_id
cols <- colnames(txdb_dump$transcripts)

txdb_dump$transcripts <- merge(txdb_dump$transcripts, var[,c('qry_id', 'rowname', 'extra', 'ref_gene_id', 'qry_gene_id')], by.x='tx_name', by.y='qry_id', all.x=T)
txdb_dump$transcripts <- txdb_dump$transcripts[order(txdb_dump$transcripts$tx_id),]
rownames(txdb_dump$transcripts) <- 1:nrow(txdb_dump$transcripts)
# replace fields 
txdb_dump$transcripts$rowname <- ifelse(is.na(txdb_dump$transcripts$rowname), txdb_dump$transcripts$tx_name, txdb_dump$transcripts$rowname)
txdb_dump$transcripts$tx_name <- txdb_dump$transcripts$rowname
txdb_dump$transcripts$extra <- ifelse(is.na(txdb_dump$transcripts$extra), txdb_dump$transcripts$tx_type, txdb_dump$transcripts$extra)
txdb_dump$transcripts$tx_type <- txdb_dump$transcripts$extra
# add gene info
txdb_dump$transcripts$gene_id <- txdb_dump$transcripts$ref_gene_id
txdb_dump$transcripts$gene_id <- ifelse(is.na(txdb_dump$transcripts$gene_id), txdb_dump$genes$gene_id, txdb_dump$transcripts$gene_id)
txdb_dump$transcripts <- txdb_dump$transcripts[,c(cols, 'gene_id')]
# genes should not be supplied if transcripts has a gene_id column
txdb_dump$genes <- NULL

txdb <- do.call(makeTxDb, txdb_dump)


saveDb(txdb, file=file.path(parent, 'Nanopore', 'dtu', db_file))

