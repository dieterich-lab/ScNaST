#! /usr/bin/env python3


import os
import sys
import logging
import argparse

import pandas as pd
import numpy as np
import scanpy as sc

import scvi
from scvi.external import RNAStereoscope, SpatialStereoscope

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Cell type deconvolution with Stereoscope""")

    parser.add_argument('dirloc', help="Output directory")
    parser.add_argument('output', help="Output file name")
    parser.add_argument('spatial', help="Input file name (pre-processed Visium data)")
    parser.add_argument('normalization', help="Gene length normalization (path to file)")
    args = parser.parse_args()
    
    # contains some hard coded fields!
    # read in single-cell data from Tabula Muris
    
    logger.info("Reading Tabula Muris")
    
    TM_droplet = sc.read("TM_droplet.h5ad",
                         backup_url="https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_mat.h5ad")
    TM_droplet.obs = pd.read_csv("https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_metadata.csv")
    TM_facs = sc.read("TM_facs.h5ad",
                      backup_url="https://s3.amazonaws.com/czbiohub-tabula-muris/TM_facs_mat.h5ad")
    TM_facs.obs = pd.read_csv("https://s3.amazonaws.com/czbiohub-tabula-muris/TM_facs_metadata.csv")

    TM_droplet = TM_droplet[(TM_droplet.obs.tissue == "Heart_and_Aorta") & (~TM_droplet.obs.cell_ontology_class.isna())].copy()
    TM_facs = TM_facs[(TM_facs.obs.tissue == "Heart") & (~TM_facs.obs.cell_ontology_class.isna())].copy()

    TM_droplet.obs["library"] = "10x"
    TM_facs.obs["library"] = "SS2"
    
    # for Smart-Seq2, apply gene-length normalization 
    logger.info("Applying gene length normalization for Smart-Seq2")
    gene_len = pd.read_csv(args.normalization,
                           delimiter=" ",
                           header=None,
                           index_col=0)
    gene_len = gene_len.reindex(TM_facs.var.index).dropna()
    TM_facs = TM_facs[:, gene_len.index]
    assert (TM_facs.var.index == gene_len.index).sum() == TM_facs.shape[1]
    TM_facs.X = TM_facs.X / gene_len[1].values * np.median(gene_len[1].values)
    # round to integer
    TM_facs.X = np.rint(TM_facs.X)
    
    # concatenate dataset
    # no need to integrate them, we use the counts and cell type labels
    TM = TM_droplet.concatenate(TM_facs)
    TM.layers["counts"] = TM.X.copy()
    sc.pp.normalize_total(TM, 
                          target_sum=1e4, 
                          exclude_highly_expressed=True)
    sc.pp.log1p(TM)
    TM.raw = TM 
    sc.pp.highly_variable_genes(TM, 
                                layer="counts", 
                                batch_key="library",
                                subset=True)
    
    logger.info(f'Cell types Tabula Muris (Heart) {TM.obs["cell_ontology_class"].value_counts()}')
    
    
    # now read in our visium data
    # this has already been pre-processed, and we can read it as h5ad (contains spatial info)
    spatial = sc.read_h5ad(args.spatial)
    
    
    # now learn cell-type expression from scRNA-seq
    # first filter genes to be the same on the spatial data
    intersect = np.intersect1d(TM.var_names, spatial.var_names)
    adata_spatial = spatial[:, intersect].copy()
    adata = TM[:, intersect].copy()
    
    logger.info(f"Using {len(intersect)} common genes!")
    
    RNAStereoscope.setup_anndata(adata, 
                                 layer="counts", 
                                 labels_key="cell_ontology_class",
                                 batch_key='library')
    sc_model = RNAStereoscope(adata)
    sc_model.train(max_epochs=200)
    sc_model.save(os.path.join(args.dirloc, 'scmodel'), overwrite=True)

    # now infer proportionm for spatial data
    
    # ISSUE: SpatialStereoscope does not seem to register layers...
    X = adata_spatial.X
    adata_spatial.X = adata_spatial.layers['counts']
    
    SpatialStereoscope.setup_anndata(adata_spatial, 
                                    # layer="counts", 
                                     batch_key='library_id')
    spatial_model = SpatialStereoscope.from_rna_model(adata_spatial, sc_model)
    spatial_model.train(max_epochs=2000)
    spatial_model.save(os.path.join(args.dirloc, 'spatialmodel'), overwrite=True)
    
    # infer proportions on each Visium spot for each cell type in the single cell reference dataset
    adata_spatial.obsm["deconvolution"] = spatial_model.get_proportions()
    # also copy as single field in the anndata for visualization
    for ct in adata_spatial.obsm["deconvolution"].columns:
        adata_spatial.obs[ct] = adata_spatial.obsm["deconvolution"][ct]

    # see above
    adata_spatial.X = X
    
    filen = os.path.join(args.dirloc, f"{args.output}_stereoscope.h5ad")
    logger.info(f"Writing to {filen}")
    adata_spatial.write_h5ad(filen)
    

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, filename="logfile.txt", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    main()
