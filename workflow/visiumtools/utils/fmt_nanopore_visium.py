#! /usr/bin/env python3


"""Helper script to link Visium spatial information and format Nanopore 
   feature-barcode matrices to h5, to use e.g. with Scanpy read_visium or
   Seurat Load10X_Spatial.
"""

import os
import sys
import logging
import argparse
import gzip
import yaml
import tables

from pathlib import Path

import numpy as np
import scipy.io as io

def module_ff(module, path):
    import importlib.util
    spec = importlib.util.spec_from_file_location(module, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

# should be called from visiumtools!
utils = module_ff("utils", "../sc-nanopore-utils/utils/utils.py")

logger = logging.getLogger(__name__)


# feature annotation (gffcompare)
final_names = [
    'var_names', 
    'transcript_id', 
    'cmp_ref', 
    'gene_id', 
    'ref_gene_id', 
    'gene_name', 
    'class_code']


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Link Visium spatial directory, and
                                     format data to h5.""")

    parser.add_argument('config', help="""The config file (yaml).""")

    parser.add_argument('--do-not-append-barcodes', help="""If this flag is set, leave the
                        Nanopore barcodes as is, otherwise '-1' is appended to the barcodes
                        to match the Illumina barcodes e.g. tissue_positions_list.csv.""", 
                        action='store_true')
    
    parser.add_argument('--feature-type', help="""feature_type field for h5 file.""", 
                        default='Transcript Expression', type=str)
    parser.add_argument('--genome', help="""genome field for h5 file.""", 
                        default='', type=str)
    
    parser.add_argument('--chemistry', help="""chemistry_description metadata field for h5 file.""", 
                        default='Spatial (Nanopore)', type=str)
    parser.add_argument('--filetype', help="""filetype metadata field for h5 file.""", 
                        default='matrix', type=str)
    parser.add_argument('--ids', help="""library_ids metadata field for h5 file.""", 
                        default='UTF-8', type=str)
    
    utils.add_logging_options(parser) 
    args = parser.parse_args()
    utils.update_logging(args)
    
    msg = "[fmt_nanopore_visium]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    required_keys = [
        'samples',
        'samples_cellranger']
    utils.check_keys_exist(config, required_keys)
    
    appd = '-1'
    if args.do_not_append_barcodes:
        appd = ''

    for sample, bam in config['samples'].items():
        parent = str(Path(bam).parents[1])
        bc_dir = os.path.join(parent, 'analysis', 'feature_bc_matrix')
        space_dir = os.path.join(parent, 'analysis', 'spatial')

        # symlink existing Visium spatial directory
        src = os.path.join(config['samples_cellranger'].get(sample), 'spatial')
        dst = os.path.join(space_dir, 'spatial')
        utils.create_symlink(src, dst, remove=False, create=True)
        
        # format data
        filen = os.path.join(bc_dir, 'matrix.mtx.gz')
        mat = io.mmread(filen).T.tocsr()
        
        filen = os.path.join(bc_dir, 'features_annotated.tsv.gz')
        var = [l.decode('utf8').strip('\n').split('\t')[0] for l in gzip.open(filen, 'rb').readlines()]
        # ids are gene_ids
        # read_visium ignores additional information, so we don't bother to add it
        ids = [l.decode('utf8').strip('\n').split('\t')[3] for l in gzip.open(filen, 'rb').readlines()]
        
        filen = os.path.join(bc_dir, 'barcodes.tsv.gz')
        obs = [l.decode('utf8').strip('\n') + appd for l in gzip.open(filen, 'rb').readlines()]
        
        filen = os.path.join(space_dir, 'filtered_feature_bc_matrix.h5')
        h5file = tables.open_file(filen, mode="w", title="")

        group = h5file.create_group("/", 'matrix', '')
        h5file.create_carray(group, 'barcodes', obj=np.asarray(obs))
        h5file.create_earray(group, 'data', obj=mat.data)
        h5file.create_earray(group, 'indices', obj=mat.indices)
        h5file.create_earray(group, 'indptr', obj=mat.indptr)
        h5file.create_earray(group, 'shape', obj=np.asarray(mat.shape[::-1])) # transpose

        features = h5file.create_group(group, 'features', '')
        h5file.create_array(features, '_all_tag_keys', np.asarray(['genome']), "")
        h5file.create_carray(features, 'feature_type', obj=np.asarray([args.feature_type]*mat.shape[1]))
        h5file.create_carray(features, 'genome', obj=np.asarray([args.genome]*mat.shape[1]))
        h5file.create_carray(features, 'id', obj=np.asarray(ids))
        h5file.create_carray(features, 'name', obj=np.asarray(var))

        h5file.root._v_attrs.chemistry_description = "Spatial 3' v1 (Nanopore)"
        h5file.root._v_attrs.filetype = "matrix"
        h5file.root._v_attrs.library_ids = np.asarray([sample.encode('UTF-8')])

        h5file.close()
        
    
if __name__ == '__main__':
    main()
     
