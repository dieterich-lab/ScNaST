#! /usr/bin/env python3


"""Format data to CellRanger-like output
"""

import os
import sys
import logging
import argparse
import shutil
import gzip

import pandas as pd
import scipy.sparse as sparse
import scipy.io as io

# for Slurm sbatch 
try:
    sys.path.append(os.getcwd())
    import utils.utils as utils
except ModuleNotFoundError:
    import utils as utils
    
logger = logging.getLogger(__name__)

default_num_cpus = 1
default_mem = '80G'


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
                                     description="""Convert to CellRanger-like files.""")

    parser.add_argument('sf', help="""The input salmon sf reformated file from all barcodes.""")

    parser.add_argument('outdir', help="""The output directory.""")

    parser.add_argument('--gtf', help="""The GTF file to annotate features, if present. Note:
                        we use default gffcompare field attributes. Var_names is 
                        transcript_id:gene_name, with gene_name replaced by gene_id if 
                        unassigned.""", default=None)
    
    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)
    
    msg = "[reformat]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    # generate standard CellRanger-like output
    data = pd.read_csv(args.sf,
                       sep='\t')
    trx_names = data.Name
    idx = 1
    cell_names = list(data)[idx:]
    if trx_names is not None:
        with gzip.open(os.path.join(args.outdir, "features.tsv.gz"), "wt") as handle:
            for name in trx_names:
                handle.write("{}\n".format(name))
    if cell_names is not None:
        with gzip.open(os.path.join(args.outdir, "barcodes.tsv.gz"), "wt") as handle:
            for name in cell_names:
                handle.write("{}\n".format(name))

    data = sparse.coo_matrix(data.iloc[:, idx:])
    io.mmwrite(os.path.join(args.outdir, "matrix.mtx"), data) 
    with open(os.path.join(args.outdir, 'matrix.mtx'),'rb') as mtx_in:
        with gzip.open(os.path.join(args.outdir, 'matrix.mtx.gz'),'wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(os.path.join(args.outdir, 'matrix.mtx')) 
    
    # re-write features to include gene annotation and class codes (gffcompare)
    if args.gtf is not None:
        msg = "Annotating features!"
        logger.info(msg)
        
        gtf = pd.read_csv(args.gtf, 
                        sep='\t', 
                        header=None,
                        names=utils.gtf_field_names)
        gtf = utils.parse_all_gtf_attributes(gtf, num_cpus=args.num_cpus)
        gtf = gtf[gtf.feature=='transcript']
        gtf.gene_name.fillna(gtf.cmp_ref_gene, inplace=True)
        
        
        df = pd.read_csv(os.path.join(args.outdir, "features.tsv.gz"), 
                         sep='\t', 
                         header=None, 
                         names=['transcript_id'])
        df = pd.merge(df, gtf[final_names[1:]], 
                      on='transcript_id', 
                      how='left')
        df.gene_name.fillna(df.gene_id, inplace=True)
        df['var_names'] = df['transcript_id'] + ':' + df['gene_name']
        df = df[final_names]
        df.to_csv(os.path.join(args.outdir, "features_annotated.tsv.gz"),
                  sep='\t', 
                  compression='gzip', 
                  index=False, 
                  header=False)
    
if __name__ == '__main__':
    main()
    
