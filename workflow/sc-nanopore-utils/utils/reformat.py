#! /usr/bin/env python3

import os
import sys
import logging
import argparse
import shutil
import gzip

import pandas as pd
import scipy.sparse as sparse
import scipy.io as io

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Convert to CellRanger-like files.""")

    parser.add_argument('sf', help="""The input salmon sf reformated file from all barcodes.""")

    parser.add_argument('outdir', help="""The output directory.""")
    
    args = parser.parse_args()
    
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
    
    
if __name__ == '__main__':
    main()
    
