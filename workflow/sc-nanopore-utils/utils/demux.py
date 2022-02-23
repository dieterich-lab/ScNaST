#! /usr/bin/env python3


"""Write BAM files one per barcode using the scNapBar results.
"""

import os
import sys
import logging
import argparse
import shutil

import pandas as pd
import pysam as ps

from contextlib import ExitStack
 
# for Slurm sbatch 
try:
    sys.path.append(os.getcwd())
    import utils.utils as utils
except ModuleNotFoundError:
    import utils as utils

logger = logging.getLogger(__name__)

default_num_cpus = 1
default_mem = '80G'


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Write BAM files one per barcode 
                                     using the scNapBar results.""")

    parser.add_argument('bam', help="""The input long read BAM file (full path).""")

    parser.add_argument('outdir', help="""The output directory (BAM files, transcriptomes, etc.).""")
    
    parser.add_argument('assignment', help="""The scNapBar results, typically 'real.label', with
                        at least 2 columns: read_id, barcode, <score>. No header, ignore comments.""")
    
    parser.add_argument('-s', '--score', help="""Post hoc filtering of scNapBar scores""", 
                        type=int, default=None)
    
    parser.add_argument('--keep-bams', help="""Use this flag to keep BAM files, otherwise they are 
                        removed.""", action='store_true')
    
    parser.add_argument('--overwrite', help="""If this flag is present, existing files 
                        will be overwritten.""", action='store_true')
    
    parser.add_argument('-t', '--tmp', help="""Optional argument: where to write 
                        temporary files. If not specified, programs-specific tmp will be used.""", 
                        type=str, default=None)

    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    args.keep_intermediate_files = False # deletes unsorted bam file
    
    msg = "[demux]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    # requires a lot of memory however...
    bam = ps.AlignmentFile(args.bam, "rb")
    qname_index = ps.IndexedReads(bam)
    qname_index.build()
    
    grouping_key = 'barcode'
    header = ['read_id', 'barcode']
    if args.score:
        header.extend(['score'])
    assignment = pd.read_csv(args.assignment, 
                             sep='\t', 
                             header=None, 
                             names=header, 
                             usecols=range(len(header)),
                             comment='#')
    if args.score:
        assignment = assignment[assignment.score>=args.score]
    barcodes = assignment.groupby(grouping_key)
    for barcode, reads in barcodes:
        filen = os.path.join(args.outdir, f"{barcode}.bam")
        with ExitStack() as stack:
            barcode_bam = stack.enter_context(ps.AlignmentFile(filen, "wb", template=bam))
            for qname in reads.read_id.values:
                alignments = qname_index.find(qname)
                for a in alignments:
                    barcode_bam.write(a)
            
        filen_sorted = os.path.join(args.outdir, f"{barcode}_sorted.bam")
        args.name = str(barcode)
        utils.sort_bam_file(filen, filen_sorted, args)
        utils.index_bam_file(filen_sorted, args)
    
    bam.close()
    
    # cleaning
    if not args.keep_bams:
        msg = "Removing temporary BAM files."
        logger.warning(msg)
        shutil.rmtree(args.outdir)
        


if __name__ == '__main__':
    main()
    
