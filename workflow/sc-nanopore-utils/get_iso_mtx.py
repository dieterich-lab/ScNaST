#! /usr/bin/env python3


"""Performs Nanopore read summarization based on the scNapBar barcode assignment
at the isoform level, using Salmon.
Use conda activate scNapBar followed by module load salmon/1.5.2
"""

import os
import sys
import logging
import argparse
import shutil
import yaml

from pathlib import Path

from paths import PROJECT_ROOT
SCRIPTS_PATH = PROJECT_ROOT / 'scripts'
UTILS_PATH = PROJECT_ROOT / 'utils'

import utils.utils as utils


logger = logging.getLogger(__name__)

default_num_cpus = 6
default_mem = '80G'


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Performs Nanopore read summarization based on 
                                     the scNapBar barcode assignment at isoform level.""")

    parser.add_argument('config', help="""The config file (yaml).""")
    
    parser.add_argument('-s', '--score', help="""Post hoc filtering of scNapBar scores""", type=int, default=None)
    
    parser.add_argument('--overwrite', help="""If this flag is present, existing files 
                        will be overwritten (limited, if directories exist or no, no files check).
                        The last step (feature-barcode matrix) is always executed and overwrites
                        existing files.""", action='store_true')

    parser.add_argument('--keep-all', help="""By default raw files, BAM files (split and minimap)
                        are removed. If this flag is set, keep these files. The salmon output is 
                        not removed.""", action='store_true')

    parser.add_argument('-t', '--tmp', help="""Optional argument: where to write 
        temporary files. If not specified, programs-specific tmp will be used.""", type=str, default=None)
    
    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)
    
    msg = "[get_iso_mtx]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    required_keys = [
        'samples',
        'strg_base',
        'idx_base',
        'idx_name']
    utils.check_keys_exist(config, required_keys)
    
    # construct path to input files/directories
    minimap_index = os.path.join(config['idx_base'], f"{config['idx_name']}.mmi")
    salmon_targets = os.path.join(config['idx_base'], f"{config['idx_name']}.fa")
    salmon_gtf = os.path.join(config['strg_base'], 'gffcompare', 'gffcmp.multi_exons.annotated.gtf')
    
    parent_dir = {}
    result_dir = {}
    for sample, bam in config['samples'].items():
        parent = str(Path(bam).parents[1])
        parent_dir[sample] = parent
        result = str(Path(bam).parents[0])
        result_dir[sample] = result

    # check that all files exist before calling 
    input_files = [minimap_index, salmon_targets, salmon_gtf]
    for sample, bam in config['samples'].items():
        input_files.append(bam)
        input_files.append(os.path.join(result_dir[sample], 'real.label'))
    exist = utils.check_files_exist(input_files, 
                                    raise_on_error=True, 
                                    logger=logger)
    
    # TODO: the overwrite flag is not well handled
    exist_ok = False
    if args.overwrite:
        exist_ok = True
    
    # first split bam files according to scNapBar results
    num_cpus = args.num_cpus
    args.num_cpus = 1
    mapping_dir = {}
    job_ids_split = {}
    all_opts = []
    if args.overwrite:
        all_opts.append('--overwrite')
    if args.tmp is not None:
        all_opts.append(f"-t {args.tmp}")
    if args.score:
        all_opts.append(f"-s {args.score}")
    all_opts_str = ' '.join(all_opts)
    for sample, bam in config['samples'].items():
        mapping_dir[sample] = os.path.join(parent_dir[sample], 'analysis', 'mapping')
        try:
            os.makedirs(mapping_dir[sample], exist_ok=exist_ok)
        except:
            msg = f"Skipping {mapping_dir[sample]} as directory exists!"
            logger.info(msg)
            continue
        cmd = f"{os.path.join(UTILS_PATH, 'demux.py')} {bam} {mapping_dir[sample]} " \
              f"{os.path.join(result_dir[sample], 'real.label')} " \
              f"--keep-bams {all_opts_str}"
        job_id = utils.check_sbatch(cmd, args=args)
        job_ids_split[sample] = job_id
    args.num_cpus = num_cpus

    # convert to fastq
    num_cpus = f"-@{args.num_cpus}"
    raw_dir = {}
    job_ids_bam2fastq = {}
    for sample in config['samples'].keys():
        raw_dir[sample] = os.path.join(parent_dir[sample], 'analysis', 'raw')
        try:
            os.makedirs(raw_dir[sample], exist_ok=exist_ok)
        except:
            msg = f"Skipping {raw_dir[sample]} as directory exists!"
            logger.info(msg)
            continue
        cmd = f"{os.path.join(SCRIPTS_PATH, 'run_bam2fastq')} {mapping_dir[sample]} {raw_dir[sample]} {num_cpus}"
        try:
            job_id = [job_ids_split[sample]]
            job_id_bam2fastq = utils.check_sbatch(cmd, args=args, dependencies=job_id)
            job_ids_bam2fastq[sample] = job_id_bam2fastq
        except:
            job_id_bam2fastq = utils.check_sbatch(cmd, args=args)
            job_ids_bam2fastq[sample] = job_id_bam2fastq

    # transcriptome mapping
    minimap_opts = f'''"{config['minimap_opts']}"'''
    num_cpus = f"-@{args.num_cpus}"
    minimap_dir = {}
    job_ids_minimap = {}
    for sample in config['samples'].keys():
        minimap_dir[sample] = os.path.join(parent_dir[sample], 'analysis', 'minimap')
        try:
            os.makedirs(minimap_dir[sample], exist_ok=exist_ok)
        except:
            msg = f"Skipping {minimap_dir[sample]} as directory exists!"
            logger.info(msg)
            continue
        cmd = f"{os.path.join(SCRIPTS_PATH, 'run_minimap2')} {raw_dir[sample]} {minimap_index} {minimap_dir[sample]} " \
              f"{minimap_opts} {args.num_cpus} {num_cpus} {args.tmp}"
        try:
            job_id = [job_ids_bam2fastq[sample]]
            job_id_minimap = utils.check_sbatch(cmd, args=args, dependencies=job_id)
            job_ids_minimap[sample] = job_id_minimap
        except:
            job_id_minimap = utils.check_sbatch(cmd, args=args)
            job_ids_minimap[sample] = job_id_minimap
            
    # salmon
    if config['salmon_env']:
        if shutil.which('salmon') != config['salmon_env']:
            msg = f"Salmon version is not compatible with that specified in the configuration file! Terminating!"
            logger.error(msg)
    
    salmon_dir = {}
    job_ids_salmon = {}
    salmon_opts = f'''"{config['salmon_opts']}"'''
    for sample in config['samples'].keys():
        salmon_dir[sample] = os.path.join(parent_dir[sample], 'analysis', 'salmon')
        try:
            os.makedirs(salmon_dir[sample], exist_ok=exist_ok)
        except:
            msg = f"Skipping {salmon_dir[sample]} as directory exists!"
            logger.info(msg)
            continue
        cmd = f"{os.path.join(SCRIPTS_PATH, 'run_salmon')} {minimap_dir[sample]} {salmon_dir[sample]} " \
              f"{salmon_targets} {salmon_gtf} {salmon_opts} {args.num_cpus}"
        try:
            job_id = [job_ids_minimap[sample]]
            job_id_salmon = utils.check_sbatch(cmd, args=args, dependencies=job_id)
            job_ids_salmon[sample] = job_id_salmon
        except:
            job_id_salmon = utils.check_sbatch(cmd, args=args)
            job_ids_salmon[sample] = job_id_salmon

    # transform salmon results into standard CellRanger-like output
    job_ids_reformat = {}
    for sample in config['samples'].keys():
        wk_dir = os.path.join(parent_dir[sample], 'analysis', 'wkdir')
        os.makedirs(wk_dir, exist_ok=True)
        cmd = f"{os.path.join(SCRIPTS_PATH, 'reformat_sf')} {salmon_dir[sample]} {wk_dir}"
        try:
            job_id = [job_ids_salmon[sample]]
            job_id_reformat = utils.check_sbatch(cmd, args=args, dependencies=job_id)
            job_ids_reformat[sample] = job_id_reformat
        except:
            job_id_reformat = utils.check_sbatch(cmd, args=args)
            job_ids_reformat[sample] = job_id_reformat

    for sample in config['samples'].keys():
        bc_dir = os.path.join(parent_dir[sample], 'analysis', 'feature_bc_matrix')
        os.makedirs(bc_dir, exist_ok=True)
        cmd = f"{os.path.join(UTILS_PATH, 'reformat.py')} {os.path.join(salmon_dir[sample], 'quant.sf')} {bc_dir}"
        job_id = [job_ids_reformat[sample]]
        utils.check_sbatch(cmd, args=args, dependencies=job_id)

    # enforce cleaning - except salmon output
    if not args.keep_all:
        for sample in config['samples'].keys():
            cmd = f"{os.path.join(SCRIPTS_PATH, 'run_clean')} {mapping_dir[sample]} {raw_dir[sample]} {minimap_dir[sample]}"
            job_id = [job_ids_salmon[sample]]
            utils.check_sbatch(cmd, args=args, dependencies=job_id)
    
    
if __name__ == '__main__':
    main()
    
