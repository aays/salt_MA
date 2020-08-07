'''
mutation_snapshots.py - create images of candidate mutations
via server instance of IGV

takes in 'table format' output file from filter_candidate_muts.py

requires that an Xvfb 'X server' is running to perform all
operations on the server and avoid having to use X11 forwarding,
as well as IGV shell script installation

to run: ./bin/Xvfb :1 -nolisten tcp -fp /path/to/X11/fonts/misc
'''

import os
import sys
import csv
import re
import time
import subprocess
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='create images of candidate mutations with IGV', 
        usage='python3.5 mutation_snapshots.py [options]')

    parser.add_argument('-f', '--fname', required=True,
                        type=str, help='Table format file with mutations')
    parser.add_argument('-i', '--igv', required=True,
                        type=str, help='Path to IGV executable')
    parser.add_argument('-r', '--reference', required=True,
                        type=str, help='Reference fasta')
    parser.add_argument('-b', '--bam_files', required=True,
                        type=str, nargs='+', help='BAM(s) to use')
    parser.add_argument('-l', '--flank_size', required=True,
                        type=int, help='bp of sequence to include flanking mutation')
    parser.add_argument('-o', '--outdir', required=True,
                        type=str, help='Dir to write outfile to')
    parser.add_argument('-s', '--snps_only', required=False,
                        action='store_true', help='Output SNPs only [optional]')
    parser.add_argument('-n', '--indels_only', required=False,
                        action='store_true', help='Output indels only [optional]')

    args = parser.parse_args()

    return args.fname, args.igv, args.reference, args.bam_files, args.flank_size, \
        args.outdir, args.snps_only, args.indels_only


def create_batch_script(fname, reference, bam_files, flank_size,
                        outdir, snps_only=None, indels_only=None):
    '''
    writes temporary batch script for IGV - this script is later
    deleted after use in generate_snapshots()

    see documentation for generate_snapshots() for arg descriptions
    '''
    fname_temp = os.path.basename(fname).rstrip('.txt') + '.batch.temp'
    with open(fname_temp, 'w') as f:
        f.write('new\n')
        f.write('snapshotDirectory {}\n'.format(outdir))
        f.write('genome {}\n'.format(reference))
        for bam_file in bam_files:
            f.write('load {}\n'.format(bam_file))
        with open(fname, 'r', newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            print('[saltMA] writing variants to batch script...')
            for record in tqdm(reader):
                chrom = record['chrom']
                start = int(record['pos']) - flank_size
                end = int(record['pos']) + flank_size
                try: # make sure we don't eval something fishy
                    pattern = "^\[[0-9\., ]+\]+"
                    match = re.search(pattern, record['gt_quals'])
                    assert match
                    pattern = "^\[['ACGT, ]+\]$"
                    match = re.search(pattern, record['alt'])
                    assert match
                except:
                    print('[saltMA] either a GT qual or an alt column in the input table is bonked')
                    print('[saltMA] The record was {}:{}'.format(record['chrom'], record['pos']))
                    raise ValueError()

                # assign snp or indel
                alt_column = eval(record['alt'])
                max_alt_length = max(set([len(s) for s in alt_column]))
                if max_alt_length > 1 or len(record['ref']) > 1:
                    variant_type = 'indel'
                elif max_alt_length == 1 and len(record['ref']) == 1:
                    variant_type = 'snp'
                else:
                    raise ValueError('[saltMA] ref/alt lengths are being weird...')

                # filter if needed
                if snps_only and variant_type == 'indel':
                    continue
                elif indels_only and variant_type == 'snp':
                    continue
                quals = eval(record['gt_quals'])
                quals = 'GQ' + '_'.join([str(int(q)) for q in quals])
                f.write('goto {}:{}-{}\n'.format(chrom, start, end))
                fname_out = [record['fname'].rstrip('_samples'),
                        record['chrom'], record['pos'], quals, variant_type]
                f.write('snapshot {}_{}_{}_{}_{}.png\n'.format(*fname_out))
        f.write('exit')


def generate_snapshots(fname, igv, reference, bam_files, flank_size, outdir,
                       snps_only=None, indels_only=None):
    '''
    use temporary batch script to run IGV

    fname : str
        name of input table formatted filename
    igv : str
        path to IGV executable (igv.sh)
    reference : str
        FASTA file containing reference genome
    bam_files : list
        list of paths to input bam file(s)
    flank_size : int
        size of flanking sequence surrounding candidate mutation
        to be included
    outdir : str
        dir to write to
    snps_only : bool
        whether to only keep SNPs [optional]
    indels_only : bool
        whether to only keep indels [optional]
    '''
    print('[saltMA] mutation screenshots for {}'.format(fname))
    print('[saltMA] generating batch script...')
    create_batch_script(fname, reference, bam_files, flank_size, outdir,
            snps_only=snps_only, indels_only=indels_only)
    fname_batch = os.path.basename(fname).rstrip('.txt') + '.batch.temp'
    try:
        assert os.path.isfile(fname_batch)
    except:
        print('[saltMA] ERROR: the batch script does not exist')
        print('[saltMA] something has gone terribly wrong. Exiting...')
        raise FileNotFoundError()

    print('[saltMA] running IGV with batch script...')
    if not igv.startswith('./'):
        igv = './' + igv
    t0 = time.time()
    proc = subprocess.Popen([igv, '--batch', fname_batch],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()

    t1 = time.time()
    if t1 - t0 < 5:
        print('[saltMA] IGV finished awful quick -')
        print('[saltMA] chances are you do not have an Xvfb server running.')
        print('[saltMA] double check to see output has been generated')
    print('[saltMA] screenshotting complete - deleting temp batch script...')
    os.remove(fname_batch)
    print('[saltMA] results written to {}. exiting...'.format(outdir))


def main():
    generate_snapshots(*args())

if __name__ == '__main__':
    main()

        

