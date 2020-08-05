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

    return args.filename, args.igv, args.reference, args.bam_files, args.flank_size, \
        args.outdir, args.snps_only, args.indels_only


def create_batch_script(fname, reference, bam_files, flank_size,
        outdir, snps_only=None, indels_only=None):
    '''
    writes temporary batch script for IGV - will need to be deleted separately

    '''
    fname_temp = os.path.basename(fname).rstrip('.txt') + '.batch.temp'
    with open(fname_temp, 'w') as f:
        f.write('new\n')
        f.write('snapshotDirectory {}\n',format(outdir))
        f.write('genome {}\n'.format(reference))
        for bam_file in bam_files:
            f.write('load {}\n'.format(bam_file))
        with open(fname, 'r', newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for record in reader:
                chrom = record['chrom']
                start = int(record['pos']) - flank_size
                end = int(record['pos']) + flank_size
                if snps_only or indels_only:
                    try:
                        pattern = "^\[['ACGT, ]+\]$"
                        match = re.search(pattern, record['alt'])
                        assert match
                    except:
                        print('[saltMA] An alt column in the input table is bonked')
                        print('[saltMA] The record was {}:{}'.format(record['chrom'], record['pos']))
                        print('[saltMA] The offending ALT column was {}'.format(record['alt']))
                        raise ValueError()
                    alt_column = eval(record['alt'])
                    max_alt_length = max(set([len(s) for s in alt_column]))
                    if max_alt_length > 1 and snps_only:
                        continue
                    elif max_alt_length == 1 and indels_only:
                        continue
                    else:
                        pass
                f.write('goto {}:{}-{}\n'.format(chrom, start, end))
                fname_out = [record['fname'].rstrip('_samples'),
                        record['chrom'], record['pos']]
                f.write('snapshot {}_{}_{}.png\n'.format(*fname_out))
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
        print('[saltMA] ERROR: the batch script has not been written')
        print('[saltMA] Something has gone terribly wrong. Exiting...')
        raise FileNotFoundError()

    print('[saltMA] running IGV with batch script...')
    print('[saltMA] you will be hearing from IGV instead of me for a bit.')
    if not igv.startswith('./'):
        igv = './' + igv
    subprocess.Popen('{} --batch {}'.format(igv, fname_batch), 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print('[saltMA] back to me - if this message showed up immediately...')
    print('[saltMA] ...chances are you do not have an Xvfb server running.')
    print('[saltMA] deleting temp batch script...')
    os.remove(fname_batch)
    print('[saltMA] done! results written to {}. exiting...'.format(outdir))


def main():
    generate_snapshots(*args())

if __name__ == '__main__':
    main()

        

