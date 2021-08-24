"""
indel_calls_redone.py - redo sketchy indel calls
to try and resolve which is the mutant sample
"""

import sys
import csv
import subprocess
from glob import glob
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='', 
        usage='python indel_calls_redone.py [options]')

    parser.add_argument('-f', '--mut_file', required=True,
        type=str, help='Mut table with indels')
    parser.add_argument('-d', '--bam_dir', required=True,
        type=str, help='Dir containing sample bams')
    parser.add_argument('-s', '--sample_bam_dirs', required=True,
        type=str, nargs='+', help='Dirs containing older bams to include')
    parser.add_argument('-r', '--flank_region', required=False, 
        type=int, default=200, help='Size of surrounding region (default 200 bp flank)')
    parser.add_argument('-o', '--outdir', required=True,
        type=str, help='Dir to write output VCF to')

    args = parser.parse_args()

    return args.mut_file, args.bam_dir, args.sample_bam_dirs, \
        args.flank_region, args.outdir

def haplotype_caller(fname, chrom, pos, bam_dir, 
    sample_bam_dirs, flank_region, outdir):
    """call variants in region surrounding possible indel

    Parameters
    ----------
    fname : str
        fname column in mut table, corresponds to sample
    chrom : str
        chromosome that indel is on
    pos : int
        position of indel in original call
    bam_dir : str
        dir containing bams for samples
    sample_bam_dirs : str
        dirs containing other bams from older samples to add to HC call
    flank_region : int
        size of region to call variants in
    outdir : str
        dir to write output VCF to

    Returns
    -------
    None
    """
    # construct cmd
    cmd = '''
    time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta
    '''

    # add bams and region to cmd
    bams = []
    if not bam_dir.endswith('/'):
        bam_dir += '/'

    sample = fname.replace('_samples', '')
    if not '_4' in sample:
        bams.extend([f"{bam_dir}{sample}_{i}.bam" for i in [0, 5]])
    elif '_4' in sample:
        prefixes = fname.lstrip('DL').split('_')
        bams.append(f"{bam_dir}DL{prefixes[0]}_0.bam")
        bams.append(f"{bam_dir}DL{prefixes[1]}_5.bam")
    for sample_bam_dir in sample_bam_dirs:
        if not sample_bam_dir.endswith('/'):
            sample_bam_dir += '/'
        bams.extend(glob(f'{sample_bam_dir}*bam'))
    for bam_fname in bams:
        cmd += f'-I {bam_fname} '
    cmd += f'-L {chrom}:{pos-flank_region}-{pos+flank_region} '

    # add other fixed parameters
    cmd += '''
    -ploidy 2 --output_mode EMIT_ALL_SITES \
    --heterozygosity 0.02 --indel_heterozygosity 0.002
    '''

    # add output info
    if not outdir.endswith('/'):
        outdir += '/'
    cmd += f'-o {outdir}{sample}_{chrom}_{pos}.vcf'
    cmd = cmd.replace('\n', '')

    cmd_current = [w for w in cmd.split() if w]
    subprocess.check_call(
        cmd_current, stdout=sys.stdout, stderr=sys.stderr)


def parse_muts(mut_file, bam_dir, sample_bam_dirs, flank_region, outdir):
    """parse through muts in mut table and call indels for each

    Parameters
    ----------
    mut_file : str
        path to mut table file containing indels
    bam_dir : str
        dir containing bams for samples
    sample_bam_dirs : str
        other bams from older samples to add to HC call
    flank_region : int
        size of region to call variants in
    outdir : str
        dir to write output VCF to

    Returns
    -------
    None
    """

    with open(mut_file, 'r') as f:
        line_count = sum([1 for line in f])

    with open(mut_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for line in tqdm(reader, desc='indels', total=line_count):
            chrom, pos = line['chrom'], int(line['pos'])
            tqdm.write(str([line['fname'], chrom, pos]))
            haplotype_caller(
                line['fname'], chrom, pos, bam_dir,
                sample_bam_dirs, flank_region, outdir)
            

def main():
    parse_muts(*args())

if __name__ == '__main__':
    main()

        

