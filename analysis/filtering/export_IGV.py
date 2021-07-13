'''
export_IGV.py - export filtered bams for offline viewing in IGV
'''

import csv
import subprocess
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='export filtered bams for offline viewing in IGV', 
        usage='python export_IGV.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='Tab-separated file containing mutation calls')
    parser.add_argument('-b', '--bam_dir', required=True,
        type=str, help='Dir containing input bams')
    parser.add_argument('-r', '--region_size', required=True,
        type=int, help='Size of region around call to keep')
    parser.add_argument('-o', '--outdir', required=True,
        type=str, help='Dir to save filtered bams to')

    args = parser.parse_args()

    return args.fname, args.bam_dir, args.region_size, args.outdir

def parse_record(line_dict, region_size):
    """
    Parse dict from single line of intervals file

    Parameters
    ----------
    line_dict : collections.OrderedDict
        dict yielded by `csv.DictReader` 
    region_size : int
        desired size of region flanking central pos

    Returns
    -------
    sample : str
        sample of interest
    region : str
        samtools-formatted region 
    """
    sample = line_dict['fname'].replace('_samples', '').replace('_all_salt', '')
    chrom, pos = line_dict['chrom'], int(line_dict['pos'])
    flank_size = int(region_size / 2)
    start, end = pos - flank_size, pos + flank_size 
    region = f'{chrom}:{start}-{end}'
    return sample, region

def export_bams(fname, bam_dir, region_size, outdir):
    """
    Iterate through intervals in interval file and
    generate filtered bams for offline download + viewing

    Parameters
    ----------
    fname : str
        name of tab-separated file containing intervals
    bam_dir : str
        directory containing original bams
    region_size : int
        size of desired region around bam - passed to parse_record
    outdir : str
        dir to save output files to (bams are automatically named)

    Returns
    -------
    None
    """
    print(f'[saltMA] reading intervals from {fname}')
    with open(fname, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        cmd = 'samtools view {} {} -o {}'
        for line in tqdm(reader):
            sample, region = parse_record(line, region_size)
            print(f'[saltMA] processing {sample} at {region}')
            if sample in ['DL41_46', 'DL46_41']:
                first, second = sample.split('_')
                last = 'DL' + second
                bams = [f'{first}_0.bam', f'DL{second}_5.bam']
            elif 'SL26' in sample:
                bams = ['SL26_5.bam', 'SL27_0.bam'] # hardcoding... 
            else:
                bams = [f'{sample}_0.bam', f'{sample}_5.bam']
            pos = int(line['pos'])
            outfile = '{sample}_{chrom}_{pos}_{val}.sam'
            for bam in bams:
                outfile_current = outfile.format(
                    sample=sample, chrom=line['chrom'], pos=pos,
                    val=bam.replace('.bam', '')[-1]) # lol
                print(f'[saltMA] writing to {outfile_current}')
                cmd_format = cmd.format(bam_dir + bam, region, outdir + outfile_current)
                proc = subprocess.run(cmd_format.split(' '))
                proc.check_returncode()

def main():
    fname, bam_dir, region_size, outdir = args()
    if not bam_dir.endswith('/'):
        bam_dir += '/'
    if not outdir.endswith('/'):
        outdir += '/'
    export_bams(fname, bam_dir, region_size, outdir)
    print('[saltMA] done.')

if __name__ == '__main__':
    main()

        

