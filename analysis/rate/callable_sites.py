'''
callable_sites.py - create giant table of callable sites
from gvcf files
'''

import os
import re
import csv
import argparse
import glob
from tqdm import tqdm
from cyvcf2 import VCF

def args():
    parser = argparse.ArgumentParser(
        description='create callable site table', 
        usage='python callable_sites.py [options]')

    parser.add_argument('-d', '--gvcf_dir', required=True,
        type=str, help='Dir with GVCF files')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.gvcf_dir, args.out

def get_chr_lengths(vcf_reader) -> dict:
    """parse VCF header for chrom lengths

    uses given cyvcf2 VCF object to create a lookup
    dict with contig lengths - works with input header

    Parameters
    ----------
    vcf_reader : cyvcf2.cyvcf2.VCF
        cyvcf2 VCF object to use

    Returns
    -------
    lengths : dict
        dict with contig names as keys and lengths as values
    """
    raw_header = vcf_reader.raw_header
    regions = ['chromosome_[0-9]+', 'scaffold_[0-9]+', 'cpDNA', 'mtDNA']
    keyvals = []
    for region in regions:
        pattern = f'({region},length=[0-9]+)'
        keyvals.extend([m.group() for m in re.finditer(pattern, raw_header)])
    keyvals = [l.split(',length=') for l in keyvals]
    lengths = {l[0]: int(l[1]) for l in keyvals}
    return lengths

def get_site_calls(samples, gvcf_dir, chrom, pos) -> dict:
    """get called status of each sample at input site

    checks VCFs across all given samples at a given site
    to see whether calls exist

    returns a lookup with 1 as value if site called and 0 if not

    Parameters
    ----------
    samples : list
        list of sample names (same as VCF filename prefixes)
    gvcf_dir : str
        path to input gVCFs
    chrom : str
        contig of interest
    pos : int
        position on contig

    Returns
    -------
    call_dict : dict
        dict with call info for input to csv.DictWriter object
    """
    call_dict = {'chrom': chrom, 'pos': pos}
    for sample in samples:
        fname = f'{gvcf_dir}{sample}.g.vcf.gz'
        reader = VCF(fname)
        records = [rec for rec in reader(f'{chrom}:{pos}-{pos}')]
        if len(records) == 1:
            record = records[0]
        else:
            raise Exception(f'record bonked, {sample} {chrom} {pos}')
        # cyvcf2 will return the nearest record rounding down...
        # ... if no record is found in the desired range
        if record.POS != pos:
            call_dict[sample] = 0
            continue # no call at this site
        else:
            call_dict[sample] = 1
    return call_dict

def parse_calls(gvcf_dir, out):
    """iterate through all sites and write per-sample called status to file

    main parsing function that iterates through all sites
    across all given gVCFs

    calls on earlier functions to get sample names and contig lengths
    before iterating through all possible sites and writing to file

    Parameters
    ---------
    gvcf_dir : str
        path to input gVCFs
    out : str
        file to write to

    Returns
    -------
    None
    """
    # get chrs and lengths from first VCF
    vcfs = glob.glob(f'{gvcf_dir}*.g.vcf.gz')
    lengths = get_chr_lengths(VCF(vcfs[0]))
    samples = sorted([os.path.basename(f).rstrip('.g.vcf.gz') for f in vcfs])
    fieldnames = ['chrom', 'pos']
    fieldnames.extend(samples)

    with open(out, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for chrom in lengths: # use dict keys here
            print(f'[saltMA] processing {chrom}')
            for pos in tqdm(range(1, lengths[chrom] + 1)):
                call_dict = get_site_calls(samples, gvcf_dir, chrom, pos)
                writer.writerow(call_dict)

def main():
    gvcf_dir, out = args()
    if not gvcf_dir.endswith('/'):
        gvcf_dir += '/'
    parse_calls(gvcf_dir, out)

if __name__ == '__main__':
    main()

