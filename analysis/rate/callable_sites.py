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
        type=str, help='Dir with gVCF files')
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

def get_lookup(samples, gvcf_dir, chrom, start, end):
    """create lookup strings for input region

    given a region L, creates a dict of length-L strings
    denoting whether a given site is callable or not for a
    given sample

    Parameters
    ----------
    samples : list
        sample names, corresponding to gVCF names
    gvcf_dir : str
        path to input gVCFs
    chrom : str
        name of desired chromosome
    start : int
        start coordinate of region (1-indexed)
    end : int
        end coordinate of region (1-indexed)

    Returns
    -------
    lookup : dict
        dict containing lookup strings for each sample
    """
    lookup = {}
    tqdm.write(f'[saltMA] creating lookups for {chrom}:{start}-{end}')
    for sample in tqdm(samples):
        fname = f'{gvcf_dir}{sample}.g.vcf.gz'
        reader = VCF(fname)
        lookup[sample] = ''
        current_pos = start + 1
        for record in reader(f'{chrom}:{start}-{end}'):
            if record.POS == current_pos:
                lookup[sample] += '1'
                current_pos += 1
            elif record.POS > current_pos:
                to_add = record.POS - current_pos # goes to 0 if pos = curr + 1
                str_to_add = '0' * to_add
                lookup[sample] += str_to_add
                lookup[sample] += '1'
                current_pos = record.POS + 1
        if current_pos < end: # add 'outer' 0s if stopped short of window end
            to_add = end - current_pos + 1
            str_to_add = '0' * to_add
            lookup[sample] += str_to_add
        elif current_pos == end:
            lookup[sample] += '0'
    return lookup
            

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
            total = lengths[chrom]
            for window_start in tqdm(range(0, total, int(1e5))):
                lookup = {}
                if window_start + 1e5 > total:
                    window_end = total
                else:
                    window_end = window_start + int(1e5) - 1
                lookup = get_lookup(samples, gvcf_dir, chrom, window_start, window_end)
                # print(window_start, window_end)
                # print([[k, len(v)] for k, v in lookup.items()])
                for i in range(window_end - window_start):
                    call_dict = {k: v[i] for k, v in lookup.items()}
                    call_dict['chrom'] = chrom
                    call_dict['pos'] = window_start + i + 1
                    writer.writerow(call_dict)

def main():
    gvcf_dir, out = args()
    if not gvcf_dir.endswith('/'):
        gvcf_dir += '/'
    parse_calls(gvcf_dir, out)

if __name__ == '__main__':
    main()

