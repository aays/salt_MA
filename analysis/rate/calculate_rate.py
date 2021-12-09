"""
calculate_rate.py - get mutation rate summaries per sample

if generation_file not provided, will set all generations to 1

otherwise, expects a lookup space separated file with the format

sample_1 30
sample_2 5
sample_3 20
"""

import re
import csv
import subprocess
import tabix
import argparse
from tqdm import tqdm
from cyvcf2 import VCF

def args():
    parser = argparse.ArgumentParser(
        description='get mutation rates for saltMA samples', 
        usage='python calculate_rate.py [options]')

    parser.add_argument('-m', '--mut_table', required=True,
        type=str, help='Mut describer outfile')
    parser.add_argument('-c', '--callables_table', required=True,
        type=str, help='Callables lookup table')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF for chromosome length lookup')
    parser.add_argument('--generation_count', required=False,
        type=int, help='Generation count to be applied to all samples [optional]')
    parser.add_argument('--generation_file', required=False,
        type=str, help='Lookup file of generation times [optional]')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.mut_table, args.callables_table, args.vcf, \
        args.generation_count, args.generation_file, args.out

def get_chr_lengths(vcf_reader) -> dict:
    """parse VCF header for chrom lengths
    uses given cyvcf2 VCF object to create a lookup
    dict with contig lengths - works with input header

    I wish there was a better way to do this that didn't
    involve counting sites in the callables table - this is really
    just because pytabix can't accept just a chromosome name as input

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
    regions = ['chromosome_[0-9]+', 'scaffold_[0-9]+', 'cpDNA', 'mtDNA', 'mtMinus']
    keyvals = []
    for region in regions:
        pattern = f'({region},length=[0-9]+)'
        keyvals.extend([m.group() for m in re.finditer(pattern, raw_header)])
    keyvals = [l.split(',length=') for l in keyvals]
    lengths = {l[0]: int(l[1]) for l in keyvals}
    return lengths


def get_callables_header(callables_table) -> list:
    """get header from callables table

    helper function for get_callable_count

    Parameters
    ----------
    callables_table : str
        path to tabixed callables table file

    Returns
    -------
    header : list
        list containing column names in order
    """
    cmd = f'tabix -H {callables_table}'
    proc = subprocess.run(cmd.split(' '), stdin=subprocess.PIPE,
        stdout=subprocess.PIPE)
    header = proc.stdout.decode('utf-8').lstrip('#').rstrip('\n').split('\t') # lol
    return header

def get_callables_count(callables_table, sample, chrom, start, end) -> int:
    """get count of callable sites for specified region in sample

    uses tabix to extract lookup strings for specified region
    and returns count of callable sites for specified sample

    Parameters
    ----------
    callables_table : str
        path to tabixed callables table file
    sample : str
        name of sample
    chrom : str
        chromosome name
    start : int
        start of region (half open)
    end : int
        end of region (half open)

    Returns
    -------
    callable_count : int
        number of callable sites in region for sample

    """
    # get callable table column idx
    header = get_callables_header(callables_table)
    sample_idx = header.index(sample)

    # pull region from callables table
    reader = tabix.open(callables_table)
    region = reader.query(chrom, start, end)

    # get region lookup string for sample
    lookup_str = ''.join([line[sample_idx] for line in region])
    callable_count = lookup_str.count('1')

    return callable_count


def summarise_muts(mut_table) -> dict:
    """return lookup dict of mut counts from mut describer output

    currently stores all mut counts in memory in nested
    dict structure (samples outer keys, chromosomes
    inner keys)

    Parameters
    ----------
    mut_table : str
        path to mut describer output

    Returns
    -------
    mut_counts : dict
        nested dict containing mutation counts per chrom per samples
    """
    mut_counts = {}
    with open(mut_table, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for record in reader:
            sample = record['mutant_sample']
            chrom = record['chromosome']
            if sample not in mut_counts:
                mut_counts[sample] = {}
            if chrom not in mut_counts[sample]:
                mut_counts[sample][chrom] = 1
            else:
                mut_counts[sample][chrom] += 1
    return mut_counts


def parse_gen_file(generation_file):
    """parse generation time file

    parses space separated generation time file, if provided

    file must be of the format

    sample_1 30
    sample_2 5
    sample_3 20

    space-separated, without a header

    Parameters
    ----------
    generation_file : str
        path to space-separated generation file

    Returns
    -------
    generations : dict
        sample lookup dict with generation values
    """
    generations = {}
    with open(generation_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            sample, gens = line
            gens = int(round(float(gens.rstrip())))
            generations[sample] = gens
    return generations

def get_rates(mut_table, callables_table, vcf, generation_count, generation_file, out):
    """get and report rates per chromosome per sample

    loads mut counts into memory and cross references each sample
    chromosome's callable sites + generations (if provided)
    before computing mut rate and writing to file

    if generations file or value not provided, generations set to 1

    if generations file and count provided, will raise an error

    Parameters
    ----------
    mut_table : str
        path to mut describer output
    callables_table : str
        path to tabixed callables table file
    vcf : str
        any VCF containing chromosome lengths in header
    generation_count : int
        number of generations to set all samples to
    generation_file : str
        path to space-separated generation file
    out : str
        file to write to

    Returns
    -------
    None
    """
    print('[saltMA] creating mutation table lookup')
    mut_counts = summarise_muts(mut_table)
    vcf_reader = VCF(vcf)
    lengths = get_chr_lengths(vcf_reader)
    fieldnames = ['sample', 'chrom', 'mutations', 'callable_sites',
        'generations', 'mut_rate']

    if generation_count and generation_file:
        raise InputError('both gen count and file provided - pick one!')
    elif generation_file and not generation_count:
        print('[saltMA] obtaining generation times')
        generations = parse_gen_file(generation_file)
    elif generation_count and not generation_file:
        generations = generation_count
    else:
        print('[saltMA] no gen file provided, setting gen = 1 for all samples')
        generations = None

    with open(out, 'w') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for sample, sample_dict in tqdm(mut_counts.items()):
            tqdm.write(f'[saltMA] processing {sample}')
            for chrom in tqdm(sample_dict, desc=sample):
                d = {}
                d['sample'] = sample
                d['chrom'] = chrom
                d['callable_sites'] = get_callables_count(
                    callables_table=callables_table, sample=sample, 
                    chrom=chrom, start=1, end=lengths[chrom])
                d['mutations'] = sample_dict[chrom]
                if isinstance(generations, int): # single value provided
                    d['generations'] = generations
                elif isinstance(generations, dict): # lookup table provided
                    d['generations'] = generations[sample]
                else:
                    d['generations'] = 1
                d['mut_rate'] = d['mutations'] / (d['callable_sites'] * d['generations'])
                writer.writerow(d)

def main():
    get_rates(*args())    

if __name__ == '__main__':
    main()

