'''
filter_candidate_muts.py - filter candidate mutations in UG VCF output

main criteria:
    1. GQ >= 30
    2. all lines are 'homozygous'
    3. mutated allele is only present in one of 0 or 5
    4. sample with mutated site has <2 reads of non-mut allele
    5. mut allele not present in ancestral lines (TODO)
'''

import os
import sys
from cyvcf2 import VCF
from cyvcf2 import Writer
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='filter SNMs in UG VCF output', 
        usage='python3.5 filter_candidate_muts.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='VCF to filter (.vcf.gz)')
    parser.add_argument('-g', '--gq', required=False, default=30,
                        type=int, help='GQ threshold (default 30)')
    parser.add_argument('-f', '--out_format', required=False, default='vcf',
                        type=str, help='VCF format or tabular format? \
                        ([table|vcf] - default VCF)') 
    parser.add_argument('-t', '--vcf_type', required=True,
                        type=str, help="VCF type ('combined' or 'pairs')")
    parser.add_argument('-l', '--verbose_level', required=False, default=2,
                        type=int, help='verbose level (0=none, 1=low, 2=all)')
    parser.add_argument('-p', '--purity_filter', required=False, action='store_true', 
                        help='enable alt allele purity filter')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.vcf, args.gq, args.out_format, args.vcf_type, \
        args.verbose_level, args.purity_filter, args.out

def check_record(record, vcf_type, sample_lookup, gq, purity_filter=False):
    '''

    Function checks for mutations that fulfill filters:
    1) GQ > 30
    2) no missing genotypes
    3) sample with mutated site has <2 reads of non mut allele
    (purity filter)

    Parameters
    -------
    record : cyvcf2.cyvcf2.Variant
        record to be checked
    vcf_type : str
        [combined|pairs] whether the VCF just contains the 0/5 pairs
        or other samples as well
    gq : int
        GQ threshold - will exclude records where _any_
        call GQ value is below
    purity_filter : bool
        if enabled, will filter out records where mutated site
        has >2 reads of non mut allele

    Returns
    ------
    bool
        True if record passes all filters.

    '''
    # unpack sample indices
    idx_0, idx_5 = sample_lookup

    # check if variant site
    if not len(record.ALT) > 0:
        return False

    # check no missing genotypes
    pair_calls = [record.genotypes[i] for i in [idx_0, idx_5]]
    if any([-1 in call for call in pair_calls]):
        return False

    # check no heterozygous calls
    if record.num_het != 0:
        return False

    # check GQ > threshold
    if not all(record.gt_quals >= gq):
        return False

    # check that sample with mutated site has <2 reads of non mut allele
    # only checked for if purity filter enabled
    if purity_filter:
        if not min(record.gt_depths - record.gt_alt_depths) < 2:
            return False

    # check for different alleles (obv) and no heterozygous calls
    if vcf_type == 'pairs':
        if record.genotypes[idx_0] == record.genotypes[idx_5]:
            return False
        if record.num_het != 0:
            return False
    elif vcf_type == 'combined':
        if record.genotypes[idx_0] != record.genotypes[idx_5]:
            sample_0_gt_count = record.genotypes.count(record_genotypes[idx_0])
            sample_5_gt_count = record.genotypes.count(record_genotypes[idx_5])
            if min(sample_0_gt_count, sample_5_gt_count) > 1:
                return 'doublemut'

    # only runs if all checks passed, unless double mutation
    return True

def parse_records(vcf, gq, out_format, vcf_type, verbose_level, purity_filter, out):
    '''
    Iterates through VCF and writes records passing above
    filters to file.

    Parameters
    -------
    vcf : str
        File path to input VCF
    gq : int
        GQ threshold (records >gq kept)
    verbose_level : int
        How much to print to console.
        If 0 - print no progress info besides tqdm bar
        If 1 - print counter at each chromosome completion
        If 2 - print all candidate information as they are found
    vcf_type : str
        [combined|pairs] whether the VCF just contains the 0/5 pairs
        or other samples as well
    out_format : str
        [table|vcf] - whether to write as new VCF or
        as a tab-separated file
    purity_filter : bool
        if enabled, will filter out records where mutated site
        has >2 reads of non mut allele (see check_record)
    out : str
        Name of file to write to

    Returns
    -------
    None
        Writes to specified file.
    '''
    vcf_in = VCF(vcf)
    sample_names = vcf_in.samples
    pair_sample_names = sorted([item for item in sample_names if item.endswith('_0')
            or item.endswith('_5')])

    # get samples
    try:
        sample_lookup = sample_names.index(pair_sample_names[0]), sample_names.index(pair_sample_names[1])
    except IndexError as e:
        print('[saltMA] ERROR: Samples seem incorrect. '
              'Ensure you have a 0 and 5 sample in the VCF.')
        print('[saltMA] Exiting...')
        sys.exit()

    print('[saltMA] initiating filtering for {}...'.format(os.path.basename(vcf)))
    counter = 0
    total_count = 0
    prev_chr = None

    if out_format == 'vcf':
        outfile = Writer(out, vcf_in)
        outfile.write_header()
    elif out_format == 'table':
        recs = [] # store records in memory - shouldn't be too many
        f = open(out, 'w')
        header_string = '\t'.join(['fname', 'chrom', 'pos', 'ref', 'alt',
        'gt_bases', 'gt_quals', 'gt_depths'])
        f.write(header_string + '\n')
        basename = os.path.basename(vcf).replace('.vcf.gz', '')

    for record in tqdm(vcf_in):
        total_count += 1
        check = check_record(record, vcf_type, sample_lookup, gq=gq, purity_filter=purity_filter)
        if check == 'doublemut':
            tqdm.write('[saltMA] doublemut at {}:{}'.format(record.CHROM, record.POS))
            continue
        elif check:
            counter += 1
            if verbose_level == 1:
                if not prev_chr:
                    prev_chr = record.CHROM
                    continue
                elif prev_chr != record.CHROM:
                    tqdm.write('[saltMA] {} completed.'.format(prev_chr))
                    tqdm.write('[saltMA] current count is {}'.format(counter))
                    prev_chr = record.CHROM
            elif verbose_level == 2:
                tqdm.write('[saltMA] candidate mut found at {}'.format(record.__repr__()))
                tqdm.write('[saltMA] current count is {}'.format(counter))
            if out_format == 'vcf':
                outfile.write_record(record)
            elif out_format == 'table':
                out_string = '\t'.join([basename, record.CHROM, str(record.POS),
                        record.REF, str(record.ALT), str(list(record.gt_bases)),
                        str(list(record.gt_quals)), str(list(record.gt_depths)) + '\n'])
                f.write(out_string)

    if out_format == 'table':
        f.close()

    print('[saltMA] completed search for {}'.format(os.path.basename(vcf)))
    print('[saltMA] found {} matches over {} sites.'.format(counter,
        total_count))


def main():
    parse_records(*args())

if __name__ == '__main__':
    main()

        

