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
        usage='python filter_candidate_muts.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF to filter (.vcf.gz)')
    parser.add_argument('-g', '--gq', required=False, default=30,
        type=int, help='GQ threshold (default 30)')
    parser.add_argument('-f', '--out_format', required=False, default='vcf',
        type=str, help='VCF format or tabular format? ([table|vcf] - default VCF)') 
    parser.add_argument('-t', '--vcf_type', required=True,
        type=str, help="VCF type ('combined' or 'pairs')")
    parser.add_argument('-l', '--verbose_level', required=False, default=2,
        type=int, help='verbose level (0=none, 1=low, 2=all)')
    parser.add_argument('-p', '--purity_filter', required=False, 
        action='store_true', help='enable alt allele purity filter')
    parser.add_argument('-s', '--single_sample', required=False,
        type=str, help='which sample to compare 0 or 5 against if only one') 
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.vcf, args.gq, args.out_format, args.vcf_type, \
        args.verbose_level, args.purity_filter, args.single_sample, args.out

def check_record(record, vcf_type, sample_lookup, gq, 
    single_sample, purity_filter=False):
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
    single_sample : str
        if enabled, will assume there is only one sample of interest
        in the VCF (e.g. just one of 0 and 5) and use provided sample as
        the 'pair'

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
            if not single_sample:
                sample_0_gt_count = record.genotypes.count(record.genotypes[idx_0])
                sample_5_gt_count = record.genotypes.count(record.genotypes[idx_5])
                if min(sample_0_gt_count, sample_5_gt_count) > 1:
                    return 'doublemut'
            elif single_sample:
                sample_5_gt_count = record.genotypes.count(record.genotypes[idx_5])
                if sample_5_gt_count > 1: # only check 0/5 sample
                    return 'doublemut'
            if len(set(record.gt_bases)) == 1: # skip if all samples have same call
                return False
        elif record.genotypes[idx_0] == record.genotypes[idx_5]:
            return False
    """
    elif single_sample:
        if vcf_type == 'pairs':
            raise Exception('pairs mode selected with single sample - exiting')
        sample_gt_count = record.genotypes.count(record.genotypes[idx])
        if sample_gt_count > 2:
            return 'doublemut'
    """

    # only runs if all checks passed, unless double mutation
    return True

def parse_records(vcf, gq, out_format, vcf_type, verbose_level, 
    purity_filter, single_sample, out):
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
    single_sample : bool
        if enabled, will assume there is only one sample of interest
        in the VCF (e.g. just one of 0 and 5) - also assumes 'combined' mode
        
    out : str
        Name of file to write to

    Returns
    -------
    None
        Writes to specified file.
    '''
    print(f'the verbose level is {verbose_level}')
    vcf_in = VCF(vcf)
    sample_names = vcf_in.samples
    pair_sample_names = sorted([item for item in sample_names if item.endswith('_0')
            or item.endswith('_5')])
    if single_sample:
        try:
            assert len(pair_sample_names) == 1
        except AssertionError as e:
            print(f'[saltMA] ERROR: samples bonked - {pair_sample_names}')
            print('[saltMA] Exiting...')
            sys.exit()


    # get samples
    try:
        if not single_sample:
            sample_lookup = sample_names.index(pair_sample_names[0]), \
                sample_names.index(pair_sample_names[1])
        elif single_sample:
            shortlist = [s for s in sample_names if single_sample in s
                        and s != pair_sample_names[0]]
            sample_lookup = sample_names.index(shortlist[0]), \
                sample_names.index(pair_sample_names[0])
            print(f'[saltMA] selected samples are {pair_sample_names[0]} and {shortlist[0]}')
    except IndexError as e:
        print('[saltMA] ERROR: Samples seem incorrect. '
              'Ensure you have a 0 and 5 sample in the VCF, '
              'unless --single_sample has been selected.')
        print('[saltMA] Exiting...')
        sys.exit()

    print(f'[saltMA] initiating filtering for {os.path.basename(vcf)}...')
    counter = 0
    total_count = 0
    doublemut_count = 0
    prev_chr = None

    if out_format == 'vcf':
        outfile = Writer(out, vcf_in)
        outfile.write_header()
    elif out_format == 'table':
        f = open(out, 'w')
        header_string = '\t'.join(['fname', 'chrom', 'pos', 'ref', 'alt',
        'gt_bases', 'gt_quals', 'gt_depths'])
        f.write(header_string + '\n')
        basename = os.path.basename(vcf).replace('.vcf.gz', '')

    for record in tqdm(vcf_in):
        total_count += 1
        check = check_record(record, vcf_type, sample_lookup, 
            gq=gq, purity_filter=purity_filter, single_sample=single_sample)
        if check == 'doublemut':
            if verbose_level == 2:
                tqdm.write(f'[saltMA] doublemut at {record.CHROM}:{record.POS}')
                doublemut_count += 1
                continue
            else:
                doublemut_count += 1
                continue
        elif check:
            counter += 1
            if verbose_level == 1:
                if not prev_chr:
                    prev_chr = record.CHROM
                    tqdm.write(f'[saltMA] first chrom with detected mut is {prev_chr}')
                    continue
                elif prev_chr != record.CHROM:
                    tqdm.write(f'[saltMA] {prev_chr} completed.')
                    tqdm.write(f'[saltMA] current count is {counter}')
                    prev_chr = record.CHROM
            elif verbose_level == 2:
                tqdm.write(f'[saltMA] candidate mut found at {record.__repr__()}')
                tqdm.write(f'[saltMA] current count is {counter}')
                tqdm.write(f'[saltMA] doublemut count is {doublemut_count}')
            if out_format == 'vcf':
                outfile.write_record(record)
            elif out_format == 'table':
                out_string = '\t'.join([basename, record.CHROM, str(record.POS),
                        record.REF, str(record.ALT), str(list(record.gt_bases)),
                        str(list(record.gt_quals)), str(list(record.gt_depths)) + '\n'])
                f.write(out_string)

    if out_format == 'table':
        f.close()

    print(f'[saltMA] completed search for {os.path.basename(vcf)}')
    print(f'[saltMA] found {counter} matches over {total_count} sites.')


def main():
    parse_records(*args())

if __name__ == '__main__':
    main()

        

