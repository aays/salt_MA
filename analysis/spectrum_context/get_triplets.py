'''
get_triplets.py - get mutation triplet context
'''

import csv
import gzip
import argparse
from tqdm import tqdm
from cyvcf2 import VCF

def args():
    parser = argparse.ArgumentParser(
        description='get mutation triplet context', 
        usage='python get_triplets.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='Mut describer outfile (tsv)')
    parser.add_argument('-v', '--vcf', required=True,
        type=str, help='VCF containing samples of interest')
    parser.add_argument('-r', '--most_recent', required=False, action='store_true', 
        help="whether to isolate 'most recent' sample in VCF")
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.vcf, args.most_recent, args.out

def get_samples(vcf, most_recent=False) -> list:
    """
    get samples in VCF of interest

    if this is an MA VCF, `samples` should be a list of size 1
    since those VCFs contain the same sample at different gen times.
    can get the 'most recent gen' (eg CC2931_13) if most_recent=True

    if this is an adaptation VCF, `samples` will be a list with
    multiple items since these VCFs contain many strains - having
    `most_recent` set to True will break this

    Parameters
    ----------
    vcf : str
        all-sites VCF containing samples of interest
    most_recent : bool
        whether to find and only use the most recent gen for MA

    Returns
    -------
    samples : list
        list of samples in VCF to be used
    """
    vcf_reader = VCF(vcf)
    if not most_recent:
        samples = [s.split('_')[0] for s in vcf_reader.samples]
        samples = list(set(samples))
    elif most_recent:
        sample_base = vcf_reader.samples[0].split('_')[0]
        sample_vals = [int(s.split('_')[1]) for s in vcf_reader.samples]
        max_val = max(sample_vals)
        samples = [f'{sample_base}_{max_val}']
    return samples


def parse_line(line_dict, vcf) -> dict:
    """
    function that parses a single line of mut describer output

    will skip over any het calls in dataset

    Parameters
    ----------
    line_dict : dict
        dict generated by csv.DictReader from mut describer output line
    vcf : str
        all-sites VCF containing samples of interest

    Returns
    -------
    out_dict : dict
        output dict passed to csv.DictWriter with triplet info
    """
    out_dict = {}
    chrom, pos = line_dict['chromosome'], line_dict['position']
    mutant_sample = line_dict['mutant_sample']
    out_dict['sample'] = line_dict['mutant_sample']
    out_dict['chromosome'] = chrom
    out_dict['position'] = pos
    out_dict['mutation'] = line_dict['mutation']

    # get VCF records for triplet
    vcf_reader = VCF(vcf)
    records = [rec for rec in vcf_reader(f'{chrom}:{int(pos)-1}-{int(pos)+1}')]
    if any([rec.is_indel for rec in records]):
        tqdm.write(f'[saltMA] WARNING: skipping {mutant_sample} {chrom} {pos} due to adj indel')
        # tqdm.write(','.join([rec.__repr__() for rec in records]))
        return None
    sample_idx = vcf_reader.samples.index(mutant_sample)
    bases = [rec.gt_bases[sample_idx] for rec in records]
    try:
        assert len(bases) == 3
    except:
        tqdm.write(f'[saltMA] WARNING: skipping unparseable triplet at {mutant_sample} {chrom} {pos}')
        return None # could not parse triplet

    # assemble triplet strings
    triplet_new = ''
    for base in bases:
        if base[0] == base[-1]: # not het, e.g. G/G
            triplet_new += base[0]
        else:
            tqdm.write(f'[saltMA] WARNING: skipping het site at {mutant_sample} {chrom} {pos}')
            return None
    triplet_old = triplet_new[0] + line_dict['mutation'][0] + triplet_new[-1]

    out_dict['triplet_old'] = triplet_old
    out_dict['triplet_new'] = triplet_new

    return out_dict
    
def parse_file(fname, vcf, most_recent, out):
    """
    main function to parse mut describer muts and return triplets

    Parameters
    ----------
    fname : str
        path to mut describer output - can be uncompressed or gzipped
    vcf : str
        all-sites VCF containing samples of interest
    most_recent : bool
        arg passed to get_samples, will only use most recent gen for MA
    out : str
        file to write to

    Returns
    -------
    None
    """
    # it would admittedly be faster to use tabix but only the MA muts are tabixed
    # and I don't have write perms for the adaptation files

    # get samples of interest
    samples = get_samples(vcf, most_recent=most_recent)
    print(f'[saltMA] processing mutation triplets for {samples}')

    # determine infile compression 
    if fname.endswith('.gz'):
        compressed = True
    else:
        compressed = False

    # parse
    with open(out, 'w') as f_out:
        fieldnames = ['sample', 'chromosome', 'position', 
                      'mutation', 'triplet_old', 'triplet_new']
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        with open(fname, 'rb' if compressed else 'rt') as f_in:
            if compressed:
                reader = (line.decode('utf-8') for line in gzip.GzipFile(fileobj=f_in))
            else:
                reader = f_in
            for line in tqdm(csv.DictReader(reader, delimiter='\t')):
                if line['mutant_sample'].split('_')[0] in samples:
                    out_dict = parse_line(line, vcf)
                    if out_dict: # parse_line returns None for indels
                        writer.writerow(out_dict)

    print(f'[saltMA] triplets written to {out}')

def main():
    parse_file(*args())

if __name__ == '__main__':
    main()

