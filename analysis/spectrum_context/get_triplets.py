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
    """
    # TODO: the logic for strain selection will have to be different for salt MA lines 
    # since they're labelled _0 and _5 and that info needs to be preserved
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
    out_dict = {}
    chrom, pos = line_dict['chromosome'], line_dict['position']
    mutant_sample = line_dict['mutant_sample']
    if '_' in mutant_sample:
        out_dict['sample'] = line_dict['mutant_sample'].split('_')[0]
    else:
        out_dict['sample'] = line_dict['mutant_sample']
    out_dict['chromosome'] = chrom
    out_dict['position'] = pos
    out_dict['mutation'] = line_dict['mutation']

    # get VCF records for triplet
    vcf_reader = VCF(vcf)
    records = [rec for rec in vcf_reader(f'{chrom}:{int(pos)-1}-{int(pos)+1}')]
    if any([rec.is_indel for rec in records]):
        return None
    sample_idx = vcf_reader.samples.index(mutant_sample)
    bases = [rec.gt_bases[sample_idx] for rec in records]
    try:
        assert len(bases) == 3
    except:
        print(bases, chrom, pos)
        raise Exception()

    # assemble triplet strings
    triplet_new = ''
    for base in bases:
        if base[0] == base[-1]: # not het, e.g. G/G
            triplet_new += base[0]
        else:
            tqdm.write(f'warning: skipping het site at {mutant_sample} {chrom} {pos}')
            return None
    triplet_old = triplet_new[0] + line_dict['mutation'][0] + triplet_new[-1]

    out_dict['triplet_old'] = triplet_old
    out_dict['triplet_new'] = triplet_new

    return out_dict
    
def parse_file(fname, vcf, most_recent, out):
    # TODO: check samples in VCF and ignore all lines in infile without them
    # it would admittedly be faster to use tabix but only the MA muts are tabixed
    # and I don't have write perms for the adaptation files

    samples = get_samples(vcf, most_recent=most_recent)
    print(samples)

    if fname.endswith('.gz'):
        compressed = True
    else:
        compressed = False
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

def main():
    parse_file(*args())

if __name__ == '__main__':
    main()

