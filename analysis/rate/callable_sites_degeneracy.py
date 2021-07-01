"""
callable_sites_degeneracy.py - get number of synonymous
and nonsynonymous callable sites
"""

import csv
import ant # custom annotation table parser

import pysam
import argparse

from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='count S and NS callable sites', 
        usage='python callable_sites_degeneracy.py [options]')

    parser.add_argument('-c', '--callables_table', required=True,
        type=str, help='Path to callables table')
    parser.add_argument('-t', '--annotation_table', required=True,
        type=str, help='Path to annotation table (.txt.gz)')
    parser.add_argument('-o', '--outname', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.callables_table, args.annotation_table, args.outname

def count_sites(callables_table, annotation_table, outname):

    callables_reader = pysam.TabixFile(callables_table)
    # ant_reader = ant.Reader(annotation_table)
    ant_reader = pysam.TabixFile(annotation_table)

    scaffolds = callables_reader.contigs
    header = callables_reader.header[0].split('\t')
    samples = [colname for colname in header
               if colname.endswith('0') or colname.endswith('5')]

    pos_idx = ant_reader.header[-1].split('\t').index('position')
    degen_idx = ant_reader.header[-1].split('\t').index('degen')

    # prep writer
    with open(outname, 'w') as f:
        fieldnames = ['sample', 'scaffold', 'fold0', 'fold2',
                      'fold3', 'fold4', 'nonsyn', 'syn']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for scaffold in tqdm(scaffolds, desc='scaffolds'):
            tqdm.write(f'[saltMA] starting {scaffold}')

            # create dict to store counts
            d = {} 
            for sample in samples:
                d[sample] = {'0': 0, '2': 0, '3': 0, '4': 0}

            # iterate through sites in scaffold
            for record in tqdm(ant_reader.fetch(scaffold)):
                record = record.split('\t')
                pos = int(record[pos_idx])
                degen = record[degen_idx]
                if degen == '.':
                    continue
                if len(degen) > 1:
                    if degen[0] != degen[-1]:
                        continue
                    else:
                        degen = degen[0]

                # get indices of samples with 1s
                callables_iter = callables_reader.fetch(scaffold, pos-1, pos)
                rec = [r for r in callables_iter][0]
                sample_indices = [
                    i for i, value in enumerate(rec.split('\t')) if value == '1']

                # increment counts for those samples based off callable file header
                for idx in sample_indices:
                    d[header[idx]][degen] += 1

            # assemble dict to write to file
            out_dict = {}
            for sample in samples:
                out_dict['sample'], out_dict['scaffold'] = sample, scaffold
                out_dict['fold0'] = d[sample]['0']
                out_dict['fold2'] = d[sample]['2']
                out_dict['fold3'] = d[sample]['3']
                out_dict['fold4'] = d[sample]['4']
                out_dict['nonsyn'] = d[sample]['0'] + (d[sample]['3'] * 1/3) + \
                  (d[sample]['2'] * 2/3)
                out_dict['syn'] = d[sample]['4'] + (d[sample]['3'] * 2/3) + \
                    (d[sample]['2'] * 1/3)
                writer.writerow(out_dict)
        
def main():
    count_sites(*args())    

if __name__ == '__main__':
    main()

        

