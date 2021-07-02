"""
callable_sites_degeneracy.py - get number of synonymous
and nonsynonymous callable sites

usage:
python callable_sites_degeneracy.py \
--callables_table path/to/callables_table.tsv.gz \
--annotation_table path/to/ant_table.txt.gz \
--outname path/to/outfile.tsv

both the callables file and annotation file need to be
bgzipped and tabixed
"""

import csv

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


def create_degen_lookup(annotation_table, scaffold) -> str:
    """create degeneracy lookup for a given chromosome

    I bit the bullet here and did the whole chromosome -
    the max lookup string length is ~9e8, so be careful
    where you run this I suppose

    given a chromosome of size L, returns a length L
    string with a degeneracy value [0234] or a missing
    value [.]

    if sites have multiple degen values - if same value,
    will reduce to that value, but if there are differing
    values that site will just be listed as missing

    Parameters
    ----------
    annotation_table : str
        path to annotation table
    scaffold : str
        scaffold of interest

    Returns
    -------
    lookup : str
        length-L string listing degeneracy of scaffold

    """
    ant_reader = pysam.TabixFile(annotation_table)
    ant_header = ant_reader.header[-1].split('\t')
    degen_idx = ant_header.index('degen')

    lookup = ''

    for record in tqdm(ant_reader.fetch(scaffold), desc=f'{scaffold} lookup'):
        degen = record.split('\t')[degen_idx]
        if len(degen) > 1:
            if degen[0] == degen[-1]:
                degen = degen[0]
            else:
                degen = '.'
        lookup += degen

    return lookup


def count_sites(callables_table, annotation_table, outname):
    """main function to count callable sites by degeneracy

    does one scaffold at a time - uses lookups generated via 
    create_degen_lookup and iterates through callables table
    to fill a nested dict of counts

    values from each scaffold across all samples written once each 
    scaffold is complete - uses csv.DictWriter to create outfile

    both callables table and ant file must be tabixed

    Parameters
    ----------

    callables_table : str
        path to callables table generated by callable_sites.py
    annotation_table : str
        path to annotation table
    outname : str
        file to write to

    Returns
    -------
    None

    """

    callables_reader = pysam.TabixFile(callables_table)

    scaffolds = callables_reader.contigs
    header = callables_reader.header[0].split('\t')
    samples = [colname for colname in header
               if colname.endswith('0') or colname.endswith('5')]

    # prep writer
    with open(outname, 'w') as f:
        fieldnames = ['sample', 'scaffold', 'fold0', 'fold2',
                      'fold3', 'fold4', 'nonsyn', 'syn']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for scaffold in tqdm(scaffolds, desc='scaffolds'):
            tqdm.write(f'[saltMA] starting {scaffold}')
            lookup = create_degen_lookup(annotation_table, scaffold)

            # create dict to store counts
            d = {} 
            for sample in samples:
                d[sample] = {'0': 0, '2': 0, '3': 0, '4': 0}

            # iterate through sites in scaffold and count callable sites per sample
            for record in tqdm(callables_reader.fetch(scaffold), desc=scaffold):
                record = record.split('\t')
                degen = lookup[int(record[1])-1]
                if degen == '.':
                    continue

                sample_indices = [i for i, value in enumerate(record) 
                    if value == '1' and i != 1] # second check to prevent pos from being counted!

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

        
