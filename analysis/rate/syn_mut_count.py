"""
syn_mut_count.py - get counts of S and NS mutations
from mut describer output

can optionally take in a gene list and only return sites
in those genes

gene list can be of any tsv format as long as
the gene names themselves are the first column

usage:
python syn_mut_count.py \
--mut_table path/to/mut_describer_file.tsv \
--gene_list path/to/gene_file.tsv \ 
--out path/to/outfile.tsv
"""

import csv
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='get counts of S and NS mutations', 
        usage='python syn_mut_count.py [options]')

    parser.add_argument('-m', '--mut_table', required=True,
        type=str, help='Mut describer output')
    parser.add_argument('-g', '--gene_file', required=False,
        type=str, help='List of genes [optional]')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.mut_table, args.gene_file, args.out

def parse_gene_file(gene_file) -> list:
    """docstring

    Parameters
    ----------
    gene_file : str
        path to file containing list of genes

    Returns
    -------
    gene_list : list
        list of genes
    """
    with open(gene_file, 'r', newline='') as f:
        gene_list = [l[0] for l in csv.reader(f, delimiter='\t')]
    
    # remove any common first column names if they're caught in
    if gene_list[0] in ['gene', 'x']: 
        _ = gene_list.pop(0)

    return gene_list


def parse_muts_gene_list(mut_table, gene_file, out):
    """parse muts given gene list

    intended for use with a provided gene list - outfile
    will look different (will list gene name instead of
    sample/scaffold)

    list of genes can be any tsv as long as the
    gene names are stored in the first column

    Parameters
    ----------
    mut_table : str
        path to mut describer outfile
    gene_file : str
        path to list of genes
    out : str
        file to write to

    Returns
    -------
    None
    """

    print('[saltMA] parsing gene file...')
    gene_list = parse_gene_file(gene_file)

    gene_dict = {}

    # parse mut describer file
    with open(mut_table, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for record in tqdm(reader):
            syn_nonsyn = record['nonsyn_v_syn']
            feature_names = eval(record['feature_names'])
            
            # pass if not syn/nonsyn
            if not syn_nonsyn:
                continue

            # gene list checks
            gene_list_membership = [gene in gene_list for gene in feature_names]
            if not any(gene_list_membership):
                continue
            if all(gene_list_membership):
                raise Exception(f'gene in gene list multiple times: {feature_names}')

            # iterate through feature names and increment counts
            for gene_name in feature_names:
                if gene_name in gene_list:
                    if gene_name not in gene_dict:
                        gene_dict[gene_name] = {'nonsynonymous': 0, 'synonymous': 0}
                    gene_dict[gene_name][syn_nonsyn] += 1
    
    # write to file
    with open(out, 'w', newline='') as f:
        print(f'[saltMA] writing to {out}...')
        fieldnames = ['gene_name', 'nonsyn_muts', 'syn_muts']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for gene_name, count_dict in tqdm(gene_dict.items(), desc='genes w/ muts'):
            d = {}
            d['gene_name'] = gene_name
            d['nonsyn_muts'] = count_dict['nonsynonymous']
            d['syn_muts'] = count_dict['synonymous']
            writer.writerow(d)
        
        # handling genes that have no muts
        gene_list_set = set(gene_list)
        found_genes_set = set(list(gene_dict.keys()))
        assert gene_list_set.issuperset(found_genes_set)
        
        missing_genes = gene_list_set.difference(found_genes_set)
        for gene_name in tqdm(missing_genes, desc='genes w/o muts'):
            d = {}
            d['gene_name'] = gene_name
            d['nonsyn_muts'] = 0
            d['syn_muts'] = 0
            writer.writerow(d)


def summarise_muts_genomewide(mut_table) -> dict:
    """return lookup dict of mut counts from mut describer output

    modified from calculate_rate

    currently stores all mut counts in memory in nested
    dict structure (samples outer keys, chromosomes
    next level keys, syn/nonsyn counts innermost)

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
            syn_nonsyn = record['nonsyn_v_syn']

            if syn_nonsyn:
                if sample not in mut_counts:
                    mut_counts[sample] = {}
                if chrom not in mut_counts[sample]:
                    mut_counts[sample][chrom] = {'nonsynonymous': 0, 'synonymous': 0}
                mut_counts[sample][chrom][syn_nonsyn] += 1

    return mut_counts
    
def parse_muts(mut_table, out):
    """summarise and write genomewide S and NS muts to outfile

    if gene list not provided - will just work off mut describer file and
    summarise counts of S and NS mutations for each scaffold in each sample

    Parameters
    ----------
    mut_table : str
        path to mut describer output
    out : str
        file to write to

    Returns
    -------
    None
    """
    fieldnames = ['sample', 'scaffold', 'nonsyn_muts', 'syn_muts']
    mut_counts = summarise_muts_genomewide(mut_table)

    with open(out, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        for sample, sample_dict in tqdm(mut_counts.items()):
            for chrom in sample_dict:
                nonsyn = sample_dict[chrom]['nonsynonymous']
                syn = sample_dict[chrom]['synonymous']

                d = {}
                d['sample'] = sample
                d['scaffold'] = chrom
                d['nonsyn_muts'] = nonsyn
                d['syn_muts'] = syn
                writer.writerow(d)


def main():
    mut_table, gene_list, out = args()
    if gene_list:
        print('[saltMA] gene list provided')
        parse_muts_gene_list(mut_table, gene_list, out)
    elif not gene_list:
        print('[saltMA] no gene list provided - calculating genomewide counts')
        parse_muts(mut_table, out)
    print('[saltMA] complete.')

if __name__ == '__main__':
    main()

        

