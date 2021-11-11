"""
clustering_coefficients.py - get clustering coefficients
given a network and a provided list of genes
"""

import csv
import networkx as nx
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='', 
        usage='python script.py [options]')

    parser.add_argument(
        '-g', '--gene_list', required=True,
        type=str, help='File containing list of genes')
    parser.add_argument(
        '-a', '--all_genes', required=False, action='store_true', 
        help='Use all genes for calculations (default: create subgraph from gene list')
    parser.add_argument(
        '-n', '--network', required=True,
        type=str, help='Network gml file')
    parser.add_argument(
        '-t', '--tag', required=True,
        type=str, help='Dataset type (for outfile)')
    parser.add_argument(
        '-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.gene_list, args.all_genes, args.network, args.tag, args.out


def parse_gene_list(gene_list, network):
    with open(gene_list, 'r') as f:
        genes = [line.rstrip('\n') for line in f]

    # count to see how many are included/missing in the network and prune list
    found_genes = []
    counter, not_found = 0, 0
    for gene in genes:
        counter += 1
        if gene not in network:
            not_found += 1
        else:
            found_genes.append(gene)

    print(f'[saltMA] {counter-not_found} genes found of {counter}')

    return found_genes


def get_coefficients(genes, all_genes, network, tag, out):
    """
    Get clustering coefficients for input genes.

    `nx.clustering` returns a single float provided a node
    """

    with open(out, 'w') as f:
        fieldnames = [
            'd', 'gene', 'clustering_coefficient', 'triangles', 'degree']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        if not all_genes:
            # create subgraph with only mutated genes
            subgraph = network.subgraph(genes)
        elif all_genes:
            subgraph = network

        for gene in tqdm(genes):
            clustering_coefficient = nx.clustering(subgraph, nodes=gene)
            triangle_count = nx.triangles(subgraph, gene)
            gene_degree = subgraph.degree(gene)

            out_dict = {
                'd': tag, 'gene': gene,
                'clustering_coefficient': clustering_coefficient,
                'triangles': triangle_count,
                'degree': gene_degree
                }
            writer.writerow(out_dict)


def main():
    gene_list, all_genes, network, tag, out = args()

    print(f'[saltMA] reading in {network}')
    G = nx.read_gml(network)

    print('[saltMA] parsing gene list...')
    found_genes = parse_gene_list(gene_list, G)

    print('[saltMA] computing clustering coefficients')
    get_coefficients(found_genes, all_genes, G, tag, out)

    print('[saltMA] complete.')

if __name__ == '__main__':
    main()

        

