"""
shortest_path_resamples.py - resample nodes and compute shortest paths
"""

import random
import csv
import argparse
import networkx as nx
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='', 
        usage='python shortest_path_resamples.py [options]')

    parser.add_argument('-f', '--fname', required=False,
        type=str, help='List of nodes to resample from [optional, default=full network]')
    parser.add_argument('-t', '--treatment', required=True,
        type=str, help='List of treatment nodes (for resample size)')
    parser.add_argument('-g', '--gml', required=True,
        type=str, help='Network in gml format')
    parser.add_argument('-r', '--replicates', required=True,
        type=int, help='Number of replicates')
    parser.add_argument('--no_replacement', required=False, action='store_true', 
        help='Resample without replacement (default: with) [optional]')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.treatment, args.gml, \
        args.replicates, args.no_replacement, args.out


def parse_gene_lists(fname, treatment, G):
    """Get number of genes to draw

    Parameters
    ----------
    fname : str
        gene list to resample from [optional]
    treatment : str
        gene list used to determine m (# of nodes to resample)
    G : networkx.classes.graph.Graph
        network to sample from

    Returns
    -------
    list(treatment_genes) : list
        list of non-unique genes from treatment fname
    list(set(gene_sample_set)) : list
        list of unique genes to resample from
    """
 
    with open(treatment, 'r') as f:
        treatment_genes = [line.rstrip('\n') for line in f]
    
    # filter treatment genes
    treatment_genes = [
        gene for gene in treatment_genes
        if gene in G]

    # get genes to sample from if provided
    if fname:
        with open(fname, 'r') as f:
            gene_sample_set = [line.rstrip('\n') for line in f]
        
        # only keep genes actually in network
        gene_sample_set = [
            gene for gene in gene_sample_set
            if gene in G]
    elif not fname:
        gene_sample_set = list(G.nodes)

    print(
        f'[saltMA] resampling {len(treatment_genes)} nodes from {len(gene_sample_set)} total genes')

    return list(treatment_genes), list(set(gene_sample_set))

def resample_nodes(G, gene_sample_set, replicate_size, no_replacement, iteration):
    """For a given iteration, resample nodes and get min path per node

    Parameters
    ----------
    G : networkx.classes.graph.Graph
        network to sample from
    gene_sample_set : list
        set of genes to sample from
    replicate_size : int
        number of nodes to draw
    no_replacement : bool
        if True, sample nodes without replacement (default: with replacement)
    iteration : int
        iteration count

    Returns
    -------
    gene1 : str
        gene name
    sample_closest_gene : str
        a given gene whose shortest path w gene 1 == min gene1 shortest path
    min_shortest_path : int
        the min shortest path gene 1 has with any other resampled node
    sample_shortest_path : list
        the path to sample_closest_gene, provided for reference
    """
    # sample nodes
    if no_replacement:
        subsampled_nodes = random.sample(gene_sample_set, replicate_size)
    elif not no_replacement:
        subsampled_nodes = random.choices(gene_sample_set, k=replicate_size)

    # get shortest paths
    for i, gene1 in tqdm(enumerate(subsampled_nodes), desc='finding min shortest paths'):
        if subsampled_nodes.count(gene1) > 1:
            # shortest path will always be itself if node is 'mutated' multiple times
            # will never execute if run in --no_replacement mode
            shortest_path_dict = {gene1: [gene1]}
        elif subsampled_nodes.count(gene1) == 1:
            shortest_path_dict = {
                gene2: nx.shortest_path(G, gene1, gene2)
                for gene2 in subsampled_nodes if gene1 != gene2
                }
        else:
            raise Exception("this probably shouldn't happen")
        shortest_path_values = [len(path) for path in shortest_path_dict.values()]
        min_shortest_path = min(shortest_path_values)

        # only keep min path
        for gene, shortest_path in shortest_path_dict.items():
            if len(shortest_path) == min_shortest_path:
                sample_closest_gene = gene
                sample_shortest_path = shortest_path

        yield gene1, sample_closest_gene, min_shortest_path, sample_shortest_path


def perform_draws(fname, treatment, gml, replicates, no_replacement, out):
    """Draw m nodes from network n times and calculate min paths

    Parameters
    ----------
    fname : str
        gene list to resample from [optional]
    treatment : str
        gene list used to determine m (# of nodes to resample)
    gml : str
        path to network file (gml)
    replicates : int
        number of iterations n
    no_replacement : bool
        whether to sample without replacement (default: with)
    out : str
        file to write to


    Returns
    -------
    None
    """

    print('[saltMA] loading in network')
    G = nx.read_gml(gml)

    print('[saltMA] parsing gene lists')
    if not fname:
        print('[saltMA] no gene list provided - resampling from full network')
    treatment_genes, gene_sample_set = parse_gene_lists(fname, treatment, G)
    replicate_size = len(treatment_genes)
    print(f'[saltMA] will be resampling {replicate_size} genes {replicates} times')

    with open(out, 'w') as f:
        fieldnames = [
            'iteration', 'gene1', 'gene2', 'min_path', 'sample_path']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        print('[saltMA] resampling...')
        for iteration in tqdm(range(replicates)):
            i = iteration + 1 # start counting at 1

            iterator = resample_nodes(
                G, gene_sample_set, replicate_size, no_replacement, i)
            for gene1, gene2, min_path, sample_path in iterator:
                out_dict = {
                    'iteration': i,
                    'gene1': gene1,
                    'gene2': gene2,
                    'min_path': min_path,
                    'sample_path': sample_path
                    }
                writer.writerow(out_dict)

                
def main():
    fname, treatment, gml, replicates, no_replacement, out = args()
    perform_draws(*args())

if __name__ == '__main__':
    main()

        

