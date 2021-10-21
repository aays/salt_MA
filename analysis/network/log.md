
## 23/9/2021

getting started with the network analysis after getting Sara's
script and installing networkx

```
$ python analysis/network/network_draft_code.py --help

usage: network_draft_code.py [-h] [-c] [-f FILE]
                             Network nodes_file num_iterations

positional arguments:
  Network          network file in gml format
  nodes_file       list of nodes in text file
  num_iterations   number of iterations

optional arguments:
  -h, --help       show this help message and exit
  -c, --connected  check if network is fully connected
  -f FILE          random nodes taken from this file rather than the entire
                   network
```

need to prep a nodes file - e.g. a file containing just
gene names on separate lines - seems the gene names
to use are the first ones in the `feature_names` column
output by mut describer

```python
import csv

with open('data/network/saltMA_gene_list.txt', 'w') as f_out:
    with open('data/mutations/mut_describer/muts_described.final.tsv', 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=['gene'])
        for record in tqdm(reader):
            feature_names = eval(record['feature_names'])
            if feature_names and not feature_names[0] == 'n/a':
                line = {'gene': feature_names[0]}
                writer.writerow(line)

# and then repeated in 'a' write mode with indels_described
```

naive run:

```bash
python analysis/network/network_draft_code.py \
data/network/chlamyNET.gml data/network/saltMA_gene_list.txt 10
```

doesn't seem like the output is written to file:

```
reading in gml file
sys.argv is ['analysis/network/network_draft_code.py', 'data/network/chlamyNET.gml', 'data/network/saltMA_gene_list.txt', '10']
first node in list Cre05.g242650
length of all mutated nodes 326
name of network file data/network/chlamyNET.gml
nodes found in the network  162
nodes NOT found in network 164
shortest path for mutated nodes  111194
i is 0
shortest path for random nodes 107698
i is 1
shortest path for random nodes 112057
i is 2
shortest path for random nodes 111376
i is 3
shortest path for random nodes 107579
i is 4
shortest path for random nodes 112650
i is 5
shortest path for random nodes 108968
i is 6
shortest path for random nodes 112331
i is 7
shortest path for random nodes 111940
i is 8
shortest path for random nodes 107139
i is 9
shortest path for random nodes 109532
shorter found is  5 out of total  10
```

need to update this so that it writes to file instead
and then repeat with more iterations (0 alone, 5 alone, all combined)

## 24/9/2021

giving this a go:

```bash
time python analysis/network/network_draft_code.py \
data/network/chlamyNET.gml data/network/saltMA_gene_list.txt 10 \
-o test_network.tsv
```

looks good - queuing this up for the overall gene list (will need to make 0 and 5 files afterwards)

```bash
time python analysis/network/network_draft_code.py \
data/network/chlamyNET.gml data/network/saltMA_gene_list.txt 100 \
-o data/network/salt_all_distance.tsv
```

## 25/9/2021

looks good, and done in two hours - now to repeat for 0 and 5 separately

creating the gene lists:

```python
import csv
from tqdm import tqdm

# 0 muts
with open('data/network/saltMA_0_gene_list.txt', 'w') as f_out:
    with open('data/mutations/mut_describer/muts_described.final.tsv', 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=['gene'])
        for record in tqdm(reader):
            if record['mutant_sample'].endswith('0'):
                feature_names = eval(record['feature_names'])
                if feature_names and not feature_names[0] == 'n/a':
                    line = {'gene': feature_names[0]}
                    writer.writerow(line)

# 0 indels
with open('data/network/saltMA_0_gene_list.txt', 'a') as f_out:
    with open('data/mutations/mut_describer/indels_described.tsv', 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=['gene'])
        for record in tqdm(reader):
            if record['mutant_sample'].endswith('0'):
                feature_names = eval(record['feature_names'])
                if feature_names and not feature_names[0] == 'n/a':
                    line = {'gene': feature_names[0]}
                    writer.writerow(line)

# and then repeated for 5
```

and now for the network stuff:

```bash
time for i in 0 5; do
    python analysis/network/network_draft_code.py \
    data/network/chlamyNET.gml data/network/saltMA_${i}_gene_list.txt 100 \
    -o data/network/salt_${i}_all_distance.tsv;
done
```

## 16/10/2021

today - rejigging the network script to return the all x all gene distances
it uses to compute shortest paths

testing:

```bash
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/saltMA_0_gene_list.txt 1 \
--distances -o test_gene_paths.tsv
```

finally got this working - also added a duplicate removal line, although
I need to check whether that should be kept for that summed analysis

next steps

1. clean up/organize the `data/network` folder
2. create gene lists for MA and HS
3. rerun this for MA and HS
4. look at distributions of shortest paths
5. redo bootstrapping as well

re: the gene lists - I might have to do some cross checking of gene names
against `chlamyNET` - genes sometimes have multiple names and I don't know
yet whether the first gene listed in mut describer output is the gene name
in the network as well

also - reading Levin and Levina 2019, it is worth looking into
whether these bootstrapped networks have similar distributions to 
the original network - need to read more on that, but there's an example
of bootstrapping for shortest path length in there although it seems
on first glance the problem of having a 'good' means of comparing
distributions for this is inconclusive? not sure

update: the more I read on this, the more I'm convinced this method is straight up wrong

concerns include -

- are genes mutated multiple times still only represented once in the network?
- resampling with replacement is flawed since drawing the same node twice will create a 0-distance edge
    - this will reduce the sum of shortest paths, and in a way reduce the size of the network
- since the mutated genes are a subset of the full network, why not work with permutation tests instead?
    - if the mutated gene set is of size L, draw L random nodes for n iterations and get their shortest path
    - could generate CIs this way
- bootstrapping networks itself is more mathematically complicated than resampling nodes, and not just for
  the second reason above - the resampled networks must have similarly distributed adjacency matrices
    - ensuring this is difficult and there's no existing clean implementation for it that I can find! 
- based on Levin and Levina 2019, our network size (<100) is likely to yield massive CIs

## 19/10/2021

going to revisit the above soon - for now, need to generate all x all shortest paths for MA and HS

```bash
# make matrices for salt combined
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/saltMA_gene_list.txt 1 \
--distances -o data/network/matrices/saltMA.tsv

```

generate gene sets:

```python
import csv
from tqdm import tqdm

# MA orig
with open('data/network/gene_lists/MA_orig_gene_list.txt', 'w') as f_out:
    with open('data/rate/mut_describer/ma.gene_sets.tsv', 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=['gene'])
        for record in tqdm(reader):
            feature_names = eval(record['feature_names'])
            if feature_names and not feature_names[0] == 'n/a':
                line = {'gene': feature_names[0]}
                writer.writerow(line)

# adaptation
with open('data/network/gene_lists/adaptation_gene_list.txt', 'w') as f_out:
    with open('data/rate/mut_describer/adaptation.gene_sets.tsv', 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=['gene'])
        for record in tqdm(reader):
            feature_names = str(record['Gene.primaryIdentifier']) # different handling!
            if feature_names and not feature_names == 'n/a':
                line = {'gene': feature_names}
                writer.writerow(line)
```

make matrices:

```bash
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/adaptation_gene_list.txt 1 \
--distances -o data/network/matrices/adaptation.tsv

# this is an upper triangle all x all with 2025 nodes - that's 4.1 million comparisons!!! 
# going to run and leave anyways, but I'm not certain we'll need all these...
# ended up taking 1.5 hours to complete
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/MA_orig_gene_list.txt 1 \
--distances -o data/network/matrices/MA_orig.tsv

```

## 21/10/2021

going to make a new version of the `MA_orig` file without the full paths listed to save
space (it's currently a 340 MB file!) 

```bash
cut -f1-3 data/network/matrices/MA_orig.tsv > data/network/matrices/MA_orig_light.tsv
```

down to 52 MB!

might as well also create matrices for salt 0 and 5 just in case:

```bash
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/saltMA_0_gene_list.txt 1 \
--distances -o data/network/matrices/saltMA_0.tsv

time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/saltMA_5_gene_list.txt 1 \
--distances -o data/network/matrices/saltMA_5.tsv
```
