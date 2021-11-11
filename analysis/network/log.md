
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

## 27/10/2021

today - starting on the betweenness centrality + degree analyses

gene lists are in `data/network/gene_lists` - network stats to be symlinked into
`data/network` (19530 line file, 5.5 Mb - so could work with it offline if needed)

seems the best way forward is to use the gene lists to just output filtered
versions of that file - or better yet, just get the gene lists locally
and do the merging etc in R directly

## 28/10/2021

so post meeting, I'm realizing I need to make sure the shortest paths
are always calculated using the full network - that way, I could take
the min shortest path for any HS gene (for example) and work with that

also need to look into clustering coefficients, how they relate to shortest
paths, and if those could be used for within-network cluster discovery vs
comparing between networks

Rob also has a bootstrap significance test I could use to compare distributions
of shortest paths, but first I need to make sure they're calculated correctly
and then those all x all values need to be transformed so we're just working with
the min path values

## 3/11/2021

update - after a multi day dive into the literature - going to be using the clustering
coefficient instead of shortest paths, since I think it gets at the question much
more directly

need to create a script that'll take in a list of genes and return their clustering
coefficients to a tsv file

trying this out for saltMA - 

```bash
time python analysis/network/clustering_coefficients.py \
--gene_list data/network/gene_lists/saltMA_gene_list.txt \
--network data/network/chlamyNET.gml \
--tag saltMA --out cluster_test.tsv
```

looks good - but just about half the genes were actually found in chlamyNET:

```
[saltMA] reading in data/network/chlamyNET.gml
[saltMA] parsing gene list...
[saltMA] 162 genes found of 326
[saltMA] computing clustering coefficients
100%|█████████████████████████████████████████████████████████████████████████████| 162/162 [00:00<00:00, 2310.50it/s]
[saltMA] complete.
```

and I should probably keep track of these numbers moving forwards

doing HS + MA + 0 and 5 separately:

```bash
for d in adaptation MA_orig saltMA_0 saltMA_5; do
    time python analysis/network/clustering_coefficients.py \
    --gene_list data/network/gene_lists/${d}_gene_list.txt \
    --network data/network/chlamyNET.gml \
    --tag ${d} --out data/network/clustering/${d}_clustering.tsv;
done
```

done in half a minute, but:

```
[saltMA] reading in data/network/chlamyNET.gml
[saltMA] parsing gene list...
[saltMA] 127 genes found of 251
[saltMA] computing clustering coefficients
100%|█████████████████████████████████████████████████████████████████████████████| 127/127 [00:00<00:00, 1638.68it/s]
[saltMA] complete.

real    0m27.393s
user    0m32.160s
sys     0m7.210s
[saltMA] reading in data/network/chlamyNET.gml
[saltMA] parsing gene list...
[saltMA] 2737 genes found of 6172
[saltMA] computing clustering coefficients
100%|███████████████████████████████████████████████████████████████████████████| 2737/2737 [00:01<00:00, 2087.03it/s]
[saltMA] complete.

real    0m28.662s
user    0m33.779s
sys     0m6.905s
[saltMA] reading in data/network/chlamyNET.gml
[saltMA] parsing gene list...
[saltMA] 88 genes found of 162
[saltMA] computing clustering coefficients
100%|███████████████████████████████████████████████████████████████████████████████| 88/88 [00:00<00:00, 1092.85it/s]
[saltMA] complete.

real    0m27.861s
user    0m33.402s
sys     0m6.412s
[saltMA] reading in data/network/chlamyNET.gml
[saltMA] parsing gene list...
[saltMA] 74 genes found of 164
[saltMA] computing clustering coefficients
100%|███████████████████████████████████████████████████████████████████████████████| 74/74 [00:00<00:00, 1480.45it/s]
[saltMA] complete.
```

looks like a good number of genes are missing here - it may be worth revisiting some missing genes
and seeing if they're listed in chlamyNET under other names

a manual search of the first few not-found-genes in the adaptation and saltMA sets yields nothing -
using the PACid doesn't change the fact that they're not in the network by the looks of it

wait - this should be done with subgraphs, since for a given mutated node we aren't concerned
with whether nearby unmutated genes are clustered

## 4/11/2021

today - the clustering coefficient seems to be a little too limited to proximate nodes,
and while I could look into an 'extended' coefficient I'm going to try and detect clusters
directly instead of computing per-node metrics 

seems the Louvain algorithm is the method of choice here, after which I could implement
the qs test (Kojaku + Masuda 2018) which compares communities to randomly generated networks
but controls for differences in community sizes

installing python packages `python-louvain`/`community` and `qstest` - looks
like the latter is on pypi only

```bash
conda install -c conda-forge python-louvain
pip install qstest # this fails - need to get from repo and run python setup.py install
```

permutation test - create random draws from MA (without replacement)
at the same size as the HS - do this ~100 times - compare with the distribution
of min shortest paths from HS 

could stack histograms of MA draws - the distribution should be less
0 focused this way

could do bin specific too - 'how many times is the 0-1 bin in the MA draws
smaller than the HS dataset'? 

## 5/11/2021

today - actually implementing the permutation test

this needs to take in both the adaptation and MA gene lists, and generate
M subsamples of size n (where n = size of the adaptation data set) 
before outputting all x all min shortest paths

it does occur to me that we specifically care about min shortest paths - it
might save space to _just_ report those, since otherwise the outfile will be
millions of lines long

## 7/11/2021

done now - giving this a go:

```bash
python analysis/network/shortest_path_resamples.py \
--fname data/network/gene_lists/MA_orig_gene_list.txt \
--treatment data/network/gene_lists/adaptation_gene_list.txt \
--gml data/network/chlamyNET.gml \
--replicates 10 \
--out resample_test.tsv
```

finally debugged and looking good - took 7 min for 10 iterations, which is
on the longer side but not too bad

here goes:

```bash
mkdir -p data/network/resamples

python analysis/network/shortest_path_resamples.py \
--fname data/network/gene_lists/MA_orig_gene_list.txt \
--treatment data/network/gene_lists/adaptation_gene_list.txt \
--gml data/network/chlamyNET.gml \
--replicates 100 \
--out data/network/resamples/ma_resample_100.tsv
```

although I just realized something - the original all x all shortest paths
seemed to all have been generated using a subgraph! which isn't quite right - and so I need
to rerun that as well

```bash
# updated to not create a subgraph and just get shortest paths from the full graph
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/adaptation_gene_list.txt 1 \
--distances -o data/network/matrices/adaptation_full.tsv

# this is an upper triangle all x all with 2025 nodes
# will filter once again later
time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/MA_orig_gene_list.txt 1 \
--distances -o data/network/matrices/MA_orig_full.tsv

cut -f1-3 data/network/matrices/MA_orig_full.tsv > data/network/matrices/MA_orig_full_light.tsv

time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/saltMA_0_gene_list.txt 1 \
--distances -o data/network/matrices/saltMA_0_full.tsv

time python analysis/network/shortest_path.py \
data/network/chlamyNET.gml data/network/gene_lists/saltMA_5_gene_list.txt 1 \
--distances -o data/network/matrices/saltMA_5_full.tsv
```

strange - these seem to be the exact same outfiles as before, and I've double
checked that the subgraph bit has been disabled. did some testing and it seems
shortest paths are supposed to raise errors if the subgraph doesn't contain
the requisite in between nodes -

```python
>>> for gene1 in tqdm(genes):
...     for gene2 in genes:
...         if gene1 != gene2:
...             s = nx.shortest_path(G, gene1, gene2)
...             for gene in s:
...                 if gene not in genes:
...                     counter += 1
100%|████████████████████████████████████████████████████████████████████████████████| 127/127 [00:52<00:00,  2.65it/s]
>>> counter
92790

>>> H = G.subgraph(genes)
>>> counter2 = 0
... for gene1 in tqdm(genes):
...     for gene2 in genes:
...         if gene1 != gene2:
...             s = nx.shortest_path(H, gene1, gene2)
...             for gene in s:
...                 if gene not in genes:
...                     counter2 += 1
  0%|                                                                                          | 0/127 [00:00<?, ?it/s]
Traceback (most recent call last):
  File "<stdin>", line 5, in <module>
  File "/home/hasans11/.conda/env/work/lib/python3.6/site-packages/networkx/algorithms/shortest_paths/generic.py", line 160, in shortest_path
    paths = nx.bidirectional_shortest_path(G, source, target)
  File "/home/hasans11/.conda/env/work/lib/python3.6/site-packages/networkx/algorithms/shortest_paths/unweighted.py", line 224, in bidirectional_shortest_path
    results = _bidirectional_pred_succ(G, source, target)
  File "/home/hasans11/.conda/env/work/lib/python3.6/site-packages/networkx/algorithms/shortest_paths/unweighted.py", line 292, in _bidirectional_pred_succ
    raise nx.NetworkXNoPath(f"No path between {source} and {target}.")
networkx.exception.NetworkXNoPath: No path between Cre01.g012700 and Cre01.g012900.

No path between Cre01.g012700 and Cre01.g012900.
>>> 'Cre01.g012900' in genes
True

```

so I'm not sure why these are exactly the same, or for that matter why errors weren't raised in
the earlier code when the subgraph was being made

WAIT - I misunderstood Sara's code - the subgraph is only made if a list of nodes is provided,
so ultimately this was still an all x all run originally - meaning I can remove the 'full' versions
of these files

alright - now to export the resampled distributions (took 70 min as expected)
into RStudio and look at this distribution

## 8/11/2021

update: need to also do resamples from the _entire_ network - modifying the script accordingly:

```bash
# running without providing --fname arg
python analysis/network/shortest_path_resamples.py \
--treatment data/network/gene_lists/adaptation_gene_list.txt \
--gml data/network/chlamyNET.gml \
--replicates 100 \
--out data/network/resamples/chlamynet_full_resample_100.tsv
```

## 9/11/2021

today - downloading the metabolic network and repeating some of these analyses with it

```
# Imam et al 2015, the plant journal - 'A refined reconstruction of Chlamydomonas metabolism...'
wget https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Ftpj.13059&file=tpj13059-sup-0017-DataS2.tar
```

is it meaningful to use degree and BC here? degree basically tells us how many
metabolites were involved in a reaction, while BC makes little sense to me here
given that species point to reactions (e.g. there are no longer paths past that)

```R
> x = d %>% filter(str_detect(transcript, 'Cre17.g698000'))
> x
# A tibble: 1 x 30
  transcript SFS0  gene  diffs0    k4 Cincerta_transc…     pi4 SFS4
  <chr>      <chr> <chr>  <dbl> <dbl> <chr>              <dbl> <chr>
1 Cre17.g69… [111… Cre1…     13 0.145 g11475.t1        0.00400 [315…
# … with 22 more variables: aligned_sites <dbl>, sites4 <dbl>, sites0 <dbl>,
#   pi0 <dbl>, ness_ID <dbl>, k0 <dbl>, PID <dbl>, diffs4 <dbl>,
#   DegreechlamynetUnDirected <dbl>, OutdegreechlamynetDirected <dbl>,
#   IndegreechlamynetUnDirected <dbl>,
#   BetweennessCentralitychlamynetDirected <dbl>,
#   BetweennessCentralitychlamynetUnDirected <dbl>,
#   OutdegreechlamynetUnDirected <dbl>, IndegreechlamynetDirected <dbl>,
#   maxDegree <dbl>, maxBC <dbl>, meanBC <dbl>, sum_product_outdegrees <dbl>,
#   mean_product_outdegrees <dbl>, meanDegree <dbl>, transcript55 <chr>
> x$maxDegree
[1] 7
> x$meanDegree
[1] 6.5
# one rxn has 6 metabolites, the other has 7... is that meaningful? 
```

one possible takeaway could be to get genes present in all three datasets (categorized as
'essential' in the paper, since they're used with and without light + with and without
external nutrients) and see whether there's more mutations in these than expected? I suppose
that overlaps with the Perrineau/expressed analysis

## 10/11/2021

will revisit this earlier stuff about the metabolic network later, but for now - need
to update the clustering coefficient script - have it output degree for a given gene
as well as the number of triangles, so that we can investigate that spike in the adaptation
dataset better

```bash
for d in adaptation MA_orig saltMA_0 saltMA_5 saltMA; do
    time python analysis/network/clustering_coefficients.py \
    --gene_list data/network/gene_lists/${d}_gene_list.txt --all_genes \
    --network data/network/chlamyNET.gml --tag ${d} \
    --out data/network/clustering/${d}_clustering_all.tsv;
done
```
