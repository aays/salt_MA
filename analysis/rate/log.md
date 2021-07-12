
## 11/6/2021

rate analysis! 

first step - need to get the callable sites for each sample

I think the easiest way to do this is to create a giant annotation table-like
file with all of the samples that can then be tabixed and references

first up - this is going to take a hot minute, but I need to bgzip and tabix
each of the gvcfs, since I'm likely going to write to the final callables table
one chromosome at a time

```bash
# in data/alignments/genotyping/gvcfs
time for fname in *.vcf; do
    echo "currently on ${fname}"
    time bgzip ${fname}
    tabix -p vcf ${fname}.gz
done
```

while this runs, a first pass at a script to do this could be something like:

```
sample_names = [...]
fieldnames = ['chrom', 'pos']
fieldnames.extend(sample_names)
writer = csv.DictWriter(..., fieldnames=fieldnames)

def get_calls(chrom, pos):
    # iterate through samples and get 'called' status as 1 or 0
    # save as dict vals
    # return dict

for chrom in chromosomes: # some sort of length lookup
    for pos in chromosome:
        d = get_calls(chrom, pos)
        writer.writerow(d)
```

## 14/6/2021

finally getting on this! 

going to start on a script called `callable_sites.py` with the above structure

alright, time for a first pass:

```bash
time python analysis/rate/callable_sites.py \
--gvcf_dir data/alignments/genotyping/gvcfs \
--out data/rate/all_callable.tsv
```

so this is VERY slow - about to take 4 hours just for chr1! 

I think I need to create lookups for each sample - here's a fun little algorithm
to efficiently do that:

```python
l = ''
current_pos = 1
for record in tqdm(v('chromosome_1:1-8033585')):
    if record.POS == current_pos:
        l += '1'
        current_pos += 1
    elif record.POS > current_pos:
        to_add = record.POS - current_pos
        str_to_add = '0' * to_add 
        l += str_to_add
        l += '1'
        current_pos = record.POS + 1
```

this runs in 10 seconds for all of chromosome 1! I do worry about holding all this in
memory though - perhaps it can be done in 'chunks' of 100k sites

## 15/6/2021

trying this with a quick lookup implementation:

```bash
time python analysis/rate/callable_sites.py \
--gvcf_dir data/alignments/genotyping/gvcfs/ \
--out data/rate/all_callable.tsv
```

## 16/6/2021

looks good - wrapped up the entire genome in just two hours! 

```bash
bgzip data/rate/all_callable.tsv
tabix -p vcf data/rate/all_callable.tsv.gz
```

second command is breaking cause the header needs to be commmented - doing
a bit of 'surgery' here

```bash
# in data/rate
bgzip -d all_callable.tsv.gz
head -n 1 all_callable.tsv > header.temp
vim header.temp # add a hash at the start to indicate header
bgzip all_callable.tsv
tabix -p vcf all_callable.tsv.gz
```

but also - somehow I didn't notice this earlier - all the positions have been dropped
by 1? so the first position in `chromosome_1` is 0 and the last is 8033584

I don't want to rerun the `callable_sites.py` script after editing it again - could
probably do this more 'manually' after unzipping the file

```python
import csv
from tqdm import tqdm
with open('all_callable.tsv', 'r') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    with open('all_callable_new.tsv', 'w') as f_out:
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for record in tqdm(reader):
            record['pos'] = str(int(record['pos']) + 1)
            writer.writerow(record)
```

## 18/6/2021

finally wrapping this up - file looks good and was done in 46 min - now to
bgzip and tabix

renaming the new file to `all_callable.tsv` first though

tabix works as intended! 

```bash
tabix all_callable.tsv.gz chromosome_2:1-2
```

alright, next steps: rate calculations! from Rob's notes:

- estimate rate overall 
- estimate rate of SaltMA, MA, SaltAdaptation

we already have the MA rates from Ness 2015 - I might have to do salt adaptation
line calculations myself based on the previous mutation tables and 'callables' reports,
while the saltMA rate should be doable once I have the generation time data from Rob

going to put this on hold since the generation time data aren't ready yet - back
to working on spectrum/context

## 24/6/2021

finally getting back on this after dealing with SL26 - managed to squeeze one good SNM out of that

going to write a script that takes in a mut describer file and a callables lookup
file to report per site rates - generations should also be a parameter there but I'll set
that to 1 as default - though since generation time might be sample specific I might have
to update that script to take in a lookup file

the output format should look something like

```
sample chrom mutations callable_sites generations mut_rate
```

where 
- `mutations` = absolute number of muts 
- `callable_sites` = absolute number of callable sites
- `generations` = number of generations muts were accumulated for
- `mut_rate` = mutations / (callable sites * generations)

the original MA files where Rob did something similar are instead structured like

```
window sample_1_muts sample_1_sitegens sample_1_mutrate sample_2_muts sample_2_sitegens ...
```

where `window` is a single column formatted like `chromosome_1:1-200000` - but given the saltMA
mutation counts are far lower I think I can get away with entire-chromosome values, after which
processing the data in R will let me convert the file to wide format if I need

let's get this going - will call this script `calculate_rate.py` - though this will likely
only work for the saltMA files on first pass

maybe I'll need to create a new callables table off of the adaptation VCFs?? a problem
for future me to figure out...
 
wait - looks like the callable table wasn't made correctly - I'm missing every x * 1e5 value
(eg 200000, 300000)

found the bug I think - line 143 didn't need to `- 1` - rerunning and taking a break
since this will take 2 hours! 

```bash
time python analysis/rate/callable_sites.py \
--gvcf_dir data/alignments/genotyping/gvcfs/ \
--out data/rate/all_callable_fix.tsv
```

## 25/6/2021

took two hours again yesterday, but file looks good! although I forgot to correct the header again...

```bash
head -n 1 all_callable_fix.tsv > header.txt
vim header.txt # add hash before chrom
tail -n +2 all_callable_fix.tsv > all_callable_partial.tsv # took 30 sec
cat header.txt all_callable_partial.tsv > all_callable.tsv
```

back to getting rate calculations - the columns should again look like

```
sample chrom mutations callable_sites generations mut_rate
```

looks good, I think? no errors in the test run, so here goes -

```bash
time python analysis/rate/calculate_rate.py \
--mut_table data/rate/mut_describer/muts_described.final.tsv \
--callables_table data/rate/all_callable.tsv.gz \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--out rate_test.tsv
```

of course I forgot about SL29...damn it 

## 28/6/2021

alright, SL29 VCFs made - now to remake the callables file

```bash
bgzip data/alignments/genotyping/gvcfs/SL29_0.g.vcf
tabix data/alignments/genotyping/gvcfs/SL29_0.g.vcf.gz
bgzip data/alignments/genotyping/gvcfs/SL29_5.g.vcf
tabix data/alignments/genotyping/gvcfs/SL29_5.g.vcf.gz

time python analysis/rate/callable_sites.py \
--gvcf_dir data/alignments/genotyping/gvcfs/ \
--out data/rate/all_callable_fix.tsv
```

and now this can finally be run again:

```bash
time python analysis/rate/calculate_rate.py \
--mut_table data/rate/mut_describer/muts_described.final.tsv \
--callables_table data/rate/all_callable.tsv.gz \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--out rate_test.tsv
```

## 29/6/2021

this will need to be rerun on indels once I've determined which to keep

need to get adaptation rates next - but in the meantime, Rob
mentions that I should set generations to 250

I *could* technically do this in R working off of `rate_test.tsv` but let's
just rerun the script for simplicity after adding the ability to give
a single constant for all generations

```bash
time python analysis/rate/calculate_rate.py \
--mut_table data/rate/mut_describer/muts_described.final.tsv \
--callables_table data/rate/all_callable.tsv.gz \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--generation_count 250 \
--out data/rate/saltMA_SNM_rate_250.tsv
```

hold on - I haven't been filtering the gvcfs for MQ > 24 and DP > 10 when making
the callables table... time to remake this a third time?

quick check:

```python
>>> from tqdm import tqdm
>>> from cyvcf2 import VCF
>>> reader = VCF('data/alignments/genotyping/gvcfs/CC1373_0.g.vcf.gz')
>>> count, failed = 0, 0
>>> for record in tqdm(reader):
...     count += 1
...     if record.INFO.get('MQ'):
...         if record.INFO.get('MQ') < 24:
...             failed += 1
...     elif record.INFO.get('DP'):
...         if record.INFO.get('DP') < 10:
...             failed += 1
26254928it [02:10, 201821.35it/s]
>>> failed, count
(3553, 27210555)
```

around 0.01% of sites, but I should still redo this... sigh

```bash
mkdir -p data/alignments/genotyping/filtered_gvcfs
```

and here's a quick and dirty script to filter the gvcfs:

```python
import sys
import os
import glob
from cyvcf2 import VCF
from cyvcf2 import Writer
from tqdm import tqdm

def filter_gvcf(vcf_name):
    reader = VCF(vcf_name):
    writer = Writer(vcf_name.replace('/gvcfs', '/filtered_gvcfs').rstrip('.gz'), reader)
    writer.write_header()
    count, failed = 0, 0
    for record in tqdm(reader, desc=os.path.basename(vcf_name)):
        count += 1
        mq, dp = record.INFO.get('MQ'), record.INFO.get('DP')
        passing = True
        elif mq:
            if mq < 24:
                passing = False
        elif dp:
            if dp < 10:
                passing = False
        elif record.num_het != 0.0:
            passing = False
        if passing:
            writer.write_record(record)
        elif not passing:
            failed += 1
    print(vcf_fname, count, failed, failed / count)

fnames = glob.glob('data/alignments/genotyping/gvcfs/*.gz')

for fname in tqdm(fnames):
    filter_gvcf(fname)
```

callable table v3:

```bash
time python analysis/rate/callable_sites.py \
--gvcf_dir data/alignments/genotyping/gvcfs/ \
--out data/rate/all_callable_fix.tsv
```

## 30/6/2021

once more, after replacing the old callables file:

```bash
time python analysis/rate/calculate_rate.py \
--mut_table data/rate/mut_describer/muts_described.final.tsv \
--callables_table data/rate/all_callable.tsv.gz \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--generation_count 250 \
--out data/rate/saltMA_SNM_rate_250.tsv
```

## 1/7/2021

today - Ka/Ks calculations for saltMA

need to get Ka/Ks genome wide, salt genes, and expressed genes

before that even - I'll need to actually get lists of these salt genes and expressed genes
to use in whatever script I end up writing

the salt genes look to be in `salt_lines/annotations/expression/` - `perrineau.geneIDs.txt`
is the list of genes while `fpkm_per_transcript.gz` contains fpkm values (where anything >= 1 
is considered an 'expressed gene')

symlinking these files into `data/rate/gene_lists` and creating a list of just expressed genes:

```python
import gzip
from tqdm import tqdm

with open('expressed_genes.tsv', 'w') as f_out:
    f_out.write('gene\tfpkm\n')
    with open('fpkm_per_transcript.gz', 'rb') as f:
        reader = gzip.GzipFile(fileobj=f)
        for line in tqdm(reader):
            gene, fpkm = line.decode('utf-8').rstrip('\n').split('\t')
            if float(fpkm) >= 1:
                f_out.write(line.decode('utf-8'))
```

need to get, for each strain:

- number of 3D sites
- number of 2D sites
- number of 4D sites
- number of NS callable sites (n0D + 1/3 n3D + 2/3 n2D)
- number of S callable sites (n4D + 2/3 n3D + 1/3 n2D)

the first three stats are technically 'constant' across strains - the only
difference maker is callability of each site - going to need to iterate through
both the annotation table and the callables table to get these values for
the saltMA lines

script should probably loop over the annotation table in chunks
and keep counts - something like this:

```
d = # nested dict containing all samples as keys
# nested dicts should have degeneracy (0D, 2D, etc) as keys
for record in ant.Reader(chrom, window_start, window_end):
    if record # check degeneracy
        # increment counts in nested dict structure above if callable site exists
        # (sample wise, of course)
```

trying it out:

```bash
time python analysis/rate/callable_sites_degeneracy.py \
--callables_table data/rate/all_callable.tsv.gz \
--annotation_table data/references/annotation_table.txt.gz \
--outname degen_callables.tsv
```

this is RIDICULOUSLY slow - crawls down to 300 it per second - redoing with the classic lookup strat
(this is the most I've ever refactored a simple script, good lord)

```bash
time python analysis/rate/callable_sites_degeneracy.py \
--callables_table data/rate/all_callable.tsv.gz \
--annotation_table data/references/annotation_table.txt.gz \
--outname data/rate/degen_callables_lookup.tsv
```

this spends about ~30 sec max making a lookup and then ROCKETS through
each chr in seconds flat - holy cow

did the whole genome in 10 minutes flat! 

## 4/7/2021

today - actually calculating Ka/Ks, and getting the bootstrap
values too if we have time

for each chromosome on each sample - need to get the number of 
NS muts and S muts, and then combine that with `degen_callables_lookup.tsv`
from earlier to get Ka/Ks 

after getting this genomewide, need to do this for candidate genes and expressed
genes - that'll have to be a check using the mut describer output against
the previous list of Perrineau + expressed genes in `data/rate/gene_lists/`

it occurs to me that the calculation for the genes likely
requires a different set of callable site counts (e.g. from just the genes)

it also occurs to me that I have to redo the triplets to separate 0 mutations
from 5 mutations... but that's for a different log at a different time

for now - first order of business - getting the number of S and NS 
mutations per chromosome - going to create a script for this that
can optionally take in a list of genes to work with

calling it `syn_mut_count.py` for lack of a better name

giving this a go:

```bash
mkdir -p data/rate/ka_ks
mv -v data/rate/degen_callables_lookup.tsv data/rate/ka_ks

time python analysis/rate/syn_mut_count.py \
--mut_table data/mutations/mut_describer/muts_described.final.tsv \
--out data/rate/ka_ks/syn_nonsyn_counts_genomewide.tsv
```

looks good - now for the Perrineau genes:

```bash
time python analysis/rate/syn_mut_count.py \
--mut_table data/mutations/mut_describer/muts_described.final.tsv \
--gene_file data/rate/gene_lists/perrineau_genes_only.txt \
--out data/rate/ka_ks/syn_nonsyn_perrineau_counts.tsv
```

also good - though there's only two genes that apply! 

finally, for the expressed genes - 

```bash
time python analysis/rate/syn_mut_count.py \
--mut_table data/mutations/mut_describer/muts_described.final.tsv \
--gene_file data/rate/gene_lists/expressed_genes.tsv \
--out data/rate/ka_ks/syn_nonsyn_expressed_counts.tsv
```

looks good! next up - getting callable sites for these genes, before estimating
Ka/Ks

## 5/7/2021

today - getting callable site counts for the two gene sets (Perrineau and expressed) 

the genomewide callable sites table (`all_callables.tsv.gz`) should suffice for genomewide Ka/Ks

going to make a modified version of `callable_sites_degeneracy` that takes in a 
gene list and returns counts of callable S and NS sites for each

the main limiting factor here is how to efficiently look up the positions of these
genes - I think I need a file that sorts them in order by chromosome too so that
I could leverage `create_degen_lookup` from the earlier script

```bash
# after creating symlink to final.strict.GFF3 in data/references
grep 'mRNA' data/references/final.strict.GFF3 > data/references/mRNA.GFF3
grep 'gene' data/references/final.strict.GFF3 | grep -v 'mRNA' > data/references/genes.GFF3
```

possibly the greediest way to work with this is to store the list of genes of interest in
memory and loop through the GFF, only keeping lines where a 'hit' is found

I could likely do this in an interpreter:

```python
import csv
from tqdm import tqdm

# store gene list lookup in memory
with open('data/rate/gene_lists/perrineau_genes_only.txt', 'r', newline='') as f:
    gene_list = [l[0] for l in csv.reader(f, delimiter='\t')]

with open('data/rate/gene_lists/perrineau_genes.gff3', 'w') as f_out:
    writer = csv.writer(f_out, delimiter='\t')
    with open('data/references/mRNA.GFF3', 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in tqdm(reader):
            start, end = int(line[3]), int(line[4])
            name_dict = {l.split('=')[0]: l.split('=')[1] for l in line[-1].split(';')} # nice
            if any([gene_name in gene_list for gene_name in name_dict.values()]):
                writer.writerow(line)
```
            
and then again for the expressed genes:

```python
with open('data/rate/gene_lists/expressed_genes.tsv', 'r', newline='') as f:
    gene_list = [l[0] for l in csv.reader(f, delimiter='\t')]

with open('data/rate/gene_lists/expressed_genes.gff3', 'w') as f_out:
    writer = csv.writer(f_out, delimiter='\t')
    with open('data/references/mRNA.GFF3', 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in tqdm(reader):
            start, end = int(line[3]), int(line[4])
            name_dict = {l.split('=')[0]: l.split('=')[1] for l in line[-1].split(';')} # nice
            if any([gene_name in gene_list for gene_name in name_dict.values()]):
                writer.writerow(line)

```

looks good - now to get callable sites for each of these

though it just occurs to me that this won't work for the perrineau genes, since
only two actually have mutations and in both cases these are NS muts - meaning Ks would be 0! 

calculating callable sites from these GFFs nonetheless - time for another quick
script that'll take in a GFF and output the number of S and NS callable sites

the output should look like:

```
gene_name start end total_syn total_nonsyn \
sample_1_fold0 sample_1_fold2 sample_1_syn sample_1_nonsyn # etc
```

will call this `callable_genes_degeneracy.py` - there really needs to be a better way to name these...

## 7/7/2021

script is ready to go - mostly riffs off `callable_sites_degeneracy.py`

though the other day I forgot I should just have `sample` as a separate column instead
of having sample-specific columns... hooray for long format! 

```bash
time python analysis/rate/callable_genes_degeneracy.py \
--callables_table data/rate/all_callable.tsv.gz \
--gff data/rate/gene_lists/expressed_genes.gff3 \
--annotation_table data/references/annotation_table.txt.gz \
--outname data/rate/ka_ks/syn_nonsyn_expressed_callables.tsv
```

looks good - done in 8 min - though the outfile is 420k lines long! going to
perhaps have to filter this down using grep or something similar

now for the Perrineau for completeness' sake:

```bash
time python analysis/rate/callable_genes_degeneracy.py \
--callables_table data/rate/all_callable.tsv.gz \
--gff data/rate/gene_lists/perrineau_genes.gff3 \
--annotation_table data/references/annotation_table.txt.gz \
--outname data/rate/ka_ks/syn_nonsyn_perrineau_callables.tsv
```

also worked like a charm! 

## 12/7/2021

today - actually getting Ka/Ks values for these three gene sets

probably the best way to go at this is just create 'combined' files that contain
the following columns:

for genome wide:

```
sample scaffold nonsyn_muts syn_muts nonsyn_callables syn_callables
```

for gene sets:

```
gene_name nonsyn_muts syn_muts nonsyn_callables syn_callables
```

in the case of gene sets, the gene counts are _not_ sample specific - so the
callable values will need to be summed (since the gene counts are also summed
across samples) 

the genome one is small enough to join quickly in R:

```R
# in data/rate/ka_ks
library(tidyverse)

counts = read_tsv('syn_nonsyn_counts_genome.tsv', col_types = cols())
callables = read_tsv('degen_callables_lookup.tsv', col_types = cols())

final = left_join(counts, callables, by = c('sample', 'scaffold'))

write_tsv(final, 'all_genome.tsv')
```

looks good - now to do this for the gene sets. given that the counts aren't sample specific,
there needs to be some way to combine counts in some sensible way as opposed to
storing 400k+ lines in memory

I can think of a quick and dirty way come to think of it:

```python
import csv
from tqdm import tqdm
import numpy as np

callables = {}

with open('data/rate/ka_ks/syn_nonsyn_perrineau_callables.tsv', 'r', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    last_gene = None
    keys = ['fold0', 'fold2', 'fold3', 'fold4', 'nonsyn', 'syn']
    for line in tqdm(reader):
        current_gene = line['gene_name']
        if current_gene != last_gene:
            last_gene = current_gene
            callables[current_gene] = np.array([float(line[k]) for k in keys])
        elif current_gene == last_gene:
            current_vals = np.array([float(line[k]) for k in keys])
            callables[current_gene] = np.add(callables[current_gene], current_vals) # done in a second!

with open('data/rate/ka_ks/all_perrineau.tsv', 'w') as f_out:
    with open('data/rate/ka_ks/syn_nonsyn_perrineau_counts.tsv', 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = ['gene_name', 'nonsyn_muts', 'syn_muts']
        fieldnames.extend(keys)
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for line in tqdm(reader):
            full_key = [k for k in callables.keys() if line['gene_name'] in k][0]
            values = callables[full_key]
            out_dict = {keys[i]: values[i] for i, colname in enumerate(keys)}
            for key, value in line.items():
                out_dict[key] = value
            writer.writerow(out_dict) # also done in a second! 

# redoing for expressed genes
callables = {}

with open('data/rate/ka_ks/syn_nonsyn_expressed_callables.tsv', 'r', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    last_gene = None
    keys = ['fold0', 'fold2', 'fold3', 'fold4', 'nonsyn', 'syn']
    for line in tqdm(reader):
        current_gene = line['gene_name']
        if current_gene != last_gene:
            last_gene = current_gene
            callables[current_gene] = np.array([float(line[k]) for k in keys])
        elif current_gene == last_gene:
            current_vals = np.array([float(line[k]) for k in keys])
            callables[current_gene] = np.add(callables[current_gene], current_vals)

with open('data/rate/ka_ks/all_expressed.tsv', 'w') as f_out:
    with open('data/rate/ka_ks/syn_nonsyn_expressed_counts.tsv', 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = ['gene_name', 'nonsyn_muts', 'syn_muts']
        fieldnames.extend(keys)
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for line in tqdm(reader):
            full_key = [k for k in callables.keys() if line['gene_name'] in k][0]
            values = callables[full_key]
            out_dict = {keys[i]: values[i] for i, colname in enumerate(keys)}
            for key, value in line.items():
                out_dict[key] = value
            writer.writerow(out_dict)

```

update - this is breaking for the expressed gene since `data/references/mRNA.gff3` doesn't have
any of the NP genes for whatever reason

```bash
cd data/references/
grep -c 'NP_' genes.GFF3 # 163
grep -c 'NP_' mRNA.GFF3 # 0
```

going to do this quick and dirty addition:

```bash
# in data/references
grep 'NP_' genes.GFF3 | grep -v 'CDS' | grep -v 'exon' | cat mRNA.GFF3 - > mRNA_w_organelles.GFF3
grep '^mtMinus' genes.GFF3 | grep -v 'CDS' | grep -v 'exon' >> mRNA_w_organelles.GFF3
```

recreating the filtered GFF first:

```python
import csv
from tqdm import tqdm

with open('data/rate/gene_lists/expressed_genes.tsv', 'r', newline='') as f:
    gene_list = [l[0] for l in csv.reader(f, delimiter='\t')]

with open('data/rate/gene_lists/expressed_genes.gff3', 'w') as f_out:
    writer = csv.writer(f_out, delimiter='\t')
    with open('data/references/mRNA_w_organelles.GFF3', 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in tqdm(reader):
            start, end = int(line[3]), int(line[4])
            name_dict = {l.split('=')[0]: l.split('=')[1] for l in line[-1].split(';')} 
            if any([gene_name in gene_list for gene_name in name_dict.values()]):
                writer.writerow(line)
```

and finally, creating the lookup again:

```bash
time python analysis/rate/callable_genes_degeneracy.py \
--callables_table data/rate/all_callable.tsv.gz \
--gff data/rate/gene_lists/expressed_genes.gff3 \
--annotation_table data/references/annotation_table.txt.gz \
--outname data/rate/ka_ks/syn_nonsyn_expressed_callables.tsv
```

going to redo the combining somewhat blindly:

```python
import csv
from tqdm import tqdm
import numpy as np

callables = {}

with open('data/rate/ka_ks/syn_nonsyn_expressed_callables.tsv', 'r', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    last_gene = None
    keys = ['fold0', 'fold2', 'fold3', 'fold4', 'nonsyn', 'syn']
    for line in tqdm(reader):
        current_gene = line['gene_name']
        if current_gene != last_gene:
            last_gene = current_gene
            callables[current_gene] = np.array([float(line[k]) for k in keys])
        elif current_gene == last_gene:
            current_vals = np.array([float(line[k]) for k in keys])
            callables[current_gene] = np.add(callables[current_gene], current_vals)

with open('data/rate/ka_ks/all_expressed.tsv', 'w') as f_out:
    with open('data/rate/ka_ks/syn_nonsyn_expressed_counts.tsv', 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = ['gene_name', 'nonsyn_muts', 'syn_muts']
        fieldnames.extend(keys)
        writer = csv.DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for line in tqdm(reader):
            try:
                full_key = [k for k in callables.keys() if line['gene_name'] in k][0]
            except:
                tqdm.write(line['gene_name']) # should only print out <30 ADF genes
                continue
            values = callables[full_key]
            out_dict = {keys[i]: values[i] for i, colname in enumerate(keys)}
            for key, value in line.items():
                out_dict[key] = value
            writer.writerow(out_dict)
```
    
broke cause apparently the callables table doesn't have mtMinus on it...

given that no expressed gene with any actual mutations is actually on mtMinus,
I'm going to get the code to outright skip mtMinus with a quick modification in 
`callable_genes_degeneracy` 

rerunning and looks like we're good! 














