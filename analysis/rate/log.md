
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


