
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

alright, next steps: rate calculations! 


