
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

