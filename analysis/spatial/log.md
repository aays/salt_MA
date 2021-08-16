
## 18/7/2021

two broad to-do items here:

1. get intermutation distances for saltMA lines
2. get proportion of new mutations in genome annotations (Fig 2)

intermutation distances - going to do multiple versions of these files:

- stratified by sample group (CC/DL/SL) and line (0/5)
- stratified by line (0/5) irrespective of sample
- not stratified at all

I think I can make a catch all script that computes distances for all three
by caching the mut describer output - still not thrilled about holding it in
memory, but it's not super big fortunately

these also need to be stratified into silent and nonsilent sites regardless
of how I shake it - need to make sure that's baked in - could potentially
do two lookups, one entirely for silent and one for non-silent

in the ms, silent sites were defined as 

- intergenic DNA
- introns
- 4D sites in CDS

wait - I don't need silent sites here - these are MA lines! that dramatically
simplifies the script at least

first attempt:

```bash
time python analysis/spatial/intermutation_distances.py \
--mut_table data/mutations/mut_describer/muts_described.final.tsv \
--stratify line \
--out intermut_test.tsv
```

wait (again) - this script is working, but it's super clunky
and I could literally just do this in R using `dplyr::lag()` 
after some arranging/grouping on the mut describer output - going
to do that instead since the output of this script would have headed 
there regardless

## 31/7/2021

next up - need to generate annotation-specific counts of callables
across all samples

cols should be:

```
annotation genome_count callable_total sample_count # maybe also frac_callable avg_callable
```

quick python script:

```python
import csv
from tqdm import tqdm
import pysam
from datetime import datetime

annotations = ['genic', 'exonic', 'intronic', 'intergenic', 
    'utr5', 'utr3', 'fold0', 'fold4', 'fold2', 'fold3', 'CDS', 'mRNA']
scaffolds = [f'chromosome_{i}' for i in range(1, 18)]
scaffolds.extend([f'scaffold_{i}' for i in range(18, 55)])
scaffolds.extend(['mtDNA', 'cpDNA'])

with open('data/rate/annotation_callables.tsv', 'w') as f:
    fieldnames = ['scaffold', 'annotation', 'genome_count', 'callable_total', 'sample_count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    

    for ant in annotations:
        for scaffold in scaffolds:
            out_dict = {}

            ant_reader = pysam.TabixFile('data/references/annotation_table.txt.gz')
            header = ant_reader.header[-1].split('\t')
            idx = header.index(ant)

            lookup = ''
            chrom_reader = ant_reader.fetch(scaffold)
            for record in tqdm(chrom_reader, desc=f'{scaffold}_{ant}'):
                if record.split('\t')[idx] == '1':
                    lookup += '1'
                else:
                    lookup += '0'

            out_dict['scaffold'] = scaffold
            out_dict['annotation'] = ant
            out_dict['genome_count'] = lookup.count('1')

            callables_reader = pysam.TabixFile('data/rate/all_callable.tsv.gz')
            out_dict['sample_count'] = len(callables_reader.header[0].split('\t')) - 2 # ignore chrom/pos cols
            callables_chrom = callables_reader.fetch(scaffold)

            callables_count = 0
            for i, record in tqdm(enumerate(callables_chrom), desc=f'{scaffold}_{ant}_callables'):
                if lookup[i] == '1':
                    record = record.split('\t')[2:]
                    callables_count += record.count('1')
                else:
                    continue
            out_dict['callable_total'] = callables_count
            writer.writerow(out_dict)
                
```

## 1/8/2021

ran overnight for nearly two hours but the deed is done! 

once more, but for specific samples - going to do long format, since I can
always transpose that in R 

```python
with open('data/rate/annotation_callables_samples.tsv', 'w') as f:
    fieldnames = ['sample', 'scaffold', 'annotation', 'genome_count', 'callable_total', 'sample_count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    samples = pysam.TabixFile('data/rate/all_callable.tsv.gz').header[0].split('\t')[2:]

    for ant in annotations:
        for scaffold in scaffolds:
            out_dict = {}

            ant_reader = pysam.TabixFile('data/references/annotation_table.txt.gz')
            header = ant_reader.header[-1].split('\t')
            idx = header.index(ant)

            lookup = ''
            chrom_reader = ant_reader.fetch(scaffold)
            for record in tqdm(chrom_reader, desc=f'{scaffold}_{ant}'):
                if record.split('\t')[idx] == '1':
                    lookup += '1'
                else:
                    lookup += '0'

            for sample in tqdm(samples, desc=f'samples_{scaffold}_{ant}'):
                out_dict['sample'] = sample
                out_dict['scaffold'] = scaffold
                out_dict['annotation'] = ant
                out_dict['genome_count'] = lookup.count('1')

                callables_reader = pysam.TabixFile('data/rate/all_callable.tsv.gz')
                out_dict['sample_count'] = len(callables_reader.header[0].split('\t')) - 2 # ignore chrom/pos cols
                callables_chrom = callables_reader.fetch(scaffold)
                sample_idx = samples.index(sample)

                callables_count = 0
                for i, record in tqdm(enumerate(callables_chrom), desc=f'{scaffold}_{ant}_callables'):
                    if lookup[i] == '1':
                        record = record.split('\t')[2:][sample_idx]
                        callables_count += record.count('1')
                    else:
                        continue
                out_dict['callable_total'] = callables_count
                writer.writerow(out_dict)
                print(datetime.now())
                    
```

## 2/8/2021

so that took a full day! 

but I forgot to include total sites in there, which is needed for the
denominator in the distribution calculations... going to do that
per sample

```python
import csv
from tqdm import tqdm
import pysam
from datetime import datetime

scaffolds = [f'chromosome_{i}' for i in range(1, 18)]
scaffolds.extend([f'scaffold_{i}' for i in range(18, 55)])
scaffolds.extend(['mtDNA', 'cpDNA'])

with open('data/rate/total_callables_samples.tsv', 'w') as f:
    fieldnames = ['sample', 'scaffold', 'annotation', 'genome_count', 'callable_total', 'sample_count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    samples = pysam.TabixFile('data/rate/all_callable.tsv.gz').header[0].split('\t')[2:]
    ant = 'total'

    for scaffold in scaffolds:
        out_dict = {}

        for sample in tqdm(samples, desc=f'samples_{scaffold}_{ant}'):
            out_dict['sample'] = sample
            out_dict['scaffold'] = scaffold
            out_dict['annotation'] = ant # total
            out_dict['genome_count'] = 0

            callables_reader = pysam.TabixFile('data/rate/all_callable.tsv.gz')
            out_dict['sample_count'] = len(callables_reader.header[0].split('\t')) - 2 # ignore chrom/pos cols
            callables_chrom = callables_reader.fetch(scaffold)
            sample_idx = samples.index(sample)

            callables_count = 0
            for i, record in tqdm(enumerate(callables_chrom), desc=f'{scaffold}_{ant}_callables'):
                out_dict['genome_count'] += 1
                record = record.split('\t')[2:][sample_idx]
                callables_count += record.count('1')
            out_dict['callable_total'] = callables_count
            writer.writerow(out_dict)
            print(datetime.now())
```
                
## 3/8/2021

today - generate annotation callable tables for Perrineau and expressed gene sets

need to also do this for the MA dataset - looks like there's a summarised
sitegens table in `bgi_full_MA/callable_sites/total_site_gens.txt.gz`

gonna be grouping by samples here to speed things up

```python
import csv
import gffutils
from tqdm import tqdm
import pysam

annotations = ['intronic', 'utr5', 'utr3', 'CDS'] # don't need intergenic here

# expressed genes
db = gffutils.FeatureDB('data/rate/gene_lists/expressed_genes.db', keep_order=True)
total_genes = len([gene for gene in db.all_features()])

with open('data/rate/gene_lists/expressed_annotation_callables.tsv', 'w') as f:
    fieldnames = ['sample_factor', 'gene_name', 'scaffold', 'start',
        'end', 'annotation', 'count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    ant_reader = pysam.TabixFile('data/references/annotation_table.txt.gz')
    annotation_header = ant_reader.header[-1].split('\t')

    callables_table = pysam.TabixFile('data/rate/all_callable.tsv.gz')
    samples = callables_table.header[0].split('\t')[2:]
    sample_groups = {
        'CC_0': [i for i, s in enumerate(samples) if s.startswith('CC') and s.endswith('0')],
        'CC_5': [i for i, s in enumerate(samples) if s.startswith('CC') and s.endswith('5')],
        'DL_0': [i for i, s in enumerate(samples) if s.startswith('DL') and s.endswith('0')],
        'DL_5': [i for i, s in enumerate(samples) if s.startswith('DL') and s.endswith('5')],
        'SL_0': [i for i, s in enumerate(samples) if s.startswith('SL') and s.endswith('0')],
        'SL_5': [i for i, s in enumerate(samples) if s.startswith('SL') and s.endswith('5')]
        }

    for gene in tqdm(db.all_features(), desc='genes', total=total_genes):
        scaffold, start, end = gene.chrom, gene.start, gene.end

        for ant in annotations:
            ant_idx = annotation_header.index(ant)
            lookup = ''
            for record in ant_reader.fetch(scaffold, start, end):
                if record.split('\t')[ant_idx] == '1':
                    lookup += '1'
                else:
                    lookup += '0'
            
            for sample_group in sample_groups.keys():
                idx_list = sample_groups[sample_group]
                callables_count = 0
                callables_region = callables_table.fetch(scaffold, start, end)
                for i, record in enumerate(callables_region):
                    if lookup[i] == '1':
                        record = [c for i, c in enumerate(record.split('\t')[2:]) if i in idx_list]
                        callables_count += record.count('1')
                    else:
                        continue
                out_dict = {
                    'sample_factor': sample_group, 'gene_name': gene.attributes['ID'],
                    'scaffold': scaffold, 'start': start, 'end': end,
                    'annotation': ant, 'count': callables_count}
                writer.writerow(out_dict)
```

ran for 1 hr 45 min (longer than it needed to given I accidentally left 
intergenic in there...) - repeating for Perrineau genes by switching the gff and outfile

update: Perrineau's done in 15 minutes - now to transform these before exporting:

```R
library(tidyverse)

# in data/rate/gene_lists
expressed = read_tsv('expressed_annotation_callables.tsv', col_types = cols())
perrineau = read_tsv('perrineau_annotation_callables.tsv', col_types = cols())

expressed_summarised = expressed %>%
    group_by(sample_factor, annotation) %>%
    summarise(callables_total = sum(count))

perrineau_summarised = perrineau %>%
    group_by(sample_factor, annotation) %>%
    summarise(callables_total = sum(count))

write_tsv(expressed_summarised, 'expressed_callables_summarised.tsv')
write_tsv(perrineau_summarised, 'perrineau_callables_summarised.tsv')

```

getting MA callables - symlinked `bgi_full_MA/callable_sites/total_site_gens.txt.gz`
into `data/rate`

```python
import csv
import gffutils
from tqdm import tqdm
import pysam

annotations = ['intronic', 'utr5', 'utr3', 'CDS'] # don't need intergenic here

# expressed genes
db = gffutils.FeatureDB('data/rate/gene_lists/expressed_genes.db', keep_order=True)
total_genes = len([gene for gene in db.all_features()])

with open('data/rate/gene_lists/expressed_ma_callables.tsv', 'w') as f:
    fieldnames = ['sample_factor', 'gene_name', 'scaffold', 'start',
        'end', 'annotation', 'count']
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    ant_reader = pysam.TabixFile('data/references/annotation_table.txt.gz')
    annotation_header = ant_reader.header[-1].split('\t')

    callables_table = pysam.TabixFile('data/rate/total_site_gens.txt.gz')

    for gene in tqdm(db.all_features(), desc='genes', total=total_genes):
        scaffold, start, end = gene.chrom, gene.start, gene.end

        for ant in annotations:
            ant_idx = annotation_header.index(ant)
            lookup = ''
            for record in ant_reader.fetch(scaffold, start, end):
                if record.split('\t')[ant_idx] == '1':
                    lookup += '1'
                else:
                    lookup += '0'
            
            callables_count = 0
            callables_region = callables_table.fetch(scaffold, start, end)
            for i, record in enumerate(callables_region):
                if lookup[i] == '1':
                    callables_count += float(record.split('\t')[-1])
                else:
                    continue
            out_dict = {
                'sample_factor': 'ma_original', 'gene_name': gene.attributes['ID'],
                'scaffold': scaffold, 'start': start, 'end': end,
                'annotation': ant, 'count': callables_count}
            writer.writerow(out_dict)
```

took 20 minutes for expressed, and then just 3 when repeated for Perrineau

quick transformations:

```R
library(tidyverse)

# in data/rate/gene_lists
expressed = read_tsv('expressed_ma_callables.tsv', col_types = cols())
perrineau = read_tsv('perrineau_ma_callables.tsv', col_types = cols())

expressed_summarised = expressed %>%
    group_by(sample_factor, annotation) %>%
    summarise(callables_total = sum(count))

perrineau_summarised = perrineau %>%
    group_by(sample_factor, annotation) %>%
    summarise(callables_total = sum(count))

write_tsv(expressed_summarised, 'expressed_ma_callables_summarised.tsv')
write_tsv(perrineau_summarised, 'perrineau_ma_callables_summarised.tsv')

```
