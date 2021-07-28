
## 24/5/2021

today - getting started on spectrum analyses for MA and salt adaptation,
since these datasets are already done

6 mut rates for MA and saltAd are already done, but need to get
the triplet mutation rates

relevant files:

```
# annotations/output
SL1_haplotype_shared_salt_muts.annotated.txt
SL1_haplotype_shared_salt_indels.annotated.txt
SL2_haplotype_shared_salt_muts.annotated.txt
SL2_haplotype_shared_salt_indels.annotated.txt
unique_salt_muts.annotated.txt
unique_salt_indels.annotated.txt

# VCFs in genotyping/unified_genotyper/UG_diploid/
# the HC ones don't have invariant sites

# /research/projects/chlamydomonas/bgi_full_MA/mutation_calls
final.curated_muts.coord_sorted.txt.gz

# VCFs in bgi_full_MA/vcfs/ with strain specific folders
```

so for each line, need to generate a file along the lines of:

```
strain   chromosome   position   mut   orig_triplet   new_triplet
CC1373   chromosome_2 200        C>A   GCA            GAA 
```

which means all I need from the mut describer columns is the strain,
chr:pos, and mut (`orig>derived`) - and then the surrounding two
sites I can get from the VCF 

in the case of the original MA lines, the VCFs have samples numbered
by generation, and I need to get the 'latest' generation to be sure
(e.g. in case of doublemuts) - although if doublemuts happen in separate
generations that would represent two separate triplets...

probably going to run under the assumption that doublemuts happen in the
same gen, but after the fact I can check for adjacent mutations - there
are probably so few that I can then manually sort those

```bash
# getting data ready
mkdir -p data/spectrum_context/MA
mkdir -p data/spectrum_context/adaptation

# in MA folder
ln -sv /research/projects/chlamydomonas/bgi_full_MA/mutation_calls/final.curated_muts.coord_sorted.txt.gz* .
for sample in 1373 1952 2342 2344 2931 2937; do
    ln -sv /research/projects/chlamydomonas/bgi_full_MA/vcfs/CC${sample}/CC${sample}.vcf.gz* .
done

mkdir -p data/spectrum_context/triplets/
```

getting started on a script called `get_triplets.py`

## 26/5/2021

first pass finally done - let's give it a whirl on the MA files:

```bash
mkdir -p data/spectrum_context/triplets/

time python analysis/spectrum_context/get_triplets.py \
--fname data/spectrum_context/MA/final.curated_muts.coord_sorted.txt \
--vcf data/spectrum_context/MA/CC1373.vcf.gz \
--out data/spectrum_context/triplets/CC1373_MA.tsv
```

looks good - repeating for other samples:

```bash
for sample in 1952 2342 2344 2931 2937; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/spectrum_context/MA/final.curated_muts.coord_sorted.txt \
    --vcf data/spectrum_context/MA/CC${sample}.vcf.gz \
    --out data/spectrum_context/triplets/CC${sample}_MA.tsv;
done
```

seems something broke in 2937 specifically - let's try that again

```bash
time python analysis/spectrum_context/get_triplets.py \
--fname data/spectrum_context/MA/final.curated_muts.coord_sorted.txt \
--vcf data/spectrum_context/MA/CC2937.vcf.gz \
--out data/spectrum_context/triplets/CC2937_MA.tsv
```

looks like there was one strange site with a 'double record' in the VCF -
just modified `get_triplets` to skip instead of breaking

next up - the adaptation lines

## 31/5/2021

time for triplets for the adaptation lines - will need the mut describer
outputs and their corresponding VCFs

from earlier in the log:

```
# annotations/output
SL1_haplotype_shared_salt_muts.annotated.txt
SL1_haplotype_shared_salt_indels.annotated.txt
SL2_haplotype_shared_salt_muts.annotated.txt
SL2_haplotype_shared_salt_indels.annotated.txt
unique_salt_muts.annotated.txt
unique_salt_indels.annotated.txt

# VCFs in alignments/genotyping/unified_genotyper/UG_diploid/
# the HC ones don't have invariant sites
```

looks like only the uniques will work with the current script, since
the shared muts are present in more than one sample - going to need
to handle these differently since their mut describer output doesn't
list a mutant sample to then get surrounding sequence from

```bash
# after symlinking files into data/spectrum_context/adaptation
# had to make a copy and remove metadata from the first line

time python analysis/spectrum_context/get_triplets.py \
--fname data/spectrum_context/adaptation/unique_salt_muts.txt \
--vcf data/spectrum_context/adaptation/salt_aligned_salt_UG_diploid.vcf.gz \
--out data/spectrum_context/triplets/unique_salt_muts.tsv
```

## 9/6/2021

to do on the triplets front:
- saltMA triplets
- are we considering shared mutations? 

point 1 is probably easier to deal with so let's start with that -
need to update `get_triplets.py` to handle the sample setup for the saltMA
lines correctly

by the looks of the code, this should work as is for the pairs files - let's
give it a go with CC1373 -

```bash
time python analysis/spectrum_context/get_triplets.py \
--fname data/mutations/mut_describer/muts_described.final.tsv \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--out trip_test.tsv
```

looks good! going to have to run a loop for this though, and I think
I'll have to handle `DL41_46` separately

let's do all the others first:

```bash
for sample in 1373 1952 2342 2344 2931 2935 2937; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/mutations/mut_describer/muts_described.final.tsv \
    --vcf data/alignments/genotyping/UG/pairs/CC${sample}_samples.vcf.gz \
    --out data/spectrum_context/triplets/CC${sample}_saltMA.tsv;
done
```

looks good - now for the DL lines:

```bash
for sample in 40 51 53 55 57 58; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/mutations/mut_describer/muts_described.final.tsv \
    --vcf data/alignments/genotyping/UG/pairs/DL${sample}_samples.vcf.gz \
    --out data/spectrum_context/triplets/DL${sample}_saltMA.tsv;
done
```

the two SL:

```bash
for sample in 27 29; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/mutations/mut_describer/muts_described.final.tsv \
    --vcf data/alignments/genotyping/UG/pairs/SL${sample}_samples.vcf.gz \
    --out data/spectrum_context/triplets/SL${sample}_saltMA.tsv;
done
```

combining the MA and saltMA files into two separate 'master' files:

```R
# in data/spectrum_context/triplets/
library(tidyverse)
library(fs)

salt_fnames = dir_ls('.', regexp = '\\w+saltMA\\.tsv')
MA_fnames = dir_ls('.', regexp = '\\w+_MA\\.tsv')

d_salt = map_dfr(salt_fnames, read_tsv, col_types = cols())
d_MA = map_dfr(MA_fnames, read_tsv, col_types = cols())

write_tsv(d_salt, 'salt_all.tsv')
write_tsv(d_MA, 'MA_all.tsv')
```

and I'm going to move everything except these two files and `unique_salt_muts.tsv`
into a new subfolder, `triplets_indiv`

```bash
mkdir -p triplets_indiv/ # in same folder
mv -v *MA.tsv triplets_indiv/
```

## 18/6/2021

to do:

- estimate and plot 6 mutation rates (AC, AG, AT, CA, CG, CT)
    - MA
    - saltMA
    - adaptation lines
- estimate and plot 32 triplet rates
    - MA
    - saltMA
    - adaptation lines

the triplet data are all good to go, so I need to review how best to present
that - but in the meantime, the 6 mutation rates is doable as is

need to see if I can find the data used to generate the original Fig 5 since that'd
make life a lot easier, but if not that might need to be recalculated from the mut tables

update: no dice - let's do this the long way

need to calculate obs/exp, where 

exp subs = base composition * prob(base changing to mut base) [33%] * num mutations

going to do this in an Rmd file - `mut_base_spectrum.Rmd`

## 18/7/2021

coming back to this a month later - just need to update the triplet
files to include 0 and 5 in the same names

rerunning - starting with SL so I can check whether it worked quickly

```bash
# SL
for sample in 27 29; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/mutations/mut_describer/muts_described.final.tsv \
    --vcf data/alignments/genotyping/UG/pairs/SL${sample}_samples.vcf.gz \
    --out data/spectrum_context/triplets/triplets_indiv/SL${sample}_saltMA.tsv;
done

# also adding SL26 since I have that now
time python analysis/spectrum_context/get_triplets.py \
--fname data/mutations/mut_describer/muts_described.final.tsv \
--vcf data/alignments/genotyping/UG/pairs/SL26_samples.vcf.gz \
--out data/spectrum_context/triplets/triplets_indiv/SL26_saltMA.tsv;

# looks good! now for CC and DL
# CC
for sample in 1373 1952 2342 2344 2931 2935 2937; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/mutations/mut_describer/muts_described.final.tsv \
    --vcf data/alignments/genotyping/UG/pairs/CC${sample}_samples.vcf.gz \
    --out data/spectrum_context/triplets/triplets_indiv/CC${sample}_saltMA.tsv;
done

# DL
for sample in 40 51 53 55 57 58; do
    time python analysis/spectrum_context/get_triplets.py \
    --fname data/mutations/mut_describer/muts_described.final.tsv \
    --vcf data/alignments/genotyping/UG/pairs/DL${sample}_samples.vcf.gz \
    --out data/spectrum_context/triplets/triplets_indiv/DL${sample}_saltMA.tsv;
done
```

combining the files:

```R
# in data/spectrum_context/triplets/
library(tidyverse)
library(fs)

salt_fnames = dir_ls('triplets_indiv/', regexp = '\\w+saltMA\\.tsv')
d_salt = map_dfr(salt_fnames, read_tsv, col_types = cols())
write_tsv(d_salt, 'salt_all.tsv')
```

## 28/7/2021

regenerating the triplet dataset for the updated adaptation dataset - first,
need to convert csv to tsv

```bash
# in data/prev/
sed 's/,/\t/g' all_mutations_w_shared_hmmIBD_corrected_FPKM.csv > adaptation_muts.tsv
```

and now for triplets:

```bash
time python analysis/spectrum_context/get_triplets.py \
--fname data/prev/adaptation_muts.tsv \
--vcf data/spectrum_context/adaptation/salt_aligned_salt_UG_diploid.vcf.gz \
--out data/spectrum_context/triplets/adaptation_all.tsv
```











