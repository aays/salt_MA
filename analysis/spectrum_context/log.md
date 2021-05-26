
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
--most_recent --out data/spectrum_context/triplets/CC1373_MA.tsv
```


