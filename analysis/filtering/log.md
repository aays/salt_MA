
## 9/7/2020

today: writing a script for first pass mutation filtering

will be working with UG files for SNMs

first pass:

```bash
time python3.5 analysis/filtering/filter_candidate_muts.py \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--out_format vcf \
--out test_out.vcf
```

looks good, and runs reasonably fast - now to create tables
of mutations:

```bash
mkdir -p data/mutations
mkdir -p data/mutations/first_pass_pairs

time for fname in data/alignments/genotyping/UG/pairs/*vcf.gz; do
    base=$(basename ${fname} _samples.vcf.gz)
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf ${fname} \
    --out_format table \
    --out data/mutations/first_pass_pairs/${base}_UG_pairs.txt;
done
```

done in less than two hours - not bad - time for some very quick counts:

```bash
wc -l data/mutations/first_pass_pairs/*
```

what in the world??

```
     15 data/mutations/first_pass_pairs/CC1373_UG_pairs.txt
      3 data/mutations/first_pass_pairs/CC1952_UG_pairs.txt
      6 data/mutations/first_pass_pairs/CC2342_UG_pairs.txt
     13 data/mutations/first_pass_pairs/CC2344_UG_pairs.txt
     12 data/mutations/first_pass_pairs/CC2931_UG_pairs.txt
      5 data/mutations/first_pass_pairs/CC2935_UG_pairs.txt
     14 data/mutations/first_pass_pairs/CC2937_UG_pairs.txt
      4 data/mutations/first_pass_pairs/DL40_UG_pairs.txt
   9569 data/mutations/first_pass_pairs/DL41_UG_pairs.txt
   7949 data/mutations/first_pass_pairs/DL46_UG_pairs.txt
     53 data/mutations/first_pass_pairs/DL51_UG_pairs.txt
     10 data/mutations/first_pass_pairs/DL53_UG_pairs.txt
     19 data/mutations/first_pass_pairs/DL55_UG_pairs.txt
     23 data/mutations/first_pass_pairs/DL57_UG_pairs.txt
     13 data/mutations/first_pass_pairs/DL58_UG_pairs.txt
     44 data/mutations/first_pass_pairs/SL27_UG_pairs.txt
      1 data/mutations/first_pass_pairs/SL29_UG_pairs.txt
  17753 total
```

by the eye test, it seems that DL41 and DL46 share a lot of mutations,
but with the 0 and 5 alleles reversed - this means there could have potentially
been a lab switch

let's see if this still holds for the HC files:

```bash
for sample in DL41 DL46; do
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf data/alignments/genotyping/HC/pairs/${sample}_samples.vcf.gz \
    --out_format table \
    --out data/mutations/first_pass_pairs/${sample}_HC_pairs.txt;
done
```

same story - look at these first few lines:

```
fname   chrom   pos     ref     alt     gt_bases        gt_quals        gt_depths
DL41_samples    chromosome_1    93952   G       ['T']   ['T/T', 'G/G']  [51.0, 60.0]    [17, 20]
DL41_samples    chromosome_1    340161  C       ['G']   ['G/G', 'C/C']  [90.0, 63.0]    [30, 21]
DL41_samples    chromosome_1    469649  G       ['A']   ['A/A', 'G/G']  [33.0, 42.0]    [11, 14]
DL41_samples    chromosome_1    607236  G       ['A']   ['A/A', 'G/G']  [75.0, 35.0]    [25, 12]

fname   chrom   pos     ref     alt     gt_bases        gt_quals        gt_depths
DL46_samples    chromosome_1    93952   G       ['T']   ['G/G', 'T/T']  [48.0, 42.0]    [16, 14]
DL46_samples    chromosome_1    340161  C       ['G']   ['C/C', 'G/G']  [54.0, 30.0]    [18, 10]
DL46_samples    chromosome_1    386230  C       ['CG']  ['CG/CG', 'C/C']        [42.0, 39.0]    [14, 13]
DL46_samples    chromosome_1    607236  G       ['A']   ['G/G', 'A/A']  [54.0, 42.0]    [18, 14]
```

## 14/7/2020

next up - need to get distribution of quality scores
and see how many we're losing < GQ30

also going to be doing a bit of 'pair switching' to figure
out this DL41/46 issue - see alignment log for details

first - create a 'database' of mutations without the GQ filter
and the purity filter (<2 reads of opposing allele) - output
as VCF

then, iterate through these VCFs to create bins of GQ and get 
distribution of GQ scores

## 18/7/2020

okay, actually getting started on this

will need to potentially use the `DL41_46` type switched-pairs files (see
alignment log) for those files, since we've established that as they are
I can't really trust the candidate mutations in those files

going to augment the filtering script to allow for exclusion of certain filters,
and then do a test run:

```bash
time python3.5 analysis/filtering/filter_candidate_muts.py \
--vcf data/alignments/genotyping/UG/pairs/CC1373_samples.vcf.gz \
--gq 0 \
--verbose_level 1 \
--out_format table \
--out test_gq.txt
```

seems this led to a lot of matches! 

```
[saltMA] completed search for CC1373_samples.vcf.gz
[saltMA] found 9469 matches over 107863739 sites.
```

before repeating this over all samples, I'm going to move the DL41 + 46
pairs into the UG pairs folder, and move the current ones into `DL_test`

```bash
mv -v data/alignments/genotyping/UG/DL_test/*vcf data/alignments/genotyping/UG/pairs # adding _samples at the end of these names
mv -v data/alignments/genotyping/UG/pairs/DL41* data/alignments/genotyping/UG/DL_test
mv -v data/alignments/genotyping/UG/pairs/DL46* data/alignments/genotyping/UG/DL_test

mkdir -p data/mutations/gq_tests
mkdir -p data/mutations/gq_tests/vcfs
mkdir -p data/mutations/gq_tests/all_candidate

time for fname in data/alignments/genotyping/UG/pairs/*vcf.gz; do
    base=$(basename ${fname} _samples.vcf.gz)
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf ${fname} --gq 0 --verbose_level 1 --out_format vcf \
    --out data/mutations/gq_tests/vcfs/${base}_UG_all.vcf
done
```

alright - now to create files from this that contain site, ref, alt, gt quals, and gt depths

I could have done this directly with the `table` format, but I think having these
'all candidate VCFs' might be handy in the long run

just a quick script to generate the above text files:

```python
import os
import sys
import csv
from tqdm import tqdm
from cyvcf2 import VCF

fname = sys.argv[-2]
outname = sys.argv[-1]

with open(outname, 'w') as f:
    fieldnames = ['fname', 'chrom', 'pos', 'ref', 'alt', 'qual',
    'gt1', 'gt2', 'gt1_qual', 'gt2_qual', 'gt1_depth', 'gt2_depth']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    v = VCF(fname)
    for record in tqdm(v):
        out_dict = dict.fromkeys(fieldnames)
        out_dict['fname'] = os.path.basename(fname).replace('_UG_all.vcf.gz', '')
        out_dict['chrom'] = record.CHROM
        out_dict['pos'] = record.POS
        out_dict['ref'], out_dict['alt'] = record.REF, record.ALT
        out_dict['qual'] = record.QUAL
        out_dict['gt1'], out_dict['gt2'] = record.gt_bases
        out_dict['gt1_qual'], out_dict['gt2_qual'] = record.gt_quals
        out_dict['gt1_depth'], out_dict['gt2_depth'] = record.gt_depths
        writer.writerow(out_dict)


```

running on one file as a test:

```bash
time python3.5 analysis/filtering/gq_test.py \
data/mutations/gq_tests/vcfs/CC1373_UG_all.vcf.gz \
data/mutations/gq_tests/all_candidate/CC1373_all.txt
```

looks good - now to run on all:

```bash
time for fname in data/mutations/gq_tests/vcfs/*.vcf.gz; do
    base=$(basename ${fname} _UG_all.vcf.gz)
    echo ${base}
    time python3.5 analysis/filtering/gq_test.py \
    ${fname} \
    data/mutations/gq_tests/all_candidate/${base}_all.txt
done
```

joining all the files together:

```bash
touch data/mutations/gq_tests/all_candidate/gq_all_samples.txt
grep '^fname' data/mutations/gq_tests/all_candidate/CC1373_all.txt >> gq_all_samples.txt # header
for fname in data/mutations/gq_tests/all_candidate/*_all.txt; do
    grep -v '^fname' ${fname} >> data/mutations/gq_tests/all_candidate/gq_all_samples.txt;
done
```

## 19/7/2020

today:

- export `gq_all_samples.txt` and analyse distribution of GQs in an Rmd


















