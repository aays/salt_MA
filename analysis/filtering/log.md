
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

