
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
