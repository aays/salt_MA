
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

## 22/7/2020

today:

start adapting mutation filtering script to take into account 'ancestral'
mutation state

let's find one from the CC1373 first pass mutations - the first
one is at `chromosome_2:6902447`:

```python
>>> from cyvcf2 import VCF
>>> fname = 'data/alignments/genotyping/UG/combined/CC1373_combined.vcf.gz'
>>> v = VCF(fname)
>>> recs = [rec for rec in v('chromosome_2:6902447-6902448')]
>>> recs
[Variant(chromosome_2:6902447 G/A), Variant(chromosome_2:6902448 C/)]
>>> recs[0]
Variant(chromosome_2:6902447 G/A)
>>> rec = recs[0]
>>> rec.genotypes
[[0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [1, 1, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False], [0, 0, False]]
>>> rec.__str__()
'chromosome_2\t6902447\t.\tG\tA\t201.41\t.\tAC=2;AF=0.043;AN=46;BaseQRankSum=0.657;DP=812;Dels=0;ExcessHet=0.7918;FS=7.993;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=59.62;MQ0=0;MQRankSum=-0.152;QD=29.97;ReadPosRankSum=0.859;SOR=0.612;set=variant2\tGT:AD:DP:GQ:PL\t0/0:.:12:.:.\t0/0:.:15:.:.\t0/0:.:31:.:.\t0/0:.:21:.:.\t0/0:.:23:.:.\t0/0:.:27:.:.\t0/0:.:13:.:.\t0/0:.:21:.:.\t0/0:.:23:.:.\t0/0:.:32:.:.\t0/0:.:24:.:.\t0/0:.:27:.:.\t0/0:.:33:.:.\t0/0:.:25:.:.\t0/0:.:32:.:.\t0/0:.:31:.:.\t0/0:.:63:.:.\t1/1:0,14:14:42:560,42,0\t0/0:21,0:21:63:0,63,840\t0/0:.:71:.:.\t0/0:.:93:.:.\t0/0:.:76:.:.\t0/0:.:84:.:.\n'
```

this looks like a bona fide new mutation - nice

going to be hard to code this up - here's what needs to be considered:

1. identify the indices at which we have our 0 and 5 samples in `genotypes`
(since they're not always at the end of the file) 
2. for now, assume that all other samples should differ from the novel allele between the 0 and 5 
(should maybe write code to specifically check for this)

it's key that the samples be allowed to differ from one another as long as none of them
have the same novel allele between the pairs - although might be good to have
a counter for these because it shouldn't happen *too* often (that or the differences
should be across groups of lines and not within them)

alright, done a first pass - going to generate table formatted files and then compare
against the first pass mutations

```bash
mkdir -p data/mutations/ancestral

# test run
time python3.5 analysis/filtering/filter_candidate_muts.py \
--vcf data/alignments/genotyping/UG/combined/CC1373_combined.vcf.gz \
--gq 30 \
--out_format table \
--vcf_type combined \
--verbose_level 1 \
--purity_filter \ # for consistency
--out data/mutations/ancestral/CC1373_UG_combined.txt
```

## 24/7/2020

so the above looks good, but I'm going to put this on hold for a bit - will need to
revisit once the specific logic for each of the sample groups is figured out

for now - need to generate tables of mutations that are `>GQ20` without the purity filter:

```bash
mkdir -p data/mutations/gq_tests/pairs_20/

time for fname in data/alignments/genotyping/UG/pairs/*vcf.gz; do
    base=$(basename ${fname} _samples.vcf.gz)
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf ${fname} --gq 20 --vcf_type pairs --verbose_level 1 \
    --out_format table --out data/mutations/gq_tests/pairs_20/${base}_GQ20.txt
done
```

## 4/8/2020

finally getting back on this - next up, need to install IGV on the server in some way
to queue up screenshots for each of the mutations of interest 

installing igvtools (all platform compatible):

```bash
cd ~/apps
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.9.zip
time unzip IGV_2.8.9.zip
```

update - need a java 8 compatible version - java 11 not supported on server

trying IGV 2.4.10 instead (I have 2.4.0 on my Mac, and apparently
2.4 is java 8 compatible):

```bash
wget https://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.19.zip
wget https://data.broadinstitute.org/igv/projects/downloads/2.4/igvtools_2.4.19.zip
time unzip IGV_2.4.19.zip
time unzip IGVTools_2.4.19
```

so apparently `igv.sh` will open the GUI to do whatever operation
is asked of it, even if a batch script is provided - need to install
`xvfb` to run it in command line only mode (since we can't use X11 otherwise)

installing xvfb:

```bash
wget http://ftp.xfree86.org/pub/XFree86/4.6.0/binaries/Linux-ix86-glibc20/Xvfb.tgz
tar zxvf Xvfb.tgz
```

update: this is totally not working for weird C reasons that are above my pay grade and
cognitive ability to figure out

instead - going to try this with X11 forwarding - this works, although it's a bit buggy and laggy,
and could be worth giving a shot to

going to set up a symlink to 'IGV.sh' in the project folder before giving that a shot
with the `--batch` setting

test batch script (after making a new `screenshots` dir):

```
new
snapshotDirectory data/mutations/screenshots
load data/alignments/bam/CC1373_0.bam index=data/alignments/bam/CC1373_0.bam
load data/alignments/bam/CC1373_5.bam index=data/alignments/bam/CC1373_5.bam
genome data/references/chlamy.5.3.w_organelles_mtMinus.fasta
maxPanelHeight 500
goto chromosome_2:6944200-6944300
snapshot CC1373_chr2_6944200-6944300.png
exit
```

here goes nothing:

```bash
./bin/igv.sh --batch=test.batch
```

so this works... technically. it's agonizingly slow (X11 forwarding apparently
generally is) and will run for >10 minutes on a single bam, let alone two

wait - there's an xvfb solution in [this anonymous blog post](http://thedusseldorfer.blogspot.com/2013/09/remotely-plotting-with-igv-even-without.html)!

here goes:

```bash
cd ~/apps
wget ftp://ftp.xfree86.org/pub/XFree86/4.8.0/binaries/Linux-x86_64-glibc23/Xfnts.tgz
wget ftp://ftp.xfree86.org/pub/XFree86/4.8.0/binaries/Linux-x86_64-glibc23/Xvfb.tgz
tar zxvf Xfnts.tgz
tar zxvf Xvfb.tgz
cd bin/
./Xvfb :1 -nolisten tcp -fp /home/hasans11/apps/lib/X11/fonts/misc/ # starts X server
export DISPLAY=:1
```

and this works! 

some notes

- display port 0 was taken, but 1 worked just fine
- there are some warnings about opening a security policy file, but overall it still works fine
- the absolute path to the fonts dir needs to be given to the Xvfb command
    - doesn't work with a relative path
- I did change `apps/lib/X11/xorg.conf` based on the manual solution found [here](https://access.redhat.com/solutions/320563)
    - this was prior to the absolute/relative path fix above, though, so I don't know if it actually helped much

next up - write a Python script that will use the table output files to create temp IGV batch scripts

things to figure out:

- file naming convention
- how much room to give around the candidate mutation

otherwise, the operation seems pretty straightforward - for a given input file,
write the batch script to a temp text file and use that to run IGV on the command line
using `subprocess` - this will require `DISPLAY` to be set and an Xvfb server to be
running of course

will output both SNMs and indels, but these need to be named differently and there
should be an option to just output one and not the other (ie an option that modifies
the batch script)

alright, script (`mutation_snapshots.py`) looks good - test this tomorrow morning
(lol it's 1 am... I got really into this)

## 5/8/2020

giving the script a shot:

```bash
mkdir -p data/mutations/screenshots/CC1373

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_20/CC1373_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/CC1373_0.bam data/alignments/bam/CC1373_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots/CC1373
```

looks good! 

queuing this for all of them:

```bash
time while read sample; do
    mkdir -p data/mutations/screenshots/${sample}
    time python3.5 analysis/filtering/mutation_snapshots.py \
    --fname data/mutations/gq_tests/pairs_20/${sample}_GQ20.txt \
    --igv bin/igv.sh \
    --reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --bam_files data/alignments/bam/${sample}_0.bam data/alignments/bam/${sample}_5.bam \
    --flank_size 50 \
    --outdir data/mutations/screenshots/${sample}
done < data/alignments/fastq/symlinks/samples.txt
```

all done - took 47 minutes too, which isn't bad at all

some issues - need to get the DL41/46 ones going manually, since they weren't
included in the `samples` file, while SL26 (predictably) failed since there
aren't any muts given the missing line

```bash
mkdir -p data/mutations/screenshots/DL41_46
mkdir -p data/mutations/screenshots/DL46_41

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_20/DL41_46_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL41_0.bam data/alignments/bam/DL46_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots/DL41_46

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_20/DL46_41_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL46_0.bam data/alignments/bam/DL41_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots/DL46_41
```

looks good - next up, need to make a csv with all the mutations in the
`pairs_20` dir and make... a spreadsheet out of them?? something to make
comments (or even some 'approved' mark on) for later reference

creating the csv:

```bash
head -n 1 $(ls -tr | head -n 1) > mutations_GQ20.tsv
for fname in *txt; do
    tail -n +2 ${fname} >> mutations_GQ20.tsv;
done
```

making snp only version:

```python
>>> import csv
>>> fname = 'mutations_GQ20.tsv'
>>> outname = 'mutations_GQ20_snps.tsv'
>>> from tqdm import tqdm
  2 with open(fname, 'r', newline='') as f_in:
  3     reader = csv.DictReader(f_in, delimiter='\t')
  4     fieldnames = reader.fieldnames
  5     with open(outname, 'w', newline='') as f_out:
  6         writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
  7         writer.writeheader()
  8         for record in tqdm(reader):
  9             alt = eval(record['alt'])
 10             max_alt_length = max([len(s) for s in alt])
 11             if max_alt_length == 1:
 12                 writer.writerow(record)
```

all done - now to go through these screenshots one by one for any SNPs where 20 < min GQ < 30

to make this a bit easier, going to actually split the GQ values into two columns:

```python
>>> import csv
>>> fname = 'mutations_GQ20_snps.tsv'
>>> outname = 'mutations_GQ20_snps_split.tsv'
>>> from tqdm import tqdm
>>> with open(fname, 'r', newline='') as f_in:
  2     reader = csv.DictReader(f_in, delimiter='\t')
  3     fieldnames = reader.fieldnames
  4     fieldnames.extend(['GQ1', 'GQ2'])
  5     with open(outname, 'w', newline='') as f_out:
  6         writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
  7         writer.writeheader()
  8         for record in tqdm(reader):
  9             quals = eval(record['gt_quals'])
 10             gt1_qual, gt2_qual = quals
 11             line_out = record
 12             removed = line_out.pop('gt_quals')
 13             line_out['GQ1'] = gt1_qual
 14             line_out['GQ2'] = gt2_qual
 15             writer.writerow(record)
```

wait - the above code (and the IGV screenshots) didn't consider
that some indel records have `len(ref)` > 1, while the alts are just
single bases - goddamn it

```python
import csv
fname = 'mutations_GQ20.tsv'
outname = 'mutations_GQ20_snps.tsv'
from tqdm import tqdm
with open(fname, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    fieldnames.extend(['GQ1', 'GQ2'])
    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for record in tqdm(reader):
            alt = eval(record['alt'])
            max_alt_length = max([len(s) for s in alt])
            if max_alt_length == 1 and len(record['ref']) == 1:
                line_out = record
                line_out['GQ1'], line_out['GQ2'] = eval(record['gt_quals'])
                writer.writerow(line_out)

```

running the IGV script again after fixing the snp bug:

```bash
export DISPLAY=:1
./bin/Xvfb :1 -nolisten tcp -fp /home/hasans11/apps/lib/X11/fonts/misc

# in separate shell with above server running
time while read sample; do
    if [ ${sample} != "DL41" ] && [ ${sample} != "DL46" ]; then
        mkdir -p data/mutations/screenshots/${sample}
        time python3.5 analysis/filtering/mutation_snapshots.py \
        --fname data/mutations/gq_tests/pairs_20/${sample}_GQ20.txt \
        --igv bin/igv.sh \
        --reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --bam_files data/alignments/bam/${sample}_0.bam data/alignments/bam/${sample}_5.bam \
        --flank_size 50 \
        --outdir data/mutations/screenshots/${sample}
    fi
done < data/alignments/fastq/symlinks/samples.txt

mkdir -p data/mutations/screenshots/DL41_46
mkdir -p data/mutations/screenshots/DL46_41

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_20/DL41_46_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL41_0.bam data/alignments/bam/DL46_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots/DL41_46

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_20/DL46_41_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL46_0.bam data/alignments/bam/DL41_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots/DL46_41
```

## 7/8/2020

now that this is done - will open the tsv in Excel (shudder)
and start manually inspecting + giving a pass/fail to
all mutations where 20 < min(GQ) < 30

## 16/8/2020

so one week, one dead laptop, and one new laptop later - time to do this
for indels using the HC VCFs

first, need to create a similar tsv to the ones in `pairs_20` containing
locations of indels and other info:

```bash
mkdir -p data/mutations/gq_tests/HC_pairs_20

time for fname in data/alignments/genotyping/HC/pairs/*vcf.gz; do
    if [ "${sample}" != "DL41" ] && [ "${sample}" != "DL46" ]; then
        base=$(basename ${fname} _samples.vcf.gz)
        time python3.5 analysis/filtering/filter_candidate_muts.py \
        --vcf ${fname} --gq 20 --vcf_type pairs --verbose_level 1 \
        --out_format table --out data/mutations/gq_tests/HC_pairs_20/${base}_GQ20.txt;
    fi;
done
```

need to also create `DL41_46` and `DL46_41` VCFs using HC:

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL41_0.bam \
-I data/alignments/bam/DL46_5.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/pairs/DL41_46.vcf

time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL46_0.bam \
-I data/alignments/bam/DL41_5.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/pairs/DL46_41.vcf

bgzip data/alignments/genotyping/HC/pairs/DL41_46.vcf
bgzip data/alignments/genotyping/HC/pairs/DL46_41.vcf
tabix -p vcf data/alignments/genotyping/HC/pairs/DL41_46.vcf.gz
tabix -p vcf data/alignments/genotyping/HC/pairs/DL46_41.vcf.gz

# including the weird DLs
time python3.5 analysis/filtering/filter_candidate_muts.py \
--vcf data/alignments/genotyping/HC/pairs/DL41_46.vcf.gz \
--gq 20 --vcf_type pairs --verbose_level 1 \
--out_format table --out data/mutations/gq_tests/HC_pairs_20/DL41_46_GQ20.txt

time python3.5 analysis/filtering/filter_candidate_muts.py \
--vcf data/alignments/genotyping/HC/pairs/DL46_41.vcf.gz \
--gq 20 --vcf_type pairs --verbose_level 1 \
--out_format table --out data/mutations/gq_tests/HC_pairs_20/DL46_41_GQ20.txt
```

once this is done, prepping a tsv for the spreadsheet:

```bash
cd data/mutations/gq_tests/HC_pairs_20/
head -n 1 $(ls -tr | head -n 1) > mutations_HC_GQ20.tsv
for fname in *txt; do
    tail -n +2 ${fname} >> mutations_HC_GQ20.tsv;
done
```

splitting the GQ values into two columns:

```python
import csv
from tqdm import tqdm
fname = 'mutations_HC_GQ20.tsv'
outname = 'mutations_HC_GQ20_split.tsv'
with open(fname, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    fieldnames.extend(['GQ1', 'GQ2'])
    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for record in tqdm(reader):
            quals = eval(record['gt_quals'])
            gt1_qual, gt2_qual = quals
            line_out = record
            removed = line_out.pop('gt_quals')
            line_out['GQ1'], line_out['GQ2'] = gt1_qual, gt2_qual
            writer.writerow(record)
```

going to leave both snps and indels in here - it'll be useful
to actually compare this with the UG snps to look for overlaps/differences

let's make an indels version as well though:

```python
import csv
from tqdm import tqdm
fname = 'mutations_HC_GQ20_split.tsv'
outname = 'mutations_HC_GQ20_indels.tsv'
with open(fname, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for record in tqdm(reader):
            alt = eval(record['alt'])
            max_alt_length = max([len(s) for s in alt])
            if max_alt_length > 1 or len(record['ref']) > 1:
                line_out = record
                writer.writerow(line_out)
```

and now for the screenshotting:

```bash
export DISPLAY=:1
./bin/Xvfb :1 -nolisten tcp -fp /home/hasans11/apps/lib/X11/fonts/misc

# in separate shell with above server running
mkdir -p data/mutations/HC_screenshots

time while read sample; do
    if [ "${sample}" != "DL41" ] && [ "${sample}" != "DL46" ]; then
        mkdir -p data/mutations/HC_screenshots/${sample}
        time python3.5 analysis/filtering/mutation_snapshots.py \
        --fname data/mutations/gq_tests/HC_pairs_20/${sample}_GQ20.txt \
        --igv bin/igv.sh \
        --reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --bam_files data/alignments/bam/${sample}_0.bam data/alignments/bam/${sample}_5.bam \
        --flank_size 50 --outdir data/mutations/HC_screenshots/${sample}
    fi
done < data/alignments/fastq/symlinks/samples.txt

# weird DL samples
mkdir -p data/mutations/HC_screenshots/DL41_46
mkdir -p data/mutations/HC_screenshots/DL46_41

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/HC_pairs_20/DL41_46_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL41_0.bam data/alignments/bam/DL46_5.bam \
--flank_size 50 --outdir data/mutations/HC_screenshots/DL41_46

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/HC_pairs_20/DL46_41_GQ20.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL46_0.bam data/alignments/bam/DL41_5.bam \
--flank_size 50 --outdir data/mutations/HC_screenshots/DL46_41
```

tomorrow: 

- import spreadsheet into google sheets and then get to inspecting indels
- do quick check of overlap b/w UG and HC variants (overall numbers, and then shared sites etc)

## 17/8/2020

checking overlap b/w UG and HC variants

```bash
$ wc -l data/mutations/gq_tests/pairs_20/*txt
   22 data/mutations/gq_tests/pairs_20/CC1373_GQ20.txt
   12 data/mutations/gq_tests/pairs_20/CC1952_GQ20.txt
    7 data/mutations/gq_tests/pairs_20/CC2342_GQ20.txt
   14 data/mutations/gq_tests/pairs_20/CC2344_GQ20.txt
   24 data/mutations/gq_tests/pairs_20/CC2931_GQ20.txt
   13 data/mutations/gq_tests/pairs_20/CC2935_GQ20.txt
   17 data/mutations/gq_tests/pairs_20/CC2937_GQ20.txt
   82 data/mutations/gq_tests/pairs_20/DL40_GQ20.txt
   27 data/mutations/gq_tests/pairs_20/DL41_46_GQ20.txt
   29 data/mutations/gq_tests/pairs_20/DL46_41_GQ20.txt
   86 data/mutations/gq_tests/pairs_20/DL51_GQ20.txt
   25 data/mutations/gq_tests/pairs_20/DL53_GQ20.txt
   31 data/mutations/gq_tests/pairs_20/DL55_GQ20.txt
   39 data/mutations/gq_tests/pairs_20/DL57_GQ20.txt
   41 data/mutations/gq_tests/pairs_20/DL58_GQ20.txt
  117 data/mutations/gq_tests/pairs_20/SL27_GQ20.txt
    1 data/mutations/gq_tests/pairs_20/SL29_GQ20.txt
  587 total

$ wc -l data/mutations/gq_tests/HC_pairs_20/*txt
   31 data/mutations/gq_tests/HC_pairs_20/CC1373_GQ20.txt
   12 data/mutations/gq_tests/HC_pairs_20/CC1952_GQ20.txt
   12 data/mutations/gq_tests/HC_pairs_20/CC2342_GQ20.txt
   44 data/mutations/gq_tests/HC_pairs_20/CC2344_GQ20.txt
   46 data/mutations/gq_tests/HC_pairs_20/CC2931_GQ20.txt
   15 data/mutations/gq_tests/HC_pairs_20/CC2935_GQ20.txt
   36 data/mutations/gq_tests/HC_pairs_20/CC2937_GQ20.txt
   14 data/mutations/gq_tests/HC_pairs_20/DL40_GQ20.txt
   31 data/mutations/gq_tests/HC_pairs_20/DL41_46_GQ20.txt
   36 data/mutations/gq_tests/HC_pairs_20/DL46_41_GQ20.txt
  112 data/mutations/gq_tests/HC_pairs_20/DL51_GQ20.txt
   32 data/mutations/gq_tests/HC_pairs_20/DL53_GQ20.txt
   34 data/mutations/gq_tests/HC_pairs_20/DL55_GQ20.txt
   51 data/mutations/gq_tests/HC_pairs_20/DL57_GQ20.txt
   47 data/mutations/gq_tests/HC_pairs_20/DL58_GQ20.txt
   20 data/mutations/gq_tests/HC_pairs_20/SL27_GQ20.txt
    3 data/mutations/gq_tests/HC_pairs_20/SL29_GQ20.txt
  576 total
```

immediate observations -
- there's some minor variation b/w the two, but HC has less total
- SL27 - which has huge stretches of mutations in the UG dataset - has much fewer in the HC dataset

need to do a more in depth investigation of overlap in python later - but onto
looking at screenshots for now

nts: accidentally let DL41 and DL46 into the table formatted file - had to remove these
with grep (since `DL41_46` and `DL46_41` should be used instead) 

## 17/9/2020

been a while! 

let's start with a quick overlap analysis:

```r
library(tidyverse)
d_ug <- read_tsv('pairs_20/mutations_GQ20_snps.tsv', col_types = cols())
d_hc <- read_tsv('HC_pairs_20/mutations_HC_GQ20_snps.tsv', col_types = cols()) %>%
    mutate(fname = ifelse(
        str_detect(fname, 'DL4[0-9]_4[0-9]'), paste0(fname, '_samples'), fname))
shared <- d_ug %>% select(fname, chrom, pos) %>%
    inner_join(select(d_hc, fname, chrom, pos))
dim(shared) # 137 rows
dim(d_ug); dim(d_hc) # 340 and 194 respectively

# checking variation across files:
shared %>% group_by(fname) %>% summarise(n = n())
# A tibble: 16 x 2
   fname               n
   <chr>           <int>
 1 CC1373_samples     10
 2 CC1952_samples      5
 3 CC2342_samples      5
 4 CC2344_samples      7
 5 CC2931_samples      5
 6 CC2935_samples      8
 7 CC2937_samples      7
 8 DL40_samples        6
 9 DL41_46_samples     9
10 DL46_41_samples     2
11 DL51_samples       29
12 DL53_samples        9
13 DL55_samples       11
14 DL57_samples        9
15 DL58_samples       10
16 SL27_samples        6
```

seems all samples that had at least one mutation had a shared mutation
between both samples - which is something

checking this against how many mutations there were in total per
filename:

```r
shared %>% group_by(fname) %>% summarise(n = n()) %>%
left_join(d_ug %>% group_by(fname) %>% summarise(n_ug = n())) %>%
left_join(d_hc %>% group_by(fname) %>% summarise(n_hc = n()))

   fname               n  n_ug  n_hc
   <chr>           <int> <int> <int>
 1 CC1373_samples     10    10    30
 2 CC1952_samples      5     7    11
 3 CC2342_samples      5     5    11
 4 CC2344_samples      7     8    43
 5 CC2931_samples      5     6    45
 6 CC2935_samples      8     8    14
 7 CC2937_samples      7    10    35
 8 DL40_samples        6    77    13
 9 DL41_46_samples     9    10    30
10 DL46_41_samples     2     3    35
11 DL51_samples       29    33   111
12 DL53_samples        9    10    31
13 DL55_samples       11    12    33
14 DL57_samples        9    11    50
15 DL58_samples       10    15    46
16 SL27_samples        6   115    19

```

seems that save for the mess that is SL27, HC has more mutations
per sample overall and that there are at least a few samples
(esp in the case of the CC strains) where all UG mutations
were also called by HC

wait shit - this HC file isn't just SNPs like the UG file is

gotta make that real quick - code from earlier:

```python
import csv
fname = 'mutations_HC_GQ20.tsv'
outname = 'mutations_HC_GQ20_snps.tsv'
from tqdm import tqdm
with open(fname, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    fieldnames.extend(['GQ1', 'GQ2'])
    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for record in tqdm(reader):
            alt = eval(record['alt'])
            max_alt_length = max([len(s) for s in alt])
            if max_alt_length == 1 and len(record['ref']) == 1:
                line_out = record
                line_out['GQ1'], line_out['GQ2'] = eval(record['gt_quals'])
                writer.writerow(line_out)
```

updated the dimensions in the bit above, but here's the new filename
report:

```r
shared %>% group_by(fname) %>% count()
# A tibble: 16 x 2
# Groups:   fname [16]
   fname               n
   <chr>           <int>
 1 CC1373_samples     10
 2 CC1952_samples      5
 3 CC2342_samples      5
 4 CC2344_samples      7
 5 CC2931_samples      5
 6 CC2935_samples      7
 7 CC2937_samples      7
 8 DL40_samples        6
 9 DL41_46_samples     9
10 DL46_41_samples     2
11 DL51_samples       29
12 DL53_samples        9
13 DL55_samples       11
14 DL57_samples        9
15 DL58_samples       10
16 SL27_samples        6
```

seems a CC2935 mutation is now no longer shared? weird that that
was eliminated when I removed indels... will have to look into that one again

here's the file comparison again:

```r
# A tibble: 16 x 4
# Groups:   fname [16]
   fname               n  n_ug  n_hc
   <chr>           <int> <int> <int>
 1 CC1373_samples     10    10    13
 2 CC1952_samples      5     7     5
 3 CC2342_samples      5     5     6
 4 CC2344_samples      7     8    14
 5 CC2931_samples      5     6    14
 6 CC2935_samples      7     8     8
 7 CC2937_samples      7    10    14
 8 DL40_samples        6    77     7
 9 DL41_46_samples     9    10     9
10 DL46_41_samples     2     3     3
11 DL51_samples       29    33    37
12 DL53_samples        9    10    11
13 DL55_samples       11    12    13
14 DL57_samples        9    11    12
15 DL58_samples       10    15    13
16 SL27_samples        6   115    15
```

overall numbers look a lot closer - which is really promising! especially
considering I'm looking for exact matches here

gonna have to look into this mystery CC2935 mutation next...

```r
> left_join(d_ug %>% select(fname, chrom, pos), d_hc_full %>% select(fname, chrom, pos)) %>%
+ filter(fname == 'CC2935_samples')
# A tibble: 8 x 3
  fname          chrom             pos
  <chr>          <chr>           <dbl>
1 CC2935_samples chromosome_6  6423831
2 CC2935_samples chromosome_8  4202087
3 CC2935_samples chromosome_9  7390928
4 CC2935_samples chromosome_11 3562572
5 CC2935_samples chromosome_13 1650911
6 CC2935_samples chromosome_13 3465246
7 CC2935_samples mtDNA             732
8 CC2935_samples mtDNA            8627
> shared %>% filter(fname == 'CC2935_samples')
# A tibble: 7 x 3
  fname          chrom             pos
  <chr>          <chr>           <dbl>
1 CC2935_samples chromosome_8  4202087
2 CC2935_samples chromosome_9  7390928
3 CC2935_samples chromosome_11 3562572
4 CC2935_samples chromosome_13 1650911
5 CC2935_samples chromosome_13 3465246
6 CC2935_samples mtDNA             732
7 CC2935_samples mtDNA            8627
```

figured it out - the 'mutation' at `chromosome_6:6423831` is a SNP
in the UG file, but an indel in the HC file - since I was only looking
at fname, chrom, and pos, R found an overlap 

overall, 137 of the 194 HC mutations being supported by UG is a good sign - 
there are more mutations in UG (340), but I wonder how much that number is
inflated by the SL27 mutations

## 11/10/2020

today: need to add screenshots + create a spreedsheet for mutations
w/ GQ b/w 10 and 20 (well, 19)

given how the filtering script works, I'll have to create a full
dataset with GQ >= 10 and then filter mutations where both samples
are above GQ20 after the fact - ideally I'd keep these all together,
but it doesn't make a lot of sense to re-screenshot a bunch of
files that have already been screenshotted + evaluated

```bash
mkdir data/mutations/gq_tests/pairs_10
mkdir data/mutations/gq_tests/pairs_10/filtered

# this includes the weird DL lines
time for fname in data/alignments/genotyping/UG/pairs/*vcf.gz; do
    base=$(basename ${fname} _samples.vcf.gz)
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf ${fname} --gq 10 --vcf_type pairs --verbose_level 1 \
    --out_format table --out data/mutations/gq_tests/pairs_10/${base}_GQ10.txt
done # takes 2 hours

mkdir -p data/mutations/screenshots_GQ10
```

## 13/10/2020

filtering the GQ10 files to only include new mutations where GQ is b/w 10-19 for one or both
of the calls:

```python
import csv
from tqdm import tqdm
from glob import glob
fnames = glob('data/mutations/gq_tests/pairs_10/*txt')

def file_filter(fname):
    with open(fname, 'r') as f_in:
        reader = csv.DictReader(f_in, delimiter='\t')
        fieldnames = reader.fieldnames
        outname = fname.replace('_10/', '_10/filtered/')
        with open(outname, 'w') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            print(fname)
            for line in tqdm(reader):
                quals = eval(line['gt_quals'])
                if any([qual < 20 for qual in quals]):
                    writer.writerow(line)

for fname in tqdm(fnames):
    file_filter(fname)
```

looks good - screenshot time:

```bash
export DISPLAY=:1
./bin/Xvfb :1 -nolisten tcp -fp /home/hasans11/apps/lib/X11/fonts/misc

# in separate shell with above server running
time while read sample; do
    if [ ${sample} != "DL41" ] && [ ${sample} != "DL46" ]; then
        mkdir -p data/mutations/screenshots_GQ10/${sample}
        time python3.5 analysis/filtering/mutation_snapshots.py \
        --fname data/mutations/gq_tests/pairs_10/filtered/${sample}_GQ10.txt \
        --igv bin/igv.sh \
        --reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --bam_files data/alignments/bam/${sample}_0.bam data/alignments/bam/${sample}_5.bam \
        --flank_size 50 \
        --outdir data/mutations/screenshots_GQ10/${sample}
    fi
done < data/alignments/fastq/symlinks/samples.txt

mkdir -p data/mutations/screenshots_GQ10/DL41_46
mkdir -p data/mutations/screenshots_GQ10/DL46_41

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_10/DL41_46_GQ10.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL41_0.bam data/alignments/bam/DL46_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots_GQ10/DL41_46

time python3.5 analysis/filtering/mutation_snapshots.py \
--fname data/mutations/gq_tests/pairs_10/DL46_41_GQ10.txt \
--igv bin/igv.sh \
--reference data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--bam_files data/alignments/bam/DL46_0.bam data/alignments/bam/DL41_5.bam \
--flank_size 50 \
--outdir data/mutations/screenshots_GQ10/DL46_41
```

done in an hour - now to make a spreadsheet for these mutations:

```bash
# in data/mutations/gq_tests/pairs_10/filtered
head -n 1 $(ls -tr | head -n 1) > mutations_GQ10.tsv
for fname in *txt; do
    tail -n +2 ${fname} >> mutations_GQ10.tsv;
done
```

SNP only version:

```python
import csv
fname = 'mutations_GQ20.tsv'
outname = 'mutations_GQ20_snps.tsv'
from tqdm import tqdm
with open(fname, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    fieldnames.extend(['GQ1', 'GQ2'])
    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for record in tqdm(reader):
            alt = eval(record['alt'])
            max_alt_length = max([len(s) for s in alt])
            if max_alt_length == 1 and len(record['ref']) == 1:
                line_out = record
                line_out['GQ1'], line_out['GQ2'] = eval(record['gt_quals'])
                writer.writerow(line_out)
```

## 23/3/2021

finally back on this! 

first order of business - creating a mutation table from the ones
that we know are good to go - e.g. all the ones above GQ20 and
those that Rob has okayed in the GQ10-20 spreadsheet

there are three sets here:

1. mutations where both were >GQ30 (autopassed)
2. mutations where one mut is 10-20 (most passed, need to go through spreadsheet)
3. mutations where both are 10-20 (only a few passed, need to go through spreadsheet)

for category 2, there are lots in DL40 and SL27 in particular where one mut is
99 and the other is very low - Rob has labelled these are likely being not
real, though a handful are also supported by HC - might need to revisit that
later on

for now though, need to assemble a table of 1, passing muts in 2 (export
from Drive spreadsheet), and passing muts in 3 (export from spreadsheet
as well)

there are some in the spreadsheet for 3 that likely require realignment - 
need to talk to Rob about best way to approach this (getting HC to specifically
realign these regions?) also need to go through ms methods once more to
refresh memory on how these were handled in previous datasets - consistency 
is also good and important! 

## 25/3/2021

today: getting the updated >GQ20 spreadsheet on the server and extracting
the 'confirmed' mutations

just leaving the file (`mutations_GQ20_snps_reviewed.tsv`) in the parent data directory for now

```bash
mkdir -p data/mutations/filtered
```

## 13/4/2021

so I dropped the ball on that - will get back to it at some point, but first I need
to get a count of >GQ30 mutations for my committee report

I need to add the `DL41_46` and `46_41` SNPs to `first_pass_pairs` though - that
hasn't been updated yet

here's what I ran back then:

```bash
time for fname in data/alignments/genotyping/UG/pairs/*vcf.gz; do
    base=$(basename ${fname} _samples.vcf.gz)
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf ${fname} \
    --out_format table \
    --out data/mutations/first_pass_pairs/${base}_UG_pairs.txt;
done
```

might actually rerun this entirely given that I believe these
were done with the purity filter on (and there were a few 
other args added since)

```bash
mkdir -p data/mutations/gq_tests/pairs_30
time for fname in data/alignments/genotyping/UG/pairs/*vcf.gz; do
    base=$(basename ${fname} _samples.vcf.gz)
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf ${fname} --gq 30 --vcf_type pairs --verbose_level 1 \
    --out_format table --out data/mutations/gq_tests/pairs_30/${base}_GQ30.txt
done
```

probably need to also generate 'combined' (eg all ancestral samples + the given pair)
files for DL41/46 - holding off on this for now though, since after mutation
filtering we'll probably only be interested in a subset (and mutation filtering
actually needs to be completed first before I end up focusing on other things...)

done in 2 hours - let's get a count of strains and of high quality mutations:

```bash
wc -l data/mutations/gq_tests/pairs_30/*
```
 
the other thing I need to do is get trimming stats out (to help debug weird
mutation clusters) - since these weren't saved the first time (!) I might need
to rerun that again with the same command I used earlier, but this time
actually explicitly saving the logs...

## 26/4/2021

today - listing out regions that I need to manually realign to try and 'solve'

these are from the variants labelled '2' in `mutations_GQ20_snps.tsv`

```
CC2937 chromosome_5 377471
CC2937 chromosome_9 2314036
DL51 chromosome_11 2143294
DL51 chromosome_11 2143301
DL51 chromosome_11 2143306
DL58 chromosome_2 4729387
```

a few of the '2' variants are actually called correctly in the HC data - will
need to manually amend final dataset to include these

```
CC2931 chromosome_1 915925 # 'solved' with HC indel call
DL58 chromosome_16 4380945 # missing in UG data
```

alright - time for some local realignments - this seems like it'd belong in
the alignment log, but I guess it's still part of filtering overall

from Josianne's draft:

"To verify that candidate mutations in regions with problematic alignment were
not the outcome of misalignment, we realigned the 100 bases around the site
using GATK HaplotypeCaller with the ‘--activeRegionIn’ option"

so let's do that up

```bash
mkdir -p data/mutations/realignment
ln -sv ~/apps/GenomeAnalysisTK-3.3/GenomeAnalysisTK.jar bin/GenomeAnalysisTK.jar

time /usr/bin/java -jar ./bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/CC2937_0.bam \
    -I data/alignments/bam/CC2937_5.bam \
    -L chromosome_5 \
    -ploidy 2 \
    --activeRegionIn chromosome_5:377421-377521 \
    --output_mode EMIT_ALL_SITES \
    --heterozygosity 0.02 \
    --indel_heterozygosity 0.002 \
    -o data/mutations/realignment/test.vcf
```

## 27/4/2021

looks like this made the variant 'disappear' altogether - and there's no variant at all
that actually is different between the two parents in this region at all

there isn't even a proximate indel that (in this case) explains what seemed to have 'created'
the origin variants

let's try automating this with a quick script - `realign_HC.py`

test run:

```
# interval file
CC2937 chromosome_5 377471
```

```bash
time python analysis/filtering/realign_HC.py \
--fname data/mutations/realignment/intervals.txt \
--outdir data/mutations/realignment
```

looks good - will look in the surrounding 4 kb, with the surrounding 100 bp defined
as the active region

updating interval list:

```
CC2937 chromosome_5 377471
CC2937 chromosome_9 2314036
DL51 chromosome_11 2143294
DL51 chromosome_11 2143301
DL51 chromosome_11 2143306
DL58 chromosome_2 4729387
```

with notes taken on the google sheets spreadsheet - long story short,
only the last indel actually shows in the HC data, but it's a heterozygote call

next up - need to download BAMs for review in real IGV for the iffy GQ10 muts
that Rob has marked - also need to get them back on the server to run with this
above script just to see if those indels explain anything

in the file `data/mutations/mutations_GQ10_snps_reviewed.tsv`, all the 
mutations where `pass == 2` need to have bams downloaded for variant review

for the latter case, would be good to create temp filtered bams for just the surrounding 2k
each time for easy download and local review - might need another script for this

also - need to create a 'master list' of mutations (probably in
`data/filtered`) containing

1. all the GQ30+ mutations
2. all the GQ20 mutations that were manually passed
3. all the GQ10 mutations that were manually passed

and then need to google sheets `data/mutations/gq_tests/HC_pairs_20/mutations_HC_GQ20_indels.tsv`
and then manually review + score the short indels using IGV

and then need to finish reviewing the short indels in `mutations_HC_GQ20_indels.tsv` (in
`data/mutations/gq_tests/HC_pairs_20` and also on google sheets) - it seems I started on
this but didn't get especially far

outstanding question - for the long stretches of mutations in DL40 and SL27, do we also
exclude the ones that even HC called (which are much fewer?)

still also need to regenerate the trimming stats and actually save the logs this time! 

## 28/4/2021

today:

1. update 'master list' of mutations containing manually passed mutations
2. create and download filtered bams from questionable GQ10 mutations for IGV review

for point 2 - need to get ~4000 bp surrounding mutation of interest - probably
overkill but it can't hurt

```bash
mkdir -p data/mutations/realignment/filtered_igv

# can use samtools view like so:
samtools view data/alignments/bam/bamfile.bam "chr:start-end" > out.sam
```

## 29/4/2021

so I ended up losing most of yesterday to NSERC stuff... time to get back on this

first, need to create a filtered version of the GQ10 spreadsheet with just
the questionable mutations:

```bash
awk -F '\t' '(NR == 1) || ($11 == 2)' mutations_GQ10_snps_reviewed.tsv > GQ10_intervals.tsv
```

and then as per custom let's write a quick and dirty python script
`export_IGV.py` to generate filtered bams based on this:

```python
time python analysis/filtering/export_IGV.py \
--fname data/mutations/GQ10_intervals.tsv \
--bam_dir data/alignments/bam \
--region_size 4000 \
--outdir data/mutations/realignment/filtered_igv
```

now to download these and have a look at them on IGV locally

## 1/5/2021

today: in addition to continuing the IGV adventures, need to create 
GQ10-30 spreadsheets for HC SNPs as well, and compare against the 'full' UG
datasets - basically asking whether all UG mutations we've accepted are
also in the HC dataset

first - assembling a sheet of 'confirmed' SNPs mutations thus far:

1. all the GQ30+ mutations (`data/mutations/gq_tests/pairs_30`)
2. all the GQ20 mutations that were manually passed (`data/mutations/mutations_GQ20_snps_reviewed.tsv`)
3. all the GQ10 mutations that were manually passed (`data/mutations/mutations_GQ10_snps_reviewed.tsv`)

the third set here is still being manually filtered and will need to be updated, 
I can do a HC check here

```R
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(fs)

gq30_fnames <- dir_ls('data/mutations/gq_tests/pairs_30')
d_30 <- map_dfr(gq30_fnames, read_delim, delim = '\t', col_types = cols())
d_20 <- read_delim('data/mutations/mutations_GQ20_snps_reviewed.tsv') %>%
    filter(pass == 1) %>%
    select(fname, chrom, pos, ref, alt, gt_bases, gt_quals, gt_depths)
all_muts <- bind_rows(d_30, d_20) %>%
    arrange(fname, chrom, pos)
write_tsv(all_muts, 'data/mutations/all_mutations.tsv')

# is every one of these muts in the HC data as well?

d_hc <- read_tsv('data/mutations/gq_tests/HC_pairs_20/mutations_HC_GQ20_snps.tsv', col_types = cols()) %>%
    mutate(fname = ifelse(
        str_detect(fname, 'DL4[0-9]_4[0-9]'), paste0(fname, '_samples'), fname)) %>%
    select(-GQ1, -GQ2)

```

might need to port these files offline for a closer analysis in RStudio
    

## 2/5/2021

wait - I never removed indels from the `pairs_30` dataset...

```bash
cd data/mutations/gq_tests/pairs_30
head -n 1 $(ls -tr | head -n 1) > mutations_GQ30.tsv
for fname in *txt; do 
    tail -n +2 ${fname} >> mutations_GQ30.tsv;
done
```

SNPs only:

```bash
import csv
from tqdm import tqdm
fname = 'mutations_GQ30.tsv'
outname = 'mutations_GQ30_snps.tsv'
with open(fname, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    fieldnames.extend(['GQ1', 'GQ2'])
    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for record in tqdm(reader):
            alt = eval(record['alt'])
            max_alt_length = max([len(s) for s in alt])
            if max_alt_length == 1 and len(record['ref']) == 1:
                line_out = record
                line_out['GQ1'], line_out['GQ2'] = eval(record['gt_quals'])
                writer.writerow(line_out)
```

regenerating the all mutations file:

```R
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(fs)

d_30 <- read_tsv('data/mutations/gq_tests/pairs_30/mutations_GQ30_snps.tsv', col_types = cols())
d_20 <- read_tsv('data/mutations/mutations_GQ20_snps_reviewed.tsv') %>%
    filter(pass == 1) %>%
    select(fname, chrom, pos, ref, alt, gt_bases, gt_quals, gt_depths)
all_muts <- bind_rows(d_30, d_20) %>%
    arrange(fname, chrom, pos)
write_tsv(all_muts, 'data/mutations/all_mutations.tsv')
```

wait - I've confused myself - the GQ20 sheet contains >GQ30 mutations as well - so the real
all mutations file would just be the passed ones in that file

just checked in R as well and the only ones in the 30 file not in the 20 are the ones where
pass != 1

```R
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(fs)

# coulda done this in awk lol
d_20 <- read_tsv('data/mutations/mutations_GQ20_snps_reviewed.tsv') %>%
    filter(pass == 1) %>%
    select(fname, chrom, pos, ref, alt, gt_bases, gt_quals, gt_depths, GQ1, GQ2)
write_tsv(d_20, 'data/mutations/all_mutations.tsv')
```

## 5/5/2021

today - making sense of the original 'mutation describer' and
recreating its functionality

the original mutation table format can be seen in `annotations/mut_list` (in the
main folder) and the mutation describer py2.7 script is at
`annotations/scripts/JL_mut_describer.py`

going to go through this script and 'pseudocode out' what it does:

```python
# args
annotation_table_file
reference_fasta
population_vcf_file # should find example of which VCFs were used and their samples
num_mutated # this is never referenced in the script again
mut_file # ostensibly the tsvs in mut_list
vcf_file # check how this differs from population_vcf_file - latter used for pi calc

# constants
num_unmutated = 1 # why assign this? 
mut_sample_count = collections.defaultdict(int) # dict that defaults to int key
h = 0 # used to ensure header prints only once apparently
ref_dict = SeqIO.to_dict(SeqIO.parse(reference_fasta, 'fasta'))

# main loop
for line in open(mut_file):
    sp = line.split('\t')
    chrom = sp[0]
    pos = int(sp[1])
    mut = mutation.mutation(chrom, pos) # USES ROB'S CUSTOM MUTATION CLASS
    # will need to get mut class working with py3.6 even though it was written in 2.7...

    # populating mut object attributes, apparently - 
    mut.ref = sp[2]
    mut.alt = sp[3]
    mut.qual = float(sp[4])
    mut.MQ = float(sp[5]) # mapping quality
    mut.call_method = sp[8] # idgi - this is just 'QUAL' in the entire file I think?

    # more complex (?) mut object attributes
    samples = mut.samples(vcf_file) # why get samples via mut class instead of pyvcf instance?
    mutation_event = mut.mutation(vcf_file, mut_quality=mut_quality) # exposes VCF info to mut class?
    mut_type = mut.mutation_type(vcf_file) # for 'ref>alt' style formatting I believe
    sample = mut.mutant_sample(vcf_file, mut_quality) # no idea what this does - seems very roundabout

    # calculations and filtering
    purity = mut.purity(vcf_file) # apparently this variable is never used
    if len(mut.purity(vcf_file)) > 0: # why not use purity variable??
        min_purity = min(mut.purity(vcf_file))
    else:
        min_purity = 'n/a'

    # variant stats by way of mut object
    het_count = mut.het_count(vcf_file)
    number_mutants = mut.number_mutants(vcf_file)
    number_genotypes = mut.number_genotypes(vcf_file)
    GTs, ADs = mut.GTs(vcf_file), mut.ADs(vcf_file) 
    DPs, GQs = mut.DPs(vcf_file), mut.GQs(vcf_file)

    # more 'genomic' stats
    gc20 = mut.gc_content(ref_dict, 10) # this should be using sample consensus sequence - I hope
    gc2000 = mut.gc_content(ref_dict, 1000)
    pi1000 = mut.pi(population_vcf_file, 500, ref_dict) # so this is what pop vcf file is for
    annotation_pos = mut.annotation_position(annotation_table_file)
    if annotation_pos.CDS == '1' and mut_type == 'SNP':
        mut.nonsyn_syn = mut.nonsynonymous(annotation_table_file, vcf_file)
        nonsyn_syn = mut.nonsynonymous(annotation_table_file, vcf_file) # same thing but as a local?
    else:
        mut.nonsyn_syn = ('n/a', 'n/a', 'n/a')
        nonsyn_syn = ('n/a', 'n/a', 'n/a') # again w the confusing local assignment
    nonsense = mut.nonsense(vcf_file, annotation_table_file)

    # prepping outfile
    annotation_header = '\t'.join([]) # has all the ant column names
    header = '\t'.join([
        'chromosome', 'position', 'ref', 'alt', 'qual', 'MQ',
        'mutation', 'type', 'mutant_sample', 'call_method', 
        'HQ_calls', 'min_purity', 'het_count', 'number_mutants',
        'number_genotypes', 'gc20', 'gc2000', 'pi1000',
        annotation_header, # ?!?!?!
        'Nonsyn_V_Syn', 'anc_codon', 'mut_codon', 'nonsense',
        "\t".join([str(j)+"_GT" for j in samples]), \
        "\t".join([str(j)+"_AD" for j in samples]), \
        "\t".join([str(j)+"_DP" for j in samples]), \
        "\t".join([str(j)+"_GQ" for j in samples])
    ])

    if h == 0: # ah so this is so that the header only prints the first time
        print(header)
        h += 1

    sys.stdout.write('\t'.join([
        mut.chromosome, mut.position, mut.ref, mut.alt, mut.qual,
        mut.MQ, '>'.join(mutation_event), # ah - so this is a C>T type deal
        mut_type, sample, mut.call_method, num_HQ_calls, # HQ calls was never defined...
        min_purity, het_count, number_mutants, number_genotypes,
        gc20, gc2000, pi1000, 
        annotation_position.output_line(), # corresponding to the header above - feels like overkill
        '\t'.join([str(x) for x in nonsyn_syn]),
        nonsense,
        '\t'.join([str(x) for x in GTs]),
        '\t'.join([str(x) for x in ADs]),
        '\t'.join([str(x) for x in DPs]),
        '\t'.join([str(x) for x in GQs])]) + '\n')

```

this is (mostly) well and good, but before I can proceed with this, need to get
the mutation and the annotation table scripts working with py3.6

`mutation.py` and `annotation_table.py` are both in `/reseach/repos/annotation` -
off the bat - looks like this needs to be converted to py3.6 altogether - breaks pretty quick! 

going to make a copy (NOT symlink) of both `mutation.py` and `annotation_table.py` in 
`.conda/env/work/lib/python3.6/site-packages/`
 
looking like `annotation_table.py` functions as is - figures since it's a fairly
lightweight parser that lets `pysam.Tabixfile` do most of the work 

update: `mutation.py` _should_ be py3 compatible now - in addition to the usual `print`/`exec`
statement fare, had to remove the `string` imports (those methods are now part of 
string methods in py3) and import `ness_vcf` and `annotation_table` at the top level
since the scripts just live 'as is' in `site-packages` (hooray conda!...)

also needed to convert all tabs to spaces - thank goodness for `:retab` in vim

(I may be a little inebriated... I hope this works tomorrow all the same)

tomorrow - need to get started on a new version of mut describer! and consolidate
the SNM (incl. GQ10 manually reviewed mutations) and indels from HC calls

other things - checking on HC in the morning (see alignment log) and also
getting those trimming stats 

## 6/5/2021

today - getting the GQ10 muts and indels on the server, and then getting started
on a new mut describer script

in the indel file, I've marked cases of alignment issues leading to the appearance
of multiple indels as '4', and all other passed indels are '1'

joining muts and making indels file:

```R
library(tidyverse)

d_all <- read_tsv('all_mutations.tsv')
filtered <- read_tsv('mutations_GQ10_snps_reviewed.tsv') %>%
    filter(pass == 1) %>%
    select(fname, chrom, pos, ref, alt, gt_bases, gt_quals, gt_depths, GQ1, GQ2)

combined <- bind_rows(d_all, filtered)

# overwrite earlier all_muts file - equivalent to appending really
write_tsv(combined, path = 'all_mutations.tsv')

indels_all <- read_tsv('mutations_HC_GQ20_indels.tsv') %>%
    filter(pass == 1 | pass == 4) %>%
    select(fname, chrom, pos, ref, alt, gt_bases, gt_depths, GQ1, GQ2) # apparently I forgot the gt_quals column
write_tsv(indels_all, path = 'all_indels.tsv')
```

now for mut describer - before I actually write this script, let's just try
to read in a single mutation and access various `mutation.mutation` attributes to
see if anything breaks 

```python
import annotation_table
import ness_vcf
import mutation
import vcf
from tqdm import tqdm

with open('all_mutations.tsv', 'r') as f:
    for line in f:
        if line.startswith('fname'):
            continue
        x = line.rstrip().split('\t')
        break

# ['CC1373_samples', 'chromosome_5', '1737824', 'G', "['T']", "['T/T', 'G/G']", '[30.0, 63.0]', '[11, 21]', '30', '63']

sample = x[0]
chrom = x[1]
pos = int(x[2])

mut = mutation.mutation(chrom, pos)

# apparently the mut object is empty - need to populate attributes
# will have to write code that accounts for len(ALT) == 2
# also need to eval (ugh) the alt column to use it
```

## 8/5/2021

time to write this script - calling it `mut_describer.py`

notes while I make this -

- need to add MQ and QUAL to the `all_mutations` and `all_indels` file just
for completeness
- need to test behaviour when len(ALT) is 2 - it seems that `mutation.mutation`
is getting the directionality wrong - (if ref is G and samples are A and C, it claims it's a `C>A` mut)
    - in this case, we probably need to look at the ancestors, right...
    - might be a separate analysis involving going through the combined files to infer what
    the ancestral state was - I should count how many times this happen to see whether this is
    feasible to do manually or whether I need to automate it

## 9/5/2021

let's check how often we get more than one alt allele:

```python
>>> with open('all_mutations.tsv', 'r') as f:
...     counter = 0
...     dbl = 0
...     for line in tqdm(f):
...         if line.startswith('fname'):
...             continue
...         sp = line.split('\t')
...         alt = eval(sp[4])
...         if len(alt) > 1:
...             dbl += 1
...         counter += 1
>>> print(dbl, counter)
1 208
```

just the one! looks like this is a much more common occurrence for indels though:


```python
>>> with open('all_indels.tsv', 'r') as f:
...     counter = 0
...     dbl = 0
...     for line in tqdm(f):
...         if line.startswith('fname'):
...             continue
...         sp = line.split('\t')
...         alt = eval(sp[4])
...         if len(alt) > 1:
...             dbl += 1
...         counter += 1
219it [00:00, 54649.72it/s]
>>> print(dbl, counter)
23 218
```

so for now I'm going to write mut describer without accounting for this much -
will deal with these manually I suppose

remaining to do before mut describer is ready:

1. add annotation table headers and fields
2. add MQ and QUAL fields to `all_mutations.tsv` and `all_indels.tsv` since these
need to be in the final file
3. what do I use for `population_vcf_file`? what was used previously? 

code for task 2 above:

```python
import csv
from copy import deepcopy
from cyvcf2 import VCF
from tqdm import tqdm

with open('data/mutations/all_mutations.tsv', 'r') as f_in:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames
    fieldnames.extend(['MQ', 'QUAL'])
    with open('data/mutations/all_muts_corrected.tsv', 'w') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for mut_dict in tqdm(reader):
            out = deepcopy(mut_dict)
            v = VCF(f"data/alignments/genotyping/UG/pairs/{out['fname']}.vcf.gz")
            chrom, pos = out['chrom'], out['pos']
            record = [r for r in v(f"{chrom}:{pos}-{pos}")][0]
            out['MQ'] = record.INFO.get('MQ')
            out['QUAL'] = record.QUAL
            writer.writerow(out)
```

reran for indels as well, just switching filenames and switching VCF folder to HC

task 1 also done, so tomorrow, need to figure out pop vcf file and then it's off
to the races!

other remaining task 'mutation-side' would be to manually check ancestral alleles for
the 'double alt' muts - I expect a lot of the indels in this category to not have
any though...

## 10/5/2021

update: won't be using pi1000 values for this anyways, so making that argument optional

here goes! 

```bash
mkdir -p data/mutations/mut_tables/
mv -v data/mutations/all* data/mutations/mut_tables/

time python analysis/filtering/mut_describer.py \
--fname data/mutations/mut_tables/all_muts_corrected.tsv \
--vcf_path data/alignments/genotyping/UG/pairs/ \
--ant_file data/references/annotation_table.txt.gz \
--ref_fasta data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--outname muts_test.tsv
```

needed to fix a bug in `mutation.py` (`maketrans` is now a str method)

```python
# old
def reverse_complement(sequence):
    return str(sequence)[::-1].translate(maketrans('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB'))

# new
def reverse_complement(sequence):
    t = sequence.maketrans('ACGTNRYKMWS?X.-BDHV', 'TGCANYRMKWS?X.-VHDB')
    return str(sequence)[::-1].translate(t)
```

regenerating the annotation table tbi cause the warning about it being older
is driving me nuts (they're 'equally old', both set to the date server maintenance
reset their timestamps)

also had to modify `mutation` to lower quality to 10

script is working, but need to use combined VCFs to get correct mutation 'direction'

```bash
time python analysis/filtering/mut_describer.py \
--fname data/mutations/mut_tables/all_muts_corrected.tsv \
--vcf_path data/alignments/genotyping/UG/combined/ \
--ant_file data/references/annotation_table.txt.gz \
--ref_fasta data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--outname muts_test_combined.tsv
```

also had to change the default `num_mutated` to 5 (though I need to go back
and look manually at the 'ancestral reappearing mutations' later) as well
as adding a check to see whether the GQ field exists since a lot of the GT/AD/etc
listcomps were breaking without that check (and the prior VCFs don't have GQ on 
invariant sites) 

next up - need to do a CombineVariants run on the `DL41_46` and vice versa samples
before running mut describer on those (since I currently don't have those combined VCFs)

## 16/5/2021

finally done CombineVariants - getting this show on the road: 

```bash
time python analysis/filtering/mut_describer.py \
--fname data/mutations/mut_tables/all_muts_corrected.tsv \
--vcf_path data/alignments/genotying/UG/combined \
--ref_fasta data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--ant_file data/references/annotation_table.txt.gz \
--outname muts_test_combined.tsv
```

finally worked - had to bump up `num_mutated` in the `mutation` class to 10 to get
this to work

remaining things to deal with - 

1. that one record with multiple ALT alleles - needs to be manually examined
2. all records where `num_mutated` is more than 1 - need to see whether these
   cases are inherited muts or recurrent (?) muts! 

also - the output file has " at the start and end of the annotation table columns, in a strange
attempt to work around the `'` in the `feature_names` column

update - R is not reading these headers correctly! I think we need to unpack the annotation
table columns... fixed this with a quick loop in `describe_mut`

next up - manually examining that record with multiple ALT alleles - 
CC1952, `chromosome_7:3938585`

looks like a `C>A` mutation! `G>C` in 1952, and `C>A` in the 0 sample in the salt MA lines

how many mutations had `num_mutated` > 1?

```python
import vcf
import mutation
import csv

with open('data/mutations/mut_tables/all_muts_corrected.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    lines = [r for r in reader]

counter = 0
vcf_dir = 'data/alignments/genotyping/UG/combined/'

for mut in lines:
    x = mutation.mutation(mut['chrom'], int(mut['pos']))
    vcf_name = vcf_dir + mut['fname'].replace('samples', '') + 'combined.vcf.gz'
    gts = x.hq_GTs(vcf_name)
    if min(gts.count(gt) for gt in gts) > 1 or len(list(set(gts))) > 2:
        counter += 1
# 44
# which 44 though -

multi = []
for mut in lines:
    x = mutation.mutation(mut['chrom'], int(mut['pos']))
    vcf_name = vcf_dir + mut['fname'].replace('samples', '') + 'combined.vcf.gz'
    gts = x.hq_GTs(vcf_name)
    if min(gts.count(gt) for gt in gts) > 1 or len(list(set(gts))) > 2:
        multi.append(vcf_name, x])

def describe(vcf_path, mut):
    print(vcf_path, mut.chromosome, mut.position)
    vcf_rec = mut.vcf_record(vcf_path)
    print(vcf_rec.REF, vcf_rec.ALT)
    return [[s['GT'], s['GQ'], s.sample] for s in vcf_rec.samples]

# writing to an outfile
with open('data/mutations/mut_tables/multi_muts.tsv', 'w', newline='') as f:
    fieldnames = [
    'fname', 'chrom', 'pos', 'ref', 'alt']
    vcf_rec = multi[0][1].vcf_record(multi[0][0])
    samples = [s.sample for s in vcf_rec.samples][:-2]
    samples.extend(['fname_0', 'fname_5'])
    cols = [(f'{s}_GT', f'{s}_GQ') for s in samples]
    fieldnames.extend([colname for pair in cols for colname in pair])
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for vcf_path, mut in multi:
        vcf_rec = mut.vcf_record(vcf_path)
        d = {}
        d['fname'] = vcf_path.replace('data/alignments/genotyping/UG/combined/', '').replace('.vcf.gz', '')
        d['chrom'] = mut.chromosome
        d['pos'] = mut.position
        d['ref'] = vcf_rec.REF
        d['alt'] = vcf_rec.ALT
        for sample in vcf_rec.samples:
            if sample.sample.endswith('_0'):
                d['fname_0_GT'] = sample['GT']
                d['fname_0_GQ'] = sample['GQ']
            elif sample.sample.endswith('_5'):
                d['fname_5_GT'] = sample['GT']
                d['fname_5_GQ'] = sample['GQ']
            else:
                d[sample.sample + '_GT'] = sample['GT']
                d[sample.sample + '_GQ'] = sample['GQ']
        writer.writerow(d)
```

that worked first try? wild

going to download these and make notes in a spreadsheet

remaining tasks -

1. check mut describer calls and see if they correlate with the manual checks in the spreadsheet
2. check intersection between HC and UG calls again and get extra HC muts in there
3. run mut describer on indels
4. similar 'multi' analysis on indels

also, smaller task - find previous mut describer outputs as well, for remaking plots etc

## 20/5/2021

so apparently there are a few mutations (eg `chromosome_15:411424`, seen in DL55) that are actually
called by mut describer in an ancestral sample (CC124 in this case) - how do we want to handle these?

this and the few 'recurrent' mutations that I've found in the spreadsheet I made probably need
to be discussed with Rob - but for now I should get mut describer going on the indels and make a
similar spreadsheet

will also need to check tomorrow how well the 'clear cut' calls in the spreadsheet match up with
what mut describer decided itself

```bash
time python analysis/filtering/mut_describer.py \
--fname data/mutations/mut_tables/all_indels_corrected.tsv \
--vcf_path data/alignments/genotying/HC/combined \
--ref_fasta data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--ant_file data/references/annotation_table.txt.gz \
--outname indels_test_combined.tsv
```

shoot - also need to create the `DL41_46` type files with CombineVariants for HC

## 21/5/2021

all done and now functioning as intended - now for the multi muts check:

```bash
import vcf
import mutation
import csv

with open('data/mutations/mut_tables/all_indels_corrected.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    lines = [r for r in reader]

counter = 0
vcf_dir = 'data/alignments/genotyping/UG/combined/'

for mut in lines:
    x = mutation.mutation(mut['chrom'], int(mut['pos']))
    vcf_name = vcf_dir + mut['fname'].replace('samples', '') + 'combined.vcf.gz'
    if '_4' in vcf_name:
        vcf_name = vcf_name.replace('combined.vcf.gz', '_combined.vcf.gz')
    gts = x.hq_GTs(vcf_name)
    if min(gts.count(gt) for gt in gts) > 1 or len(list(set(gts))) > 2:
        counter += 1
# 105! of 218... oh boy

multi = []
for mut in lines:
    x = mutation.mutation(mut['chrom'], int(mut['pos']))
    vcf_name = vcf_dir + mut['fname'].replace('samples', '') + 'combined.vcf.gz'
    if '_4' in vcf_name:
        vcf_name = vcf_name.replace('combined.vcf.gz', '_combined.vcf.gz')
    gts = x.hq_GTs(vcf_name)
    if min(gts.count(gt) for gt in gts) > 1 or len(list(set(gts))) > 2:
        multi.append(vcf_name, x])

# writing to an outfile
with open('data/mutations/mut_tables/multi_indels.tsv', 'w', newline='') as f:
    fieldnames = [
    'fname', 'chrom', 'pos', 'ref', 'alt']
    vcf_rec = multi[20][1].vcf_record(multi[0][0]) # first record was weird
    samples = [s.sample for s in vcf_rec.samples][:-2]
    samples.extend(['fname_0', 'fname_5'])
    cols = [(f'{s}_GT', f'{s}_GQ') for s in samples]
    fieldnames.extend([colname for pair in cols for colname in pair])
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    for vcf_path, mut in multi:
        vcf_rec = mut.vcf_record(vcf_path)
        d = {}
        d['fname'] = vcf_path.replace('data/alignments/genotyping/UG/combined/', '').replace('.vcf.gz', '')
        d['chrom'] = mut.chromosome
        d['pos'] = mut.position
        d['ref'] = vcf_rec.REF
        d['alt'] = vcf_rec.ALT
        for sample in vcf_rec.samples:
            if sample.sample.endswith('_0'):
                d['fname_0_GT'] = sample['GT']
                d['fname_0_GQ'] = sample['GQ']
            elif sample.sample.endswith('_5'):
                d['fname_5_GT'] = sample['GT']
                d['fname_5_GQ'] = sample['GQ']
            else:
                d[sample.sample + '_GT'] = sample['GT']
                d[sample.sample + '_GQ'] = sample['GQ']
        writer.writerow(d)
```

reviewing to-do list:

1. check mut describer calls and see if they correlate with the manual checks in the spreadsheet - SORTA DONE
2. check intersection between HC and UG calls again and get extra HC muts in there
3. run mut describer on indels - DONE
4. similar 'multi' analysis on indels

eyeballing item 1 above, and a few things to note -

- the manual checks don't correlate a lot of the time - mut describer tends to pick the wt/D/S lines
when it comes to muts that appear to be 'recurrent'
- mut describer currently ignores anything with multiple alts, so we're still missing those
    - ^ this also means a lot of the indels are missing in the mut describer output!  

## 22/5/2021

it occurs to me that `mutation.py` can handle multiple alts just fine - I should
update `mut_describer` to account for this and rerun the muts at minimum to see those calls

```bash
time python analysis/filtering/mut_describer.py \
--fname data/mutations/mut_tables/all_muts_corrected.tsv \
--vcf_path data/alignments/genotying/UG/combined \
--ref_fasta data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--ant_file data/references/annotation_table.txt.gz \
--outname muts_test_combined_new.tsv
```

## 24/5/2021

calls for multi-alt samples still defaulting to wt/D/S strains in some cases - need to
get the 'pedigree' for these lines figured out with Rob asap before redoing this for indels

in the meantime, checked point 2 from the to-do list above and I had it the wrong way around -
a few (9) muts are present in UG but missing in HC - but looking at IGV for each of these,
they look like legitimate calls, so keeping them in there

so that leaves the 'multi muts/indels' analysis (muts mostly checked, save for a few 
particularly confusing ones) and calculating the fraction of callable sites before
getting started with the remaining questions/analyses:

- does salt stress affect de novo mutation?
    - overall mut rate (SNM rate, indel rate) - can get this sooner (e.g. before 'multi' analysis)
        - need to do this for salt MA, MA, salt adaptation, and report
        - indel:SNM ratio
    - mut context
    - base spectrum
        - 6 mutation rates for MA, saltMA, adaptation
        - 32 triplet rates for MA, saltMA, adaptation 
            - could I start this for the previous files? is in the ms already?
- effect of strain on salt MA
    - compare de novo mut rate between salt and non-salt conditions


looking through the ms, currently we have -

- context for SNMs + indels in MA and adaptation lines
- distribution for SNMs + indels in MA and adaptation lines
    - incl. inter-mutation distance
- base spectrum for MA and adaptation
- eq. GC for MA and adaptation

so for now I can do -

- 32 triplet rates for MA and adaptation lines

first - let's find the files with the original mut describer output to
calculate this from, as well as the relevant VCFs - ideally we want to grab the
surrounding bases from the same line instead of from the reference in case they
differ from the ref strain

```
# annotations/output
SL1_haplotype_shared_salt_muts.annotated.txt 
SL1_haplotype_shared_salt_indels.annotated.txt 
SL2_haplotype_shared_salt_muts.annotated.txt 
SL2_haplotype_shared_salt_indels.annotated.txt 
unique_salt_muts.annotated.txt
unique_salt_indels.annotated.txt

# /research/projects/chlamydomonas/bgi_full_MA/mutation_calls
final.curated_muts.coord_sorted.txt.gz
```

starting a new analysis folder for this! `analysis/spectrum_context`



