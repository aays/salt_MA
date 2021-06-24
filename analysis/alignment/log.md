
## 18/5/2020

today:

- downloading all these sweet, sweet fastq files

command in private notion project page, for obvious reasons - needed to use
`wget --no-remove-listing` to get the actual file listings and then `wget -r` to
ensure recursive downloading

took 30 minutes but it's looks like we're good to go! 

checking checksums:

```python
import subprocess
from tqdm import tqdm
import glob
import os
import re

counter = 0
pass_counter = 0
dirnames = os.listdir()
for d in tqdm(dirnames):
    if os.path.isdir(d) and re.search('[A-Z]{2}[0-9]{2,4}_[05]', d):
        fastq_files = glob.glob(d + '/*fq.gz')
        checksum_file = d + '/MD5.txt'
        with open(checksum_file, 'r') as f:
            lines = [line.rstrip('\n').split('  ') for line in f] # two spaces
            checksums = {d + '/' + line[1]: line[0] for line in lines} # fname: checksum
        own_checksums = dict.fromkeys(fastq_files)
        for fname in fastq_files:
            proc = subprocess.Popen(['md5sum', fname],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            ch_split = out.decode('utf-8').rstrip('\n').split('  ')
            if checksums[fname] == ch_split[0]:
                pass_counter += 1
            counter += 1

print(counter, pass_counter)
# 82 82
```

finally, changing permissions for `data/fastq` to 544 just to be safe (dirs need
to be set to r-x to even be openable)

## 20/5/2020

today: 

- QC and trimming reads
- alignment with bwa mem

getting the reference fasta prepped:

```bash
mkdir -p data/references
mkdir -p data/bam
mkdir -p data/bam/sam
cd data/references/
ln -sv /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta* . 
```

making SAM files - going to make temp symlinks in `data/fastq/symlinks` to speed this up

```bash
cd data/fastq
mkdir -p symlinks
cd symlinks

for fname in $(find ../ -name '*fq.gz'); do
    ln -sv ${fname} symlinks/;
done

# writing to initial file
ls *fq.gz > fnames.txt
```

cleaning the symlink names:

```python
import subprocess 
import re

cmd = 'mv -v {} {}'
pattern = '([CDLS]{2}[0-9]{2,4}_[05])[A-Za-z0-9_\-]+(_[12]\.fq\.gz)'

with open('fnames.txt', 'r') as f:
    for line in f:
        line = line.rstrip()
        match = re.search(pattern, line)
        name_start, name_end = match.group(1), match.group(2)
        outname = name_start + name_end
        subprocess.call(cmd.format(line, outname), shell=True)
```

file with just prefixes:

```bash
# in data/fastq/symlinks
touch prefixes.txt
for read in 1 2; do
    for fname in *_${read}.fq.gz; do
        basename ${fname} _${read}.fq.gz >> prefixes.txt;
    done
done

sort prefixes.txt | uniq | sponge prefixes.txt
```

fastqc - this will take over two hours

```bash
mkdir -p data/fastq/fastqc
fastqc --kmers 7 --outdir data/fastq/fastqc \
data/fastq/symlinks/*fq.gz
```

## 21/5/2020

looking through the fastqc data:

```bash
cd data/fastq/fastqc
mkdir raw
mv -v *zip raw/
mv -v *html raw/
for fname in raw/*zip; do
    unzip ${fname};
done
```

some of these are failing the per base sequence quality check:

```bash
cd data/fastq/fastqc
grep -P 'FAIL\tPer base sequence quality' */summary.txt
```

prepping trimmomatic:

```bash
mkdir -p bin/
ln -sv /scratch/research/tmp_apps/Trimmomatic-0.36/trimmomatic-0.36.jar bin/trimmomatic
```

based on the sequencing report, we've got TruSeq2-PE adapters:

```bash
ln -sv /scratch/research/tmp_apps/Trimmomatic-0.36/adapters/TruSeq2-PE.fa bin/TruSeq2-PE.fa
```

trimming reads:

```bash
mkdir data/fastq_trim

# starting on just one
java -jar bin/trimmomatic PE -threads 4 -phred33 \
data/fastq/symlinks/CC1373_0_1.fq.gz \
data/fastq/symlinks/CC1373_0_2.fq.gz \
data/fastq_trim/CC1373_0_trim_1.fq.gz \
data/fastq_trim/CC1373_0_trim_unpaired_1.fq.gz \
data/fastq_trim/CC1373_0_trim_2.fq.gz \
data/fastq_trim/CC1373_0_trim_unpaired_2.fq.gz \
ILLUMINACLIP:bin/TruSeq2-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20
```

looks good - took 6 minutes with 4 threads

let's just do the CC files and then rerun fastqc on them to see whether that helps

```bash
time while read fname; do
    if [[ ${fname} =~ "CC" ]]; then
        time java -jar bin/trimmomatic PE -threads 4 -phred33 \
        data/fastq/symlinks/${fname}_1.fq.gz \
        data/fastq/symlinks/${fname}_2.fq.gz \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_unpaired_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz \
        data/fastq_trim/${fname}_trim_unpaired_2.fq.gz \
        ILLUMINACLIP:bin/TruSeq2-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20;
    fi;
    sleep 3;
done < data/fastq/symlinks/prefixes.txt
```

## 23/5/2020

so that just took 1.5 hours total

let's fastqc these:

```bash
mkdir data/fastq_trim/fastqc
mkdir data/fastq_trim/fastqc/raw

fastqc --kmers 7 --outdir data/fastq_trim/fastqc/raw \
data/fastq_trim/CC*fq.gz
```

## 25/5/2020

apparently fastqc crashed 5 minutes in on the first unpaired file,
though it worked fine on the paired trim files

wait - I should just be using the paired files here anyways (unpaired spans <
0.1% of the total read count) - though double check with Rob about this later

```bash
fastqc --kmers 7 --outdir data/fastq_trim/fastqc/raw \
data/fastq_trim/CC*_?_trim_?.fq.gz
```

## 27/5/2020

checking the fastq results:

```bash
cd data/fastq_trim/fastqc
for fname in raw/*zip; do
    unzip ${fname};
done

grep -P 'FAIL\tPer base sequence quality' */summary.txt
grep -P 'FAIL' */summary.txt 
```

the only failing tests are per base sequence content, which is expected
given the ridiculously high GC in the chlamy genome

now to trim the SL reads:

```bash
time while read fname; do
    if [[ ${fname} =~ "SL" ]]; then
        time java -jar bin/trimmomatic PE -threads 4 -phred33 \
        data/fastq/symlinks/${fname}_1.fq.gz \
        data/fastq/symlinks/${fname}_2.fq.gz \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_unpaired_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz \
        data/fastq_trim/${fname}_trim_unpaired_2.fq.gz \
        ILLUMINACLIP:bin/TruSeq2-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20
    fi;
    sleep 3;
done < data/fastq/symlinks/prefixes.txt
```

and the DL ones:

```bash
time while read fname; do
    if [[ ${fname} =~ "DL" ]]; then
        time java -jar bin/trimmomatic PE -threads 4 -phred33 \
        data/fastq/symlinks/${fname}_1.fq.gz \
        data/fastq/symlinks/${fname}_2.fq.gz \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_unpaired_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz \
        data/fastq_trim/${fname}_trim_unpaired_2.fq.gz \
        ILLUMINACLIP:bin/TruSeq2-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20
    fi;
    sleep 3;
done < data/fastq/symlinks/prefixes.txt
```

all done - the SL took 28 min, and DL 2 hours

now to fastq them:

```bash
fastqc --kmers 7 --outdir data/fastq_trim/fastqc/raw \
data/fastq_trim/DL*_?_trim_?.fq.gz

fastqc --kmers 7 --outdir data/fastq_trim/fastqc/raw \
data/fastq_trim/SL*_?_trim_?.fq.gz
```

I could glob both with `?L*_` but this way both can be run at once

while these are running, might as well get started with alignment for the CC strains:

```bash
# from project root
time while read fname; do
    ref=data/references/chlamy.5.3.w_organelles_mtMinus.fasta
    if [[ ${fname} =~ "CC" ]]; then
        time bwa mem -t 12 $ref \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz | \
        samtools sort -@4 -T ${fname}.sorting.tmp -O bam -o data/bam/${fname}.sorted.bam
    fi;
done < data/fastq/symlinks/prefixes.txt
```

## 28/5/2020

so the alignments finished just short of 3 hours

checking on the fastqc results for the SL and DL strains:

```bash
cd data/fastq_trim/fastqc
for fname in raw/?L*zip; do
    unzip ${fname};
done

grep -P 'FAIL\tPer base sequence quality' */summary.txt # no failing files
grep -P 'FAIL' */summary.txt # all for GC content + sequence content
```

now for alignment with these new trimmed reads - although bash regex is a bit
unfriendly, and I had to google around a fair bit to get this pattern right

```bash
# from project root
time while read fname; do
    ref=data/references/chlamy.5.3.w_organelles_mtMinus.fasta
    if [[ ${fname} = *@(DL|SL)*" ]]; then
        time bwa mem -t 12 $ref \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz | \
        samtools sort -@4 -T ${fname}.sorting.tmp -O bam -o data/bam/${fname}.sorted.bam
    fi;
done < data/fastq/symlinks/prefixes.txt
```

apparently samtools disappeared from `opt` overnight... I'll have to use the
one in `/scratch/research/tmp_apps` for now

```bash
cd bin/
ln -sv /scratch/research/tmp_apps/samtools-1.9/samtools .
```

## 29/5/2020

trying again:

```bash
# from project root
time while read fname; do
    ref=data/references/chlamy.5.3.w_organelles_mtMinus.fasta
    if [[ ${fname} = *@(DL|SL)* ]]; then
        time bwa mem -t 12 $ref \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz | \
        ./bin/samtools sort -@5 -T ${fname}.sorting.tmp -O bam -o data/bam/${fname}.sorted.bam
    fi;
done < data/fastq/symlinks/prefixes.txt
```

post bwa - fixing mates and then adding read groups:

```bash
time while read fname; do
    if [[ ${fname} =~ "CC" ]]; then
        time java -jar /opt/picard-tools/picard.jar FixMateInformation \
        I=data/bam/${fname}.sorted.bam \
        O=data/bam/${fname}.fixMate.bam \
        VALIDATION_STRINGENCY=LENIENT;
    fi;
done < data/fastq/symlinks/prefixes.txt
```

## 1/6/2020

fixing mate pair info for the remaining bams:

```bash
time while read fname; do
    if [[ ${fname} = *@(DL|SL)* ]]; then
        time java -jar /opt/picard-tools/picard.jar FixMateInformation \
        I=data/bam/${fname}.sorted.bam \
        O=data/bam/${fname}.fixMate.bam \
        VALIDATION_STRINGENCY=LENIENT;
    fi;
done < data/fastq/symlinks/prefixes.txt
```

## 2/6/2020

next up - adding read groups:

```bash
time while read fname; do
    time java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups \
    I=data/bam/${fname}.fixMate.bam \
    O=data/bam/${fname}.RG.bam \
    RGID=${fname} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=${fname} \
    VALIDATION_STRINGENCY=LENIENT;
    sleep 3;
done < data/fastq/symlinks/prefixes.txt

```

## 5/6/2020

finally - using MarkDuplicates:

```bash
mkdir -p data/bam_final
time while read fname; do
    time java -jar /opt/picard-tools/picard.jar MarkDuplicates \
    I=data/bam/${fname}.RG.bam \
    O=data/bam_final/${fname}.bam \
    METRICS_FILE=analysis/dup_metrics.txt
    sleep 3;
done < data/fastq/symlinks/prefixes.txt # took 3.5 hours
```

## 9/6/2020

creating bam index files:

```bash
cd data/bam_final
time for fname in *.bam; do
    samtools index ${fname}
done
```

reorganizing data folder to match Josianne's file structure:

```bash
rm -rfv data/bam
mv -v data/bam_final data/bam
mkdir -p data/alignments
mv -v data/bam/ data/alignments/bam
mv -v data/fastq data/alignments/fastq
mv -v data/fastq_trim data/alignments/fastq_trim
```

next up: starting out with genotyping

need to make VCFs with both UG and HC, both in diploid mode for consistency
with how it was done previously

I think it makes most sense to group by strain (CC, SL, DL), and within each
strain group, split by treatment (0 and 5)

will need to check the notebook commands and try to match them as closely
as possible 

## 13/6/2020

so the original UnifiedGenotyper and HaplotypeCaller commands don't seem to
be anywhere in the project folder - will just try to match the methods as closely
as possible then

need to keep the 0 and 5 treatments separate for each of
our sample groups

will start with UnifiedGenotyper for the SL lines - from the methods:

- set ploidy to diploid
- `heterozygosity 0.02`
- `indel_heterozygosity 0.002`
- `--output_mode EMIT_ALL_SITES`
- `--genotype_likelihoods_model BOTH`

going to start with just chromosome 1 for debugging purposes

```bash
mkdir -p data/alignments/genotyping
mkdir -p data/alignments/genotyping/UG
mkdir -p data/alignments/genotyping/HC

time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL26_5.bam \
-I data/alignments/bam/SL27_5.bam \
-I data/alignments/bam/SL29_5.bam \
-L chromosome_1 \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/SL_chr1.vcf
```

looks good, and was done in 11 min - let's queue this up for the full genome,
across both 0 and 5

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL26_5.bam \
-I data/alignments/bam/SL27_5.bam \
-I data/alignments/bam/SL29_5.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/SL_5.vcf

sleep 3

time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL27_0.bam \ # no SL26_0
-I data/alignments/bam/SL29_0.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/SL_0.vcf
```

haplotypecaller for these samples - going to be using the same settings:

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL26_5.bam \
-I data/alignments/bam/SL27_5.bam \
-I data/alignments/bam/SL29_5.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/SL_5.vcf

sleep 3

time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL27_0.bam \
-I data/alignments/bam/SL29_0.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/SL_0.vcf
```

## 14/6/2020

running HC commands today

## 16/6/2020

that script (one UG, two HC) took a full day to run

queuing up UG commands for DL lines next:

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL40_5.bam \
-I data/alignments/bam/DL41_5.bam \
-I data/alignments/bam/DL46_5.bam \
-I data/alignments/bam/DL51_5.bam \
-I data/alignments/bam/DL53_5.bam \
-I data/alignments/bam/DL55_5.bam \
-I data/alignments/bam/DL57_5.bam \
-I data/alignments/bam/DL58_5.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/DL_5.vcf

sleep 3

time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL40_0.bam \
-I data/alignments/bam/DL41_0.bam \
-I data/alignments/bam/DL46_0.bam \
-I data/alignments/bam/DL51_0.bam \
-I data/alignments/bam/DL53_0.bam \
-I data/alignments/bam/DL55_0.bam \
-I data/alignments/bam/DL57_0.bam \
-I data/alignments/bam/DL58_0.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/DL_0.vcf
```

done in 10 hours - I guess UG is a bit faster than HC

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL40_5.bam \
-I data/alignments/bam/DL41_5.bam \
-I data/alignments/bam/DL46_5.bam \
-I data/alignments/bam/DL51_5.bam \
-I data/alignments/bam/DL53_5.bam \
-I data/alignments/bam/DL55_5.bam \
-I data/alignments/bam/DL57_5.bam \
-I data/alignments/bam/DL58_5.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/DL_5.vcf

sleep 3

time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL40_0.bam \
-I data/alignments/bam/DL41_0.bam \
-I data/alignments/bam/DL46_0.bam \
-I data/alignments/bam/DL51_0.bam \
-I data/alignments/bam/DL53_0.bam \
-I data/alignments/bam/DL55_0.bam \
-I data/alignments/bam/DL57_0.bam \
-I data/alignments/bam/DL58_0.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/DL_0.vcf
```

## 19/6/2020

I've been doing this incorrectly - it makes a lot more sense to
be calling pairs of the same sample (eg both treatments) - once
this is done, will need to add in other 'reference' strains
(eg wild type, earlier SLs, earlier DLs) using GATK

going to zip the existing VCFs and put them in new directories called `grouped`

also - HC `--output-mode` should be `EMIT_ALL_CONFIDENT_SITES`, not `EMIT_ALL_SITES`
(which is a UG command)

```bash
mkdir -p data/alignments/genotyping/HC/grouped
mkdir -p data/alignments/genotyping/UG/grouped
```

getting started with UG - can get the exact command that was run by
grepping for `GATKCommandLine`:

```bash
# in parent project dir
zgrep -m 1 'GATKCommandLine' alignments/genotyping/unified_genotyper/CC2935.UG.vcf.gz
```

modding the output in vim with `s/,/\r/g` followed by `s/ /\r/g` on the final line

going to delete anything I'm absolutely sure won't be needed

```
ID=UnifiedGenotyper
Version=3.5-0-g36282e4
Date="Thu Aug 04 15:05:17 EDT 2016"
Epoch=1470337517036
CommandLineOptions="analysis_type=UnifiedGenotyper
input_file=[./CC2935//CC2935.realigned.bam]
downsampling_type=BY_SAMPLE
downsample_to_fraction=null
downsample_to_coverage=250
validation_strictness=SILENT
sites_only=false
genotype_likelihoods_model=BOTH
min_base_quality_score=1
max_deletion_fraction=0.05
min_indel_count_for_genotyping=5
min_indel_fraction_per_sample=0.25
indelGapContinuationPenalty=10
indelGapOpenPenalty=45
indelHaplotypeSize=80
indelDebug=false
min_quality_score=1
max_quality_score=40
site_quality_prior=20
min_power_threshold_for_calling=0.95
heterozygosity=0.02
indel_heterozygosity=0.002
standard_min_confidence_threshold_for_calling=0.0
standard_min_confidence_threshold_for_emitting=0.0
max_alternate_alleles=6
sample_ploidy=1 # this should be 2 for these samples
genotyping_mode=DISCOVERY
output_mode=EMIT_ALL_CONFIDENT_SITES
allSitePLs=false
out=/scratch/research/projects/chlamydomonas/salt_lines/alignments/genotyping/unified_genotyper/CC2935.UG.vcf
```

going to create a text file called `sample_names.txt` similar to `prefixes.txt` from earlier

```python
fname = 'data/alignments/fastq/symlinks/prefixes.txt'
outname = 'data/alignments/fastq/symlinks/samples.txt'
with open(fname, 'r') as f:
    lines = [l.split('_') for l in f.readlines()]
samples = [l[0] for l in lines]
samples = list(set(samples)) # get unique sample names
samples.sort()
with open(outname, 'w') as f:
    for sample in samples:
        f.write(sample + '\n')
```

and now to use this file:

```bash
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T UnifiedGenotyper \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        -I data/alignments/bam/${sample}_0.bam \
        -I data/alignments/bam/${sample}_5.bam \
        -glm BOTH \
        -ploidy 2 \
        --output_mode EMIT_ALL_SITES \ 
        --heterozygosity 0.02 \
        --indel_heterozygosity 0.002 \
        -o data/alignments/genotyping/UG/${sample}_samples.vcf;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

(will have to handle SL26 differently since there's only one sample there)

## 21/6/2020

whole thing took just under a day - 3 hrs/run

queuing up the DL files:

```bash
time while read sample; do
    if [[ ${sample} =~ "DL" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T UnifiedGenotyper \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        -I data/alignments/bam/${sample}_0.bam \
        -I data/alignments/bam/${sample}_5.bam \
        -glm BOTH \
        -ploidy 2 \
        --output_mode EMIT_ALL_SITES \ 
        --heterozygosity 0.02 \
        --indel_heterozygosity 0.002 \
        -o data/alignments/genotyping/UG/${sample}_samples.vcf;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

while this runs, let's test `CombineVariants`:

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T CombineVariants \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--variant data/alignments/genotyping/UG/CC1373_samples.vcf \
--variant data/alignments/genotyping/UG/CC1952_samples.vcf \
-L chromosome_1 \
-o test_combine.vcf \
-genotypeMergeOptions UNSORTED
```

this works well - the merge setting could also be UNIQUIFY if there are sample
name issues, but there aren't any right now, so UNSORTED works fine for now

trying this with one of the older VCFs 

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
-T CombineVariants \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
--variant ../alignments/genotyping/unified_genotyper/UG_diploid/salt/salt_aligned_salt_UG.diploid.vcf.gz \
--variant data/alignments/genotyping/UG/CC1373_samples.vcf \
--variant data/alignments/genotyping/UG/CC1952_samples.vcf \
-L chromosome_1 \
-o test_combine_salt.vcf \
-genotypeMergeOptions UNSORTED
```

works great! will be revisiting this once all the pairs have been called. 

## 22/6/2020

today: UG for SL lines, and getting HC queued up

given that there are 2 (well, 2 and a half) pairs here, I should just type these up manually

```bash
time for i in 27 29; do
    time java -jar ./bin/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/SL${i}_0.bam \
    -I data/alignments/bam/SL${i}_5.bam \
    -glm BOTH \
    -ploidy 2 \
    --output_mode EMIT_ALL_SITES \ 
    --heterozygosity 0.02 \
    --indel_heterozygosity 0.002 \
    -o data/alignments/genotyping/UG/SL${i}_samples.vcf;
done

time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL26_5.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \ 
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/SL26_samples.vcf
```

## 23/6/2020

today: 

- bgzip yesterday's VCFs
- queue up HC commands for pairs
- combine UG calls with earlier files

bgzipping:

```bash
time for fname in data/alignments/genotyping/UG/*vcf; do
    time bgzip ${fname}
    echo "done ${fname}"
done
```

getting HC command from previous runs:

```bash
# in parent project dir
zgrep -m 1 'GATKCommandLine' alignments/genotyping/HC_diploid/all_SL_HC/allSL.diploid.HC.variants.vcf.gz
```

output - copying over in tmux and then using `s/ /\r/g` 

```
##GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller,CommandLineOptions="analysis_type=HaplotypeCaller
input_file=[S21/RG.bam, S23/RG.bam, S24/RG.bam, S25/RG.bam,
S26/RG.bam, S27/RG.bam, S29/RG.bam, S30/RG.bam, S31/RG.bam, S33/RG.bam]
showFullBamList=false
read_buffer_size=null
phone_home=AWS
gatk_key=null
intervals=[/scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/interval_list/10.interval_list]
excludeIntervals=null
interval_set_rule=UNION
interval_merging=ALL
interval_padding=0
reference_sequence=/scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta
downsampling_type=BY_SAMPLE
downsample_to_fraction=null
downsample_to_coverage=500
fix_misencoded_quality_scores=false
allow_potentially_misencoded_quality_scores=false
useOriginalQualities=false
defaultBaseQualities=-1
performanceLog=null
disable_indel_quals=false
emit_original_quals=false
validation_strictness=SILENT
sites_only=false
never_trim_vcf_format_field=false
bcf=false
bam_compression=null
simplifyBAM=false
disable_bam_indexing=false
generate_md5=false
pedigreeValidationType=STRICT
allow_intervals_with_unindexed_bam=false
generateShadowBCF=false
variant_index_type=DYNAMIC_SEEK
variant_index_parameter=-1
reference_window_stop=0
logging_level=INFO
log_to_file=null
help=false
version=false
out=/scratch/research/projects/chlamydomonas/salt_lines/alignments/./all_SL_HC/allSL.10.diploid.HC.variants.vcf
likelihoodCalculationEngine=PairHMM
heterogeneousKmerSizeResolution=COMBO_MIN
dontTrimActiveRegions=false
maxDiscARExtension=25
maxGGAARExtension=300
paddingAroundIndels=150
paddingAroundSNPs=20
emitRefConfidence=NONE
bamWriterType=CALLED_HAPLOTYPES
disableOptimizations=false
annotateNDA=false
heterozygosity=0.02
indel_heterozygosity=0.002
standard_min_confidence_threshold_for_calling=0.0
standard_min_confidence_threshold_for_emitting=0.0
max_alternate_alleles=6
input_prior=[]
sample_ploidy=2
genotyping_mode=DISCOVERY
contamination_fraction_to_filter=0.0
contamination_fraction_per_sample_file=null
exactcallslog=null
output_mode=EMIT_VARIANTS_ONLY ### variants only! 
allSitePLs=false
gcpHMM=10
pair_hmm_implementation=VECTOR_LOGLESS_CACHING
pair_hmm_sub_implementation=ENABLE_ALL
always_load_vector_logless_PairHMM_lib=false
phredScaledGlobalReadMismappingRate=45
GVCFGQBands=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
70, 80, 90, 99]
errorCorrectReads=false
pcr_indel_model=CONSERVATIVE
maxReadsInRegionPerSample=10000
minReadsPerAlignmentStart=10
maxProbPropagationDistance=50
activeProbabilityThreshold=0.002
min_mapping_quality_score=20
filter_reads_with_N_cigar=false
filter_mismatching_base_and_quals=false
filter_bases_not_stored=false"
Date="Sat Oct 01 07:00:08 EDT 2016",Epoch=1475319608804,Version=3.5-0-g36282e4>
```

getting started with CC samples:

```bash
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        -I data/alignments/bam/${sample}_0.bam \
        -I data/alignments/bam/${sample}_5.bam \
        -ploidy 2 \
        --output_mode EMIT_VARIANTS_ONLY \ 
        --heterozygosity 0.02 \
        --indel_heterozygosity 0.002 \
        -o data/alignments/genotyping/HC/${sample}_samples.vcf;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

while this runs, can get started on using CombineVariants. the SL
lines make sense as first pass candidates - these should include the
dark ancestors, CC2935 (diploid), the previous SL lines, and finally the pairs

will need to tabix and bgzip the files in `data/alignments/genotyping/UG` for compatibility
with CombineVariants

```bash
mkdir -p data/alignments/genotyping/UG/pairs
mkdir -p data/alignments/genotyping/UG/combined
mv -v data/alignments/genotyping/UG/*vcf* data/alignments/genotyping/UG/pairs

# doing all sequentially - no reason to split these up
time while read sample; do
    time java -jar ./bin/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/dark_ancestors/dark_aligned_dark_UG_diploid.vcf.gz \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/salt/salt_aligned_salt_UG.diploid.vcf.gz \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/CC2935_diploid.UG.vcf.gz \
    --variant data/alignments/genotyping/UG/pairs/${sample}_samples.vcf.gz \
    -o data/alignments/genotyping/UG/combined/${sample}_combined.vcf \
    -genotypeMergeOptions UNSORTED;
done < data/alignments/fastq/symlinks/samples.txt
```

## 25/6/2020

today:

- cancel HC and modify:
    - emit all sites
    - use interval files to parallelize over chromosomes
    - also include scaffolds + organelles as intervals in separate file
- do first pass mutation calling check in UG files using cyvcf2
    - diff b/w 0 and 5 lines, GQ > 30, homozygous calls

combining UG files is done, took a full day (~90 min per VCF)

new HC code - going to make a temp script for this (`HC_temp.sh`) - will
start with CC samples

```bash
chr=$1
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        -I data/alignments/bam/${sample}_0.bam \
        -I data/alignments/bam/${sample}_5.bam \
        -L chromosome_${chr}
        -ploidy 2 \
        --output_mode EMIT_ALL_ACTIVE_SITES \ 
        --heterozygosity 0.02 \
        --indel_heterozygosity 0.002 \
        -o data/alignments/genotyping/HC/${sample}_samples.vcf;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

running this in parallel - just the first five chrs, as a test:

```bash
parallel -j 5 -i sh -c 'time bash HC_temp.sh {}' -- {1..5}
```

## 26/6/2020

queuing up remaining chrs:

```bash
parallel -j 2 -i sh -c 'time bash HC_temp.sh {}' -- 6 7
```

while this is happening, going to:

- bgzip and tabix the UG combined files 
- do a first pass check for muts that differ between 0 and 5 in the UG files

```bash
time for fname in data/alignments/genotyping/UG/combined/*vcf; do
    time bgzip ${fname}
    tabix -p vcf ${fname}.gz
    echo "done ${fname}"
done
```

meanwhile, working off the pairs file to look for some variants:

```python
>>> from cyvcf2 import VCF
>>> fname = 'CC1373_samples.vcf.gz'
>>> v = VCF(fname)
>>> snps = []
>>> v = VCF(fname)
>>> snps = []
>>> for rec in tqdm(v):
  2     if len(rec.ALT) > 0 and rec.num_het == 0 and not rec.is_indel:
  3         if rec.genotypes[0] != rec.genotypes[1]:
  4             if not any([-1 in call for call in rec.genotypes]):
  5                 if all(rec.gt_quals > 30):
  6                     snps.append(rec)
  7                     print('found', rec)
  8     if len(snps) == 10:
  9         break
>>> snps
[Variant(chromosome_2:6902447 G/A), Variant(chromosome_7:5244681 C/A), Variant(chromosome_12:6997391 G/T), 
Variant(chromosome_17:4899631 G/C), Variant(mtDNA:573 T/C), Variant(mtDNA:3693 T/G)]
```

oh shit! those are some solid first pass SNPs - let's do the same for CC1952:

```bash
>>> fname = 'CC1952_samples.vcf.gz'
>>> snps_1952 = []
>>> v = VCF(fname)
  2 for rec in tqdm(v):
  3     if len(rec.ALT) > 0 and rec.num_het == 0 and not rec.is_indel:
  4         if rec.genotypes[0] != rec.genotypes[1]:
  5             if not any([-1 in call for call in rec.genotypes]):
  6                 if all(rec.gt_quals > 30):
  7                     snps_1952.append(rec)
  8                     print('found', rec)
  9     if len(snps_1952) == 10:
 10         break
found chromosome_12    343437  .       C       T       418.61  .       AC=2;AF=0.5;AN=4;BaseQRankSum=0.423;DP=28;Dels=0;ExcessHet=0.7918;FS=3.556;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=60;MQ0=0;MQRankSum=0.047;QD=28.1;ReadPosRankSum=-1.035;SOR=0.531     GT:AD:DP:GQ:PL  1/1:0,11:11:33:440,33,0 0/0:17,0:17:51:0,51,670
```

just the one this time - moving on to the next (2342)

```bash
>>> snps_2342 = [] # I should have used a dict, lol
>>> fname = 'CC2342_samples.vcf.gz'
>>> v = VCF(fname)
>>> v = VCF(fname)
  2 for rec in tqdm(v):
  3     if len(rec.ALT) > 0 and rec.num_het == 0 and not rec.is_indel:
  4         if rec.genotypes[0] != rec.genotypes[1]:
  5             if not any([-1 in call for call in rec.genotypes]):
  6                 if all(rec.gt_quals > 30):
  7                     snps_2342.append(rec)
  8                     print('found', rec)
  9     if len(snps_2342) == 10:
 10         break
found chromosome_1       509598  .       C       T       568.61  .       AC=2;AF=0.5;AN=4;BaseQRankSum=0.37;DP=27;Dels=0;ExcessHet=0.7918;FS=1.657;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=60;MQ0=0;MQRankSum=0.025;QD=25.19;ReadPosRankSum=-0.666;SOR=0.269     GT:AD:DP:GQ:PL  0/0:11,0:11:33:0,33,440 1/1:0,16:16:45:590,45,0

found chromosome_7     667812  .       C       T       532.61  .       AC=2;AF=0.5;AN=4;BaseQRankSum=-0.487;DP=28;Dels=0;ExcessHet=0.7918;FS=1.606;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=59.31;MQ0=0;MQRankSum=-1.416;QD=33.29;ReadPosRankSum=0.627;SOR=1.127        GT:AD:DP:GQ:PL  0/0:12,0:12:36:0,36,480 1/1:0,16:16:42:554,42,0

found chromosome_11    83491   .       A       C       703.6   .       AC=2;AF=0.5;AN=4;BaseQRankSum=-0.48;DP=35;Dels=0;ExcessHet=0.7918;FS=8.193;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=60;MQ0=0;MQRankSum=-0.414;QD=26.32;ReadPosRankSum=-0.381;SOR=0.094   GT:AD:DP:GQ:PL  0/0:16,0:16:48:0,48,640 1/1:0,19:19:54:725,54,0

found chromosome_11    1661849 .       G       A       408.61  .       AC=2;AF=0.5;AN=4;BaseQRankSum=0.406;DP=24;Dels=0;ExcessHet=0.7918;FS=1.569;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=59.84;MQ0=0;MQRankSum=0.753;QD=27.52;ReadPosRankSum=1.043;SOR=1.085  GT:AD:DP:GQ:PL  0/0:13,0:13:39:0,39,520 1/1:0,11:11:33:430,33,0

found chromosome_13    1323388 .       G       T       448.6   .       AC=2;AF=0.5;AN=4;BaseQRankSum=1.546;DP=35;Dels=0;ExcessHet=0.7918;FS=0;HaplotypeScore=0;MLEAC=2;MLEAF=0.5;MQ=60;MQ0=0;MQRankSum=0.99;QD=33.76;ReadPosRankSum=1.477;SOR=0.922  GT:AD:DP:GQ:PL  1/1:0,12:12:36:470,36,0 0/0:23,0:23:57:0,57,779
```

## 28/6/2020

queuing up next HC step:

```bash
parallel -j 3 -i sh -c 'time bash HC_temp.sh {}' -- 8 9 10
```

and now that bgzipping + tabix is done, going to make variant only versions
of the combined VCFs with bcftools:

```bash
mkdir -p data/alignments/genotyping/UG/variants
cd bin/
ln -sv ~/apps/bcftools/bcftools

# quick test
time ./bin/bcftools filter -i 'TYPE!="ref"' \
data/alignments/genotyping/UG/combined/CC1373_combined.vcf.gz \
test.vcf # took ~8 min

# looks good - doing this over all VCFs
time for fname in data/alignments/genotyping/UG/combined/*.vcf.gz; do
    echo "currently on ${fname}"
    base=$(basename ${fname} .vcf.gz)
    time ./bin/bcftools filter -i 'TYPE!="ref"' ${fname} > \
    data/alignments/genotyping/UG/variants/${base}.variants.vcf;
    echo "done ${fname}"
done
```

## 30/6/2020

next HC step:

```bash
parallel -j 3 -i sh -c 'time bash HC_temp.sh {}' -- 11 12 13
```

### 1/7/2020

final HC step (for CC strains...)

```bash
parallel -j 4 -i sh -c 'time bash HC_temp.sh {}' -- {14..17}
```

## 2/7/2020

CC strains done! will have to combine chromosomal VCFs into full ones,
bgzip + tabix those, as well as make combined versions with prev strains
after this (although the prev strains for these will be different - doesn't
make sense to stick older SL lines with CCs)

wait - I did do that with the UG CC strains - might want to revisit those
and see whether I need to make different combined VCFs

now to do this with DL strains (swapping out CC in the regex match for DL)

```bash
parallel -j 4 -i sh -c 'time bash HC_temp.sh {}' -- {1..4}
```

combining all the CC haplotypecaller chromosomal VCFs:

```bash
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        echo "${sample}"
        grep '^#' data/alignments/genotyping/HC/${sample}_samples_chr1.vcf > \
            data/alignments/genotyping/HC/${sample}_samples.vcf
        for i in {1..17}; do
            grep -v '^#' data/alignments/genotyping/HC/${sample}_samples_chr${i}.vcf >> \
            data/alignments/genotyping/HC/${sample}_samples.vcf
            echo "added chromosome ${i} for ${sample}";
        done
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

manually checked that samples look alright - clearing out files:

```bash
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        rm -v data/alignments/genotyping/HC/${sample}_samples_chr*.vcf*
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

bgzipping and tabixing full files:

```bash
time for fname in data/alignments/genotyping/HC/CC*; do
    bgzip ${fname};
    tabix -p vcf ${fname}.gz;
done
```

now to use CombineVariants to add prev variants:

```bash
mkdir -p data/alignments/genotyping/HC/pairs
mkdir -p data/alignments/genotyping/HC/combined
# mv -v data/alignments/genotyping/HC/CC* data/alignments/genotyping/HC/pairs

time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
        --variant data/alignments/genotyping/HC/pairs/${sample}_samples.vcf.gz \
        -o data/alignments/genotyping/HC/combined/${sample}_combined.vcf \
        -genotypeMergeOptions UNSORTED;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

## 3/7/2020

doing organelles (cpDNA, mtDNA) as well - defining intervals in
`organelles.intervals` and modding `HC_temp.sh` as needed

shoot - I think I forgot to do these for the CC strains - will
queue that up next

after organelles are ready, DL strains are done - 
combining, bgzipping, and tabixing:

```bash
time while read sample; do
    if [[ ${sample} =~ "DL" ]]; then
        echo "${sample}"
        grep '^#' data/alignments/genotyping/HC/${sample}_samples_chr1.vcf > \
            data/alignments/genotyping/HC/${sample}_samples.vcf
        for i in {1..17}; do
            grep -v '^#' data/alignments/genotyping/HC/${sample}_samples_chr${i}.vcf >> \
            data/alignments/genotyping/HC/${sample}_samples.vcf
            echo "added chromosome ${i} for ${sample}";
        done
        grep -v '^#' data/alignments/genotyping/HC/${sample}_samples_organelles.vcf >> \
            data/alignments/genotyping/HC/${sample}_samples.vcf
    fi;
done < data/alignments/fastq/symlinks/samples.txt

# clearing chromosomal vcfs
time while read sample; do
    if [[ ${sample} =~ "DL" ]]; then
        rm -v data/alignments/genotyping/HC/${sample}_samples_chr*.vcf*
        rm -v data/alignments/genotyping/HC/${sample}_samples_organelles.vcf*
    fi;
done < data/alignments/fastq/symlinks/samples.txt

# bgzip and tabix
time for fname in data/alignments/genotyping/HC/DL*; do
    echo "starting ${fname}"
    time bgzip ${fname};
    tabix -p vcf ${fname}.gz;
    echo "done ${fname}"
done
```

## 5/7/2020

doing organelles for CC strains:

```bash
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        -I data/alignments/bam/${sample}_0.bam \
        -I data/alignments/bam/${sample}_5.bam \
        --intervals organelles.intervals \
        -ploidy 2 \
        --output_mode EMIT_ALL_SITES \
        --heterozygosity 0.02 \
        --indel_heterozygosity 0.002 \
        -o data/alignments/genotyping/HC/${sample}_samples_organelles.vcf;
    fi;
done < data/alignments/fastq/symlinks/samples.txt

date
```

while this runs, have to use CombineVariants to create combined versions
of the DL files

```bash
time while read sample; do
    if [[ ${sample} =~ "DL" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
        --variant ../alignments/genotyping/HC_diploid/all_dark_HC/all_dark_diploid.HC.variants.vcf.gz \
        --variant data/alignments/genotyping/HC/pairs/${sample}_samples.vcf.gz \
        -o data/alignments/genotyping/HC/combined/${sample}_combined.vcf \
        -genotypeMergeOptions UNSORTED;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

once organelles for CC strains are done, will need to:
- append those to paired VCFs
- do CombineVariants (on *just* the organelle files)
- add the CombineVariants output to the existing combined files
    - make sure that the samples are in the correct order before merging

## 6/7/2020

need to unzip paired CC VCFs before adding organelles, and then bgzip again

```bash
rm -v data/alignments/genotyping/HC/pairs/CC*tbi # clear index files
for fname in data/alignments/genotyping/HC/pairs/CC*; do
    echo ${fname}
    time bgzip -d ${fname}
done

# adding organelles
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        grep -v '^#' data/alignments/genotyping/HC/${sample}_samples_organelles.vcf >> \
        data/alignments/genotyping/HC/pairs/${sample}_samples.vcf
    fi;
done < data/alignments/fastq/symlinks/samples.txt

# rezipping and tabixing
for fname in data/alignments/genotyping/HC/pairs/CC*; do
    echo ${fname}
    time bgzip ${fname}
    tabix -p vcf ${fname}.gz
done

# combine variants for temp organelle files - to add to combined files
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --intervals organelles.intervals \
        --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
        --variant data/alignments/genotyping/HC/${sample}_samples_organelles.vcf \
        -o data/alignments/genotyping/HC/${sample}_combined_organelles.vcf \
        -genotypeMergeOptions UNSORTED;
    fi;
done < data/alignments/fastq/symlinks/samples.txt

# adding combined organelle files to combined full files
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        grep -v '^#' data/alignments/genotyping/HC/${sample}_combined_organelles.vcf >> \
        data/alignments/genotyping/HC/combined/${sample}_combined.vcf
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

...wait I forgot the scaffolds and mtMinus

*headdesk*

remaining to-do:

- HC for SL lines (chromosomes, scaffolds, organelles)
- combined files for SL lines
- HC for scaffolds in CC and DL lines
- update combined HC files for CC and DL lines with scaffolds

cleaning up organelle files now that they've been combined:

```bash
rm -v data/alignments/genotyping/HC/*organelles*
```

creating scaffold intervals file:

```bash
# scaffolds range from 18 to 54
touch scaffolds.intervals
for i in {18..54}; do
    echo "scaffold_${i}" >> scaffolds.intervals
done
echo "mtMinus" >> scaffolds.intervals
```

gonna try a different parallel strategy here - will parallelize samples over
all of these missing regions

```bash
# HC_scaffolds.sh
sample=$1

time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/${sample}_0.bam \
-I data/alignments/bam/${sample}_5.bam \
--intervals scaffolds.intervals \
-ploidy 2
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/${sample}_samples_scaffolds.vcf

date
```

test run:

```bash
bash HC_scaffolds.sh CC1373
```

looks good, was done in 16 min

next up:

```bash
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- CC1952 CC2342
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- CC2344 CC2931
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- CC2935 CC2937
```

all done - adding these to the pairs files:

```bash
rm -v data/alignments/genotyping/HC/pairs/CC*tbi # clear index files
for fname in data/alignments/genotyping/HC/pairs/CC*; do
    echo ${fname}
    time bgzip -d ${fname}
done

# adding scaffolds + mtMinus
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        grep -v '^#' data/alignments/genotyping/HC/${sample}_samples_scaffolds.vcf >> \
        data/alignments/genotyping/HC/pairs/${sample}_samples.vcf
    fi;
done < data/alignments/fastq/symlinks/samples.txt

# rezipping and tabixing
time for fname in data/alignments/genotyping/HC/pairs/CC*; do
    echo ${fname}
    time bgzip ${fname}
    tabix -p vcf ${fname}.gz
done

# combine variants for temp scaffold files - to add to combined files
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --intervals scaffolds.intervals \
        --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
        --variant data/alignments/genotyping/HC/${sample}_samples_scaffolds.vcf \
        -o data/alignments/genotyping/HC/${sample}_combined_scaffolds.vcf \
        -genotypeMergeOptions UNSORTED;
    fi;
done < data/alignments/fastq/symlinks/samples.txt

# adding combined scaffold files to combined full files
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        grep -v '^#' data/alignments/genotyping/HC/${sample}_combined_scaffolds.vcf >> \
        data/alignments/genotyping/HC/combined/${sample}_combined.vcf
    fi;
done < data/alignments/fastq/symlinks/samples.txt

```

## 7/7/2020

remaining to do:

- do scaffold runs for DL files and regenerate combined files
- bgzip and tabix DL and CC combined files
- HC for SL lines (chromosomes, scaffolds, organelles)
- create combined HC files for SL lines

for DL scaffold runs - just rerunning the code from just above but with
the regex matching DL instead

```bash
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- DL40 DL46
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- DL51 DL53
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- DL55 DL57
parallel -j 2 -i sh -c 'bash HC_scaffolds.sh {}' -- DL41 DL58
```

once this is done, rerun yesterday's code to clear existing combined files,
update pairs with scaffolds and mtMinus

might need to regen combined files from scratch in full - also do this
for just one of the CC files to make sure that those were also concatenated
correctly (bc I'm no longer sure) 

shit - the contigs need to be reordered - this might be an issue for the CC files as well

I can try fixing this with grep:

```bash
grep '^#' DL40_samples.vcf > DL40_samples_sorted.vcf
grep -v '^#' DL40_samples.vcf | grep '^chromosome' >> DL40_samples_sorted.vcf
grep -v '^#' DL40_samples.vcf | grep '^scaffold' >> DL40_samples_sorted.vcf
grep -v '^#' DL40_samples.vcf | grep '^cpDNA' >> DL40_samples_sorted.vcf
grep -v '^#' DL40_samples.vcf | grep '^mt[MD]' >> DL40_samples_sorted.vcf
```

wait - this actually doubled up the scaffold variants, and now the sorting
above actually made it *harder* to fix that

here's the fix:

```bash
sample=$1
grep '^#' ${sample}_samples.vcf > ${sample}_samples_fixed.vcf
grep '^chromosome' ${sample}_samples.vcf >> ${sample}_samples_fixed.vcf
grep '^scaffold' ${sample}_samples.vcf | sort -k1,1 -k2,2n | uniq >> ${sample}_samples_fixed.vcf
grep -v '^#' ${sample}_samples.vcf | grep '^[cm][pt]D' >> ${sample}_samples_fixed.vcf
grep '^mtMinus' ${sample}_samples.vcf | sort -k1,1 -k2,2n | uniq >> ${sample}_samples_fixed.vcf
```

meanwhile, the CC files (combined and pairs) are out of order, and can
be fixed with the grep code above:

```bash
sample=$1
grep '^#' ${sample}_samples.vcf > ${sample}_samples_sorted.vcf
grep -v '^#' ${sample}_samples.vcf | grep '^chromosome' >> ${sample}_samples_sorted.vcf
grep -v '^#' ${sample}_samples.vcf | grep '^scaffold' >> ${sample}_samples_sorted.vcf
grep -v '^#' ${sample}_samples.vcf | grep '^cpDNA' >> ${sample}_samples_sorted.vcf
grep -v '^#' ${sample}_samples.vcf | grep '^mt[MD]' >> ${sample}_samples_sorted.vcf
```

and now, just going to rerun CombineVariants in full:

```bash
# not specifying intervals
time while read sample; do
    if [[ ${sample} =~ "CC" ]]; then
        time java -jar ./bin/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
        --variant data/alignments/genotyping/HC/${sample}_samples_scaffolds.vcf \
        -o data/alignments/genotyping/HC/${sample}_combined_scaffolds.vcf \
        -genotypeMergeOptions UNSORTED;
    fi;
done < data/alignments/fastq/symlinks/samples.txt
```

## 8/7/2020

today: 
- HC for SL files
- CombineVariants for SL files
- bgzip and tabix pairs and combined vcfs

SL HC:

```bash
time for i in 27 29; do
    time java -jar ./bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/SL${i}_0.bam \
    -I data/alignments/bam/SL${i}_5.bam \
    -ploidy 2 \
    --output_mode EMIT_ALL_SITES \
    --heterozygosity 0.02 \
    --indel_heterozygosity 0.002 \
    -o data/alignments/genotyping/HC/pairs/SL${i}_samples.vcf;
done

time java -jar ./bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/SL26_5.bam \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/HC/pairs/SL26_samples.vcf
```

bgzipping and tabixing the pairs + combined vcfs:

```bash
# combined files
time for fname in data/alignments/genotyping/HC/combined/*combined.vcf; do
    echo "${fname}"
    bgzip ${fname}
    tabix -p vcf ${fname}.gz
    echo "done ${fname}"
done

# pairs
time for fname in data/alignments/genotyping/HC/pairs/*samples.vcf; do
    echo "${fname}"
    bgzip ${fname}
    tabix -p vcf ${fname}.gz
    echo "done ${fname}"
done
```


## 9/7/2020

- bgzip and tabix the SL files (once SL26 is completed)
- create combined versions of the SL files
- get started on script for creating filtered vcfs containing only potential mutations

getting started with filtering:

```bash
mkdir -p analysis/filtering
touch analysis/filtering/log.md
```

## 10/7/2020

- bgzip and tabix the SL files (once SL26 is completed)
- create combined versions of the SL files

creating combined versions:


```bash
for i in 26 27 29; do
    time java -jar ./bin/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
    --variant ../alignments/genotyping/HC_diploid/all_dark_HC/all_dark_diploid.HC.variants.vcf.gz \
    --variant ../alignments/genotyping/HC_diploid/all_SL_HC/allSL.diploid.HC.variants.vcf.gz \
    --variant data/alignments/genotyping/HC/pairs/SL${i}_samples.vcf.gz \
    -o data/alignments/genotyping/HC/combined/SL${i}_combined.vcf \
    -genotypeMergeOptions UNSORTED;
done
```

## 12/7/2020

oops - I missed a VCF with wild type strains when creating the combined UG files:

```bash
mv -v data/alignments/genotyping/UG/combined data/alignments/genotyping/UG/combined_old
mkdir -p data/alignments/genotyping/UG/combined

time while read sample; do
    time java -jar ./bin/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/wt_ancestors/wt_diploid.UG.vcf.gz \
    --variant data/alignments/genotyping/UG/combined_old/${sample}_combined.vcf.gz \
    -o data/alignments/genotyping/UG/combined/${sample}_combined.vcf \
    -genotypeMergeOptions UNSORTED;
done < data/alignments/fastq/symlinks/samples.txt
```

## 13/7/2020

doing the same for the HC lines:

```bash
mv -v data/alignments/genotyping/HC/combined data/alignments/genotyping/HC/combined_old
mkdir -p data/alignments/genotyping/HC/combined

time while read sample; do
    time java -jar ./bin/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
    --variant ../alignments/genotyping/HC_diploid/all_dark_HC/all_dark_diploid.HC.variants.vcf.gz \
    --variant ../alignments/genotyping/HC_diploid/all_SL_HC/allSL.diploid.HC.variants.vcf.gz \
    --variant data/alignments/genotyping/HC/pairs/${sample}_samples.vcf.gz \
    --o data/alignments/genotyping/HC/combined/${sample}_combined.vcf \
    -genotypeMergeOptions UNSORTED;
done < data/alignments/fastq/symlinks/samples.txt
```

looks good - now to just bgzip and tabix

next up - recreating the DL41 and 46 but with pairs
switched to see if that fixes the issue

should also check combined VCFs to see if the observed
variants are persistent across the ancestral pooled DL samples

## 15/7/2020

queueing up DL 'test' pairs

```bash
mkdir -p data/alignments/genotyping/UG/DL_test

time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL41_0.bam \
-I data/alignments/bam/DL46_5.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/DL_test/DL41_46.vcf

time java -jar ./bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
-I data/alignments/bam/DL46_0.bam \
-I data/alignments/bam/DL41_5.bam \
-glm BOTH \
-ploidy 2 \
--output_mode EMIT_ALL_SITES \
--heterozygosity 0.02 \
--indel_heterozygosity 0.002 \
-o data/alignments/genotyping/UG/DL_test/DL46_41.vcf
```

## 16/7/2020

next up - bgzip, tabix, and then check for mutation counts:

```bash
for test in 41_46 46_41; do
    time python3.5 analysis/filtering/filter_candidate_muts.py \
    --vcf data/alignments/genotyping/UG/DL_test/DL${test}.vcf.gz \
    --out_format table \
    --out data/alignments/genotyping/UG/DL_test/DL${test}_muts.txt
done
```

would you look at that:

```
$ wc -l data/alignments/genotyping/UG/DL_test/*txt
  15 data/alignments/genotyping/UG/DL_test/DL41_46_muts.txt
  16 data/alignments/genotyping/UG/DL_test/DL46_41_muts.txt
  31 total
```


## 20/7/2020

trying haplotypecaller again with the `--emitRefConfidence` option to see if that
affects anything - though it seems to only work on one sample at a time

testing:

```bash
time java -jar ./bin/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
        -I data/alignments/bam/CC1373_0.bam \
        -L chromosome_1 \
        -ploidy 2 \
        --emitRefConfidence BP_RESOLUTION \
        --output_mode EMIT_ALL_SITES \ 
        --heterozygosity 0.02 \
        --indel_heterozygosity 0.002 \
        -o test.vcf
```

## 5/4/2021

picking this up many moons later - need to get these invariant sites for
our 'denominator' (e.g. the number of callable sites) - and
GenotypeGVCFs is probably the actual way to do it! 

```bash
mkdir -p data/alignments/gvcfs/

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/CC1373_0.bam \
    -L chromosome_1:1-1000000 \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \ # required due to different gvcf indexing method
    -variant_index_parameter 128000 \
    -o test.g.vcf

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/CC1373_5.bam \
    -L chromosome_1:1-1000000 \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \ # required due to different gvcf indexing method
    -variant_index_parameter 128000 \
    -o test_2.g.vcf

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant test.g.vcf \
    --variant test_2.g.vcf \
    -L chromosome_1:1-1000000 \
    --includeNonVariantSites \
    -o combined_test.vcf
```

so far so good - although it seems the combined file contains all 1m sites - going to have
to keep an eye out to see how much 'dropoff' happens as we call more of each sample

```bash
mkdir -p data/alignments/genotyping/HC_invariant

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/CC1373_0.bam \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o data/alignments/genotyping/gvcfs/CC1373_0.g.vcf

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/CC1373_5.bam \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o data/alignments/genotyping/gvcfs/CC1373_5.g.vcf
```

I could parallelize this over regions like I did earlier, but it feels pretty
extraneous right now - this mostly just needs to run in the background while
I do other things

quick temp script (should work for all except the DL41/46 samples)

```bash
sample=$1
time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/${sample}_0.bam \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o data/alignments/genotyping/gvcfs/${sample}_0.g.vcf

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/${sample}_5.bam \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o data/alignments/genotyping/gvcfs/${sample}_5.g.vcf

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant data/alignments/genotyping/gvcfs/${sample}_0.g.vcf \
    --variant data/alignments/genotyping/gvcfs/${sample}_5.g.vcf \
    --includeNonVariantSites \
    -o data/alignments/genotyping/HC_invariant/${sample}.vcf
```

but I should wait for the 1373 ones to finish to make sure that they 
actually look as intended before running this over all the others! 

## 6/5/2021

HC is done - took between 13-16 hours between the two samples

now for the final step:

```bash
time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant data/alignments/genotyping/gvcfs/CC1373_0.g.vcf \
    --variant data/alignments/genotyping/gvcfs/CC1373_5.g.vcf \
    --includeNonVariantSites \
    -o data/alignments/genotyping/HC_invariant/CC1373.vcf
```

done in 30 min - looks like this 'contains all sites', but the INFO
column lists sites without calls (`NCC=2`, or `NCC=1` if one missing)

bgzipping and tabixing the combined file (it's 7.4 GB!) before continuing:

I think I'll need the script to take in 0 or 5 as an arg - that way I
can have two running at all times

`gvcf_temp.sh`:

```bash
sample=$1
line=$2

time /usr/bin/java -jar bin/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    -I data/alignments/bam/${sample}_${line}.bam \
    -ploidy 2 \
    --emitRefConfidence GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -o data/alignments/genotyping/gvcfs/${sample}_${line}.g.vcf
```

followed by:

```bash
time bash gvcf_temp.sh CC1952 0
time bash gvcf_temp.sh CC1952 5 # in separate shell
```

## 8/5/2021

so apparently there was a typo..... I wrote `${lines}` instead of
`${line}` in the final `-o` argument, which meant both the 0 and 5 files
were written to the same gvcf...

here we go again:

```bash
time bash gvcf_temp.sh CC1952 0
```

the server is being hammered right now, so I might not queue up the 5 file just yet

## 9/5/2021

queueing up 5 and 2342 0:

```bash
time bash gvcf_temp.sh CC1952 5
time bash gvcf_temp.sh CC2342 0
```

## 10/5/2021

next up:

```bash
time bash gvcf_temp.sh CC2342 5
time bash gvcf_temp.sh CC2344 0
```

## 12/5/2021

happy post-vax day!

```bash
time bash gvcf_temp.sh CC2344 5
time bash gvcf_temp.sh CC2931 0
```

## 14/5/2021

```bash
time bash gvcf_temp.sh CC2931 5
time bash gvcf_temp.sh CC2935 0
```

## 15/5/2021

need to run CombineVariants for `DL41_46` and vice versa

```bash
for sample in DL41_46 DL46_41; do
    time /usr/bin/java -jar ./bin/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/dark_ancestors/dark_aligned_dark_UG_diploid.vcf.gz \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/salt/salt_aligned_salt_UG.diploid.vcf.gz \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/CC2935_diploid.UG.vcf.gz \
    --variant ../alignments/genotyping/unified_genotyper/UG_diploid/wt_ancestors/wt_diploid.UG.vcf.gz \
    --variant data/alignments/genotyping/UG/pairs/${sample}_samples.vcf.gz \
    --o data/alignments/genotyping/UG/combined/${sample}_combined.vcf \
    -genotypeMergeOptions UNSORTED;
done
```

## 16/5/2021

moving on:

```bash
# bgzip and tabix
bgzip data/alignments/genotyping/UG/combined/DL41_46_combined.vcf
bgzip data/alignments/genotyping/UG/combined/DL46_41_combined.vcf
tabix -p vcf data/alignments/genotyping/UG/combined/DL41_46_combined.vcf.gz
tabix -p vcf data/alignments/genotyping/UG/combined/DL46_41_combined.vcf.gz

# and then the gvcf-ing continues
time bash gvcf_temp.sh CC2935 5
time bash gvcf_temp.sh CC2937 0
```

## 18/5/2021

```bash
time bash gvcf_temp.sh CC2937 5
time bash gvcf_temp.sh DL40 0
```

## 19/5/2021

DL40 was done in an hour and 40 minutes?? file looks legit,
so I'm surprised it went so fast 

```bash
time bash gvcf_temp.sh DL40 5
time bash gvcf_temp.sh DL41 0
```

## 20/5/2021

seems the DLs are continuing the theme of being done in 1.5-2 hours

```bash
time bash gvcf_temp.sh DL41 5
time bash gvcf_temp.sh DL46 0
time bash gvcf_temp.sh DL46 5
```

shoot - I also have to CombineVariants the `HC/combined` files

```bash
for sample in DL41_46 DL46_41; do
    time /usr/bin/java -jar ./bin/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R data/references/chlamy.5.3.w_organelles_mtMinus.fasta \
    --variant ../alignments/genotyping/HC_diploid/anc_wt_HC/anc_wt_diploid.HC.variants.vcf.gz \
    --variant ../alignments/genotyping/HC_diploid/all_dark_HC/all_dark_diploid.HC.variants.vcf.gz \
    --variant ../alignments/genotyping/HC_diploid/all_SL_HC/allSL.diploid.HC.variants.vcf.gz \
    --variant data/alignments/genotyping/HC/pairs/${sample}_samples.vcf.gz \
    --o data/alignments/genotyping/HC/combined/${sample}_combined.vcf \
    -genotypeMergeOptions UNSORTED;
done
```

done in 7 min! 

## 21/5/2021

```bash
time bash gvcf_temp.sh DL51 0
time bash gvcf_temp.sh DL51 5
```

## 22/5/2021

```bash
time bash gvcf_temp.sh DL53 0
time bash gvcf_temp.sh DL53 5
```

## 23/5/2021

```bash
time bash gvcf_temp.sh DL55 0
time bash gvcf_temp.sh DL55 5
time bash gvcf_temp.sh DL57 0
time bash gvcf_temp.sh DL57 5
```

## 24/5/2021

```bash
time bash gvcf_temp.sh DL58 0
time bash gvcf_temp.sh DL58 5
time bash gvcf_temp.sh SL26 5 # no 0
time bash gvcf_temp.sh SL27 0
time bash gvcf_temp.sh SL27 5
```

seems the SL lines took longer, interestingly enough (10 hrs for 26,
4.5 hrs for SL27 0)





