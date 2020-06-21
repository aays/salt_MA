
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




















