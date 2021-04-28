# alignment.sh - QC and trim fastqs + align to reference
# run from project root

### DIR SETUP ###

mkdir -p data/references
mkdir -p data/alignments/bam
mkdir -p data/alignments/fastq
mkdir -p data/alignments/fastq/symlinks
mkdir -p data/alignments/fastq/fastqc
mkdir -p data/alignments/fastq/fastqc/raw
mkdir -p data/alignments/fastq_trim
mkdir -p data/alignments/fastq_trim/fastqc
mkdir -p data/alignments/fastq_trim/fastqc/raw

# download fastq files to data/alignments/fastq

for fname in $(find data/alignments/fastq -name '*fq.gz'); do
    ln -sv ${fname} data/alignments/fastq/symlinks/
done

# see log 20/5 for symlink name cleaning with python

touch data/alignments/fastq/symlinks
for read in 1 2; do
    for fname in data/alignments/fastq/symlinks/*_${read}.fq.gz; do 
        basename ${fname} _${read}.fq.gz;
    done
done

### READ QC AND TRIMMING ###

fastqc --kmers 7 --outdir data/alignments/fastq/fastqc/raw \
    data/alignments/fastq/symlinks/*fq.gz

for fname in data/alignments/fastq/fastqc/raw/*zip; do
    unzip ${fname} -d data/alignments/fastq/fastqc;
done

# trimming with trimmomatic
time while read fname; do
    time java -jar bin/trimmomatic PE -threads 4 -phred33 \
        data/alignments/fastq/symlinks/${fname}_1.fq.gz \
        data/alignments/fastq/symlinks/${fname}_2.fq.gz \
        data/alignments/fastq_trim/${fname}_trim_1.fq.gz \
        data/alignments/fastq_trim/${fname}_trim_unpaired_1.fq.gz \
        data/alignments/fastq_trim/${fname}_trim_2.fq.gz \
        data/alignments/fastq_trim/${fname}_trim_unpaired_2.fq.gz \
        ILLUMINACLIP:bin/TruSeq2-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20;
    sleep 3;
done < data/alignments/fastq/symlinks/prefixes.txt

# fastqc trimmed reads

fastqc --kmers 7 --outdir data/alignments/fastq_trim/fastqc/raw \
    data/alignments/fastq_trim/??*_?_trim_?.fq.gz

for fname in data/alignments/fastq_trim/fastqc/raw/*zip; do
    unzip ${fname} -d data/alignments/fastq_trim/fastqc;
done

### ALIGNMENT ###

# bwa
time while read fname; do
    ref=data/references/chlamy.5.3.w_organelles_mtMinus.fasta
    time bwa mem -t 12 $ref \
        data/fastq_trim/${fname}_trim_1.fq.gz \
        data/fastq_trim/${fname}_trim_2.fq.gz | \
        samtools sort -@4 -T ${fname}.sorting.tmp -O bam \
        -o data/alignments/bam/${fname}.sorted.bam
done < data/alignments/fastq/symlinks/prefixes.txt

# fix mates
time while read fname; do
    time java -jar /opt/picard-tools/picard.jar FixMateInformation \
        I=data/alignments/bam/${fname}.sorted.bam \
        O=data/alignments/bam/${fname}.fixMate.bam \
        VALIDATION_STRINGENCY=LENIENT;
done < data/alignments/fastq/symlinks/prefixes.txt

# add read groups
time while read fname; do
    time java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups \
        I=data/alignments/bam/${fname}.fixMate.bam \
        O=data/alignments/bam/${fname}.RG.bam \
        RGID=${fname} \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=${fname} \
        VALIDATION_STRINGENCY=LENIENT;
    sleep 3;
done < data/alignments/fastq/symlinks/prefixes.txt

mkdir -p data/alignments/bam_final # temp dir

# mark duplicates
time while read fname; do
    time java -jar /opt/picard-tools/picard.jar MarkDuplicates \
        I=data/alignments/bam/${fname}.RG.bam \
        O=data/alignments/bam_final/${fname}.bam \
        METRICS_FILE=analysis/dup_metrics.txt
    sleep 3;
done < data/alignments/fastq/symlinks/prefixes.txt

# create index files
time for fname in data/alignments/bam_final/*.bam; do
    samtools index ${fname} data/alignments/bam_final/${fname}.bai;
done < data/alignments/fastq/symlinks/prefixes.txt

# reorganize directories and remove intermediate bam files
rm -rf data/alignments/bam/*
mv -v data/alignments/bam_final data/alignments/bam





