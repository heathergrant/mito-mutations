#!/bin/bash
#SBATCH --job-name=extract-mito
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-21
#SBATCH --output=logs/make_new_mitobams_%J

module load GATK/4.1.9.0
conda activate gatk4.1.9.0

# Extract information from samplelist.txt
filelist=sample_list_ibd.txt
bamlist=all_bam_list_ibd.txt
TUM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" $filelist | awk '{print $1}')"
INPUTBAM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" $bamlist | awk '{print $1}')"

#step 1 Subset Bam to CHrm
gatk --java-options -Xmx2G PrintReads -L chrM -I ${INPUTBAM} -O mitobams/mito_${TUM}.bam

#step 2 revert sam
module load picard-tools/2.23.8

java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar RevertSam INPUT=mitobams/mito_${TUM}.bam OUTPUT_BY_READGROUP=false OUTPUT=mitobams/${TUM}.revertsam

rm mitobams/mito_${TUM}.bam
rm mitobams/mito_${TUM}.bai

#step 3 sam to fastq
java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar SamToFastq INPUT=mitobams/${TUM}.revertsam FASTQ=mitobams/${TUM}_mito.fastq INTERLEAVE=true NON_PF=true VALIDATION_STRINGENCY=SILENT

rm mitobams/${TUM}.revertsam

#step4 bwa mem and pipe to new bam file
module load BWA/0.7.17
module load SAMtools/1.11
bwa mem -K 100000000 -M -R "@RG\tID:${TUM}\tCN:ICR_TPU\tPU:1\tSM:${TUM}\tLB:${TUM}\tPL:ILLUMINA" -p -t 8 \
chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta mitobams/${TUM}_mito.fastq | samtools view --threads 8 -o mitobams/${TUM}_mt_ref.bam

#step 5 bwa mem and pipe to new bam file WITH SHIFTED MT REFERENC E  # no interval flag

bwa mem -K 100000000 -M -R "@RG\tID:${TUM}\tCN:ICR_TPU\tPU:1\tSM:${TUM}\tLB:${TUM}\tPL:ILLUMINA" -p -t 8 \
chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta mitobams/${TUM}_mito.fastq | samtools view --threads 8 -o mitobams/${TUM}_mt_shifted_ref.bam

rm mitobams/${TUM}_mito.fastq

#step 7 Mark duplicates
java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar  MarkDuplicates \
      INPUT=mitobams/${TUM}_mt_ref.bam   OUTPUT=mitobams/${TUM}_mt_ref_md.bam \
      METRICS_FILE=logs/${TUM}.mark.duplicates.log \
      VALIDATION_STRINGENCY=SILENT  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname"  CLEAR_DT="false" ADD_PG_TAG_TO_READS=false

rm mitobams/${TUM}_mt_ref.bam

java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar  MarkDuplicates \
      INPUT=mitobams/${TUM}_mt_shifted_ref.bam   OUTPUT=mitobams/${TUM}_mt_shifted_ref_md.bam \
      METRICS_FILE=logs/${TUM}.shifted_mark.duplicates.log \
      VALIDATION_STRINGENCY=SILENT  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname"  CLEAR_DT="false" ADD_PG_TAG_TO_READS=false

rm mitobams/${TUM}_mt_shifted_ref.bam

# step 8 Sort Sam
java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar SortSam \
INPUT=mitobams/${TUM}_mt_ref_md.bam OUTPUT=mitobams/${TUM}_sorted.bam SORT_ORDER='coordinate' CREATE_INDEX=true

rm mitobams/${TUM}_mt_ref_md.bam

java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar SortSam \
INPUT=mitobams/${TUM}_mt_shifted_ref_md.bam OUTPUT=mitobams/${TUM}_shifted_sorted.bam SORT_ORDER='coordinate' CREATE_INDEX=true

rm mitobams/${TUM}_mt_shifted_ref_md.bam
