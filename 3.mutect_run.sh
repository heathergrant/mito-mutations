#!/bin/bash
#SBATCH --job-name=mutec
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-21
#SBATCH --output=logs/mutec_%J

module load GATK/4.1.9.0
conda activate gatk4.1.9.0

# Extract information from samplelist.txt
filelist=sample_list_ibd.txt
bamlist=all_bam_list_ibd.txt
TUM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" $filelist | awk '{print $1}')"
INPUTBAM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" $bamlist | awk '{print $1}')"
NORM=IBD002_B_S1

# make blood bams see (1.make_bloods.sh)
#steps 1 to 8 make the bam files (2.make_mitobams.sh)

#step 9 call mutec2 with mitochondria mode on

gatk --java-options -Xmx4G Mutect2 --reference chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta \
	--input mitobams/${TUM}_sorted.bam \
       	--input mitobams/${NORM}_sorted.bam \
	--tumor-sample ${TUM} --normal-sample ${NORM} --output vcfs/${TUM}_raw.vcf \
	--read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \
	--mitochondria-mode -L chrM:576-16024

gatk --java-options -Xmx4G Mutect2 --reference chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
	--input mitobams/${TUM}_shifted_sorted.bam \
        --input mitobams/${NORM}_shifted_sorted.bam \
        --tumor-sample ${TUM} --normal-sample ${NORM} --output vcfs/${TUM}_shifted_raw.vcf\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \
        --mitochondria-mode -L chrM:8025-9144

#step 10 gatk filter mutec calls

gatk --java-options -Xmx4G FilterMutectCalls --variant vcfs/${TUM}_raw.vcf --output vcfs/${TUM}_filter.vcf \
	--reference chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta --mitochondria-mode

gatk --java-options -Xmx4G FilterMutectCalls --variant vcfs/${TUM}_shifted_raw.vcf --output vcfs/${TUM}_shifted_filter.vcf \
        --reference chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta --mitochondria-mode

#step 11 variant filtration - blacklisted sites from file

gatk VariantFiltration -V vcfs/${TUM}_filter.vcf \
        -O vcfs/${TUM}_vf.vcf \
        --apply-allele-specific-filters \
        --mask chrM_refs/hg38_v0_chrM_blacklist_sites.hg38.chrM.bed \
        --mask-name "blacklisted_site"

gatk VariantFiltration -V vcfs/${TUM}_shifted_filter.vcf \
        -O vcfs/${TUM}_shifted_vf.vcf \
        --apply-allele-specific-filters \
        --mask chrM_refs/hg38_v0_chrM_blacklist_sites.hg38.chrM.bed \
        --mask-name "blacklisted_site"

#step 12 Left align and trim to simplify the indel calls

gatk LeftAlignAndTrimVariants \
      -R chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta \
      -V vcfs/${TUM}_vf.vcf \
      -O vcfs/${TUM}_split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

gatk LeftAlignAndTrimVariants \
      -R chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta \
      -V vcfs/${TUM}_shifted_vf.vcf \
      -O vcfs/${TUM}_shifted_split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

#step 13 merge the first reference vcfs with the shifted reference vcfs

module load picard-tools/2.23.8

java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar LiftoverVcf I=vcfs/${TUM}_shifted_split.vcf \
	O=vcfs/${TUM}_SHIFTED_BACK.vcf\
	R=chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta \
	CHAIN=chrM_refs/hg38_v0_chrM_ShiftBack.chain\
	REJECT=vcfs/${TUM}_rejected.vcf

#merge the shifted back vcf and the regular vcf
java -Xmx4g -jar /opt/software/applications/picard-tools/2.23.8/picard.jar MergeVcfs \
	I=vcfs/${TUM}_SHIFTED_BACK.vcf\
	I=vcfs/${TUM}_split.vcf\
	O=vcfs/${TUM}_MERGED.vcf

#step 14 bzip

module load SAMtools/1.11
bgzip -c vcfs/${TUM}_MERGED.vcf > vcfs/${TUM}_MERGED.vcf.gz
tabix vcfs/${TUM}_MERGED.vcf.gz

