#!/bin/bash
#SBATCH --job-name=filter
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-21
#SBATCH --output=filter.log

module load SAMtools/1.11
module load GATK/4.1.9.0
conda activate variant_calling

# Extract information from samplelist.txt
filelist=sample_list_ibd.txt
bamlist=all_bam_list_ibd.txt
TUM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" $filelist | awk '{print $1}')"
INPUTBAM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" $bamlist | awk '{print $1}')"
NORM=IBD002_B_S1

#  run custom python filter
python3 filter_mutec_mt.py vcfs/${TUM}_MERGED.vcf vcfs/${TUM}_PY_filtered.vcf

bgzip -c vcfs/${TUM}_PY_filtered.vcf > vcfs/${TUM}_PY_filtered.vcf.gz
tabix vcfs/${TUM}_PY_filtered.vcf.gz

#Do last step outside of the slurm array
module load BCFtools/1.11
#step 15 merge with bcftools
bcftools merge --force-samples vcfs/*_PY_filtered.vcf.gz | bgzip -c > vep/MT_IBD002_ALL_SAMPLES.vcf.gz
tabix -p vcf vep/MT_IBD002_ALL_SAMPLES.vcf.gz
