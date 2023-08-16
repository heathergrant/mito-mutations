#!/bin/bash
#SBATCH --job-name=vep
#SBATCH --time=5:00:00
#### #SBATCH --cpus-per-task=1
#SBATCH --output=vep.log
#SBATCH --mem=32G

#module load SAMtools/1.11
#module load GATK/4.1.9.0
#conda activate variant_calling

conda activate ensembl-vep101.0

#step 15 - annotation with VEP
vep --offline -i vep/MT_IBD002_ALL_SAMPLES.vcf.gz --format vcf --vcf -o vep/MT_IBD002_AL_SAMPLES_VEP.vcf.gz\
	--dir /data/scratch/DMP/UCEC/GENEVOD/repository/vep --assembly GRCh38 --fasta chrM_refs/hg38_v0_chrM_Homo_sapiens_assembly38.chrM.fasta --everything -v


gatk VariantsToTable \
     -V vep/MT_IBD002_ALL_SAMPLES.vcf.gz  \
     -F CHROM -F POS -F REF -F ALT -F TYPE -GF AD -GF AF -GF DP \
     -O vep/MT_IBD002_ALL_SAMPLES.table
