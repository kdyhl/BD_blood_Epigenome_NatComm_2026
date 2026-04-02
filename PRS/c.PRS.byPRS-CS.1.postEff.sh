#!/bin/bash
#$ -S /bin/bash
#$ -V -cwd
#$ -e ./error.$JOB_NAME.$JOB_ID
#$ -o ./outpt.$JOB_NAME.$JOB_ID
#$ -l h_vmem=20g
#$ -l h_rt=5:00:00
##$ -pe smp 4
source ~/.bashrc
source ~/.my.bashrc


# export LD_LIBRARY_PATH=$LIBRARY_PATH:/broad/uge/8.4.3/lib/lx-amd64 
source ~/lhou.compbio/software/anaconda3_2023/etc/profile.d/conda.sh
conda activate prscs_env


# GWAS summary input for PRS-CS
#zcat ./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.bed.gz |awk 'BEGIN{OFS="\t"; print "SNP", "A1", "A2", "BETA", "P"} NR>1 {print $4, $5, $6, $8, $7}' > ./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.forPRS-CS.txt


#calculate effective sample size for each SNP and median value for the input
#zcat ./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.bed.gz | awk 'NR>1 {print 4 / (1/$13 + 1/$14)}' |sort -n |awk ' { a[i++] = $1 } END { print a[int(i/2)] }'
#147632

for chr in {1..22}; do
   python ~/lhou.compbio/software/PRScs/PRScs.py \
   --ref_dir=/home/unix/leihou/lhou.compbio/data/LD/1K_Genome_byPRS-CS/ldblk_1kg_eur \
   --bim_prefix=./b_Mayo.BP_PGC.leaveOutMayo.GWAS_MAF0.05.locID_V1.1_noSwapping_hg38_QC/plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC \
   --sst_file=./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.forPRS-CS.txt \
   --n_gwas=147632 \
   --chrom=$chr \
   --out_dir=c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/BD_Mayo
done

echo "job done"
