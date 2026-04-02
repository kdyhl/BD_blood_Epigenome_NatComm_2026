#include three different kinds of annotation set
#include only limited number of GWAS
#ARE grp only coming from 3 active marks


#remove the filter snp step, since I rerun the cacluation of l2 for baseline model, and now snps set matching problem should be solved
#use GWAS atlas data, some redundant summary statistics


#for the baseline they get different number of SNPs as I got.
#need to use the overlapping snps, ignore possible change to .M, .M_5_50
#should guarantee all the new annotations share the same set of SNPs


CHRs = as.character(1:22)


sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                #"#$ -P compbio_lab", 
                "#$ -binding linear:1", 
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=20g",
                "#$ -l h_rt=10:00:00", 
                "##$ -pe smp 5",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "export PATH=$HOME/lhou.compbio/software/anaconda2/bin:$PATH",
                "source activate ldsc",
                "export PATH=$HOME/lhou.compbio/software/ldsc/:$PATH",
                sep="\n")


# dir.in.plink = "~/lhou.compbio/data/LDSC/1000G_EUR_Phase3_plink/"
# dir.in.snp = "~/lhou.compbio/data/LDSC/hapmap3_snps/"
file.in.baseline.pref = "~/lhou.compbio/data/LDSC/1000G_Phase3_baselineLD_v2.1_ldscores/baselineLD."
file.in.weight.pref = "~/lhou.compbio/data/LDSC/weights_hm3_no_hla/weights."
file.in.freq.pref = "~/lhou.compbio/data/LDSC/1000G_Phase3_frq/1000G.EUR.QC."


files.in.GWAS=c()
files.in.GWAS["PGC_SCZ"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_pgc.scz2.bed.gz"
files.in.GWAS["RA_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1464_self_reported_rheumatoid_arthritis.bed.gz"
files.in.GWAS["IBD_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1461_self_reported_inflammatory_bowel_disease.bed.gz"
files.in.GWAS["MS_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1261_self_reported_multiple_sclerosis.bed.gz"
files.in.GWAS["PD_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1262_self_reported_parkinsons_disease.bed.gz"
files.in.GWAS["T1D_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1222_self_reported_type_1_diabetes.bed.gz"
files.in.GWAS["T2D_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1223_self_reported_type_2_diabetes.bed.gz"
files.in.GWAS["Neurot_UKB"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20127_Neuroticism_score.bed.gz"
files.in.GWAS["PGC_BP"]="~/lhou.compbio/data/GWAS/pgc_bipolar/daner_PGC_BIP32b_mds7a_0416a.hg38.sumstats.gz"
files.in.GWAS["AD_NG2022"]="~/lhou.compbio/data/GWAS/AD/AD_Bellenguez_NG2022/GCST90027158_buildGRCh38.sumstats.gz"
files.in.GWAS["HighBloodPressure"]="~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/highBloodPressure_6150_4_logistic.EUR.sumstats.gz"
files.in.GWAS["SmokeInit"]="~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/SmokingInitiation.txt.hg38.sumstats.gz"
files.in.GWAS["DrinkPerWeek"]="~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/DrinksPerWeek.txt.hg38.sumstats.gz"
files.in.GWAS["COVID19.HG"] = "~/lhou.compbio/data/GWAS/COVID19_hg/COVID19_HGI_ANA_C2_V2.sumstats.gz"

files.in.GWAS["height_GIANT"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_GIANT_HEIGHT.bed.gz"
files.in.GWAS["MDD_PGC"]= "~/lhou.compbio/data/GWAS/PGC/PGC_UKB_MDD_2019_forLDSC.sumstats.gz"
files.in.GWAS["crossDisorder_PGC"]= "~/lhou.compbio/data/GWAS/PGC/pgc_cdg2_meta_no23andMe_oct2019_v2.daner.forLDSC.sumstats.gz"
files.in.GWAS["CRP"]= "~/lhou.compbio/data/GWAS/GWAS.Catalog/CRP_GCST90029070_buildGRCh37.pval_floored.forLDSC.sumstats.gz"
files.in.GWAS["neutrophil"]= "~/lhou.compbio/data/GWAS/GWAS.Catalog/neutropCounts_GCST90691823.h.forLDSC.sumstats.gz"



# files.in.gwasSummary["FEV1"] = "~/lhou.compbio/data/GWAS/Lung.function_fromXushen/Shrine_30804560_FEV1_meta-analysis.parse.txt.gz"
# files.in.gwasSummary["FVC"] = "~/lhou.compbio/data/GWAS/Lung.function_fromXushen/Shrine_30804560_FVC_meta-analysis.parse.txt.gz"
# files.in.gwasSummary["PEF"] = "~/lhou.compbio/data/GWAS/Lung.function_fromXushen/Shrine_30804560_PEF_meta-analysis.parse.txt.gz"




dirs.in.annot = c("a_annotations_bulkDP/limma_V1.4.1/", "a_annotations_gARE_emp.p.fdr.cutoff0.05/", "a_annotations_3actHM_AREGrp/")

# dir.tmp.annot = "~/hptmp/eGTEx-H3K27ac/ldsc_a_2_haQTL_V1.2/" #annot_ovlpWithBaseline
# if(!dir.exists(dir.tmp.annot))
# {
#   dir.create(dir.tmp.annot)
# }
# dir.tmp.annot.baseline = paste0(dir.tmp.annot, "baseline/")
# if(!dir.exists(dir.tmp.annot.baseline))
# {
#   dir.create(dir.tmp.annot.baseline)
# }


dir.ou = paste0("a_2_heritPartition_byLDSC_AREGrp.bulkDP.gCRE_slctGWAS/")
dir.create(dir.ou, showWarnings = F, recursive = T)

#annotation of interest by regular expression



for(d.i in dirs.in.annot)
{
  annot.nms = levels(factor(sapply(dir(d.i, paste0("annot.gz"), full.names=F), #"haQTL.*\\.1.annot.gz"  "diffPeak.*\\.1.annot.gz""
  FUN=function(f.i)
  {
   return(gsub(".\\d+.annot.gz", "", f.i, perl=T))
  })))
  
  
  
  
  for(annot.nm in annot.nms)#  #[1:8]  [25:29] [9:16]#[17:24]
  {
    print(annot.nm)
  #LDSC --h2
    # f.o.sh = paste0("a_2_LDSCPartit_", annot.nm, ".", TISSUE, ".sh")
    # write(sh.head, f.o.sh)
    # 
    f.o.sh = paste0("a_2_LDSCPartit_", annot.nm, ".sh")
    cmds=""
     
    for(pheno in names(files.in.GWAS)) 
    {
      #d.tmp.annot.nm = paste0(dir.tmp.annot, annot.nm, "/")
      f.o.res = paste0(dir.ou, pheno,"_", annot.nm)
      f.o.log = paste0(f.o.res, ".log")
      if(file.exists(f.o.log) &&
         length(system(paste0("grep \"Analysis finished\" ", f.o.log), intern = T))==1 )
      {
        next
      }
        
      
      
      cmd = paste0("ldsc.py",
                 " --h2 ", files.in.GWAS[pheno],
                 #" --ref-ld-chr ", dir.tmp.annot.baseline, "baselineLD.,", d.tmp.annot.nm, annot.nm, ".", 
                  " --ref-ld-chr ", file.in.baseline.pref, ",", d.i, annot.nm, ".", 
                 " --w-ld-chr ", file.in.weight.pref,
                 " --overlap-annot",
                 " --frqfile-chr ", file.in.freq.pref,
                 " --out ", f.o.res
                  )
                 
      cmds = paste(cmds, cmd, sep = "\n")
      
      
    }
    if(cmds!="")
    {
      write(sh.head, f.o.sh)
      write(cmds, f.o.sh, append = T)    
      write("echo done!", f.o.sh, append = T)    
      system(paste("qsub " , f.o.sh, sep=""))
    } 
      
  }

  
}
