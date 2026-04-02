#including all the peaks

#args = commandArgs(trailingOnly=TRUE)
#DISEASE =  args[2] # "CD", "UC", "RA" "AD"

# DIS2GRP=c(CD="IBD",
#           UC="IBD",
#           RA="IBD",
#           AD="AD")
# GRP =DIS2GRP[DISEASE]           


# HM="H3K27ac"
# 
# REFs.SLCT = c("E040","E039", "E048", "E047", "E029", "E032", "E046", "BPNeutrop")

# DISEASEs = list()
# DISEASEs$AD="AD"
# DISEASEs$IBD=c("CD", "UC", "RA")#c("CD", "UC", "RA", "AD") #for disease specific haQTL

CHRs = 1:22

#TYPE="Q0.1"


dir.in.plink = "~/lhou.compbio/data/LDSC/1000G_EUR_Phase3_plink/"
dir.in.snp = "~/lhou.compbio/data/LDSC/hapmap3_snps/"

files.in.bed = c()

dir.tmp = paste0("~/hptmp/MayoBipolar/ldsc_a_annot_CREGrpAndbulkDP/")
dir.create(dir.tmp, showWarnings = F, recursive = T)

file.ou.sh.pref = "a.getAnnotScore."
dir.ou.annot = paste0("a_annotations_CREGrpAndbulkDP/")
dir.create(dir.ou.annot, showWarnings = F, recursive = T)

#file.ou.RDS = paste(dir.ou.annot, "totalAnnots.RDS", sep="")

sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                #"#$ -P compbio_lab", 
                #"#$ -binding linear:1", 
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=10g",
               # "#$ -l h_rt=15:00:00", 
                "##$ -pe smp 5",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "source activate ldsc",
                sep="\n")




# ###differential peaks#######################################
#PEAK.WIN =2000

files.in.DP.bed = system(paste0("find ../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/*/*.diseasedisease.Q0.05.*.bed"), intern = T) #paste0("../deconv/d_2_diffPeak_deconv_real_rmBE_rmBulkDP_V1.6.7.1_addCov_withIntrcpt_", DISEASE.GRP, "_summary/diffPeak.") #"../untilBatch8.2.2018/c_2_diffPeak_deconv_real_V1.6.1_rmBE_AD_summary/diffPeak." #"../untilBatch8.2.2018/c_2_diffPeak_deconv_real_V1.6.5_rmBE_lmWithIntrcpt_allPks_IBD_summary/diffPeak."
names(files.in.DP.bed)= sapply(files.in.DP.bed,
FUN=function(f.i.bed)
{
  gsub(".diffPeak.diseasedisease.Q0.05|.bed", "", basename(f.i.bed))
})
files.in.CREGrp.bed= system(paste0("find ../peakVariationAcrossTiss/a_3_AREModulesInEpimap/H3K27ac/AREGrp.*.bed"), intern = T) #paste0("../deconv/d_2_diffPeak_deconv_real_rmBE_rmBulkDP_V1.6.7.1_addCov_withIntrcpt_", DISEASE.GRP, "_summary/diffPeak.") #"../untilBatch8.2.2018/c_2_diffPeak_deconv_real_V1.6.1_rmBE_AD_summary/diffPeak." #"../untilBatch8.2.2018/c_2_diffPeak_deconv_real_V1.6.5_rmBE_lmWithIntrcpt_allPks_IBD_summary/diffPeak."
names(files.in.CREGrp.bed)= sapply(files.in.CREGrp.bed,
FUN=function(f.i.bed)
{
  gsub(".bed", "", basename(f.i.bed))
})


files.in.bed=c(files.in.DP.bed, files.in.CREGrp.bed)
# 
# files.in.bed=c()
# for(nm in names(files.in.all))
# {
#   annot.nm = paste(nm, ".",  PEAK.WIN/1000, "kbp",  sep="")
#   f.tmp = paste(dir.tmp, annot.nm, ".bed", sep="")
#   
#   
#   cmd = paste("awk  '(NF>=3){print $1 \"\\t\" ($2-", PEAK.WIN, "<0?0:$2-", PEAK.WIN, ") \"\\t\" $3+", PEAK.WIN, " }' ", files.in.all[nm], " >", f.tmp, sep="") #remove those empty files
#   system(cmd)
#   
#   peak.NO = as.numeric(unlist(strsplit(system(paste0("wc -l ", f.tmp), intern = T), split=" "))[1])
#   if(peak.NO >=10)
#   {
#     files.in.bed[paste0(annot.nm, "_", peak.NO)] = f.tmp
#   }
# 
# }
#   






#calculate

for(nm in names(files.in.bed))
{

  for(chr in CHRs[1:22])  #need to make sure the total job submitted is less than 1000
  {
    f.o.sh = paste(file.ou.sh.pref, nm, ".", chr, ".sh", sep="")
    write(sh.head, f.o.sh)
    
    
    f.o.annot = paste(dir.ou.annot, nm, ".", chr, ".annot.gz", sep="")
    f.o.l2 = paste(dir.ou.annot, nm, ".", chr, sep="")
 
    cmd= paste("make_annot.py --bed-file ", files.in.bed[nm], 
               " --bimfile ", dir.in.plink, "/1000G.EUR.QC.", chr, ".bim", 
               " --annot-file ", f.o.annot,
               sep="")
    write(cmd, f.o.sh, append = T)    
    
    cmd= paste("ldsc.py --l2 --bfile ", dir.in.plink, "/1000G.EUR.QC.", chr, 
               " --ld-wind-cm 1 --annot ", f.o.annot, 
               " --thin-annot --out ", f.o.l2,
               " --print-snps ", dir.in.snp, "/hm.", chr, ".snp",
               sep="")
    write(cmd, f.o.sh, append = T)    
    
    cmd = "echo done!"
    write(cmd, f.o.sh, append = T)    
    
    system(paste("qsub " , f.o.sh, sep=""))

  }  
  
}

