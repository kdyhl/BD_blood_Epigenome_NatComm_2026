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
CUTOFF="emp.p.fdr.cutoff0.05"

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


#TYPE="Q0.1"
HMs=c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3")

dir.in.plink = "~/lhou.compbio/data/LDSC/1000G_EUR_Phase3_plink/"
dir.in.snp = "~/lhou.compbio/data/LDSC/hapmap3_snps/"

files.in.haQTLPeak.RData = paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData")
names(files.in.haQTLPeak.RData) = HMs
files.in.pk.hg38tohg19.RDS = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".pk.hg38tohg19.RDS")
names(files.in.pk.hg38tohg19.RDS) = HMs


dir.tmp=paste0("~/hptmp/MayoBipolar/ldsc_a_annot_gARE/")
dir.create(dir.tmp, showWarnings = F, recursive = T)

file.ou.sh.pref = "a.getAnnotScore."
dir.ou.annot = paste0("a_annotations_gARE_", CUTOFF, "/")
dir.create(dir.ou.annot, showWarnings = F, recursive = T)

#file.ou.RDS = paste(dir.ou.annot, "totalAnnots.RDS", sep="")




files.tmp.bed=c()
for(hm in HMs)
{
  print(hm)
  files.tmp.bed[paste0(hm, ".gARE.hg19")]=paste0(dir.tmp, hm, ".gARE.hg19.bed")
  load(files.in.haQTLPeak.RData[hm])
  pks.hg38tohg19=readRDS(files.in.pk.hg38tohg19.RDS[hm])
  
  pks.hg38=hQTLPeaks.list[[CUTOFF]]
  pks.hg19=pks.hg38tohg19[pks.hg38]
  
  pks.hg19.bed=gsub(":|-", "\t", pks.hg19)
  write(pks.hg19.bed, file=files.tmp.bed[paste0(hm, ".gARE.hg19")])
  
}




#calculate

for(nm in names(files.tmp.bed))
{

  for(chr in CHRs[1:22])  #need to make sure the total job submitted is less than 1000
  {
    f.o.sh = paste(file.ou.sh.pref, nm, ".", chr, ".sh", sep="")
    write(sh.head, f.o.sh)
    
    
    f.o.annot = paste(dir.ou.annot, nm, ".", chr, ".annot.gz", sep="")
    f.o.l2 = paste(dir.ou.annot, nm, ".", chr, sep="")
 
    cmd= paste("make_annot.py --bed-file ", files.tmp.bed[nm], 
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

