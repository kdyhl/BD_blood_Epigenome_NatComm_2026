#V1.1
#a.testBulkHaQTLVSGWAS_byMR_V1.1_allGWAS.R

#shared LD structure since GWAS is imputed based on GTEx EUR cohort

args = commandArgs(trailingOnly=TRUE)
HM = args[1] #H3K27ac, H3K4me1 ...
# #print(HM)
GWAS.NM=args[2]

if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))




options(scipen = 999)

#library(data.table)
#library(pheatmap)
#library("RColorBrewer")
# library(grid)
# library(gridExtra)
library(MendelianRandomization)
library(ggplot2)
library(RColorBrewer)
#cols <- rev(brewer.pal(11, 'RdYlBu'))
MY.PALETTE <- colorRampPalette(c("royalblue2", "green1", "orange2", "darkred")) #rev(brewer.pal(11, "Spectral"))
SC_COLOR = scale_colour_gradientn(colours = MY.PALETTE(8), limits=c(0,1)) 

source("~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R")




PK.WIND = 100000
GWAS.P.CUTOFF = 1e-6
GWAS.Z.CUTOFF = abs(qnorm(GWAS.P.CUTOFF/2))
hQTL.EMP.FDR.CUTOFF=0.2
#
SNP.PRUNE.R2.CUTOFF = 0.6
SNP.MIN.NUM = 3
haQTL.QVAL.CUTOFF =0.2


HM2sampN = c(H3K27ac=176,
             H3K4me1=149,
             H3K4me3=165,
             H3K36me3=165,
             H3K27me3=178
               )

HM2COVNum = c(H3K27ac = 10,
              H3K36me3 = 10,
              H3K4me1 = 10,
              H3K27me3 = 10,
              H3K4me3 = 10)

#haQTL.PVAL.CUTOFF =0.001

#GWAS info
files.in.GWAS.gz = c(PGC.BP = "~/lhou.compbio/data/GWAS/PGC_leaveMayoBDout/BDlooMAYO_PGC3_hg38.sorted.bed.gz", #"~/lhou.compbio/data/GWAS/pgc_bipolar/pgc_bip_hg38.sorted.bed.gz",
                     GWASAtlas_ukb2.highBloodPressure="~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/highBloodPressure_6150_4_logistic.EUR.hg38.sorted.bed.gz",
                     GWASAtlas_ukb2.diabetes="~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/diabetes_f.2443.0.0_logistic.EUR.hg38.sorted.bed.gz",
                     DrinksPerWeek = "~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/DrinksPerWeek.txt.hg38.sorted.bed.gz",
                     SmokingInitiation = "~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/SmokingInitiation.txt.hg38.sorted.bed.gz")


#zcat ~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/diabetes_f.2443.0.0_logistic.EUR.hg38.sorted.bed.gz|awk '($11!="N"){if(max<$11){max=$11}}END{print max}'
GWAS2N = c(PGC.BP=411505, #20352+31358,
           GWASAtlas_ukb2.highBloodPressure=385699,
           GWASAtlas_ukb2.diabetes=385420,
           DrinksPerWeek=537349,
           SmokingInitiation=632802) #https://www.nature.com/articles/s41588-019-0397-8

GWAS2type=c(PGC.BP="cc",
            GWASAtlas_ukb2.highBloodPressure="quant",
            GWASAtlas_ukb2.diabetes="cc",
            DrinksPerWeek="quant",
            SmokingInitiation="quant")
GWAS2CaseProp=c(PGC.BP=40802/411505,
                GWASAtlas_ukb2.diabetes=18483/385420
)

GWAS2beta=c(PGC.BP="BETA",
            GWASAtlas_ukb2.highBloodPressure="BETA",
            GWASAtlas_ukb2.diabetes="BETA",
            DrinksPerWeek="BETA",
            SmokingInitiation="BETA")

GWAS2se=c(PGC.BP="SE",
          GWASAtlas_ukb2.highBloodPressure="SE",
          GWASAtlas_ukb2.diabetes="SE",
          DrinksPerWeek="SE",
          SmokingInitiation="SE")
GWAS2P=c(PGC.BP="P",
         GWASAtlas_ukb2.highBloodPressure="PVAL",
         GWASAtlas_ukb2.diabetes="PVAL",
         DrinksPerWeek="PVAL",
         SmokingInitiation="PVAL")
GWAS2P.ID = c(PGC.BP=7,
              GWASAtlas_ukb2.highBloodPressure=10,
              GWASAtlas_ukb2.diabetes=10,
              DrinksPerWeek=9,
              SmokingInitiation=9)



files.in.hQTL.bed =  system(paste0("find ../haQTL/a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, "_nominal/*sorted.bed.gz"), intern=T) #paste0("
names(files.in.hQTL.bed)= sapply(files.in.hQTL.bed,
FUN=function(f.i.hQTL)
{
  unlist(strsplit(basename(f.i.hQTL), split=".", fixed = T))[1]
})
#file.in.haQTLPeak.bed = paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, ".haQTLPeak.fdr0.05.hg38.bed")
file.in.hQTLPeak.RData = paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, ".haQTLPeak.cutoffs.hg38.RData")

dir.in.Mayo.genotype.plink = "~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/"
file.in.Mayo.genotype.vcf="~/lhou.compbio/data/Mayo_Bipolar/WGS/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID.vcf.gz"
#dir.in.1KG.genotype.plink ="/broad/compbio/data/1KG_phase3/eur/"

dir.tmp=paste0("~/hptmp/MayoBipolar/a_bulkhaQTLvsGWAS_byMR_V1.1/", GWAS.NM, "/", HM, "/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.hQTLPk.bed = paste0(dir.tmp, HM, ".hQTLPk.hg38.bed")
file.tmp.GWAS.hits.bed = paste0(dir.tmp, "GWAS.hits.hg38.bed")
file.tmp.qtlPk.ovlp.GWAS.hits.bed = paste0(dir.tmp, HM, ".QTLPK.ovlp.GWAS.hits.hg38.bed")

#file.tmp.eQTL.chr.RDS = paste0(dir.tmp, HM, ".", CHR, ".eQTL.RDS")
#file.tmp.haQTL.chr.RDS = paste0(dir.tmp, HM, ".", CHR, ".haQTL.RDS")

#file.tmp.peak.bed = paste0(dir.tmp, HM, ".", CHR, ".pk.hg38.bed")
#file.tmp.gene.bed = paste0(dir.tmp, HM, ".", CHR,".eGene.hg38.bed")

# file.tmp.haQTL.hg38.bed = paste0(dir.tmp, HM, ".", CHR,".haQTL.hg38.bed")
# file.tmp.haQTL.unmap.bed = paste0(dir.tmp, HM, ".", CHR,".haQTL.unmap.bed")

dir.ou = paste0("a_bulkhaQTLvsGWAS_byMR_Egger_V1.1_allGWAS/", GWAS.NM, "/", HM, "_100k/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.pdf =paste0(dir.ou, HM, ".bulkHaQTLvsGWAS.pdf")
file.ou.RData =paste0(dir.ou, HM, ".bulkHaQTLvsGWAS.RData")
#

load(file.in.hQTLPeak.RData)
#focus on those peak overlapped with strong signals in GWAS
#load(file.in.haQTL.RData)

hQTLPk.hg38= hQTLPeaks.list[[paste0("emp.p.fdr.cutoff", hQTL.EMP.FDR.CUTOFF)]]#rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.P"] <= hQTL.EMP.PVAL]
hQTLPk.hg38.bed = paste(gsub(":|-", "\t", hQTLPk.hg38), hQTLPk.hg38, sep="\t")
write(hQTLPk.hg38.bed, file.tmp.hQTLPk.bed)


cmd = paste0("zcat ", files.in.GWAS.gz[GWAS.NM], "|awk '($", GWAS2P.ID[GWAS.NM], "<", GWAS.P.CUTOFF, "){print }'>", file.tmp.GWAS.hits.bed)
system(cmd)

cmd = paste0("bedtools window -w ", PK.WIND, " -a ", file.tmp.hQTLPk.bed, " -b ", file.tmp.GWAS.hits.bed, ">", file.tmp.qtlPk.ovlp.GWAS.hits.bed)
system(cmd)
buf = read.table(file.tmp.qtlPk.ovlp.GWAS.hits.bed, "\t", header=F, row.names = NULL)
hQTLPeaks.ovlpGWAS = levels(factor(buf[,4])) #paste0(buf[,1], ":", buf[,2], "-", buf[,3])))



#prune function
prunedSortedSNPWithR2 = function(snps.sorted,  #sorted with increasing p-value
                                 corr.r2,
                                 r2.cutoff)
{
  snps.pruned = snps.sorted[1]
  
  for(snp.test in snps.sorted[-1])
  {
    if(all(!is.na(corr.r2[snp.test, snps.pruned])) && all(corr.r2[snp.test, snps.pruned]<= r2.cutoff))
    {
      snps.pruned = c(snps.pruned, snp.test)
    }
  }
  
  return(snps.pruned)
}


#hQTLpeak2signifQTL = hQTLpeak2signifQTL.bycutoffs[[paste0("emp.p.fdr.cutoff", hQTL.EMP.FDR.CUTOFF)]]
pdf(file.ou.pdf, width=12, height=8)
#layout(matrix(1:9, ncol=3))
hQTLPeaks.QTLvsGWAS=lapply(hQTLPeaks.ovlpGWAS,
FUN=function(hQTLPk)
{
  print(hQTLPk)
  f.tmp.ovlp.snps =paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".", PK.WIND/1000, "k.ovlpSNP.txt")
  f.tmp.ld.pre = paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".ovlpSNP")
  
  # f.tmp.haQTL.z = paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".haQTL.z.txt")
  # f.tmp.gwas.z = paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".gwas.z.txt")
  # 
  #f.o.eCaviar = paste0(d.o.gwas, GWAS.nm, ".", HM, ".", gsub(":", "-", hQTLPk), ".txt")
  
  #SNPs from haQTL
  buf = unlist(strsplit(hQTLPk, ":|-", perl=T))
  pk.chr= buf[1]
  pk.start = as.numeric(buf[2])
  pk.end = as.numeric(buf[3])
  query = paste0(pk.chr, ":", max(1, pk.start-PK.WIND), "-", pk.end+PK.WIND)
  hQTL.tab = seqminer::tabix.read.table(files.in.hQTL.bed[pk.chr], query) 
  hQTL.tab.filt = hQTL.tab[hQTL.tab$peak==hQTLPk,]
  hQTL.tab.filt$hqtl.qval = p.adjust(hQTL.tab.filt$pval, method="BH")
  rownames(hQTL.tab.filt) = hQTL.tab.filt$rs
  hQTL.tab.filt = hQTL.tab.filt[hQTL.tab.filt$hqtl.qval<=haQTL.QVAL.CUTOFF,]
  #hQTL.tab.filt$varbeta=(-abs(hQTL.tab.filt[, "beta"])/qt(hQTL.tab.filt[, "pval"]/2, df=HM2sampN[HM]-2-HM2COVNum[HM]))^2

  
  #SNPs from GWAS
  gwas.tab = seqminer::tabix.read.table(files.in.GWAS.gz[GWAS.NM], query) 
  rownames(gwas.tab) = paste(gwas.tab[,1], gwas.tab[,3], gwas.tab[,5], gwas.tab[,6], sep="_")
  gwas.tab.flip = gwas.tab
  gwas.tab.flip[,5]= gwas.tab[,6]
  gwas.tab.flip[,6]= gwas.tab[,5]
  rownames(gwas.tab.flip)  = paste(gwas.tab.flip[,1], gwas.tab.flip[,3], gwas.tab.flip[,5], gwas.tab.flip[,6], sep="_")
  gwas.tab.all =rbind(gwas.tab, gwas.tab.flip)
  

  snp.loc.ovlp = intersect(rownames(hQTL.tab.filt), rownames(gwas.tab.all))
  
  ################################################################
  ##r2 extracted
  if(length(snp.loc.ovlp)==0) return(NA)
  print(paste0(length(snp.loc.ovlp), " of ", nrow(hQTL.tab.filt), " haQTL is with GWAS zscore"))
  write(snp.loc.ovlp, f.tmp.ovlp.snps)
  
  # plot(haQTL.z[snps.locID.ovlp],
  #     gwas.z[snps.locID.ovlp], 
  #     xlab="haQTL zscore",
  #     ylab="GWAS zscore",
  #     main=paste(HM, GWAS.nm, hQTLPk, sep=" ")
  #    )
  
  #ld of overlapping SNP from eGTEx
  cmd = paste0("plink --r2 square  --bfile ", dir.in.Mayo.genotype.plink, pk.chr, " --extract ",
                 f.tmp.ovlp.snps,  " --write-snplist  --out ", f.tmp.ld.pre)
  system(cmd)

  snp.loc.ovlp.kept=rownames(read.table(paste0(f.tmp.ld.pre, ".snplist"), sep="\t", row.names = 1, header = F))
  ld.matrix=as.matrix(read.table(paste0(f.tmp.ld.pre, ".ld"), sep="\t", row.names = NULL, header = F))
  rownames(ld.matrix)=snp.loc.ovlp.kept
  colnames(ld.matrix)=snp.loc.ovlp.kept
  #snp.loc.ovlp.loc =sapply(snp.loc.ovlp, FUN=function(x) as.numeric(unlist(strsplit(x, "_"))[2]))
  snp.loc.ovlp=snp.loc.ovlp.kept
  #snp.loc.ovlp = intersect(hQTLpeak2signifQTL[[hQTLPk]]$V2, snp.loc.ovlp.kept)
  if(length(snp.loc.ovlp)==0) return(NA)
  
  
  #################################################################
  ##r2 extracted
  #need to filter snps whose r2 is not available
  #prunings snps step by step on those with strong correlation with lead haQTLs
  
  
  
  #prune haQTLs snps
  snp.loc.ovlp.sorted = snp.loc.ovlp[order(hQTL.tab.filt[snp.loc.ovlp, "pval"], decreasing = F)]
  snp.loc.ovlp.sorted.pruned = prunedSortedSNPWithR2(snp.loc.ovlp.sorted, 
                                                         ld.matrix[snp.loc.ovlp.sorted, snp.loc.ovlp.sorted],
                                                         SNP.PRUNE.R2.CUTOFF)
  if(length(snp.loc.ovlp.sorted.pruned) < SNP.MIN.NUM)
  {
    return(NA)
  }
  
  
  b_se = cbind(bx = hQTL.tab.filt[snp.loc.ovlp.sorted.pruned, "beta"],
               bxse = abs(hQTL.tab.filt[snp.loc.ovlp.sorted.pruned, "beta"]/qnorm(hQTL.tab.filt[snp.loc.ovlp.sorted.pruned, "pval"]/2)),
               by = gwas.tab.all[snp.loc.ovlp.sorted.pruned, GWAS2beta[GWAS.NM]],
               byse = gwas.tab.all[snp.loc.ovlp.sorted.pruned, GWAS2se[GWAS.NM]]
               )
  rownames(b_se) = snp.loc.ovlp.sorted.pruned
  b_se.filt = b_se[apply(b_se, 1, FUN=function(x){sum(is.na(x))==0}), ]
  snp.loc.ovlp.sorted.pruned = rownames(b_se.filt)
    
  
  MRInputObject.cor <- mr_input(bx=b_se.filt[,"bx"],
                                  bxse=b_se.filt[,"bxse"],
                                  by=b_se.filt[,"by"],
                                  byse=b_se.filt[,"byse"],
                                  snps = snp.loc.ovlp.sorted.pruned,
                                  corr = ld.matrix[snp.loc.ovlp.sorted.pruned, snp.loc.ovlp.sorted.pruned]
                               )
  
  result = tryCatch(
  {
    mr_egger(MRInputObject.cor,
                      correl = TRUE,
                      distribution = "normal",
                      alpha = 0.05)
    
    # mr_ivw(MRInputObject.cor,
    #                   correl = TRUE,
    #                   distribution = "normal",
    #                   alpha = 0.05) 
  }, error = function(e) {
      print(e)
    return(NA)
  })
  
  if(is.na(result)) {return(NA)}
  print(result@Causal.pval)
  #return(result)
    
  #if(!any(abs(gwas.z[snp.loc.ovlp]) >= GWAS.Z.CUTOFF)) next
  
  #snp.topGWAS = snp.loc.ovlp[which.min(gwas.tab.all[snp.loc.ovlp.kept, "p"])[1]]
  
  if(!is.na(result@Causal.pval) && result@Causal.pval <=0.01)
  {
  
    r2.topGWAS = ld.matrix[snp.loc.ovlp.sorted.pruned[1], snp.loc.ovlp.sorted.pruned]
    
    df <- data.frame(x = b_se[, "bx"],
                     y = b_se[, "by"],
                     ymin = b_se[, "by"] - b_se[, "byse"],
                     ymax = b_se[, "by"] + b_se[, "byse"],
                     xmin = b_se[, "bx"] - b_se[, "bxse"],
                     xmax = b_se[, "bx"] + b_se[, "bxse"],
                     r2 = r2.topGWAS)
  
    p=ggplot(data = df,aes(x = x,y = y)) + 
        geom_point(aes(col= r2)) + 
        geom_errorbar(aes(ymin = ymin,ymax = ymax, col= r2), width=0) + 
        geom_errorbarh(aes(xmin = xmin,xmax = xmax, col= r2), height=0) +
        SC_COLOR +
        geom_hline(yintercept=0)+
        geom_vline(xintercept=0)+
        geom_abline(slope=result@Estimate, intercept = result@Intercept, linetype="dashed") +
        ggtitle(paste0(HM, " bipolar disorder ", hQTLPk, "\ncausal pval: ", signif(result@Causal.pval)))
  
    print(p)
  }

  
  
  
  return(list(qtl=b_se, MR=result))
})
names(hQTLPeaks.QTLvsGWAS) = hQTLPeaks.ovlpGWAS
dev.off()


save(hQTLPeaks.QTLvsGWAS, 
     file = file.ou.RData)
