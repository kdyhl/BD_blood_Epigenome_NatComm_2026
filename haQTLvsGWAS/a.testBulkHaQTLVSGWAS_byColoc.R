#V1.1
#a.testBulkHaQTLVSGWAS_byColoc_V1.1_allGWAS.R
#all other GWAS

#shared LD structure since GWAS is imputed based on GTEx EUR cohort

args = commandArgs(trailingOnly=TRUE)
HM = args[1] #H3K27ac, H3K4me1 ...
GWAS.NM=args[2]
# #print(HM)

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
library(ggplot2)
library(RColorBrewer)
MY.PALETTE <- colorRampPalette(c("royalblue2", "green1", "orange2", "darkred")) #rev(brewer.pal(11, "Spectral"))
SC_COLOR = scale_colour_gradientn(colours = MY.PALETTE(8), limits=c(0,1)) 


library(coloc)
source("~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R")




PK.WIND = 100000
GWAS.P.CUTOFF = 1e-6
GWAS.Z.CUTOFF = abs(qnorm(GWAS.P.CUTOFF/2))
hQTL.EMP.FDR.CUTOFF=0.2


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

file.in.hQTLPeak.RData = paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, ".haQTLPeak.cutoffs.hg38.RData")

dir.in.Mayo.genotype.plink = "~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/"
file.in.Mayo.genotype.vcf="~/lhou.compbio/data/Mayo_Bipolar/WGS/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID.vcf.gz"



#dir.in.1KG.genotype.plink ="/broad/compbio/data/1KG_phase3/eur/"

dir.tmp=paste0("~/hptmp/MayoBipolar/a_bulkhQTLvsGWAS_bycoloc_V1.1/", GWAS.NM, "/", HM, "/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.hQTLPk.bed = paste0(dir.tmp, HM, ".hQTLPk.hg38.bed")
file.tmp.GWAS.hits.bed = paste0(dir.tmp, "GWAS.hits.hg38.bed")
file.tmp.qtlPk.ovlp.GWAS.hits.bed = paste0(dir.tmp, HM, ".QTLPK.ovlp.GWAS.hits.hg38.bed")

#file.tmp.eQTL.chr.RDS = paste0(dir.tmp, HM, ".", CHR, ".eQTL.RDS")
#file.tmp.hQTL.chr.RDS = paste0(dir.tmp, HM, ".", CHR, ".hQTL.RDS")

#file.tmp.peak.bed = paste0(dir.tmp, HM, ".", CHR, ".pk.hg38.bed")
#file.tmp.gene.bed = paste0(dir.tmp, HM, ".", CHR,".eGene.hg38.bed")

# file.tmp.hQTL.hg38.bed = paste0(dir.tmp, HM, ".", CHR,".hQTL.hg38.bed")
# file.tmp.hQTL.unmap.bed = paste0(dir.tmp, HM, ".", CHR,".hQTL.unmap.bed")

dir.ou = paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/", GWAS.NM, "/", HM, "_100k/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.pdf =paste0(dir.ou, HM, ".bulkhQTLvsGWAS.pdf")
file.ou.RData =paste0(dir.ou, HM, ".bulkhQTLvsGWAS.RData")
#


#focus on those peak overlapped with strong signals in GWAS
load(file.in.hQTLPeak.RData)
hQTLPk.hg38= hQTLPeaks.list[[paste0("emp.p.fdr.cutoff", hQTL.EMP.FDR.CUTOFF)]]#rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.P"] <= hQTL.EMP.PVAL]
hQTLPk.hg38.bed = paste(gsub(":|-", "\t", hQTLPk.hg38), hQTLPk.hg38, sep="\t")
write(hQTLPk.hg38.bed, file.tmp.hQTLPk.bed)

cmd = paste0("zcat ", files.in.GWAS.gz[GWAS.NM], "|awk '($", GWAS2P.ID[GWAS.NM], "<", GWAS.P.CUTOFF, "){print }'>", file.tmp.GWAS.hits.bed)
system(cmd)

cmd = paste0("bedtools window -w ", PK.WIND, " -a ", file.tmp.hQTLPk.bed, " -b ", file.tmp.GWAS.hits.bed, ">", file.tmp.qtlPk.ovlp.GWAS.hits.bed)
system(cmd)
buf = read.table(file.tmp.qtlPk.ovlp.GWAS.hits.bed, "\t", header=F, row.names = NULL)
hQTLPeaks.ovlpGWAS = levels(factor(buf[,4])) #paste0(buf[,1], ":", buf[,2], "-", buf[,3])))


pdf(file.ou.pdf, width=12, height=8)
#layout(matrix(1:9, ncol=3))
hQTLPeaks.QTLvsGWAS=lapply(hQTLPeaks.ovlpGWAS,
FUN=function(hQTLPk)
{
  print(hQTLPk)
  f.tmp.ovlp.snps =paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".", PK.WIND/1000, "k.ovlpSNP.txt")
  f.tmp.ld.pre = paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".ovlpSNP")
  
  f.tmp.hQTL.z = paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".hQTL.z.txt")
  f.tmp.gwas.z = paste0(dir.tmp, HM, ".", gsub(":", "-", hQTLPk), ".gwas.z.txt")
  
  #f.o.eCaviar = paste0(d.o.gwas, GWAS.nm, ".", HM, ".", gsub(":", "-", hQTLPk), ".txt")
  
  #SNPs from hQTL
  buf = unlist(strsplit(hQTLPk, ":|-", perl=T))
  pk.chr= buf[1]
  pk.start = as.numeric(buf[2])
  pk.end = as.numeric(buf[3])
  query = paste0(pk.chr, ":", max(1, pk.start-PK.WIND), "-", pk.end+PK.WIND)
  hQTL.tab = seqminer::tabix.read.table(files.in.hQTL.bed[pk.chr], query) 
  hQTL.tab.filt = hQTL.tab[hQTL.tab$peak==hQTLPk,]
  rownames(hQTL.tab.filt) = hQTL.tab.filt$rs
  hQTL.tab.filt$varbeta=(-abs(hQTL.tab.filt[, "beta"])/qt(hQTL.tab.filt[, "pval"]/2, df=HM2sampN[HM]-2-HM2COVNum[HM]))^2

  
  #SNPs from GWAS
  gwas.tab = seqminer::tabix.read.table(files.in.GWAS.gz[GWAS.NM], query) 
  rownames(gwas.tab) = paste(gwas.tab[,1], gwas.tab[,3], gwas.tab[,5], gwas.tab[,6], sep="_")
  gwas.tab.flip = gwas.tab
  gwas.tab.flip[,5]= gwas.tab[,6]
  gwas.tab.flip[,6]= gwas.tab[,5]
  rownames(gwas.tab.flip)  = paste(gwas.tab.flip[,1], gwas.tab.flip[,3], gwas.tab.flip[,5], gwas.tab.flip[,6], sep="_")
  gwas.tab.all =rbind(gwas.tab, gwas.tab.flip)
  
  snp.loc.ovlp = intersect(rownames(hQTL.tab.filt), rownames(gwas.tab.all))
  
  
  #SNP for MAF
  snp.tab = seqminer::tabix.read.table(file.in.Mayo.genotype.vcf, query)
  snp2maf = sapply(snp.tab$INFO,
  FUN=function(x)
  {
    as.numeric(gsub(".*AF=(0.\\d+);.*", "\\1", x, perl=T))
  })
  names(snp2maf) = snp.tab$ID
   
  
  dataset.gwas= list(beta=gwas.tab.all[snp.loc.ovlp, GWAS2beta[GWAS.NM]], 
                        varbeta=gwas.tab.all[snp.loc.ovlp, GWAS2se[GWAS.NM]]^2,
                        N=GWAS2N[GWAS.NM], type=GWAS2type[GWAS.NM])
  if(GWAS2type[GWAS.NM]=="cc")
  {
    dataset.gwas$s=GWAS2CaseProp[GWAS.NM]  
  }
  
  
  dataset.hQTL = list(beta=hQTL.tab.filt[snp.loc.ovlp, "beta"], 
                        varbeta=hQTL.tab.filt[snp.loc.ovlp, "varbeta"],
                        N=HM2sampN[HM], type="quant")
    
  my.res <- coloc.abf(dataset1=dataset.gwas,
                  dataset2=dataset.hQTL,
                  MAF=snp2maf[snp.loc.ovlp])#,
    
  
  
  #overlaped SNPs and flapping
  if(length(snp.loc.ovlp)!=0 && 
     (!is.na(my.res$summary["PP.H4.abf"])) &&
     my.res$summary["PP.H4.abf"]>=0.5)
  {
    print(paste0(length(snp.loc.ovlp), " of ", nrow(hQTL.tab.filt), " hQTL is with GWAS zscore"))
    write(snp.loc.ovlp, f.tmp.ovlp.snps)
  
    # plot(hQTL.z[snps.locID.ovlp],
    #     gwas.z[snps.locID.ovlp], 
    #     xlab="hQTL zscore",
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
    snp.loc.ovlp.loc =sapply(snp.loc.ovlp, FUN=function(x) as.numeric(unlist(strsplit(x, "_"))[2]))
  
    if(length(snp.loc.ovlp)==0) next
    #if(!any(abs(gwas.z[snp.loc.ovlp]) >= GWAS.Z.CUTOFF)) next
    
    snp.topGWAS = snp.loc.ovlp[which.min(gwas.tab.all[snp.loc.ovlp.kept, GWAS2P[GWAS.NM]])[1]]
    
    df =data.frame(log10p= c(-log10(gwas.tab.all[snp.loc.ovlp.kept, GWAS2P[GWAS.NM]]),
                             -log10(hQTL.tab.filt[snp.loc.ovlp.kept, "pval"])),
                   type= c(rep("GWAS", length(snp.loc.ovlp.kept)),
                           rep("hQTL", length(snp.loc.ovlp.kept))),
                   r2.withTopGWASHit = rep(ld.matrix[snp.topGWAS,], 2),
                   loc = rep(snp.loc.ovlp.loc[snp.loc.ovlp.kept], 2))
    p=ggplot(df, aes(x=loc, y=log10p)) +
      geom_point(aes(col=r2.withTopGWASHit)) +
      SC_COLOR + #scale_colour_gradientn(colours = cols) +
      ggtitle(paste0(HM, " bipolar disorder ", hQTLPk, "\nH4 ABF: ", signif(my.res$summary["PP.H4.abf"],3))) +
      facet_grid(type ~ .,scales="free_y") +
      theme_bw()
    print(p)
  }
  # if(length(snp.loc.ovlp)==0) next
  # write(paste(snp.loc.ovlp, hQTL.z[snp.loc.ovlp], sep="\t"), file=f.tmp.hQTL.z)
  # write(paste(snp.loc.ovlp, gwas.z[snp.loc.ovlp], sep="\t"), file=f.tmp.gwas.z)


  
  # system("gunzip ", paste0(f.tmp.ld, ".ld.gz") )
  # 
  # #run eCaviar
  # 
  # cmd = paste0(eCAVIAR, " -o ", f.o.eCaviar,
  #              " -l ", paste0(f.tmp.ld.pre, ".ld"),
  #              " -l ", paste0(f.tmp.ld.pre, ".ld"),
  #              " -z ", f.tmp.gwas.z,
  #              " -z ", f.tmp.hQTL.z,
  #              #" -r ", 0.01,
  #              " -c 5"
  #              )
  # 
  # write(cmd, file.ou.sh, append = T)
  return(list(qtl=df, coloc=my.res))
})
names(hQTLPeaks.QTLvsGWAS) = hQTLPeaks.ovlpGWAS
dev.off()


save(hQTLPeaks.QTLvsGWAS, 
     file = file.ou.RData)
