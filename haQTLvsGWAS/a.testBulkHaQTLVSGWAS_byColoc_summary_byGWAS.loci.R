#V1.1
#a.testBulkHaQTLVSGWAS_byColoc_V1.1_allGWAS_summary_byGWAS.loci.R
#summary
#per GWAS loci rather than mCREs


#shared LD structure since GWAS is imputed based on GTEx EUR cohort



if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))




options(scipen = 999)

#library(data.table)
library(pheatmap)
#library(rGREAT)
#library("RColorBrewer")
# library(grid)
# library(gridExtra)

COLOC.H4.CUTOFF=0.5
COLOC.WIND=100000
gARE.cutoff="emp.p.fdr.cutoff0.2"
HMs=c("H3K27ac", "H3K36me3", "H3K4me1", "H3K27me3", "H3K4me3")

dir.in = paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/PGC.BP/")
files.in.coloc.RData=system(paste0("find ", dir.in, "*/*RData"), intern = T)

file.in.BD.GWAS.clump.bed="a_GWAS_clump/BD_PGC_bip/BDlooMAYO_PGC3_hg38.P0.000001.R20.2.dis250kb.hg38.bed"

files.in.gARE=paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData")
names(files.in.gARE) = HMs

dir.tmp="~/hptmp/MayoBipolar/hQTLvsGWAS/coloc.summary.byLocus/"
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.hm.gAREs.bed=paste0(dir.tmp, "hm.colcoc.gAREs.bed")
file.tmp.hm.gAREsOvlpGWAS.bed=paste0(dir.tmp, "hm.colcoc.gAREs.ovlp.GWAS.bed")

dir.ou=paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS_summary_byGWAS.loci/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.tsv =paste0(dir.ou, "bulkhQTLvsGWAS.coloc.tsv")
file.ou.RData =paste0(dir.ou, "bulkhQTLvsGWAS.coloc.RData")
#

#
hm.coloc.hg38 = lapply(files.in.coloc.RData,
FUN=function(f.i.RData)
{
  print(f.i.RData)
  load(f.i.RData)
  H4s=sapply(hQTLPeaks.QTLvsGWAS,
  FUN=function(res)
  {
    res$coloc$summary["PP.H4.abf"]
  })
  names(H4s)= gsub(".PP.H4.abf", "", names(H4s))
  return(H4s)
})
names(hm.coloc.hg38)= gsub(".bulkhQTLvsGWAS.RData", "", basename(files.in.coloc.RData))


#
BD.GWAS.clump.df=read.table(file.in.BD.GWAS.clump.bed, head=T, row.names="SNP.loc_ID", sep="\t", comment.char="", stringsAsFactors=F)

#
GWAS.clump2gARE.coloc.list=lapply(names(hm.coloc.hg38),
FUN=function(hm)
{
  print(hm)
  
  load(files.in.gARE[hm])
  gARE2hQTLs=hQTLpeak2signifQTL.bycutoffs[[gARE.cutoff]]

  gAREs.coloc=names(hm.coloc.hg38[[hm]])[hm.coloc.hg38[[hm]]>=COLOC.H4.CUTOFF]
  gAREs.coloc=gAREs.coloc[!is.na(gAREs.coloc)]
  gAREs.coloc.bed=gsub(":|-", "\t", gAREs.coloc)

  write(paste(gAREs.coloc.bed, gAREs.coloc, sep="\t"), file.tmp.hm.gAREs.bed)

  cmd=paste0("bedtools window -w ", COLOC.WIND, " -a ", file.tmp.hm.gAREs.bed, " -b ", file.in.BD.GWAS.clump.bed, " > ", file.tmp.hm.gAREsOvlpGWAS.bed)
  system(cmd)

  d=read.table(file.tmp.hm.gAREsOvlpGWAS.bed, head=F, row.names=NULL, sep="\t", stringsAsFactors=F)
  d.GWASishQTL=sapply(1:nrow(d),
  FUN=function(i)
  {
    gARE=d[i, 4]
    snp=d[i, 8]
    snp %in% gARE2hQTLs[[gARE]][,2]
  })
  d.filt=d[d.GWASishQTL, ]

  #GWAS.leadSNP2gARE=split(d.filt[,4], d.filt[,10])
  df=data.frame(lead.SNP=d.filt[,10],
                gARE=d.filt[,4],
                histoneMark=hm,
                coloc.H4=hm.coloc.hg38[[hm]][d.filt[,4]],
                stringsAsFactors=F)
  df.filt=df[!duplicated(paste(df$lead.SNP, df$gARE, sep=";")), ]

})
GWAS.clump2gARE.coloc.df=do.call(rbind, GWAS.clump2gARE.coloc.list)


GWAS.clump2histone.counts.matrix = as.matrix(with(GWAS.clump2gARE.coloc.df, table(lead.SNP, histoneMark)))

GWAS.clump2histone.counts.matrix=cbind(GWAS.clump2histone.counts.matrix, 
                                        pval=-log10(BD.GWAS.clump.df[rownames(GWAS.clump2histone.counts.matrix),"P"]))

save(hm.coloc.hg38,
    GWAS.clump2gARE.coloc.list,
    GWAS.clump2histone.counts.matrix,

    file=file.ou.RData)


write.table(GWAS.clump2histone.counts.matrix, file=file.ou.tsv, sep="\t", quote=F)

# #
# d=read.table(file.in.BD.GWAS.clump.bed, sep="\t", header=T, row.names=NULL, comment.char="", stringsAsFactors=F)


# # load(file.in.mARE.modules.RData)
# # names(mARE2ModGrpName) = gsub(":|-", "\t", names(mARE2ModGrpName))
# # 

# pdf(file.ou.pdf, height=12, width=6)
# # 

# dev.off()




