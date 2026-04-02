#V1.1
#a.testBulkHaQTLVSGWAS_byColoc_V1.1_allGWAS_summary.R
#summary
#great


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
library(rGREAT)
#library("RColorBrewer")
# library(grid)
# library(gridExtra)

GREAT.FDR.CUTOFF=0.1
H4.cutoff=0.5
gARE.cutoff="emp.p.fdr.cutoff0.2"
HMs=c("H3K27ac", "H3K36me3", "H3K4me1", "H3K27me3", "H3K4me3")

dir.in = paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/")
files.in.RData=system(paste0("find ", dir.in, "*/*/*RData"), intern = T)

file.in.mAREs.hg19.RData="../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2HMPeak.hg19.RData"
file.in.mARE.modules.RData="../peakVariationAcrossTiss/a_4_modules_annotations_mergedPk5HMs_Epimap_CSandeGTExActivity/mARE2ModuleGrpName.RData"
files.in.ARE.hg38toHg19.RDS=paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".pk.hg38tohg19.RDS")
names(files.in.ARE.hg38toHg19.RDS) = HMs

files.in.gARE=paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData")
names(files.in.gARE) = HMs

dir.ou=paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS_summary/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.pdf =paste0(dir.ou, "bulkhQTLvsGWAS.coloc.H4", H4.cutoff, ".pdf")
file.ou.tsv =paste0(dir.ou, "bulkhQTLvsGWAS.coloc.H4", H4.cutoff, ".tsv")
file.ou.RData =paste0(dir.ou, "bulkhQTLvsGWAS.coloc.H4", H4.cutoff, ".RData")
file.ou.enrich.RData =paste0(dir.ou, "bulkhQTLvsGWAS.BPcoloc.enrich.H4", H4.cutoff, ".RData")
file.ou.enrich.pdf =paste0(dir.ou, "bulkhQTLvsGWAS.BPcoloc.enrich.H4", H4.cutoff, ".pdf")
#

#
load(file.in.mAREs.hg19.RData)

hm.hg38tohg19=lapply(files.in.ARE.hg38toHg19.RDS,
FUN=function(f.i)
{
  readRDS(f.i)
})

mAREs.gARE.hg19=lapply(HMs,
FUN=function(hm)
{
  load(files.in.gARE[hm])
  hQTLPk.hg38= hQTLPeaks.list[[gARE.cutoff]]
  hQTLPk.hg19 = hm.hg38tohg19[[hm]][hQTLPk.hg38]
  hm.peak2mergedPeak[[hm]][hQTLPk.hg19]
})
mAREs.gARE.hg19 = levels(as.factor(unlist(mAREs.gARE.hg19)))

write(names(mergedPeak2hmPeaks), file=paste0("../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2HMPeak.GREATbg.hg19.bed"))
write(mAREs.gARE.hg19, file=paste0("../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2HMPeak.gARE.GREATbg.hg19.bed"))
mAREs.bg=list()
mAREs.bg$all=names(mergedPeak2hmPeaks)
mAREs.bg$gAREs=mAREs.gARE.hg19
#
hm.gwas.coloc.hg38 = lapply(files.in.RData,
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
hm.gwas.coloc2HM= sapply(files.in.RData,
FUN=function(f.nm)
{
  hm=unlist(strsplit(f.nm, split="/", fixed = T))[3]
  hm=gsub("_100k", "", hm)
})
hm.gwas.coloc2GWAS= sapply(files.in.RData,
FUN=function(f.nm)
{
  unlist(strsplit(f.nm, split="/", fixed = T))[2]
})
hm2gwas2coloc.hg38.list = list()
for(i in 1:length(hm.gwas.coloc.hg38))
{
  hm2gwas2coloc.hg38.list[[hm.gwas.coloc2HM[[i]]]][[hm.gwas.coloc2GWAS[i]]]=hm.gwas.coloc.hg38[[i]]
}

#some haQTL peaks are not included in the final peaks
mARE.hm.gwas.coloc=lapply(1:length(hm.gwas.coloc.hg38),
FUN=function(i)
{
  h4s=hm.gwas.coloc.hg38[[i]]
  hm=hm.gwas.coloc2HM[i]
  pk.hg38= names(h4s)
  pk.hg19=hm.hg38tohg19[[hm]][pk.hg38]
  names(h4s)=hm.peak2mergedPeak[[hm]][pk.hg19]
  h4s.filt=h4s[!is.na(names(h4s))]
  
  h4s2H4 = split(h4s.filt, names(h4s.filt))
  h4s2maxH4= sapply(h4s2H4, max)
  return(h4s2maxH4)
})

mAREs.all = levels(factor(unlist(lapply(mARE.hm.gwas.coloc, names))))

mARE.hm_gwas.coloc.H4=matrix(NA, ncol=length(mARE.hm.gwas.coloc), nrow=length(mAREs.all))
colnames(mARE.hm_gwas.coloc.H4) = paste(hm.gwas.coloc2HM, hm.gwas.coloc2GWAS, sep=";")
rownames(mARE.hm_gwas.coloc.H4)=mAREs.all
for(i in 1:length(mARE.hm.gwas.coloc))
{
  h4s=mARE.hm.gwas.coloc[[i]]
  mARE.hm_gwas.coloc.H4[names(h4s), i] = h4s
}

mARE.hm_gwas.coloc.H4.filt =mARE.hm_gwas.coloc.H4[apply(mARE.hm_gwas.coloc.H4>H4.cutoff, 1, any, na.rm=T),]
#mARE.hm_gwas.coloc.H4.filt[is.na(mARE.hm_gwas.coloc.H4.filt)]=0


is.coloc.BP=apply(mARE.hm_gwas.coloc.H4.filt[, grepl("PGC.BP", colnames(mARE.hm_gwas.coloc.H4.filt))]>=H4.cutoff, 1 , any, na.rm=T)
mARE.hm_gwas.coloc.H4.filt.BP=mARE.hm_gwas.coloc.H4.filt[is.coloc.BP, ]
mARE.hm_gwas.coloc.H4.filt.BP[mARE.hm_gwas.coloc.H4.filt.BP<0.1]=0

gwas2cond=split(colnames(mARE.hm_gwas.coloc.H4), hm.gwas.coloc2GWAS)
mARE.gwas.colocMaxH4.aggregHMs=sapply(gwas2cond,
FUN=function(nms)
{
  apply(mARE.hm_gwas.coloc.H4[, nms], 1, max, na.rm=T)
})

mARE.gwas.coloc.shareProp= sapply(names(gwas2cond),
FUN=function(gwas.1)
{
  sapply(names(gwas2cond),
  FUN=function(gwas.2)
  {
    mARE.gwas.1 = rownames(mARE.gwas.colocMaxH4.aggregHMs)[mARE.gwas.colocMaxH4.aggregHMs[, gwas.1]>=H4.cutoff]
    mARE.gwas.2 = rownames(mARE.gwas.colocMaxH4.aggregHMs)[mARE.gwas.colocMaxH4.aggregHMs[, gwas.2]>=0.1]
    mARE.gwas.ovlp= intersect(mARE.gwas.1, mARE.gwas.2)
    return(length(mARE.gwas.ovlp)/length(mARE.gwas.1))
  })
})

save(hm.gwas.coloc.hg38,
     hm2gwas2coloc.hg38.list,
     hm.gwas.coloc2HM,
     hm.gwas.coloc2GWAS,
     mARE.hm_gwas.coloc.H4,
     mARE.hm_gwas.coloc.H4.filt,
     mARE.hm_gwas.coloc.H4.filt.BP,
     mARE.gwas.colocMaxH4.aggregHMs,
     mARE.gwas.coloc.shareProp,
     file=file.ou.RData)



# load(file.in.mARE.modules.RData)
# names(mARE2ModGrpName) = gsub(":|-", "\t", names(mARE2ModGrpName))
# 

pdf(file.ou.pdf, height=12, width=6)
# 
buf=mARE.hm_gwas.coloc.H4.filt
buf[is.na(buf)]=0
hclu=hclust(dist(buf, method = "euclidean"))
pheatmap(mARE.hm_gwas.coloc.H4.filt,
         color= colorRampPalette(c("white", "darkred"))(10),
         cluster_cols=F,
         cluster_row=hclu,
         show_rownames=F,
         na_col="grey",
         main="haQTL vs GWAS")


buf=mARE.hm_gwas.coloc.H4.filt.BP[, gwas2cond$PGC.BP]
buf[is.na(buf)]=0
hm.topColoc=colnames(buf)[(apply(buf, 1, which.max))]
names(hm.topColoc) = rownames(mARE.hm_gwas.coloc.H4.filt.BP)
clu2mAREs=split(names(hm.topColoc), hm.topColoc)
files.ou.clu=c()
for(hm in names(clu2mAREs))
{
  files.ou.clu[hm]=paste0(dir.ou, hm, ".topBPcoloc.mARE")
  write(clu2mAREs[[hm]], file=files.ou.clu[hm])
}
mARE.sorted= unlist(clu2mAREs)
# hclu=hclust(dist(buf, method = "euclidean"))
# mARE.BP.coloc.clu=as.character(cutree(hclu, k = 5))
# names(mARE.BP.coloc.clu)=rownames(mARE.hm_gwas.coloc.H4.filt.BP)
# other.GWAS.coloc=mARE.hm_gwas.coloc.H4.filt.BP[, setdiff(colnames(mARE.hm_gwas.coloc.H4.filt.BP), gwas2cond$PGC.BP)]
# other.GWAS.coloc.filt=other.GWAS.coloc[,apply(other.GWAS.coloc!=0, 2, any)]
pheatmap(mARE.hm_gwas.coloc.H4.filt.BP[mARE.sorted, gwas2cond$PGC.BP],
         color= colorRampPalette(c("white", "darkred"))(10),
         #annotation_row= data.frame(clu=mARE.BP.coloc.clu),
         annotation_row= data.frame(clu=hm.topColoc[mARE.sorted]),
         #annotation_row =as.data.frame(CREGrp=mARE2ModGrpName[rownames(mARE.hm_gwas.coloc.H4.filt.BP)]),
         #annotation_colors = "green",
         cluster_row=F,#hclu,
         show_rownames=F,
         na_col="grey",
         main="haQTL vs GWAS\ncoloc in BP")


# 
# 
# diag(mARE.gwas.coloc.shareProp) = max(c(mARE.gwas.coloc.shareProp[upper.tri(mARE.gwas.coloc.shareProp)],
#                                         mARE.gwas.coloc.shareProp[lower.tri(mARE.gwas.coloc.shareProp)]
#                                         ))+0.1
# pheatmap(mARE.gwas.coloc.shareProp,
#          cluster_rows=F,
#          cluster_cols=F,
#          main="proportion of shared coloc")

dev.off()

write.table(mARE.hm_gwas.coloc.H4.filt.BP[mARE.sorted, gwas2cond$PGC.BP], file=file.ou.tsv, sep="\t", quote=F)


#enrichment

pks2bed.df = function(pks)
{
  res= t(sapply(pks, FUN=function(pk) unlist(strsplit(pk, split="\t"))))
  data.frame(chr=res[,1],
             start=as.numeric(res[,2]),
             end=as.numeric(res[,3]),
             stringsAsFactors = F)
}

mAREs.df.bg = lapply(mAREs.bg, pks2bed.df)
clu.df.mARE= lapply(clu2mAREs, pks2bed.df)

enrich.res=lapply(mAREs.df.bg,
FUN=function(bg)
{
  lapply(clu.df.mARE,
  FUN=function(clu.mARE)
  {
    print(length(clu.mARE))
    job = submitGreatJob(clu.mARE, bg, species = "hg19")
    tb = getEnrichmentTables(job, category=c("GO", "Phenotype", "Genes"))
    lapply(tb,
    FUN=function(r)
    {
     r[r$Hyper_Adjp_BH<=GREAT.FDR.CUTOFF, c("name", "Hyper_Fold_Enrichment", "Hyper_Region_Set_Coverage", "Hyper_Raw_PValue")]
    })
    
  })
    
})

pdf(file.ou.enrich.pdf, width=12, height=15)
par(mar=c(3,22,4,4))
for(hm in names(clu.df.mARE))
{
  
  for(nm in names(enrich.res))
  {
   
    layout(matrix(1:8, ncol=2))
    res=enrich.res[[nm]][[hm]]
    for(type in names(res))
    {
      print(paste(hm, nm, type))
      if(nrow(res[[type]])==0) next
        
      term.logP=-log10(res[[type]]$Hyper_Raw_PValue)
      names(term.logP)=res[[type]]$name
      term.logP.sorted=sort(term.logP, decreasing = T)
      
      barplot(rev(term.logP.sorted[1:min(15, length(term.logP.sorted))]),
              col="red",
              horiz=T,
              las=2,
              main=paste(hm, nm, type))  
    }
    
  }
}


dev.off()

save(enrich.res,
     file=file.ou.enrich.RData)
