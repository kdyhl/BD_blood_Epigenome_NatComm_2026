#V1.1
#a.2.linkGenePeak_MR-Egger_V1.1.1_gARE_summary.R
#include enhancer 1MB around gene
#QTL around 100kb around enhancers
#filter hQTL based on multiple test for each enhancer
#further filter based on R2
#only check enhancers 10kb near eQTL hits

##use MR-Egger regression to prediction link between enhancer and gene 
#work on all genes 
#data loading  part is the same from PRS code
  
if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



args = commandArgs(trailingOnly=TRUE)
TISS = args[1] #Brain Heart Muscle Lung
HM= args[2]
hQTL.WIND = as.numeric(args[3])

                

# library(data.table)
# library(pheatmap)
# library("RColorBrewer")
# library(grid)
# library(gridExtra)
# library(ggplot2)
library(MendelianRandomization)
#library(coloc)
#library(ggforce)
#library(segclust2d)
#library("RSpectra")
library(PRROC)


tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
                  Lung = "Lung",
                  Muscle = "Muscle_Skeletal",
                  Heart = "Heart_Left_Ventricle",
                  Blood = "Whole_Blood")

# COND2HiCTISS = c(Brain = "Dorsolateral_Prefrontal_Cortex", #Hippocampus|
#                   Lung = "Lung",
#                   #Muscle = "Muscle_Skeletal",
#                   Heart = "Left_Ventricle")


files.in.MR.RData = system(paste0("find ~/hptmp/MayoBipolar/GeneEnhLink_MR-egger_V1.1.1/", TISS, "_", HM, "_", hQTL.WIND/1000, "k/", TISS, "*geneEnhLink_*.RData"), intern = T)
#file.in.eQTL2pk.RData=paste0("~/hptmp/eGTEx-H3K27ac/GeneEnhLink_MR-egger_V1.1/", TISS, "/", TISS, ".gene2peaks.hg38.RData")
file.in.g2pk.RData = paste0("./a_gAREsNeareGene_emp.p.fdr.cutoff0.2/", TISS, "_1000k/", HM, ".eGene2gARE.hg38.RData")
file.in.eGene.info = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.egenes.txt.gz")

#file.in.g2pk.RData = paste0(dir.tmp, TISS, ".gene2peaks.hg38.RData")

#file.in.PromHiC.hg38.bed= "~/lhou.compbio/data/HiC/PromCapture_HiC_GSE86189_hg19_2019/Promoter-other.signif.allTissues.hg38.bed"
# file.in.peak.hg382hg19.bed = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", TISS, ".pks.hg382hg19Uniq.RDS")

# 
# dir.tmp= paste0("~/hptmp/eGTEx-H3K27ac/GeneEnhLink_MR-egger_V1.1/", TISS, "/")
# file.tmp.eQTL.pks.hg38.bed = paste0(dir.tmp, TISS, ".eQTL.pks.bed")

dir.ou="./a_2_linkGenePeak_MR_Egger_V1.1.1_summary/"
dir.create(dir.ou, recursive = T, showWarnings = F)
file.ou.RData=paste0(dir.ou, TISS, "_", HM, "_", hQTL.WIND/1000, "k.gene2peaks_MR.RData")
file.ou.pdf=paste0(dir.ou, TISS, "_", HM, "_", hQTL.WIND/1000, "k.gene2peaks_MR.pdf")
#
gene2pks.MR = list()
for(f.i.RData in files.in.MR.RData)
{
  load(f.i.RData)  
  gene2pks.MR=c(gene2pks.MR, gene2enh.MR_Egger)
  
}

gene2pks.MR.pvals = lapply(gene2pks.MR,
FUN=function(pk2MR)
{
  pvals = lapply(pk2MR,
  FUN=function(res)
  {
    if(is.na(res))  
    {
      return(NA)  
    }else
    {
      return(res@Pvalue.Est)
    }
  })
  
  #pvals[is.na(pvals)]=1
  #qvals=rep(1, length(pvals))
  
  return(pvals)#[!is.na(pvals)]
  
  
})
#gene2pks.MR.pvals= gene2pks.MR.pvals[sapply(gene2pks.MR.pvals, length)!=0]

gene2pks.MR.qvals = lapply(gene2pks.MR.pvals,
FUN=function(pvals)
{
  #pvals
  qvals = p.adjust(pvals, method="BH")
  qvals[is.na(qvals)]=1
  #qvals=rep(1, length(pvals))
  
  return(qvals)
  
  
})


# 
# #overlapping with HiC
# eQTL.pks = levels(factor(unlist(lapply(gene2pks.MR, names))))
# write(paste(gsub(":|-", "\t", eQTL.pks), eQTL.pks, sep="\t"),
#       file= file.tmp.eQTL.pks.hg38.bed)
# 
# 
# cmd = paste0("intersectBed -wa -wb -a ", file.tmp.eQTL.pks.hg38.bed, " -b ", file.in.PromHiC.hg38.bed)
# buf =system(cmd, intern = T)
# pksAndGeneAndTiss.byHiC = t(sapply(buf,
# FUN=function(b)
# {
#   buff = unlist(strsplit(b, split="\t"))
#   pk = buff[4] #paste0(buff[1], ":", buff[2], "-", buff[3])
#   gene = buff[8]
#   tiss = buff[9]
#   
#   return(c(pk, gene, tiss))
# }))
# colnames(pksAndGeneAndTiss.byHiC) = c("ARE", "gene", "tissue")
# 
# pksAndGeneAndTiss.byHiC.tiss = pksAndGeneAndTiss.byHiC[pksAndGeneAndTiss.byHiC[,"tissue"]==COND2HiCTISS[TISS],1:2]
# tissHiC.g2pk = split(pksAndGeneAndTiss.byHiC.tiss[,"ARE"], 
#                   pksAndGeneAndTiss.byHiC.tiss[,"gene"])
# 
# gene.ENSG2name=as.matrix(read.table(file.in.eGene.info, sep="\t", header = T, row.names = 1, stringsAsFactors = F))[,1]

# 
# #scores
# load(file.in.eQTL2pk.RData)
# 
# 
# gene2Pk.dis = lapply(names(gene2pks.MR),
# FUN=function(g)
# {
#   pks = names(gene2pks.MR[[g]])
#   pks.pos = sapply(pks,
#   FUN=function(pk)
#   {
#     buf = unlist(strsplit(pk, split=":|-"))  
#     mean(as.numeric(buf[2:3]))
#   })
#   
#   tss= gene2TSS[g, "gene.TSS"]
#   dis=abs(pks.pos-tss)  
#   names(dis)= pks
#   
#   return(dis)
# })
# names(gene2Pk.dis) = names(gene2pks.MR)
# 
# 
# genes.ovlp=names(gene2pks.MR.qvals)[gene.ENSG2name[names(gene2pks.MR.qvals)] %in% names(tissHiC.g2pk)]
# 
# gAndpeakPair.scores=list()
# gAndpeakPair.scores$MR.qval = -unlist(gene2pks.MR.qvals[genes.ovlp])
# gAndpeakPair.scores$MR.pval = -unlist(gene2pks.MR.pvals[genes.ovlp])
# gAndpeakPair.scores$MR.pval[is.na(gAndpeakPair.scores$MR.pval)]=-1
# gAndpeakPair.scores$distance = -unlist(gene2Pk.dis[genes.ovlp])
# 
# 
# gAndpeakPair.class = unlist(lapply(genes.ovlp,
# FUN=function(g)
# {
#   pks = names(gene2pks.MR[[g]])
#   classLab=rep("bg", length(pks))
#   names(classLab)=pks
#   
#   classLab[intersect(pks, tissHiC.g2pk[[gene.ENSG2name[g]]])] = "fg"
#   return(classLab)
# }))
# 
# 

# #visualize
# PRs.res<-lapply(gAndpeakPair.scores,
# FUN=function(scores)
# {
#   pr.curve(scores.class0 = scores[gAndpeakPair.class=="fg" & gAndpeakPair.scores$MR.pval != -1], 
#            scores.class1 = scores[gAndpeakPair.class=="bg" & gAndpeakPair.scores$MR.pval != -1], 
#            curve=T)
# })
# 
# 
pdf(file.ou.pdf, height=8, width = 8)
hist(unlist(gene2pks.MR.pvals), xlab="pvalue")
# plot(c(0,1), c(0,1), 
#      xlab="recall", 
#      ylab= "precision",
#      type="n",
#      main=TISS)
# cols = rainbow(length(PRs.res))
# names(cols) = names(PRs.res)
# for(nm  in names(PRs.res))
# {
#   lines(PRs.res[[nm]]$curve[,1:2], col=cols[nm], lwd=2)
# }
# legend("topright", 
#        legend=paste(names(PRs.res), "AUPRC:", signif(sapply(PRs.res, FUN=function(x) x$auc.integral),3)),
#        lwd=2,
#        col= cols
#        )
dev.off()



save(gene2pks.MR,
     gene2pks.MR.pvals,
     gene2pks.MR.qvals,
     # genes.ovlp,
     # pksAndGeneAndTiss.byHiC,
     # tissHiC.g2pk,
     # gAndpeakPair.scores,
     # PRs.res,
     file=file.ou.RData)
