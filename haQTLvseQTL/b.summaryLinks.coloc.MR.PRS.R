#V1.3
#b.summaryLinks.coloc.MR.PRS_V1.3_gARE.R

#compare link prediction from different resources by cutoffs
#adding fm-eQTL too

args = commandArgs(trailingOnly=TRUE)
TISS = args[1]
HM= args[2]

tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
                  Lung = "Lung",
                  Muscle = "Muscle_Skeletal",
                  Heart = "Heart_Left_Ventricle",
                  Blood = "Whole_Blood")

SCORE2COL = c(FMeQTL2ARE.DistanceInv = "turquoise3",
              #eQTL2ARE.minP ="steelblue",
              coloc.H4.ABF="red2", 
              coloc.H4vsH3.ratio="orangered3", 
              MR.qval="blueviolet",  
              #PGS.pcor.qval="darkorange1",#, #"magenta3",
              #scores.comb = "magenta3",
              gene2ARE.DistanceInv="gold1"
              # EpiMap="khaki3",
              # ABCscore="goldenrod3"
              
)
# COND2HiCTISS = c(Brain = "Dorsolateral_Prefrontal_Cortex", #Hippocampus|
#                   Lung = "Lung",
#                   #Muscle = "Muscle_Skeletal",
#                   Heart = "Left_Ventricle")

# #GWAS.PVAL.CUTOFF = 1e-6
coloc.H4.ABF.CUTOFF.forRatio=0.1
# #GSP.TOP.Perc= 0.01
# coloc.H4.ABF.CUTOFF=0.5
# coloc.H4_H3.Ratio.ABF.CUTOFF=5
# MR.Q.CUTOFF = -0.05
# #PGS.pcor.log10PVAL.CUTOFF=-log10(0.05)
# FMeQTL.DistInV.CUTOFF=1/2000


GSP.cutoffs = c("coloc.H4.ABF" =0.5,
                  "coloc.H4vsH3.ratio" =5, 
                  "MR.qval" =-0.05, 
                  "FMeQTL2ARE.DistanceInv" = 1/2000)
PROMT.WIND = 2000

#library(Sushi)
#library(igraph)
library(PRROC)

file.in.eGene.info = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.egenes.txt.gz")
file.in.g2pk.RData = paste0("./a_gAREsNeareGene_emp.p.fdr.cutoff0.2/", TISS, "_1000k/", HM, ".eGene2gARE.hg38.RData")

#file.in.ARE.modules.RData = "../variationAcrossTissueandIndiv/a_3_AREModulesInEpimap_V1.4_mergedPk4Tiss/4TissMergedPeaks.inEpimap.H3K27acSignal.RData"
#file.in.mARE.RData = "../peakMergeNormal/c_mergePeaksV1.3FromTissues_hg19/4Tiss.mergedPeak2tissuePeak.hg19.RData"
#file.in.hg38to2hg19.RDS= paste0("../peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", TISS, ".pks.hg382hg19Uniq.RDS")
file.in.hg38to2hg19.RDS= paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".pk.hg38tohg19.RDS") 


file.in.g2pks.coloc.RData = paste0("./a_2_coloc_haQTLvseQTL_bycoloc_V1.3_summary/", TISS, "_", HM, "_100k.haQTLvseQTL.RData")
file.in.g2pks.MR.RData = paste0("./a_2_linkGenePeak_MR_Egger_V1.1.1_summary/", TISS, "_", HM, "_100k.gene2peaks_MR.RData")
#file.in.g2pks.PRC.pcor.RData = paste0("./a_PRSExpvsPeak_V1.3_corWithCovar/", TISS, "/", TISS, ".gene2Peak.10keQTL.pCor.RData")
file.in.g2pks.FmeQTL.RData = paste0("./a_genePeakLink_byFMeQTL_V1.1_gARE/", TISS, "_10k/", TISS, ".", HM, ".gen2peaks.byFMeQTL.RData")
#file.in.PromHiC.hg38.bed= "~/lhou.compbio/data/HiC/PromCapture_HiC_GSE86189_hg19_2019/Promoter-other.signif.allTissues.hg38.bed"
# file.in.peak.hg382hg19.bed = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", TISS, ".pks.hg382hg19Uniq.RDS")



# dir.tmp= paste0("~/hptmp/eGTEx-H3K27ac/b_cmp_links/", TISS, "/")
# dir.create(dir.tmp, showWarnings = F, recursive = T)
# file.tmp.eQTL.pks.hg38.bed = paste0(dir.tmp, TISS, ".eQTL.pks.bed")
# file.tmp.TSS.hg38.bed = paste0(dir.tmp, TISS, ".TSS.hg38.bed")
# file.tmp.AREs.hg38.bed = paste0(dir.tmp, TISS, ".ARE.hg38.bed")
# file.tmp.AREsOVLPTSS.hg38.bed = paste0(dir.tmp, TISS, ".AREOVLPTSS.hg38.bed")

dir.ou="b_summary_Links.coloc.MR.FMeQTL_V1.3_gARE/"
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.scores.RData = paste(dir.ou, TISS, "_", HM, ".linkFromColoc.MR.RData", sep="")
file.ou.cmp.pdf = paste(dir.ou, TISS, "_", HM, ".cmp.linkFromColoc.MR.noPromoter.pdf", sep="")
#file.ou.haQTLDiffPks2GWAS.RData = paste(file.out.pref, ".haQTLDiffPks2GWAS.RData", sep="")
#file.ou.geneNeighb = paste(file.out.pref, ".txt", sep="")


######
buf=as.matrix(read.table(file.in.eGene.info, sep="\t", header = T, row.names = NULL, stringsAsFactors = F))[,1:2]
gene.gname2ESG = split(buf[,1], buf[,2])
gene.gname2ESG = gene.gname2ESG[sapply(gene.gname2ESG, length)==1]
gene.gname2ESG = unlist(gene.gname2ESG)
#coloc link ########################################################################
load(file.in.g2pks.coloc.RData)
geneAndPeak.coloc= lapply(names(gene2peak.H4ABF),
FUN=function(g)
{
  pks = names(gene2peak.H4ABF[[g]])
  data.frame(gene=g, peak= pks, H4.ABF=gene2peak.H4ABF[[g]], H4vsH3.ratio = gene2peak.H4vsH3Ratio[[g]][pks] , stringsAsFactors = F)
})
geneAndPeak.coloc.df = do.call(rbind, geneAndPeak.coloc)
rownames(geneAndPeak.coloc.df) = paste(geneAndPeak.coloc.df[,1], geneAndPeak.coloc.df[,2], sep=";")

# geneAndPeak.coloc.filt.df = geneAndPeak.coloc.df[(!is.na(geneAndPeak.coloc.df$H4vsH3.ABF)) & geneAndPeak.coloc.df$H4vsH3.ABF>=coloc.H4vsH3ABF.CUTOFF,]
# rownames(geneAndPeak.coloc.filt.df) = paste(geneAndPeak.coloc.filt.df[,1], geneAndPeak.coloc.filt.df[,2], sep=";")
#MR link ########################################################################
load(file.in.g2pks.MR.RData)

geneAndPeak.MR_Pval= lapply(names(gene2pks.MR.pvals),
FUN=function(g)
{
  data.frame(gene=g, peak= names(gene2pks.MR.pvals[[g]]), MR.pval=unlist(gene2pks.MR.pvals[[g]]), stringsAsFactors = F)
})
geneAndPeak.MR_P.df = do.call(rbind, geneAndPeak.MR_Pval)
rownames(geneAndPeak.MR_P.df) = paste(geneAndPeak.MR_P.df[,1], geneAndPeak.MR_P.df[,2], sep=";")
# geneAndPeak.MR.filt.df = geneAndPeak.MR_Q.df[geneAndPeak.MR_Q.df$MR.qval<=MR.QVAL.CUTOFF,]
# rownames(geneAndPeak.MR.filt.df) = paste(geneAndPeak.MR.filt.df[,1], geneAndPeak.MR.filt.df[,2], sep=";")

# #PRC pcor link ########################################################################
# load(file.in.g2pks.PRC.pcor.RData)
# geneAndPeak.PRS_pcor.Pval= lapply(names(g2peaks.PRS.pcor.pval),
# FUN=function(g)
# {
#   #print(g)
#   data.frame(gene=g, peak= names(g2peaks.PRS.pcor.pval[[g]]), PGS.pcor.pval=g2peaks.PRS.pcor.pval[[g]], stringsAsFactors = F)
# })
# geneAndPeak.PRS_pcor.Pval.df = do.call(rbind, geneAndPeak.PRS_pcor.Pval)
# rownames(geneAndPeak.PRS_pcor.Pval.df) = paste(geneAndPeak.PRS_pcor.Pval.df[,1], geneAndPeak.PRS_pcor.Pval.df[,2], sep=";")


#FMeQTL
load(file.in.g2pks.FmeQTL.RData)


#distance ########################################################################
load(file.in.g2pk.RData)
gene2Pk.dist = lapply(names(eGene2gAREs.hg38),
FUN=function(g)
{
  #print(g)
  pks = eGene2gAREs.hg38[[g]]
  pks.pos = sapply(pks,
  FUN=function(pk)
  {
    buf = unlist(strsplit(pk, split=":|-"))  
    mean(as.numeric(buf[2:3]))
  })
  
  tss= eGene2TSS[g, "gene.tss"]
  dis=abs(pks.pos-tss)  
  names(dis)= pks
  
  data.frame(gene=g, peak= pks, distance=dis, stringsAsFactors = F)
})
gene2Pk.dis.df = do.call(rbind, gene2Pk.dist)
rownames(gene2Pk.dis.df) = paste(gene2Pk.dis.df[,1], gene2Pk.dis.df[,2], sep=";")

#

scores.all = matrix(NA, ncol=5, nrow=nrow(gene2Pk.dis.df))
rownames(scores.all) = rownames(gene2Pk.dis.df)
colnames(scores.all) = c("coloc.H4.ABF", "coloc.H4vsH3.ratio", "MR.qval", 
                         "FMeQTL2ARE.DistanceInv", "gene2ARE.DistanceInv")

scores.all[rownames(geneAndPeak.coloc.df),"coloc.H4.ABF"] = geneAndPeak.coloc.df[,"H4.ABF"]
scores.all[rownames(geneAndPeak.coloc.df),"coloc.H4vsH3.ratio"] = geneAndPeak.coloc.df[,"H4vsH3.ratio"]
scores.all[rownames(geneAndPeak.MR_P.df),"MR.qval"] = -p.adjust(geneAndPeak.MR_P.df[, "MR.pval"], method="BH")
scores.all[rownames(geneAndPeak.FMeQTL_DistInv.df),"FMeQTL2ARE.DistanceInv"] = geneAndPeak.FMeQTL_DistInv.df[, "FMeQTL.distInv"]
scores.all[rownames(gene2Pk.dis.df),"gene2ARE.DistanceInv"] = 1/gene2Pk.dis.df[, "distance"]


scores.all[,"MR.qval"][is.na(scores.all[,"MR.qval"])] = -1
#scores.all[,"PGS.pcor.qval"][is.na(scores.all[,"PGS.pcor.qval"])]= -1
scores.all[, "coloc.H4vsH3.ratio"][scores.all[, "coloc.H4.ABF"]<coloc.H4.ABF.CUTOFF.forRatio]=0
scores.all[is.na(scores.all)]=0

save(scores.all,
     file=file.ou.scores.RData)

##################################################################
#with different measurement as golden standard
#PRC curve

PRC.res.amongEach = lapply(names(GSP.cutoffs),
FUN=function(GSP.nm)
{
  classLab=rep(NA, nrow(scores.all))
  names(classLab)=rownames(scores.all)
  
  classLab[scores.all[, GSP.nm] >= GSP.cutoffs[GSP.nm]] = "fg"
  classLab[scores.all[, GSP.nm] < GSP.cutoffs[GSP.nm]] = "bg"
  
  res= lapply(setdiff(colnames(scores.all), GSP.nm),
  FUN=function(classifier.nm)
  {
    samps.ovlp = rownames(scores.all)[(!is.na(scores.all[, classifier.nm])) &
                                        (!is.na(scores.all[, GSP.nm])) & 
                                        scores.all[, "gene2ARE.DistanceInv"] <1/PROMT.WIND]
    
    prc = pr.curve(scores.class0 = scores.all[samps.ovlp,classifier.nm][classLab[samps.ovlp]=="fg"],
           scores.class1 = scores.all[samps.ovlp,classifier.nm][classLab[samps.ovlp]=="bg"], 
           curve=T)
    fg.prop = sum(classLab[samps.ovlp]=="fg")/length(samps.ovlp)
    
    return(list(prc= prc,
                fg.prop = fg.prop))
  })
  names(res) = setdiff(colnames(scores.all), GSP.nm)
  
  return(res)
})
names(PRC.res.amongEach) = names(GSP.cutoffs)

save(scores.all,
     PRC.res.amongEach,
     file=file.ou.scores.RData)

#visualization

pdf(file.ou.cmp.pdf, height=6, width = 24)
layout(matrix(1:4, ncol=4, byrow = T))
for(GS.nm in names(PRC.res.amongEach))
{
  plot(c(0,1), c(0,1), 
     xlab="recall", 
     ylab= "precision",
     type="n",
     main=paste(TISS, GS.nm, "as Gold Standard Positive"))
  
  # cols = rainbow(length(PRC.res.amongEach[[GS.nm]]))
  # names(cols) = names(PRC.res.amongEach[[GS.nm]])  
  nms.test=names(PRC.res.amongEach[[GS.nm]])
  for(nm  in nms.test)
  {
    lines(PRC.res.amongEach[[GS.nm]][[nm]]$prc$curve[,1:2], col=SCORE2COL[nm], lwd=2, lty=1)
    segments(x0=0,
             y0=PRC.res.amongEach[[GS.nm]][[nm]]$fg.prop,
             x1=1,
             y1=PRC.res.amongEach[[GS.nm]][[nm]]$fg.prop, col=SCORE2COL[nm], lwd=2, lty=2)
  }
  legend("topright", 
         legend=c(paste(names(PRC.res.amongEach[[GS.nm]]), "AUPRC:", signif(sapply(PRC.res.amongEach[[GS.nm]], FUN=function(x) x$prc$auc.integral),3)), 
                  paste(names(PRC.res.amongEach[[GS.nm]]), "baseline AUPRC:", signif(sapply(PRC.res.amongEach[[GS.nm]], FUN=function(x) x$fg.prop),3))),
         lwd=2,
         lty = rep(1:2, each= length(PRC.res.amongEach[[GS.nm]])),
         col= c(SCORE2COL[nms.test], SCORE2COL[nms.test])
         )
    
}
dev.off()

