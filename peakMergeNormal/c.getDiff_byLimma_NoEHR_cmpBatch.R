#V1.4.1.2
#c.getDiff_byLimma_V1.4.1.2_sva_NoEHR_cmpBatch.R

#use sva to consider unwanted variables
#compare samples of each batch to all control batches and all case batch separately
#leveraging the sv inferred before
args = commandArgs(trailingOnly=TRUE)
HM = args[1]
BATCH.i= as.numeric(args[2])

#DRUG.N.CUTOFF=5
COUNT.CUTOFF=5
PRES.PROP = 0.2
#TOP.NUM=1000


# HM2DIFF.Q= c(H3K27ac=0.01,
#              H3K4me1=0.01,
#              H3K36me3=0.1)
# POIS.DIFF.Q=HM2DIFF.Q[HM]

library(limma)
library(data.table)
#library("SmartSVA")
library(sva)
library(pheatmap)
# library(clusterProfiler)
# library(AnnotationHub)
# library(org.Hs.eg.db)
#

TYPE2COL=c("grey","red")
names(TYPE2COL)=c("control", "case") #c("IBDRA_Control", "CD", "UC", "RA")

file.in.meta.RData = paste0("a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.goodSamples.meta.RData")
#file.in.norm.RData = paste0("a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.heightsMean.depthGCNormed.RData")
#file.in.ref2pkOvlprefSlct.combBP= paste0("a_2_peakHeightsNormalized_V1.4/", HM, ".ref2pkOvlpWithTest.refSlct.combBP.RDS")
file.in.covarite.csv = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
file.in.sampInfo.RData = "../sampleMeta/a_sampleInfo_lithiumAndTobacco_V1.1_researchMeta/samp.info.RData"

file.in.DP.RData=paste0("c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HM, "/", HM, ".diffPeak.RData")

#file.in.deconv.RData =  paste0("../deconv/b_real_estFract_V1.4.3_nnPoisR_peakNorm_colinearity_6cellType/", HM, "/estFraction.results.nnPoisR.FiltQ", POIS.DIFF.Q, ".RData")
file.in.batch="~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/c_metaAndQC/QCed.Samples.vs.libBatch.txt"
#file.in.deconv.RData =  paste0("../deconv/b_real_estFract_V1.4.4.1_nnPoisR_noIntcpt_peakNorm_6cellType_autoSig/", HM, "/estFraction.results.nnPoisR.FiltQ", POIS.DIFF.Q, ".RData")


dir.ou=paste0("./c_bulkDiffPeak_limma_V1.4.1.2_sva_noEHR_compBatch/", HM, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.diffPeak.RData = paste0(dir.ou, HM, ".diffPeak.batch.", BATCH.i, ".RData")
#file.ou.diffPeak.pdf = paste(dir.ou, HM, ".diffPeak.pdf", sep="")
# file.ou.diffPeak.sva.diag.pdf = paste(dir.ou, HM, ".diffPeak.sva.diag.pdf", sep="")
# file.ou.diffPeak.bed.pref = paste(dir.ou, HM, ".diffPeak.", sep="")
# file.ou.diffPeak.BG.bed = paste(dir.ou, HM, ".diffPeak.BG.bed", sep="")



#load(file.in.norm.RData)
#ref2pkOvlpWithTest.refSlct = get('ref2pkOvlpWithTest.refSlct', e1)
#ref2pkOvlpWithTest.refSlct.combBP = readRDS(file.in.ref2pkOvlprefSlct.combBP)
#load(file.in.deconv.RData)
load(file.in.meta.RData)
load(file.in.sampInfo.RData)
load(file.in.DP.RData)
#deconv#######
# load(file.in.deconv.RData)

# Grp.slct.props.unNorm.mean = sapply(colnames(Grp.slct.props.sampling.unNorm[[1]]),
# FUN=function(ref)
# {
#   ref.props= sapply(Grp.slct.props.sampling.unNorm, FUN=function(x) x[,ref])
#   return(apply(ref.props, 1, mean))
# })



##############
covars = read.table(file.in.covarite.csv, sep=",", header=T, row.names=1)
rownames(covars)=paste("id", rownames(covars), sep="_")

batch.df=read.table(file.in.batch, sep="\t", header=T, row.names=1, stringsAsFactor=F)
batch2grp=lapply(split(batch.df$group, batch.df$batch.exp), FUN=function(x)x[!duplicated(x)])
batch.uniqGrp=names(batch2grp)[sapply(batch2grp, length)==1]
batch.uniqGrp.df=batch.df[batch.df$batch.exp %in% batch.uniqGrp, ]
grp2batch=lapply(split(batch.uniqGrp.df$batch.exp, batch.uniqGrp.df$group), FUN=function(x)x[!duplicated(x)])
batch2samp.id=split(rownames(batch.uniqGrp.df), batch.uniqGrp.df$batch.exp)
#
if(BATCH.i > length(batch.uniqGrp)) q("no")
print(BATCH.i)
batch.slct=batch.uniqGrp[BATCH.i]
batch.slct.grp=batch2grp[batch.slct]
batch.slct.sampIDs=batch2samp.id[[batch.slct]]
rest.list.sampIDs=list()
rest.list.sampIDs$control=unlist(batch2samp.id[setdiff(grp2batch$control, batch.slct)])
rest.list.sampIDs$case=unlist(batch2samp.id[setdiff(grp2batch$case, batch.slct)])

ids.all=unlist(dis2samps[c("control", "case")])


#
# pks.slct = rownames(ref2pkOvlpWithTest.refSlct.combBP)[apply(ref2pkOvlpWithTest.refSlct.combBP,1,sum)>=0.5]
# samps.counts.disease = round(samp2HeightsMean$depth.GCnormTeng[pks.slct,unlist(dis2samps[c("control", "case")])])
# pks.filt = pks.slct[apply(samps.counts.disease>=COUNT.CUTOFF, 1, sum)/ncol(samps.counts.disease)>=PRES.PROP]


case.id=intersect(batch.slct.sampIDs, ids.all)
if(length(case.id)==0) q("no")

res=lapply(rest.list.sampIDs,
FUN=function(rest.ids)
{
  ctrl.id=intersect(rest.ids, ids.all)

  meta.filt.df = data.frame(disease=c(rep("disease", length(case.id)), #including sv inferred for previous DP analysis
                                  rep("control", length(ctrl.id))),
                            design.matrix.mod2[c(case.id, ctrl.id), -(1:2)])
  counts.filt=samps.counts.disease[pks.filt, c(case.id, ctrl.id)]
  # # #smartSVA to estimate the surrogate variates
  # exp.r <- t(resid(lm(t(exp) ~., data=meta.df)))
  # # ## Add one extra dimension to compensate potential loss of 1 degree of freedom
  # # ## in confounded scenarios (very important)
  # n.sv <- EstDimRMT(exp.r, FALSE)$dim + 1
  # smartSVA.obj <- smartsva.cpp(as.matrix(exp), mod=design.matrix.mod, mod0=design.matrix.mod0, n.sv=n.sv)


  #voom + sva 
  lib.size = rep(1e6, nrow(meta.filt.df))
  
  design.matrix.mod2=model.matrix(~., meta.filt.df)
  v = voom(counts.filt, lib.size=lib.size, design =design.matrix.mod2, plot=F)
  fit <- lmFit(v, design.matrix.mod2)
  fit <- eBayes(fit)

  res=topTable(fit, coef="diseasedisease", adjust="BH", number=nrow(counts.filt))
})

save(
  batch.slct,
  batch.slct.grp,
  res, 
  file=file.ou.diffPeak.RData)




