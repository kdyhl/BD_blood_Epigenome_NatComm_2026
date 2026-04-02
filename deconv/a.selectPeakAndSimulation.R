#V1.7.2
#automatically select peaks rather than manually select signature peaks
#a.selectPeakAndSimulation_V1.7.2_colineariyAuto_rowNorm_6Celltypes.R

#V1.7.1
#row normalize 
#6 cell types, CD4, CD8, B, NK, MONO, Neutrophil


#V1.7 
#1) do not distinguish CD4 and CD8, Neutrophil and monocyte 
#2) besides selecting peaks based on DP peaks
#select based on colinearity  [https://www.nature.com/articles/s41467-019-09990-5]
#not remove from bulk differential data
#remove simulation to a separate script

#V1.6
#remove differential peaks from bulk analysis
#Pois different qvalue 0.05

#V1.5.5
#use dirichlet distrubiton for cell fraction


#V1.5.4
#put more constrict for the counts of peaks for deconvolution
#

#V1.5.1
#use neutrophil
#output merged peaks overlapped with reference peaks in cell type

#V1.5
#use non-negative poisson regression for fraction estimation
#

#V1.4
#use log10 (x+1), before deconvolution

#V1.3
#use un normalized fraction
#use residules and estimated fraction to diagnose

#V1.2
#plan to add DNase senstivity sites to filter, not enough data, maybe updated later
#output peak markers for annotation
#also test them on blueprint samples
#also test a way to better scale before regression
#test another way to rank and select markers
#only show qood qulity samples

#Version 1.1
#using features from peak length normalized samples
#use a different way to generate tissue specific peak list
# #use those peaks  whose max counts across reference samples overlapping with a peak in the same tissue
# #the max one should be larger enough than the second max one

#model changed to 
#a) mean of constitutents for each sample is the count from reference
#b) mean of constitutents for each sample is from  a poisson distribution with mean of the count from reference
#c) all constitutents mean modeling as b), except monocytes which coming from a random pick from blueprint sample of monocytes


#Version 1
#simulation could choose three models
#a) mean of constitutents for each sample is the count from reference
#b) mean of constitutents for each sample is from  a negative binomial distribution with mean of the count from reference
#c) all constitutents mean modeling as b), except monocytes which coming from a random pick from blueprint sample of monocytes

#test two algorithm
#robust linear regression and support vector regression

#also deconvolve the real data too

args = commandArgs(trailingOnly=TRUE)
HM = args[1] #"H3K27ac" #"H3K36me3"#"H3K4me1"#"H3K27ac"


R.LIBS.MYPATHs = c(V3.6 ="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/",
                   V3.5="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.5/",
                   V3.4="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.4/",
                   V3.3="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.3/")
R.LIBS.PATH = .libPaths()


if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}
.libPaths(c(R.LIBS.MYPATHs["V3.6"], R.LIBS.PATH))


#
# library(affy)
library(MASS)
#library(addreg)
library(MCMCpack)
# library(preprocessCore)
# library(EpiDISH)
#library(gcapc)
library("RColorBrewer")
# library("devtools")
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(pheatmap)
library(umap)
library(ggplot2)
library(Rtsne)
library(FNN)

#source("~/lhou.compbio/codes/mylib/RCode/H3K27ac.functions.R")


HM2Neutroph = c(H3K27ac = "BPNeutrop",
                H3K4me1 = "E030",
                H3K36me3 = "E030")


#REFs.SELECT= c("E040", "E039", "E048", "E047", "E029", "E032", "E046", "BPNeutrop") #"E037", "E038",#, "E037", "Neutrop")#, "E030")#"E044", "E041", "E042", "E043", "E048", "E034", "E039", "E043", "E045", 
#REF2COL=brewer.pal(n = 8, name = "Set2")
#REF2COL = rainbow(length(REFs.SELECT))
#names(REF2COL) = REFs.SELECT
GRP2REFs.SLCT = c(CD4T= "E040|E039|E037|E038",
                  CD8T= "E048|E047",
                  BCell = "E032",
                  NKCell = "E046",
                  Monocyte = "E029",
                  Neutrophil = HM2Neutroph[HM]) #BPNeutrop
names(GRP2REFs.SLCT)=c("CD4T", "CD8T", "BCell", "NKCell", "Monocyte", "Neutrophil")
GRP2COL = brewer.pal(n = length(GRP2REFs.SLCT), name = "Set2")
names(GRP2COL) = names(GRP2REFs.SLCT)
#
BULK.DP.PVAL=0.01

POIS.DIFF.Q=c(H3K27ac=0.01,
              H3K4me1=0.01,
              H3K36me3=0.1)[HM]
#DIFF.PEAK.COUNT.UB = 40 #V1.5.4
DIFF.PEAK.COUNT.LB = 10
CLUSTER.CENTER.PCC =0.2
#DIFF.PEAK.COUNT.LB = 5

#MEAN.NB.DISPERS = 0.1
SAMPLE.NB.DISPERS = 1 #c(0.1, 0.5, 1, 2, 5)# c(20, 10, 5, 2, 1)#, 0.5, 0.2, 0.1, 0.05)
POISSON.MEAN = 5
INDIVs =rep("health", 40)
SIM.MODEL = c("mean=Ref", "mean=sampleFromRef","mean=sampleFromRef+BPmono")[2]
EXPM.NO = 20
SAMP.FRACT.DIRICHLET.ALPHA = c(30, 15, 10, 15, 10, 40)



#

file.in.norm.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.heightsMean.depthGCNormed.RData")
file.in.meta.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.goodSamples.meta.RData")
#file.in.testMeta="~/lhou.compbio/data/pfizer_aim2_5.2017/H3K27ac_3.27.2018/sample.metaData.until3.27.2018.txt"
file.in.ref2pkOvlpTest.RDS = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".ref2pkOvlpWithTest.refSlct.combBP.RDS")
#file.in.testMeta="~/lhou.compbio/data/pfizer_aim2_5.2017/H3K27ac_3.27.2018/sample.metaData.until3.27.2018.txt"
#file.in.DNase = "/broad/compbio/anshul/projects/roadmap/peaks/consolidated/narrowPeak/E000-DNase.macs2.narrowPeak.gz"
files.in.bulkDP.RData=paste0("../peakMergeNormal/c_bulkDiffPeak_deseq2_V1.2_addCov/", HM, ".diffPeak.RData")
names(files.in.bulkDP.RData) =c("BP")
#dir.tmp = "~/hptmp/pfizer2/"

dir.ou = paste0("a_selectPeaksAndSimulation_V1.7.2_colinearAuto_rowNorm_6cellType/", HM, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
script.nm = sub(".*=", "", commandArgs()[4])
system(paste0("cp ", script.nm, " ", dir.ou))

#PEKAS.SLCT.NUM = 2400
file.ou.cellSpePeak.pref = paste0(dir.ou,  "cellTypeSpecPeaks.")
file.ou.cellDiffPeak.pref = paste0(dir.ou,  "cellTypeDiffPeaks.")
#file.ou.real.pdf = paste(dir.ou, "estFractionForPfizer.results.RPC.FiltQ", POIS.DIFF.Q, ".pdf", sep="")
file.ou.RData = paste0(dir.ou,  "signaturePeaksAndSimulation.FiltQ", POIS.DIFF.Q, ".above", DIFF.PEAK.COUNT.LB,  ".RData")

#
load(file.in.norm.RData) #save(samp2HeightsMean, ref2pkOvlpWithTest, pk2refOvlpTest.Weightsum, file=file.ou.RData)
#merged.peaks = rownames(samp2HeightsMean[[1]])
#merged.peaks.inGrp = merged.peaks[pk2refOvlpTest.Weightsum>=0.5]
load(file.in.meta.RData)
samps= unlist(dis2samps)
samp.dis = unlist(lapply(names(dis2samps), FUN=function(nm) {rep(nm, length(dis2samps[[nm]]))}))
names(samp.dis) = samps

colnames(samp2HeightsMean$depth.GCnormTeng) = gsub("BP-", "BP", colnames(samp2HeightsMean$depth.GCnormTeng) )
BP.monoSamp = colnames(samp2HeightsMean$depth.GCnormTeng)
BP.monoSamp = BP.monoSamp[grepl("BPmono", BP.monoSamp)]
#colnames(ref2pkOvlpWithTest) = gsub("BP-", "BP",colnames(ref2pkOvlpWithTest) )

ref2pkOvlpWithTest.refSlct.combBP = readRDS(file.in.ref2pkOvlpTest.RDS)

TEST.SAMP.SIZE= nrow(meta.goodsamples)
samp2HeightsMean.depth.GCnormTeng = samp2HeightsMean$depth.GCnormTeng[,1:TEST.SAMP.SIZE]

samp2HeightsMean.depth.GCnormTeng.medians = apply(samp2HeightsMean.depth.GCnormTeng, 1, median)
peaks.all = rownames(samp2HeightsMean.depth.GCnormTeng)[samp2HeightsMean.depth.GCnormTeng.medians>=1] #filter those peaks too small in bulk level
samp2HeightsMean.depth.GCnormTeng = samp2HeightsMean.depth.GCnormTeng[peaks.all, ]
###############################################
grp.counts = sapply(GRP2REFs.SLCT,
FUN=function(ref)
{
  if(ref %in% colnames(samp2HeightsMean$depth.GCnormTeng))
  {
    return(samp2HeightsMean$depth.GCnormTeng[peaks.all,ref])
  }else
  {
    samps = colnames(samp2HeightsMean$depth.GCnormTeng)
    samps = samps[grepl(ref, samps)]
    apply(samp2HeightsMean$depth.GCnormTeng[peaks.all,samps], 1, median)
  }
})
colnames(grp.counts)= names(GRP2REFs.SLCT)

samp2HeightsMean.depth.GCnormTeng.medians = apply(samp2HeightsMean.depth.GCnormTeng, 1, median)
samp2HeightsMean.depth.GCnormTeng.medians.median =  median(samp2HeightsMean.depth.GCnormTeng.medians)
samp2HeightsMean.depth.GCnormTeng.medians.factor = samp2HeightsMean.depth.GCnormTeng.medians.median/samp2HeightsMean.depth.GCnormTeng.medians
samp2HeightsMean.depth.GCnormTeng.rowNorm = sweep(samp2HeightsMean.depth.GCnormTeng, 1, 
                                                       samp2HeightsMean.depth.GCnormTeng.medians.factor, "*" )
grp.counts.rowNorm = sweep(grp.counts, 1, samp2HeightsMean.depth.GCnormTeng.medians.factor, "*" )


# ref.counts = sapply(REFs.SELECT,
# FUN=function(ref)
# {
#   if(ref %in% colnames(samp2HeightsMean$depth.GCnormTeng))
#   {
#     return(samp2HeightsMean$depth.GCnormTeng[,ref])
#   }else
#   {
#     samps = colnames(samp2HeightsMean$depth.GCnormTeng)
#     samps = samps[grepl(ref, samps)]
#     apply(samp2HeightsMean$depth.GCnormTeng[,samps], 1, median)
#   }
# }) 

# 
# bulkDPs=lapply(files.in.bulkDP.RData,
# FUN=function(f.i.RData)
# {
#   load(f.i.RData)
#   dps=rownames(res)[res@listData$pvalue<=BULK.DP.PVAL]
#   dps[!is.na(dps)]
# })
# bulkDPs = levels(factor(unlist(bulkDPs)))
# 


#
print("feature selection")##########################################################

#feature selections
pk2refOvlpTestNum.GrpRef =  sapply(GRP2REFs.SLCT,
FUN=function(ref)
{
  if(ref %in% colnames(ref2pkOvlpWithTest.refSlct.combBP))
  {
    return(ref2pkOvlpWithTest.refSlct.combBP[peaks.all,ref])
  }else
  {
    samps = colnames(ref2pkOvlpWithTest.refSlct.combBP)
    samps = samps[grepl(ref, samps)]
    apply(cbind(ref2pkOvlpWithTest.refSlct.combBP[peaks.all,samps]), 1, mean)
  }
})
names(pk2refOvlpTestNum.GrpRef) = names(GRP2REFs.SLCT)
saveRDS(pk2refOvlpTestNum.GrpRef, file=paste0(dir.ou, "pk2refOvlpTestNum.GrpRef.RDS"))  
#peaks in each cell type
for(GRP in names(GRP2REFs.SLCT))
{
  f.o.bed = paste(dir.ou, "AllPeaks.", GRP, ".bed", sep="")
  write(gsub(":|-", "\t", rownames(pk2refOvlpTestNum.GrpRef)[pk2refOvlpTestNum.GrpRef[,GRP]>=0.5]), file=f.o.bed)
}



#cell type specific peak
peaks.cellSpecif = rownames(pk2refOvlpTestNum.GrpRef)[apply(pk2refOvlpTestNum.GrpRef, 1, FUN=function(x){sum(x>=0.5)==1})]
peaks.cellSpecif.top2max = t(apply(grp.counts[peaks.cellSpecif, ],1, FUN=function(x){sort(x,decreasing = T)[1:2]}))
peaks.cellSpecif.maxVs2ndMax = peaks.cellSpecif.top2max[,1]/(peaks.cellSpecif.top2max[,2]+1)
peaks.cellSpecif.MaxWithPeak.DiffQ = p.adjust(apply(round(peaks.cellSpecif.top2max[peaks.cellSpecif,]),1,
FUN=function(x)
{
  ppois(x[2],x[1])
}), method="BH")

peaks.cellSpecif.isT =  pk2refOvlpTestNum.GrpRef[peaks.cellSpecif, ]>=0.5 
peaks.cellSpecif.counts = t(grp.counts[peaks.cellSpecif, ])[t(peaks.cellSpecif.isT)]
names(peaks.cellSpecif.counts) = peaks.cellSpecif
print(paste0(sum(peaks.cellSpecif.counts==peaks.cellSpecif.top2max[,1]), " out of ", length(peaks.cellSpecif)))

peaks.cellSpecif.filt = peaks.cellSpecif[peaks.cellSpecif.counts==peaks.cellSpecif.top2max[,1] 
                                         & peaks.cellSpecif.MaxWithPeak.DiffQ<= POIS.DIFF.Q
                                         & peaks.cellSpecif.top2max[,1] >= DIFF.PEAK.COUNT.LB
                                        #  #& peaks.cellSpecif.top2max[,1]<=100
                                        #& peaks.cellSpecif.top2max[,1] <=DIFF.PEAK.COUNT.UB
                                        #    & peaks.cellSpecif.maxVs2ndMax >=2
                                          ]




peaks.cellSpecif.filt.list = lapply(names(GRP2REFs.SLCT),
FUN=function(grp)
{
  print(grp)
  pks=peaks.cellSpecif.filt[pk2refOvlpTestNum.GrpRef[peaks.cellSpecif.filt,grp]==1]
  
  if(length(pks)==0)
  {
    return(NA)
  }
  plot(peaks.cellSpecif.top2max[pks,], main=grp, xlog=T, ylog=T)
  
  
  f.o.pk = paste(file.ou.cellSpePeak.pref, "cellSpecPeak.", grp, ".bed", sep="")
  pks.bed = gsub(":|-", "\t", pks)
  write(pks.bed, f.o.pk)
  return(pks)
  #sample(peaks.cellSpecif.filt[ref2pkOvlpWithTest.refSlct.combBP[peaks.cellSpecif.filt,ref]==1],300)
})
names(peaks.cellSpecif.filt.list) = names(GRP2REFs.SLCT)

peaks.cellSpecif.filt.all = unlist(peaks.cellSpecif.filt.list)
peaks.cellSpecif.filt.grp = unlist(sapply(names(GRP2REFs.SLCT), FUN=function(grp){rep(grp, length(peaks.cellSpecif.filt.list[[grp]]))}))
names(peaks.cellSpecif.filt.grp) =peaks.cellSpecif.filt.all



#######################
#differential peaks
peaks.inGrp = rownames(pk2refOvlpTestNum.GrpRef)[apply(pk2refOvlpTestNum.GrpRef, 1, FUN=function(x) any(x>=0.5))]
peaks.inGrp.top2max = t(apply(grp.counts[peaks.inGrp, ],1, FUN=function(x){sort(x,decreasing = T)[1:2]}))
#peaks.inGrp.maxVs2ndMax = peaks.inGrp.top2max[,1]/(peaks.inGrp.top2max[,2]+1)
peaks.inGrp.isT = pk2refOvlpTestNum.GrpRef[peaks.inGrp, ]>=0.5   
peaks.inGrp.isMax = grp.counts[peaks.inGrp, ] == peaks.inGrp.top2max[,1]
peaks.inGrp.MaxWithPeak = peaks.inGrp[apply(peaks.inGrp.isT & peaks.inGrp.isMax, 1, any)]
peaks.inGrp.MaxWithPeak.DiffQ = p.adjust(apply(round(peaks.inGrp.top2max[peaks.inGrp.MaxWithPeak,]),1,
FUN=function(x)
{
  ppois(x[2],x[1])
}), method="BH")
peaks.inGrp.diff = peaks.inGrp.MaxWithPeak[peaks.inGrp.MaxWithPeak.DiffQ<=POIS.DIFF.Q]

peaks.inGrp.diff.list = lapply(names(GRP2REFs.SLCT),
FUN=function(grp)
{
  #print(ref)
  pks=peaks.inGrp.diff[pk2refOvlpTestNum.GrpRef[peaks.inGrp.diff,grp]>=0.5 & 
                         grp.counts[peaks.inGrp.diff, grp] == peaks.inGrp.top2max[peaks.inGrp.diff,1]]
  plot(peaks.inGrp.top2max[pks,], main=grp, xlog=T, ylog=T)
  
  return(pks)
  #sample(peaks.inGrp.filt[pk2refOvlpTestNum.GrpRef[peaks.inGrp.filt,ref]==1],300)
})
names(peaks.inGrp.diff.list) = names(GRP2REFs.SLCT)

dev.off()
# 
# 
# lapply(REFs.SELECT,
# FUN=function(ref)
# {
#   #pri
#   f.o.pk = paste(file.ou.cellSpePeak.pref, "diffPeak.FiltQ", POIS.DIFF.Q, ".", ref, ".bed", sep="")
#   pks.bed = gsub(":|-", "\t", peaks.inGrp.diff.list[[ref]])
#   write(pks.bed, f.o.pk)
#   
# })

peaks.inGrp.diff.filt.list = lapply(names(GRP2REFs.SLCT),
FUN=function(grp)
{
  #print(ref)
  pks=peaks.inGrp.diff.list[[grp]]
    
  pks.filt=pks[ #grp.counts[pks, grp] <= DIFF.PEAK.COUNT.UB&
                  grp.counts[pks, grp] >= DIFF.PEAK.COUNT.LB
                ]
  #pks.filt =  setdiff(pks.filt, bulkDPs)
  #
  # f.o.pk = paste(file.ou.cellSpePeak.pref, "diffPeak.FiltQ", POIS.DIFF.Q, ".rmBulkDP.p", BULK.DP.PVAL, ".", DIFF.PEAK.COUNT.LB, "-", DIFF.PEAK.COUNT.UB, ".", grp, ".bed", sep="") #, 
  # pks.bed = gsub(":|-", "\t", pks.filt)
  # write(pks.bed, f.o.pk)
  # 
  return(pks.filt)
  #sample(peaks.inGrp.filt[pk2refOvlpTestNum.GrpRef[peaks.inGrp.filt,ref]==1],300)
})
names(peaks.inGrp.diff.filt.list) = names(GRP2REFs.SLCT)

# cellDiff.pks.minNum = min(sapply(peaks.inGrp.diff.filt.list,length))
# cellDiff.pks.minNum.grp = names(GRP2REFs.SLCT)[sapply(peaks.inGrp.diff.filt.list,length)==cellDiff.pks.minNum]
# cellDiff.pks.minNum.grp.maxCounts = max(grp.counts[peaks.inGrp.diff.filt.list[[cellDiff.pks.minNum.grp]],cellDiff.pks.minNum.grp])

#####
#explore  colinearity
#

#samp2HeightsMean.depth.GCnormTeng = samp2HeightsMean$depth.GCnormTeng
peaks.inGrp.diff.filt.all = unlist(peaks.inGrp.diff.filt.list)
peaks.inGrp.diff.filt.grp = unlist(sapply(names(GRP2REFs.SLCT), FUN=function(grp){rep(grp, length(peaks.inGrp.diff.filt.list[[grp]]))}))
names(peaks.inGrp.diff.filt.grp )=peaks.inGrp.diff.filt.all
#peaks.inGrp.diff.filt.cors= cor(t(samp2HeightsMean.depth.GCnormTeng[peaks.inGrp.diff.filt.all,]))

peaks.inGrp.diff.filt.Rtsne = Rtsne(samp2HeightsMean.depth.GCnormTeng.rowNorm[peaks.inGrp.diff.filt.all,])
colnames(peaks.inGrp.diff.filt.Rtsne$Y) =c("x", "y")
peaks.inGrp.diff.filt.Rtsne.df = data.frame(peaks.inGrp.diff.filt.Rtsne$Y,
                                  cellType=peaks.inGrp.diff.filt.grp)#,
                                 # isBulkDP=peaks.inGrp.diff.filt.all %in% bulkDPs)
rownames(peaks.inGrp.diff.filt.Rtsne.df) = peaks.inGrp.diff.filt.all

# #remove  subtype of peak from NK cell type sitting with CD peaks from tsne
# #
# print("need to manually find the right formula based on visualization of peaks.inGrp.diff.filt.Rtsne.df below") #c(-10, 25) c(10, 0)
# peaks.inGrp.diff.rmSubsets.list=peaks.inGrp.diff.filt.list
# # pks.CD4T = peaks.inGrp.diff.filt.list$CD4T
# # peaks.inGrp.diff.rmSubsets.list$CD4T = pks.CD4T[peaks.inGrp.diff.filt.Rtsne.df[pks.CD4T, "y"] <= -8] 
# #
# pks.NKCell = peaks.inGrp.diff.filt.list$NKCell
# peaks.inGrp.diff.rmSubsets.list$NKCell = pks.NKCell[#peaks.inGrp.diff.filt.Rtsne.df[pks.NKCell, "y"] >= 5 &
#                                                       peaks.inGrp.diff.filt.Rtsne.df[pks.NKCell, "x"] >=-17] 
# #
# pks.BCell = peaks.inGrp.diff.filt.list$BCell
# peaks.inGrp.diff.rmSubsets.list$BCell = pks.BCell[#peaks.inGrp.diff.filt.Rtsne.df[pks.BCell,"y"]>= 10 &
#                                                     peaks.inGrp.diff.filt.Rtsne.df[pks.BCell,"x"] <= -15 ]
# #
# pks.Monocyte = peaks.inGrp.diff.filt.list$Monocyte
# peaks.inGrp.diff.rmSubsets.list$Monocyte = pks.Monocyte[peaks.inGrp.diff.filt.Rtsne.df[pks.Monocyte,"x"]<=0 &
#                                                           #peaks.inGrp.diff.filt.Rtsne.df[pks.Monocyte,"x"]>=0 &
#                                                           peaks.inGrp.diff.filt.Rtsne.df[pks.Monocyte,"y"]<=-10 #&
#                                                           #peaks.inGrp.diff.filt.Rtsne.df[pks.Monocyte,"y"]>= -10
#                                                           ]
# 
# pks.Neutrophil = peaks.inGrp.diff.filt.list$Neutrophil
# peaks.inGrp.diff.rmSubsets.list$Neutrophil = pks.Neutrophil[peaks.inGrp.diff.filt.Rtsne.df[pks.Neutrophil,"y"]>=5 |
#                                                               peaks.inGrp.diff.filt.Rtsne.df[pks.Neutrophil,"x"]>=5 
#                                                               ]

extractSignaturePeaks=function(peaks.inGrp.diff.filt.Rtsne.df, 
                               K=100, 
                               prop.cutoff=0.8, 
                               cellGrp.sigPeak.minSize=500, # proportion cutoff in the nearest neighbors
                               cellGrp.filt.minSize=100 #smaller than this size no filter will be done in this step
                               )
{
  #dis.matrix = parDist(as.matrix(peaks.inGrp.diff.filt.Rtsne.df[,1:2]), method="euclidean") #parDist ##keep the in dist format after dist function in function getConstraintDis
  knn.res=get.knn(as.matrix(peaks.inGrp.diff.filt.Rtsne.df[,1:2]), k=K)
  
  #peak to neighbor proportion
  cellGrp2pks = split(rownames(peaks.inGrp.diff.filt.Rtsne.df), as.character(peaks.inGrp.diff.filt.Rtsne.df[,3]))
  cellGrp2size = sapply(cellGrp2pks, length)
  cellGrp2sigPkMinSize = sapply(cellGrp2size, FUN=function(x) return(min(x, cellGrp.sigPeak.minSize)))
  cellGrps = names(cellGrp2pks)
  pk2neighbCellGrpCounts = t(apply(knn.res$nn.index, 1,
  FUN=function(ids)
  {
    summary(factor(as.character(peaks.inGrp.diff.filt.Rtsne.df[ids, 3]), levels=cellGrps))
  }))
  rownames(pk2neighbCellGrpCounts) = rownames(peaks.inGrp.diff.filt.Rtsne.df)
  
  pk2cellGrpPropMax = apply(pk2neighbCellGrpCounts/K, 1, max)
  pk2cellGrp.Max = colnames(pk2neighbCellGrpCounts)[apply(pk2neighbCellGrpCounts/K, 1, which.max)]
  
  pks.isSig = rownames(pk2neighbCellGrpCounts)[pk2cellGrpPropMax>=prop.cutoff & pk2cellGrp.Max==peaks.inGrp.diff.filt.Rtsne.df[,3]] 
  cellGrp2sigPks =  split(pks.isSig, 
                          as.character(peaks.inGrp.diff.filt.Rtsne.df[pks.isSig, 3]))
  
  cellGrps.enoughSigPks = names(cellGrp2sigPks)[sapply(cellGrp2sigPks, length) >= cellGrp2sigPkMinSize[names(cellGrp2sigPks)]]
  if(length(cellGrps.enoughSigPks) == length(cellGrps))
  {
    return(cellGrp2sigPks)
  }
  cellGrp2sigPks = cellGrp2sigPks[cellGrps.enoughSigPks]
  
  #second round fiter
  cellGrps.left = setdiff(cellGrps,cellGrps.enoughSigPks)
  pks.left= unlist(cellGrp2pks[cellGrps.left])
  pk2neighbCellGrpCounts.left = pk2neighbCellGrpCounts[pks.left, cellGrps.left]
  pks.left= pks.left[apply(pk2neighbCellGrpCounts.left, 1, sum)!=0]
  pk2cellGrpPropMax.left = apply(pk2neighbCellGrpCounts.left[pks.left,], 1, max)/apply(pk2neighbCellGrpCounts.left[pks.left, ], 1, sum)
  pk2cellGrp.Max.left = colnames(pk2neighbCellGrpCounts.left)[apply(pk2neighbCellGrpCounts.left[pks.left, ], 1, which.max)]
  
  pksLeft.isSig = pks.left[pk2cellGrpPropMax.left>=prop.cutoff & pk2cellGrp.Max.left==peaks.inGrp.diff.filt.Rtsne.df[pks.left,3]]
  cellGrp2sigPks.left =  split(pksLeft.isSig, 
                          as.character(peaks.inGrp.diff.filt.Rtsne.df[pksLeft.isSig, 3]))
  
  
  cellGrp2sigPks =c(cellGrp2sigPks, cellGrp2sigPks.left)
  
  for(grp in names(cellGrp2pks))
  {
    if(cellGrp2size[grp] <=cellGrp.filt.minSize)
    {
      cellGrp2sigPks[[grp]]=cellGrp2pks[[grp]]
    }
  }
  
  
  return(cellGrp2sigPks)
}

#peaks.inGrp.diff.subset.list = extractSignaturePeaks(peaks.inGrp.diff.filt.Rtsne.df, K=100, prop.cutoff=0.9, cellGrp.sigPeak.minSize=500)
peaks.inGrp.diff.subset.list = extractSignaturePeaks(peaks.inGrp.diff.filt.Rtsne.df, K=50, prop.cutoff=0.8, cellGrp.sigPeak.minSize=200)

#peak center
print("peak center")


peaks.inGrp.diff.rmSubsets.refCenter= sapply(peaks.inGrp.diff.subset.list,
FUN=function(pks)
{
  apply(samp2HeightsMean.depth.GCnormTeng.rowNorm[pks,], 2, median)
})
peaks.inGrp.diff.deconv.list =lapply(names(peaks.inGrp.diff.subset.list),
FUN=function(ref)
{
  cors2center=cor(peaks.inGrp.diff.rmSubsets.refCenter[,ref], t(samp2HeightsMean.depth.GCnormTeng.rowNorm[peaks.inGrp.diff.subset.list[[ref]],]))
  colnames(cors2center)[cors2center>=CLUSTER.CENTER.PCC]
})
names(peaks.inGrp.diff.deconv.list) = names(peaks.inGrp.diff.subset.list)
#
peaks.inGrp.diff.deconv.list.all = unlist(peaks.inGrp.diff.deconv.list)
peaks.inGrp.diff.deconv.list.Rtsne = Rtsne(samp2HeightsMean.depth.GCnormTeng.rowNorm[peaks.inGrp.diff.deconv.list.all,])
colnames(peaks.inGrp.diff.deconv.list.Rtsne$Y) =c("x", "y")
peaks.inGrp.diff.deconv.list.Rtsne.df = data.frame(peaks.inGrp.diff.deconv.list.Rtsne$Y,
                                  cellType=peaks.inGrp.diff.filt.grp[peaks.inGrp.diff.deconv.list.all])
#




#
pdf(paste0(dir.ou, "Peaks.deconv.tsne.DiffQ", POIS.DIFF.Q, ".diag.pdf"), height=15, width=15)
# p=ggplot(peaks.inGrp.diff.filt.Rtsne.df, aes(x=x, y=y))  +  
#     geom_point(aes(col=isBulkDP), size=1) +
#     #scale_colour_manual(values=GRP2COL) +
#     ggtitle("peaks.inGrp.diff.filt")
#   
# print(p)

p=ggplot(peaks.inGrp.diff.filt.Rtsne.df, aes(x=x, y=y))  +  
    geom_point(aes(col=cellType), size=1) +
    scale_colour_manual(values=GRP2COL) +
    ggtitle("peaks.inGrp.diff.filt")
  
print(p)

p=ggplot(peaks.inGrp.diff.deconv.list.Rtsne.df, aes(x=x, y=y))  +  
  geom_point(aes(col=cellType), size=2) +
  scale_colour_manual(values=GRP2COL) +
  ggtitle("peaks.inGrp.diff.deconv")

print(p)
dev.off()



for(grp in names(GRP2REFs.SLCT))
{
  #print(ref)
  pks=peaks.inGrp.diff.deconv.list[[grp]]
    
  #
  f.o.pk = paste(dir.ou, "diffPeak.deconv.FiltQ", POIS.DIFF.Q, ".above", DIFF.PEAK.COUNT.LB, ".cluCenterPCC", CLUSTER.CENTER.PCC, ".", grp, ".bed", sep="")
  pks.bed = gsub(":|-", "\t", pks)
  write(pks.bed, f.o.pk)
  
}

#
pdf(paste0(dir.ou, "Peaks.deconv.signatures.pdf"), height=9, width=9)
pheatmap(log10(grp.counts[peaks.inGrp.diff.deconv.list.all,]+1),
         cluster_row=F, cluster_col=F,
         show_rownames=F,
         show_colnames=T,
         main="log10 signal of signature peaks in reference samples"
         )

peaks.inGrp.diff.deconv.list.testMedian=lapply(peaks.inGrp.diff.deconv.list,
FUN=function(pks)
{
  apply(samp2HeightsMean.depth.GCnormTeng[pks, ], 1, median)
})
boxplot(peaks.inGrp.diff.deconv.list.testMedian,
        col=GRP2COL,
        ylab="median signals of signature peaks across test samples")

dev.off()





print("simulation")########################################################

peaks.sim.all = unlist(peaks.inGrp.diff.deconv.list)
#pdf(file.ou.log.pdf, height=9, weight=12)
sim.counts=lapply(1:EXPM.NO,
FUN=function(i)
{
  print(i)
  # ind2fractions = sapply(1:length(INDIVs),
  # FUN=function(i)
  # {
  #   fracts.i = runif(ncol(ref.counts))
  #   fracts.i = fracts.i/sum(fracts.i)
  # })
  ind2fractions = t(rdirichlet(length(INDIVs), SAMP.FRACT.DIRICHLET.ALPHA))
  rownames(ind2fractions) = colnames(grp.counts)
  
  sim.counts.wholeBlood = sapply(1:length(INDIVs),
  FUN=function(i)
  {
    grp.counts.i = sapply(colnames(grp.counts),
    FUN=function(ref)
    {
      negbin.mean = grp.counts[peaks.sim.all,ref]
      
      if(SIM.MODEL=="mean=Ref")
      {
        return(rnegbin(length(peaks.sim.all), mu=negbin.mean, theta=1/SAMPLE.NB.DISPERS) + rpois(length(peaks.sim.all), POISSON.MEAN))
      }else if(SIM.MODEL=="mean=sampleFromRef")
      {
        return(rnegbin(length(peaks.sim.all), mu=rpois(length(peaks.sim.all), negbin.mean), theta=1/SAMPLE.NB.DISPERS) + rpois(length(peaks.sim.all), POISSON.MEAN))
      }else if(SIM.MODEL=="mean=sampleFromRef+BPmono")
      {
        if(ref=="E029")
        {
          return(samp2HeightsMean$depth.GCnormTeng[peaks.sim.all, sample(BP.monoSamp,1)])   
        }else
        {
          return(rnegbin(length(peaks.sim.all), mu=rpois(length(peaks.sim.all), negbin.mean), theta=1/SAMPLE.NB.DISPERS) + rpois(length(peaks.sim.all), POISSON.MEAN))
        }
      }
      
    })
    grp.counts.i %*% cbind(ind2fractions[,i])
  })
  rownames(sim.counts.wholeBlood) = peaks.sim.all
  colnames(sim.counts.wholeBlood) = 1:length(INDIVs)
  
  return(list(ind2fractions=ind2fractions, sim.counts.wholeBlood=sim.counts.wholeBlood))
})





save(
     ref2pkOvlpWithTest.refSlct.combBP,
     peaks.cellSpecif.filt.list, 
     peaks.inGrp.diff.list, 
     peaks.inGrp.diff.filt.list, 
     peaks.inGrp.diff.subset.list, 
     peaks.inGrp.diff.deconv.list,# finally the right set to use
     peaks.inGrp.diff.filt.Rtsne.df,
     peaks.inGrp.diff.deconv.list.Rtsne.df,
     samp2HeightsMean.depth.GCnormTeng, 
     samp2HeightsMean.depth.GCnormTeng.rowNorm,
     grp.counts,
     grp.counts.rowNorm,
     sim.counts,
     file=file.ou.RData)



