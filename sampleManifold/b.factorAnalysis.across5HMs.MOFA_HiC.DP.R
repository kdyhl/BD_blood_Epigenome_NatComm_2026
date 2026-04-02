#V1.2 
#b.factorAnalysis.across5HMs.MOFA_HiC.DP_V1.2_pvalCutoff.R
#instead of top 10000
#select number of peaks based on pvalue


#instead of cluster with each histone marks, we cluster samples with all of them together 
#instead of using PCs, using peaks
#under a NMF framework, called MOFA
#
# R.LIBS.MYPATHs = c(V3.5="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.5/",
#                 V3.4="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.4/",
#                 V3.3="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.3/")
#args = commandArgs(trailingOnly=TRUE)
#TYPE= args[1]  #DPs HiC.DPs topCVPks

R.LIBS.PATH = .libPaths()
.libPaths(c("~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/", R.LIBS.PATH))


library(reticulate)
library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)
library(pheatmap)
use_python("/broad/compbio/lhou/software/anaconda2/bin/python", required = TRUE)


TYPE = "HiC.DPs"#"topCVPks"

DP.PVAL.CUTOFF=0.01
#PEAKS.TOP.N = 10000 #c(50, 40, 30, 20, 10 ,5)
HMs= c("H3K27ac", "H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")

HiCBlock.DPEnrich.Q.CUTOFF=0.1

MOFA.FACTORS.N=40
MOFA.PROP.VAR.CUTOFF = 0.005
MOFA.SEED =666


#
file.in.norm.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".MayoWithRef.heightsMean.depthGCNormed.RData")
names(file.in.norm.RData) = HMs

# files.in.peakHeights.RData =paste0("../deconv/a_selectPeaksAndSimulation_V1.7.2_colinearAuto_rowNorm_6cellType/", HMs, "/signaturePeaksAndSimulation.FiltQ", HM2DIFF.Q[HMs], ".above10.RData")
# names(files.in.peakHeights.RData) = HMs

#files.in.DP.RData = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.0_addLiMed/", HMs, ".diffPeak.RData")
files.in.DP.RData = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
names(files.in.DP.RData) = HMs

file.in.HiC.block.RData="../hiC.DiffPeakandDEGandTF/a_2_cmpBulkDPAcrossHMs_hiCBlock/DP.V1.4.1_sva_noEHR/Bulk.hicBlocks.DPstat.RData"
file.in.HiCRegion2Pk.RData="../hiC.DiffPeakandDEGandTF/a_hiCLinks_withPeakAnnot_bulk/hiCLinks.RData"


file.in.meta1 = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
file.in.meta2 = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"

dir.ou=paste0("b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/", TYPE, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData = paste0(dir.ou, "ind.MOFA.DP-pval", DP.PVAL.CUTOFF, ".varProp", MOFA.PROP.VAR.CUTOFF, ".RData")
file.ou.pdf = paste0(dir.ou, "ind.MOFA.DP-pval", DP.PVAL.CUTOFF, ".varProp", MOFA.PROP.VAR.CUTOFF, ".pdf")
#file.ou.cluCmp.pdf = paste0(dir.ou, "patientHeter.cmpAcrossHMs.pdf")
#
#

meta1 = read.table(file.in.meta1, sep=",", header=T, row.names=1)
rownames(meta1)=paste("id", rownames(meta1), sep="_")
meta2 = read.table(file.in.meta2, sep=",", header=T, row.names=1)
rownames(meta2)=paste("id", rownames(meta2), sep="_")

meta.all = cbind(meta2, group=meta1[rownames(meta2), "group"])
#samps.case = rownames(meta)[meta$group=="case"]


#

#check meta and samples

# for(hm in HMs)
# {
#   print(hm)
#   load(paste0("a_2_peakHeightsNormalized_V1.4/", hm, ".MayoWithRef.goodSamples.meta.RData"))
#   load(paste0("a_2_peakHeightsNormalized_V1.4/", hm, ".MayoWithRef.heightsMean.depthGCNormed.RData"))

#   print(nrow(meta.goodsamples))
#   samps.ovlp = intersect(rownames(meta.all), colnames(samp2HeightsMean$depth.GCnormTeng))
#   print(length(samps.ovlp))
# }


#




load(file.in.HiC.block.RData)
#load(file.in.HiCRegion2DP.RData)
load(file.in.HiCRegion2Pk.RData)

hiCBlock.slct = rownames(hiCBlock.DP.summary.adjp)[apply(hiCBlock.DP.summary.adjp<=HiCBlock.DPEnrich.Q.CUTOFF,1, any)] #select differential peaks in side DP.enriched hicBlocks
hiCRegion.slct = unlist(hicBlock2hicRegs[hiCBlock.slct])

HM.matrix.list = lapply(HMs,
FUN=function(hm)
{
  print(hm)
  load(file.in.norm.RData[hm])
  load(files.in.DP.RData[hm])

  
  hicReg2pks=lapply(hicReg.df[[hm]],
  FUN=function(x)
  {
    unlist(strsplit(x, split=","))
  })
  names(hicReg2pks) = rownames(hicReg.df)
  
  hic.pks.slct = levels(factor(unlist(hicReg2pks[hiCRegion.slct])))
  
  pks.ovlp =intersect(hic.pks.slct, intersect(rownames(res$diseasedisease), rownames(samp2HeightsMean$depth.GCnormTeng)))
  
  samps.ovlp = intersect(rownames(meta.all), colnames(samp2HeightsMean$depth.GCnormTeng))
  
  
  
  dps = pks.ovlp[res$diseasedisease[pks.ovlp, "P.Value"] <= DP.PVAL.CUTOFF] 
 
  
  # peak2CVs = apply(samp2HeightsMean.depth.GCnormTeng[dps,], 1,
  # FUN=function(x)
  # {
  #   sd(x)/mean(x)
  # })
  # dp.topCVs = names(peak2CVs)[order(peak2CVs, decreasing = T)[1:min(PEAKS.TOP.N, length(peak2CVs))]]
  
  samp2peakHeights.topPks.scaled = t(scale(t(samp2HeightsMean$depth.GCnormTeng[dps,samps.ovlp]), center = T, scale = T ))



  hm.allSamp = matrix(NA, ncol=nrow(meta.all), nrow=length(dps))
  colnames(hm.allSamp) = rownames(meta.all)
  rownames(hm.allSamp) = dps
  
  hm.allSamp[dps, colnames(samp2peakHeights.topPks.scaled)] = samp2peakHeights.topPks.scaled
  
  return(hm.allSamp)
})
names(HM.matrix.list) = HMs

samps.union = levels(factor(unlist(lapply(HM.matrix.list, FUN=function(x){ colnames(x)[!is.na(x[1,])]}))))
HM.matrix.list = lapply(HM.matrix.list, FUN=function(x) x[,samps.union])


# Build the MOFA object
MOFAobject <- createMOFAobject(MultiAssayExperiment(experiments=HM.matrix.list, colData = meta.all[samps.union,]))
MOFAobject


#
DataOptions <- getDefaultDataOptions()
#
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- MOFA.FACTORS.N
#
TrainOptions <- getDefaultTrainOptions()
# Automatically drop factors that explain less than 2% of variance in all omics
TrainOptions$DropFactorThreshold <- MOFA.PROP.VAR.CUTOFF
TrainOptions$seed <- MOFA.SEED
TrainOptions
#

MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

MOFAobject <- regressCovariates(
  object = MOFAobject,
  views = HMs,
  covariates = data.frame(MOFAobject@InputData@colData@listData$Gender, MOFAobject@InputData@colData@listData$age_atsamp)
)
  


MOFAobject <- runMOFA(MOFAobject)

save(MOFAobject,
     samps.union,
     file=file.ou.RData)



#downstream 
pdf(file.ou.pdf, heigh=12, width=15)
# for(hm in HMs)
# {
#   pheatmap(HM.matrix.list[[hm]][, !is.na(HM.matrix.list[[hm]][1,])], main = hm)
# }

plotDataOverview(MOFAobject)

#variance
r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total
plotVarianceExplained(MOFAobject)

f.o.tsv=paste0(dir.ou, "ind.MOFA.DP-pval", DP.PVAL.CUTOFF, ".varProp", MOFA.PROP.VAR.CUTOFF, ".var.tsv")
write.table(r2$R2Total, file=f.o.tsv, sep="\t", quote=F)
write.table(r2$R2PerFactor, file=f.o.tsv, sep="\t", quote=F, append=T)

#
for(hm in HMs)
{
  plotWeightsHeatmap(
    MOFAobject, 
    view = hm,
    #factors = 1:FACTORS.N,
    show_colnames = FALSE,
    main=hm)
}


plotFactorScatters(
  MOFAobject,
  #factors = 1:3,
  color_by = meta.all[samps.union, "group"]
)

annotation_row = data.frame(disease=meta.all[samps.union, "group"])
rownames(annotation_row) = samps.union
pheatmap(MOFAobject@Expectations$Z,
         main="factors for each samples",
         annotation_row = annotation_row)
  

dev.off()
