#V1.4
#a.2.normalizePeakCounts.V1.4_other.R
#peak heights calc. are moved to the previous step
#to generally normalized total reads in peaks around 5,000,000 before use size factor to normalize


#V1.3 with reference
#overlap with reference are quantified and visualzied
#peak counts are extraced from bigwig by bigWigAverageOverBed
#since blue print bw is from different process, they are not included in calculation of geoMeans for size factor 

#V1.2
#add reads length normalization after depth normalization
#normalize to 1kb

#V1.1 to better evaluate the normalization 
#samples from CD4+ naive and monocyte form blueprint are also included
#add depth normalization, and all the following normalization is based on it

#normalize test and reference samples together 
#by different means and test consistence across replicates



args = commandArgs(trailingOnly=TRUE)
HM = args[1]





R.LIBS.MYPATHs = c(V3.5="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.5/",
                V3.4="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.4/",
                V3.3="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.3/")
R.LIBS.PATH = .libPaths()

if(!grepl("R version 3.5", R.version$version.string))
{
  stop("should use R-3.5") 
}
.libPaths(c(R.LIBS.MYPATHs["V3.5"], R.LIBS.PATH))
 
# TYPE2COL=c("grey","red", "purple", "green")
# names(TYPE2COL)=c("IBDRA_Control", "CD", "UC", "RA")

GOOD.TOTALREADs.CUTOFF=2e7
GOOD.RSC.CUTOFF=1

TYPE2PCH=c(4, 2)
names(TYPE2PCH)=c("case", "control")

REFs.SELECT= c("E040","E039", "E048", "E047", "E029", "E032", "E046", "E030", "E037", "E038")#"BPNeutrop", "E037", "Eneutrop")#, "E030")#"E044", "E041", "E042", "E043", "E048", "E034", "E039", "E043", "E045", 
# REF2COL = rainbow(length(REFs.SELECT))
# names(REF2COL) = REFs.SELECT

REF2CELLTYPE= 
  c(
    E062="Primary mononuclear cells from peripheral blood",
    E034="Primary T cells from peripheral blood",
    E045="Primary T cells effector/memory enriched from peripheral blood",
    #E033="Primary T cells from cord blood",
    E044="Primary T regulatory cells",# from peripheral blood",
    E043="Primary T helper cells",# from peripheral blood",
    E041="Primary T helper cells PMA-I stimulated",
    E042="Primary T helper 17 cells PMA-I stimulated",
    E040="Primary T helper memory cells1",#  from peripheral blood 1",
    E037="Primary T helper memory cells2",# from peripheral blood 2",
    E038="Primary T helper naive cells",#  from peripheral blood",
    E039="Primary T helper naive cells",# from peripheral blood",
    E048="Primary T CD8+ memory cells",#  from peripheral blood",
    E047="Primary T CD8+ naive cells",#  from peripheral blood",   
    E029="Primary monocytes",#  from peripheral blood",
    #E031="Primary B cells from cord blood",
    # E035="Primary hematopoietic stem cells",
    # E051="Primary hematopoietic stem cells G-CSF-mobilized Male",
    #E050="Primary hematopoietic stem cells G-CSF-mobilized Female",
    # E036="Primary hematopoietic stem cells short term culture",
    E032="Primary B cells",#  from peripheral blood",
    E046="Primary Natural Killer cells",#,#  from peripheral blood",
    E030="Primary neutrophils",#  from peripheral blood"
    BPNeutrop="neutrophil"
    )

blueprint.nm2celltypes = c("mature_neutrophil", "CD4-positive_alpha-beta_T_cell", "CD14-positive_CD16-negative_classical_monocyte")
names(blueprint.nm2celltypes) = c("BP-Neutrophil", "BP-CD4NaiveT", "BP-monocyte")

#
library(affy)
library(preprocessCore)
library(MASS)
library(splines)
library(matrixStats)
library(DESeq2)
library(BSgenome)
#library(EpiDISH)
#library(gcapc)
library(pheatmap)
# library("RColorBrewer")
# library("devtools")
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

source("~/lhou.compbio/codes/mylib/RCode/H3K27ac.functions.R")


#
file.in.RData = paste0("a_mergePeak2readsCounts_V4.2_filteredMergedPeaks/", HM, ".hg19.narrowPeak.RSC2.Tot3e+07.logQ2.merged.filt2GoodSamps.heightsMean.RData")

file.in.mergedPeak = paste0("a_mergePeak2readsCounts_V4.2_filteredMergedPeaks/", HM, ".hg19.narrowPeak.RSC2.Tot3e+07.logQ2.merged.filt2GoodSamps.bed")

#file.in.testMeta="~/lhou.compbio/data/Mayo_aim2_5.2017/H3K27ac_3.27.2018/sample.metaData.until3.27.2018.txt"


dir.ou = "a_2_peakHeightsNormalized_V1.4/"
if(!dir.exists(dir.ou))
{
  dir.create(dir.ou)
}

file.ou.GCTeng.pdf = paste(dir.ou, HM, ".diagnose.GCTeng.pdf", sep="")
file.ou.pdf = paste(dir.ou, HM, ".MayoWithRef.heightsMean.depthGCNormed.pdf", sep="")
file.ou.RData = paste(dir.ou, HM, ".MayoWithRef.heightsMean.depthGCNormed.RData", sep="")
file.ou.meta.RData = paste(dir.ou, HM, ".MayoWithRef.goodSamples.meta.RData", sep="")
file.ou.ref2pkOvlprefSlct.combBP= paste0(dir.ou, HM, ".ref2pkOvlpWithTest.refSlct.combBP.RDS")

#
load(file.in.RData)

meta.goodsamples = qc.info[qc.info$totalReads>=GOOD.TOTALREADs.CUTOFF & qc.info$RSC>=GOOD.RSC.CUTOFF & qc.info$isBad.byMismatch==0,]
#samp2meta = read.table(file.in.testMeta, sep="\t", header=T, row.names=1)
dis2samps = lapply(names(TYPE2PCH),
FUN=function(dis)
{
  rownames(meta.goodsamples)[meta.goodsamples$Phenotype==dis]
})
names(dis2samps)=names(TYPE2PCH)

save(meta.goodsamples, dis2samps, file=file.ou.meta.RData)


#to valdiate the merged.peaks.inRef makes sense
# for(i in 1:60)
# {
#   print(i)
#   print(paste(sum(samp2HeightsMean$rawCounts[merged.peaks.inRef,i]>=10)/length(merged.peaks.inRef),
#   sum(samp2HeightsMean$rawCounts[setdiff(merged.peaks, merged.peaks.inRef),i]>=10)/(length(merged.peaks)-length(merged.peaks.inRef))))
# }


#samp2HeightsMean$ROADMAP= samp2HeightsMean$ROADMAP[,!is.na(samp2HeightsMean$ROADMAP[1,])]



###################################################
print("normalize for library sequencing depth")
merged.peaks = names(merged.peaks.len)
merged.peaks.inRef = merged.peaks[pk2refOvlpTest.Weightsum>=0.5]

samp2HeightsMean.upto30M = apply(cbind(samp2HeightsMean$Mayo[merged.peaks.inRef,], samp2HeightsMean$ROADMAP[merged.peaks.inRef,]),2,
FUN=function(x)
{
  total.reads = rbind(x) %*% cbind(merged.peaks.len[merged.peaks.inRef])/200
  x = x* 3*10^7/total.reads 
  #x = x* 5*10^6/(sum(x)*500/200) #500 is average peak length, 200 is average fragment length
})


samp2HeightsMean$rawCounts = cbind(samp2HeightsMean.upto30M, samp2HeightsMean$blueprint[merged.peaks.inRef,])


samp2HeightsMean.NOBlueprint.geoMean = apply(samp2HeightsMean.upto30M, 1,
FUN=function(x)
{
  exp(sum(log(x[x > 0])) / length(x))
})
samp2HeightsMean.sizeFactor = estimateSizeFactorsForMatrix(samp2HeightsMean$rawCounts, 
                                                              geoMeans =samp2HeightsMean.NOBlueprint.geoMean)
samp2HeightsMean$depthNorm = t(t(samp2HeightsMean$rawCounts)/samp2HeightsMean.sizeFactor)

# ###################################################
# print("normalize for library peak length")
# 
# samp2HeightsMean$depth.pkLen.Norm = samp2HeightsMean$depthNorm/merged.peaks.len[rownames(samp2HeightsMean$depthNorm)]*1000

# ##################################################
# print("GC correction based on reads counts for each bin by myself")
# 
# peak2GC = read.table(file.in.pk2GC, sep="\t", row.names=NULL, header=F)[,1]
# names(peak2GC) = merged.peaks
# 
# samp2HeightsMean$depthNorm.quantile = normalize.quantiles(samp2HeightsMean$depthNorm)
# colnames(samp2HeightsMean$depthNorm.quantile)=colnames(samp2HeightsMean$depthNorm)
# rownames(samp2HeightsMean$depthNorm.quantile)=rownames(samp2HeightsMean$depthNorm)
# 
#   
# samp2HeightsMean$depthNorm.GCnormMyself = getGCnormalizedPeakHeight(peak2GC[merged.peaks.inRef], samp2HeightsMean$depthNorm, GC.binNum)
# save(samp2HeightsMean, ref2pkOvlpWithTest, pk2refOvlpTest.Weightsum, merged.peaks.inRef, merged.peaks, file=file.ou.RData)
# 
# 
# samp2HeightsMean$depthNorm.GCnormMyself.quantile = normalize.quantiles(samp2HeightsMean$depthNorm.GCnormMyself)
# colnames(samp2HeightsMean$depthNorm.GCnormMyself.quantile)=colnames(samp2HeightsMean$depthNorm.GCnormMyself)
# rownames(samp2HeightsMean$depthNorm.GCnormMyself.quantile)=rownames(samp2HeightsMean$depthNorm.GCnormMyself)

##################################################
print("GC correction based on mixture model adapted from Mingxiang Teng")

pdf(file.ou.GCTeng.pdf)
samp2HeightsMean$depth.GCnormTeng=refineSites(round(samp2HeightsMean$depthNorm[merged.peaks.inRef,]), flank=0, 
                                           sites=makeGRangesFromDataFrame(merged.peaks.df[merged.peaks.inRef,], 
                                           seqnames.field="chr", start.field="start", end.field="end"),
                                           gcrange=c(0,1),
                                           model="poisson", gctype="uniform", emtrace=T, mu1=10)
rownames(samp2HeightsMean$depth.GCnormTeng)=merged.peaks.inRef
dev.off()




# samp2HeightsMean$depth.GCnormTeng.quant = normalize.quantiles(samp2HeightsMean$depth.GCnormTeng)
# colnames(samp2HeightsMean$depth.GCnormTeng.quant)=colnames(samp2HeightsMean$depth.GCnormTeng)
# rownames(samp2HeightsMean$depth.GCnormTeng.quant)=rownames(samp2HeightsMean$depth.GCnormTeng)

# samp2HeightsMean$GCnormTeng.quantSeparate = cbind(normalize.quantiles(samp2HeightsMean$GCnormTeng[,!grepl("^E|^BP",colnames(samp2HeightsMean$GCnormTeng))]),
#                                                       normalize.quantiles(samp2HeightsMean$GCnormTeng[,grepl("^E|^BP",colnames(samp2HeightsMean$GCnormTeng))]))
# colnames(samp2HeightsMean$GCnormTeng.quantSeparate)=c(colnames(samp2HeightsMean$GCnormMyself)[!grepl("^E|^BP",colnames(samp2HeightsMean$GCnormTeng))],
#                                                            colnames(samp2HeightsMean$GCnormMyself)[grepl("^E|^BP",colnames(samp2HeightsMean$GCnormTeng))])
# rownames(samp2HeightsMean$GCnormTeng.quantSeparate)=rownames(samp2HeightsMean$GCnormTeng)


samp2HeightsMean$Mayo=NULL
samp2HeightsMean$ROADMAP=NULL
samp2HeightsMean$blueprint=NULL
save(samp2HeightsMean,  merged.peaks.inRef, merged.peaks, merged.peaks.len, file=file.ou.RData)
print("overall sense of how normalization worked")##################################################


# merged.peaks.cellSpecif1 = merged.peaks.inRef[pk2refOvlpTest.Weightsum[merged.peaks.inRef]==1]
# merged.peaks.cellSpecif2 = merged.peaks.inRef[pk2refOvlpTest.Weightsum[merged.peaks.inRef]<=5 & pk2refOvlpTest.Weightsum[merged.peaks.inRef]>=2]
#colnames(samp2HeightsMean$depth.GCnormTeng) = gsub("BP-", "BP", colnames(samp2HeightsMean$depth.GCnormTeng) )
colnames(ref2pkOvlpWithTest) = gsub("BP-", "BP",colnames(ref2pkOvlpWithTest) )

ref2pkOvlpWithTest.refSlct = ref2pkOvlpWithTest[,grepl(paste(REFs.SELECT, collapse="|"), colnames(ref2pkOvlpWithTest))]
weights.refSlct = rep(1, ncol(ref2pkOvlpWithTest.refSlct))
names(weights.refSlct) = colnames(ref2pkOvlpWithTest.refSlct)
weights.refSlct[grepl(paste("^BPNeutrophil", sep=""), colnames(ref2pkOvlpWithTest.refSlct))] = 1/sum(grepl(paste("^BPNeutrophil", sep=""), colnames(ref2pkOvlpWithTest.refSlct))) #from the same cell type, share the same weights.refSlct  


ref2pkOvlpWithTest.refSlct.combBP = sapply(REFs.SELECT,
FUN=function(ref)
{
  if(ref %in% colnames(ref2pkOvlpWithTest.refSlct))
  {
    return(ref2pkOvlpWithTest.refSlct[,ref])
  }else
  {
    samps = colnames(ref2pkOvlpWithTest.refSlct)
    samps = samps[grepl(ref, samps)]
    apply(ref2pkOvlpWithTest.refSlct[,samps], 1, mean)
  }
  
})
saveRDS(ref2pkOvlpWithTest.refSlct.combBP, 
        file= file.ou.ref2pkOvlprefSlct.combBP)


pk2refOvlpTest.refSlct.Weightsum = apply(ref2pkOvlpWithTest.refSlct,1,
FUN=function(x)
{
  sum(weights.refSlct[x==1])
})
#pk2refOvlpTestNum.refSlct = apply(ref2pkOvlpWithTest[,REFs.SELECT[-8]],1,sum)
peaks.cellSpecif1 = names(pk2refOvlpTest.refSlct.Weightsum)[pk2refOvlpTest.refSlct.Weightsum>=1]
peaks.cellSpecif14 = names(pk2refOvlpTest.refSlct.Weightsum)[pk2refOvlpTest.refSlct.Weightsum>=1 &pk2refOvlpTest.refSlct.Weightsum<=4]

peaks.cellSpecif =  peaks.cellSpecif14

save(samp2HeightsMean, ref2pkOvlpWithTest.refSlct, pk2refOvlpTest.refSlct.Weightsum, merged.peaks.inRef, merged.peaks, merged.peaks.len, file=file.ou.RData)



pdf(file.ou.pdf, width=20, height=12)
layout(matrix(1:2, ncol=1))
for(i in 1:length(samp2HeightsMean))
{
  print(i)
  
  samples = colnames(samp2HeightsMean[[i]])
  samp2type = sapply(samples,
  FUN=function(s)
  {
    if(grepl("^E\\d+$", s))
    {
      return("ROADMAP")
    }else if(grepl("^BP-Neutrophil", s))
    {
      return("BP-Neutrophil")
    }else if(grepl("^BP-CD4NaiveT", s))
    {
      return("BP-CD4NaiveT")
    }else if(grepl("^BP-monocyte", s))
    {
      return("BP-monocyte")
    }else
    {
      return("Mayo")
    }
  })
  type2col = c("red", "blue", "purple", "orange", "grey")
  names(type2col) = c("ROADMAP", "BP-Neutrophil", "BP-CD4NaiveT", "BP-monocyte", "Mayo")
  
  samp2pch = rep(4, length(samples))
  names(samp2pch)=samples
  for(dis in names(dis2samps))
  {
    samp2pch[dis2samps[[dis]]] = TYPE2PCH[dis]
  }

  counts.log.scaled = scale(t(log10(samp2HeightsMean[[i]][peaks.cellSpecif,]+1)))
  plot(hclust(as.dist(1-cor(t(counts.log.scaled)))), main= paste("log10", names(samp2HeightsMean)[i]), cex=0.7)
  
  pca = prcomp(counts.log.scaled)  #merged.peaks.inRef, merged.peaks.cellSpecif
  pca.imptnce = summary(pca)$importance
  plot(pca$x[,"PC1"], pca$x[,"PC2"], #xlab="PC1", ylab="PC2",
       xlab=paste("PC1: ", pca.imptnce["Proportion of Variance","PC1"]*100, "%", sep=""),
       ylab=paste("PC2: ", pca.imptnce["Proportion of Variance","PC2"]*100, "%", sep=""),
       main= paste("log10", names(samp2HeightsMean)[i]),
       col=type2col[samp2type],
       pch=samp2pch
       )
  y.range = max(pca$x[,"PC2"])-min(pca$x[,"PC2"])
  text(x=pca$x[samp2type=="ROADMAP","PC1"], 
       y=pca$x[samp2type=="ROADMAP","PC2"]+y.range*0.02,
       cex=0.6,
       labels = samples[samp2type=="ROADMAP"])
  legend("bottomleft", legend = c(names(type2col), names(dis2samps)), 
         col=c(type2col,rep("grey", length(dis2samps))), 
         pch=c(rep(4,length(type2col)), TYPE2PCH[names(dis2samps)])
         )
}

# pk2refOvlpTestNum.refSlct = apply(ref2pkOvlpWithTest[,REFs.SELECT],1,sum)
# peaks.cellSpecif = names(pk2refOvlpTestNum.refSlct)[pk2refOvlpTestNum.refSlct==1]
# prctage.testLessThanRef=lapply(samp2HeightsMean,
# FUN=function(samp2HeightsMean)
# {
#   peaks.cellSpecif.max = apply(samp2HeightsMean[peaks.cellSpecif,REFs.SELECT], 1,max)
#   sapply(1:60,
#   FUN=function(i){
#     sum(samp2HeightsMean[peaks.cellSpecif,i] <= peaks.cellSpecif.max)/length(peaks.cellSpecif)
#   })
# 
# })
# boxplot(prctage.testLessThanRef, ylab="percentage of cell specific peaks in test samples with counts<= max in reference sample", las=2)


dev.off()



pdf(paste(dir.ou, HM, ".overallPropertyOfPeaksAfterNorm.pdf",sep=""), width=12, height=20)
merged.peaks.anyRef = merged.peaks[pk2refOvlpTest.refSlct.Weightsum>=0.5]

boxplot(merged.peaks.len[merged.peaks.anyRef]~as.factor(round(pk2refOvlpTest.refSlct.Weightsum[merged.peaks.anyRef])),
        xlab="presence in Ref samples",ylab="length of peak")

layout(matrix(1:14,ncol=2, byrow=T))
for(ref in REFs.SELECT)
{
  merged.peaks.ref = merged.peaks[ref2pkOvlpWithTest.refSlct.combBP[,ref]>=0.5]
  
  
  if(ref != "BPNeutrop")
  {
    plot(merged.peaks.len[merged.peaks.ref], samp2HeightsMean$depth.GCnormTeng[merged.peaks.ref, ref]+1, log="y", cex=0.2,
         xlab="length of peak", ylab="peak mean heights", main=ref)  
    boxplot((samp2HeightsMean$depth.GCnormTeng[merged.peaks.ref, ref]+1)~as.factor(round(pk2refOvlpTest.refSlct.Weightsum[merged.peaks.ref])), log="y",
         xlab="presence in Ref samples", ylab="peak mean heights", main=ref)

   # boxplot((samp2HeightsMean$depth.pkLen.GC.Norm[merged.peaks.ref, ref]+1)/merged.peaks.len[merged.peaks.ref]*1000~as.factor(pk2refOvlpTestNum.refSlct[merged.peaks.ref]), log="y",
   #     xlab="occurence in Ref samples", ylab="normalized counts RPK", main=ref)
  }else
  {
    merged.peaks.ref.counts = apply(samp2HeightsMean$depth.GCnormTeng[merged.peaks.ref, grepl("BP-Neutrophil", colnames(samp2HeightsMean$depth.GCnormTeng))],1,mean)
    plot(merged.peaks.len[merged.peaks.ref], 
         merged.peaks.ref.counts+1, log="y", cex=0.2,
         xlab="length of peak", ylab="peak mean heights", main=ref)  
    boxplot((merged.peaks.ref.counts+1)~as.factor(round(pk2refOvlpTest.refSlct.Weightsum[merged.peaks.ref])), log="y",
         xlab="presence in Ref samples", ylab="peak mean heights", main=ref)
  }
}


dev.off()

#
