#V1.4.4
#b.real.estFract_V1.4.4_nnPoisR_noIntercept_6CellType.R
#use nnPois without intercept

#V1.4.3
#row normalized need to work for each individual rather than as bNMF for each row
#so that median of signature peaks for each cell type is the same within one individual


#V1.4.2
#peaks for deconvolution are further filter by their counts in 
#a.selectPeakAndSimulation_V1.5.5_dirichlet.R

#V1.4.1
#also add neutrophil

#V1.4
#use non-negative poisson regression

#V1.3
#use log10 transform

#V1.2
#use only unnormalized cell fraction estimation
#use R-3.3
#use differential peaks


#V1.1
#apply peak markers and method from the b.estFrac.simulation*R
#add checking differential peaks for these marker peaks [should not change healthy samples cell fraction much, but may affect unhealthy ones]
#


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
#library(affy)
library(MASS)
library(addreg)
library("RColorBrewer")
#library(preprocessCore)
#library(EpiDISH)
#library(gcapc)

#source("~/lhou.compbio/codes/mylib/RCode/H3K27ac.functions.R")

args = commandArgs(trailingOnly=TRUE)
HM=args[1]
  
  
HM2PK.NUM = c(H3K27ac=30,
              H3K4me1=30,
              H3K36me3=15)
HM2DIFF.Q= c(H3K27ac=0.01,
              H3K4me1=0.01,
              H3K36me3=0.1)


TYPE2COL=c("grey","red") #, "purple", "green", "blue", "cyan")
names(TYPE2COL)= c("control", "case")#c("IBDRA_Control", "CD", "UC", "RA", "AD_age_matched_control", "AD") #c("IBDRA_Co


GRP2REFs.SLCT = c("CD4T"= "E040|E039|E037|E038",
                  "CD8T"= "E048|E047",
                  "BCell" = "E032",
                  "NKCell" = "E046",
                  "Monocyte" = "E029",
                  "Neutrophil" = "BPNeutrop") #, "E037", "Neutrop")#, "E030")#"E044", "E041", "E042", "E043", "E048", "E034", "E039", "E043", "E045",  #"E047",
GRP2COL = brewer.pal(n = length(GRP2REFs.SLCT), name = "Set2")
names(GRP2COL) = names(GRP2REFs.SLCT)


POIS.DIFF.Q=HM2DIFF.Q[HM]
DIFF.PEAK.COUNT.LB = 10

SAMP.N = 100
PKS.SLCT.REF.NUM=HM2PK.NUM[HM]


#
file.in.RData =paste0("a_selectPeaksAndSimulation_V1.7.2_colinearAuto_rowNorm_6cellType/", HM, "/signaturePeaksAndSimulation.FiltQ", POIS.DIFF.Q, ".above10.RData") #"a_peakSelecti
# file.in.RData =paste0("a_simulation_withRMandBP_V1.5.5_dirichlet/estFractionForPfizer.FiltQ", POIS.DIFF.Q, ".", DIFF.PEAK.COUNT.LB, "-", DIFF.PEAK.COUNT.UB, ".RData")
#file.in.norm.RData = "../untilBatch8.2.2018/a_3_norm_rmBatchEffect/PfizerWithRef.heightsMean.depthGCNormed.rmBE.RData"
file.in.meta.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.goodSamples.meta.RData")

#file.in.testMeta="~/lhou.compbio/data/pfizer_aim2_5.2017/H3K27ac_3.27.2018/sample.metaData.until3.27.2018.txt"
#file.in.DNase = "/broad/compbio/anshul/projects/roadmap/peaks/consolidated/narrowPeak/E000-DNase.macs2.narrowPeak.gz"


dir.ou = paste0("b_real_estFract_V1.4.4.1_nnPoisR_noIntcpt_peakNorm_6cellType_autoSig/", HM, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.real.pdf = paste(dir.ou, "estFraction.results.nnPoisR.FiltQ", POIS.DIFF.Q, ".pdf", sep="")
file.ou.RData=paste(dir.ou, "estFraction.results.nnPoisR.FiltQ", POIS.DIFF.Q, ".RData", sep="")

#

load(file.in.RData) #save(samp2Counts.withRef, grp2pkOvlpWithTest, pk2grpOvlpTest.Weightsum, file=file.ou.RData)
#load(file.in.norm.RData)
load(file.in.meta.RData)



TEST.SAMP.SIZE= nrow(meta.goodsamples)#sum(grepl("^id", colnames(samp2HeightsMean.depth.GCnormTeng), perl=T))

print("real data with differential sites")###########################


peaks.deconv.diff.sampling = sapply(1:SAMP.N,
FUN=function(i)
{
  unlist(lapply(names(GRP2REFs.SLCT),
  FUN=function(grp)
  {
    #print(grp)
    sample(peaks.inGrp.diff.deconv.list[[grp]], PKS.SLCT.REF.NUM)
  }))
})

fml = as.formula(paste0("y~", paste0(names(GRP2REFs.SLCT), collapse="+"), "-1"))
Grp.slct.props.sampling.unNorm=lapply(1:ncol(peaks.deconv.diff.sampling), 
FUN=function(i)
{
  print(i)
  #pks=peaks.deconv[[nm]]
  pks=peaks.deconv.diff.sampling[,i]
  test.pks =samp2HeightsMean.depth.GCnormTeng[pks,c(1:TEST.SAMP.SIZE)]
  grp.pks = grp.counts[pks,]
  estF = t(apply(test.pks, 2,
  FUN=function(y)
  {
    #median of the peak from the same cell type keep same
    #print(y[1])
    y.matrix= matrix(y, nrow=PKS.SLCT.REF.NUM)
    y.cellType.mean = apply(y.matrix, 2, mean)
    y.cellType.mean.median=median(y.cellType.mean)
    cellFract.init = y.cellType.mean/sum(y.cellType.mean)
    y.cellType.mean.factor=rep(y.cellType.mean.median/y.cellType.mean, each=PKS.SLCT.REF.NUM)
    y.norm=y*y.cellType.mean.factor
    grp.pks.norm = sweep(grp.pks, 1, y.cellType.mean.factor, "*" )
    
    mat.mod= model.matrix(fml, data.frame(y=round(y.norm), grp.pks.norm))
    #mod=addreg(fml, data = data.frame(y=round(y.norm), grp.pks.norm), mono=names(GRP2REFs.SLCT), family="poisson", method="em", accelerate="squarem")  
    #mod=addreg(fml, data = df.mod, mono=names(GRP2REFs.SLCT), family="poisson", method="em", accelerate="squarem")  
    mod = nnpois(round(y.norm), mat.mod, standard=1, offset=0, start = cellFract.init)
    #mod=addreg(fml, data = df.mod, mono=names(GRP2REFs.SLCT), family="poisson", method="em", accelerate="squarem")  
    
    # rlm.o <- rlm(log10(y+1) ~ log10(grp.pks+1), maxit = 150)  
    # coef.v <- summary(rlm.o)$coef[2:(ncol(grp.pks) + 1), 1]
    # coef.v[which(coef.v < 0)] <- 0
    
    return(mod$coef)
  }))
  colnames(estF) = colnames(grp.pks)
  rownames(estF) = colnames(test.pks)

  return(estF)
})


Grp.slct.props.sampling= lapply(Grp.slct.props.sampling.unNorm,
FUN=function(estF)
{
  t(apply(estF, 1, 
  FUN=function(x)
  {
    x/sum(x)
  }))
})


Grp.slct.props.unNorm.mean = sapply(names(GRP2REFs.SLCT),
FUN=function(grp)
{
  grp.props= sapply(Grp.slct.props.sampling.unNorm, FUN=function(x) x[,grp])
  return(apply(grp.props, 1, mean))
})



Grp.slct.props.mean = sapply(names(GRP2REFs.SLCT),
FUN=function(grp)
{
  grp.props= sapply(Grp.slct.props.sampling, FUN=function(x) x[,grp])
  return(apply(grp.props, 1, mean))
})




print("visualization")##############################

conditions=list()
conditions$all =  c("control", "case")



pdf(file.ou.real.pdf, width=12, height=6)

#mean estimate fraction
layout(matrix(1:length(dis2samps), ncol=1))
par(mar=c(4,7,5,5))
for(dis in conditions$all)
{
  samps=dis2samps[[dis]]
  #order by neutrophil
  samps.ordered=order(Grp.slct.props.mean[samps,"Neutrophil"])

  f.o.tsv=paste(dir.ou, "estFraction.results.nnPoisR.FiltQ", POIS.DIFF.Q, ".", dis, ".tsv", sep="")
  write.table(t(Grp.slct.props.mean[samps[samps.ordered],]), file=f.o.tsv, sep="\t", quote=F)

  barplot(t(Grp.slct.props.mean[samps[samps.ordered],]), col=GRP2COL, las=2,
          xlim=c(0,130),
          ylab="cell fraction estimated",
          main= dis,
          legend.text = paste(colnames(Grp.slct.props.mean)),
          axisnames=F,
          cex.main=3,
          cex.axis=1.5,
          cex.lab=2)
          #cex.names=0.3)

}


#Grp.slct.props.mean.rmNeutr=Grp.slct.props.mean[,1:7]/apply(Grp.slct.props.mean[,1:7], 1, sum)

lapply(conditions,
FUN=function(conds)
{
  layout(matrix(1:6, ncol=3))
  par(mar=c(6,10,5,5))
  
  samps= unlist(dis2samps[conds])
  for(j in 1:ncol(grp.counts))
  {
    # if(REF2CELLTYPE[colnames(grp.counts)[j]]=="Neutrophil")
    # {
    #   boxplot(Grp.slct.props.mean[,j]~factor(meta.goodsamples[rownames(Grp.slct.props.mean),"Phenotype"],levels=names(TYPE2COL)),
    #         xlab="conditions", ylab="estimated fraction",
    #         main=paste(colnames(Grp.slct.props.mean)[j], REF2CELLTYPE[colnames(Grp.slct.props.mean)[j]], sep="\n"),
    #         cex.main=2,
    #         cex.axis =1.5,
    #         cex.lab = 2,
    #         outline=F,
    #         col= TYPE2COL[names(TYPE2COL)])#levels(factor(meta.goodsamples[rownames(Grp.slct.props.mean), "Phenotype"]))])
    # 
    # 
    #   stripchart(Grp.slct.props.mean[,j]~factor(meta.goodsamples[rownames(Grp.slct.props.mean), "Phenotype"],levels=names(TYPE2COL)),
    #            vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black', cex=2)
    # 
    # }else
    # {
      mod=lm(Grp.slct.props.mean[samps,j]~factor(meta.goodsamples[samps,"Phenotype"]))
      boxplot(Grp.slct.props.mean[samps,j]~factor(meta.goodsamples[samps,"Phenotype"],levels=conds),
            xlab="conditions", ylab="estimated fraction in buffy coat\n",#[after removal of neutrophil]
            main=paste0(colnames(Grp.slct.props.mean)[j], " P:", signif(summary(mod)$coef[2, "Pr(>|t|)"],2)),
            cex.main=2,
            cex.axis =1.5,
            cex.lab = 2,
            outline=F,
            col= TYPE2COL[1:length(conds)]) #col= TYPE2COL[levels(factor(meta.goodsamples[rownames(Grp.slct.props.mean.rmNeutr), "Phenotype"]))])
      stripchart(Grp.slct.props.mean[samps,j]~factor(meta.goodsamples[samps,"Phenotype"],levels=conds),
               vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black', cex=1)
    # }
    
  }
})

# #after removing diff peak
# layout(matrix(1:length(dis2samps), ncol=1))
# for(dis in c("control", "case"))
# {
# 
#   barplot(t(Grp.slct.props.sampling.rmDP.case.mean[dis2samps[[dis]],]), col=GRP2COL, las=2,
#           xlim=c(0,40),
#           main= dis,
#           legend.text = paste(colnames(Grp.slct.props.sampling.rmDP.case.mean), REF2CELLTYPE[colnames(Grp.slct.props.sampling.rmDP.case.mean)]),
#           cex.names=0.3)
# 
# }
# 
# layout(matrix(1:9, ncol=3))
# for(j in 1:ncol(grp.counts))
# {
#   boxplot(Grp.slct.props.sampling.rmDP.case.mean[,j]~factor(samp2meta[rownames(Grp.slct.props.sampling.rmDP.case.mean),"phenotype"],levels=names(TYPE2COL)),
#           xlab="conditions", ylab="estimated fraction",
#           main=paste(colnames(Grp.slct.props.sampling.rmDP.case.mean)[j], REF2CELLTYPE[colnames(Grp.slct.props.sampling.rmDP.case.mean)[j]], sep="\n"),
#           col= TYPE2COL[levels(factor(samp2meta[rownames(Grp.slct.props.sampling.rmDP.case.mean), "phenotype"]))])
# 
# 
#   stripchart(Grp.slct.props.sampling.rmDP.case.mean[,j]~factor(samp2meta[rownames(Grp.slct.props.sampling.rmDP.case.mean), "phenotype"],levels=names(TYPE2COL)),
#              vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black', cex=2)
# 
# }

#RMSE
layout(matrix(1:6, ncol=2))
Grp.slct.props.RMSEs = lapply(names(dis2samps),
FUN=function(dis)
{
  samps = dis2samps[[dis]]
  RMSEs = sapply(names(GRP2REFs.SLCT),
  FUN=function(grp)
  {
    grps.props = sapply(Grp.slct.props.sampling,
    FUN=function(props)
    {
      props[samps, grp]
    })
    
    grps.props.mean =Grp.slct.props.mean[samps,grp]
    
    rmse = apply(grps.props,2,
    FUN=function(x)
    {
      (mean((x-grps.props.mean)^2))^0.5
    })
    
    
  })
  boxplot(RMSEs, main=dis, ylim=c(0,0.5), col=TYPE2COL)
  return(RMSEs)
})
names(Grp.slct.props.RMSEs) = names(dis2samps)

# #differential peaks 1st round
# layout(matrix(1:1, ncol=1))
# boxplot(-log10(diffPeaks.pvals), xlab="different samplings", ylab="differential Peak -log10 pvalue", main="before rm Diff. Peaks")

dev.off()


#samp2HeightsMean.depth.GCnormTeng.rmBE=samp2HeightsMean$depth.GCnormTeng.rmBE
save(meta.goodsamples, 
     samp2HeightsMean.depth.GCnormTeng, 
     grp.counts,
     peaks.deconv.diff.sampling,
     Grp.slct.props.sampling.unNorm, 
     Grp.slct.props.sampling, 
     Grp.slct.props.unNorm.mean,
     Grp.slct.props.mean,
     #Grp.slct.props.sampling.rmDP.case, Grp.slct.props.sampling.rmDP.case.mean, diffPeaks.pvals.rmDP.case,
     file=file.ou.RData)




#
file.ou.sum.pdf = paste0(dir.ou, "cellFraction.sum.FiltQ", POIS.DIFF.Q, ".pdf")

fract.sum =t(sapply(Grp.slct.props.sampling.unNorm, FUN=function(x){apply(x,1,sum)}))
dis2samp = lapply(split(rownames(meta.goodsamples), as.character(meta.goodsamples$Phenotype)),
FUN=function(samps)
{
  samps[order(apply(fract.sum[,samps],2,mean))]
})
pdf(file.ou.sum.pdf, width=30, height = 10)
samples.dis = unlist(dis2samp[names(TYPE2COL)])
boxplot(fract.sum[,samples.dis],
        ylab="sum of fraction",
        xlab="samples",
        col = TYPE2COL[as.character(meta.goodsamples[samples.dis, "Phenotype"])],
        #col=c("grey","red")[as.logical(qc.info[colnames(fract.sum),"GoodQuality"])+1], 
        las=2)
abline(h=1, lty=2, lwd=2)
legend("topright", legend = names(TYPE2COL), fill=TYPE2COL)

dev.off()




#
file.ou.cor.pdf = paste(dir.ou, "cellFraction.cor.FiltQ", POIS.DIFF.Q, ".pdf", sep="")
pdf(file.ou.cor.pdf, width=20, height = 10)
for(dis in names(TYPE2COL))
{
  samps = rownames(meta.goodsamples)[meta.goodsamples$Phenotype==dis]
  layout(matrix(1:28, ncol=7))
  print(length(samps))
  for(i in 1:(length(names(GRP2REFs.SLCT))-1))
  {
    grp.i=names(GRP2REFs.SLCT)[i]
    for(j in (i+1):length(names(GRP2REFs.SLCT)))
    {
      grp.j = names(GRP2REFs.SLCT)[j]
      pcc = cor(Grp.slct.props.mean[samps,grp.i], Grp.slct.props.mean[samps,grp.j])
      plot(x = Grp.slct.props.mean[samps,grp.i],
           y = Grp.slct.props.mean[samps,grp.j],
           xlab= paste(grp.i),
           ylab= paste(grp.j),
           main=paste(dis, " pcc=", signif(pcc,2), sep=""))
    }
  }
}
dev.off()





#estimate of mixture for each different cell type

file.ou.mixture.pdf = paste(dir.ou, "mixtureEst.FiltQ", POIS.DIFF.Q, ".pdf", sep="")
pdf(file.ou.mixture.pdf, width=15, height = 8)
layout(matrix(1:12, ncol=4))
for(i in 1: ncol(peaks.deconv.diff.sampling))
{
  
  pks =  peaks.deconv.diff.sampling[,i]
  
  estR = grp.counts.rowNorm[pks,]
  estF = Grp.slct.props.sampling.unNorm[[i]]
  estMix = estR %*% t(estF)
  
  
  for(grp in names(GRP2REFs.SLCT))
  {
    grp.peaks = intersect(pks, peaks.inGrp.diff.deconv.list[[grp]])
    plot(samp2HeightsMean.depth.GCnormTeng.rowNorm[grp.peaks, ],
         estMix[grp.peaks, ],
         xlab="mixture counts from test samples",
         ylab="mixture counts from estimation",
         col=GRP2COL[grp],
         main=paste(i, grp))
    abline(a=0, b=1, lty=2)
  }
}

dev.off()
