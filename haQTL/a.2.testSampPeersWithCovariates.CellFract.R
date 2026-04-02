#for the same overlapped peak
#check correlation from two tissues across same individuals 
#partial correlation by sex

args = commandArgs(trailingOnly=TRUE)
HM = args[1]
#TISS=gsub("(.*?)\\..*", "\\1", COND)

#
library(pheatmap)
library(RColorBrewer)
#cols <- brewer.pal(9, 'YlOrRd')
HM2DIFF.Q= c(H3K27ac=0.01,
             H3K4me1=0.01,
             H3K36me3=0.1)
POIS.DIFF.Q=HM2DIFF.Q[HM]


CellTYPEs =c("CD4T", "CD8T", "BCell", "NKCell", "Monocyte", "Neutrophil" )

colfunc <- colorRampPalette(c("grey99", "darkred"))
cols <- colfunc(10)
# library(preprocessCore)
file.in.deconv.RData =  paste0("../deconv/b_real_estFract_V1.4.4.1_nnPoisR_noIntcpt_peakNorm_6cellType_autoSig/", HM, "/estFraction.results.nnPoisR.FiltQ", POIS.DIFF.Q, ".RData")
#file.in.mergedPk2TissPk.RData="../peakMergeNormal/c_mergePeaksV1.3FromTissues_hg19/4Tiss.mergedPeak2tissuePeak.hg19.RData"
file.in.peer= paste0("a_haQTLCalling_covariate_V1.2.3_pkFiltByMedianCV_gINT_Peer/", HM, ".top50.PEER_covariates.txt")
file.in.meta.csv = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
file.in.wgs.sampleInfo="~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"

file.in.hm.meta = paste0("~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/c_metaAndQC/", HM, ".sample.metaDataAndQC.txt")
  
dir.ou="a_2_PeersByKnownCovarAndCellFractions/"
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.overall.pdf = paste0(dir.ou, HM, ".VarByKnownCovar.pdf")
file.ou.perPair.pdf = paste0(dir.ou, HM, ".VarByKnownCovar.perPair.pdf")
file.ou.varProp.pdf =  paste0(dir.ou, HM, ".VarProportion.pdf")
file.ou.RData = paste0(dir.ou, HM, ".VarByKnownCovar.RData")

#

indID2info = read.table(file.in.wgs.sampleInfo, sep=",", header=T, row.names = 1, comment.char = "", quote="", stringsAsFactors = F)
IndID2sampID  = paste0("id_", rownames(indID2info))
names(IndID2sampID)= indID2info$HA.Sample.Name

load(file.in.deconv.RData)
peers = as.matrix(read.table(file.in.peer, sep="\t", header=T, row.names = 1, check.names = F))
colnames(peers) = IndID2sampID[colnames(peers)]

factors.est.list = list(peers = t(peers))

#
meta.samp= read.table(file.in.meta.csv, sep=",", header=T, row.names = 1, comment.char = "", quote="", stringsAsFactors = F)
rownames(meta.samp) = paste0("id_", rownames(meta.samp))

meta.HM= read.table(file.in.hm.meta, sep="\t", header=T, row.names = 1, comment.char = "", quote="", stringsAsFactors = F)


hm.samps = intersect(colnames(peers), rownames(Grp.slct.props.unNorm.mean))

hm.knownCovar = data.frame(meta.samp[hm.samps, c("Gender", "age_atsamp")],
                        meta.HM[hm.samps, c("FRiP", "NSC", "RSC", "totalReads", "Phenotype")],
                        Grp.slct.props.unNorm.mean[hm.samps,],
                        stringsAsFactors = F   )

  
save(
     factors.est.list,
     hm.knownCovar,
     file=file.ou.RData) 
#

pdf(file.ou.perPair.pdf, width=15, height=12)

R2.PCvsKnownCovar= lapply(names(factors.est.list),
FUN=function(cond)
{
  layout(matrix(1:9, ncol=3))
  
  
  factors.est = factors.est.list[[cond]]
  res=sapply(colnames(factors.est),
  FUN=function(pc.id)
  {
    pc =factors.est[,pc.id]
    sapply(names(hm.knownCovar), 
    FUN=function(nm)
    {
      #print(nm)
      covar= hm.knownCovar[hm.samps, nm]
      mod = lm(pc~covar)
      r2=summary(mod)$r.squared
      
      if(r2>=0.2)
      {
        if(class(covar)=="numeric")
        {
          plot(pc~covar, 
               xlab=nm, 
               ylab=paste0("peak Heights ", pc.id), 
               main=paste(HM, cond))
        }else
        {
          boxplot(pc~covar, 
                  xlab=nm, 
                  ylab=paste0("peak Heights ", pc.id), 
                  main=paste(HM, cond))
        }
      }
      return(r2)
    })
  })
  colnames(res) = colnames(factors.est)
  return(res)
}) 
names(R2.PCvsKnownCovar) = names(factors.est.list)
dev.off()


#correlation among each covar
hm.covarNM.num = colnames(hm.knownCovar)[sapply(colnames(hm.knownCovar), FUN=function(nm) class(hm.knownCovar[[nm]])!="character")]
hm.knownCovar.adjR2 = sapply(hm.covarNM.num,
FUN=function(nm)
{
  y=hm.knownCovar[[nm]]
  sapply(colnames(hm.knownCovar),
  FUN=function(nm2)
  {
    if (nm == nm2)
    {
      return(1)
    }
     #print(paste(nm, nm2))
    x=hm.knownCovar[[nm2]]
    
    if(length(levels(factor(x[!is.na(y)])))==1)
    {
      return(NA)
    }
    if(sum((!is.na(y)) & (!is.na(x))) <=10)
    {
      return(NA)
    }
    
    mod=lm(y~x)
    adjR2 = summary(mod)$adj.r.squared
    
    if(!is.na(adjR2) && adjR2<0)
      adjR2=0
    return(adjR2)
  })
})
colnames(hm.knownCovar.adjR2) = hm.covarNM.num
rownames(hm.knownCovar.adjR2) = colnames(hm.knownCovar)
hm.knownCovar.adjR2 = hm.knownCovar.adjR2[apply(hm.knownCovar.adjR2, 1, FUN=function(x) sum(!is.na(x))>1), ]



pdf(file.ou.overall.pdf, height=8, width=12)
#layout(matrix(1:4, ncol=2))

# smoothScatter(pks.median, pks.CV,
#               main=HM)
# smoothScatter(pks.median, pks.sd,
#               main=HM)

# cor.amongPCs = cor(PCA.filt.results$topCV$x, PCA.filt.results$topMedian$x)
# rownames(cor.amongPCs)= paste0("topCV_",rownames(cor.amongPCs))
# colnames(cor.amongPCs)= paste0("topMedian_",colnames(cor.amongPCs))

# pheatmap(cor.amongPCs,
#          #color = cols,
#          cluster_rows = F,
#          cluster_cols = F)
# 

corrs = cor(t(hm.knownCovar.adjR2), use="pairwise.complete.obs")
hclu.res = hclust(as.dist(1-corrs))
hm.knownCovar.sorted = rownames(hm.knownCovar.adjR2)[hclu.res$order]

pheatmap(hm.knownCovar.adjR2[hm.knownCovar.sorted, intersect(hm.knownCovar.sorted, colnames(hm.knownCovar.adjR2))],
         #color = cols,
         cluster_rows = F,
         cluster_cols = F,
         main="adjust.R2 among known covariates of individuals")


for(nm in names(R2.PCvsKnownCovar))
{
  R2.nm = R2.PCvsKnownCovar[[nm]]
  R2.nm.filt = R2.nm[apply(R2.nm, 1, max)>0.1, ]
  
  
  knownCov.sorted= rownames(R2.nm.filt)[order(apply(R2.nm.filt, 1, FUN=function(x)which.max(abs(x))))]
  #pval.notes=matrix("", ncol=ncol(R2.PCvsKnownCovar[[nm]]), nrow=nrow(R2.PCvsKnownCovar[[nm]]))
  cellFract.filt = intersect(knownCov.sorted, CellTYPEs)
  
  
  pheatmap(R2.nm.filt[knownCov.sorted, ],
           color = cols,
           cluster_cols = F,
           cluster_row = F,
           show_colnames = T,
           show_rownames = T,
           #notes = 
           main=paste(HM, nm))
  

  #strongest R2
  layout(matrix(1:4, ncol=1))
  R2.max=apply(R2.nm.filt[knownCov.sorted, 1:10], 1, max)
  R2.max.mean = mean(R2.max[cellFract.filt])
  R2.max.sd= sd(R2.max[cellFract.filt])
  R2.max.max = max(R2.max[cellFract.filt])
  R2.max.min = min(R2.max[cellFract.filt])
  barplot(R2.max,
          las=2,
          main= paste(nm, "strongest R2 with top 10 PCs\n cell fraction", 
                      R2.max.mean, R2.max.sd,
                      R2.max.max, R2.max.min))
}

# write.table(R2.nm.filt[knownCov.sorted, ],
#             file=paste0(dir.ou, "Fig.S2g.Brain.PCvscovariates.tsv"),
#             sep="\t",
#             quote=F)
# write.table(R2.nm.filt[knownCov.sorted, ],
#             file=paste0(dir.ou, "Fig.S3a.Brain.Peervscovariates.tsv"),
#             sep="\t",
#             quote=F)

# 
dev.off()

#
f.o.tsv=paste0(dir.ou, HM, ".VarByKnownCovar.tsv")
nm="peers"
R2.nm = R2.PCvsKnownCovar[[nm]]
R2.nm.filt = R2.nm[apply(R2.nm, 1, max)>0.1, ]
knownCov.sorted= rownames(R2.nm.filt)[order(apply(R2.nm.filt, 1, FUN=function(x)which.max(abs(x))))]
write.table(R2.nm.filt[knownCov.sorted, ], file=f.o.tsv, sep="\t", quote=F)



# 
# pdf(file.ou.varProp.pdf, width = 20, height=6)
# 
# for(nm in names(PCA.filt.results))
# {
#   pca.summary = summary(PCA.filt.results[[nm]])
#   pc.varProp = pca.summary$importance["Proportion of Variance", ]
#   
#   R2.nm = R2.PCvsKnownCovar[[nm]]
#   R2.nm.filt = R2.nm[apply(R2.nm, 1, max)>0.1, ]
#   
#   
#   knownCov.varProp = R2.nm %*% cbind(pc.varProp)
#   knownCov.sorted= rownames(R2.nm.filt)[order(apply(R2.nm.filt, 1, FUN=function(x)which.max(abs(x))))]
#   #knownCov.sorted.nm = knownCov.sorted
#   #knownCov.sorted.nm[knownCov.sorted.nm %in% names(cellType2nm)]= paste0("prop. ", cellType2nm[knownCov.sorted.nm[knownCov.sorted.nm %in% names(cellType2nm)]])
#   cellFract.filt = knownCov.sorted[grepl("Cell Fraction", knownCov.sorted)]
#   cellFract.varProp.max= max(knownCov.varProp[cellFract.filt,])
#   
#   barplot(pc.varProp, main= paste0(nm, " porportion of Variance for each PC"), las=2, cex.names = 0.4)
#   barplot(knownCov.varProp[knownCov.sorted,1], 
#           names.arg=knownCov.sorted,
#           main= paste0(nm, " porportion of Variance for each covariate\nmax prop by cell fraction: ", cellFract.varProp.max), 
#           las=2, cex.names = 0.5)
#   
#   
#   
# }
# 
# dev.off()


save(#pks.median,
     #pks.CV,
     hm.knownCovar,
     hm.knownCovar.adjR2,
     factors.est.list,
     R2.PCvsKnownCovar,
     file=file.ou.RData)


