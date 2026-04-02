#V1.5.1
#c.getDiff_byLimma_V1.5.1_inflamSubgrps.R
#two groups, subgroup1&2 inflammation, subgrp3&4&5, non-inflammation

#V1.5
#get differential peaks for each subgroup vs control
#focus on immune peaks
#also binned by GWAS signals

#use sva to consider unwanted variables
#


args = commandArgs(trailingOnly=TRUE)
HM = args[1]



library(limma)
library(data.table)
#library("SmartSVA")
library(sva)
library(pheatmap)
# library(clusterProfiler)
# library(AnnotationHub)
# library(org.Hs.eg.db)
#
file.in.model.RData = paste0("./c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HM, "/", HM, ".diffPeak.RData")

file.in.subgrp.RData="../sampleManifold/b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.sampClu_V1.1.RData"


dir.ou=paste0("./c_bulkDiffPeak_limma_V1.5.1_inflamSubgrp/", HM, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.diffPeak.RData = paste(dir.ou, HM, ".diffPeak.RData", sep="")
file.ou.diffPeak.pdf = paste(dir.ou, HM, ".diffPeak.pdf", sep="")
file.ou.diffPeak.BG.bed = paste(dir.ou, HM, ".diffPeak.BG.bed", sep="")
file.ou.diffPeak.bed.pref = paste(dir.ou, HM, ".diffPeak.", sep="")




load(file.in.model.RData)
load(file.in.subgrp.RData)

#

subgroup2samps=list()
subgroup2samps$inflamSubgrps=rownames(samp2grp)[samp2grp[,"group"]=="control" | samp2grp[,"clu.corHclust5"]==1|samp2grp[,"clu.corHclust5"]==2]
subgroup2samps$nonInflamSubgrps=rownames(samp2grp)[samp2grp[,"group"]=="control" | (samp2grp[,"clu.corHclust5"]!=1 &samp2grp[,"clu.corHclust5"]!=2)]

#
lib.size = rep(1e6, nrow(design.matrix.mod2))
counts.filt=samps.counts.disease[pks.filt, ]
v = voom(counts.filt, lib.size=lib.size, design =design.matrix.mod2, plot=TRUE)


res=lapply(subgroup2samps,
FUN=function(samps)
{
  samps.ovlp=intersect(samps, rownames(design.matrix.mod2))
  v.filt=v[,samps.ovlp]

  fit=lmFit(v.filt, design.matrix.mod2[samps.ovlp,])
  fit <- eBayes(fit)
  topTable(fit, coef="diseasedisease", adjust="BH", number=nrow(v))
})

res[["all"]]=topTable(fit, coef="diseasedisease", adjust="BH", number=nrow(v))



save(subgroup2samps,
  design.matrix.mod2,
  v,
  counts.filt, 
  res,
  file=file.ou.diffPeak.RData)

#
pdf(file.ou.diffPeak.pdf, width=10, height=10)
layout(matrix(1:4, ncol=2))

pks.filt=rownames(counts.filt)
smoothScatter(res$inflamSubgrps[pks.filt, "t"], res$nonInflamSubgrps[pks.filt, "t"],
              xlab="t (inflamSubgrps vs control)",
              ylab="t (nonInflamSubgrps vs control)",
              main=HM)

smoothScatter(res$all[pks.filt, "t"], res$inflamSubgrps[pks.filt, "t"],
              xlab="t (all vs control)",
              ylab="t (inflamSubgrps vs control)",
              main=HM)


smoothScatter(res$all[pks.filt, "t"], res$nonInflamSubgrps[pks.filt, "t"],
              xlab="t (all vs control)",
              ylab="t (nonInflamSubgrps vs control)",
              main=HM)


dev.off()


#
pks.filt.bed= gsub(":|-", "\t", pks.filt)
write(pks.filt.bed, file=file.ou.diffPeak.BG.bed)

#FDR
for(sbgrp in names(res))
{
  for(cutoff in c(0.01, 0.05, 0.1, 0.2))
  {
    diffPk.up = rownames(res[[sbgrp]])[res[[sbgrp]]$adj.P.Val<=cutoff & res[[sbgrp]]$t >0]
    diffPk.up= diffPk.up[!is.na(diffPk.up)]
    write(gsub(":|-", "\t", diffPk.up), file=paste(file.ou.diffPeak.bed.pref, sbgrp, ".Q", cutoff, ".up.bed", sep=""))
    
    diffPk.dw = rownames(res[[sbgrp]])[res[[sbgrp]]$adj.P.Val<=cutoff & res[[sbgrp]]$t <0]
    diffPk.dw= diffPk.dw[!is.na(diffPk.dw)]
    write(gsub(":|-", "\t", diffPk.dw), file=paste(file.ou.diffPeak.bed.pref, sbgrp, ".Q", cutoff, ".dw.bed", sep=""))
  }
}
