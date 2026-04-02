#V1.2
#c.findMarkerDPforCluster_V1.2_uniq.R
#since strong overlap of markers between different clusters for each histone marks from V1.1
#remove overlapped markers witth each other
#add genetic colocalization with different GWAS

#V1.1
#find marker DP based on clusters from all samples
#compare cluster against others for cases only taking into consideration of age and sex

#find marker DP that distinguish patient cluster from the rest
args = commandArgs(trailingOnly=TRUE)
HM = args[1]


if(!grepl("R version 4.0", R.version$version.string))
{
  stop("use R-4.0")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

colors = brewer.pal(10, "Set3")

HMs= c("H3K27ac", "H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")
GWAS.NMS=c("DrinksPerWeek", "SmokingInitiation", "GWASAtlas_ukb2.diabetes", "GWASAtlas_ukb2.highBloodPressure", "PGC.BP")
HiCBlock.DPEnrich.Q.CUTOFF=0.1
DP.PVAL.CUTOFF=0.01
DP.MARKER.Q.CUTOFF = 0.05


CLU.RES = "clu.corHclust5"

file.in.norm.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".MayoWithRef.heightsMean.depthGCNormed.RData")
names(file.in.norm.RData) = HMs
files.in.DP.RData = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
names(files.in.DP.RData) = HMs

file.in.HiC.block.RData="../hiC.DiffPeakandDEGandTF/a_2_cmpBulkDP.limma_addLitium.AcrossHMs_hiCBlock/Bulk.hicBlocks.DPstat.RData"
file.in.HiCRegion2Pk.RData="../hiC.DiffPeakandDEGandTF/a_hiCLinks_withPeakAnnot_bulk/hiCLinks.RData"


file.in.clu.RData=paste0("./b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.sampClu_V1.1.RData")
file.in.meta="~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"

file.in.hm.gwas.coloc.RData="../haQTLvsGWAS/a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS_summary/bulkhQTLvsGWAS.coloc.RData"

file.in.ARE.hg38toHg19.RDS=paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".pk.hg38tohg19.RDS")


dir.ou=paste0("c_markerDPForCluster_V1.2_uniq/HiC.DPs/", CLU.RES, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData=paste0(dir.ou, HM, ".markerDP.", CLU.RES, ".RData")
file.ou.pdf=paste0(dir.ou, HM, ".markerDP.", CLU.RES, ".pdf")
file.ou.markervsBD_DP.pdf=paste0(dir.ou, HM, ".markerDP.", CLU.RES, ".vsBD_DP.pdf")

#
load(file.in.HiC.block.RData)
#load(file.in.HiCRegion2DP.RData)
load(file.in.HiCRegion2Pk.RData)

hiCBlock.slct = rownames(hiCBlock.DP.summary.adjp)[apply(hiCBlock.DP.summary.adjp<=HiCBlock.DPEnrich.Q.CUTOFF,1, any)] #select differential peaks in side DP.enriched hicBlocks
hiCRegion.slct = unlist(hicBlock2hicRegs[hiCBlock.slct])


#
load(file.in.clu.RData)
# clu2samp=names(samp2clu)
# names(clu2samp)=paste0("clu_", samp2clu)
# buf =read.table(file.in.meta, sep=",", row.names = NULL, header = T, stringsAsFactors = F)
# samp2group=buf$group
# names(samp2group)=paste0("id_", rownames(buf))

samp.case=rownames(samp2grp)[samp2grp$group=="case"]
  
clus.all = levels(factor(samp2grp[samp.case, CLU.RES]))

clu2Rest.samp =lapply(clus.all,
FUN=function(clu.i)
{
  res=data.frame(clu.grp=rep("control.other", nrow(samp2grp)), stringsAsFactors = F)
  rownames(res)=rownames(samp2grp)
  res$clu.grp[samp2grp$group=="case" & samp2grp[[CLU.RES]]==clu.i] = "case.cluster"
  res$clu.grp[samp2grp$group=="control" & samp2grp[[CLU.RES]]==clu.i] = "control.cluster"
  res$clu.grp[samp2grp$group=="case" & samp2grp[[CLU.RES]]!=clu.i] = "case.other"
  res$clu.grp =factor(res$clu.grp, levels=c("case.cluster", "case.other", "control.cluster", "control.other"))
  return(res)
})
names(clu2Rest.samp) = paste0("clu_", clus.all)                         


# clu2Rest.samp = lapply(clu.all,
# FUN=function(clu)
# {
#   res=samp2clu
#   res[samp2clu!=clu]="rest"
#   res[samp2clu==clu]="clu"
#   return(res)
# })
# names(clu2Rest.samp) = paste0("clu_", clu.all)

##marker
print(HM)
load(file.in.norm.RData[HM])
load(files.in.DP.RData[HM])

samps.ovlp = intersect(rownames(samp2grp), colnames(samp2HeightsMean$depth.GCnormTeng))

  
hicReg2pks=lapply(hicReg.df[[HM]],
FUN=function(x)
{
  unlist(strsplit(x, split=","))
})
names(hicReg2pks) = rownames(hicReg.df)

hic.pks.slct = levels(factor(unlist(hicReg2pks[hiCRegion.slct])))
DPs = rownames(res$diseasedisease)[res$diseasedisease[, "P.Value"] <= DP.PVAL.CUTOFF]
hic.DPs.slct =intersect(hic.pks.slct, DPs) #intersect(rownames(res$disease), rownames(samp2HeightsMean$depth.GCnormTeng))

HM.markDP.hg19.res=lapply(clu2Rest.samp,
FUN=function(clu2Rest)
{
  
  cases.ovlp= intersect(samp.case, samps.ovlp)
  
  
  print(length(hic.DPs.slct))
  t(sapply(hic.DPs.slct,
  FUN=function(pk)
  {
    df =data.frame(pk=samp2HeightsMean$depth.GCnormTeng[pk, cases.ovlp],
                   sex =design.matrix.mod2[cases.ovlp, "GenderM"],
                   age =design.matrix.mod2[cases.ovlp, "age_atsamp"],
                   clu=clu2Rest[cases.ovlp, "clu.grp"])
    
    mod.summary=summary(lm(pk~clu+sex+age, data=df))
    mod.t = mod.summary$coef[2, "t value"]
    mod.p =mod.summary$coef[2, "Pr(>|t|)"]
    
    return(c(t=mod.t, 
             p=mod.p,
             p.adj=p.adjust(mod.p, method = "BH")))
  }))
})

HM.markDP.hg19.list=lapply(HM.markDP.hg19.res,
FUN=function(dps)
{
  rownames(dps)[dps[,"p.adj"] <= DP.MARKER.Q.CUTOFF]
})
HM.markDP.uniq.hg19.list=lapply(names(HM.markDP.hg19.list),
FUN=function(clu)
{
  setdiff(HM.markDP.hg19.list[[clu]], unlist(HM.markDP.hg19.list[setdiff(names(HM.markDP.hg19.list), clu)]))
})
names(HM.markDP.uniq.hg19.list) = names(HM.markDP.hg19.list)

save(clu2Rest.samp,
     samp2grp,
     HM.markDP.hg19.res,
     HM.markDP.uniq.hg19.list,
     file=file.ou.RData)

#

#
load(file.in.hm.gwas.coloc.RData)
pk.hg38tohg19=readRDS(file.in.ARE.hg38toHg19.RDS) 
col_fun.coloc = colorRamp2(c(0, 0.1, 0.5), c("grey", "white", "red"))


pdf(file.ou.pdf, width=12, height=12)
layout(matrix(1:20, ncol=4, byrow = T))
for(clu in names(HM.markDP.hg19.list))
{
  hist(HM.markDP.hg19.res[[clu]][,"p"], 
       main=paste(HM, clu),
       xlab="pval")
}

colors.corHclust =colors[1:length(levels(factor(samp2grp[[CLU.RES]])))]
names(colors.corHclust) = as.character(1:length(colors.corHclust))

for(clu in names(HM.markDP.uniq.hg19.list))
{
  
  samps.ovlp= intersect(rownames(clu2Rest.samp[[clu]]), colnames(samp2HeightsMean$depth.GCnormTeng))
    
  f.o.pk= paste0(dir.ou, HM, ".", clu, ".markDP.uniq.hg19.bed")
  dp.clu.marker = HM.markDP.uniq.hg19.list[[clu]]

  if(length(dp.clu.marker)!=0)
  {
    dp.clu.marker.GWASColoc=matrix(0, nrow=length(dp.clu.marker), ncol=length(hm2gwas2coloc.hg38.list[[HM]]))
    rownames(dp.clu.marker.GWASColoc)=dp.clu.marker
    colnames(dp.clu.marker.GWASColoc)=GWAS.NMS
    for(gwas in GWAS.NMS)
    {
      pk2coloc.hg19=hm2gwas2coloc.hg38.list[[HM]][[gwas]]
      names(pk2coloc.hg19) = pk.hg38tohg19[names(pk2coloc.hg19)]
  
      pks.ovlp=intersect(dp.clu.marker, names(pk2coloc.hg19))
      dp.clu.marker.GWASColoc[pks.ovlp, gwas] = pk2coloc.hg19[pks.ovlp]
    }
    
    write.table(cbind(gsub(":|-", "\t", dp.clu.marker), dp.clu.marker), file=f.o.pk, sep="\t", quote=F, col.names = F, row.names = F)
    
    pk.norm.signal = t(apply(samp2HeightsMean$depth.GCnormTeng[dp.clu.marker, samps.ovlp], 1, FUN=function(x){(x-mean(x))/sd(x)}))
    
    ha.top=HeatmapAnnotation(clu=samp2grp[samps.ovlp, CLU.RES],
                         col=list(clu=colors.corHclust))
    ha.left=rowAnnotation(coloc=dp.clu.marker.GWASColoc,
                          col=list(coloc=col_fun.coloc))
    
    ht= Heatmap(pk.norm.signal,
                column_split= clu2Rest.samp[[clu]][samps.ovlp, "clu.grp"],
                cluster_column_slices=F,
                top_annotation = ha.top,
                left_annotation = ha.left,
                row_title = paste(HM, clu),
                show_row_names=F)
    draw(ht)
  }
  
}


dev.off()
  
  
#
pdf(file.ou.markervsBD_DP.pdf, width=12, height=9)
layout(matrix(1:6, nrow=2))
for(grp in names(HM.markDP.hg19.res))
{
  pks= rownames(HM.markDP.hg19.res[[grp]])

  smoothScatter(res$diseasedisease[pks, "t"],
                HM.markDP.hg19.res[[grp]][pks, "t"],
                xlab="BD differential t",
                ylab="markerdifferential t",
                main=grp)
  abline(h=0, col="red", lty=2)
  abline(v=0, col="red", lty=2)
}
dev.off()

