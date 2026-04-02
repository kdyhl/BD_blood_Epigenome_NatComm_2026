#V1.2.1
#a.geneBodyDP.meta_V1.2.1_allGenes.R
#extend analysis to all genes

#V1.2
#filter by distance
#remove those negative MR-z 
#differential peak from 

#mapping peaks to gene body based on hg38 coordinates
#do meta-analysis of differential signal of peaks taking into consideration of sign between peak and gene expression
#including peaks


#link H3K36me3 to gene body
#check eQTL colocalizaiton with H3K36me3
#function of those genes

if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

library(meta)
library(pheatmap)
library(limma)
library(enrichR)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ggplot2)
# library(DOSE)
# source("~/lhou.compbio/codes/mylib/RCode/clusterProfiler_auxiliary_function.R")

#KEGG.FDR=0.05
#
HMs= c("H3K27ac","H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")
#MR.PVAL.CUTOFF=0.05 #correction pvalue threshold
GENEBODY.UPPER=100000 #based on the distribution from both ./haQTLvseQTL/b_summary_Links.coloc.MR.PRS_V1.2_addGeneBody and ../peakMergeNormal/d_CorrAcrossIndivBtwHMs_V1.2_promVSgeneBody
GENEBODY.LOWER= 1000
#file.in.gene.hg19.bed = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.annotation.OnlyGene.gtf

files.in.diffPeak.RData = paste0("~/lhou.compbio/codes/MayoBipolar_2019_10/peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
names(files.in.diffPeak.RData) = HMs

# files.in.haQTLPeak.hg38.bed = paste("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.2_pkfilt_gINT_peer_rmSexAge_100k/", HMs, ".haQTLPeak.fdr0.05.hg38.bed", sep="")
# names(files.in.haQTLPeak.hg38.bed)= HMs

files.in.peak2geneBody.RData= paste0("../peakMergeNormal/b_bgPeak_annotByGeneBody_hg19/", HMs, ".geneBody.RData")
names(files.in.peak2geneBody.RData) = HMs

file.in.eGenes = "/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/Whole_Blood.v8.egenes.txt.gz"

file.in.DEG.RData="~/lhou.compbio/data/expression/GSE124326_BP_blood_exp/BP_blood_limmaRes.RData"

# dir.tmp="~/hptmp/MayoBipolar/peakToGeneBody/"
# dir.create(dir.tmp, showWarnings = F, recursive = T)
# file.tmp.geneBody.hg19.bed = paste0(dir.tmp, "geneBody.hg19.bed")
# file.tmp.pks.hg19.bed = paste0(dir.tmp, "pks.hg19.bed")
# file.tmp.geneBodyVSpk.hg19.bed = paste0(dir.tmp, "geneBodyVspk.ovlp.hg19.bed")

dir.ou="a_geneBody_DP_meta_V1.2.1_filtByDist_allGenes/DP.V1.4.1_disease/"
dir.create(dir.ou, showWarnings = F, recursive = T)
#file.ou.pkAndGeneBody.RDS = paste0(dir.ou, "pkAndGeneBody.RDS")
file.ou.RData=paste0(dir.ou, "genebody.DP.meta.RData")
file.ou.pdf=paste0(dir.ou, "genebody.DP.meta.allGenes.pdf")
#file.ou.KEGG.pdf=paste0(dir.ou, "geneBody2Peak.DP.meta.KEGG.pdf")

#
eGenes.info = read.table(file.in.eGenes, sep="\t", row.names = 1, header = T, stringsAsFactors = F)
esmGid.v2gname=eGenes.info[, "gene_name"]
names(esmGid.v2gname) = rownames(eGenes.info)


HM.DP.score.hg19=lapply(HMs,
FUN=function(HM)
{
  load(files.in.diffPeak.RData[HM])

  dp.signals=cbind(log2FoldChange= res$diseasedisease$logFC, 
                   pvalue=res$diseasedisease$P.Value, 
                   z = sign(res$diseasedisease$logFC)* abs(qnorm(res$diseasedisease$P.Value/2)),
                   padj = res$diseasedisease$adj.P.Val)
  rownames(dp.signals) = rownames(res$diseasedisease)

  return(dp.signals)
})
names(HM.DP.score.hg19)=HMs


# HM.linkScores = lapply(HMs,
# FUN=function(HM)
# {
#   load(files.in.peak2gene.geneticScore.RData[HM])
#   scores.filt=scores.all[scores.all$geneBody & scores.all$distance<=GENEBODY.DIST & scores.all$distance>0, ]
#   gene2scores=split(scores.filt, scores.filt$gene)
# })
# names(HM.linkScores)=HMs

#

HM.geneBody2DP.meta = list()

for(hm in HMs)
{
  print(hm)
  load(files.in.peak2geneBody.RData[hm])
  pkAndGenebody.filt=pkAndGenebody[abs(pkAndGenebody$pk2tss.dist)<GENEBODY.UPPER & abs(pkAndGenebody$pk2tss.dist)>GENEBODY.LOWER,]
  
  #g2scores.filt = g2scores[sapply(g2scores, FUN=function(x) sum(x$geneBody, na.rm = T)>=1)]    
  g2gbodyPks=split(pkAndGenebody.filt, pkAndGenebody.filt$gene)    
  
  HM.geneBody2DP.meta[[hm]]=lapply(g2gbodyPks,
  FUN=function(df)
  {
    gb.pks.hg19 = df$pk
    gb.pks.hg19.ovlp= intersect(rownames(HM.DP.score.hg19[[hm]]), gb.pks.hg19)
    if(length(gb.pks.hg19.ovlp)==0) return(NA)
    gb.pks.hg19.ovlp.DP.log2fc = HM.DP.score.hg19[[hm]][gb.pks.hg19.ovlp,"log2FoldChange"] #* gb.pks.hg38.ovlp.MR.crct
    gb.pks.hg19.ovlp.DP.pval = HM.DP.score.hg19[[hm]][gb.pks.hg19.ovlp, "pvalue"]
    gb.pks.hg19.ovlp.DP.log2fc.ste = abs(gb.pks.hg19.ovlp.DP.log2fc/qnorm(gb.pks.hg19.ovlp.DP.pval/2))
    
    m = metagen(TE=gb.pks.hg19.ovlp.DP.log2fc, 
                seTE=gb.pks.hg19.ovlp.DP.log2fc.ste, 
                studlab = gb.pks.hg19.ovlp,
                comb.fixed = F, 
                comb.random=T, 
                sm = "SMD")
    # gb.pks.hg19.ovlp.num = length(gb.pks.hg19.ovlp)
    # 
    # c(pks.num=gb.pks.hg19.ovlp.num, 
    #   effect.meta=m$TE.random,
    #   pval.meta=m$pval.random
    #   ) 
  })
  save(
      HM.DP.score.hg19,
      HM.geneBody2DP.meta,
      file=file.ou.RData)
  
}


HM.geneBody2DP.meta.summary = lapply(HM.geneBody2DP.meta,
FUN=function(meta.res)
{
  res=t(sapply(meta.res,
  FUN=function(m)
  {
    if(is.na(m)) return(c(pks.num=NA, effect.meta=NA, pval.meta=NA))
    c(pks.num=m$k, 
      effect.meta=m$TE.random,
      pval.meta=m$pval.random
      ) 
  }))
  rownames(res) = names(meta.res)
  
  return(res)
})

for(hm in HMs)
{
  zscores = abs(qnorm(HM.geneBody2DP.meta.summary[[hm]][,"pval.meta"]/2)) * sign(HM.geneBody2DP.meta.summary[[hm]][,"effect.meta"])
  HM.geneBody2DP.meta.summary[[hm]] =cbind(HM.geneBody2DP.meta.summary[[hm]],
                             zscore.meta= zscores,
                             qvalue.meta = p.adjust(HM.geneBody2DP.meta.summary[[hm]][,"pval.meta"], method = "BH"))
  
}


genes.all = levels(factor(unlist(lapply(HM.geneBody2DP.meta.summary, rownames))))
#
HM.geneBody2DP.meta.z.matrix = sapply(HMs,
FUN=function(hm)
{
  zs = rep(NA, length(genes.all))
  names(zs)= genes.all
  zs[rownames(HM.geneBody2DP.meta.summary[[hm]])]=HM.geneBody2DP.meta.summary[[hm]][,"zscore.meta"]
  return(zs)
})
colnames(HM.geneBody2DP.meta.z.matrix) =HMs


#
HM.geneBody2DP.meta.TE.matrix = sapply(HMs,
FUN=function(hm)
{
  zs = rep(NA, length(genes.all))
  names(zs)= genes.all
  zs[rownames(HM.geneBody2DP.meta.summary[[hm]])]=HM.geneBody2DP.meta.summary[[hm]][,"effect.meta"]
  return(zs)
})
colnames(HM.geneBody2DP.meta.TE.matrix) =HMs


save(
  HM.DP.score.hg19,
  HM.geneBody2DP.meta,
  HM.geneBody2DP.meta.summary,
  HM.geneBody2DP.meta.z.matrix,
  HM.geneBody2DP.meta.TE.matrix,
  file=file.ou.RData
)

#
load(file.in.DEG.RData)
#rownames(diag.info) =gsub("\\.\\d+", "", rownames(diag.info))

esmGid2esmGid.v = genes.all
names(esmGid2esmGid.v) = gsub("\\.\\d+", "", genes.all)
esmGid.v2esmGid = names(esmGid2esmGid.v)
names(esmGid.v2esmGid) = genes.all



pdf(file.ou.pdf, width=12, height=12)
Lab.palette <- colorRampPalette(c("blue4", "orange", "red"), space = "Lab")
layout(matrix(1:12, nrow=4))
for(hm in HMs)
{
  hist(HM.geneBody2DP.meta.summary[[hm]][,"pval.meta"], main=hm, xlab="meta pvalue of DP on the gene body")
  smoothScatter(HM.geneBody2DP.meta.summary[[hm]][,"pks.num"],
                -log10(HM.geneBody2DP.meta.summary[[hm]][,"pval.meta"]),
                nbin=200,
                colramp = Lab.palette,
                xlab="number of peaks on the same gene",
                ylab="-log10 (meta pvalue) for the gene")
}


layout(matrix(1:25, nrow=5))
for(hm.a in HMs)
{
  
  for(hm.b in HMs)
  {
    
    gs.ovlp = intersect(rownames(HM.geneBody2DP.meta.summary[[hm.a]]),
                        rownames(HM.geneBody2DP.meta.summary[[hm.b]]))
    cor=cor(HM.geneBody2DP.meta.summary[[hm.a]][gs.ovlp,"zscore.meta"],
            HM.geneBody2DP.meta.summary[[hm.b]][gs.ovlp,"zscore.meta"],
            use="pairwise.complete.obs"
            )
    cor.test=cor.test(HM.geneBody2DP.meta.summary[[hm.a]][gs.ovlp,"zscore.meta"],
                HM.geneBody2DP.meta.summary[[hm.b]][gs.ovlp,"zscore.meta"],
                use="pairwise.complete.obs")
    smoothScatter(HM.geneBody2DP.meta.summary[[hm.a]][gs.ovlp,"zscore.meta"],
                HM.geneBody2DP.meta.summary[[hm.b]][gs.ovlp,"zscore.meta"],
                colramp = Lab.palette,
                nbin=200,
                xlab = paste(hm.a, "meta z score"),
                ylab = paste(hm.b, "meta z score"),
                main=paste0("PCC = ", signif(cor,2), 
                            "\npval=", signif(cor.test$p.value,2),
                            "\nN=", length(gs.ovlp)))
    abline(h=0, lty=2)
    abline(v=0, lty=2)
  }
}

layout(matrix(1:25, nrow=5))

for(cond in c("diagBP", "lithium", "tobacco"))
{
  for(hm in HMs)
  {
    deg.info = topTable(fit, coef=cond, number=nrow(coefficients(fit)))
    rownames(deg.info)=gsub("\\.\\d+", "", rownames(deg.info))
    #vs expression
    esmGid.ovlp = intersect(rownames(deg.info), esmGid.v2esmGid[rownames(HM.geneBody2DP.meta.summary[[hm]])])
    cor=cor(HM.geneBody2DP.meta.summary[[hm]][esmGid2esmGid.v[esmGid.ovlp],"zscore.meta"],
            deg.info[esmGid.ovlp,"t"],
            use="pairwise.complete.obs"
              )
    cor.test=cor.test(HM.geneBody2DP.meta.summary[[hm]][esmGid2esmGid.v[esmGid.ovlp],"zscore.meta"],
                      deg.info[esmGid.ovlp,"t"],
                      use="pairwise.complete.obs")
    smoothScatter(HM.geneBody2DP.meta.summary[[hm]][esmGid2esmGid.v[esmGid.ovlp],"zscore.meta"],
                  deg.info[esmGid.ovlp,"t"],
                  colramp = Lab.palette,
                  nbin=200,
                  xlab = paste(hm, "meta z score"),
                  ylab = paste(cond, "exp t value"),
                  main=paste0("PCC = ", signif(cor,2), 
                            "\npval=", signif(cor.test$p.value,2),
                            "\nN=", length(esmGid.ovlp)))
    abline(h=0, lty=2)
    abline(v=0, lty=2)
  
  }
}

#dev.off()


#

DEG.q.cutoff=0.01
gENS.list=list()
gENS.list$bg = rownames(HM.geneBody2DP.meta.z.matrix)
for(hm in HMs)
{
  genes=rownames(HM.geneBody2DP.meta.summary[[hm]])
  degs.up = genes[HM.geneBody2DP.meta.summary[[hm]][,"effect.meta"]>0 & HM.geneBody2DP.meta.summary[[hm]][,"qvalue.meta"]<=DEG.q.cutoff]
  degs.dw = genes[HM.geneBody2DP.meta.summary[[hm]][,"effect.meta"]<0 & HM.geneBody2DP.meta.summary[[hm]][,"qvalue.meta"]<=DEG.q.cutoff]
  gENS.list[[paste0(hm, ".up")]]=degs.up[!is.na(degs.up)]
  gENS.list[[paste0(hm, ".dw")]]=degs.dw[!is.na(degs.dw)]
}

for(nm in names(gENS.list))
{
  gName=esmGid.v2gname[gENS.list[[nm]]]
  gName=gName[!is.na(gName)]
  write(gName, file=paste0(dir.ou, nm, ".q", DEG.q.cutoff, ".geneSymbol.txt"))
  write(esmGid.v2esmGid[gENS.list[[nm]]], file=paste0(dir.ou, nm, ".q", DEG.q.cutoff, ".ENSG.txt"))
}

#library(HM.geneBody2DP.meta.z.matrix)
HM.geneBody2DP.meta.z.matrix[HM.geneBody2DP.meta.z.matrix>=5]=5
HM.geneBody2DP.meta.z.matrix[HM.geneBody2DP.meta.z.matrix<=-5]=-5



pheatmap(HM.geneBody2DP.meta.z.matrix[order(HM.geneBody2DP.meta.z.matrix[,"H3K36me3"]), c("H3K36me3", "H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3")],
         cluster_rows = F,
         cluster_cols = F)


dev.off()

# 
# 
# #visualize meta results
# hm="H3K36me3"
# pdf(paste0(hm,".meta.genebygene.pdf"), height=12, width=15)
# for(g in c("ENSG00000178852.15", "ENSG00000214176.9", "ENSG00000168778.11", "ENSG00000123610.4"))
# {
#   forest(HM.geneBody2DP.meta[[hm]][[g]], main=g)
# }
# dev.off()
# # 
# # 
# for(hm in HMs)
# {
#   print(hm)
#   file.ou.pdf=paste0(dir.ou, hm, ".DP2geneticLink.pdf")
# 
#   pdf(file.ou.pdf, width=12, height=12)
# 
#   layout(matrix(1:12, ncol=3))
# 
#   for(g in sample(names(HM.linkScores[[hm]]), 200))
#   {
# 
#     pks.hg38 = HM.linkScores[[hm]][[g]]$peak
#     pks.hg19 = HM.hg38tohg19[[hm]][pks.hg38]
#     pks.hg19.ovlp = intersect(pks.hg19, rownames(HM.DP.score.hg19[[hm]]))
#     pks.hg19.isOvlp = which(pks.hg19%in% rownames(HM.DP.score.hg19[[hm]]))
#     pks.hg19.ovlp.DP.z = HM.DP.score.hg19[[hm]][pks.hg19.ovlp,"z"]
#     if(sum(HM.DP.score.hg19[[hm]][pks.hg19.ovlp,"padj"]<=0.05)>=1 &
#        sum(!is.na(HM.linkScores[[hm]][[g]]$MR.z[pks.hg19.isOvlp]))>=2)
#     {
#       print(g)
#       plot(HM.linkScores[[hm]][[g]]$distance[pks.hg19.isOvlp],
#          pks.hg19.ovlp.DP.z,# * sign(HM.linkScores[[hm]][[g]]$MR.z[pks.hg19.isOvlp]),
#          xlab="distance to TSS",
#          ylab="DP.z",
#          col= c("blue", "red")[(sign(HM.linkScores[[hm]][[g]]$MR.z[pks.hg19.isOvlp])==1)+1],
#          pch = c(1,16)[HM.linkScores[[hm]][[g]]$geneBody[pks.hg19.isOvlp]+1],
#          main=g)
#       abline(h=0, lty=2)
#       abline(v=0, lty=2)
#       legend("topright",
#              col=c("red", "blue", "black", "black"),
#              legend=c("postive MR.z", "negative MR.z", "genebody", "rest"),
#              pch=c(1 ,1, 16, 1))
#     }
#   }
# 
# 
# 
#   dev.off()
# 
# 
# }

  

