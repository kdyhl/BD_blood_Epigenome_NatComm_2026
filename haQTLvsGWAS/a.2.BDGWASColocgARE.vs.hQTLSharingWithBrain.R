#V1.1
#summary
#great


#shared LD structure since GWAS is imputed based on GTEx EUR cohort



if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))




options(scipen = 999)

#library(data.table)
library(pheatmap)
#library(rGREAT)
#library("RColorBrewer")
# library(grid)
# library(gridExtra)
gARE.cutoff="emp.p.fdr.cutoff0.2"
HM="H3K27ac"
Brain.shared.p.cutoff=0.01


file.in.coloc.blood.gARE.RData=paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/PGC.BP/", HM, "_100k/H3K27ac.bulkhQTLvsGWAS.RData")
file.in.coloc.brain.gARE.RData="a_eGTExhQTLvsGWAS_bycoloc_V1.1_allGWAS/PGC.BP/Brain_100k/Brain.bulkhQTLvsGWAS.RData"

#file.in.mAREs.hg19.RData="../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/3activeHMs.mergedPeak2HMPeak.hg19.RData"
#file.in.mARE.modules.RData="../peakVariationAcrossTiss/a_4_modules_annotations_mergedPk5HMs_Epimap_CSandeGTExActivity/mARE2ModuleGrpName.RData"
#file.in.ARE.hg38toHg19.RDS=paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".pk.hg38tohg19.RDS")

file.in.gARE=paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, ".haQTLPeak.cutoffs.hg38.RData")

file.in.haQTLSharedWithBrain.RData="../haQTL/c_compHaQTLWitheGTExBrain_V1.1_sharedgAREs/H3K27ac.vseGTEx_Brain_100k_BPhaQTLemp.p.fdr.cutoff0.2.RData"


dir.ou=paste0("a_2_GWAS.coloc.gARE_VS_haQTLSharedWithBrain/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData =paste0(dir.ou, HM, ".BD.coloc.gARE.QTLSharing.RData")
file.ou.pdf =paste0(dir.ou, HM, ".BD.coloc.gARE.QTLSharing.pdf")
#

#
#pk.hg38tohg19=readRDS(file.in.ARE.hg38toHg19.RDS)
load(file.in.coloc.blood.gARE.RData)

load(file.in.haQTLSharedWithBrain.RData)
load(file.in.gARE)

gAREs = hQTLPeaks.list[[gARE.cutoff]]

blood.gARE.hg38.H4=sapply(hQTLPeaks.QTLvsGWAS,
FUN=function(res)
{
  res$coloc$summary["PP.H4.abf"]
})
names(blood.gARE.hg38.H4)= gsub(".PP.H4.abf", "", names(blood.gARE.hg38.H4))

load(file.in.coloc.brain.gARE.RData)
brain.gARE.hg38.H4 = sapply(gwasvshaQTL.coloc,
FUN=function(x)
{
  x$coloc["PP.H4.abf"]
})
names(brain.gARE.hg38.H4) = names(gwasvshaQTL.coloc)
#


gAREs.sharedWithBrain = levels(factor(topBPhaQTL.ovlpPeaks.filt.df$BP.pk)) #all BD gARE detected in brain
gAREs.withBrain.ovlpGWAS= intersect(names(blood.gARE.hg38.H4), gAREs.sharedWithBrain) #gAREs near GWAS

gARE2bestBrain.p=sapply(split(topBPhaQTL.ovlpPeaks.filt.df, topBPhaQTL.ovlpPeaks.filt.df$BP.pk),
FUN=function(df)
{
  min(df$eGTEx.pval)
})

gARE2bestBrain.hQTL.H4=sapply(split(topBPhaQTL.ovlpPeaks.filt.df, topBPhaQTL.ovlpPeaks.filt.df$BP.pk),
FUN=function(df)
{
  brain.gAREs=levels(factor(df$eGTEx.pk))
  max(brain.gARE.hg38.H4[brain.gAREs], na.rm=T)
  
})

gARE2DirectConsist=sapply(split(topBPhaQTL.ovlpPeaks.filt.df, topBPhaQTL.ovlpPeaks.filt.df$BP.pk),
FUN=function(df)
{
  any(sign(df$BP.beta) == sign(df$eGTEx.beta))
})

#

gAREs.ovlpBrain.EnrichGWAS.p=fisher.test(matrix(c(length(gAREs.withBrain.ovlpGWAS),
                                               length(setdiff(gAREs.sharedWithBrain, gAREs.withBrain.ovlpGWAS)),
                                               length(setdiff(names(blood.gARE.hg38.H4), gAREs.withBrain.ovlpGWAS)),
                                               length(setdiff(gAREs, union(gAREs.sharedWithBrain, names(blood.gARE.hg38.H4))))),ncol=2),
                                       alternative="greater"
                                       )$p.val

#

gAREs.coloc=names(blood.gARE.hg38.H4)[blood.gARE.hg38.H4>=0.5]
gAREs.withBrain.ovlpGWAS.coloc=intersect(gAREs.withBrain.ovlpGWAS, gAREs.coloc)#gAREs.withBrain.ovlpGWAS[blood.gARE.hg38.H4[gAREs.withBrain.ovlpGWAS]>=0.5]
gAREs.ovlpBrain.EnrichColoc.p=fisher.test(matrix(c(length(gAREs.withBrain.ovlpGWAS.coloc),
                                                  length(setdiff(gAREs.withBrain.ovlpGWAS, gAREs.withBrain.ovlpGWAS.coloc)),
                                                  length(setdiff(gAREs.coloc, gAREs.withBrain.ovlpGWAS.coloc)),
                                                  length(setdiff(names(blood.gARE.hg38.H4), union(gAREs.coloc, gAREs.withBrain.ovlpGWAS)))),ncol=2),
                                         alternative="greater"
                                         )$p.val

#
#gAREs.withBrain.ovlpGWAS.directConsist = gAREs.withBrain.ovlpGWAS[gARE2DirectConsist[gAREs.withBrain.ovlpGWAS]]
gAREs.QTLSharedwithBrain = names(gARE2DirectConsist)[gARE2DirectConsist
                                                            & gARE2bestBrain.p<=Brain.shared.p.cutoff]

gAREs.QTLSharedwithBrain.ovlpGWAS = intersect(gAREs.QTLSharedwithBrain, blood.gARE.hg38.H4)
gAREs.QTLSharedAndcoloc = intersect(gAREs.coloc, gAREs.QTLSharedwithBrain )
gAREs.QTLShared.EnrichColoc.p=fisher.test(matrix(c(length(gAREs.QTLSharedAndcoloc),
                                                   length(setdiff(gAREs.coloc, gAREs.QTLSharedAndcoloc)),
                                                   length(setdiff(gAREs.QTLSharedwithBrain.ovlpGWAS, gAREs.QTLSharedAndcoloc)),
                                                   length(setdiff(blood.gARE.hg38.H4, union(gAREs.coloc, gAREs.QTLSharedwithBrain )))),ncol=2),
                                          alternative="greater"
)$p.val


pdf(file.ou.pdf, width=9, height=6)
plot(x=-log10(gARE2bestBrain.p[gAREs.withBrain.ovlpGWAS]),
     y=blood.gARE.hg38.H4[gAREs.withBrain.ovlpGWAS],
     col=c("grey", "red")[gARE2DirectConsist[gAREs.withBrain.ovlpGWAS]+1],
     pch=c(1, 4)[(gARE2bestBrain.hQTL.H4[gAREs.withBrain.ovlpGWAS]>0.2) +1],
     cex=c(.5, 2)[(gARE2bestBrain.hQTL.H4[gAREs.withBrain.ovlpGWAS]>0.2) +1],
     lwd=c(.5, 3)[(gARE2bestBrain.hQTL.H4[gAREs.withBrain.ovlpGWAS]>0.2) +1])
     #cex=0.3,
     # main=paste0("gAREs.ovlpBrain.EnrichGWAS.p:", signif(gAREs.ovlpBrain.EnrichGWAS.p, 2), "\n",
     #             "gAREs.ovlpBrain.EnrichColoc.p:", signif(gAREs.ovlpBrain.EnrichColoc.p, 2), "\n",
     #             "gAREs.directConsist.EnrichColoc.p:", signif(gAREs.directConsist.EnrichColoc.p, 2)))

legend("topright",
       legend = c("consistentDirect", "nonConsistent", "overlapped brain gARE coloc >.2"),
       col=c("red", "grey", "black"),
       pch=c(1, 1, 4))
dev.off()



save(gAREs.withBrain.ovlpGWAS,
     gARE2bestBrain.p,
     blood.gARE.hg38.H4,
     gARE2bestBrain.hQTL.H4,
     gAREs.ovlpBrain.EnrichGWAS.p,
     gAREs.ovlpBrain.EnrichColoc.p,
     gAREs.directConsist.EnrichColoc.p,
     file=file.ou.RData
     )