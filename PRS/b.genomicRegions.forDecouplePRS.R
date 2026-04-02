#generated BD associated genomic regions related to immune component, and other background regions, use them to decouple PRS
#potential regions
#a) all DPs regions
#b) permissive DPs but filter by enriched signal in Hi-C neighborhood
#c) genetically supported regions, hQTLs
#d) GWAS supported regions, GWAS-hQTL colocalizated signal
#e) ARE groups based on eGTEx
#f) ARE groups baesd on this study


LIFTOVER = "~/lhou.compbio/software/ucscUtil/liftOver"
FILE.CHAIN.hg19to38 = "~/lhou.compbio/software/ucscUtil/hg19ToHg38.over.chain"
FILE.CHAIN.hg38to19 = "~/lhou.compbio/software/ucscUtil/hg38ToHg19.over.chain"

DP.adjP.CUTOFF=0.1
DP.PVAL.CUTOFF=0.01
BD.COLOC.CUTOFF=0.1
BD.MR.CUTOFF=0.1
hQTL.CUTOFF="emp.p.fdr.cutoff0.2"
#PEAKS.TOP.N = 10000 #c(50, 40, 30, 20, 10 ,5)
HMs= c("H3K27ac", "H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")

HiCBlock.DPEnrich.Q.CUTOFF=0.1


#
# file.in.norm.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".MayoWithRef.heightsMean.depthGCNormed.RData")
# names(file.in.norm.RData) = HMs

files.in.DP.RData = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
names(files.in.DP.RData) = HMs

file.in.HiC.block.RData="../hiC.DiffPeakandDEGandTF/a_2_cmpBulkDPAcrossHMs_hiCBlock/DP.V1.4.1_sva_noEHR/Bulk.hicBlocks.DPstat.RData"
file.in.HiCRegion2Pk.RData= "../hiC.DiffPeakandDEGandTF/a_hiCLinks_withPeakAnnot_bulk_DP.BG_limma_V1.4.1//hiCLinks.RData"#"./a_hiCLinks_withPeakAnnot_bulk/hiCLinks.RData"

files.in.hQTL.RData=paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData")
names(files.in.hQTL.RData)=HMs

files.in.coloc.RData=paste0("../haQTLvsGWAS/a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/PGC.BP/", HMs, "_100k/", HMs, ".bulkhQTLvsGWAS.RData")
names(files.in.coloc.RData)=HMs

files.in.MR.RData=paste0("../haQTLvsGWAS/a_bulkhaQTLvsGWAS_byMR_Egger_V1.1_allGWAS/PGC.BP/", HMs, "_100k/", HMs, ".bulkHaQTLvsGWAS.RData")
names(files.in.coloc.RData)=HMs


file.in.AREgrps.eGTEx.RData="../../eGTEX-H3K27ac/variationAcrossTissueandIndiv/a_3_AREModulesInEpimap_V1.6_multiTiss/4TissMergedPeaks.inEpimap.H3K27acSignal.RData"

file.in.AREgrps.BDBlood.RData="../peakVariationAcrossTiss/a_4_modules_annotations_mergedPk3actHMs_Epimap/mARE2ModuleGrpName.RData"
#file.in.mARE2AREs.RData
dir.tmp="~/hptmp/MayoBipolar/b_genomicRegions_forPRS/"
dir.create(dir.tmp, showWarnings = F, recursive = T)


dir.ou=paste0("b_genomicRegions_forPRS/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RDS=paste0(dir.ou, "genomicRegions.RDS")
#file.ou.bed.pref = paste0(dir.ou, "BD_peaks.")
#file.ou.pdf = paste0(dir.ou, "ind.MOFA.DP-pval", DP.PVAL.CUTOFF, ".varProp", MOFA.PROP.VAR.CUTOFF, ".pdf")
#file.ou.cluCmp.pdf = paste0(dir.ou, "patientHeter.cmpAcrossHMs.pdf")
#

pks.list=list()

#BD differential peaks #################################
hm2DPs.list=lapply(HMs,
FUN=function(hm)
{
  load(files.in.DP.RData[hm])
  res$diseasedisease
})
names(hm2DPs.list)=HMs

pks.list$BD_DPs.hg19=unlist(lapply(hm2DPs.list,
FUN=function(res)
{
  rownames(res)[res$adj.P.Val<=DP.adjP.CUTOFF]
}))
#BD associated peaks filter by HiC
load(file.in.HiC.block.RData)
#load(file.in.HiCRegion2DP.RData)
load(file.in.HiCRegion2Pk.RData)
hiCBlock.slct = rownames(hiCBlock.DP.summary.adjp)[apply(hiCBlock.DP.summary.adjp<=HiCBlock.DPEnrich.Q.CUTOFF,1, any)] #select differential peaks in side DP.enriched hicBlocks
hiCRegion.slct = unlist(hicBlock2hicRegs[hiCBlock.slct])
hm2DPs.filtByHiC.list = lapply(HMs,
FUN=function(hm)
{
  print(hm)

  hicReg2pks=lapply(hicReg.df[[hm]],
  FUN=function(x)
  {
    unlist(strsplit(x, split=","))
  })
  names(hicReg2pks) = rownames(hicReg.df)
  
  hic.pks.slct = levels(factor(unlist(hicReg2pks[hiCRegion.slct])))
  
  pks.ovlp =intersect(hic.pks.slct, rownames(hm2DPs.list[[hm]]))
  
  #samps.ovlp = intersect(rownames(meta.all), colnames(samp2HeightsMean$depth.GCnormTeng))
  dps = pks.ovlp[hm2DPs.list[[hm]][pks.ovlp, "P.Value"] <= DP.PVAL.CUTOFF] 
 
  return(dps)
})
pks.list$BD_DPs_filtByHiC.hg19=unlist(hm2DPs.filtByHiC.list)


#QTL and coloc ########################################
hQTLPks.list=lapply(files.in.hQTL.RData,
FUN=function(f.i.RData)
{
  load(f.i.RData)
  hQTLPeaks.list[[hQTL.CUTOFF]]
})
pks.list$blood_gAREs.hg38=unlist(hQTLPks.list)


BD.coloc.pks.list=lapply(files.in.coloc.RData,
FUN=function(f.i.RData)
{
  load(f.i.RData)
  H4s=sapply(hQTLPeaks.QTLvsGWAS,
  FUN=function(res)
  {
    res$coloc$summary["PP.H4.abf"]
  })
  names(H4s)=names(hQTLPeaks.QTLvsGWAS)

  return(names(H4s)[H4s>=BD.COLOC.CUTOFF])
})
pks.list$blood_gAREs_BDColoc.hg38=unlist(BD.coloc.pks.list)
pks.list$blood_gAREs_BDColoc.hg38=pks.list$blood_gAREs_BDColoc.hg38[!is.na(pks.list$blood_gAREs_BDColoc.hg38)]



BD.MR.pks.list=lapply(files.in.MR.RData,
FUN=function(f.i.RData)
{
  load(f.i.RData)
  
  peaks.MR=t(sapply(hQTLPeaks.QTLvsGWAS,
  FUN=function(res)
  {
      if(is.na(res[[1]])) return(c(p=NA, est=NA))

      return(c(p=res$MR@Pvalue.Est, est=res$MR@Estimate))
  }))

  rownames(peaks.MR)[!is.na(peaks.MR[,"p"]) & peaks.MR[,"p"]<=BD.MR.CUTOFF]
})
pks.list$blood_gAREs_BD.MR.hg38=unlist(BD.MR.pks.list)
#pks.list$blood_gAREs_BD.MR.hg38=pks.list$blood_gAREs_BD.MR.hg38[!is.na(pks.list$blood_gAREs_BDColoc.hg38)]



#ARE groups ########################################
new_env <- new.env()
load(file.in.AREgrps.eGTEx.RData, envir = new_env)
grp2ARE.list=lapply(split(new_env$ARE.modName.df$id, new_env$ARE.modName.df$ARE.module.groups),
FUN=function(clus)
{
  unlist(new_env$cluster2ARE[clus])
})
names(grp2ARE.list)=gsub("/| |-", "_", names(grp2ARE.list))
names(grp2ARE.list)=paste0(names(grp2ARE.list), "_eGTEx.hg19")
pks.list=c(pks.list, grp2ARE.list)
#
new_env <- new.env()
load(file.in.AREgrps.BDBlood.RData, envir = new_env)
grp2ARE.list=split(names(new_env$mARE2ModGrpName), new_env$mARE2ModGrpName)
names(grp2ARE.list)=gsub("/| |-", "_", names(grp2ARE.list))
names(grp2ARE.list)=paste0(names(grp2ARE.list), "_BDBlood.hg19")
pks.list=c(pks.list, grp2ARE.list)

saveRDS(pks.list, file=file.ou.RDS)




#output in hg38
for (nm in names(pks.list))
{
  pks=pks.list[[nm]]
  pks.bed=gsub(":|-", "\t", pks)
  file.ou.hg38.bed=paste0(dir.ou, gsub("hg19", "hg38", nm), ".bed")

  print(file.ou.hg38.bed)

  if(grepl("hg38", nm))
  {
    write(paste(pks.bed, pks, sep="\t"), file.ou.hg38.bed)
  }else
  {
    file.tmp.hg19.bed=paste0(dir.tmp, "pks.hg19.bed")
    file.tmp.unmap.bed=paste0(dir.tmp, "pks.unmap.bed")

    write(paste(pks.bed, pks, sep="\t"), file.tmp.hg19.bed)
    cmd = paste(LIFTOVER, "-bedPlus=3", file.tmp.hg19.bed, FILE.CHAIN.hg19to38, file.ou.hg38.bed, file.tmp.unmap.bed, sep=" ")
    system(cmd)
  }
}

#output in hg19
for (nm in names(pks.list))
{
  pks=pks.list[[nm]]
  pks.bed=gsub(":|-", "\t", pks)
  file.ou.hg19.bed=paste0(dir.ou, gsub("hg38", "hg19", nm), ".bed")

  print(file.ou.hg19.bed)

  if(grepl("hg19", nm))
  {
    write(paste(pks.bed, pks, sep="\t"), file.ou.hg19.bed)
  }else
  {
    file.tmp.hg38.bed=paste0(dir.tmp, "pks.hg38.bed")
    file.tmp.unmap.bed=paste0(dir.tmp, "pks.unmap.bed")

    write(paste(pks.bed, pks, sep="\t"), file.tmp.hg38.bed)
    cmd = paste(LIFTOVER, "-bedPlus=3", file.tmp.hg38.bed, FILE.CHAIN.hg38to19, file.ou.hg19.bed, file.tmp.unmap.bed, sep=" ")
    system(cmd)
  }
}
