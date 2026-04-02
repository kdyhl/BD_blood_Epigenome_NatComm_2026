#V1.4
#separate cluster visualization based on up or down for each disease

#V1.3
#add background peaks when comparing fractions
#and cluster annotation together while visualization

#V1.2 add neutrophil

#v1.1 for nnPoisR results
#annotate by ChIPseeker



library(rGREAT)
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# library(clusterProfiler)
#library(corrplot)
library(pheatmap)

HMs = c("H3K27ac", "H3K4me1", "H3K36me3", "H3K4me3", "H3K27me3")
RANGE=10000
FDR = 0.05


CLU.RES = "clu.corHclust5"

dir.in.dp = "../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/"
#dir.in= paste0("c_markerDPForCluster_V1.2_uniq/V1.2.2_allSampGrp/", CLU.RES, "/") 
dir.in= paste0("c_markerDPForCluster_V1.2_uniq/HiC.DPs/", CLU.RES, "/") 
files.in.markerDP.bed = dir(dir.in, ".hg19.bed", full.names = T)

#dir.in=paste0("c_markerDPForCluster_V1.3_vsMatchedControl/HiC.DPs/", CLU.RES, "/")
#files.in.markerDP.bed = dir(dir.in, paste0(CLU.RES, ".hg19.bed"), full.names = T)

#"d_bulkDiffPeakV1.2_annotByrGREAT/"
#dir.create(dir.ou, showWarnings = F, recursive = T)

files.in.allPeaks = lapply(HMs,
FUN=function(hm)
{
  fs = c()
  fs["all.BG"]=paste(dir.in.dp, hm, "/", hm, ".diffPeak.BG.bed", sep="")  #paste(ref, REF2CELLTYPE[ref], "BG", sep="\n") 
  
  fs.hm = files.in.markerDP.bed[grepl(hm, files.in.markerDP.bed)]
  names(fs.hm)= sapply(fs.hm,
  FUN=function(f)
  {
    gsub(".*(clu_\\d+).*", "\\1",basename(f))
  })
  fs=c(fs, fs.hm)
  return(fs)
})
names(files.in.allPeaks) = HMs

files.in.allPeaks = files.in.allPeaks[sapply(files.in.allPeaks, length)!=1]


file.ou.RData = paste0(dir.in, "HMs.subtype.markerDP.fdr", FDR, ".", CLU.RES, ".annotByrGREAT.RData")
file.ou.pdf = paste0(dir.in, "HMs.subtype.markerDP.fdr", FDR, ".", CLU.RES, ".annotByrGREAT.pdf")


#

enrichRes.HM=list()

for(hm in names(files.in.allPeaks))
{
  print(hm)
  fs = files.in.allPeaks[[hm]]
  f.bg = fs["all.BG"]
  bg.bed = read.table(f.bg, sep="\t", row.names = NULL, header=F)[,1:3]
  colnames(bg.bed) = c("chr", "start", "end")
  
  enrichRes.HM[[hm]]=list()
  for(nm in names(fs[-1]))
  {
    print(paste0(" ", nm))
    dp.bed =  read.table(fs[nm], sep="\t", row.names = NULL, header=F)
    colnames(dp.bed) = c("chr", "start", "end")
    
    job = submitGreatJob(dp.bed, bg.bed, species = "hg19")
    tb = getEnrichmentTables(job, category=c("GO", "Phenotype", "Genes"))
    tb.filt = lapply(tb,
    FUN=function(r)
    {
      r[r$Hyper_Adjp_BH<=FDR, c("name", "Hyper_Fold_Enrichment", "Hyper_Region_Set_Coverage", "Hyper_Raw_PValue", "Total_Genes_Annotated")]
    })
    enrichRes.HM[[hm]][[nm]]=tb.filt
    save(enrichRes.HM, file=file.ou.RData)
  }
}





enrichRes.eachOT = list()
for(HM in names(enrichRes.HM))
{
  for(nm in names(enrichRes.HM[[HM]]))
  {
    for(ot in names(enrichRes.HM[[HM]][[nm]]))
    {
      if(nrow(enrichRes.HM[[HM]][[nm]][[ot]] )>1)
        enrichRes.eachOT[[ot]][[paste(HM, nm, sep=".")]] = enrichRes.HM[[HM]][[nm]][[ot]]  
    }
  }
}

enrichRes.HMs.merge=lapply(enrichRes.eachOT,
FUN=function(res.ots)
{
  ot.names=levels(factor(unlist(lapply(res.ots,
  FUN=function(r)
  {
    r$name[order(r$Hyper_Raw_PValue)[1:min(5, nrow(r))]]
  }))))
  
  enrichMatrix = matrix(0, ncol=length(res.ots), nrow=length(ot.names))
  rownames(enrichMatrix) = ot.names
  colnames(enrichMatrix)=names(res.ots)
  for(nm in names(res.ots))
  {
    ot.nms = res.ots[[nm]][,"name"]
    ot.nms.ovlp = intersect(ot.nms, ot.names)
    enrichMatrix[ot.nms.ovlp, nm] = res.ots[[nm]][res.ots[[nm]][,"name"] %in% ot.nms.ovlp,"Hyper_Fold_Enrichment"]
  }
  
  return(enrichMatrix)
})
  

pdf(file.ou.pdf, width=12, height=12)

for(nm in names(enrichRes.HMs.merge))
{
  pheatmap(log10(enrichRes.HMs.merge[[nm]]+1),
           main=nm)
}



dev.off()


#
f.o.heatmap.tsv=paste0(dir.in, "HMs.subtype.markerDP.fdr", FDR, ".", CLU.RES, ".annotByrGREAT.heatmap.tsv")
nm="GO Biological Process"
write.table(log10(enrichRes.HMs.merge[[nm]]+1), f.o.heatmap.tsv, sep="\t", quote=F)


# for(ot.nm in names(enrichRes.HMs.merge))
# {
#   print(paste(ot.nm))
#   matrix.log10enrich = log10(enrichRes.HMs.merge[[ot.nm]]+1)
#   matrix.log10enrich.filt = matrix.log10enrich[apply(matrix.log10enrich>=log10(2+1), 1, any), ]
#   matrix.no0 = matrix(0, 
#                       nrow=nrow(matrix.log10enrich),
#                       ncol=ncol(matrix.log10enrich))
#   matrix.no0[matrix.log10enrich!=0]=1
#   if(nrow(matrix.no0)<=2 || ncol(matrix.no0)<=2 )
#     corrplot(matrix.log10enrich, 
#        #p.mat=10^(-ldsc.log10ps)[row.hclust$order, col.hclust$order],
#        is.corr=F, #insig = "blank", sig.level = .05,
#        method = "square",
#        #, order = "hclust"
#        title = paste(ot.nm)
#        )
#   else
#   {
#     # row.hclust.order=hclust(dist(matrix.no0, method = "euclidean"))$order
#     # col.hclust.order=hclust(dist(t(matrix.no0), method = "euclidean"))$order 
#     # corrplot(matrix.log10enrich[row.hclust.order, col.hclust.order], 
#     #    #p.mat=10^(-ldsc.log10ps)[row.hclust$order, col.hclust$order],
#     #    is.corr=F, #insig = "blank", sig.level = .05,
#     #    method = "square",
#     #    #, order = "hclust"
#     #    title = paste(ot.nm)
#     #    )
#     pheatmap(matrix.log10enrich,
#              fontsize_row=3,
#              main=ot.nm)
#   }
  
#}



save(enrichRes.HM,
     enrichRes.HMs.merge,
     file=file.ou.RData)






