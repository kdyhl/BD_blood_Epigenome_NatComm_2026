#V1.4
#d.annotateDiffPeak_byrGREAT_V1.4.1_sva_NoEHR.R
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
library("RColorBrewer")

HMs = c("H3K27ac", "H3K4me1", "H3K36me3", "H3K4me3", "H3K27me3")
CONDs.sorted= c("H3K27me3.up", "H3K27ac.dw", "H3K4me1.dw", "H3K36me3.dw", "H3K36me3.dw", 
        "H3K27ac.up", "H3K36me3.up", "H3K4me3.up", "H3K27me3.dw", "H3K4me1.up")
RANGE=10000
FDR = 0.05
FC.CUTOFF=1.5
TERMs.TOPNum=15
DP.TYPE="diseasedisease.Q0.05"
#DP.TYPE="top1000"


#
dir.in = "c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/" #"c_bulkDiffPeak_deseq2_V1.2_addCov/"
dir.ou = "d_bulkDP_limma_V1.4.1_annotByrGREAT/" #"d_bulkDiffPeakV1.2_annotByrGREAT/"
dir.create(dir.ou, showWarnings = F, recursive = T)



files.in.allPeaks = lapply(HMs,
FUN=function(hm)
{
  fs = c()
  fs["all.BG"]=paste(dir.in, hm, "/", hm, ".diffPeak.BG.bed", sep="")  #paste(ref, REF2CELLTYPE[ref], "BG", sep="\n") 
  for(type in c("up", "dw"))
  {
    f.i = paste(dir.in, hm, "/", hm, ".diffPeak.", DP.TYPE, ".", type, ".bed", sep="")  
    if(file.exists(f.i))
    {
      buf = unlist(strsplit(system(paste("wc -l ", f.i, sep=""), intern = T), split=" "))
      print(paste(f.i, buf[1]))
      if(as.numeric(buf[1])>10)
      {
        fs[type] = f.i
      }
    }
  }
  
  return(fs)
})
names(files.in.allPeaks) = HMs

files.in.allPeaks = files.in.allPeaks[sapply(files.in.allPeaks, length)!=1]


file.ou.pdf = paste(dir.ou, "diffPeak.", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.pdf", sep="")
file.ou.RData = paste(dir.ou, "diffPeak.", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.RData", sep="")


#

enrichRes.HM=lapply(names(files.in.allPeaks),
FUN=function(hm)
{
  print(hm)
  fs = files.in.allPeaks[[hm]]
  f.bg = fs["all.BG"]
  bg.bed = read.table(f.bg, sep="\t", row.names = NULL, header=F)[,1:3]
  colnames(bg.bed) = c("chr", "start", "end")
  
  res=list()
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
      r[r$Hyper_Adjp_BH<=FDR, c("name", "Hyper_Fold_Enrichment", "Hyper_Region_Set_Coverage", "Hyper_Raw_PValue",  "Total_Genes_Annotated")]
    })
    res[[nm]]=tb.filt
  }
  return(res)
})
names(enrichRes.HM) = names(files.in.allPeaks)

save(enrichRes.HM, file=file.ou.RData)


enrichRes.HM.filt=enrichRes.HM
for(HM in names(enrichRes.HM.filt))
{
  for(nm in names(enrichRes.HM.filt[[HM]]))
  {
    for(ot in names(enrichRes.HM.filt[[HM]][[nm]]))
    {
      x=enrichRes.HM.filt[[HM]][[nm]][[ot]]
      
      if(nrow(x)>1)
        enrichRes.HM.filt[[HM]][[nm]][[ot]]=x[x[,"Total_Genes_Annotated"]<=500 & 
                                              x[,"Total_Genes_Annotated"]>=10, ] # 
        
    }
  }
}

enrichRes.eachOT = list()
for(HM in names(enrichRes.HM.filt))
{
  for(nm in names(enrichRes.HM.filt[[HM]]))
  {
    for(ot in names(enrichRes.HM.filt[[HM]][[nm]]))
    {
      if(nrow(enrichRes.HM.filt[[HM]][[nm]][[ot]] )>1)
        enrichRes.eachOT[[ot]][[paste(HM, nm, sep=".")]] = enrichRes.HM.filt[[HM]][[nm]][[ot]]  
    }
  }
}

enrichRes.HMs.merge=lapply(enrichRes.eachOT,
FUN=function(res.ots)
{
  ot.names=levels(factor(unlist(lapply(res.ots,
  FUN=function(r)
  {
    r$name
  }))))
  
  enrichMatrix = matrix(0, ncol=length(res.ots), nrow=length(ot.names))
  rownames(enrichMatrix) = ot.names
  colnames(enrichMatrix)=names(res.ots)
  for(nm in names(res.ots))
  {
    ot.nms = res.ots[[nm]][,"name"]
    
    enrichMatrix[ot.nms, nm] = res.ots[[nm]][,"Hyper_Fold_Enrichment"]
  }
  
  return(enrichMatrix)
})


enrichRes.HMs.merge.filter=lapply(enrichRes.eachOT,
FUN=function(res.ots)
{
  ot.names=levels(factor(unlist(lapply(res.ots,
  FUN=function(r)
  {
    r$name
  }))))
  
  enrichMatrix = matrix(0, ncol=length(res.ots), nrow=length(ot.names))
  rownames(enrichMatrix) = ot.names
  colnames(enrichMatrix)=names(res.ots)
  for(nm in names(res.ots))
  {
    #ot.nms = res.ots[[nm]][,"name"]

    rownames(res.ots[[nm]]) = res.ots[[nm]][,"name"]
    ot.nms = rownames(res.ots[[nm]])[res.ots[[nm]][,"Hyper_Fold_Enrichment"]>=FC.CUTOFF]
    if(length(ot.nms)>TERMs.TOPNum)
    {
      ot.nms = ot.nms[order(res.ots[[nm]][ot.nms,"Hyper_Raw_PValue"], decreasing = F)][1:TERMs.TOPNum]
    }
    enrichMatrix[ot.nms, nm] = res.ots[[nm]][ot.nms,"Hyper_Fold_Enrichment"]
  }
  enrichMatrix.filt= enrichMatrix[apply(enrichMatrix, 1, FUN=function(x) any(x!=0)), ]
  return(enrichMatrix.filt)
})

pdf(file.ou.pdf, width=12, height=20)

for(ot.nm in names(enrichRes.HMs.merge.filter))
{
  print(paste(ot.nm))
  matrix.log2enrich = log2(enrichRes.HMs.merge.filter[[ot.nm]]+1)
  #matrix.log2enrich.filt = matrix.log2enrich[apply(matrix.log2enrich>=log2(2+1), 1, any), ]
  matrix.no0 = matrix(0, 
                      nrow=nrow(matrix.log2enrich),
                      ncol=ncol(matrix.log2enrich))
  matrix.no0[matrix.log2enrich!=0]=1
  if(nrow(matrix.no0)<=2 || ncol(matrix.no0)<=2 ) next
    # corrplot(matrix.log2enrich, 
    #    #p.mat=10^(-ldsc.log10ps)[row.hclust$order, col.hclust$order],
    #    is.corr=F, #insig = "blank", sig.level = .05,
    #    method = "square",
    #    #, order = "hclust"
    #    title = paste(ot.nm)
    #    )
  else
  {
    # row.hclust.order=hclust(dist(matrix.no0, method = "euclidean"))$order
    # col.hclust.order=hclust(dist(t(matrix.no0), method = "euclidean"))$order 
    # corrplot(matrix.log10enrich[row.hclust.order, col.hclust.order], 
    #    #p.mat=10^(-ldsc.log10ps)[row.hclust$order, col.hclust$order],
    #    is.corr=F, #insig = "blank", sig.level = .05,
    #    method = "square",
    #    #, order = "hclust"
    #    title = paste(ot.nm)
    #    )
    # pheatmap(matrix.log10enrich,
    #          fontsize_row=2,
    #          main=ot.nm)
    conds.ovlp=intersect(CONDs.sorted, colnames(matrix.log2enrich))
    pheatmap(matrix.log2enrich[, conds.ovlp],
            col=colorRampPalette(c("white", "purple"))(100),
            cluster_cols=F,
            fontsize_row=8,
            main=ot.nm)
  }
}


dev.off()


save(enrichRes.HM,
     enrichRes.HMs.merge.filter,
     file=file.ou.RData)




#data of heatmap for visualization
ot.nm="GO Biological Process"
file.ou.tsv=paste0(dir.ou, "diffPeak.", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.", ot.nm, ".tsv")

matrix.log2enrich = log2(enrichRes.HMs.merge.filter[[ot.nm]]+1)
conds.ovlp=intersect(CONDs.sorted, colnames(matrix.log2enrich))
# matrix.log10enrich.filt = matrix.log10enrich[apply(matrix.log10enrich>=log10(2+1), 1, any), ]
# matrix.no0 = matrix(0, 
#                     nrow=nrow(matrix.log10enrich),
#                     ncol=ncol(matrix.log10enrich))
# matrix.no0[matrix.log10enrich!=0]=1

write.table(matrix.log2enrich[, conds.ovlp], file=file.ou.tsv, sep="\t", quote=F)



