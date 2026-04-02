#V1.4
#d.annotateDiffPeak_byrGREAT_V1.4_sva.R
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
library(tidyr)
library(ggplot2)

colfunc <- colorRampPalette(c("white", "purple"))
cols= colfunc(10)


HMs = c("H3K27ac", "H3K4me1", "H3K36me3", "H3K4me3", "H3K27me3")
DP.NUM.CUTOFF=50
RANGE=10000
FDR = 0.05
FE.CUTOFF=1.5
DP.TYPE="_current.Q0.2"
#DP.TYPE="top1000"


#
dir.in = "/" #"c_bulkDiffPeak_deseq2_V1.2_addCov/"
dir.ou = "d_bulkDP_limma_V1.4_annotByrGREAT/" #"d_bulkDiffPeakV1.2_annotByrGREAT/"
dir.create(dir.ou, showWarnings = F, recursive = T)



files.in.allPeaks = lapply(HMs,
FUN=function(hm)
{
  fs = c()
  fs["all.BG"]=paste(dir.in, hm, "/", hm, ".diffPeak.BG.bed", sep="")  #paste(ref, REF2CELLTYPE[ref], "BG", sep="\n") 
  
  fs.drug.dp =system(paste0("find ", dir.in, hm, "/*", DP.TYPE, "*bed"), intern = T)
  names(fs.drug.dp)= sapply(fs.drug.dp,
  FUN=function(f.i)
  {
    nm=sub(paste0(hm,".diffPeak."), "", basename(f.i))
    nm=sub(".bed", "", nm)
  })
  fs
  for(nm in names(fs.drug.dp))
  {
    f.i = fs.drug.dp[nm]
    
    buf = unlist(strsplit(system(paste("wc -l ", f.i, sep=""), intern = T), split=" "))
    print(paste(f.i, buf[1]))
    if(as.numeric(buf[1])>=DP.NUM.CUTOFF)
    {
      fs[nm] = f.i
    }
    
  }
  
  return(fs)
})
names(files.in.allPeaks) = HMs

files.in.allPeaks = files.in.allPeaks[sapply(files.in.allPeaks, length)!=1]


file.ou.pdf = paste(dir.ou, "diffPeak.drug", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.pdf", sep="")
#file.ou.plot.pdf = paste(dir.ou, "diffPeak.drug", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.plot.pdf", sep="")
file.ou.RData = paste(dir.ou, "diffPeak.drug", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.RData", sep="")
file.ou.tsv.pref=paste0(dir.ou, "diffPeak.drug", DP.TYPE, ".fdr", FDR, ".annot.byrGREAT.") 

#total mARE
#file.in.5HM.mARE.RData = "../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2HMPeak.hg19.RData"
#files.in = system(paste0("find c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/*/*diseasedisease.Q0.05*bed"), intern=T)
#
#load(file.in.5HM.mARE.RData)
# mAREs.dp= lapply(files.in, 
# FUN=function(f.i)
# {
#   hm=unlist(strsplit(f.i, split="/", fixed = T))[2]
#   buf=read.table(f.i, sep="\t", header=F, row.names=NULL, stringsAsFactors=F)
#   dps=paste0(buf[,1], ":", buf[,2], "-", buf[,3])
#   hm.peak2mergedPeak[[hm]][dps]
# })
# print(length(levels(factor(unlist(mAREs.dp))))
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
      print(sum(r$Hyper_Adjp_BH<=0.1))
      r[r$Hyper_Adjp_BH<=FDR, c("name", "Hyper_Fold_Enrichment", "Hyper_Region_Set_Coverage", "Hyper_Raw_PValue", "Total_Genes_Annotated")]
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

  

pdf(file.ou.pdf, width=15, height=20)

par(mar=c(3,30,4,5))

types.all = levels(factor(unlist(lapply(enrichRes.HM.filt, names))))
# for(ot in c("GO Molecular Function", "GO Biological Process", "GO Cellular Component", "Mouse Phenotype", 
#             "Mouse Phenotype Single KO", "Human Phenotype", "Ensembl Genes"))
#   
for(ot in c("GO Biological Process"))
{
  layout(matrix(1:10, ncol=2))
  
  drug.hm.term.p.df=data.frame()
  all.terms.slct=c()
  for (type in types.all)
  {
    drug=sub("(.*?)_.*", "\\1", type)
    for(hm in names(enrichRes.HM.filt))
    {
      if(!(type %in% names(enrichRes.HM.filt[[hm]])) || !(ot %in% names(enrichRes.HM.filt[[hm]][[type]]))|| nrow(enrichRes.HM.filt[[hm]][[type]][[ot]])==0 ) next
        
      
      res = enrichRes.HM.filt[[hm]][[type]][[ot]]
      rownames(res) = res$name
      res.filt=res[res$Hyper_Fold_Enrichment>=FE.CUTOFF,]
      #res= res[res$Hyper_Fold_Enrichment>1.2, ]
      if(nrow(res)==0) next
      term2log10p=-log10(res.filt$Hyper_Raw_PValue)
      names(term2log10p) = res.filt$name
      
      term2log10p.slct=rev(sort(term2log10p, decreasing = T)[1:min(length(term2log10p), 20)])
      barplot(term2log10p.slct,
              horiz=T,
              main=paste(hm, type, ot),
              col="red",
              las=2,
              xlab="-log10P")
      
      #terms.slct= res$name[res$Hyper_Fold_Enrichment>=FE.CUTOFF]
      #print(paste(paste0(hm, ".", type), names(term2log10p.slct)))
      all.terms.slct=c(all.terms.slct, names(term2log10p.slct))
      drug.hm.term.p.df=rbind(drug.hm.term.p.df,
                        data.frame(drug=drug,
                                   hm=hm,
                                   term=res$name,
                                   p=res$Hyper_Raw_PValue,
                                   foldEnrich=res$Hyper_Fold_Enrichment,
                                   stringsAsFactors = F))
    }
  } 
  
  all.terms.slct=levels(factor(all.terms.slct))
  terms.slct.list=list()
  terms.slct.list$immune = all.terms.slct[grepl("lymphocyte|cytokine|immun|leukocyte|myeloid|chemokine|neutrophil|T cell|B cell|cytotoxicity|inflammatory|antigen|toll-like receptor|interferon|macrophage|interleukin|granulocyte|T follicular helper cell", all.terms.slct)]
  terms.slct.list$rest = setdiff(all.terms.slct,terms.slct.list$immune)
    
    
  drug2term.p= lapply(split(drug.hm.term.p.df, drug.hm.term.p.df$drug),
  FUN=function(drug.term.ps.df)
  {
    term2minP=sapply(split(drug.term.ps.df, drug.term.ps.df$term), 
    FUN=function(df)
    {
      min(df$p)
    })
    data.frame(drug=drug.term.ps.df$drug[1],
               term=names(term2minP),
               p=term2minP,
               stringsAsFactors = F
               )
  })
  drug2term.p.df =do.call(rbind, drug2term.p)
    
  buf=spread(drug2term.p.df, term, p, fill=1)
  drug.term.p.matrix=as.matrix(buf[, -1])
  rownames(drug.term.p.matrix)=buf[,1]
  
  for(nm in names(terms.slct.list))
  {
    ps = -log10(t(drug.term.p.matrix[,terms.slct.list[[nm]]]))
    ps.sorted = ps[, order(apply(ps, 2, max), decreasing = T)] 
    
  
    #
    p=pheatmap(ps.sorted,
           color = cols,
           cluster_row=T,
           cluster_col=F,
           main=paste(ot, nm))
    
    write.table(ps.sorted[p$tree_row$order,],
                file=paste0(file.ou.tsv.pref, nm, ".tsv"),
                sep="\t",
                quote=F)

    #linegraph
    ps.rowColSorted=ps.sorted[rev(p$tree_row$order),]

    p.list=lapply(1:ncol(ps.rowColSorted), 
    FUN=function(i)
    {
      data.frame(term=rownames(ps.rowColSorted),
                logP=ps.rowColSorted[,i],
                drug=colnames(ps.rowColSorted)[i],
                stringsAsFactors=F)   
    })
    p.df=do.call(rbind, p.list)
    p.df$term=factor(p.df$term, levels=rownames(ps.rowColSorted))

    g=ggplot(p.df, aes(x = term, y = logP, group = drug, color = drug)) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      labs(title = "Enrichment P-Value Across drug-associated differential signals",
           x = "GO biological processes",
           y = "-log10 P-Value")
    print(g)
    
  }
    # if(nrow(hm.term.foldEnrich.df)<=1) next
    # 
    # #print(c(ot, type, hm, nrow(hm.term.foldEnrich.df)))
    # buf=spread(hm.term.foldEnrich.df[, c("cond", "term", "p")], term, p, fill=1)
    # hm.term.p.matrix=as.matrix(cbind(buf[, -1]))
    # rownames(hm.term.p.matrix)=buf[,1]
    # colnames(hm.term.p.matrix)=colnames(buf)[-1]
    # pheatmap(-log10(t(hm.term.p.matrix[,terms.slct])), 
    #          color = cols,
    #          cluster_row=T,
    #          cluster_col=F,
    #          main=ot)
    
 
  
  
  
}

dev.off()


#pdf(file.ou.pdf, width=15, height=20)

# 
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
#   
# }



save(enrichRes.HM,
     enrichRes.HMs.merge,
     file=file.ou.RData)






