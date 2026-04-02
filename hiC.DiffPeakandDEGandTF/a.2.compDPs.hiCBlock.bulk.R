#bulk level for both DP and HiC (all links found)
#compare DPs within each hiC block
#based on network of hiC regions
#peaks from different kinds or mapped to this regions as annotation
#list genes asscciated with DP in each HiC-block
#annotate hiC Block with genetic signal



if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



options(scipen = 999)
library(igraph)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(DOSE)
DP.Q.CUTOFF=0.05
HIC.DIS.CUTOFF =1e6
HiC.REG.LENGTH.BIN=cbind(min=c(0, 1000, 2000, 5000, 10000),
                         max=c(1000, 2000, 5000, 10000, Inf))
HiCBlock.DP.ENRICH.ADJP.CUTOFF=0.05
PERM.NUM = 200

KEGG.FDR=0.05

HMs= c("H3K27ac", "H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")

DP.V = "DP.V1.4.1_sva_noEHR"
#
files.in.DP.RData = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
names(files.in.DP.RData)= HMs
file.in.hiC.RData="a_hiCLinks_withPeakAnnot_bulk_DP.BG_limma_V1.4.1/hiCLinks.RData"
file.in.gNeighbor.annot.RData = "a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.annot.RData"

dir.ou = paste0("a_2_cmpBulkDPAcrossHMs_hiCBlock/", DP.V, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.hicReg.RData = paste0(dir.ou, "Bulk.hicRegion.RData")  
file.ou.hicBlock.RData = paste0(dir.ou, "Bulk.hicBlocks.DPstat.RData")  
file.ou.hiCReg.pdf = paste0(dir.ou, "Bulk.hicReg.DPcomp.pdf")  
file.ou.hiCBlock.pdf = paste0(dir.ou, "Bulk.hicBlock.DPcomp.pdf") 
file.ou.hiCBlock.annot.pdf = paste0(dir.ou, "Bulk.hicBlock.annot.pdf")
#
# HM2DPStat = lapply(files.in.DP.RData,
# FUN=function(f.i.RData)
# {
#   load(f.i.RData)
#   pk.info = cbind(zscore=res@listData$stat, padj=res@listData$padj)
#   rownames(pk.info) = res@rownames
#   
#   return(pk.info)
# })

HM2DPStat = lapply(files.in.DP.RData,
FUN=function(f.i.RData)
{
 load(f.i.RData)
 pk.info = cbind(zscore= res$diseasedisease$t/abs(res$diseasedisease$t) * (-qnorm(res$diseasedisease$P.Value/2)), 
                 padj=res$diseasedisease$adj.P.Val)
 rownames(pk.info) = rownames(res$diseasedisease)
 
 return(pk.info)
})



#
load(file.in.hiC.RData)
load(file.in.gNeighbor.annot.RData)


#
pks.covered = sapply(HMs,
FUN=function(hm)
{
  pks = lapply(hicReg.df[, hm], FUN=function(x){unlist(strsplit(x, split=","))})
  levels(factor(unlist(pks)))
})
names(pks.covered)=HMs
barplot(signif(sapply(pks.covered,length)/sapply(HM2DPStat, nrow),2), 
        ylab="proportion covered by Hi-C region")
dev.off()

print("hm.hicReg2pks")
if(file.exists(file.ou.hicReg.RData)) load(file.ou.hicReg.RData)
  
if(!exists("hm.hicReg2pks"))
{
  hm.hicReg2pks=lapply(HMs,
  FUN=function(hm)
  {
    print(hm)
    hicReg2pks.list = lapply((1:nrow(hicReg.df))[hicReg.df[[hm]]!=""],
    FUN=function(i)
    {
      #print(i)
      pks=unlist(strsplit(hicReg.df[[hm]][i], split=","))
      res=data.frame(pks,
                 rownames(hicReg.df)[i],
                 stringsAsFactors = F)
      
    })
    hicReg2pks.df=do.call(rbind, hicReg2pks.list)
    
    hicReg2pks.df = hicReg2pks.df[hicReg2pks.df[,1] %in% rownames(HM2DPStat[[hm]]), ]
    hicReg2pks = split(hicReg2pks.df[,1], hicReg2pks.df[,2])
    
    #hicReg2pks=hicReg2pks[sapply(hicReg2pks, length)!=0]
    
  })
  names(hm.hicReg2pks)=HMs
  
  save(HM2DPStat,
       hm.hicReg2pks,
       file=file.ou.hicReg.RData)
}



print("hicReg.DPStat")
if(!exists("hicReg.DPStat"))
{
  hicReg.DPStat=sapply(HMs,
  FUN=function(hm)
  {
    print(hm)
    dp.stat = sapply(hm.hicReg2pks[[hm]],
    FUN=function(pks)
    {
      #print(Sys.time()-start_time)
      #pks.ovlp=intersect(pks, rownames(HM2DPStat[[hm]]))
      if(length(pks)==0)
      {
        return(NA)
      }
      #print(Sys.time()-start_time)
      zscores=HM2DPStat[[hm]][pks, "zscore"]
      zscores[which.max(abs(zscores))[1]]
      #print(Sys.time()-start_time)
    })  
    names(dp.stat) = names(hm.hicReg2pks[[hm]])
    dp.stat.all=rep(NA, nrow(hicReg.df))
    names(dp.stat.all) = rownames(hicReg.df)
    dp.stat.all[names(dp.stat)] = dp.stat
    
    return(dp.stat.all)
  })
  save(HM2DPStat,
       hm.hicReg2pks,
       hicReg.DPStat,
       file=file.ou.hicReg.RData)
}

print("hicReg.DP.adjp.list")
if(!exists("hicReg.DP.adjp"))
{
  hicReg.DP.adjp.list=lapply(HMs,
  FUN=function(hm)
  {
    print(hm)
    dp.adjp = t(sapply(hm.hicReg2pks[[hm]],
    FUN=function(pks)
    {
      pks.ovlp=intersect(pks, rownames(HM2DPStat[[hm]]))
      padj.min=c(UP=NA, DW=NA)
      
      if(length(pks.ovlp)==0)
      {
        return(padj.min)
      }
      zscores=HM2DPStat[[hm]][pks, "zscore"]
      if(any(zscores>0))
      {
        padj.min["UP"] = min(HM2DPStat[[hm]][pks[zscores>0], "padj"])
      }
      if(any(zscores<0))
      {
        padj.min["DW"] = min(HM2DPStat[[hm]][pks[zscores<0], "padj"])
      }
      return(padj.min)
      
    }))  
    colnames(dp.adjp) = paste0(hm, ".DP.", colnames(dp.adjp))
    
    dp.adjp.all=matrix(NA, ncol=2, nrow=nrow(hicReg.df))
    rownames(dp.adjp.all) = rownames(hicReg.df)
    colnames(dp.adjp.all) = colnames(dp.adjp)
    
    dp.adjp.all[rownames(dp.adjp), ] = dp.adjp
    
    
    return(dp.adjp.all)
  })
  hicReg.DP.adjp = do.call(cbind, hicReg.DP.adjp.list)

  save(
       HM2DPStat,
       hm.hicReg2pks,
       hicReg.DPStat,
       hicReg.DP.adjp,
       file=file.ou.hicReg.RData)
}

pdf(file.ou.hiCReg.pdf, width=20, height=20)
layout(matrix(1:25, ncol=5))
for(i in 1:nrow(HiC.REG.LENGTH.BIN))
{
  len.min=HiC.REG.LENGTH.BIN[i, "min"]
  len.max=HiC.REG.LENGTH.BIN[i, "max"]
  for(hm.a in HMs)
  {
    for(hm.b in HMs)
    {
      regs.ovlp = rownames(hicReg.DPStat)[(!is.na(hicReg.DPStat[,hm.a])) &
                                           (!is.na(hicReg.DPStat[,hm.b])) &
                                            hicReg.df$length>=len.min &
                                            hicReg.df$length<len.max 
                                           ]
      pcc = signif(cor(hicReg.DPStat[regs.ovlp,hm.a], hicReg.DPStat[regs.ovlp, hm.b]),3)
      smoothScatter(hicReg.DPStat[regs.ovlp,hm.a], hicReg.DPStat[regs.ovlp, hm.b], 
                    xlab=paste0(hm.a, " max bulk DP zscore"),
                    ylab=paste0(hm.b, " max bulk DP zscore"),
                    main=paste("bin: ", len.min, "-", len.max, " bp\n", length(regs.ovlp), " overlapped Hi-C regions\nPCC:", pcc)
                    )
    }
  }
}

dev.off()

# 
regNet.filt = subgraph.edges(regNet, eids =E(regNet)[E(regNet)$distance<=HIC.DIS.CUTOFF], delete.vertices = T)
#print(paste(dis.cutoff, length(V(regNet.sub)), length(E(regNet.sub)), sep=" "))
regNet.filt.cmpnts = components(regNet.filt)
hicBlock2hicRegs = split(names(regNet.filt.cmpnts$membership), paste0("hicBlock",regNet.filt.cmpnts$membership))

hiCBlock.topol.summary = t(sapply(hicBlock2hicRegs,
FUN=function(hicRegs)
{
  #hicRegs = names(regNet.filt.cmpnts$membership)[regNet.filt.cmpnts$membership==i]
  regNet.filt.sub = induced_subgraph(regNet.filt, hicRegs)
  nodes.N= length(hicRegs)
  edges.N = length(E(regNet.filt.sub))
  edge.prob = edges.N/(nodes.N*(nodes.N-1)/2)
  return(c(nodes.N= nodes.N,
           edges.N = edges.N,
           edge.prob = edge.prob
           ))
   
}))
#rownames(hiCBlock.topol.summary) = paste0("hicBlock", 1:regNet.filt.cmpnts$no)


# 
HM2DP.numAndRatio = sapply(HMs,
FUN=function(hm)
{
  pk.all.no = nrow(HM2DPStat[[hm]])
  dp.up.no = sum(HM2DPStat[[hm]][, "zscore"]>0 & HM2DPStat[[hm]][, "padj"]<=DP.Q.CUTOFF, na.rm=T)
  dp.dw.no = sum(HM2DPStat[[hm]][, "zscore"]<0 & HM2DPStat[[hm]][, "padj"]<=DP.Q.CUTOFF, na.rm=T)
  return(c(pk.all.no=pk.all.no,
           dp.up.no = dp.up.no,
           dp.dw.no = dp.dw.no,
           dp.up.ratio = dp.up.no/pk.all.no, 
           dp.dw.ratio= dp.dw.no/pk.all.no))
})
colnames(HM2DP.numAndRatio) = HMs

hiCBlock.DP.summary = t(sapply(hicBlock2hicRegs,
FUN=function(hicRegs)
{
  #print(i)
  #hicRegs = names(regNet.filt.cmpnts$membership)[regNet.filt.cmpnts$membership==i]
  # hm2pks = lapply(HMs,
  # FUN=function(hm)
  # {
  #   unlist(lapply(hicReg.df[hicRegs, hm], FUN=function(x) unlist(strsplit(x, split=","))))
  # })
  # names(hm2pks)=HMs
  hm2pks = lapply(HMs,
  FUN=function(hm)
  {
    hicRegs.ovlp= intersect(hicRegs, names(hm.hicReg2pks[[hm]]))
    
    if(length(hicRegs.ovlp)==0) return(NULL)
    
    levels(factor(unlist(hm.hicReg2pks[[hm]][hicRegs.ovlp])))
  })
  names(hm2pks)=HMs
  
  dp.summary = lapply(HMs,
  FUN=function(hm)
  {
    pks = hm2pks[[hm]]
    dp.up.no = sum(HM2DPStat[[hm]][pks, "zscore"]>0 & HM2DPStat[[hm]][pks, "padj"]<=DP.Q.CUTOFF, na.rm=T)
    dp.dw.no= sum(HM2DPStat[[hm]][pks, "zscore"]<0 & HM2DPStat[[hm]][pks, "padj"]<=DP.Q.CUTOFF, na.rm=T)
    dp.up.fishp = fisher.test(matrix(c(dp.up.no, 
                                       HM2DP.numAndRatio["dp.up.no", hm]-dp.up.no,
                                       length(pks)-dp.up.no, 
                                       HM2DP.numAndRatio["pk.all.no", hm]-length(pks)-HM2DP.numAndRatio["dp.up.no", hm]+dp.up.no), ncol=2))$p
    dp.dw.fishp = fisher.test(matrix(c(dp.dw.no, 
                                       HM2DP.numAndRatio["dp.dw.no", hm]-dp.dw.no,
                                       length(pks)-dp.dw.no, 
                                       HM2DP.numAndRatio["pk.all.no", hm]-length(pks)-HM2DP.numAndRatio["dp.dw.no", hm]+dp.dw.no), ncol=2))$p
    return(c(total=length(pks), 
             DP.UP.NO = dp.up.no, DP.DW.NO=dp.dw.no,
             DP.UP.enrichFC = dp.up.no/length(pks)/HM2DP.numAndRatio["dp.up.ratio", hm], 
             DP.UP.fisherP = dp.up.fishp, 
             DP.DW.enrichFC = dp.dw.no/length(pks)/HM2DP.numAndRatio["dp.dw.ratio", hm],
             DP.DW.fisherP = dp.dw.fishp))
  })
  names(dp.summary) = HMs
  
  return(unlist(dp.summary))
})) 
#rownames(hiCBlock.DP.summary) = paste0("hicBlock", 1:regNet.filt.cmpnts$no)

hiCBlock.DP.summary.adjp= sapply(colnames(hiCBlock.DP.summary)[grepl("fisherP", colnames(hiCBlock.DP.summary))],
FUN=function(nm)
{
  p.adjust(hiCBlock.DP.summary[,nm], method = "BH")
})
colnames(hiCBlock.DP.summary.adjp) = gsub("fisherP", "fisherAdjP", colnames(hiCBlock.DP.summary.adjp))

hiCBlock.DP.summary=cbind(hiCBlock.DP.summary,hiCBlock.DP.summary.adjp)

save(regNet.filt,
     regNet.filt.cmpnts,
     hicBlock2hicRegs,
     hiCBlock.topol.summary,
     HM2DP.numAndRatio,
     hiCBlock.DP.summary,
     hiCBlock.DP.summary.adjp,
     file=file.ou.hicBlock.RData)


#
pdf(file.ou.hiCBlock.pdf, width=15, height=20)

layout(matrix(1:10, ncol=2, byrow = T))
for(nm in colnames(hiCBlock.DP.summary)[grepl("fisherP", colnames(hiCBlock.DP.summary))])
{
  hist(hiCBlock.DP.summary[,nm],
       xlab=nm)
  
}

#
layout(matrix(1, ncol=, byrow = T))
hiCBlock.DP.FC.summary = hiCBlock.DP.summary[,grepl("FC", colnames(hiCBlock.DP.summary))]
hiCBlock.DP.FC.summary[is.na(hiCBlock.DP.FC.summary)]=0
hiCBlock.GWAS= t(sapply(rownames(hiCBlock.DP.FC.summary),
FUN=function(hicB)
{
  hicRegs = gsub(":|-", "\t", hicBlock2hicRegs[[hicB]])
  GWAS.ovlp = sum(hic.reg2HM.PK.annot.matrix[hicRegs, grepl("_GWAS", colnames(hic.reg2HM.PK.annot.matrix))])!=0
  GWAS.coloc = sum(hic.reg2HM.PK.annot.matrix[hicRegs, grepl("_coloc", colnames(hic.reg2HM.PK.annot.matrix))])!=0
  GWAS.MR = sum(hic.reg2HM.PK.annot.matrix[hicRegs, grepl("_MR", colnames(hic.reg2HM.PK.annot.matrix))])!=0
  
  return(as.numeric(c(GWAS.ovlp=GWAS.ovlp, GWAS.coloc= GWAS.coloc, GWAS.MR=GWAS.MR)))
}))
rownames(hiCBlock.GWAS) = rownames(hiCBlock.DP.FC.summary)
colnames(hiCBlock.GWAS) =c("GWAS.ovlp", "GWAS.coloc", "GWAS.MR")
pheatmap(log10(hiCBlock.DP.FC.summary+1),
         annotation_row = data.frame(hiCBlock.GWAS),
         cluster_cols = F)
#pvalue of enrichment
GWASEnrichDP.block.ps=apply(hiCBlock.GWAS, 2,
FUN=function(x)
{
  apply(hiCBlock.DP.FC.summary, 2,
  FUN=function(y)
  {
    mod=lm(log(y+1)~x)
    formatC(signif(summary(mod)$coef[2, "Pr(>|t|)"],2), format="e")
  })
})



# hiCBlock.DPEnrich.ovlp = apply(hiCBlock.DP.summary.adjp, 2,
# FUN=function(adjp.a)
# {
#   apply(hiCBlock.DP.summary.adjp, 2,
#   FUN=function(adjp.b)
#   {
#     block.a.no = sum(adjp.a<=HiCBlock.DP.ENRICH.ADJP.CUTOFF)
#     block.b.no = sum(adjp.b<=HiCBlock.DP.ENRICH.ADJP.CUTOFF)
#     block.ovlp.no = sum(adjp.a<=HiCBlock.DP.ENRICH.ADJP.CUTOFF & adjp.b<=HiCBlock.DP.ENRICH.ADJP.CUTOFF)
#     block.all.no =nrow(hiCBlock.DP.summary.adjp)
#     
#     fc = block.ovlp.no/block.a.no/(block.b.no/block.all.no)
#   })
# })
# diag(hiCBlock.DPEnrich.ovlp)= max(hiCBlock.DPEnrich.ovlp[upper.tri(hiCBlock.DPEnrich.ovlp)])
# 
# pheatmap(hiCBlock.DPEnrich.ovlp, 
#          cluster_rows = F,
#          cluster_cols = F)

#pdf(file.ou.hiCBlock.pdf, width=15, height=20)

layout(matrix(1:30, ncol = 5))
for(nm in colnames(hiCBlock.DP.FC.summary))
{
  plot(log10(hiCBlock.topol.summary[, "nodes.N"]+1),
       log10(hiCBlock.DP.FC.summary[, nm]+1),
       xlab="log10 nodes N",
       ylab= "DP log10 enrichment fold change",
       main=nm)
  plot(log10(hiCBlock.topol.summary[, "edges.N"]+1),
       log10(hiCBlock.DP.FC.summary[, nm]+1),
       xlab="log10 edges N",
       ylab= "DP log10 enrichment fold change",
       main= nm)

  #
  edge.prob=log10(hiCBlock.topol.summary[, "edge.prob"])
  log10FC = log10(hiCBlock.DP.FC.summary[, nm]+1)
  edge.prob.filt=edge.prob[log10FC!=0]
  log10FC.filt=log10FC[log10FC!=0]
  mod=lm(log10FC.filt~edge.prob.filt)

  plot(edge.prob,
       log10FC,
       #xlim=c(0, 0.2),
       xlab="log10 edge probability",
       ylab= "log10 DP enrichment fold change",
       main=paste0(nm, "\nP=", formatC(signif(summary(mod)$coeff[2, "Pr(>|t|)"],2), format="e")))
  abline(mod, col="red", lty=2)
}

dev.off()


#
file.ou.hiCBlock.DP.tsv = paste0(dir.ou, "Bulk.hicBlock.DP.tsv") 
file.ou.hiCBlock.topo.tsv = paste0(dir.ou, "Bulk.hicBlock.topo.tsv") 
file.ou.hiCBlock.GWAS.tsv = paste0(dir.ou, "Bulk.hicBlock.GWAS.tsv") 
write.table(hiCBlock.DP.FC.summary, file=file.ou.hiCBlock.DP.tsv, sep="\t", quote=F)
write.table(hiCBlock.topol.summary, file=file.ou.hiCBlock.topo.tsv, sep="\t", quote=F)
write.table(hiCBlock.GWAS, file=file.ou.hiCBlock.GWAS.tsv, sep="\t", quote=F)





#genes from enriched Hi-C block
hicBlock2DPAssGene=lapply(names(hicBlock2hicRegs),
FUN=function(hb.i)
{
  print(hb.i)
  if(all(hiCBlock.DP.summary.adjp[hb.i,]>HiCBlock.DP.ENRICH.ADJP.CUTOFF))
  {
    return(NA)
  }
  hicRegs = hicBlock2hicRegs[[hb.i]]
  
  hm.conds= colnames(hiCBlock.DP.summary.adjp)[hiCBlock.DP.summary.adjp[hb.i,]<=HiCBlock.DP.ENRICH.ADJP.CUTOFF]
  genes=lapply(hm.conds,
  FUN=function(hm.cond)
  {
    nm =gsub(".fisherAdjP","", hm.cond)   
    hicRegs.dp = hicRegs[hicReg.DP.adjp[hicRegs, nm]<=DP.Q.CUTOFF]
    hicRegs.dp = hicRegs.dp[!is.na(hicRegs.dp)]
    hicRegs.neighb = neighbors(regNet.filt, hicRegs.dp)
    hicRegs.neighb.genes=levels(factor(unlist(lapply(hicReg.df[union(names(hicRegs.neighb), hicRegs.dp),"TSS"], FUN=function(x){unlist(strsplit(x, split=","))}))))
  })
  names(genes) = hm.conds
  
  return(genes)
  
})
names(hicBlock2DPAssGene) = names(hicBlock2hicRegs)
hicBlock2DPAssGene= hicBlock2DPAssGene[!is.na(hicBlock2DPAssGene)]

hicBlock2DPAssGene.unionHM = lapply(hicBlock2DPAssGene,
FUN=function(genes)
{
  levels(factor(unlist(genes)))
})

HM2DPAssGene.unionHicBlock =list()

HM2DPAssGene.unionHicBlock$all = lapply(colnames(hiCBlock.DP.summary.adjp),
FUN=function(cond)
{
  levels(factor(unlist(lapply(hicBlock2DPAssGene[!is.na(hicBlock2DPAssGene)],
  FUN=function(genes)
  {
    if(cond %in% names(genes))
    {
      return(genes[[cond]])
    }
  }))))
})
names(HM2DPAssGene.unionHicBlock$all) = colnames(hiCBlock.DP.summary.adjp)
HM2DPAssGene.unionHicBlock$all=HM2DPAssGene.unionHicBlock$all[sapply(HM2DPAssGene.unionHicBlock$all, length)!=0]

for(nm in colnames(hiCBlock.GWAS))
{
  hicBs.nm = rownames(hiCBlock.GWAS)[hiCBlock.GWAS[,nm]==1]
  HM2DPAssGene.unionHicBlock[[nm]] = lapply(colnames(hiCBlock.DP.summary.adjp),
  FUN=function(cond)
  {
    levels(factor(unlist(lapply(hicBlock2DPAssGene[intersect(hicBs.nm, names(hicBlock2DPAssGene))],
    FUN=function(genes)
    {
      if(cond %in% names(genes))
      {
        return(genes[[cond]])
      }
    }))))
  })
  names(HM2DPAssGene.unionHicBlock[[nm]]) = colnames(hiCBlock.DP.summary.adjp)
  HM2DPAssGene.unionHicBlock[[nm]]=HM2DPAssGene.unionHicBlock[[nm]][sapply(HM2DPAssGene.unionHicBlock[[nm]], length)!=0]
}

save(regNet.filt,
     regNet.filt.cmpnts,
     hicBlock2hicRegs,
     hiCBlock.topol.summary,
     HM2DP.numAndRatio,
     hiCBlock.DP.summary,
     hiCBlock.DP.summary.adjp,
     hiCBlock.GWAS,
     hicBlock2DPAssGene,
     hicBlock2DPAssGene.unionHM,
     HM2DPAssGene.unionHicBlock,
     file=file.ou.hicBlock.RData)



#annotation
genes.bg = levels(factor(unlist(lapply(hicReg.df[,"TSS"], FUN=function(x){unlist(strsplit(x, split=","))}))))
genes.symbolAndgid = bitr(genes.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gsymbol2gid = genes.symbolAndgid$ENTREZID
names(gsymbol2gid) = genes.symbolAndgid$SYMBOL
gids.bg =genes.symbolAndgid$ENTREZID



dotplot.cluster = function(df, x="Cluster", y = "Description" , size.by="GeneRatio", col.by="pvalue", showCategory=30, title="", font.size=12)
{
  gsets = levels(factor(df[[x]]))
  gs2category = lapply(gsets,
  FUN=function(gs)
  {
    cat2pvalues= df[df[[x]] == gs, "pvalue"] 
    names(cat2pvalues) = df[df[[x]] == gs, y] 
    cat2pvalues.sorted = sort(cat2pvalues, decreasing = F)
    return(names(cat2pvalues.sorted)[1:min(showCategory, length(cat2pvalues.sorted))])
  })
  gs2category.all= levels(factor(unlist(gs2category)))
  
  gsVScate.matrix = matrix(0, ncol=length(gsets), nrow=length(gs2category.all))
  colnames(gsVScate.matrix) = gsets
  rownames(gsVScate.matrix) = gs2category.all
  
  for(i in 1:nrow(df))
  {
    print(i)
    if(df[i,y] %in% gs2category.all &&  df[i,x] %in% gsets)
    {
      gsVScate.matrix[df[i,y],as.character(df[i,x])]=1
    }
  }
  
  if(nrow(gsVScate.matrix)<=2)
  {
    cate.order=1:nrow(gsVScate.matrix)
  }else
  {
    cate.order = hclust(dist(gsVScate.matrix, method="euclidean"))$order
  }
  
  if(ncol(gsVScate.matrix)<=2)
  {
    gset.order=1:ncol(gsVScate.matrix)
  }else
  {
    gset.order = hclust(dist(t(gsVScate.matrix), method="euclidean"))$order
  }
  
  df.filt = df[df[[y]] %in% gs2category.all,]
  df.filt[[x]]= factor(df.filt[[x]], levels=gsets[gset.order])
  df.filt[[y]]= factor(df.filt[[y]], levels=gs2category.all[cate.order])
  
  p <- ggplot(df.filt, aes_(x = as.formula(paste("~", x, sep="")), y = as.formula(paste("~", y, sep="")), size = as.formula(paste("~", size.by, sep=""))))
  p <- p + geom_point() +
            aes_string(color=col.by) +
            scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
  p <- p + xlab("") + ylab("") + ggtitle(title) +
            theme_dose(font.size)
  print(p)
  #return(p)
}


# 
# hicBlock2DPAssGid.unionHM = lapply(hicBlock2DPAssGene.unionHM,
# FUN=function(gs)
# {
#   gids = gsymbol2gid[gs]
#   gids=gids[!is.na(gids)]
# })

HM2DPAssGid.unionHicBlock = lapply(HM2DPAssGene.unionHicBlock,
FUN=function(gs.cond)
{
  lapply(gs.cond,
  FUN=function(gs)
  {
    gids = gsymbol2gid[gs]
    gids=gids[!is.na(gids)]  
  })
  
})

# compKEGG.hicBlock  = tryCatch(
# {
#   compareCluster(geneCluster   = hicBlock2DPAssGid.unionHM,
#                      fun           = "enrichKEGG",
#                      pvalueCutoff  = KEGG.FDR,
#                      universe = gids.bg,
#                      pAdjustMethod = "BH")
# 
# },
# error = function(e)
# {
#   print(e)
#   return(NA)
#   
# })
#names(compKEGG.hicBlock) = names(hicBlock2DPAssGid.unionHM)




pdf(file.ou.hiCBlock.annot.pdf, width=15, height=12)

for(nm in names(HM2DPAssGid.unionHicBlock))
{
  compKEGG.HM= tryCatch(
  {
    compareCluster(geneCluster   = HM2DPAssGid.unionHicBlock[[nm]],
                       fun           = "enrichKEGG",
                       pvalueCutoff  = KEGG.FDR,
                       universe = gids.bg,
                       pAdjustMethod = "BH")
  
  },
  error = function(e)
  {
    print(e)
    return(NA)
    
  })
  #names(compKEGG.HM) = names(HM2DPAssGid.unionHicBlock)
  
  compKEGG.HM@compareClusterResult$GeneRatio = sapply(compKEGG.HM@compareClusterResult$GeneRatio, FUN=function(x){buf=unlist(strsplit(x,split="/", fixed=T)); as.numeric(buf[1])/as.numeric(buf[2])})
  compKEGG.HM@compareClusterResult$Cluster = sapply(compKEGG.HM@compareClusterResult$Cluster,
  FUN=function(nm)
  {
    nm = gsub(".fisherAdjP", "", nm)
    nm = sub(".", "\n", nm, fixed = T)
  })
  
  
  
  
  #compKEGG.merge.sub = compKEGG.merge[grepl(paste0(dis,"."), compKEGG.merge$Cluster, fixed=T), ]
  dotplot.cluster(compKEGG.HM@compareClusterResult, 
                  showCategory = 100, 
                  title = paste0(nm, " HiC-blocks\n",
                                 "KEGG Pathway Enrichment Analysis for genes associated with DP"))
  
  
}


dev.off()
  
  
