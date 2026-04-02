#V1.2
#c.2.cmbPotentialGenes_V1.2_bulkAndcellsortedeQTL.R
#add bulk eQTL colocalization
#add DEG annotation

#
library(ggplot2)
library(ggrepel)
library(pheatmap)

library(RColorBrewer)
cols <- brewer.pal(9, 'OrRd')

COLOC.H4.CUTOFF=0.5
COLOC.H4.LOW.CUTOFF=0.1

file.in.gNeighb.gLinks.RData="c_geneNeighb_visGWASANDBulkAnnot_byGene_GeneticLink_V1.1_3HM//3HM.geneNB.peaks.annot_gLinks.RData"
file.in.gNeighb.HiC.RData= "../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.annot.RData"
#2HM
# file.in.gNeighb.gLinks.RData="c_geneNeighb_visGWASANDBulkAnnot_byGene_GeneticLink_V1.2.1_meta.2HM.DEG/3HM.geneNB.peaks.annot_gLinks.RData"
# file.in.gNeighb.HiC.RData= "../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.2.1_meta.2HM.DEG/geneNB.annot.RData"

#5cond
# file.in.gNeighb.gLinks.RData="c_geneNeighb_visGWASANDBulkAnnot_byGene_GeneticLink_V1.2.2_meta.5conds.DEG/3HM.geneNB.peaks.annot_gLinks.RData"
# file.in.gNeighb.HiC.RData= "../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.2.2_meta.5conds.DEG/geneNB.annot.RData"

file.in.meta.DEG.RData="../hiC.DiffPeakandDEGandTF/a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds/DEG.metaDP.DEGs.RData"

# file.in.Brain.eQTL.GWAS.coloc.RData = "../eQTLvsGWAS/a_GTExeaQTLvsGWAS_bycoloc/Brain_1000k/Brain.eQTLvsGWAS.RData"
# file.in.Blood.eQTL.GWAS.coloc.RData = "../eQTLvsGWAS/a_GTExeaQTLvsGWAS_bycoloc/Blood_1000k/Blood.eQTLvsGWAS.RData"
files.in.eQTL.GWAS.coloc.RData = list()
files.in.eQTL.GWAS.coloc.RData$Brain= system(paste0("find ../eQTLvsGWAS/a_BrainscRNA.eQTLvsGWAS_bycoloc/*_1000k/*eQTLvsGWAS.RData"), intern = T)
names(files.in.eQTL.GWAS.coloc.RData$Brain) = sapply(files.in.eQTL.GWAS.coloc.RData$Brain, FUN=function(x){gsub(".*/(.*?)_1000k/.*", "\\1", x)})
files.in.eQTL.GWAS.coloc.RData$Brain = c(files.in.eQTL.GWAS.coloc.RData$Brain,
                                         bulk="../eQTLvsGWAS/a_GTExeaQTLvsGWAS_bycoloc/Brain_1000k/Brain.eQTLvsGWAS.RData")
#
files.in.eQTL.GWAS.coloc.RData$Blood= system(paste0("find ../eQTLvsGWAS/a_Blueprint.eQTLvsGWAS_bycoloc/*/*.eQTLvsGWAS.RData"), intern = T)
names(files.in.eQTL.GWAS.coloc.RData$Blood) = sapply(files.in.eQTL.GWAS.coloc.RData$Blood, FUN=function(x){gsub(".*/(.*?)_1000k/.*", "\\1", x)})
files.in.eQTL.GWAS.coloc.RData$Blood= c(files.in.eQTL.GWAS.coloc.RData$Blood,
                                        bulk="../eQTLvsGWAS/a_GTExeaQTLvsGWAS_bycoloc/Blood_1000k/Blood.eQTLvsGWAS.RData")

files.in.eQTL.GWAS.coloc.RData=unlist(files.in.eQTL.GWAS.coloc.RData)


#dir.ou="c_2_targetGenes_comb_V1.1_cellsortedEQTL_meta.2HM.DEG/"
dir.ou="c_2_targetGenes_comb_V1.2_bulkAndcellsortedEQTL/"
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.txt=paste0(dir.ou, "cellsorted.candidateGenes.evidence.txt")
file.ou.driver.txt=paste0(dir.ou, "215drivergenes.cellsorted.coloc.txt")
file.ou.driver.braineQTL.txt=paste0(dir.ou, "drivergenes.cellsorted.coloc.braineQTL.pp4_0.5.txt")
file.ou.driver.bloodeQTL.txt=paste0(dir.ou, "drivergenes.cellsorted.coloc.bloodeQTL.pp4_0.5.txt")
file.ou.RData=paste0(dir.ou, "cellsorted.candidateGenes.evidence.RData")
file.ou.pdf=paste0(dir.ou, "cellsorted.candidateGenes.evidence.pdf")

#
genes.candid =list()
load(file.in.gNeighb.gLinks.RData)
genes.candid$gLinks = genes.filt.sorted

load(file.in.gNeighb.HiC.RData)
genes.candid$HiC = genes.filt.sorted

load(file.in.meta.DEG.RData)

# gESGid.noVersion2gName = gESGid2gName
# names(gESGid.noVersion2gName) = gsub("\\..*", "", names(gESGid.noVersion2gName) )

for(nm in names(files.in.eQTL.GWAS.coloc.RData))
{
  f.i.RData = files.in.eQTL.GWAS.coloc.RData[nm]
  load(f.i.RData)
  genes.candid[[nm]] = sapply(eGenes.eQTLvsGWAS, FUN=function(x){x$coloc$summary["PP.H4.abf"]})
  names(genes.candid[[nm]]) = names(eGenes.eQTLvsGWAS)
}
  

# load(file.in.Brain.eQTL.GWAS.coloc.RData)
# genes.candid$Brain.coloc = sapply(eGenes.eQTLvsGWAS, FUN=function(x){x$coloc$summary["PP.H4.abf"]})
# names(genes.candid$Brain.coloc) = names(eGenes.eQTLvsGWAS)
# 
# load(file.in.Blood.eQTL.GWAS.coloc.RData)
# genes.candid$Blood.coloc = sapply(eGenes.eQTLvsGWAS, FUN=function(x){x$coloc$summary["PP.H4.abf"]})
# names(genes.candid$Blood.coloc) = names(eGenes.eQTLvsGWAS)


candiGenes.evid = t(sapply(names(gESGid2gName),
FUN=function(gID)
{
  gID.noV = gsub("\\..*", "", gID)
  gName = gESGid2gName[gID]
  
  g.isgLinks=0
  g.isHiC =0
  g.cellType.coloc=rep(0, length(files.in.eQTL.GWAS.coloc.RData))
  names(g.cellType.coloc) = names(files.in.eQTL.GWAS.coloc.RData)
  #g.Blood.coloc=0
  
  if(gID %in% genes.candid$gLinks)
  {
    g.isgLinks=1
  }
  if(gName %in% genes.candid$HiC)
  {
    g.isHiC=1
  }
  
  for(nm in names(files.in.eQTL.GWAS.coloc.RData))
  {
    if(grepl("bulk", nm))
    {
      if(gID %in% names(genes.candid[[nm]]))
      {
        g.cellType.coloc[nm]=genes.candid[[nm]][gID]
      }
    }else
    {
      if(gID.noV %in% names(genes.candid[[nm]]))
      {
        g.cellType.coloc[nm]=genes.candid[[nm]][gID.noV]
      }  
    }
    
  }

  
  
  return(c(g.isgLinks, 
           g.isHiC,
           g.cellType.coloc
           ))
  
  
}))
colnames(candiGenes.evid) = c("gLinks", "bulk.HiC",  names(files.in.eQTL.GWAS.coloc.RData))
rownames(candiGenes.evid) = names(gESGid2gName)

candiGenes.evid.df = data.frame(gid=rownames(candiGenes.evid),
                                gName=gESGid2gName[rownames(candiGenes.evid)],
                                candiGenes.evid,
                                stringsAsFactors = F,
                                check.names = F)

candiGenes.evid.df$Brain.eQTL.coloc.max=apply(candiGenes.evid[, grepl("Brain", colnames(candiGenes.evid))], 1, max)
candiGenes.evid.df$Blood.eQTL.coloc.max=apply(candiGenes.evid[, grepl("Blood", colnames(candiGenes.evid))], 1, max)


save(candiGenes.evid.df, file=file.ou.RData)
#

propAndCounts.coloc=sapply(c("Blood", "Brain"),
FUN=function(tiss)
{
  cond = paste0(tiss, ".eQTL.coloc.max")
  
  count=c(bg=sum(candiGenes.evid.df[,cond]>=0.5),
         bulk.HiC = sum(candiGenes.evid.df[,cond]>=0.5 & candiGenes.evid.df[,"bulk.HiC"]==1),
         gLinks = sum(candiGenes.evid.df[,cond]>=0.5 & candiGenes.evid.df[,"gLinks"]==1))
  prop=100*count/c(nrow(candiGenes.evid.df),
                  sum(candiGenes.evid.df[,"bulk.HiC"]==1),
                  sum(candiGenes.evid.df[,"gLinks"]==1))
  
  c(count, prop)
})
colnames(propAndCounts.coloc) = c("Blood", "Brain")

prop.test.res=list()
prop.test.res$braineQTL.hicVSBg=prop.test(x=c(propAndCounts.coloc["bulk.HiC", "Brain"], propAndCounts.coloc["bg", "Brain"]),
                                          n=c(sum(candiGenes.evid.df[,"bulk.HiC"]==1), nrow(candiGenes.evid.df)), 
                                          alternative="greater")$p.val
prop.test.res$braineQTL.gLinkVSBg=prop.test(x=c(propAndCounts.coloc["gLinks", "Brain"], propAndCounts.coloc["bg", "Brain"]),
                                          n=c(sum(candiGenes.evid.df[,"gLinks"]==1), nrow(candiGenes.evid.df)), 
                                          alternative="greater")$p.val
prop.test.res$bloodeQTL.hicVSBg=prop.test(x=c(propAndCounts.coloc["bulk.HiC", "Blood"], propAndCounts.coloc["bg", "Blood"]),
                                          n=c(sum(candiGenes.evid.df[,"bulk.HiC"]==1), nrow(candiGenes.evid.df)), 
                                          alternative="greater")$p.val
prop.test.res$bloodeQTL.gLinkVSBg=prop.test(x=c(propAndCounts.coloc["gLinks", "Blood"], propAndCounts.coloc["bg", "Blood"]),
                                          n=c(sum(candiGenes.evid.df[,"gLinks"]==1), nrow(candiGenes.evid.df)), 
                                          alternative="greater")$p.val
    #


#
candiGenes.evid.filt.df= candiGenes.evid.df[apply(candiGenes.evid.df[, -(1:2)], 1, FUN=function(x) {any(x>=0.0001)}), ]

# candiGenes.evid.filt.df = data.frame(gid=rownames(candiGenes.evid.filt),
#                                      gName=gESGid2gName[rownames(candiGenes.evid.filt)],
#                                      candiGenes.evid.filt,
#                                      stringsAsFactors = F,
#                                      check.names = F)


write.table(candiGenes.evid.filt.df, file.ou.txt,sep="\t", quote=F, row.names = F)
#

#

genes.filt.brain.eQTL=candiGenes.evid.df$gName[candiGenes.evid.df$Brain.eQTL.coloc.max>=0.5]
genes.filt.blood.eQTL=candiGenes.evid.df$gName[candiGenes.evid.df$Blood.eQTL.coloc.max>=0.5]
write(genes.filt.brain.eQTL, file.ou.driver.braineQTL.txt)
write(genes.filt.blood.eQTL, file.ou.driver.bloodeQTL.txt)



pdf(file.ou.pdf, width=9, height=9)

#


xs=barplot(propAndCounts.coloc[4:6,],
           col=c("gray", "firebrick2", "green4"),
           ylab="percentage of BP-coloc genes",
           beside = T)
text(x=as.vector(xs),
     y=as.vector(propAndCounts.coloc[4:6,]/2),
     labels=as.vector(propAndCounts.coloc[1:3,]))
legend("topleft", 
       fill=c("gray", "firebrick2", "green4"),
       legend = c("bg", "bulk.HiC", "gLinks")
)

f.o.colocGene.tsv=paste0(dir.ou, "coloc.gene.proportion.tsv")
write.table(propAndCounts.coloc[4:6,], file=f.o.colocGene.tsv, sep="\t", quote=F)


genes.filt=rownames(candiGenes.evid.filt.df)[candiGenes.evid.filt.df$gLinks!=0 | candiGenes.evid.filt.df$bulk.HiC!=0]

# plot(candiGenes.evid.filt.df[genes.filt, "Blood.eQTL.coloc"],
#      candiGenes.evid.filt.df[genes.filt, "Brain.eQTL.coloc"],
#      xlab="Blood.eQTL.coloc",
#      ylab="Brain.eQTL.coloc")

candiGenes.evid.filt.df.filt=candiGenes.evid.filt.df[genes.filt,]
write.table(candiGenes.evid.filt.df[genes.filt,], file.ou.driver.txt, sep="\t", quote=F, row.names = F)
candiGenes.evid.filt.df.filt$type = paste(c("", "genetic")[candiGenes.evid.filt.df.filt$gLinks+1], 
                                          c("", "Hi-C")[candiGenes.evid.filt.df.filt$bulk.HiC+1])
ggplot(candiGenes.evid.filt.df.filt, 
       aes(x= Blood.eQTL.coloc.max, 
           y = Brain.eQTL.coloc.max, 
           label=gName, 
           color=type)) +
  geom_point(size = 2) +
  geom_label_repel(
          aes(label = gName),
          data = subset(candiGenes.evid.filt.df.filt, 
                        candiGenes.evid.filt.df.filt$Blood.eQTL.coloc.max>=COLOC.H4.CUTOFF |candiGenes.evid.filt.df.filt$Brain.eQTL.coloc.max>=COLOC.H4.CUTOFF),
          cex=5,
          force=12)+
  theme_bw()


f.o.candiGeneScatterPlot.tsv=paste0(dir.ou, "candiGeneScatterPlot.tsv")
write.table(candiGenes.evid.filt.df.filt, file=f.o.candiGeneScatterPlot.tsv, sep="\t", quote=F)


#
candiGenes.evid.filt.df.filt2 = candiGenes.evid.filt.df.filt[candiGenes.evid.filt.df.filt$Blood.eQTL.coloc.max>=COLOC.H4.CUTOFF,] # |candiGenes.evid.filt.df.filt$Brain.eQTL.coloc.max>=COLOC.H4.CUTOFF, ]
annotation_row=data.frame(is.meta.DEG= c(0,1)[rownames(candiGenes.evid.filt.df.filt2) %in% DEGs.filt + 1],
  mock=0)
rownames(annotation_row)=rownames(candiGenes.evid.filt.df.filt2)
gs.sharedWithBrain = rownames(candiGenes.evid.filt.df.filt2)[candiGenes.evid.filt.df.filt2[, "Brain.eQTL.coloc.max"]>=COLOC.H4.LOW.CUTOFF]



res=pheatmap(as.matrix(candiGenes.evid.filt.df.filt2[,names(files.in.eQTL.GWAS.coloc.RData)]),
         labels_row = candiGenes.evid.filt.df.filt2$gName,
         cluster_rows = T,
         cluster_cols = T,
         annotation_row=annotation_row, 
         col=cols)

gs.clu.ordered = rownames(candiGenes.evid.filt.df.filt2)[res$tree_row$order]
gs.reordered = c(intersect(gs.clu.ordered, gs.sharedWithBrain),
                setdiff(gs.clu.ordered, gs.sharedWithBrain))

candiGenes.evid.filt.ordered=as.matrix(candiGenes.evid.filt.df.filt2[,names(files.in.eQTL.GWAS.coloc.RData)])[gs.reordered, ]
rownames(candiGenes.evid.filt.ordered)=candiGenes.evid.filt.df.filt2[gs.reordered, "gName"]
pheatmap(candiGenes.evid.filt.ordered,
         labels_row = candiGenes.evid.filt.df.filt2[gs.reordered, "gName"],
         cluster_rows = F,
         cluster_cols = T,
         annotation_row=annotation_row[gs.reordered, ], 
         col=cols)

#
f.o.candiGeneheatmap.tsv=paste0(dir.ou, "candiGene.heatmap.tsv")
write.table(candiGenes.evid.filt.ordered, 
          file=f.o.candiGeneheatmap.tsv, sep="\t", quote=F)


dev.off()






