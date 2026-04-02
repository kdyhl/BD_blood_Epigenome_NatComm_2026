#V1.1
#a.2.cmpDPs.geneRegDomainAndGeneBody_V1.1_rmOvlp.R
#remove those CREs counted in the gene bodies
#use effect size not pvalue to check consistency

#compared DP signal from regulatory domain
#and DP signal at gene body

#1. meta analysis on DP regulatory domain, consistent with Hi-C interaction density
#2. meta analysis on DP regulatory domain, consistent with meta analysis on gene body




library(meta)
library(pheatmap)
library(RColorBrewer)

COND2PCH=c(
  #gb_H3K27ac=2,
  gb=2,
  grd=3
  #grd_H3K36me3=3
)
HM2COL=c(
  #gb_H3K27ac="#E5803C",
  
  H3K36me3="#5FAA37",
  H3K27me3="#4474EA",
  H3K27ac="#E5803C",
  H3K4me1="#E974FF",
  H3K4me3="#DD4545"
  )




file.in.geneRegDomain.annot.RData= "a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.annot.RData"
file.in.geneRegDomain.link.RData= "a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.link.RData"
#file.in.geneBody.meta.RData = "a_geneBody_DP_meta_V1.2_filtByDist/DP.V1.4.1_disease/genebody.DP.meta.RData"
file.in.geneBody.meta.RData = "a_geneBody_DP_meta_V1.2.1_filtByDist_allGenes/DP.V1.4.1_disease/genebody.DP.meta.RData"

file.in.eGenes = "/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/Whole_Blood.v8.egenes.txt.gz"


dir.ou="a_2_cmpDP.geneRegDomainVSgeneBody_V1.1_rmOvlpPk/DP.V1.4.1_disease.geneBody_meta_V1.2.1/"
dir.create(dir.ou, recursive = T, showWarnings = F)
file.ou.RData=paste0(dir.ou, "geneRegDomainVSgeneBody.metaDP.RData")
file.ou.pdf=paste0(dir.ou, "geneRegDomainVSgeneBody.metaDP.pdf")
#
load(file.in.geneRegDomain.annot.RData)
load(file.in.geneRegDomain.link.RData)
load(file.in.geneBody.meta.RData)


#estimate differential signal for overall regulatory domain by meta analysis

eGenes.info = read.table(file.in.eGenes, sep="\t", row.names = 1, header = T, stringsAsFactors = F)
# esmGid.v2gname=eGenes.info[, "gene_name"]
# names(esmGid.v2gname) = rownames(eGenes.info)
gname2esmGid=rownames(eGenes.info)
names(gname2esmGid)=eGenes.info[, "gene_name"]


gNames.ovlp=intersect(names(gname2esmGid), names(gene2Neighb))
                      
gene.hm.regDomainDP.meta=lapply(gNames.ovlp,                             
FUN=function(gName)
{
  hicRegs=gene2Neighb[[gName]]
  gid=gname2esmGid[gName]
  print(gName)
  #print(length(hicRegs))
  hm.hiCReg.dp.meta=lapply(names(hic.reg2HM.Pks),
  FUN=function(hm)
  {
    #print(hm)
    reg2pks = lapply(hicRegs,
    FUN=function(hicReg)
    {
      hic.reg2HM.Pks[[hm]][[hicReg]]
    })
    reg2pks.filt= reg2pks[sapply(reg2pks, length)!=0]
    if(length(reg2pks.filt)==0) return(NA)
    
    pks.filt = levels(factor(unlist(reg2pks.filt)))
    if(gid %in% names(HM.geneBody2DP.meta[[hm]]) && !is.na(HM.geneBody2DP.meta[[hm]][[gid]]))
    {
      pks.inGeneBody = HM.geneBody2DP.meta[[hm]][[gid]]$data$.studlab  
      pks.filt=setdiff(pks.filt, pks.inGeneBody)
    }
    
    if(length(pks.filt)==0) return(NA)
    reg2DP.log2fc = HM.DP.score.hg19[[hm]][pks.filt, "log2FoldChange"] #* gb.pks.hg38.ovlp.MR.crct
    reg2DP.pval = HM.DP.score.hg19[[hm]][pks.filt, "pvalue"]
    reg2DP.DP.log2fc.ste = abs(reg2DP.log2fc/qnorm(reg2DP.pval/2))
    
    
    
    m = metagen(TE=reg2DP.log2fc, 
                seTE=reg2DP.DP.log2fc.ste, 
                studlab = pks.filt,
                comb.fixed = F, 
                comb.random=T, 
                sm = "SMD")
   
    
    # if(hm=="H3K4me1" || hm=="H3K4me3")
    # {
    #   print(hm)
    #   print(paste(reg2DP.log2fc, sep=","))
    #   print(paste(reg2DP.pval, sep=","))
    #   print(pval.meta) 
    # }
    
  }) 
  names(hm.hiCReg.dp.meta) = names(hic.reg2HM.Pks)
  
  return(hm.hiCReg.dp.meta)
})
names(gene.hm.regDomainDP.meta) = gNames.ovlp
  

save(gene.hm.regDomainDP.meta,
     file=file.ou.RData)


#z
gene.hm.regDomainDP.z = t(sapply(gene.hm.regDomainDP.meta,
FUN=function(g.meta)
{
  hm.hiCReg.dp.meta.z=sapply(names(hic.reg2HM.Pks),
  FUN=function(hm)
  {
    m= g.meta[[hm]]
    if(is.na(m)) return(NA)
   
    effect.meta=m$TE.random
    pval.meta=m$pval.random
    z = abs(qnorm(pval.meta/2)) * sign(effect.meta)
    
  })
  
  names(hm.hiCReg.dp.meta.z) = names(hic.reg2HM.Pks)
  
  return(hm.hiCReg.dp.meta.z)
}))
rownames(gene.hm.regDomainDP.z) = gname2esmGid[gNames.ovlp]

#effect size
gene.hm.regDomainDP.TE = t(sapply(gene.hm.regDomainDP.meta,
FUN=function(g.meta)
{
  hm.hiCReg.dp.meta.TE=sapply(names(hic.reg2HM.Pks),
  FUN=function(hm)
  {
    m= g.meta[[hm]]
    if(is.na(m)) return(NA)
    
    effect.meta=m$TE.random
  })
 names(hm.hiCReg.dp.meta.TE) = names(hic.reg2HM.Pks)
 
 return(hm.hiCReg.dp.meta.TE)
}))
rownames(gene.hm.regDomainDP.TE) = gname2esmGid[gNames.ovlp]


save(gene.hm.regDomainDP.meta,
     gene.hm.regDomainDP.z,
     gene.hm.regDomainDP.TE,
     file=file.ou.RData)


#z
meta.z.geneBodyAndgeneRegDom = matrix(NA, ncol=10, nrow=length(gNames.ovlp))
colnames(meta.z.geneBodyAndgeneRegDom)= c(paste0("gb_", colnames(HM.geneBody2DP.meta.z.matrix)), 
                                          paste0("grd_", colnames(gene.hm.regDomainDP.z)))
rownames(meta.z.geneBodyAndgeneRegDom) = gname2esmGid[gNames.ovlp]

meta.z.geneBodyAndgeneRegDom[, paste0("grd_", colnames(gene.hm.regDomainDP.z))] = gene.hm.regDomainDP.z
gs.ovlp=intersect(rownames(HM.geneBody2DP.meta.z.matrix), rownames(meta.z.geneBodyAndgeneRegDom))
meta.z.geneBodyAndgeneRegDom[gs.ovlp, paste0("gb_", colnames(HM.geneBody2DP.meta.z.matrix))]=HM.geneBody2DP.meta.z.matrix[gs.ovlp,]
meta.z.geneBodyAndgeneRegDom[abs(meta.z.geneBodyAndgeneRegDom)==Inf]=NA


meta.z.geneBodyAndgeneRegDom.signed=meta.z.geneBodyAndgeneRegDom
is.repressive=grepl("H3K27me3", colnames(meta.z.geneBodyAndgeneRegDom.signed))
meta.z.geneBodyAndgeneRegDom.signed[, is.repressive]=-meta.z.geneBodyAndgeneRegDom.signed[, is.repressive]


z.cors.signed = cor(meta.z.geneBodyAndgeneRegDom.signed, 
                  use="pairwise.complete.obs",
                  method = "pearson")


#effect
meta.TE.geneBodyAndgeneRegDom = matrix(NA, ncol=10, nrow=length(gNames.ovlp))
colnames(meta.TE.geneBodyAndgeneRegDom)= c(paste0("gb_", colnames(HM.geneBody2DP.meta.TE.matrix)), 
                                          paste0("grd_", colnames(gene.hm.regDomainDP.TE)))
rownames(meta.TE.geneBodyAndgeneRegDom) = gname2esmGid[gNames.ovlp]

meta.TE.geneBodyAndgeneRegDom[, paste0("grd_", colnames(gene.hm.regDomainDP.TE))] = gene.hm.regDomainDP.TE
gs.ovlp=intersect(rownames(HM.geneBody2DP.meta.TE.matrix), rownames(meta.TE.geneBodyAndgeneRegDom))
meta.TE.geneBodyAndgeneRegDom[gs.ovlp, paste0("gb_", colnames(HM.geneBody2DP.meta.TE.matrix))]=HM.geneBody2DP.meta.TE.matrix[gs.ovlp,]
meta.TE.geneBodyAndgeneRegDom[abs(meta.TE.geneBodyAndgeneRegDom)==Inf]=NA


meta.TE.geneBodyAndgeneRegDom.signed=meta.TE.geneBodyAndgeneRegDom
is.repressive=grepl("H3K27me3", colnames(meta.TE.geneBodyAndgeneRegDom.signed))
meta.TE.geneBodyAndgeneRegDom.signed[, is.repressive]=-meta.TE.geneBodyAndgeneRegDom.signed[, is.repressive]


TE.cors.signed = cor(meta.TE.geneBodyAndgeneRegDom.signed, 
                    use="pairwise.complete.obs",
                    method = "pearson")



save(gene.hm.regDomainDP.meta,
     gene.hm.regDomainDP.z,
     gene.hm.regDomainDP.TE,
     meta.z.geneBodyAndgeneRegDom,
     meta.z.geneBodyAndgeneRegDom.signed,
     z.cors.signed,
     meta.TE.geneBodyAndgeneRegDom,
     meta.TE.geneBodyAndgeneRegDom.signed,
     TE.cors.signed,
     file=file.ou.RData)



#



meta.p.geneBodyAndgeneRegDom=2*pnorm(-abs(meta.z.geneBodyAndgeneRegDom.signed))

gb_H3K36me3.p=meta.p.geneBodyAndgeneRegDom[,"gb_H3K36me3"]
gs.gb_H3K36me3 = names(gb_H3K36me3.p)[!is.na(gb_H3K36me3.p)]
gs.sortedBy.gb_H3K36me3 = sort(gb_H3K36me3.p[gs.gb_H3K36me3], decreasing = T)
N=length(gs.sortedBy.gb_H3K36me3)
bin.N= 10
bin.gs_N = N/bin.N
g.bins = floor((1:N)/bin.gs_N)
g.bins[g.bins==bin.N] = bin.N-1               
g2bins.by.gb_H3K36me3.p = split(names(gs.sortedBy.gb_H3K36me3), 
                                g.bins)


pdf(file.ou.pdf, width=35, height=35)
layout(matrix(1:100, ncol=10))
g.bins2corr = sapply(colnames(meta.TE.geneBodyAndgeneRegDom.signed),
FUN=function(cond)
{
  y=meta.TE.geneBodyAndgeneRegDom.signed[, cond]
  sapply(g2bins.by.gb_H3K36me3.p,
  FUN=function(bin.gs)
  {
    smoothScatter(
      meta.TE.geneBodyAndgeneRegDom.signed[bin.gs,"gb_H3K36me3"],
      y[bin.gs],
      xlab="gb_H3K36me3",
      ylab=cond
      )
    abline(h=0, lty=2)
    abline(v=0, lty=2)
    cor(y[bin.gs],
        meta.TE.geneBodyAndgeneRegDom.signed[bin.gs,"gb_H3K36me3"],
        use="pairwise.complete.obs",
        method = "pearson")  
  })
  
})




plot(c(0,10),
     c(0,1),
     type="n",
     ylab="pearson correlation",
     xlab="bins of genes based on gb_H3K36me3 p value"
     )
cond2hm=c()
cond2type=c()
for(cond in colnames(g.bins2corr))
{
  hm = gsub("gb_|grd_", "", cond)
  type = gsub("_H.*?$", "", cond)
  lines(x=as.numeric(rownames(g.bins2corr)),
        y=g.bins2corr[, cond],
        type="b",
        pch=COND2PCH[type],
        col=HM2COL[hm])
  
  
  cond2hm[cond]=hm
  cond2type[cond]=type
}
legend("topleft",
       legend=colnames(g.bins2corr),
       lwd=1,
       col=HM2COL[cond2hm],
       pch=COND2PCH[cond2type],
       
)


f.o.bins.cors.tsv=paste0(dir.ou, "gb_H3K36me3.bin.corrs.tsv")
write.table(t(g.bins2corr),file=f.o.bins.cors.tsv, sep="\t", quote=F)




# 
# cols.pos <- brewer.pal(5, 'Reds')
# #cols.neg <- brewer.pal(3, 'Blues')
# breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
# cols=c("white", cols.pos)
# pheatmap(z.cors.signed,
#          breaks=breaks,
#          main="comparison of meta-analysis of differential signal over gene regulatory domain and gene body",
#          col=cols
# )
# 
# # 
# layout(matrix(1:25, ncol=5))
# 
# for(i in 1:(ncol(meta.z.geneBodyAndgeneRegDom.signed)-1))
# {
#   for(j in (i+1): ncol(meta.z.geneBodyAndgeneRegDom.signed))
#   {
#     smoothScatter(x=meta.z.geneBodyAndgeneRegDom.signed[,i],
#                   y=meta.z.geneBodyAndgeneRegDom.signed[,j],
#                   xlab=colnames(meta.z.geneBodyAndgeneRegDom.signed)[i],
#                   ylab=colnames(meta.z.geneBodyAndgeneRegDom.signed)[j],
#                   main=paste0(z.cors.signed[i,j], " zscore"))
#   }
# }
# 
# 
#TE
cols.pos <- brewer.pal(5, 'Reds')
#cols.neg <- brewer.pal(3, 'Blues')
breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
cols=c("white", cols.pos)
pheatmap(TE.cors.signed,
         breaks=breaks,
         main="comparison of meta-analysis of differential signal over gene regulatory domain and gene body",
         col=cols
)
# 
# # 
# layout(matrix(1:25, ncol=5))
# 
# for(i in 1:(ncol(meta.TE.geneBodyAndgeneRegDom.signed)-1))
# {
#   for(j in (i+1): ncol(meta.TE.geneBodyAndgeneRegDom.signed))
#   {
#     smoothScatter(x=meta.TE.geneBodyAndgeneRegDom.signed[,i],
#                   y=meta.TE.geneBodyAndgeneRegDom.signed[,j],
#                   xlab=colnames(meta.TE.geneBodyAndgeneRegDom.signed)[i],
#                   ylab=colnames(meta.TE.geneBodyAndgeneRegDom.signed)[j],
#                   main= paste(TE.cors.signed[i,j], " effect size"))
#   }
# }
# dev.off()
# 

dev.off()

f.o.TE.cors.tsv=paste0(dir.ou, "geneRegDomainVSgeneBody.heatmap.tsv")
write.table(TE.cors.signed,file=f.o.TE.cors.tsv, sep="\t", quote=F)
