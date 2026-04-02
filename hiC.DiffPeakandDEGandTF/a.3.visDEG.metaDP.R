#V1.2
#a.3.visDEG.metaDP_V1.2_5conds.R

#gb_H3K36me3
#gb_H3K27ac 
#gb_H3K27me3
#grd_H3K36me3
#grd_H3K27ac

#add annotated gene sets
#add UpSetR for DEG supported by different marks


#V1.1
#H3K36me3 as a main signature
#together with H3K27me3 in gene body or regulatory domain
#H3K27ac in regulatory domain
#need sign to be consisent with pvalue <=0.01

#identify final DEG list
#and visualize its DP signal in meta analysis

#re-scale gene body region



library(UpSetR)


DEG.P.permCUTOFF=0.01
DEG.Q.CUTOFF=0.05

COND2PCH=c(
          gb_H3K27ac=2,
          gb_H3K36me3=2,
          gb_H3K27me3=2,
          grd_H3K27ac=3,
          grd_H3K36me3=3
          #grd_H3K27me3=3
         
         )
COND2COL=c(
          gb_H3K27ac="#E5803C",
          gb_H3K36me3="#5FAA37",
          gb_H3K27me3="#4474EA",
          grd_H3K27ac="#E5803C",
          grd_H3K36me3="#5FAA37")
#
geneSets.Annot=list()
geneSets.Annot$NET=c("SELPLG", "TLR2", "TLR4", "TLR7", "TLR8", "ITGAM", "ITGB2", "AGER" ,"HMGB1", "CLEC6A", "CLEC7A", "SIGLEC5", "SIGLEC9", "SIGLEC14", "SIRL1", "CYBB",
            "FCGR1A", "FCGR1BP", "FCGR1CP", "FCGR3A", "FCGR2A", 
            "IRAK4", "RIPK1", "RIPK3", 
            "PIK3CD", "PIK3CA", "PIK3CG", "PIK3CB", "AKT1", "AKT2", "AKT3", "MTOR",
            "ATG7", "MPO", "PADI4", "ELANE", "CTSG" 
            ) #Venizelos Papayannopoulos, Nat. Rev. Immu., 2017

geneSets.Annot$TFs =c("FLI1", "IRF8", "JDP2", "SPI1", "RFX4", "E2F4", "ATF1", "MBD2", 
  "MECP2", "NRF1", "FEV", "EGR3", "CEBPB", "RFX2", "RORA", "ZBT14", "MYBL1", "NF1A")
geneSets.Annot$DNA.damage.repair=c("NEIL1", "MUTYH", "LIG3", "APEX1", "UNG", "OGG1", "PARP1", "XRCC1", "XRCC3", "LIG1", "MPG", "MUTYH", "NTHL1", "UNG") #from Kucuker ,2022 J Affect Disord, 2022
geneSets.Annot$Treg=c("IL10", "IL12A", "EBI3", "TGFB1", 
                      "PRF1", "GZMB", "GZMA", "GZMK", "GZMM", "GZMH", 
                      "IL2RA", "IL2RB", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "ENTPD1", "NT5E") #Romano, Front Immun, 2019
geneSets.Annot$TCellAnergy=c("LAG3", "PDCD1", "CTLA4", "CD28" ,"TRA", "IL2", "NFATC2" , "SIRT1", "IKZF1")  
geneSets.Annot$TCellExhaustion=c("LAG3", "PDCD1", "HAVCR2", "CD160", "TRA", "CD28", "CD244", "IL2", "KLRG1", "TNF", "IFNAR1", "IL10RA", "IL10RB", "CTLA4")  
geneSets.Annot$TSenescence=c("KLRG1", "B3GAT1", "CD160", "TRA", "CD28", "CDKN2A", "CDKN1A")
  
file.in.eGenes = "/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/Whole_Blood.v8.egenes.txt.gz"
file.in.hg19.gtf= "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.annotation.OnlyGene.gtf"

file.in.geneBody.meta.RData = "a_geneBody_DP_meta_V1.2.1_filtByDist_allGenes/DP.V1.4.1_disease/genebody.DP.meta.RData"
file.in.geneRegDom.meta.RData ="a_2_cmpDP.geneRegDomainVSgeneBody_V1.1_rmOvlpPk/DP.V1.4.1_disease.geneBody_meta_V1.2.1/geneRegDomainVSgeneBody.metaDP.RData"



dir.ou="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds/"
dir.create(dir.ou, recursive = T, showWarnings = F)
#file.ou.RData=paste0(dir.ou, "geneRegDomainVSgeneBody.metaDP.RData")
file.ou.pdf=paste0(dir.ou, "DEG.metaDP.byGene.hg19.pdf")
file.ou.RData=paste0(dir.ou, "DEG.metaDP.DEGs.RData")
files.ou.txt=c(up=paste0(dir.ou, "DEG.metaDP.DEGs.up.txt"),
               dn=paste0(dir.ou, "DEG.metaDP.DEGs.dn.txt"))

file.ou.upset.pdf=paste0(dir.ou, "DEG.metaDP.upset.hg19.pdf")
#
#eGenes.info = read.table(file.in.eGenes, sep="\t", row.names = 1, header = T, stringsAsFactors = F)
buf=read.table(file.in.hg19.gtf, sep="\t", header=F, row.names = NULL, stringsAsFactors = F)
buf= buf[grepl("^chr", buf[,1], perl=T),]
gENSAndName = t(sapply(buf[,5],
                       FUN=function(x)
                       {
                         gid = gsub(".*gene_id (.*?);.*", "\\1", x, perl=TRUE)
                         gid=gsub("_\\d+", "", gid, perl=T)
                         gname = gsub(".*gene_name (.*?);.*", "\\1", x, perl=TRUE)
                         return(c(gid,gname))
                       }))
gene2info=data.frame(chr=buf[,1],
                     start=buf[,2],
                     end=buf[,3],
                     strand=buf[,4],
                     geneName=gENSAndName[,2],
                     stringsAsFactors = F)
rownames(gene2info) = gENSAndName[,1]
# esmGid.v2gname=eGenes.info[, "gene_name"]
# names(esmGid.v2gname) = rownames(eGenes.info)
# gname2esmGid=rownames(eGenes.info)
# names(gname2esmGid)=eGenes.info[, "gene_name"]


#
load(file.in.geneBody.meta.RData)
load(file.in.geneRegDom.meta.RData)

#
#DEG.meta.p = pnorm(-abs(HM.geneBody2DP.meta.z.matrix))*2
#DEG.meta.q= apply(DEG.meta.p, 2, p.adjust, method="BH")

DEG.meta.p = pnorm(-abs(meta.z.geneBodyAndgeneRegDom.signed))*2
DEG.meta.q= apply(DEG.meta.p, 2, p.adjust, method="BH")

# DEGs.sign= apply(sign(meta.z.geneBodyAndgeneRegDom.signed[, names(COND2PCH)]), 1, 
# FUN=function(x)
# {
#   #return(x[1]==x[2] && x[2]!=x[3])
#   x=x[!is.na(x)]
  
#   if(length(x) >=2 && all(x==x[1])) return(x[1])
#   else return(0)
  
# })
DEGs.sign= sapply(rownames(DEG.meta.p), 
FUN=function(g)
{
  
  g.sign =sign(meta.z.geneBodyAndgeneRegDom.signed[g, names(COND2PCH)])
  
  g.sign.filt=g.sign[DEG.meta.p[g, names(COND2PCH)]<=DEG.P.permCUTOFF ] #filter by pvalue
   #return(x[1]==x[2] && x[2]!=x[3])
  g.sign.filt=g.sign.filt[!is.na(g.sign.filt)]
   
  if(length(g.sign.filt) >=2 && all(g.sign.filt==g.sign.filt[1])) return(g.sign.filt[1])
  else return(0)
 
})
names(DEGs.sign)=rownames(DEG.meta.p)
#DEGs.sign.isConsist[is.na(DEGs.sign.isConsist)] =F

DEGs.filt= rownames(DEG.meta.q)[DEGs.sign!=0 &
                                  (is.na(DEG.meta.p[, "gb_H3K36me3"]) | DEG.meta.p[, "gb_H3K36me3"] <= DEG.P.permCUTOFF)  &
                                  apply(DEG.meta.q[, names(COND2PCH)]<=DEG.Q.CUTOFF, 1, sum, na.rm=T)>=1]

DEGs.nm.filt=gene2info[DEGs.filt,  "geneName"]
write(DEGs.nm.filt[DEGs.sign[DEGs.filt]==1], file=files.ou.txt["up"])
write(DEGs.nm.filt[DEGs.sign[DEGs.filt]==-1], file=files.ou.txt["dn"])

for(nm in names(geneSets.Annot))
{
  print(paste(nm, paste(intersect(geneSets.Annot[[nm]], DEGs.nm.filt), collapse=";")))
}



#upset by different HM
DEGs.filt.byCond=lapply(names(COND2PCH),
FUN=function(cond)
{
  genes = DEGs.filt[DEG.meta.p[DEGs.filt, cond]<=DEG.P.permCUTOFF]
  gs.filt=genes[!is.na(genes)]

  write(gene2info[gs.filt,  "geneName"], 
  file=paste0(dir.ou, "DEGs.filt.supBy.", cond, ".pCUTOFF", DEG.P.permCUTOFF,".txt"))

  return(gs.filt)

})
names(DEGs.filt.byCond) = names(COND2PCH)

pdf(file.ou.upset.pdf, width=12, height=6)

upset(fromList(DEGs.filt.byCond), order.by = "freq")

dev.off()

#
f.o.upset.txt=paste0(dir.ou, "DEG.metaDP.upset.hg19.txt")
DEGs.filt=sapply(DEGs.filt.byCond, paste, collapse="\t")
DEGs.filt=paste(names(DEGs.filt.byCond), DEGs.filt, sep="\t")
write(DEGs.filt, f.o.upset.txt)


save(
  DEG.meta.p,
  DEG.meta.q,
  DEGs.sign,
  DEGs.filt,
  gene2info,
  DEGs.filt.byCond,
  file=file.ou.RData)

rescale.transform=function(x, xmin, xmax, g.start, g.end)
{
  #if(! (xmin<g.start && xmax>g.end)) return(x)
  if(xmin>g.start && xmax<g.end) return(x)
  #only rescale the gene body part
  if(x <g.start) return(x)
  if(x>=g.start && x<=g.end)
  {
    x=g.start+ (x-g.start)/(g.end-g.start)*((x.max-g.end + g.start-x.min)/2)
    return(x)
  }
  if(x>g.end)
  {
    x= g.start +(x.max-g.end + g.start-x.min)/2 + x-g.end
  }
}


DEGs.slct = intersect(unlist(geneSets.Annot), DEGs.nm.filt)
DEGs.slct.id = rownames(gene2info)[gene2info[,"geneName"] %in% DEGs.slct]

pdf(file.ou.pdf, width=10, height=40)
layout(matrix(1:12, ncol=1))
for(deg in DEGs.slct.id) #c("ENSG00000066336.11", "ENSG00000185697.16")
{
  print(deg)
  deg.nm = gene2info[deg, "geneName"]
  deg.tss=gene2info[deg, "start"]
  deg.tts=gene2info[deg, "end"]
  f.o.tsv=paste0(dir.ou, "DEG.metaDP.", gene2info[deg, "geneName"], ".hg19.tsv")
  
  if(is.na(gene2info[deg, "strand"])) next
  if(gene2info[deg, "strand"]=="-")
  {
    deg.tts=gene2info[deg, "start"]
    deg.tss=gene2info[deg, "end"]
  }
   
  cond2TEs.list=list()
  cond2pk.center.list=list()
  for(cond in names(COND2PCH))
  {
    hm=gsub("gb_|grd_", "", cond)
    if(grepl("gb", cond))
    {
      cond2TEs.list[[cond]]=HM.geneBody2DP.meta[[hm]][[deg]]$data$.TE
      pks=HM.geneBody2DP.meta[[hm]][[deg]]$data$.studlab
      cond2pk.center.list[[cond]]=sapply(pks,
      FUN=function(pk)
      {
        buf=unlist(strsplit(pk, split=":|-"))
        mean(as.numeric(buf[2:3]))
      })
    }else
    {
      if(!is.na(gene.hm.regDomainDP.meta[[deg.nm]][[hm]]))
      { 
        cond2TEs.list[[cond]]=gene.hm.regDomainDP.meta[[deg.nm]][[hm]]$data$.TE
        pks=rownames(gene.hm.regDomainDP.meta[[deg.nm]][[hm]]$data)

        cond2pk.center.list[[cond]]=sapply(pks,
        FUN=function(pk)
        {
          buf=unlist(strsplit(pk, split=":|-"))
          mean(as.numeric(buf[2:3]))
        })
      }
    }
    
  }



  
  #
  y.min= min(unlist(cond2TEs.list))
  y.max= max(unlist(cond2TEs.list))
  
  x.min=min(c(unlist(cond2pk.center.list), deg.tss, deg.tts), na.rm = T)
  x.max=max(c(unlist(cond2pk.center.list), deg.tss, deg.tts), na.rm=T)
  
  x.min.rescale=rescale.transform(x.min, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])
  x.max.rescale=rescale.transform(x.max, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])

  #
  cond.dp.df=data.frame(log2FC=unlist(cond2TEs.list),
                        chrom.pos=unlist(cond2pk.center.list),
                        COND=unlist(lapply(names(cond2TEs.list), FUN=function(x) rep(x, length(cond2TEs.list[[x]])))),
                        CRE=unlist(lapply(cond2pk.center.list, names)),
                        stringsAsFactors=F)
  cond.dp.df$x.rescale=sapply(cond.dp.df$chrom.pos, rescale.transform, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])

  x.rescale= sapply(cond2pk.center.list[[cond]], rescale.transform, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])
  write.table(cond.dp.df, file=f.o.tsv, sep="\t", quote=F)

  #
  plot(x=c(x.min.rescale, x.max.rescale)/1000000,
       y=c(y.min, y.max),
       main=paste0(deg, " ", gene2info[deg, "geneName"], "\n",
                   paste(paste(names(COND2PCH), signif(DEG.meta.p[deg, names(COND2PCH)],2) , sep = ": "), collapse=","), "\n",
                   gene2info[deg, "start"], " ", gene2info[deg, "end"], "\n",
                   x.min, " ", x.max),
       xaxt = 'n',
       xlab=paste(gene2info[deg, "chr"], "Mb"),
       ylab="log2FC",
       type="n")
  abline(h=0, lty=2)
  arrows(x0=rescale.transform(deg.tss, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])/1000000,
         x1=rescale.transform(deg.tts, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])/1000000,
         y0=0,
         y1=0,
         length=0.1,
         col="black",
         lwd=2)
  # for(cond in names(COND2PCH))
  # {
  #   if(sum(!is.na(cond2pk.center.list[[cond]]))!=0)
  #   {
  #     x.rescale= sapply(cond2pk.center.list[[cond]], rescale.transform, x.min, x.max, gene2info[deg, "start"], gene2info[deg, "end"])
  #     print(x.rescale)
  #     points(x=x.rescale/1000000,
  #          y=cond2TEs.list[[cond]],
  #          pch=COND2PCH[cond],
  #          col=COND2COL[cond])
  #   }
  # }
  points(x=cond.dp.df$x.rescale/1000000,
          y=cond.dp.df$log2FC,
           pch=COND2PCH[cond.dp.df$COND],
           col=COND2COL[cond.dp.df$COND])
  legend("bottomleft", 
         legend =  names(COND2PCH),
         pch=COND2PCH,
         col=COND2COL[names(COND2PCH)])
  
  
}  
  
dev.off()

