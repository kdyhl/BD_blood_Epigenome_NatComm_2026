#calculation the distribution of peaks around genes


HM2COL=c("H3K27ac"="#E5803C", 
         "H3K4me1"="#E974FF", 
         "H3K4me3"="#DD4545",
         "H3K36me3"="#5FAA37",
         "H3K27me3"="#4474EA"
         )

UPSTREAM.WIND=50000
DNSTREAM.WIND=50000
BIN.SIZE=500
BIN.NUM=UPSTREAM.WIND/BIN.SIZE
BIN.PROM.NUM=2000/BIN.SIZE
library(data.table)


file.in.hg19.gtf= "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.annotation.gtf"
files.in.peak=paste0("c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", names(HM2COL), "/", names(HM2COL), ".diffPeak.BG.bed")
names(files.in.peak) = names(HM2COL)


dir.tmp = "~/hptmp/b_bgPeak_homerAnnot/"
dir.create(dir.tmp, showWarnings = F, recursive = T)


dir.ou = "b_bgPeak_homerAnnot/"
dir.create(dir.ou, showWarnings = F, recursive = T)
files.ou.annot = paste0(dir.ou, names(HM2COL), ".bgPeaks.homerAnnot.txt")
names(files.ou.annot) = names(HM2COL)
file.ou.RData = paste0(dir.ou, "5HM.peakDistrib.scaledGeneBody.RData")
file.ou.pdf=paste0(dir.ou, "5HM.peakDistrib.scaledGeneBody.pdf")

#
for(hm in names(files.in.peak))
{
  print(hm)
  f.tmp= paste0(dir.tmp, hm, ".pk.tmp")
  
  if(file.exists(files.ou.annot[hm])) next
  
  buf=read.table(files.in.peak[hm], sep="\t", head=F, row.names = NULL)
  write(paste(paste0(buf[,1], ":", buf[,2], "-", buf[,3]), buf[,1], buf[,2], buf[,3], "+", sep="\t"), file=f.tmp)
  
  cmd=paste0("annotatePeaks.pl ", f.tmp, " hg19 -gtf ", file.in.hg19.gtf, ">", files.ou.annot[hm])
  
  system(cmd)
}

#
if(file.exists(file.ou.RData)) 
{
  load(file.ou.RData)
}else
{
  buf=fread(file.in.hg19.gtf, sep="\t", header=T, data.table=F, stringsAsFactors=F)
  gene.annots=buf[buf[,3]=="gene", ]
  rownames(gene.annots)=sapply(gene.annots[, 9],
  FUN=function(x)
  {
    gsub("gene_id \"(ENSG.*?)\".*", "\\1", x)
  })

  save(gene.annots,
      
       file=file.ou.RData)
}
  #
hm.distb.toGene=sapply(names(HM2COL),
FUN=function(hm)
{
  print(hm)
  buf=read.table(files.ou.annot[hm], sep="\t", header = F, row.names = 1, stringsAsFactors = F, skip=5)
  
  
  res=sapply(1:nrow(buf),
  FUN=function(i)
  {
    if(i %% 1000 ==0) print(i)
    pk.center=(buf[i,2]+buf[i,3])/2
    g=buf[i,11]
    if(g =="" || is.na(g) || !(g %in% rownames(gene.annots))) return(NA)
    
    g.start=gene.annots[g,4]
    g.end=gene.annots[g,5]
    g.strand=gene.annots[g,7]
    
    pk.bin=NA
    if(g.strand=="+")
    {
      if(pk.center<g.start)
      {
        if(g.start-pk.center>UPSTREAM.WIND) return (NA)
        pk.bin= paste0("up", round((g.start-pk.center)/BIN.SIZE))
      }else
      {
        if(pk.center>=g.start && pk.center<=g.end)
        {
          pk.bin= paste0("gb", round((pk.center-g.start)/((g.end-g.start)/BIN.NUM)))
        }else
        {
          if(pk.center-g.end>DNSTREAM.WIND) return (NA)
          pk.bin= paste0("dn", round((pk.center-g.end)/BIN.SIZE))
        }
      }
    }else
    {
      if(pk.center>g.end)
      {
        if(pk.center-g.end>UPSTREAM.WIND) return (NA)
        pk.bin= paste0("up", round((pk.center-g.end)/BIN.SIZE))
      }else
      {
        if(pk.center>=g.start && pk.center<=g.end)
        {
          pk.bin= paste0("gb", round((g.end-pk.center)/((g.end-g.start)/BIN.NUM)))
        }else
        {
          if(g.start-pk.center>DNSTREAM.WIND) return (NA)
          pk.bin= paste0("dn", round((g.start-pk.center)/BIN.SIZE))
          #pk.bin= round((g.start-pk.center)/BIN.SIZE)
        }
      }
    }
    
    return(pk.bin)
  })
  
  names(res) =rownames(buf)
  
  return(res)
})
save(gene.annots,
     hm.distb.toGene, 
     file=file.ou.RData)


pdf(file.ou.pdf, width=12, height=10)

layout(matrix(1:5, nrow=5))

barplot.list=list()

for(hm in names(HM2COL))
{
  x = hm.distb.toGene[[hm]]
  x.gbProp = sum(!is.na(x))/length(x)*100
  x = x[!is.na(x)]
  x[x=="up100"]="up99"
  x[x=="gb100"]="gb99"
  x[x=="dn100"]="dn99"
  bin2Num= summary(factor(x, levels=c(paste0("up", 99:0), paste0("gb", 0:99), paste0("dn", 0:99))), maxsum=3*BIN.NUM+1)
  
  bin.promot.avg= mean(bin2Num[c(paste0("up", 0:(BIN.PROM.NUM-1)), "gb0")])
  bin.gb.avg=mean(bin2Num[paste0("gb", 0:99)])
  
  xs=barplot(bin2Num, 
          las=2, 
          #lwd=0,
          #col=NA,
          border=HM2COL[hm],
          col=HM2COL[hm],
          xlab="genomic bins relative to genebody",
          ylab="counts of CREs",
          main=paste0(hm, " ", signif(x.gbProp,2), 
                      "% CREs within 50kb of a gene\ngeneBody to prom signal ratio: ", 
                      signif(bin.gb.avg/bin.promot.avg,2)))
  barplot.list[[hm]]=bin2Num
  abline(v=c(xs[c(101, 200)]), lty=2)
  
}
dev.off()

#data for visualization
barplot.df=do.call(rbind,barplot.list)
file.ou.tsv=paste0(dir.ou, "5HM.peakDistrib.scaledGeneBody.counts.tsv")
write.table(barplot.df, file=file.ou.tsv, sep="\t", quote=F)










# #print how much in each category
# for(hm in names(hm.distb.toGene))
# {
#   print(hm)
#   props=summary(factor(sub("\\d+", "", hm.distb.toGene[[hm]])))/length(hm.distb.toGene[[hm]])
#   print(props)
# }





