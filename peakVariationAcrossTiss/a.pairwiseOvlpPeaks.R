#pairwise overlap across histone modification peaks


library(pheatmap)
HMs = c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3",  "H3K27me3")

files.in.peaks.bed = sapply(HMs,
FUN=function(hm)
{
  f.i = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", hm, ".diffPeak.BG.bed")
})
names(files.in.peaks.bed) = HMs

dir.tmp = "~/hptmp/MayoBipolar/peaksOvlpFromHMs/"
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.peaks.ovlp =  paste0(dir.tmp, "peaks.ovlp.bed")
#file.tmp.peaks.all.sorted = paste0(dir.tmp, "allTiss.peaks.sorted.bed")
#file.tmp.peaks.merge = paste0(dir.tmp, "allTiss.peaks.merged")

dir.ou = "./a_peaksOvlpBetweenHMs_hg19/"
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData = paste0(dir.ou, "hmPair.peakOvlp.hg19.RData")
file.ou.pdf = paste0(dir.ou, "hmPair.peakOvlp.prop.pdf")


hmPair.peakOvlp=list()
for(i in 1:4)
{
  for(j in (i+1):5)
  {
    hm.i=HMs[i]
    hm.j=HMs[j]
    cmd = paste0("intersectBed -wa -wb -a ", files.in.peaks.bed[hm.i], " -b ", files.in.peaks.bed[hm.j], ">", file.tmp.peaks.ovlp)
    system(cmd)
    
    buf=read.table(file.tmp.peaks.ovlp, sep="\t", header = F, row.names = NULL)
    df=data.frame(pk.i=paste0(buf[,1], ":", buf[,2], "-", buf[,3]),
                  pk.j=paste0(buf[,4], ":", buf[,5], "-", buf[,6]),
                  stringsAsFactors = F)
    
    hmPair.peakOvlp[[paste0(hm.i,"to", hm.j)]] = split(df$pk.j, df$pk.i)
    hmPair.peakOvlp[[paste0(hm.j,"to", hm.i)]] = split(df$pk.i, df$pk.j)
    
  }
}



hm.totalPks=sapply(files.in.peaks.bed,
FUN=function(f.i)
{
  buf=as.numeric(unlist(strsplit(system(paste0("wc -l ", f.i), intern = T), " "))[1])
})

hm.pkOvlp.size = sapply(HMs,
FUN=function(hm.i)
{
  sapply(HMs,
  FUN=function(hm.j)
  {
    if(hm.i==hm.j)
    {
      hm.totalPks[hm.i]
    }else
    {
      length(hmPair.peakOvlp[[paste0(hm.i,"to", hm.j)]])
    }
  })
})
rownames(hm.pkOvlp.size) = HMs



hm.pkOvlp.prop = sweep(hm.pkOvlp.size,2,hm.totalPks,FUN="/") 



save(hmPair.peakOvlp,
     hm.totalPks,
     hm.pkOvlp.prop,
     hm.pkOvlp.size,
     file=file.ou.RData)

ovlp.annot= matrix(paste0(as.vector(round(hm.pkOvlp.size/1000)), 
                          "k\n(", 
                          as.vector(round(hm.pkOvlp.prop*100)), 
                          "%)"), 
                  ncol=5)

pdf(file.ou.pdf)
pheatmap(hm.pkOvlp.prop,
         cluster_row=F,
         cluster_col=F,
         display_numbers=ovlp.annot,
         col=colorRampPalette(c("white", "red"))(10))
dev.off()


#data for visualization tsv
file.ou.tsv = paste0(dir.ou, "hmPair.peakOvlp.prop.tsv")
write.table(hm.pkOvlp.prop, file=file.ou.tsv, sep="\t", quote=F)
