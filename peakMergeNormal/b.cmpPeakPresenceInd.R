#compare the peak presence in each individual samples


#args = commandArgs(trailingOnly=TRUE)
HMs=c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3")


N.CUTOFF=2 #need to appear in at least 10 percent of all samples
SUPER.TOTALREADs.CUTOFF=3e7
SUPER.RSC.CUTOFF=2
PEAK.LogQ.CUTOFF =2
GOOD.TOTALREADs.CUTOFF=2e7
GOOD.RSC.CUTOFF=1


dir.ou="b_cmpPeakPresenceAcrossInd/"
dir.create(dir.ou, showWarnings=F, recursive=T)
file.ou.RDS=paste0(dir.ou, "hms.peakPresence.RDS")
file.ou.pdf=paste0(dir.ou, "hms.peakPresence.pdf")


hm.peak2grp.prop.list=lapply(HMs,
FUN=function(hm)
{  
  print(hm)
  file.in.meta = paste0("~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/c_metaAndQC/", hm, ".sample.metaDataAndQC.txt")
  file.in.peaks.bg.bed=paste0("c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", hm, "/", hm, ".diffPeak.BG.bed")
  dir.in = paste0("/broad/compbio_ce/lhou/Mayo_Bipolar/AQUAS_", hm, "/")
  file.in.peak.suf = "out/peak/macs2/rep1/*nodup.tagAlign.pval0.01.500K.filt.narrowPeak.gz"

  dir.tmp = "~/hptmp/MayoBipolar/peakMerge/"
  dir.create(dir.tmp, showWarnings=F, recursive=T)
  file.tmp.ovlp.bed = paste(dir.tmp, hm, ".hg19.ovlp.bed", sep="")


  qc.info = read.table(file.in.meta, sep="\t", row.names=1, header = T)
  samps.QCed=rownames(qc.info)[qc.info$totalReads>=GOOD.TOTALREADs.CUTOFF & qc.info$RSC>=GOOD.RSC.CUTOFF & qc.info$isBad.byMismatch==0]
  grp2samps=split(samps.QCed, qc.info[samps.QCed, "Phenotype"])

  allPks=read.table(file.in.peaks.bg.bed, sep=",", header=F, row.names=NULL, stringsAsFactors=F)[,1]

  peak.pres.mat=sapply(samps.QCed,
  FUN=function(samp.id)
  {
    f.i.pk = paste0(dir.in, samp.id, "/", file.in.peak.suf)
    print(f.i.pk)
    #cmd = paste("zcat ", f.i.pk, " >>", file.tmp.allPeak)
    cmd = paste0("zcat ", f.i.pk, "|awk '($9>=", PEAK.LogQ.CUTOFF, "){print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $4 \"_\" $8}' | intersectBed -u -a ", file.in.peaks.bg.bed, " -b stdin >", file.tmp.ovlp.bed)
    system(cmd)


    pks.ovlp=read.table(file.tmp.ovlp.bed, header=F, row.names=NULL, sep=",", stringsAsFactors=F)[,1]

    pres=rep(0, length(allPks))
    names(pres)=allPks
    pres[pks.ovlp]=1
    return(pres)
  })


  peak2prop=apply(peak.pres.mat, 1, sum)/ncol(peak.pres.mat)

  peak2grp.prop=sapply(grp2samps,
  FUN=function(grp.samps)
  {
    apply(peak.pres.mat[, grp.samps], 1, sum)/length(grp.samps)
  })

  return(peak2grp.prop)

})
names(hm.peak2grp.prop.list)=HMs
saveRDS(hm.peak2grp.prop.list,
      file=file.ou.RDS)


pdf(file.ou.pdf, height=16, width=9)
layout(matrix(1:15, ncol=3, byrow=T))

for(hm in HMs)
{
  f.o.tsv=paste0(dir.ou, hm, ".peakPresence.BDvsCtrl.tsv")

  write.table(hm.peak2grp.prop.list[[hm]], file=f.o.tsv, sep="\t", quote=F)

  hist(hm.peak2grp.prop.list[[hm]][,"case"], main=paste0(hm, " case"), xlab="peak presence")
  hist(hm.peak2grp.prop.list[[hm]][,"control"], main=paste0(hm, " control"), xlab="peak presence")
  smoothScatter(hm.peak2grp.prop.list[[hm]], main=paste0(hm, " comparing peak presence\nbetween case and control"))
  abline(a=0, b=1, col="red", lty=2)



}
dev.off()
