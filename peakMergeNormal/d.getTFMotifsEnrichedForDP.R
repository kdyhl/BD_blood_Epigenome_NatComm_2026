#V1.1
#d.getTFMotifsEnrichedForDP_V1.1_EpimapMotif.R
#use motif from epimap with threshold computed
#peaks are not further extend
#motifs are split into 10 sets




args = commandArgs(trailingOnly=TRUE)
MOD = args[1]
BG= args[2] #"randomBG"  #"localBG" #"#"givenBG" ##"

PK.WIN=0 #2000

CLU.SIZE.CUTOFF=50
NTASK=10

HMs = c("H3K27ac", "H3K4me1", "H3K36me3", "H3K4me3", "H3K27me3")
DP.TYPE="diseasedisease.Q0.05"

FILE.SH = "d.getTFMotifsEnrichForDP_byHomer.sh"



qsub.pref = c()
qsub.pref["norm"]= paste0("qsub -cwd -t 1-", NTASK, " -l h_vmem=20G -l h_rt=3:00:00 -j y -r y")
qsub.pref["rerun"]= paste0("qsub -cwd -t 1-", NTASK, " -l h_vmem=30G -l h_rt=8:00:00 -j y -r y")

dir.in = "c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/"
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
      if(as.numeric(buf[1])>CLU.SIZE.CUTOFF)
      {
        fs[type] = f.i
      }
    }
  }
  
  return(fs)
})
names(files.in.allPeaks) = HMs

files.in.allPeaks = files.in.allPeaks[sapply(files.in.allPeaks, length)!=1]


files.in.motif= system("find ~/lhou.compbio/data/Epimap/motifs/motif_forHomer/collated_motifs_bgFreq_homerFormat.pval_1e-4.10*.txt", intern=T)
names(files.in.motif) = sapply(files.in.motif,
FUN=function(f.i)
{
  gsub(".*10_(\\d+)\\.txt", "\\1", f.i)
})

dir.tmp=paste0("~/hptmp/MayoBipolar/d_motif_DP.V1.4.1_byHomer_V1.1_", PK.WIN/1000, "kb_", BG, "/")
dir.create(dir.tmp, showWarnings = F, recursive = T)

#
dir.ou= paste0("./d_TFMotifForPeak_DP.V1.4.1_V1.1_EpimapMotif_", PK.WIN/1000, "kb_", BG, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData = paste0(dir.ou, "DP.V1.4.1.RData")
  

#AREs by HM
pks.all.bed = list()
for(hm in HMs)
{
  pks.all.bed[[hm]]=sapply(files.in.allPeaks[[hm]],
  FUN=function(f.i)
  {
    buf=read.table(f.i, sep="\t", row.names = NULL, head=F)
    paste(buf[,1], buf[,2], buf[,3], sep="\t")
  })

}

#peak file
files.tmp.bed = lapply(names(pks.all.bed),
FUN=function(hm)
{
  print(hm)
  files.tmp = c()
  
  #nms.filt = names(pks.all.bed[[hm]])[!grepl(paste(TRAITS.RM, collapse="|"), names(pks.all.bed[[hm]]))]
  for(nm in names(pks.all.bed[[hm]]))
  {
    d = t(sapply(pks.all.bed[[hm]][[nm]], FUN=function(x){unlist(strsplit(x, split="\t", perl=T))}))  
    
    if(nrow(d)>=CLU.SIZE.CUTOFF)
    {
      files.tmp[nm] =  paste0(dir.tmp, hm, ".", nm, ".txt")
      d.forHomer = cbind(d[,1], as.numeric(d[,2])-PK.WIN, as.numeric(d[,3])+PK.WIN, paste(d[,1], d[,2], d[,3], sep="_"), 1, "+")
      write.table(d.forHomer, sep="\t", quote=F, file=files.tmp[nm], col.names = F, row.names = F)
    }
      
  
  }

  return(files.tmp)
  
})
names(files.tmp.bed) = names(pks.all.bed)



save(pks.all.bed,
     files.tmp.bed,
     file=file.ou.RData)
#
for(hm in names(files.tmp.bed))
{
  for(nm in names(files.tmp.bed[[hm]])[-1])
  {
    
    print(nm)
    d.o = paste(dir.ou, hm, ".", nm, "/", sep="")
    
    if(BG=="givenBG")
    {
      
      cmd = paste0(qsub.pref[MOD],
                   " -N ", paste0(hm, ".", nm),
                   " -o ", paste0(dir.ou, "out/"),
                   " -v MOD=", MOD, ",fg=", files.tmp.bed[[hm]][nm], ",genome=hg19,dir_ou=", d.o, ",bg=", files.tmp.bed[[hm]]["all.BG"],
                   " ", FILE.SH)
      
      
    }else if(BG=="randomBG")
    {
      cmd = paste0(qsub.pref[MOD],
                   " -N ", paste0(hm, ".", nm),
                   " -o ", paste0(dir.ou, "out/"),
                   " -v MOD=", MOD, ",fg=", files.tmp.bed[[hm]][nm], ",genome=hg19,dir_ou=", d.o, ",bg=", BG,
                   " ", FILE.SH)
      
    }
    
    print(cmd)
    system(cmd)
  }
}















