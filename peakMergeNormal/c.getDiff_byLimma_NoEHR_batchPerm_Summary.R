#V1.4.1.1
#c.getDiff_byLimma_V1.4.1.1_sva_NoEHR_batchPerm_Summary.R

#use sva to consider unwanted variables
#use permutated batches

# args = commandArgs(trailingOnly=TRUE)
# HM = args[1]
# #PERM.N=500


HM2DP.NUM=c("H3K27ac"=9029+13338, 
          "H3K4me3"=1468+97, 
          "H3K4me1"=728+77,
          "H3K36me3"=2783+13363, 
          "H3K27me3"=321+1324)


P.ADJ.CUTOFF=0.05

dir.ou=paste0("./c_bulkDiffPeak_limma_V1.4.1.1_sva_noEHR_batchPerm/")
file.ou.pdf=paste0(dir.ou, "perm.DPs.num.pdf")  
file.ou.RDS=paste0(dir.ou, "perm.DPs.num.RDS")  

hm.permDP.list=list()

pdf(file.ou.pdf, height=4, width=20)
layout(matrix(1:5, nrow=1))
for(hm in names(HM2DP.NUM))
{
  print(hm)
  d.o.hm=paste0(dir.ou, hm, "/")
  file.in.RDS=dir(d.o.hm, "*RDS", full.names=T)
  
  

  DPs.num=sapply(file.in.RDS,
  FUN=function(f.i.RDS)
  {
    df=readRDS(f.i.RDS)
    sum(df$adj.P.Val<=P.ADJ.CUTOFF)
  })
  print(length(DPs.num))

  #
  hm.permDP.list[[hm]]=DPs.num

  hist(log10(DPs.num+1), 
          xlim=c(0, ceiling(max(c(log10(DPs.num+1), log10(HM2DP.NUM[hm]+1))))+1),
          breaks=100, 
          xlab="log10(#dCRE+1)",
          main=paste0(hm, " ", signif(sum(DPs.num>HM2DP.NUM[hm])/length(DPs.num)*100,4), "% > observed dCREs"),
          )
  abline(v=log10(HM2DP.NUM[hm]+1), col="red", lty=2)
}

dev.off()


saveRDS(hm.permDP.list,
        file=file.ou.RDS)
#
f.o.tsv=paste0(dir.ou, "hms.perm.DPs.num.tsv")
hm.permDP=sapply(hm.permDP.list, paste, collapse="\t")
write(hm.permDP, file=f.o.tsv)
