#V1.4.1.2
#c.getDiff_byLimma_V1.4.1.2_sva_NoEHR_cmpBatch_Summary.R

#use sva to consider unwanted variables
#use permutated batches

# args = commandArgs(trailingOnly=TRUE)
# HM = args[1]
# #PERM.N=500

library(ggplot2)
library(gridExtra)

HMs=c("H3K27ac", 
          "H3K4me3", 
          "H3K4me1",
          "H3K36me3", 
          "H3K27me3")


P.ADJ.CUTOFF=0.05


dir.ou=paste0("./c_bulkDiffPeak_limma_V1.4.1.2_sva_noEHR_compBatch/")
file.ou.DP.num.pdf=paste0(dir.ou, "batch.DPs.num.pdf")  
file.ou.DP.t.pdf=paste0(dir.ou, "cmp.batch.DPs.t.pdf") 
file.ou.DP.t.RData=paste0(dir.ou, "batch.DPs.t.PCC.RData")  

pdf(file.ou.DP.num.pdf, height=4, width=20)
layout(matrix(1:5, nrow=1))  
for(hm in HMs)
{
  print(hm)
  d.o.hm=paste0(dir.ou, hm, "/")
  files.in.RData=dir(d.o.hm, "*RData", full.names=T)
  

  DPs.num=sapply(files.in.RData,
  FUN=function(f.i.RData)
  {
    load(f.i.RData)
    grp=unlist(batch.slct.grp)
    dps.num=sapply(res,
    FUN=function(df)
    {
      sum(df$adj.P.Val<=P.ADJ.CUTOFF)  
    })

    dps.num.SameGrp=c("same"=0, "diff"=0)


    if(grp=="case")
    {
      dps.num.SameGrp["same"]=dps.num["case"]
      dps.num.SameGrp["diff"]=dps.num["control"]
    }else
    {
      dps.num.SameGrp["same"]=dps.num["control"]
      dps.num.SameGrp["diff"]=dps.num["case"]
    }
    

    return(dps.num.SameGrp)
 
  })
  # print(length(DPs.num))
  # hist(log10(DPs.num+1), 
  #         breaks=100, 
  #         xlab="log10(peak Num+1)",
  #         main=paste0(hm, " ", signif(sum(DPs.num>HM2DP.NUM[hm])/length(DPs.num)*100,4), "% > observed dCREs"),
  #         )
  # abline(v=log10(HM2DP.NUM[hm]+1), col="red", lty=2)

  plot(t(DPs.num), 
      main=hm)
  abline(a=0, b=1, col="red", lty=2)
}
dev.off()


#


hm.pcc.list=list()
for(hm in HMs)
{
  print(hm)
  f.i.DP.RData=paste0("c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", hm, "/", hm, ".diffPeak.RData")
  

  load(f.i.DP.RData)

  d.o.hm=paste0(dir.ou, hm, "/")
  files.in.RData=dir(d.o.hm, "*RData", full.names=T)
  
  DP.ts.list=list()
  DP.ts.list$AD.vs.ctrl=list()
  DP.ts.list$same=list()

  for(f.i.RData in files.in.RData)
  {
    print(f.i.RData)
    load(f.i.RData)
    grp=unlist(batch.slct.grp)
    
    if(grp=="case")
    {
      DP.ts.list$AD.vs.ctrl=c(DP.ts.list$AD.vs.ctrl, list(res$control[pks.filt, "t"]))
      DP.ts.list$same=c(DP.ts.list$same, list(res$case[pks.filt, "t"]))
    }else
    {
      DP.ts.list$AD.vs.ctrl=c(DP.ts.list$AD.vs.ctrl, list(-res$case[pks.filt, "t"]))
      DP.ts.list$same=c(DP.ts.list$same, list(res$control[pks.filt, "t"]))
    }
  }
  DP.ts.list=lapply(DP.ts.list, 
  FUN=function(ts.list)
  {
    do.call(cbind, ts.list)
  })

  hm.pcc.list[[hm]]=lapply(DP.ts.list,
  FUN=function(ts.mat)
  {
    pcc.mat=cor(ts.mat, method="pearson")
    pcc.mat[upper.tri(pcc.mat)]
  })
 
}
save(hm.pcc.list,
    file=file.ou.DP.t.RData)


plot.list <- list()
for(hm in HMs)
{

  DP.ts.PCC.list=hm.pcc.list[[hm]]
  pcc.df=data.frame(PCC=c(DP.ts.PCC.list$AD.vs.ctrl, DP.ts.PCC.list$same),
                    cond=c(rep("BD.vs.ctrl", length(DP.ts.PCC.list$AD.vs.ctrl)), 
                          rep("same", length(DP.ts.PCC.list$same))))
  #
  p=ggplot(pcc.df, aes(PCC, fill = cond, color = cond)) +
    #geom_histogram(aes(y = ..density..), alpha = 0.4, position = "identity") +
    geom_density(alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("BD.vs.ctrl" = "red", "same" = "blue"))+
    ggtitle(hm)+
    theme_bw()

  plot.list[[hm]] <- p
  #print(p)
}



# Arrange and output to PDF
pdf(file.ou.DP.t.pdf, height = 3, width = 18)
grid.arrange(grobs = plot.list, ncol = 5)
dev.off()
