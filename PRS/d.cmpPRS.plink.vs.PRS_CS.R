
R.LIBS.PATH = .libPaths()
.libPaths(c("~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/", R.LIBS.PATH))


options(scipen=999)


#compare PRS with EHR data and BD subgroups
library(pheatmap)


CLUMP.LD.R2.CUTOFF=0.4
CLUMP.DIS.CUTOFF=100000
PEAK.WIND=2000

SNP.PVAL.CUTOFF="pval_0.05"
set.seed(666)

subGRP2COL = brewer.pal(10, "Set3")[1:5]
names(subGRP2COL)=paste0("subgrp", 1:5)


dir.in.prs.1=paste0("c_PRS_byPlink_PGC.leaveOutMayo.GWAS_genomicRegions/ld", CLUMP.LD.R2.CUTOFF, ".clumpDis", CLUMP.DIS.CUTOFF/1000, "kb.wind", PEAK.WIND/1000, "kb/")
file.in.PRS.1.RData = paste0(dir.in.prs.1, "plink_PRS_genomicRegions.RData")

dir.in.prs.2=paste0("c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/PRS_wind", PEAK.WIND/1000, "kb/")
file.in.PRS.2.RData = paste0(dir.in.prs.2, "plink_PRS_CS_genomicRegions.RData")



dir.ou=paste0("d_PRS.plink_vs_PRS.CS/wind", PEAK.WIND/1000, "kb/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.PRS.corr.pdf = paste0(dir.ou, "PRS.corr.pdf")


load(file.in.PRS.1.RData)

plink.cutoff2PRS.list=lapply(samp2PRS.list,
FUN=function(df)
{
  df[, grepl("PRS.z", names(df))]
})



load(file.in.PRS.2.RData)
PRS_CS.PRS.df=samp2PRS.df[, grepl("PRS.z", names(samp2PRS.df))]
PRS_CS.PRS.mat=as.matrix(PRS_CS.PRS.df)
colnames(PRS_CS.PRS.mat)=paste0("PRS_CS.", colnames(PRS_CS.PRS.mat))

samps.all=rownames(PRS_CS.PRS.mat)
corr.list=lapply(plink.cutoff2PRS.list,
FUN=function(plink.prs.df)
{
  plink.prs.mat=as.matrix(plink.prs.df[samps.all,])
  colnames(plink.prs.mat)=paste0("C+T.", colnames(plink.prs.mat))

  corr=cor(plink.prs.mat, PRS_CS.PRS.mat)

})

pdf(file.ou.PRS.corr.pdf, height=10, width=10)
for(nm in names(corr.list))
{
  pheatmap(corr.list[[nm]],
          cluster_rows=F,
          cluster_cols=F,
          main=nm)

  f.o.corr.tsv = paste0(dir.ou, "PRS.corr.", nm, ".tsv")
  write.table(corr.list[[nm]], f.o.corr.tsv, sep="\t", quote=F)
}

dev.off()



