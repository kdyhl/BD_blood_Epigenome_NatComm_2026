#V1.3
#a.2.colocalizeTest.haQTLVSeQTL.coloc_V1.3.gARE_summary.R
  
if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



args = commandArgs(trailingOnly=TRUE)
TISS = args[1] #Brain Heart Muscle Lung
HM= args[2]
hQTL.WIND = as.numeric(args[3])

dir.tmp = paste0("~/hptmp/MayoBipolar/a.2_haQTLvseQTL_coloc_V1.3/", TISS, "_", HM, "/")
file.in.g2pk.RData = paste0("./a_gAREsNeareGene_emp.p.fdr.cutoff0.2/", TISS, "_1000k/", HM, ".eGene2gARE.hg38.RData")

cmd = paste0("find ", dir.tmp, TISS, ".", HM ,"*coloc_bycoloc_*.RData")
files.tmp.coloc.RData = system(cmd, intern = T)

dir.ou=paste0("a_2_coloc_haQTLvseQTL_bycoloc_V1.3_summary/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.pdf=paste0(dir.ou, TISS, "_", HM, "_", hQTL.WIND/1000, "k.haQTLvseQTL.pdf")
file.ou.RData=paste0(dir.ou, TISS, "_", HM, "_", hQTL.WIND/1000, "k.haQTLvseQTL.RData")

#load(file.tmp.g2pk.RData)

gene2peak.coloc.res=list()
for(f.t.RData in files.tmp.coloc.RData)
{
  load(f.t.RData)
  gene2peak.coloc.res= c(gene2peak.coloc.res, eQTLvshaQTL.coloc)
}


gene2peak.H4ABF = lapply(gene2peak.coloc.res,
FUN=function(peaks.abf)
{
  res=sapply(peaks.abf,
  FUN=function(pk.abf)
  {
    if(is.na(pk.abf))
    {
      return(NA)
    }else
    {
      return(pk.abf["PP.H4.abf"])
    }
  })
  names(res) = names(peaks.abf)
  return(res)
})

gene2peak.H4vsH3Ratio = lapply(gene2peak.coloc.res,
FUN=function(peaks.abf)
{
  res=sapply(peaks.abf,
  FUN=function(pk.abf)
  {
    if(is.na(pk.abf))
    {
      return(NA)
    }else
    {
      return(pk.abf["PP.H4.abf"]/pk.abf["PP.H3.abf"])
    }
  })
  names(res) = names(peaks.abf)
  return(res)
})


pdf(file.ou.pdf, width =30, height =4)
gene2peak.H4ABF.max= sapply(gene2peak.H4ABF, max, na.rm=T)
hist(gene2peak.H4ABF.max, breaks=100)

gene2peak.H4vsH3Ratio.max= sapply(gene2peak.H4vsH3Ratio, max, na.rm=T)
hist(log10(gene2peak.H4vsH3Ratio.max+ 1e-5), breaks=100)

boxplot(gene2peak.H4ABF[order(gene2peak.H4ABF.max)],
        las=2,
        ylab="H4 ABF for each gene-peak pair")
dev.off()

save(gene2peak.coloc.res,
     gene2peak.H4vsH3Ratio,
     gene2peak.H4ABF,
     file=file.ou.RData)


  