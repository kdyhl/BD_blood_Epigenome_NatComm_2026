#V1.1
#a.5.multipleTestOnQTL_forPeer_V1.1_cutoffs.R
#multiple cutoff
#version for results from peer

#based on V8 GTEx pipeline 
#cis-eQTL mapping was performed using FastQTL [16]. The mapping window was defined as 1 Mb up- and down-stream of the transcription start site (TSS), and the adaptive permutation mode was used with the setting --permute 1000 10000. The phased VCF described in Section 2.4 was used, and all variants with minor allele frequency ≥0.01 across the 838 donors were included. The same set of variants was tested in all tissues. The beta distribution-extrapolated empirical P-values from FastQTL were used to calculate gene-level q-values [17] with a fixed P-value interval for the estimation of π0 (the ‘lambda’ parameter was set to 0.85). A false discovery rate (FDR) threshold of ≤0.05 was applied to identify genes with at least one significant cis-eQTL (“eGenes”). To identify the list of all significant variant-gene pairs associated with cis-eGenes, a genome-wide empirical P-value threshold,pt, was defined as the empiricalP-value of the gene closest to the 0.05 FDR threshold. pt was then used to calculate a nominal P-value threshold for each gene based on the beta distribution model (from FastQTL) of the minimum P-value distribution f(pmin) obtained from the permutations for the gene. Specifically, the nominal threshold was calculated asF−1(pt), where F−1 is the inverse cumulative distribution. For each gene, variants with a nominalP-value below the gene-level threshold were consideredsignificant and included in the final list of variant-gene pairs



args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) 
{
  stop("need one argument for organs, any of Heart, Brain, Lung, Muscle.\n", call.=FALSE)
} 
HM = args[1]
CIS.WIND.LEN = as.numeric(args[2])

# R.LIBS.MYPATHs = c(V3.5="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.5/",
#                 V3.4="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.4/",
#                 V3.3="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.3/")
# R.LIBS.PATH = .libPaths()
# 
# if(grepl("R version 3.5", R.version$version.string))
# {
#   .libPaths(c(R.LIBS.MYPATHs["V3.5"], R.LIBS.PATH))
# }
# if(grepl("R version 3.4", R.version$version.string))
# {
#   .libPaths(c(R.LIBS.MYPATHs["V3.4"], R.LIBS.PATH))
# }
# if(grepl("R version 3.3", R.version$version.string))
# {
#   .libPaths(c(R.LIBS.MYPATHs["V3.3"], R.LIBS.PATH))
# }

library(qvalue)

options(scipen = 999)

FDR.CUTOFF=0.05
FDR.CUTOFFs=c(0.05, 0.1, 0.2)
#Emp.P.CUTOFFs=c(0.005, 0.01, 0.05)

dir.in = paste0("a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_", CIS.WIND.LEN/1000, "k/")
files.in.permute = paste0(dir.in, HM, "_permutation/*.txt.gz")
files.in.nominal = system(paste0("find ", dir.in, HM, "_nominal/*.txt.gz"), intern=T)


file.in.peakHeight.RData =paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.heightsMean.depthGCNormed.RData") ##paste0("./a_haQTLCalling_covariate_V1.2.1_filtPeaks.gINT.peer/", HM, ".H3K27ac.bed")


dir.ou = paste0("a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_", CIS.WIND.LEN/1000, "k/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.perm.all.bed = paste0(dir.ou, HM, ".permutation.all.txt")
file.ou.hQTLPeak.bed = paste0(dir.ou, HM, ".haQTLPeak.fdr", FDR.CUTOFF, ".hg38.bed")
file.ou.bgPeak.bed =  paste0(dir.ou, HM, ".bgPeak.hg38.bed")
file.ou.hQTLPeak.RData = paste0(dir.ou, HM, ".haQTLPeak.fdr", FDR.CUTOFF, ".hg38.RData")
file.ou.hQTLPeak.cutoffs.RData = paste0(dir.ou, HM, ".haQTLPeak.cutoffs.hg38.RData")
file.ou.hQTLPeak.pdf = paste0(dir.ou, HM, ".haQTLPeak.fdr", FDR.CUTOFF, ".hg38.pdf")
file.ou.hQTLandPk.Pairs.tsv = paste0(dir.ou, HM, ".haQTLandPeak.pairs.emp.0.005.hg38.tsv")
file.ou.hQTLPeak.dist.pdf = paste0(dir.ou, HM, ".haQTLPeak.dist.fdr", FDR.CUTOFF, ".hg38.pdf")


#permutation 
cmd = paste0("zcat ", files.in.permute, ">", file.ou.perm.all.bed) #|awk {'print $1 \"\\t\" $3 \"\\t\" $4 \"\\t\" $6 \"\\t\" $8 \"\\t\" $9 \"\\t\" $11'}")
system(cmd)
emp.pvals = read.table(file.ou.perm.all.bed, sep=" ", head=F, row.names = 1, stringsAsFactors = F)

#emp.pvals = t(sapply(buf, FUN=function(x) unlist(strsplit(x, split="\t"))))
peaks = rownames(emp.pvals)
peaks.len= sapply(peaks,
FUN=function(pk)
{
  buf=as.numeric(unlist(strsplit(pk, split=":|-", perl=T))[2:3])
  buf[2]-buf[1]
})
names(peaks.len) = peaks


# emp.pvals=emp.pvals[,-1]  
# emp.pvals= matrix(as.numeric(emp.pvals), ncol=4, byrow = F)
colnames(emp.pvals) =c("variant#", "shape1", "shape2", "dummy", "topVariant", "distance", "topNomial.P", "slope", "emp.P.perm", "emp.P.model")
#rownames(emp.pvals) =peaks
emp.pvals.filtNA= emp.pvals[!is.na(emp.pvals[,4]),]


emp.qval.res = qvalue(emp.pvals.filtNA[,"emp.P.model"])
emp.pvals.filtNA = cbind(emp.pvals.filtNA, emp.pval.fdr= emp.qval.res$qvalues)

#different cutoffs
hQTLPeaks.list=list()
for(fdr.cutoff in FDR.CUTOFFs)
{
  hQTLPeaks.list[[paste0("emp.p.fdr.cutoff", fdr.cutoff)]] = rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.pval.fdr"]<=fdr.cutoff]
}
# for(emp.p.cutoff in Emp.P.CUTOFFs)
# {
#   hQTLPeaks.list[[paste0("emp.p.cutoff", emp.p.cutoff)]] = rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.P.model"]<=emp.p.cutoff]
# }
#
cutoff2emp.p.cutoff = sapply(FDR.CUTOFFs,
FUN=function(fdr.cutoff)
{
  hQTLPeaks = rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.pval.fdr"]<=fdr.cutoff]
  emp.pval.cutoff = sort(emp.pvals.filtNA[,"emp.P.model"], decreasing = F)[length(hQTLPeaks)]
})
names(cutoff2emp.p.cutoff) = paste0("emp.p.fdr.cutoff",FDR.CUTOFFs)
#cutoff2emp.p.cutoff[paste0("emp.p.cutoff", Emp.P.CUTOFFs)] = Emp.P.CUTOFFs


#for each hQTL peak, calculate nominal pvalue cutoff
hQTLPeaks2nominPCutoff.bycutoffs = lapply(names(cutoff2emp.p.cutoff),
FUN=function(nm.cutoff)
{
  print(nm.cutoff)
  emp.p.cutoff=cutoff2emp.p.cutoff[nm.cutoff]
  hQTLPeaks =  hQTLPeaks.list[[nm.cutoff]]
  print(length(hQTLPeaks))
  res=sapply(hQTLPeaks,
  FUN=function(pk)
  {
    shape1 = emp.pvals.filtNA[pk, "shape1"]
    shape2 = emp.pvals.filtNA[pk, "shape2"]
    nominP.cutoff = qbeta(emp.p.cutoff, shape1, shape2)
  })
  names(res) = hQTLPeaks
  
  return(res)
})
names(hQTLPeaks2nominPCutoff.bycutoffs) = names(cutoff2emp.p.cutoff)

#significant haQTLs 
hQTLpeak2signifQTL.bycutoffs = lapply(names(hQTLPeaks2nominPCutoff.bycutoffs),
FUN=function(nm.cutoff)
{
  hQTLPeaks = hQTLPeaks.list[[nm.cutoff]]
  hQTLPeaks2nominPCutoff = hQTLPeaks2nominPCutoff.bycutoffs[[nm.cutoff]]
  
  buf = lapply(files.in.nominal,
  FUN=function(f.i)
  {
    print(f.i)
    d= read.table(gzfile(f.i), sep=" ", row.names = NULL, head=F, stringsAsFactors = F)
    d.filt = d[d[,1] %in% hQTLPeaks,]
    pk2QTLs = split(d.filt, d.filt[,1])
    pk2QTLs.filt=lapply(pk2QTLs,
    FUN=function(qtls)
    {
      pk = qtls[1,1]
      p.cutoff=hQTLPeaks2nominPCutoff[pk]
      qtls[qtls[,4]<=p.cutoff,]
    })
    
  })
  hQTLpeak2signifQTL = do.call(c, buf)
  
  print(paste0("haQTLPeak:", length(hQTLPeaks), "; haQTLs:", sum(sapply(hQTLpeak2signifQTL, nrow))))
  
  return(hQTLpeak2signifQTL)
})
names(hQTLpeak2signifQTL.bycutoffs) =names(hQTLPeaks2nominPCutoff.bycutoffs)


save(emp.pvals.filtNA,
     hQTLPeaks.list,
     hQTLPeaks2nominPCutoff.bycutoffs,
     hQTLpeak2signifQTL.bycutoffs,
     cutoff2emp.p.cutoff,
     #hQTLPeaks,
     #hQTLpeak2signifQTL, 
     file=file.ou.hQTLPeak.cutoffs.RData)


#output haQTL 
df.slct = do.call(rbind, hQTLpeak2signifQTL.bycutoffs$emp.p.fdr.cutoff0.05)
colnames(df.slct)=c("gARE", "haQTL", "distance", "nominal-p", "beta")
write.table(df.slct, file=file.ou.hQTLandPk.Pairs.tsv, sep="\t", quote=F, row.names = F)




#selected cutoff 
hQTLPeaks = rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.pval.fdr"]<=FDR.CUTOFF]
emp.pval.cutoff = sort(emp.pvals.filtNA[,"emp.P.model"], decreasing = F)[length(hQTLPeaks)]
write(gsub(":|-", "\t", hQTLPeaks), file=file.ou.hQTLPeak.bed)
write(gsub(":|-", "\t", rownames(emp.pvals.filtNA)), file=file.ou.bgPeak.bed)


#distance distribution
#hist(do.call(rbind, hQTLpeak2signifQTL.bycutoffs[[1]])[,3], xlab="all signifcant haQTL to ARE")
pdf(file.ou.hQTLPeak.pdf, height = 9, width=12)
layout(matrix(1:6, ncol=2))
hist(emp.pvals.filtNA[, "distance"], xlab="distance of top haQTL to all ARE", main=paste0(CIS.WIND.LEN/1000, "kb"))
for(nm in rev(names(cutoff2emp.p.cutoff)))
{
  hist(emp.pvals.filtNA[emp.pvals.filtNA$emp.P.model<=cutoff2emp.p.cutoff[nm], "distance"], xlab="distance of top haQTL to gARE", main=paste0(CIS.WIND.LEN/1000, "kb\n", nm))
}
#hist(do.call(rbind, hQTLpeak2signifQTL.bycutoffs[[1]])[,3], xlab="distance of all haQTLs to gARE", main=paste0(CIS.WIND.LEN/1000, "kb\nFDR cutoff ", FDR.CUTOFF))


dev.off()



#lead haQTL distance
hQTLPk.dist = emp.pvals.filtNA[hQTLPeaks.list[[paste0("emp.p.fdr.cutoff", FDR.CUTOFF)]],"distance"]
hQTLPk.dist.prop= sum(abs(hQTLPk.dist)<=20000)/length(hQTLPk.dist)

leadQTL.dist=list()
leadQTL.dist$all=abs(emp.pvals.filtNA$distance)/1000
leadQTL.dist$hQTLPk=abs(hQTLPk.dist)/1000

ks.test.p=ks.test(leadQTL.dist$hQTLPk, leadQTL.dist$all, alternatice="less")$p.val
leadQTL.dist.cdf=lapply(leadQTL.dist, ecdf)

pdf(file.ou.hQTLPeak.dist.pdf, width=12, height=6)
layout(matrix(1:2, ncol=2))
hist(hQTLPk.dist/1000, 
     xlab="distance of lead hQTL to gCRE (kb)", 
     breaks=20, 
     main= paste0( signif(hQTLPk.dist.prop*100, 1), "% lead haQTL within 20kb"))
plot(leadQTL.dist.cdf$all, 
     xlab="distance of lead hQTL to gCRE (kb)",
     ylab="cumulative probablity",
     col="grey", main=paste0("pval ", ks.test.p))
lines(leadQTL.dist.cdf$hQTLPk, col="red")
legend("bottomright",
       legend=c("all", "gCREs (hQTL peaks)"),
       col=c("grey", "red"),
       lwd=1)

dev.off()

f.o.txt=paste0(dir.ou, HM, ".haQTLPeak.dist_kb.fdr", FDR.CUTOFF, ".hg38.txt")
write(hQTLPk.dist/1000, file=f.o.txt, sep="\n")

# #
# #load(file.in.peakHeight.RData)
# load(file.in.peakHeight.RData)
# pks.hg382hg19.uniq= readRDS(file.in.hg382hg19.RDS)
# peaks.hg19 = rownames(samp2HeightsMean$depth.GCnormTeng)
# pks.hg38tohg19 = pks.hg382hg19.uniq[pks.hg382hg19.uniq %in% peaks.hg19]
# 
# Y.hg38 = samp2HeightsMean$depth.GCnormTeng[pks.hg38tohg19, ]
# rownames(Y.hg38) = names(pks.hg38tohg19)
# Y.hg38.mean = apply(Y.hg38, 1, mean)
# #
# pdf(file.ou.hQTLPeak.pdf, height = 9, width=12)
# 
# layout(matrix(1:6, ncol=2, byrow = T))
# 
# hist(emp.pvals.filtNA[,"topNomial.P"],
#      xlab="top nominal pvalue for each peak",
#      main=HM)
# hist(emp.pvals.filtNA[,"emp.P.model"],
#      xlab="empirical pvalue for each peak",
#      main=HM)
# 
# hist(peaks.len, xlab="peak length", breaks=120)
# hist(peaks.len[hQTLPeaks], xlab="hQTL peaks length", breaks=40)
# 
# hist(Y.hg38.mean, xlab="peak mean height")
# hist(Y.hg38.mean[hQTLPeaks], xlab="hQTL peaks peak mean height")
# 
# hist(sapply(hQTLpeak2signifQTL, nrow), xlab="number of significant snps per peak",
#      breaks=100)
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
