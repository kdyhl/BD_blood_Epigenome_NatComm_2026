#

library(UpSetR)


QTL.Q.CUTOFF="emp.p.fdr.cutoff0.05"

HMs = c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3",  "H3K27me3")
files.in.gCREs = paste0("a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData")
names(files.in.gCREs) = HMs
file.in.merged2Pks.RData = "a_6_mergePeaksFromHMs_hg38/5HMs.mergedPeak2HMPeak.hg38.RData"

dir.ou=paste0("b_gCRE_sharing_upset/")
dir.create(dir.ou, recursive=F, showWarnings=F)
file.ou.RData= paste0(dir.ou, "5HMs.gCRE.sharing.RData")
file.ou.pdf= paste0(dir.ou, "5HMs.gCRE.sharing.upset.pdf")

load(file.in.merged2Pks.RData)

hm.merged.gCRE=lapply(HMs,
FUN=function(hm)
{
	load(files.in.gCREs[hm])
	gCREs=hQTLPeaks.list[[QTL.Q.CUTOFF]]

	levels(factor(hm.peak2mergedPeak[[hm]][gCREs]))
})
names(hm.merged.gCRE)= HMs

save(hm.merged.gCRE,
	file=file.ou.RData)


for(hm in HMs)
{
	f.o.hm=paste0(dir.ou, hm, ".merged.gCRE.txt")
	write(hm.merged.gCRE[[hm]], f.o.hm)
}

pdf(file.ou.pdf, width=12, height=6)

upset(fromList(hm.merged.gCRE), order.by = "freq")

dev.off()