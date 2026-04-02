#V1.2.3
#filter peaks by median level first before calculate CV

#V1.2.1
#filter peaks based on signal proportion
#use INT istead of quantile

#V1.2
#use Peer instead of RUV_g


#before calling haQTL, RUV is applied to identify covriates 
#the basic idea is that those factors that shared across different chromosomes should probably not be due to genetics, and thus should be removed for QTL calling
#step 1, select the top 100 variable peaks for each chromsome
#step 2, for each chromosome, using top variable peaks from other chromsomes as negative control peaks to feed in RUV-g
#step 3, check with known variables to see whether it is consistent


args = commandArgs(trailingOnly=TRUE)
HM = args[1] #
TOPPEAK.NUM=args[2] 

TOPCOV.NUM = 50 #args[2] #50
print(HM)


SIGNAL.CUTOFF=1
SIGNAL.PROP = 0.2


PEER.SCRIPT = "~/lhou.compbio/codes/mylib/RCode/run_PEER.R"

CHRs= paste0("chr", c(1:22, "X"))
#TOPPEAKS.PERCHR = 100



sh.head = paste("#!/bin/bash",
               "#$ -binding linear:1",
                #"#$ -P compbio_lab",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=30g",
                "#$ -l h_rt=15:00:00",
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "use R-3.5",
                sep="\n")


options(scipen = 999)


#library(RUVSeq)
library(umap)
library(ggplot2)
library(pheatmap)
#library(preprocessCore)
library(RNOmni)




file.in.histon.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.heightsMean.depthGCNormed.RData")
file.in.meta.RData= paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.goodSamples.meta.RData")
file.in.sampleInfo="~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"


dir.ou="./a_haQTLCalling_covariate_V1.2.3_pkFiltByMedianCV_gINT_Peer/"
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.HM.bed = paste0(dir.ou, HM, ".bed")
file.ou.HM.forPeer.bed = paste0(dir.ou, HM, ".forPeer.bed")
file.ou.sh = paste0("a.identifCov.peer.", HM, ".sh")
file.ou.covariates.prefix = paste0(dir.ou, HM, ".top", TOPCOV.NUM)
# file.ou.pdf = paste0(dir.ou, HM, ".covariates.pdf")
# file.ou.RData = paste0(dir.ou, HM, ".covariates.RData")


#
load(file.in.histon.RData)
load(file.in.meta.RData)

samp.info=read.table(file.in.sampleInfo, sep=",", header = T, row.names = 1, stringsAsFactors = F)
id2WGSid=samp.info$HA.Sample.Name
names(id2WGSid)= paste0("id_", rownames(samp.info))

#known factor
covar.known = meta.goodsamples
rownames(covar.known) = id2WGSid[rownames(covar.known)]
covar.known = covar.known[,]

# pks.filt.hg19 = merged.peaks.filted
# pks.hg382hg19.uniq = pks.hg382hg19.uniq[pks.hg382hg19.uniq %in% pks.filt.hg19]

Y.hg38 = samp2HeightsMean$depth.GCnormTeng[pks.hg382hg19.uniq, rownames(meta.goodsamples)]
rownames(Y.hg38) = names(pks.hg382hg19.uniq)

peaks.signalProp = apply(Y.hg38>=SIGNAL.CUTOFF, 1, sum)/ncol(Y.hg38)
peaks.hg38.filt = rownames(Y.hg38)[peaks.signalProp>=SIGNAL.PROP]
peaks.hg38.filt.median= apply(Y.hg38[peaks.hg38.filt, ], 1, median)
peaks.hg38.filtByMedian = peaks.hg38.filt[peaks.hg38.filt.median>=quantile(peaks.hg38.filt.median, 0.25)]
peaks.hg38.filtByMedian.CV= apply(Y.hg38[peaks.hg38.filtByMedian, ], 1, FUN=function(x) {sd(x)/mean(x)})
peaks.hg38.filtByMedianCV = peaks.hg38.filtByMedian[order(peaks.hg38.filtByMedian.CV, decreasing = T)][1:TOPPEAK.NUM]



Y.hg38.gINT = t(apply(Y.hg38[peaks.hg38.filt, ],1, rankNorm))

#
# sampID2PartID = as.character(meta$AliasCollaboratorParticipantID)
# names(sampID2PartID) = rownames(meta)
#
Y.hg38.gINT.partID = Y.hg38.gINT
colnames(Y.hg38.gINT.partID) = id2WGSid[colnames(Y.hg38)]
rownames(Y.hg38.gINT.partID) = peaks.hg38.filt


write(paste0("#chr\tstart\tend\tID\t", paste(colnames(Y.hg38.gINT.partID), collapse="\t")),
      file=file.ou.HM.bed)
peer.input = cbind(gsub(":|-", "\t", rownames(Y.hg38.gINT.partID)), rownames(Y.hg38.gINT.partID), Y.hg38.gINT.partID)
write.table(peer.input, file=file.ou.HM.bed, sep="\t", row.names=F, col.names = F, quote=F, append = T)

#
Y.hg38.gINT.partID.forPeer = Y.hg38.gINT.partID[peaks.hg38.filtByMedianCV,]

write(paste0("#chr\tstart\tend\tID\t", paste(colnames(Y.hg38.gINT.partID.forPeer), collapse="\t")),
      file=file.ou.HM.forPeer.bed)
peer.input = cbind(gsub(":|-", "\t", rownames(Y.hg38.gINT.partID.forPeer)), rownames(Y.hg38.gINT.partID.forPeer), Y.hg38.gINT.partID.forPeer)
write.table(peer.input, file=file.ou.HM.forPeer.bed, sep="\t", row.names=F, col.names = F, quote=F, append = T)


write(sh.head, file=file.ou.sh)
cmd = paste("Rscript", PEER.SCRIPT, file.ou.HM.forPeer.bed, file.ou.covariates.prefix, TOPCOV.NUM, sep=" ")
#system(cmd)
write(cmd, file = file.ou.sh, append = T)
write("echo work done", file = file.ou.sh, append = T)
system(paste0("qsub ", file.ou.sh))


