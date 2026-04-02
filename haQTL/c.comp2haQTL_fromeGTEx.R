#V1.1
#c.comp2haQTL_fromeGTEx_V1.1_sharedgARE.R
#check haQTL discovered from blood and test direction in brain

#compare with eGTEx Brain haQTL
#a1 and a2 are defined in vcf files which are reference allele and alternative allele when mapping QTLs

source("~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R")
library(ggplot2)
options(scipen = 999)

#CHRs = "chr19"

PEAK.WIN = 100000
haQTL.BP.CUTOFF="emp.p.fdr.cutoff0.2"

#
files.in.haQTL.eGTEx = system(paste0("find ~/lhou.compbio/codes/eGTEX-H3K27ac/haQTL/a_4_haQTL_FixedFactorNum_withGTExMinCov_V1.2.2_gINT_peer_rmAgeBatch_100k/Brain_nominal/chr*sorted.bed.gz"), intern=T)
names(files.in.haQTL.eGTEx)= sapply(files.in.haQTL.eGTEx, 
FUN=function(f.i)
{
  unlist(strsplit(basename(f.i), split=".", fixed=T))[1]
})
files.in.haQTL.BP = system(paste0("find ./a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_100k/H3K27ac_nominal/chr*sorted.bed.gz"), intern=T)
names(files.in.haQTL.BP )= sapply(files.in.haQTL.BP , 
FUN=function(f.i)
{
  unlist(strsplit(basename(f.i), split=".", fixed=T))[1]
})

# file.in.pk.hg38.eGTEx ="~/lhou.compbio/codes/eGTEX-H3K27ac/peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/Brain.pks.hg382hg19Uniq.RDS"
# file.in.pk.hg38.BP="../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/H3K27ac.MayoWithRef.heightsMean.depthGCNormed.RData"


#file.in.pk.hg38.eGTEx.RData="~/lhou.compbio/codes/eGTEX-H3K27ac/haQTL/a_5_haQTL_multTest_FixedFactorNum_withGTExMinCov_V1.2.2_gINT_peer_rmAgeBatch_100k/Brain.haQTLPeak.fdr0.05.hg38.RData"

file.in.pk.hg38.eGTEx.bed="~/lhou.compbio/codes/eGTEX-H3K27ac/haQTL/a_5_haQTL_multTest_FixedFactorNum_withGTExMinCov_V1.2.2_gINT_peer_rmAgeBatch_100k/Brain.bgPeak.hg38.bed"
file.in.pk.hg38.BP.RData="../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/H3K27ac.haQTLPeak.cutoffs.hg38.RData"

file.in.mARE.RData="../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/3activeHMs.mergedPeak2HMPeak.hg19.RData"
file.in.mAREgrp.RData="../peakVariationAcrossTiss/a_4_modules_annotations_mergedPk3actHMs_Epimap/mARE2ModuleGrpName.RData"
file.in.hg38tohg19.RDS="../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/H3K27ac.pk.hg38tohg19.RDS"

dir.tmp = "~/hptmp/MayoBipolar/a_cmpHaQTLvseGTEx_V1.1/"
dir.create(dir.tmp, showWarnings = F, recursive = T)  

#file.tmp.peak.eGTEx.bed = paste(dir.tmp, "H3K27ac.eGTEx.peak.bed", sep="")
file.tmp.peak.BP.bed = paste(dir.tmp, "H3K27ac.BP.peak.bed", sep="")
file.tmp.peakOvlp.bed = paste(dir.tmp,"H3K27ac.peakOvlp.bed", sep="")



dir.ou= "c_compHaQTLWitheGTExBrain_V1.1_sharedgAREs/"
dir.create(dir.ou, showWarnings = F, recursive = T)  

file.ou.pdf = paste(dir.ou, "H3K27ac.vseGTEx_Brain_", PEAK.WIN/1000, "k_BPhaQTL.", haQTL.BP.CUTOFF, ".pdf", sep="")
file.ou.RData = paste(dir.ou, "H3K27ac.vseGTEx_Brain_", PEAK.WIN/1000, "k_BPhaQTL", haQTL.BP.CUTOFF, ".RData", sep="")

#

#generate overlap peaks
#pks.hg38tohg19.eGTEx = readRDS(file.in.pk.hg38.eGTEx)
#load(file.in.pk.hg38.BP)


load(file.in.pk.hg38.BP.RData)
pk.hg38.filtByEmpP.BP = hQTLPeaks.list[[haQTL.BP.CUTOFF]]
write.table(cbind(gsub(":|-", "\t", pk.hg38.filtByEmpP.BP), pk.hg38.filtByEmpP.BP), sep="\t", file=file.tmp.peak.BP.bed, quote=F, row.names = F, col.names = F)

# 
# load(file.in.pk.hg38.eGTEx.RData)
# pk.hg38.filtByEmpP.eGTEx = rownames(emp.pvals.filtNA)[emp.pvals.filtNA[,"emp.P"]<=EMP.P.CUTOFF]
# 

#write.table(cbind(gsub(":|-", "\t", pk.hg38.filtByEmpP.eGTEx), pk.hg38.filtByEmpP.eGTEx), sep="\t", file=file.tmp.peak.eGTEx.bed, quote=F, row.names = F, col.names = F)



system(paste0("intersectBed -a ", file.in.pk.hg38.eGTEx.bed, " -b ", file.tmp.peak.BP.bed, " -wa -wb  > ", file.tmp.peakOvlp.bed))

peaks.ovlp.eGTExAndBP = read.table(file.tmp.peakOvlp.bed, sep="\t", header = F, row.names = NULL, stringsAsFactors = F)
colnames(peaks.ovlp.eGTExAndBP)=c("eGTEx.chr", "eGTEx.start", "eGTEx.end", "BP.chr", "BP.start", "BP.end", "BP.pk")
peaks.ovlp.eGTExAndBP$eGTEx.pk = paste0(peaks.ovlp.eGTExAndBP[,1], ":", peaks.ovlp.eGTExAndBP[,2], "-", peaks.ovlp.eGTExAndBP[,3])
#peaks.ovlp.BP2Deconv = split(peaks.ovlp.eGTExAndBP[,2], peaks.ovlp.eGTExAndBP[,1])
#peaks.ovlp.BP = levels(factor(peaks.ovlp.eGTExAndBP[,1]))
#chr2peaks.ovlp.BP = split(peaks.ovlp.BP, sapply(peaks.ovlp.BP, FUN=function(x){unlist(strsplit(x, "\t"))[1]}))


#generate haQTL for overlap peaks
# peaks.BP.isSlct = which(peaks$BP %in% peaks.uniq$BP[peaks.ovlp.BP])
# peaks.BP.isSlct = paste(peaks.BP.isSlct, collapse="p;")
# peaks.BP.isSlct = paste(peaks.BP.isSlct, "q", sep="")
# buf = system(paste("zcat ", files.in.haQTL.blueprint[CELLTYPE.SLCT], " |sed -ne '", peaks.BP.isSlct, "'", sep=""), intern=T)

topBPhaQTL.ovlpPeaks=lapply(1:nrow(peaks.ovlp.eGTExAndBP),
FUN=function(i)
{
  if(i %% 1000==0) print(i)
  pk.eGTEx = peaks.ovlp.eGTExAndBP$eGTEx.pk[i]
  pk.BP = peaks.ovlp.eGTExAndBP$BP.pk[i]
  
  #buf = unlist(strsplit(pk.eGTEx, split=":|-", perl=T))
  pk.chr= peaks.ovlp.eGTExAndBP$BP.chr[i]
  if(!(pk.chr %in% names(files.in.haQTL.BP)))
  {
    return(data.frame(snpid = NA, BP.pk = NA, BP.beta = NA, BP.pval= NA, eGTEx.pk = NA, eGTEx.beta= NA, eGTEx.pval= NA))
  }
  pk.start.BP = peaks.ovlp.eGTExAndBP$BP.start[i]
  pk.end.BP = peaks.ovlp.eGTExAndBP$BP.end[i]
  
  pk.start.eGTEx = peaks.ovlp.eGTExAndBP$eGTEx.start[i]
  pk.end.eGTEx= peaks.ovlp.eGTExAndBP$eGTEx.end[i]
  
  #
  pk.eGTEx.lb = max(1, pk.start.eGTEx-PEAK.WIN)
  pk.eGTEx.ub = pk.end.eGTEx+PEAK.WIN
  query = pk.chr %&&% ':' %&&% pk.eGTEx.lb %&&% '-' %&&% pk.eGTEx.ub
  haQTL.eGTEx = seqminer::tabix.read.table(files.in.haQTL.eGTEx[pk.chr], query) 
  if(length(haQTL.eGTEx)==0)
  {
    return(data.frame(lead.snpid = NA, BP.pk = NA, BP.beta = NA, BP.pval= NA, eGTEx.pk = NA, eGTEx.beta= NA, eGTEx.pval= NA))
  }
  haQTL.eGTEx.filt = haQTL.eGTEx[haQTL.eGTEx$peak==pk.eGTEx,]
  rownames(haQTL.eGTEx.filt) = paste0(haQTL.eGTEx.filt$chr, "_", haQTL.eGTEx.filt$end, "_", haQTL.eGTEx.filt$a1, "_", haQTL.eGTEx.filt$a2)
  #
  haQTL.eGTEx.filt.flip = haQTL.eGTEx.filt
  haQTL.eGTEx.filt.flip$beta= -haQTL.eGTEx.filt$beta
  haQTL.eGTEx.filt.flip$a1= haQTL.eGTEx.filt$a2
  haQTL.eGTEx.filt.flip$a2= haQTL.eGTEx.filt$a1
  rownames(haQTL.eGTEx.filt.flip) = paste0(haQTL.eGTEx.filt.flip$chr, "_", haQTL.eGTEx.filt.flip$end, "_", haQTL.eGTEx.filt.flip$a1, "_", haQTL.eGTEx.filt.flip$a2)
  #
  haQTL.eGTEx.all = rbind(haQTL.eGTEx.filt, haQTL.eGTEx.filt.flip)
  
  #
  pk.BP.lb = max(1, pk.start.BP-PEAK.WIN)
  pk.BP.ub = pk.end.BP+PEAK.WIN
  query = pk.chr %&&% ':' %&&% pk.BP.lb %&&% '-' %&&% pk.BP.ub
  haQTL.BP = seqminer::tabix.read.table(files.in.haQTL.BP[pk.chr], query) 
  if(length(haQTL.BP)==0)
  {
    return(data.frame(lead.snpid = NA, BP.pk = NA, BP.beta = NA, BP.pval= NA, eGTEx.pk = NA, eGTEx.beta= NA, eGTEx.pval= NA))
  }
  haQTL.BP.filt = haQTL.BP[haQTL.BP$peak==pk.BP,]
  rownames(haQTL.BP.filt) = paste0(haQTL.BP.filt$chr, "_", haQTL.BP.filt$end, "_", haQTL.BP.filt$a1, "_", haQTL.BP.filt$a2)
  
  
  snps.loc.ovlp = intersect(rownames(haQTL.eGTEx.all), rownames(haQTL.BP.filt))
  
  haQTL.cmb = lapply(snps.loc.ovlp,
  FUN=function(snp.loc)
  {
    return(data.frame(lead.snpid = snp.loc, 
                      BP.pk = haQTL.BP.filt[snp.loc, "peak"], BP.beta = haQTL.BP.filt[snp.loc, "beta"], BP.pval= haQTL.BP.filt[snp.loc, "pval"],
                      eGTEx.pk = haQTL.eGTEx.all[snp.loc, "peak"], eGTEx.beta= haQTL.eGTEx.all[snp.loc, "beta"], eGTEx.pval= haQTL.eGTEx.all[snp.loc, "pval"]))
 
    
  })
  #
  haQTL.cmb.df = do.call(rbind, haQTL.cmb)
  
  
  #haQTL.cmb.df = haQTL.cmb.df[!is.na(haQTL.cmb.df$snpid), ]
  
  
  return(haQTL.cmb.df[which.min(haQTL.cmb.df$BP.pval)[1], ])
})
save(topBPhaQTL.ovlpPeaks, 
     file=file.ou.RData)


topBPhaQTL.ovlpPeaks.df =do.call(rbind, topBPhaQTL.ovlpPeaks)
topBPhaQTL.ovlpPeaks.filt.df = topBPhaQTL.ovlpPeaks.df[!is.na(topBPhaQTL.ovlpPeaks.df[,1]),]

topBPhaQTL.ovlpPeaks.filt.df$log10BP.pval = -log10(topBPhaQTL.ovlpPeaks.filt.df$BP.pval)
topBPhaQTL.ovlpPeaks.filt.df$log10eGTEx.pval = -log10(topBPhaQTL.ovlpPeaks.filt.df$eGTEx.pval)
topBPhaQTL.ovlpPeaks.filt.df$log10eGTEx.pval.bin= floor(topBPhaQTL.ovlpPeaks.filt.df$log10eGTEx.pval)
topBPhaQTL.ovlpPeaks.filt.df$log10eGTEx.pval.bin[topBPhaQTL.ovlpPeaks.filt.df$log10eGTEx.pval.bin>=5]=">=5"
topBPhaQTL.ovlpPeaks.filt.df$BPSigned.log10eGTEx.pval = topBPhaQTL.ovlpPeaks.filt.df$log10eGTEx.pval * sign(topBPhaQTL.ovlpPeaks.filt.df$BP.beta)
#
#topBPhaQTL.ovlpPeaks.filt.df$BPSigned.log10BP.pval = topBPhaQTL.ovlpPeaks.filt.df$log10BP.pval * sign(topBPhaQTL.ovlpPeaks.filt.df$BP.beta)
# topBPhaQTL.ovlpPeaks.filt.df$log10BP.pval.bin= floor(topBPhaQTL.ovlpPeaks.filt.df$log10BP.pval)
# topBPhaQTL.ovlpPeaks.filt.df$log10BP.pval.bin[topBPhaQTL.ovlpPeaks.filt.df$log10BP.pval.bin>=5]=">=5"


save(topBPhaQTL.ovlpPeaks,
     topBPhaQTL.ovlpPeaks.filt.df,
     file=file.ou.RData)
#load(file.ou.RData)


#directionality consistency
#
load(file.in.mAREgrp.RData)
load(file.in.mARE.RData)
pk.hg38tohg19=readRDS(file.in.hg38tohg19.RDS)

topBPhaQTL.ovlpPeaks.filt.df$mPk.hg19=hm.peak2mergedPeak[["H3K27ac"]][pk.hg38tohg19[topBPhaQTL.ovlpPeaks.filt.df$BP.pk]]
topBPhaQTL.ovlpPeaks.filt.df$mPk.hg19=sub("\t", ":", topBPhaQTL.ovlpPeaks.filt.df$mPk.hg19)
topBPhaQTL.ovlpPeaks.filt.df$mPk.hg19=sub("\t", "-", topBPhaQTL.ovlpPeaks.filt.df$mPk.hg19)
topBPhaQTL.ovlpPeaks.filt.df$mPk.AREGrp=mARE2ModGrpName[topBPhaQTL.ovlpPeaks.filt.df$mPk.hg19]
eGTEx.pval.cutoffs=c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001)
is.signConsist = sign(topBPhaQTL.ovlpPeaks.filt.df$BP.beta) == sign(topBPhaQTL.ovlpPeaks.filt.df$eGTEx.beta)
is.AREGrpBrain =  topBPhaQTL.ovlpPeaks.filt.df$mPk.AREGrp=="Brain/neuron"
is.Ubiq =  topBPhaQTL.ovlpPeaks.filt.df$mPk.AREGrp=="Ubiquitious"
is.Mult=  topBPhaQTL.ovlpPeaks.filt.df$mPk.AREGrp=="Multi-tissue"
directConsist.perc = t(sapply(eGTEx.pval.cutoffs,
                              FUN=function(p.cutoff)
                              {
                                is.pass= topBPhaQTL.ovlpPeaks.filt.df$eGTEx.pval<=p.cutoff
                                
                                N=sum(is.pass)
                                N.consist = sum(is.signConsist[is.pass])
                                N.Brain = sum(is.AREGrpBrain[is.pass], na.rm=T)
                                N.Ubiq = sum(is.Ubiq[is.pass], na.rm=T)
                                N.Mult = sum(is.Mult[is.pass], na.rm=T)
                                
                                return(c(N=N, 
                                         directConsist.perc=N.consist/N * 100,
                                         BrainAREGrp.perc = N.Brain/N*100,
                                         UbiqAREGrp.perc = N.Ubiq/N*100,
                                         MultAREGrp.perc = N.Mult/N*100))
                              }))
rownames(directConsist.perc) = eGTEx.pval.cutoffs

save(topBPhaQTL.ovlpPeaks,
     topBPhaQTL.ovlpPeaks.filt.df,
     eGTEx.pval.cutoffs,
     directConsist.perc,
     file=file.ou.RData)



#
f.o.tsv = paste0(dir.ou, "H3K27ac.vseGTEx_Brain_", PEAK.WIN/1000, "k_BPhaQTL", haQTL.BP.CUTOFF, ".tsv")
write.table(topBPhaQTL.ovlpPeaks.filt.df, file=f.o.tsv, sep="\t", quote=F, row.names=F)


  
pdf(file.ou.pdf, height=8, width=15)

#
p =ggplot(topBPhaQTL.ovlpPeaks.filt.df, aes(x=BPSigned.log10eGTEx.pval, y=eGTEx.beta)) +
  geom_point(aes(color=log10eGTEx.pval.bin), size=0.5) +
  scale_colour_manual(values = c("0"="#7AC5CD",
                                 "1"="#FFC0CB",
                                 "2"="#EEA9B8",
                                 "3"="#CD919E",
                                 "4"="#8B636C",
                                 ">=5"="#8B0000")) +
  theme_bw()
  
p

#
p =ggplot(topBPhaQTL.ovlpPeaks.filt.df, aes(x=BP.beta, y=eGTEx.beta)) +
  geom_point(aes(color=log10eGTEx.pval.bin), size=0.5) +
  scale_colour_manual(values = c("0"="#7AC5CD",
                                 "1"="#FFC0CB",
                                 "2"="#EEA9B8",
                                 "3"="#CD919E",
                                 "4"="#8B636C",
                                 ">=5"="#8B0000")) +
  theme_bw()

p




#
layout(matrix(1:4, ncol=2))
par(mar=c(5, 8, 5, 8))
plot(-log10(eGTEx.pval.cutoffs),
     directConsist.perc[, "directConsist.perc"],
     xlab="-log10 eGTEx pval as cutoff",
     ylab="directinality consistency",
     col="red",
     type="b")
par(new=T)
plot(-log10(eGTEx.pval.cutoffs), 
     directConsist.perc[, "N"], 
     col="blue",
     #type="b",
     xaxt="n", yaxt="n",xlab="",ylab="", type="b", bty="n", 
     lwd=2, cex= 1.5)
axis(4)
mtext("number of CREs filtered", side = 4, line = 4, col="blue")

plot(-log10(eGTEx.pval.cutoffs),
     directConsist.perc[, "BrainAREGrp.perc"],
     xlab="-log10 eGTEx pval as cutoff",
     ylab="Brain ARE group percentage",
     col="red",
     type="b")


plot(-log10(eGTEx.pval.cutoffs),
     directConsist.perc[, "UbiqAREGrp.perc"],
     xlab="-log10 eGTEx pval as cutoff",
     ylab="Ubiq ARE group percentage",
     col="red",
     type="b")

plot(-log10(eGTEx.pval.cutoffs),
     directConsist.perc[, "MultAREGrp.perc"],
     xlab="-log10 eGTEx pval as cutoff",
     ylab="Multi-tissue ARE group percentage",
     col="red",
     type="b")


dev.off()
  




