#V4.2
#a.mergePeak2ReadsCounts_V4.2_filteredMergedPeak_H3K27ac.R
#keep reads overlapping with all filtered merged peaks 
#compare the proportion of reads in reference overlapped peaks, selected reference overlapping peaks

#V4.1 
#merge peak only from samples with enough reads and good RSC
#use logq as cutoff to filter peaks

#V4
#merge peaks from samples with good quality [defined by mismatch] and only non-AD related samples would be included
#do not use wiggler to calculate peak heights, but use shifted reads coverage instead
#also caculate counts for referecen samples
#normalized to 1kb region are not used any more

#V3
#merge peaks from samples only with good quality and non-AD samples would be included
#remove those merged peaks overlapping with peak from at least two individuals
#get peak counts from normalized bw, and re-normalzied to 30,000,000 reads in total from 1*10^9


#V2
#merge peaks, get reads counts and normalize by GC content
#merge peak
# need to have overlap with peak from at least one other individuals


#for GC  normalization
#a) base-pair-level union of peak regions from all individuals as the set of active regions in our samples
#b) H3K27ac peak height was defined as the number of ChIP-seq fragments in the peak regions. For paired-end reads, the ChIP-seq fragment was defined as the entire region between the outer edges of the read pair, whereas for single-end reads we used the median fragment length across all libraries. Fragments that lay only partially within the peak region were counted according to their overlap fraction. 
#c) Peak heights were normalized for G+C content biases as in Pickrell et al
#https://www.nature.com/nature/journal/v464/n7289/full/nature08872.html


args = commandArgs(trailingOnly=TRUE)
HM = args[1]





R.LIBS.MYPATHs = c(V3.5="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.5/",
                V3.4="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.4/",
                V3.3="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.3/")
R.LIBS.PATH = .libPaths()

if(grepl("R version 3.5", R.version$version.string))
{
  .libPaths(c(R.LIBS.MYPATHs["V3.5"], R.LIBS.PATH))
}
if(grepl("R version 3.4", R.version$version.string))
{
  .libPaths(c(R.LIBS.MYPATHs["V3.4"], R.LIBS.PATH))
}
if(grepl("R version 3.3", R.version$version.string))
{
  .libPaths(c(R.LIBS.MYPATHs["V3.3"], R.LIBS.PATH))
}



N.CUTOFF=2 #need to appear in at least 10 percent of all samples
SUPER.TOTALREADs.CUTOFF=3e7
SUPER.RSC.CUTOFF=2
PEAK.LogQ.CUTOFF =2
GOOD.TOTALREADs.CUTOFF=2e7
GOOD.RSC.CUTOFF=1

TYPE2COL=c("grey","red")
names(TYPE2COL)=c("control", "case") #c("IBDRA_Control", "CD", "UC", "RA")

REFs.SELECT= c("E040","E039", "E048", "E047", "E029", "E032", "E046", "BPNeutrop")

REF2CELLTYPE= 
  c(
    E062="Primary mononuclear cells from peripheral blood",
    E034="Primary T cells from peripheral blood",
    E045="Primary T cells effector/memory enriched from peripheral blood",
    #E033="Primary T cells from cord blood",
    E044="Primary T regulatory cells",# from peripheral blood",
    E043="Primary T helper cells",# from peripheral blood",
    E041="Primary T helper cells PMA-I stimulated",
    E042="Primary T helper 17 cells PMA-I stimulated",
    E040="Primary T helper memory cells1",#  from peripheral blood 1",
    E037="Primary T helper memory cells2",# from peripheral blood 2",
    E038="Primary T helper naive cells",#  from peripheral blood",
    E039="Primary T helper naive cells",# from peripheral blood",
    E048="Primary T CD8+ memory cells",#  from peripheral blood",
    E047="Primary T CD8+ naive cells",#  from peripheral blood",   
    E029="Primary monocytes",#  from peripheral blood",
    #E031="Primary B cells from cord blood",
    # E035="Primary hematopoietic stem cells",
    # E051="Primary hematopoietic stem cells G-CSF-mobilized Male",
    #E050="Primary hematopoietic stem cells G-CSF-mobilized Female",
    # E036="Primary hematopoietic stem cells short term culture",
    E032="Primary B cells",#  from peripheral blood",
    E046="Primary Natural Killer cells"#,#  from peripheral blood",
    #E030="Primary neutrophils"#  from peripheral blood"
    )

sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                "#$ -P compbio_lab",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=30g",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                sep="\n")


##Mayo
dir.in = paste0("~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/AQUAS_", HM, "/")
file.in.meta = paste0("~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/c_metaAndQC/", HM, ".sample.metaDataAndQC.txt")
file.in.peak.suf = "out/peak/macs2/rep1/*nodup.tagAlign.pval0.01.500K.filt.narrowPeak.gz"
file.in.tagAlign.suf = "out/align/rep1/*nodup.tagAlign.gz"
# file.in.badSamp = "b_tissueSpecif/sample.QC.RSCcutoff1.txt"
# file.in.qc = "b_QCmetric/QC_metrics.txt"
# file.in.testMeta="~/lhou.compbio/data/Mayo_aim2_5.2017/H3K27ac_3.27.2018/sample.metaData.until3.27.2018.txt"

#roadmap
file.in.ref.tag=paste0("/broad/compbio/anshul/projects/roadmap/alignments/consolidated_nosubsampling/E000-", HM, ".tagAlign.gz")
file.in.ref.pk=paste0("/broad/compbio/anshul/projects/roadmap/peaks/consolidated/narrowPeak/E000-", HM, ".narrowPeak.gz")
file.in.ref.fragLen = "/broad/compbio/anshul/projects/roadmap/metadata/strandcorr/all.encode.roadmap.fraglen"

#blueprint
dir.blueprint = "/broad/compbio/data/ihec/Venous_blood_2018/Venous_blood/"
blueprint.indNames = c("S00D9Y", "S00DKC", "S00DP2", "S00E9U", "S00F9Q", "S00J9A", "S00P04", "S00R1V", "S0100M", "S01342")
blueprint.nm2celltypes = c("mature_neutrophil", "CD4-positive_alpha-beta_T_cell", "CD14-positive_CD16-negative_classical_monocyte")
names(blueprint.nm2celltypes) = c("BP-Neutrophil", "BP-CD4NaiveT", "BP-monocyte")
files.in.blueprint = list()
for (ind in blueprint.indNames)
{
  for(nm in names(blueprint.nm2celltypes))
  {
    ct = blueprint.nm2celltypes[nm]
    
    f.bw = system(paste("find ", dir.blueprint, ind, "/", ct, "/*/*/*", HM, "*.bw", sep=""), intern=T)
    f.pk = system(paste("find ", dir.blueprint, ind, "/", ct, "/*/*/*", HM, "*.bed.gz", sep=""), intern=T)
    
    if(length(f.bw)==1 & length(f.pk)==1)
    {
      files.in.blueprint[[paste(nm, ".", ind, sep="")]]=c(bw=f.bw, pk=f.pk)
    }else
    {
      print(paste("NO", ind, nm))
    }
  }
  
}


dir.tmp = "~/hptmp/MayoBipolar/peakMerge/"
if(! dir.exists(dir.tmp))
{
  dir.create(dir.tmp)
}
file.tmp.allPeak = paste(dir.tmp, HM, ".hg19.nodup.filt.narrowPeak.all.bed", sep="")
file.tmp.mergedPeak = paste(dir.tmp, HM, ".hg19.nodup.filt.narrowPeak.merged.bed", sep="")
file.tmp.filt.mergedPeak = paste(dir.tmp, HM, ".hg19.nodup.filt.narrowPeak.merged.filt", N.CUTOFF, "Samps.bed", sep="")
file.tmp.pkInRef = paste(dir.tmp, HM, ".merged.peaks.inRef.bed", sep="")


dir.ou="./a_mergePeak2readsCounts_V4.2_filteredMergedPeaks/"
if(! dir.exists(dir.ou))
{
  dir.create(dir.ou)
}
file.ou.filt.mergedPeak = paste(dir.ou, HM, ".hg19.narrowPeak.RSC", SUPER.RSC.CUTOFF, ".Tot", SUPER.TOTALREADs.CUTOFF, ".logQ", PEAK.LogQ.CUTOFF, ".merged.filt", N.CUTOFF, "GoodSamps.bed", sep="")
file.ou.RData=paste(dir.ou, HM, ".hg19.narrowPeak.RSC", SUPER.RSC.CUTOFF, ".Tot", SUPER.TOTALREADs.CUTOFF, ".logQ", PEAK.LogQ.CUTOFF, ".merged.filt", N.CUTOFF, "GoodSamps.heightsMean.RData", sep="")
#file.ou.fractMatrix = paste(dir.ou, "hg19.nodup.filt.narrowPeak.merged.No0Fraction", sep="")

file.ou.peakMerged.pdf = paste(dir.ou, HM, ".diagnose.peakMerge.pdf", sep="")
file.ou.ovlp.pdf = paste(dir.ou, HM, ".diagnose.ovlpInRef.pdf", sep="")
file.ou.pdf =   paste(dir.ou, HM, ".hg19.narrowPeak.merged.filt", N.CUTOFF, "GoodSamps.heigthsMean.pdf", sep="")
file.ou.readPerct.pdf =   paste(dir.ou, HM, ".hg19.narrowPeak.merged.filt", N.CUTOFF, "GoodSamps.readsPercent.pdf", sep="")

#
print("input meta info") ########################################################################
# qc.info=c() #sample names, peak file path, bw file path, whether good or not based RSC and comparison to ROADMAP, 
# 
# for(d.i in dirs.in)
# {
#   isBad = as.matrix(read.table(paste(d.i, file.in.badSamp, sep=""), sep="\t", header=T, row.names=1))[,"isBad.byMismatch"]
#   fragLen = as.matrix(read.table(paste(d.i, file.in.qc, sep=""), sep="\t", header=T, row.names=1))[,"FragLen"]
#   for(samp in names(isBad))
#   {
#     f.pk = system(paste("find ", d.i, "AQUAS/", samp, "/", file.in.peak.suf, sep=""), intern = T)
#     f.tag = system(paste("find ", d.i, "AQUAS/", samp, "/", file.in.tagAlign.suf, sep=""), intern = T)
#     qc.info=rbind(qc.info, c(samp, isBad[samp]==0, fragLen[samp], f.pk, f.tag))
#   }
#   
# }
# rownames(qc.info)=qc.info[,1]
# qc.info=qc.info[,2:ncol(qc.info)]
# colnames(qc.info) = c("GoodQuality", "fragLen", "Peak", "Tag")


qc.info = read.table(file.in.meta, sep="\t", row.names=1, header = T)


#
#samp2meta = read.table(file.in.testMeta, sep="\t", header=T, row.names=1)#[colnames(samp2Heights.withRef$testWithRef.GCnorm.quant)[1:60],]
# dis2samps = lapply(names(TYPE2COL),
# FUN=function(dis)
# {
#   rownames(samp2meta)[samp2meta$phenotype==dis]
# })
# names(dis2samps)=names(TYPE2COL)
#samps.nonAD = unlist(dis2samps[c("IBDRA_Control", "CD", "UC", "RA", "AD_control")])

#########################################################################
print("merged peak") 
if(file.exists(file.tmp.allPeak))
{
  file.remove(file.tmp.allPeak)
}

for(samp.id in rownames(qc.info)[qc.info$totalReads>=SUPER.TOTALREADs.CUTOFF & qc.info$RSC>=SUPER.RSC.CUTOFF & qc.info$isBad.byMismatch==0])
{
  f.i.pk = paste0(dir.in, samp.id, "/", file.in.peak.suf)
  print(f.i.pk)
  #cmd = paste("zcat ", f.i.pk, " >>", file.tmp.allPeak)
  cmd = paste("zcat ", f.i.pk, "|awk '($9>=", PEAK.LogQ.CUTOFF, "){print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $4 \"_\" $8}' >>", file.tmp.allPeak, sep="")
  system(cmd)
}

cmd = paste("sort -k1,1 -k2,2n ", file.tmp.allPeak, " | mergeBed -i stdin -c 4 -o collapse -delim \",\" >", file.tmp.mergedPeak, sep="")
system(cmd)


#diagnosis
pk2len = as.numeric(system(paste("awk '{print $3-$2}'", file.tmp.mergedPeak), intern = T)) #- as.numeric(system(paste("awk '{print $2}'", file.tmp.mergedPeak), intern = T))
buf=system(paste("awk '{print $4}'", file.tmp.mergedPeak), intern = T)
pk2pvals = lapply(buf,
FUN=function(pks)
{
  sapply(unlist(strsplit(pks, split=",")),
  FUN=function(pk)
  {
    as.numeric(unlist(strsplit(pk, split="_"))[3])
  })
})
pk2pvals.size = sapply(pk2pvals,length)
size.seq = sort(as.numeric(levels(factor(pk2pvals.size))))
size2pval = lapply(size.seq,
FUN=function(size)
{
  unlist(pk2pvals[pk2pvals.size==size])
})
names(size2pval) = as.character(size.seq)

size2SecondTopPval = lapply(size.seq,
FUN=function(size)
{
  if(size==1)
  {
    return(unlist(pk2pvals[pk2pvals.size==1]))
  }
  sapply(pk2pvals[pk2pvals.size==size],
  FUN=function(qs)
  {
    
    sort(qs, decreasing = T)[2]
  })
  
})
pdf(file.ou.peakMerged.pdf, height=9, width=12)
layout(1:4)
hist(pk2pvals.size, xlab="number of peaks overlapped with merged peaks", breaks=100)
boxplot(size2pval, xlab="number of peaks overlapped with merged peaks", ylab="pvalues")
boxplot(size2SecondTopPval, xlab="number of peaks overlapped with merged peaks", ylab="2nd largest pvalues")
boxplot(pk2len~pk2pvals.size, xlab="number of peaks overlapped with merged peaks", ylab="peak length")
layout(matrix(1:9, ncol=3))
for(i in  1:sum(qc.info$totalReads>=SUPER.TOTALREADs.CUTOFF & qc.info$RSC>=SUPER.RSC.CUTOFF & qc.info$isBad.byMismatch==0))
{
  print(i)
  hist(size2pval[[i]], main=paste("size ", i, sep=""), xlab="qvalue", breaks=200)
}
dev.off()

 


cmd = paste("awk '{if (gsub(/,/,\"\",$4)>=", N.CUTOFF-1,") print}' ", file.tmp.mergedPeak," > ", file.tmp.filt.mergedPeak, sep="")
system(cmd)

cmd= paste("awk '{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $1 \":\" $2 \"-\" $3}' ", file.tmp.filt.mergedPeak, ">", file.ou.filt.mergedPeak, sep="" )
system(cmd)

# cmd =paste("awk '{print $1 \":\" $2 \"-\" $3}' ", file.ou.filt.mergedPeak, sep="")
# merged.peaks = system(cmd, intern = T)
merged.peaks.df = read.table(file.ou.filt.mergedPeak, header=F, sep="\t", row.names = 4)
colnames(merged.peaks.df) = c("chr", "start", "end")
merged.peaks = rownames(merged.peaks.df)
merged.peaks.len = merged.peaks.df[,3]-merged.peaks.df[,2]
names(merged.peaks.len) = merged.peaks
# 

##################################################
print("get reads overlapping with filtered merged peaks for both testing samples and reference samples")
####################################################

print("mean counts in each peak for Mayo sample")
#ref region  names
 

for(samp in rownames(qc.info)[qc.info$totalReads>=GOOD.TOTALREADs.CUTOFF & qc.info$RSC>=GOOD.RSC.CUTOFF & qc.info$isBad.byMismatch==0])
{
  f.i.tag=paste(dir.in, samp, "/", file.in.tagAlign.suf, sep="")
  frag.len = qc.info[samp,"FragLen"]
  f.o.sh = paste("a.", samp, ".ovlp.merged.sh", sep="")
  f.tmp.merged = paste(dir.tmp, samp, ".", HM, ".ovlp.merged.bed", sep="")
  
  # if(unlist(strsplit(system(paste("wc -l ", f.tmp.merged, sep=""), intern = T), split=" "))[1] != 0)
  # {
  #   next
  # }
  write(sh.head, f.o.sh)
  
  cmd= paste("zcat ", f.i.tag, "|awk '{if($6==\"+\") {print $1 \"\\t\" $2 \"\\t\" $2+", frag.len, "} if($6==\"-\") {start =1>$3-", frag.len, "?1:$3-", frag.len, ";print $1 \"\\t\" start \"\\t\" $3 }}'|sortBed -i stdin| intersectBed -b stdin -a ", file.ou.filt.mergedPeak, " -wo -sorted|mergeBed -i stdin -c 4,8 -o distinct,sum >", f.tmp.merged, sep="")
  write(cmd, file=f.o.sh, append = T)
  
  write("echo done",file=f.o.sh, append = T)
  system(paste("qsub ", f.o.sh, sep=""))
}

#need to wait for the above codes running

samp2HeightsMean=list()

samp2HeightsMean$Mayo = sapply(rownames(qc.info)[qc.info$totalReads>=GOOD.TOTALREADs.CUTOFF & qc.info$RSC>=GOOD.RSC.CUTOFF & qc.info$isBad.byMismatch==0],
FUN=function(samp)
{
  print(samp)
  f.tmp.merged = paste(dir.tmp, samp, ".", HM, ".ovlp.merged.bed", sep="")

  heights  =as.matrix(read.table(f.tmp.merged, sep="\t", header=F, row.names=4, colClasses = "character"))[,4]
  pk2heightMean = rep(0, length(merged.peaks))
  names(pk2heightMean) = merged.peaks
  pk2heightMean[names(heights)] = as.numeric(heights)/merged.peaks.len[names(heights)]
  
  return(pk2heightMean)
})
colnames(samp2HeightsMean$Mayo) = rownames(qc.info)[qc.info$totalReads>=GOOD.TOTALREADs.CUTOFF & qc.info$RSC>=GOOD.RSC.CUTOFF & qc.info$isBad.byMismatch==0]

save(samp2HeightsMean, merged.peaks, merged.peaks.len, qc.info,file=file.ou.RData)

#rownames(samp2Heights$Mayo) = merged.peaks.inRef
#rownames(samp2Heigths) =merge.peaks

#################################################
print("get counts for ROADMAP samples")

file2fragLen = as.matrix(read.table(file.in.ref.fragLen, sep="\t", row.names=1, header=F))

for(samp in names(REF2CELLTYPE))
{
  f.i.tag= sub("E000", samp, file.in.ref.tag)
  frag.len = file2fragLen[basename(f.i.tag),1]
  f.o.sh = paste("a.", samp, ".ovlp.merged.sh", sep="")
  f.tmp.merged = paste(dir.tmp, samp, ".", HM, ".ovlp.merged.bed", sep="")
  
  # if(unlist(strsplit(system(paste("wc -l ", f.tmp.merged, sep=""), intern = T), split=" "))[1] != 0)
  # {
  #   next
  # }
  
  write(sh.head, f.o.sh)
  
  cmd= paste("zcat ", f.i.tag, "|awk '{if($6==\"+\") {print $1 \"\\t\" $2 \"\\t\" $2+", frag.len, "} if($6==\"-\") {start =1>$3-", frag.len, "?1:$3-", frag.len, ";print $1 \"\\t\" start \"\\t\" $3 }}'|sortBed -i stdin| intersectBed -b stdin -a ", file.ou.filt.mergedPeak, " -wo -sorted|mergeBed -i stdin -c 4,8 -o distinct,sum >", f.tmp.merged, sep="")
  write(cmd, file=f.o.sh, append = T)
  
  write("echo done",file=f.o.sh, append = T)
  system(paste("qsub ", f.o.sh, sep=""))
}

#need to wait for the above codes running
samp2HeightsMean$ROADMAP = sapply(names(REF2CELLTYPE),
FUN=function(samp)
{
  print(samp)
  f.tmp.merged =  paste(dir.tmp, samp, ".", HM, ".ovlp.merged.bed", sep="")
  
  heights  =as.matrix(read.table(f.tmp.merged, sep="\t", header=F, row.names=4, colClasses = "character"))[,4]
  pk2heightMean = rep(0, length(merged.peaks))
  names(pk2heightMean) = merged.peaks
  pk2heightMean[names(heights)] = as.numeric(heights)/merged.peaks.len[names(heights)]
  
  return(pk2heightMean)
  #return(counts)
})
colnames(samp2HeightsMean$ROADMAP) = names(REF2CELLTYPE)
#rownames(samp2Heights$ROADMAP) = merged.peaks.inRef


#####################################################
print("get counts for blueprint samples")

samp2HeightsMean$blueprint = sapply(names(files.in.blueprint),
FUN=function(samp)
{
  print(samp)
  f.i.bw = files.in.blueprint[[samp]]["bw"]
  f.tmp.heights = paste(dir.tmp, samp, ".bw2height.txt", sep="")
  
  cmd = paste("bigWigAverageOverBed ", f.i.bw, " ", file.ou.filt.mergedPeak, " ", f.tmp.heights, sep="")
  system(cmd)
  
  heights  =as.matrix(read.table(f.tmp.heights, sep="\t", header=F, row.names=1))[,4]
  return(heights)
})
colnames(samp2HeightsMean$blueprint) =  names(files.in.blueprint)
# rownames(samp2Heights$blueprint) =  merged.peaks.inRef


save(samp2HeightsMean, qc.info, merged.peaks, merged.peaks.df, merged.peaks.len, file=file.ou.RData)
#save(samp2HeightsMean, qc.info, ref2pkOvlpWithTest, pk2refOvlpTest.Weightsum, merged.peaks.inRef, merged.peaks, merged.peaks.df, merged.peaks.len, file=file.ou.RData)

##################################################
print("get peaks overlap with reference peaks")


pdf(file.ou.ovlp.pdf, height=6, width=12)
layout(matrix(1:12, ncol=4))
ref2pkOvlpWithTest = sapply(c(names(REF2CELLTYPE), names(files.in.blueprint)),
FUN=function(ref)
{
  f.i.pk = ""
  if(grepl("^E",ref))
  {
    f.i.pk=sub("E000", ref, file.in.ref.pk)
  }else if(grepl("^BP",ref))
  {
    f.i.pk = files.in.blueprint[[ref]]["pk"]  
  }
  print(f.i.pk)
  
  cmd = paste("intersectBed -b ", file.ou.filt.mergedPeak, " -a ", f.i.pk, " -wo|awk '{print $3-$2 \"\t\" $NF \"\t\" $(NF-1)  }'", sep="")
  pksOvlp.info=t(sapply(system(cmd, intern = T), FUN=function(x){unlist(strsplit(x, split="\t"))}))
  smoothScatter(as.numeric(pksOvlp.info[,1]), as.numeric(pksOvlp.info[,2])/as.numeric(pksOvlp.info[,1]), 
                xlab="peak length in ref", ylab="fraction of each ref peak\noverlapping by test",
                nbin = 500,
                main=ref) 
  
  #cmd = paste("intersectBed -a ", file.ou.filt.mergedPeak, " -b ", f.i.pk, " -u|awk '{print $1 \":\" $2 \"-\" $3}'  ", sep="")
  pks.isOvlp=rep(0, length(merged.peaks))
  names(pks.isOvlp) = merged.peaks
  
  pks.isOvlp[levels(factor(pksOvlp.info[,3]))]=1
  
  
  
  return(pks.isOvlp)
})
colnames(ref2pkOvlpWithTest) = c(names(REF2CELLTYPE), names(files.in.blueprint))
dev.off()

weights = rep(1, ncol(ref2pkOvlpWithTest))
names(weights) = colnames(ref2pkOvlpWithTest)
for(ct in names(blueprint.nm2celltypes))
{
  weights[grepl(paste("^", ct, sep=""), colnames(ref2pkOvlpWithTest))] = 1/sum(grepl(paste("^", ct, sep=""), colnames(ref2pkOvlpWithTest))) #from the same cell type, share the same weights  
}
weights["E040"] <- weights["E037"] <- 0.5
weights["E038"] <- weights["E039"] <- 0.5
pk2refOvlpTest.Weightsum = apply(ref2pkOvlpWithTest,1,
FUN=function(x)
{
  sum(weights[x==1])
})
merged.peaks.inRef = merged.peaks[pk2refOvlpTest.Weightsum>=0.5]

merged.peaks.inRef.bed = gsub(":", "\t", merged.peaks.inRef)
merged.peaks.inRef.bed = gsub("-", "\t", merged.peaks.inRef.bed)
write(paste(merged.peaks.inRef.bed, merged.peaks.inRef, sep="\t"), file.tmp.pkInRef)


save(samp2HeightsMean, qc.info, ref2pkOvlpWithTest, pk2refOvlpTest.Weightsum, merged.peaks.inRef, merged.peaks, merged.peaks.df, merged.peaks.len, file=file.ou.RData)


#proporiton of reads
colnames(ref2pkOvlpWithTest) = gsub("BP-", "BP",colnames(ref2pkOvlpWithTest) )
ref2pkOvlpWithTest.refSlct = ref2pkOvlpWithTest[,grepl(paste(REFs.SELECT, collapse="|"), colnames(ref2pkOvlpWithTest))]

ref2pkOvlpWithTest.refSlct.combBP = sapply(REFs.SELECT,
FUN=function(ref)
{
  if(ref %in% colnames(ref2pkOvlpWithTest.refSlct))
  {
    return(ref2pkOvlpWithTest.refSlct[,ref])
  }else
  {
    samps = colnames(ref2pkOvlpWithTest.refSlct)
    samps = samps[grepl(ref, samps)]
    apply(ref2pkOvlpWithTest.refSlct[,samps], 1, mean)
  }
  
})


peaks = list()
peaks$merged.2goodSample = merged.peaks
peaks$merged.2GS.inRef = merged.peaks.inRef
peaks$merged.2GS.inRef.slcted = rownames(ref2pkOvlpWithTest.refSlct.combBP)[apply(ref2pkOvlpWithTest.refSlct.combBP,1,max)>=0.5]

reads.total = lapply(samp2HeightsMean,
FUN=function(samps.heightMean)
{
  res= sapply(peaks,
  FUN=function(pks)
  {
    rbind(merged.peaks.len[pks]) %*% samps.heightMean[pks, ]
  })
  
})

reads.prop = lapply(reads.total,
FUN=function(rds.tot)
{
  rds.tot/rds.tot[,1]
})

#
library(ggplot2)
df = data.frame(reads.prop = c(reads.prop$Mayo[, "merged.2GS.inRef"], reads.prop$Mayo[, "merged.2GS.inRef.slcted"],
                                     reads.prop$ROADMAP[, "merged.2GS.inRef"], reads.prop$ROADMAP[, "merged.2GS.inRef.slcted"],
                                     reads.prop$blueprint[, "merged.2GS.inRef"], reads.prop$blueprint[, "merged.2GS.inRef.slcted"]),
                datasets =c(rep("Mayo", 2*nrow(reads.prop$Mayo)), 
                            rep("ROADMAP", 2*nrow(reads.prop$ROADMAP)),
                            rep("Blueprint", 2*nrow(reads.prop$blueprint))),
                type = c(rep("reads from peaks overlapping those from any reference samples", nrow(reads.prop$Mayo)), 
                         rep("reads from peaks overlapping those from selected reference samples", nrow(reads.prop$Mayo)),
                         rep("reads from peaks overlapping those from any reference samples", nrow(reads.prop$ROADMAP)), 
                         rep("reads from peaks overlapping those from selected reference samples", nrow(reads.prop$ROADMAP)),
                         rep("reads from peaks overlapping those from any reference samples", nrow(reads.prop$blueprint)), 
                         rep("reads from peaks overlapping those from selected reference samples", nrow(reads.prop$blueprint))
                         ),
                stringsAsFactors=F)
df$datasets=factor(df$datasets, levels=c("Mayo", "ROADMAP", "Blueprint"))
df$type=factor(df$type, levels=c("reads from peaks overlapping those from any reference samples", 
                                 "reads from peaks overlapping those from selected reference samples"))


pdf(file.ou.readPerct.pdf, width =12, height=5)

p <-ggplot(df, aes(x=datasets, y=reads.prop, fill=type)) +
  geom_violin() + 
  ggtitle("compared to reads in all filtered merged peaks") +
  xlab("datasets") + 
  ylab("reads proportion")

print(p)
# # Change the position
# p<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=supp)) +
#   geom_violin(position=position_dodge(1))
# p

dev.off()


save(samp2HeightsMean, qc.info, ref2pkOvlpWithTest, pk2refOvlpTest.Weightsum, peaks, merged.peaks.df, merged.peaks.len, reads.total, reads.prop, file=file.ou.RData)
