#V1.5
#a.2.extractRefSignalForMergedPeaksFromEpimap_V1.5_HM.R

#V1.4
#filter signal below threshoold for bw file
#so that signal extraction by bigWigAverageOverBed would be more acurate on mean column

#V1.1
#use peaks merged across all four tissues

#extract pval signal from Roadmap reference  for peaks identified from each tissue
#

options(scipen = 999)

args = commandArgs(trailingOnly=TRUE)
HM = args[1] #Brain Heart Muscle Lung
# 

#PROJECT="BiP"


sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                #"#$ -P compbio_lab", 
                "#$ -binding linear:1", 
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=10g",
                "#$ -l h_rt=2:00:00", 
                "##$ -pe smp 5",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                #"export LD_LIBRARY_PATH=$HOME/lhou.compbio/libs/:$LIBRARY_PATH:/broad/uge/8.4.3/lib/lx-amd64",
                sep="\n")

BW.CUTOFF=1
# 
# files.in.histon.RData = c(Brain ="/broad/compbio/data/eGTEx/H3K27ac/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2/Brain.heightsMean.depthGCNormed.RData",
#                           Heart="/broad/compbio/data/eGTEx/H3K27ac/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2/Heart.heightsMean.depthGCNormed.RData",
#                           Muscle="/broad/compbio/data/eGTEx/H3K27ac/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2/Muscle.heightsMean.depthGCNormed.RData",
#                           Lung = "/broad/compbio/data/eGTEx/H3K27ac/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2/Lung.heightsMean.depthGCNormed.RData")

file.in.mergedPeak.bed.gz = "a_mergePeaksFromHMs_hg19/3activeHMs.peaks.merged.hg19.bed.gz"#"a_mergePeaksFromHMs_hg19/5HMs.peaks.merged.hg19.bed.gz"#"../peakMergeNormal/c_mergePeaksFromTissues_hg19/4Tiss.peaks.merged.hg19.bed"
files.in.ref.bw = dir(paste0("~/lhou.compbio/data/Epimap/bw_hg19/", HM, "/"), "bigWig$", full.names = T)

file.in.chrom.sizes = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes"

dir.tmp =  paste0("~/hptmp/MayoBipolar/vsEpimap/", HM, "/")
dir.create(dir.tmp, showWarnings = F, recursive = T)


file.tmp.mergedPeaks.bed= paste0(dir.tmp, "3activeHMs.mergedPeaks.bed")


#

cmd = paste0("zcat ", file.in.mergedPeak.bed.gz, "|awk '{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $1 \":\" $2 \"-\" $3 }' >", file.tmp.mergedPeaks.bed)
system(cmd)

for(f.i.bw in files.in.ref.bw)
{
  f.o.sh = paste0("a.2.",basename(f.i.bw), ".exSig.sh")
  write(sh.head, f.o.sh)
  f.o = paste0(dir.tmp, basename(f.i.bw), ".Signal")
  #cmd = paste("bigWigSummaryBatch", f.i.bw, file.tmp.mergedPeaks.bed, "1 > ", f.o, sep=" ")
  
  f.o.bdg = gsub("bigWig", "bdg", f.i.bw)
  f.o.filt.bdg = gsub("bigWig", "filt.bdg", f.i.bw)
  f.o.filt.bw = gsub("bigWig", paste0("filt", BW.CUTOFF, ".bw"), f.i.bw)
  
  
  cmd= paste("bigWigToBedGraph", f.i.bw, f.o.bdg)
  write(cmd, f.o.sh, append = T)
  
  cmd= paste0("awk 'BEGIN {OFS=\"\\t\"}{if($4>", BW.CUTOFF, ") print $1,$2,$3,$4}' ", f.o.bdg, " >", f.o.filt.bdg) #filter those signal
  write(cmd, f.o.sh, append = T)
  
  
  cmd= paste("bedGraphToBigWig ", f.o.filt.bdg,  file.in.chrom.sizes, f.o.filt.bw)
  write(cmd, f.o.sh, append = T)
  
  cmd = paste("bigWigAverageOverBed", f.o.filt.bw, file.tmp.mergedPeaks.bed, f.o, sep=" ")
  write(cmd, f.o.sh, append = T)
  
  
  cmd= paste("rm ", f.o.bdg)
  write(cmd, f.o.sh, append = T)
  
  cmd= paste("rm ", f.o.filt.bdg)
  write(cmd, f.o.sh, append = T)
  
  cmd= paste("rm ", f.o.filt.bw)
  write(cmd, f.o.sh, append = T)
  
  
  write(paste0("echo done!", basename(f.i.bw)), f.o.sh, append = T)
  system(paste("qsub ", f.o.sh))

}



