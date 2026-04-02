
#merge peaks across tissues
#and map merged peaks to peaks in tissue

library(VennDiagram)

# TISS2COLOR = c(Brain ="#E5803C",
#                Heart = "#DD4545",
#                Muscle = "#5FAA37",
#                Lung = "#4474EA")

HMs = c("H3K27ac", "H3K4me1", "H3K4me3")#, "H3K36me3",  "H3K27me3")

files.in.peaks.bed = sapply(HMs,
FUN=function(hm)
{
  f.i = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", hm, ".diffPeak.BG.bed")
})
names(files.in.peaks.bed) = HMs

dir.tmp = "~/hptmp/MayoBipolar/peaksMergedFromHMs/"
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.peaks.all =  paste0(dir.tmp, "allTiss.peaks.bed")
#file.tmp.peaks.all.sorted = paste0(dir.tmp, "allTiss.peaks.sorted.bed")
#file.tmp.peaks.merge = paste0(dir.tmp, "allTiss.peaks.merged")

dir.ou = "./a_mergePeaksFromHMs_hg19/"
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.peak.merged = paste0(dir.ou, "3activeHMs.peaks.merged.hg19.bed")
file.ou.peak.merged.RData = paste0(dir.ou, "3activeHMs.mergedPeak2HMPeak.hg19.RData")
file.ou.mergedPeakProp.RDS = paste0(dir.ou, "3activeHMs.mergedPeak2HMPeakProp.hg19.RData")
file.ou.pdf = paste0(dir.ou, "3activeHMs.mergedPeak2HMPeak.hg19.pdf")
#put peak id into file
if(file.exists(file.tmp.peaks.all)) file.remove(file.tmp.peaks.all)
for(hm in HMs)
{
  cmd = paste0("awk '{print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" \"", hm, "_\" $1 \":\" $2 \"-\" $3", "  }' ", files.in.peaks.bed[hm], ">>", file.tmp.peaks.all)
  system(cmd)
}

cmd = paste0("sortBed -i ", file.tmp.peaks.all, "|mergeBed -i stdin -c 4 -o collapse -delim \",\" >", file.ou.peak.merged)
system(cmd)

cmd = paste0("bgzip ", file.ou.peak.merged, " && tabix -p bed ", file.ou.peak.merged, ".gz")
system(cmd)

#

buf = read.table(gzfile(paste0(file.ou.peak.merged, ".gz")), sep="\t", header=F, row.names=NULL, stringsAsFactors = F)
peaks.len= as.numeric(buf[,3])-as.numeric(buf[,2])

mergedPeak2hmPeaks = lapply(buf[,4],
FUN=function(x)
{
  unlist(strsplit(x, split=","))
})  
names(mergedPeak2hmPeaks) = paste(buf[,1],buf[,2],buf[,3], sep="\t")

hm.mergedPeak2peaks = lapply(HMs,
FUN=function(hm)
{
  lapply(mergedPeak2hmPeaks,
  FUN=function(pks)
  {
    hm.pks =  pks[grepl(hm, pks)]
    hm.pks = gsub(paste0(hm, "_"), "", hm.pks)
    return(hm.pks)
  })
})
names(hm.mergedPeak2peaks)= HMs

save(mergedPeak2hmPeaks, 
     hm.mergedPeak2peaks,
     #chr2composit, 
     #chr.comp2Frac,
     file=file.ou.peak.merged.RData)


# hm.peak2mergedPeak = lapply(hm.mergedPeak2peaks,
# FUN=function(mp2pks)
# {
#   mp2pks.df.list =lapply(names(mp2pks)[sapply(mp2pks, length)>0],
#   FUN=function(mp)
#   {
#     data.frame(mp,  mp2pks[mp])
#   })
#   
#   mp2pks.df = do.call(rbind, mp2pks.df.list)
#   
#   pk2mp = mp2pks.df[,1]
#   names(pk2mp) = mp2pks.df[,2]
#   
#   return(pk2mp)
# })


hm.peak2mergedPeak = lapply(HMs,
FUN=function(hm)
{
  print(hm)
  pk2mps = unlist(lapply(names(mergedPeak2hmPeaks),
  FUN=function(mp)
  {
    pks = mergedPeak2hmPeaks[[mp]]
    hm.pks =  pks[grepl(hm, pks)]
    hm.pks = gsub(paste0(hm, "_"), "", hm.pks)
    if(length(hm.pks)>0)
    {
      pk2mp=c()
      pk2mp[hm.pks] = mp
      return(pk2mp)
    }
    
  }))
  return(pk2mps)
  
})
names(hm.peak2mergedPeak)= HMs

# 
# chr2composit = lapply(split(buf[,4], buf[,1]),
# FUN=function(peaks.ids)
# {
#   hm.comp = sapply(peaks.ids,
#   FUN=function(pks)
#   {
#     tc = c()
#     for(hm in HMs)
#     {
#       if(grepl(hm, pks))
#       {
#         tc = c(tc, hm)
#       }
#     }
#     return(paste(tc, collapse=";"))
#   })
# })
# 
# chr2composit.levels = levels(factor(unlist(chr2composit)))
# 
# chr.comp2Frac = sapply(chr2composit,
# FUN=function(comps)
# {
#   comps.counts =  summary(factor(comps,levels=chr2composit.levels))  
#   comps.fract = comps.counts/sum(comps.counts)
# })


save(mergedPeak2hmPeaks, 
     hm.mergedPeak2peaks,
     hm.peak2mergedPeak,
     #chr2composit, 
     #chr.comp2Frac,
     file=file.ou.peak.merged.RData)

mPk.length=sapply(names(mergedPeak2hmPeaks), FUN=function(x) {buf=as.numeric(unlist(strsplit(x, "\t"))[2:3]); buf[2]-buf[1]})
mPk.pkNum=sapply(mergedPeak2hmPeaks, length)

pdf(file.ou.pdf, width=8, height=6)
boxplot(mPk.length~mPk.pkNum, 
        log="y",
        xlab="number of peaks merged",
        ylab="log10 peak length")
hist(log10(mPk.length),
     xlab="log10 peak length",
     main=paste0(round(sum(mPk.length<5000)/length(mPk.length)*100), "% peaks shorter than 5000bp"))
abline(v=log10(5000))
dev.off()


#proportion of peak in length for each merged peak
# 
# mergedPeak2lengthProp = sapply(names(mergedPeak2hmPeaks),
# FUN=function(m)
# {
#   buf = unlist(strsplit(m , split="\t"))
#   len = as.numeric(buf[3]) - as.numeric(buf[2])
#   
#   len.hms=sapply(HMs,
#   FUN=function(hm)
#   {
#     hm.pks = mergedPeak2hmPeaks[[m]][grepl(hm, mergedPeak2hmPeaks[[m]])]
#     if(length(hm.pks)==0)
#     {
#       return(0)
#     }else
#     {
#       sum(sapply(hm.pks,
#       FUN=function(pk)
#       {
#         buf = unlist(strsplit(pk , split=":|-", perl = T))
#         len = as.numeric(buf[3]) - as.numeric(buf[2])
#       }))
#     }
#   })
#   
#   len.hms/len
# })
# names(mergedPeak2lengthProp) = names(mergedPeak2hmPeaks)
# 
# saveRDS(mergedPeak2lengthProp, file=file.ou.mergedPeakProp.RDS)

# #
# hm.mergedPeaks = lapply(hm.mergedPeak2peaks,
# FUN=function(mPks)
# {
#   names(mPks)[sapply(mPks,length)!=0]
# })
# mPks.size=c(B = length(hm.mergedPeaks$Brain),
#             L = length(hm.mergedPeaks$Lung),
#             H = length(hm.mergedPeaks$Heart),
#             M = length(hm.mergedPeaks$Muscle),
#             
#             
#             BH = length(intersect(tiss.mergedPeaks$Brain,
#                                   tiss.mergedPeaks$Heart)),
#             BM = length(intersect(tiss.mergedPeaks$Brain,
#                                   tiss.mergedPeaks$Muscle)),
#             BL = length(intersect(tiss.mergedPeaks$Brain,
#                                   tiss.mergedPeaks$Lung)),
#             HM = length(intersect(tiss.mergedPeaks$Heart,
#                                   tiss.mergedPeaks$Muscle)),
#             LH = length(intersect(tiss.mergedPeaks$Heart,
#                                   tiss.mergedPeaks$Lung)),
#             LM = length(intersect(tiss.mergedPeaks$Muscle,
#                                   tiss.mergedPeaks$Lung)),
#             
#             BHM = length(intersect(intersect(tiss.mergedPeaks$Brain,
#                                   tiss.mergedPeaks$Heart),
#                                   tiss.mergedPeaks$Muscle)),
#             BLH = length(intersect(intersect(tiss.mergedPeaks$Brain,
#                                   tiss.mergedPeaks$Heart),
#                                   tiss.mergedPeaks$Lung)),
#             BLM = length(intersect(intersect(tiss.mergedPeaks$Brain,
#                                   tiss.mergedPeaks$Muscle),
#                                   tiss.mergedPeaks$Lung)),
#             LHM = length(intersect(intersect(tiss.mergedPeaks$Heart,
#                                   tiss.mergedPeaks$Muscle),
#                                   tiss.mergedPeaks$Lung)),
#             BLHM = length(intersect(intersect(intersect(tiss.mergedPeaks$Heart,
#                                   tiss.mergedPeaks$Muscle),
#                                   tiss.mergedPeaks$Lung),
#                                   tiss.mergedPeaks$Brain))
#             )



# 
# 
# pdf(file.ou.pdf, width=12, height=9)
# # par(xpd=TRUE, mar=c(7,5,5,12))
# # 
# # cols = rainbow(nrow(chr.comp2Frac))
# # chr.comp2Frac=chr.comp2Frac[, paste0("chr", c(1:22, "X"))]
# # barplot(chr.comp2Frac, beside = F, las=2, col=cols, ylab="fraction of peaks", main="composition of merged peaks")
# # legend("topright", inset=c(-0.2,0), fill=cols, legend = rownames(chr.comp2Frac))
# # 
# #grid.draw(venn.plot);
# venn.plot <- draw.quad.venn(
# 	area1 = mPks.size["B"],
# 	area2 = mPks.size["L"],
# 	area3 = mPks.size["H"],
# 	area4 = mPks.size["M"],
# 	n12 = mPks.size["BL"],
# 	n13 = mPks.size["BH"],
# 	n14 = mPks.size["BM"],
# 	n23 = mPks.size["LH"],
# 	n24 = mPks.size["LM"],
# 	n34 = mPks.size["HM"],
# 	n123 = mPks.size["BLH"],
# 	n124 = mPks.size["BLM"],
# 	n134 = mPks.size["BHM"],
# 	n234 = mPks.size["LHM"],
# 	n1234 = mPks.size["BLHM"],
# 	category = c("Brain", "Lung", "Heart", "Muscle"),
# 	col = TISS2COLOR[c("Brain", "Lung", "Heart", "Muscle")],
# 	lty = 1,
# 	lwd=3,
# 	cex = 2,
# 	cat.cex = 2,
# 	cat.col = TISS2COLOR[c("Brain", "Lung", "Heart", "Muscle")]
# 	);
# #
# dev.off()
# 
# 
# 
# 
