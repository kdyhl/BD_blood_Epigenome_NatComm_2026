#build link between HiC regions
#anotated with different Peak sets from different marks


#args = commandArgs(trailingOnly=TRUE)

#
#library(Sushi)
library(igraph)
library(ggplot2)
#
#GENES.CANDIDATE = c("IL23R", "ATG16L1", "IL17RB", "CD83",
                    # "RPS6KA2", "CCR6", "GPR141", "GPR65",
                    # "IL27", "NOD2", "TOX3", "ACSF3")


#DISEASES = c("CD", "UC", "RA")
# DIS2COLOR = c("brown", "purple", "orange")
# names(DIS2COLOR) = DISEASES
PRMOTER.WINDOW = 2000
HIC.DIS.CUTOFF =1e6
#HAQTL.COLs = c(rgb(1, 0, 0, 0.5), rgb(1, 1, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 1, 1, 0.5), rgb(0, 0, 1, 0.5), rgb(1, 0, 1, 0.5))

HMs= c("H3K36me3", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3")

file.in.hic = "~/lhou.compbio/data/HiC/Promoter_HiC_blood_bluePrint_2016/CHiCAGO5/allCellTypes.links.merged.uniq.sorted.bed.gz"
file.in.hic.regions = "~/lhou.compbio/data/HiC/Promoter_HiC_blood_bluePrint_2016/CHiCAGO5/allCellTypes.regions.sorted.merged.bed"

file.in.TSS.bed = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.TSSUp0Down1bp.bed"
# files.in.TssA = c(E047= "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE047TssAFlank.bed",
#                   E048="~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE048TssAFlank.bed",
#                   E039 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE039TssAFlank.bed",
#                   E040 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE040TssAFlank.bed",
#                   E029 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE029TssAFlank.bed",
#                   BPNeutrop = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE030TssAFlank.bed",#
#                   E032 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE032TssAFlank.bed"
#                   )

files.in.peak.bed = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.BG.bed")
names(files.in.peak.bed) = HMs

dir.tmp = paste0("~/hptmp/MayoBipolar/a_HiC_bulk/")
dir.create(dir.tmp, showWarnings = F, recursive = T)

file.tmp.hicRegsOvlpTss.bed = paste(dir.tmp, "hicRegionsOvlpTssA.bed", sep="")
file.tmp.hicRegsOvlpPk.bed = paste(dir.tmp, "hicRegionsOvlpTssA.bed", sep="")
  
dir.ou=paste0("a_hiCLinks_withPeakAnnot_bulk_DP.BG_limma_V1.4.1/")
dir.create(dir.ou, showWarnings = F, recursive = T)

#file.out.pref = paste(dir.ou, "hiCBlock.", REF, sep="")
file.ou.hiCBlock.RData = paste(dir.ou, "hiCLinks.RData", sep="")
file.ou.network.pdf = paste(dir.ou, "hiCNetwork.pdf", sep="")
#file.ou.haQTLDiffPks2GWAS.RData = paste(file.out.pref, ".haQTLDiffPks2GWAS.RData", sep="")
#file.ou.geneNeighb = paste(file.out.pref, ".txt", sep="")


################################################################################
print("hiC region with annotation")
#HIC regions overlap with Tss
hicRegs = rownames(read.table(gzfile(file.in.hic.regions), sep="\t", header = F, row.names = 4))

cmd = paste0("bedtools window -w ", PRMOTER.WINDOW, " -a ", file.in.TSS.bed, " -b ", file.in.hic.regions, ">", file.tmp.hicRegsOvlpTss.bed)
system(cmd)
buf=read.table(file.tmp.hicRegsOvlpTss.bed, sep="\t", header = F, row.names = NULL, stringsAsFactors = F)
buf.hicReg2TSS = sapply(split(buf[,5], buf[,12]), FUN=function(x) paste(levels(factor(x)), collapse=","))

hicReg2TSS = rep("", length(hicRegs))
names(hicReg2TSS)= hicRegs
hicReg2TSS[names(buf.hicReg2TSS)] = buf.hicReg2TSS

#HIC regions overlap with peak
hicReg2Pk.hg19 = sapply(files.in.peak.bed,
FUN=function(f.i.bed)
{
  cmd = paste0("intersectBed -wa -wb -a ", file.in.hic.regions, " -b ", f.i.bed, ">", file.tmp.hicRegsOvlpPk.bed)
  system(cmd)
  buf=read.table(file.tmp.hicRegsOvlpPk.bed, sep="\t", header = F, row.names = NULL, stringsAsFactors = F)
  buf.r2pks = sapply(split(paste0(buf[,5], ":", buf[,6], "-", buf[,7]), buf[,4]), paste, collapse=",")
  
  hicReg2Pks = rep("", length(hicRegs))
  names(hicReg2Pks)= hicRegs
  hicReg2Pks[names(buf.r2pks)] = buf.r2pks
  
  
  return(hicReg2Pks)
  
})

hicRegs.pos = t(sapply(hicRegs,
FUN=function(x)
{
  buf=unlist(strsplit(x, ":|-", perl=T))
}))

hicReg.df =data.frame(id =hicRegs,
                      TSS = hicReg2TSS, 
                      hicReg2Pk.hg19,
                      chr = hicRegs.pos[,1],
                      loc = (as.numeric(hicRegs.pos[,2])+as.numeric(hicRegs.pos[,3]))/2,
                      length= as.numeric(hicRegs.pos[,3])-as.numeric(hicRegs.pos[,2]),
                      stringsAsFactors = F)
################################################################################
print("generating HiC region network") ###############################################

hic.bed = read.table(gzfile(file.in.hic), sep="\t", header = F, row.names=NULL, stringsAsFactors=F)
hic.bed.filt = hic.bed[hic.bed[,1] == hic.bed[,5],]
hicReg.link.df = data.frame(from=hic.bed.filt[,4], 
                            to =hic.bed.filt[,8], 
                            distance= abs((hic.bed.filt[,2]+hic.bed.filt[,3])/2-(hic.bed.filt[,6]+hic.bed.filt[,7])/2 ),
                            stringsAsFactors = F)

regNet <- graph_from_data_frame(d=hicReg.link.df, vertices=hicReg.df, directed=F)
regNet=simplify(regNet, remove.multiple = TRUE, edge.attr.comb = "first")

# regNet=set_edge_attr(regNet, name="distance", index = E(regNet), 
#               value = sapply(attr(E(regNet),"vnames"),
#               FUN=function(vnames)
#               {
#                 vs = unlist(strsplit(vnames, split="|", fixed=T))
#                 return(abs(V(regNet)[vs[1]]$loc - V(regNet)[vs[2]]$loc))
#               }))


save(hicReg.df, hicReg.link.df, 
     regNet,
     file=file.ou.hiCBlock.RData)



#overall properties
pdf(file.ou.network.pdf, height=12, width=20)
layout(matrix(1:12, ncol=4))

hist(log10(hicReg.df$length))
hist(log10(E(regNet)$distance), xlab="log10 distance", main="edge distance", breaks = 100)



#remove edges based on distance and test how many components left, and median size of components

cmpsize = sapply(c(10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000),
FUN=function(dis.cutoff)
{
  print(dis.cutoff)
  #nodes = unlist(lapply(attr(E(regNet)[E(regNet)$distance<=dis.cutoff], "vnames"), FUN=function(vname){unlist(strsplit(vname, split="|", fixed=T))}))
  #regNet.sub = induced_subgraph(regNet, nodes)
  regNet.sub = subgraph.edges(regNet, eids =E(regNet)[E(regNet)$distance<=dis.cutoff], delete.vertices = T)
  print(paste(dis.cutoff, length(V(regNet.sub)), length(E(regNet.sub)), sep=" "))
  regNet.sub.cmpnts = components(regNet.sub)
  hist(regNet.sub.cmpnts$csize, 
       main = paste("distance cutoff ", dis.cutoff, ": ", length(regNet.sub.cmpnts$csize), " components\n(no singleton)", sep=""), 
       breaks=1000)
  
  return(length(regNet.sub.cmpnts$csize))
})
plot(log10(c(10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000)),
     cmpsize,
     xlab="distance cutoff",
     ylab="number of components")

dev.off()
