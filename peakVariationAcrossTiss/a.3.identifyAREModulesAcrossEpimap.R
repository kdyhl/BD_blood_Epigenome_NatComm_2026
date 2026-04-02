#V1.6.1
#a.3.identifyAREModulesAcrossEpimap_V1.6.1_merge3ActiveMarks.R
#use merged ARE from three active marks

if(!grepl("R version 3.6.", R.Version()$version.string ))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



#
options(scipen = 999)



library(pheatmap)
library(flexclust)

#library(colorRampPalette)
library(Rtsne)
library(ggplot2)

#
ARE.SIGNAL.CUTOFF=2
CLUSTER.SIZE = 2000
#SAMPLES.FILT = c("BSS01866", "BSS00333")
#
dir.in.ref =  "~/hptmp/MayoBipolar/vsEpimap/H3K27ac/"
files.in.ref=dir(path=dir.in.ref, pattern=paste0(".*Signal"), full.names=T)
names(files.in.ref)=sapply(files.in.ref,
FUN=function(f.i)
{
   gsub(".*?H3K27ac_(BSS\\d+)(\\.sub).*", "\\1", basename(f.i), perl=TRUE)
})

files.tmp.avg = sapply(files.in.ref,
FUN=function(f.i)
{
  
  f.tmp  = gsub("Signal", "k27acAvg", f.i)
  if(!file.exists(f.tmp))
  {
    cmd = paste0("awk '{print $5}' " , f.i, "> ", f.tmp)
    system(cmd)
  }
  return(f.tmp)
})


  
#file.in.mergedPeak.bed.gz = "a_mergePeaksFromHMs_hg19/5HMs.peaks.merged.hg19.bed.gz"
file.in.epimap.meta = "~/lhou.compbio/data/Epimap/Epimap_metadata_sampleMatrix.csv"
file.in.epimap.grp2col = "~/lhou.compbio/data/Epimap/Epimap_metadata_group2color.csv"

#file.in.mergedPk2TissPk.RData="a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2tissuePeak.hg19.RData"
#file.in.peakHeights.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.2/", COND, ".heightsMean.depthGCNormed.RData")


dir.ou = "a_3_AREModulesInEpimap_V1.6.1_merge3activeMarks/H3K27ac/"#"b_tissueSpecif/"
dir.create(dir.ou, showWarnings = F, recursive = T)



file.ou.pref=paste0(dir.ou, "3actHMMergedPeaks.inEpimap.H3K27acSignal.", sep="")
file.ou.matrix = paste0(file.ou.pref, "matrix.txt", sep="")
file.ou.RData = paste0(file.ou.pref, "RData", sep="")
file.ou.clu.RData = paste0(file.ou.pref, "cluster.RData", sep="")
file.ou.cluster.pdf= paste0(file.ou.pref, "cluster.pdf", sep="")
file.ou.tsne.pdf= paste0(file.ou.pref, "tsne.pdf", sep="")

file.ou.clu.bed.pref = paste0(dir.ou, "AREs.")
#


# load(file.in.mergedPk2TissPk.RData)
# 
# buf=read.table(gzfile(file.in.mergedPeak.bed.gz), sep="\t", header=F, row.names=4)
# peaks.4TissMerged=paste0(buf[,1], ":", buf[,2], "-", buf[,3])
# peaks.filter = sub("\t", ":", peaks.filter.bed)
# peaks.filter = sub("\t", "-", peaks.filter)

# getSingleFilesIntoMatrix2=function(files.in, f.o)#paste in column wise 
# {
#   write(paste(names(files.in), collapse="\t"), f.o)
#   
#   cmd = paste("paste ", paste(files.in, collapse=" "), "  | column -s \"\\t\" -t  >> ", f.o, sep="")
#   print(cmd)
#   system(cmd)
# 
#   
# }

# getSingleFilesIntoMatrix2(files.tmp.avg, file.ou.matrix)

#


if(!file.exists(file.ou.RData))
{
  signal.matrix=sapply(files.tmp.avg,
  FUN=function(f.i)
  {
    read.table(f.i, sep="\t", header = F, row.names = NULL, colClasses = "numeric")[,1]
  })
  rownames(signal.matrix) = system(paste0("awk '{print $1}' ", files.in.ref[1]), intern = T)


  signal.matrix.filter=signal.matrix[!grepl("chrY", rownames(signal.matrix)), ] #[peaks.filter,]
  #signal.matrix.filter = signal.matrix.filter[, setdiff(colnames(signal.matrix.filter), SAMPLES.FILT)]

#RPKMs.cmb = cbind(RPKMs$RNA,PKMs$m6A[rownames(RPKMs$RNA),colnames(RPKMs$RNA)])



  
  #
  signal.matrix.filter.avgAndVar = t(apply(signal.matrix.filter, 1,
  FUN=function(x)
  {
    return(c(sd(x), mean(x)))
  }))
  colnames(signal.matrix.filter.avgAndVar) = c("sd", "mean")
  
  signal.matrix.filter.binary = signal.matrix.filter
  signal.matrix.filter.binary[signal.matrix.filter>=ARE.SIGNAL.CUTOFF]=1
  signal.matrix.filter.binary[signal.matrix.filter<ARE.SIGNAL.CUTOFF]=0
  
  save(#signal.matrix,
       signal.matrix.filter,
       #peaks.filter,
       #ref.samp2color,
       signal.matrix.filter.avgAndVar,
       signal.matrix.filter.binary,
       #dists, 
       file = file.ou.RData)
}else
{
  load(file.ou.RData)
}
###########################################

ref.meta = read.table(file.in.epimap.meta, sep=",", header=T, row.names=3, comment="", stringsAsFactors = F)[colnames(signal.matrix.filter), 1:12]
ref2grp= ref.meta$Group
names(ref2grp)=rownames(ref.meta)
ref2sampInfo= ref.meta$Extended.Info
names(ref2sampInfo)=rownames(ref.meta)

ref.ID2sampName= paste0(rownames(ref.meta), "_", ref2sampInfo)
names(ref.ID2sampName) = rownames(ref.meta)

ref.grp2col = read.table(file.in.epimap.grp2col, sep=",", header = T, row.names = 1, stringsAsFactors = F, comment.char = "")
ref.samp2color = ref.grp2col[ref2grp, "color"]
names(ref.samp2color) = ref.ID2sampName

# 
print("KCCA")

pks.kept = rownames(signal.matrix.filter.binary)[apply(signal.matrix.filter.binary!=0, 1, any)]

AREs.kcca =kcca(signal.matrix.filter.binary[pks.kept, ], 
               #k=round(nrow(signal.matrix.filter.binary)/CLUSTER.SIZE), 
               k=round(length(pks.kept)/CLUSTER.SIZE), 
               family=kccaFamily("ejaccard"),
               simple=T)


save(#matrix,
     signal.matrix.filter,
     #peaks.filter,
     #ref.samp2color,
     signal.matrix.filter.avgAndVar,
     signal.matrix.filter.binary,
     pks.kept,
     AREs.kcca,
     file = file.ou.RData)


#ARE center
ARE.module.center = t(AREs.kcca@centers)
colnames(ARE.module.center)=paste0("clu_", 1:ncol(ARE.module.center))
# rownames(ARE.module.center) = paste(rownames(ARE.module.center), ref2sampInfo[rownames(ARE.module.center)], sep="_")
# 
rownames(ARE.module.center) = ref.ID2sampName[rownames(ARE.module.center)]
#
#ARE.module.center.all0 = apply(ARE.module.center, 2, sum)==0
#ARE.module.center.filt = ARE.module.center[, !ARE.module.center.all0]

#
cluster2ARE = split(names(AREs.kcca@cluster), AREs.kcca@cluster)
names(cluster2ARE) = paste0("clu_", names(cluster2ARE))
for(clu in names(cluster2ARE))
{
  file.ou.clu.bed = paste0(file.ou.clu.bed.pref, clu, ".bed")
  write(gsub("-|:", "\t", cluster2ARE[[clu]]), file.ou.clu.bed)
}
write(gsub("-|:", "\t", names(AREs.kcca@cluster)), 
        paste0(file.ou.clu.bed.pref, "all.bed"))


#
clu2size = sapply(cluster2ARE, length)

AREs.module.ubiqScore = apply(ARE.module.center, 2,
FUN=function(x)
{
  sum(x>=0.2)/length(x)
})

# AREs.module.GrpPresScore = apply(ARE.module.center, 2, 
# FUN=function(x)
# {
#   sum(sapply(cluster2samp,
#   FUN=function(samps)
#   {
#     max(x[samps])
#   })>=0.5)/length(cluster2samp)
#   # sapply(cluster2samp,
#   # FUN=function(samps)
#   # {
#   #   max(x[samps])
#   # })
# })

module.ubiq = colnames(ARE.module.center)[AREs.module.ubiqScore>=0.5]
#module.new  = colnames(ARE.module.center)[AREs.module.ubiqScore==0]
modules.multiTiss = colnames(ARE.module.center)[AREs.module.ubiqScore<0.5 & AREs.module.ubiqScore>=0.1]
modules.rest = setdiff(colnames(ARE.module.center), c(module.ubiq,  modules.multiTiss)) #module.new,


labelCol = function(x)
{
  if (is.leaf(x)) 
  {
    ## fetch label
    label <- attr(x, "label")
    attr(x, "nodePar") <- list(lab.col=ref.samp2color[label])
    #attr(x, "label") = paste(label, samp2tissue[label], sep="-")
  }
  return(x)
}

draw.hclu = function(dists, main)
{
  hc<-hclust(dists, method ="complete")
  d <- dendrapply(as.dendrogram(hc), labelCol)
  #par(cex=1)
  par(mar=c(30,5,5,5))
  plot(d ,main = main, cex.main=2, cex.axis=1)
  
  return(hc)
  #legend("topleft", legend=TISSUES, col=TISS2COL, text.col=TISS2COL, pch="o", cex=4)
}



#visualization
pdf(file.ou.cluster.pdf, width=30, height=16)
#
# dists = as.dist(1-cor(t(ARE.module.center[,modules.rest]), method = "pearson"))
# Roadmap.samp.hclu.res=draw.hclu(dists,
#           main=paste0("hierarchical struture of Roadmap samples\nbased on ARE modules from eGTEx") #, COND, " H3K27ac ChIP-seq samples")
# )
# #
# dists = as.dist(1-cor(t(ARE.module.center[,c(modules.multiTiss, modules.rest)]), method = "pearson"))
# Roadmap.samp.hclu.res=draw.hclu(dists,
#                                 main=paste0("hierarchical struture of Roadmap samples\nbased on ARE modules from eGTEx") #, COND, " H3K27ac ChIP-seq samples")
# )
# # 
dists = as.dist(1-cor(t(ARE.module.center[,c(modules.rest)]), method = "pearson"))#c(modules.multiTiss, modules.rest)
Epimap.samp.hclu.res=draw.hclu(dists, 
                                main=paste0("hierarchical struture of Roadmap samples\nbased on ARE modules from eGTEx") #, COND, " H3K27ac ChIP-seq samples")
)

# dists = as.dist(1-cor(t(ARE.module.center[,c(module.ubiq, modules.multiTiss, modules.rest)]), method = "pearson"))#c(modules.multiTiss, modules.rest)
# Epimap.samp.hclu.res=draw.hclu(dists, 
#                                main=paste0("hierarchical struture of Roadmap samples\nbased on ARE modules from eGTEx") #, COND, " H3K27ac ChIP-seq samples")
# )
# Roadmap.samp.hclu.res=draw.hclu(dists, 
#                                 main=paste0("hierarchical struture of Roadmap samples\nbased on ARE modules from eGTEx") #, COND, " H3K27ac ChIP-seq samples")
# )
#order of AREs modules
modules.multiTiss.col.order =   hclust(as.dist(1-cor(ARE.module.center[, modules.multiTiss], method="pearson")), method ="complete")$order
modules.rest.col.max = apply(ARE.module.center[Epimap.samp.hclu.res$order, modules.rest], 2, which.max)
modules.rest.col.order = order(modules.rest.col.max, decreasing = F)

ARES.modules.ordered=  ARE.module.center[, c(module.ubiq, 
                                             modules.multiTiss[modules.multiTiss.col.order], 
                                             #module.new, 
                                             modules.rest[modules.rest.col.order])]


pheatmap(ARES.modules.ordered,
         cluster_rows = Epimap.samp.hclu.res,
         cluster_cols=F
         )

dev.off()


save(#matrix,
  signal.matrix.filter,
  #peaks.filter,
  #ref.samp2color,
  signal.matrix.filter.avgAndVar,
  signal.matrix.filter.binary,
  pks.kept,
  AREs.kcca,
  cluster2ARE,
  ARE.module.center,
  clu2size,
  ARES.modules.ordered,
  file = file.ou.RData)


# manually assign group name
samp2AREModGrp = c(rep("Blood/immune", 32),#1-32
                   rep("Epithelial/Epidermal_cells", 19), #33-51
                   rep("Placenta", 16), # 52-67
                   rep("Digestive", 15),  #68-82
                   rep("Epithelial_cancer", 19), #83-101
                   rep("Organ", 32), #102-133
                   rep("Fibroblast/Aorta", 26), #134-159
                   rep("Muscle", 31),#160-190
                   rep("Spleen/Adrenal", 15),#191-205
                   rep("Brain/neuron", 18), #206-223
                   rep("ESC/iPSderive", 19)#224-242
                   
)


names(samp2AREModGrp) = rownames(ARES.modules.ordered)[Epimap.samp.hclu.res$order]
modules.rest2AREModGrp = samp2AREModGrp[modules.rest.col.max]
names(modules.rest2AREModGrp) = names(modules.rest.col.max)


ARE.modName.df =data.frame(id = colnames(ARES.modules.ordered),
                           "ARE.module.groups" = c(
                             "Ubiquitious",
                             rep("Multi-tissue", length(modules.multiTiss)),
                             #"Newly detected",
                             modules.rest2AREModGrp[modules.rest[modules.rest.col.order]]
                           ),
                           stringsAsFactors = F)
rownames(ARE.modName.df) = ARE.modName.df$id
# 
# 
# 
# 
# #luster2samp = split(names(samp2AREModGrp), samp2Cluster)
# 
# #colnames(ARES.modules.ordered) = c(module.ubiq, colnames(module.nonUbiq)[module.nonUbiq.col.order])
# 
# 
AREandCluster.df = do.call(rbind, lapply(names(cluster2ARE), FUN=function(clu) data.frame(cluster=clu, ARE=cluster2ARE[[clu]], stringsAsFactors = F)))
AREandCluster.df$modGrp = ARE.modName.df[AREandCluster.df$cluster, "ARE.module.groups"]
ARE2modGrp = AREandCluster.df$modGrp
names(ARE2modGrp)= AREandCluster.df$ARE
  



grpName2mARE = split(names(ARE2modGrp), ARE2modGrp)
for(nm in names(grpName2mARE))
{
  nm.1 = gsub("/| ", "_", nm)
  print(nm.1)
  mAREs.bed = gsub(":|-", "\t", grpName2mARE[[nm]])
  
  f.o.bed=paste0(dir.ou, "AREGrp.", nm.1, ".bed")
  write(mAREs.bed, file=f.o.bed)
}



save(pks.kept,
     ARE.module.center,
     cluster2ARE,
     clu2size,
     ARES.modules.ordered,
     samp2AREModGrp,
     ARE.modName.df,
     ARE2modGrp,
     module.ubiq,
     modules.multiTiss,
     modules.rest,
     Epimap.samp.hclu.res,
     file=file.ou.clu.RData)

#
file.ou.heatmap.tsv=paste0(file.ou.pref, "heatmap.tsv")
write.table(t(ARES.modules.ordered),
            file=file.ou.heatmap.tsv,
            sep="\t",
            quote=F)