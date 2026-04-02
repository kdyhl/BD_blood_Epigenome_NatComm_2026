#V1.2.3
#use matched tissue specific chromatin states rather than over all samples
#only distinguish proportion of enhancer or promoter


#V1.2.1
#flip 90
#add columns of novel ARE module
#add annotation with haDHS


#V1.2
#stretch Annot size
#proportion into bar
#bag of word

#V1.1 use complex heatmap

#modularity from Epimap

#annotate the modules 
#a) number of AREs for each module
#b) module intensity across each eGTEx tissues
#c) proportion of eGTEx tissue AREs
#d) median length of AREs
#e) genomic annotation. such as 3' UTR, 5'UTR, promoter or intergenic region
#f) haQTL proportion for each tissue
#g) ARE  activity in single cell ATAC in Brain
#h) genes near ARE  in scRNA in heart
 

# # for annotation
# if(!grepl("R version 3.6.", R.Version()$version.string ))
# {
#   stop("use .r-3.6.0-bioconductor")
# }
# 
# R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
# R.LIBS.PATH = .libPaths()
# .libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



# 

# if(!grepl("R version 3.6", R.version$version.string))
# {
#   stop("use R-3.6")
# }

# R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
# R.LIBS.PATH = .libPaths()
# .libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

#for final results heatmap
if(!grepl("R version 4.0", R.version$version.string))
{
  stop("use R-4.0")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))
# 
library(ComplexHeatmap)
library(annotatr)
library(ggplot2)
library(circlize)
# library(pheatmap)
# library(fgsea)
library(scales)
source("~/lhou.compbio/codes/mylib/RCode/auxfunctions_bagOfword.Heatmap_byCarles.R")
colfunc <- colorRampPalette(c("white", "darkred"))
COLs.MATRIX= colfunc(10)
#
colfunc <- colorRampPalette(c("white", "aquamarine4"))
COLs.GENOMIC= colfunc(10)

ARE2GENE.RNAGE=10000
#SCRNA.GENEs.NUM.CUTOFF=3000

# CS.GROUP=list(promoter="TssA|TssFlnk|TssFlnkU|TssFlnkD|TssBiv",
#               enhancer="EnhG1|EnhG2|EnhA1|EnhA2|EnhWk|EnhBiv")

#


  
ARE.ModGrp2COLOR = c("Ubiquitious"="#FF65AC",
                     "Multi-tissue" ="#D575FE",
                     "Blood/immune" ="#F8766D",
                     "Epithelial/Epidermal_cells" = "#24B700",
                     "Placenta" = "#D16FAA",
                     "Digestive" = "#BE9C00",
                     "Epithelial_cancer" ="#00ACFC",
                     "Organ" = "#4474EA",
                     "Fibroblast/Aorta" ="#00BE70",
                     "Muscle" = "#D85D5D",
                     "Spleen/Adrenal"="#D6D65F",
                     "Brain/neuron"="#E5803C",
                     "ESC/iPSderive" = "#8CAB00")

# 
# ARE.modName.df =data.frame(id = colnames(ARES.modules.ordered),
#                            "ARE.module.groups" = c(
#                              "Ubiquitious",
#                              rep("Multi-tissue", 18), 
#                              "Newly detected",
#                              #rep("Other", 2),
#                              rep("Brain/neuron",15),
#                              rep("ES/iPS", 6),
#                              rep("Digestive", 15),
#                              #rep("other", ),
#                              rep("Organs", 14),
#                              #rep("Lung & inner organs 2", 9),
#                              rep("Muscle/heart", 28),
#                              rep("Blood/immune", 14),
#                              rep("Fibroblast", 12),
#                              rep("Placenta", 4),
#                              rep("Epithelial", 5),
#                              rep("Lung cancer", 9)
# 
#                            ),
#                            stringsAsFactors = F)
# rownames(ARE.modName.df) = ARE.modName.df$id


HM2COLOR = c(H3K27ac ="orange",
             # H3K36me3 = "green",
             # H3K27me3 = "blue",
              H3K4me1 = "purple",
              H3K4me3 = "red")




#
dir.in.bed =paste0("a_3_AREModulesInEpimap_V1.6.1_merge3activeMarks/H3K27ac/")# paste0("a_3_AREModulesInRoadmap_V1.2_addRefPeaks/", COND, "/")

file.in.bg.peak.bed = paste0(dir.in.bed, "AREs.all.bed")
files.in.clu.bed=dir(dir.in.bed, "AREs.clu.*?bed", full.names = T)
names(files.in.clu.bed) = sapply(files.in.clu.bed,
FUN=function(f.i)
{
  gsub(paste0("AREs.|.bed"), "", basename(f.i))
})
file.in.clu.RData= paste0(dir.in.bed, "3actHMMergedPeaks.inEpimap.H3K27acSignal.cluster.RData")
#
file.in.mergedPk2HMPk.RData="a_mergePeaksFromHMs_hg19/3activeHMs.mergedPeak2HMPeak.hg19.RData"
#
# files.in.gARE.RData = paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_withGTExMinCov_V1.2.2_gINT_peer_rmAgeBatch_10k/", c("Brain", "Heart", "Muscle"), ".haQTLPeak.fdr0.05.hg38.RData")
# names(files.in.gARE.RData) = c("Brain", "Heart", "Muscle")
# files.in.hg382hg19.RDS = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", c("Brain", "Heart", "Muscle"), ".pks.hg382hg19Uniq.RDS")
# names(files.in.hg382hg19.RDS) = c("Brain", "Heart", "Muscle")
#
files.in.peakHeights.RData = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", names(HM2COLOR), ".MayoWithRef.heightsMean.depthGCNormed.RData")
names(files.in.peakHeights.RData) = names(HM2COLOR)

#file.in.mARE.annotation.RData = "../peakMergeNormal/c_2_mARE_annotation/mARE_annotation.RData"

file.in.hg19.tss.bed = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.TSSUp0Down1bp.geneName.ENSG.bed"

# file.in.CS2color.txt="~/lhou.compbio/data/Roadmap/ROADMAP_colorCodesFor18ChromStates.txt"
# 
# #file.in.haDHS = "~/lhou.compbio/data/genome/humanAcceleratedRegion/humanAcceleratedDHS_GR_2015_JoshuaAkey.hg19.bed"
# 
# file.in.ARE.Activity.RData = "a_3_EpiMapAREModulesACtivity_acrosseGTEx_merged4Tiss/ARE.eGTExInd.signal.RData"
# 
# file.in.meta.txt = "~/lhou.compbio/data/eGTEX-H3K27ac/Batch_summary_6.16.2018/sample.metaAndQC.RSC0.6.totalReads5e+06.RSC0.8OrMatchedForLung.txt"


dir.tmp= "~/hptmp/MayoBipolar/a4_annot_epimap/"
dir.create(dir.tmp, showWarnings = F, recursive = T)


dir.ou = "a_4_modules_annotations_mergedPk3actHMs_Epimap/"
dir.create(dir.ou, showWarnings = F, recursive = T)
#file.ou.clu.txt = paste0(dir.ou, "modules.term.txt")
file.ou.RData = paste0(dir.ou, "modules.annot.RData")
file.ou.size.pdf = paste0(dir.ou, "modules.size.pdf")
file.ou.otherAnnot.pdf = paste0(dir.ou, "modules.otherAnnot.pdf")
file.ou.ARE2ModGrp.RData= paste0(dir.ou, "mARE2ModuleGrpName.RData")

#
load(file.in.clu.RData)
load(file.in.mergedPk2HMPk.RData)
#load(file.in.mARE.annotation.RData)
# 
#Module name

#module size
clu2size = sapply(cluster2ARE[colnames(ARES.modules.ordered)], length)
AREs.totalNum = sum(clu2size)

grp.Num = length(levels(factor(ARE.modName.df$ARE.module.groups)))

# ARE.ModGrp2color = hue_pal()(grp.Num)
# names(ARE.ModGrp2color) = levels(factor(ARE.modName.df$ARE.module.groups))

modGrp2mAREMod = split(ARE.modName.df$id, ARE.modName.df$ARE.module.groups)
modGrp2size = sapply(modGrp2mAREMod, FUN=function(mods){sum(clu2size[mods])})


# mAREMod2GrpName = rep("Unknown", length(cluster2ARE))
# names(mAREMod2GrpName)   = names(cluster2ARE)
mAREMod2GrpName = ARE.modName.df$ARE.module.groups
names(mAREMod2GrpName) = ARE.modName.df$id

mARE2ModGrpName = lapply(names(cluster2ARE),
FUN=function(clu)
{
  mAREs=cluster2ARE[[clu]]
  ARE2GrpName = rep(mAREMod2GrpName[clu], length(mAREs))
  names(ARE2GrpName) = mAREs

  return(ARE2GrpName)
})
mARE2ModGrpName=unlist(mARE2ModGrpName)



save(cluster2ARE,
     ARE.modName.df,
     ARE.ModGrp2COLOR,
     mARE2ModGrpName,
     modGrp2size,
     file=file.ou.ARE2ModGrp.RData)

#stretch a little on the left side for better alignment
polygon.x=c()
x.right.bin = sum(clu2size)/length(clu2size)*9/11 #adjust for cancer and ubiquitious modules
x.adj = sum(clu2size)/40

#clu2size.rev = rev(clu2size)
for(i in 1:length(clu2size))
{
  x.right = x.adj+x.right.bin * c(i-1, i)
  x.left = c(sum(c(0,clu2size)[1:i]), sum(c(0,clu2size)[1:(i+1)]))

  polygon.x=rbind(polygon.x,
                  c(x.right[1], x.right[1],
                    x.left[1], x.left[1],
                    x.left[2], x.left[2],
                    x.right[2], x.right[2]
                    ))
}
polygon.x = as.vector(t(polygon.x))
#
polygon.y=rep(c(0, 1, 4, 5, 5, 4, 1, 0), length(clu2size))

datapoly = data.frame(id = rep(ARE.modName.df$id, each=8),
                 ARE.module.groups = rep(ARE.modName.df[["ARE.module.groups"]], each=8),
                 x=polygon.x,
                 y=polygon.y
                 )
datapoly[["ARE.module.groups"]] = factor(as.vector(datapoly[["ARE.module.groups"]]) )#,
                        #levels=c("ubiquitious", "other", "brain & neuron", "ESC & iPSC & ES-derive", "digestive", "inner organs 1",  "lung & inner organs 2",  "muscle & heart", "blood & immune cells", "fibroblast", "placenta & EEM", "epithelial cell", "cancer"))
pdf(file.ou.size.pdf, width=20, height=4)
p <- ggplot(datapoly, aes(x = x, y = y)) +
  geom_polygon(aes(fill = ARE.module.groups, group = id), color="black") +
  scale_fill_manual(values=ARE.ModGrp2COLOR) +
  theme(axis.text.x=element_text(size=30))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)

dev.off()




annotates = list()

#genomic region annotation

##need library(annotatr) from R-V3.6
annotates$genomicRegions = sapply(colnames(ARES.modules.ordered),
FUN=function(clu)
{
  print(clu)
  f.i.bed= files.in.clu.bed[clu]
  
  AREs.bed = read_regions(con = f.i.bed, genome = 'hg19', format = 'bed')
  annots = c('hg19_genes_intergenic', 'hg19_basicgenes')#,'hg19_genes_intronexonboundaries')
  annotations = build_annotations(genome = 'hg19', annotations = annots)
  AREs.annotated = annotate_regions(regions = AREs.bed,
                                    annotations = annotations,
                                    ignore.strand = TRUE,
                                    quiet = FALSE)
  AREs.annotated.df=data.frame(AREs.annotated, stringsAsFactors = T)
  
  buf=as.data.frame(summarize_annotations(
    annotated_regions = AREs.annotated,
    quiet = TRUE))
  AREs.annotated.summary=buf[,2]
  names(AREs.annotated.summary) =buf[,1]
  AREs.annotated.prop=AREs.annotated.summary[c("hg19_genes_intergenic",
                                               "hg19_genes_1to5kb", 
                                               "hg19_genes_promoters",
                                               "hg19_genes_5UTRs",
                                               "hg19_genes_exons",
                                               "hg19_genes_introns",
                                               "hg19_genes_3UTRs")]/length(AREs.bed@ranges@start)
  
})
annotates$genomicRegions[is.na(annotates$genomicRegions)]=0
rownames(annotates$genomicRegions) = gsub("hg19_genes_", "prop. of AREs in ", rownames(annotates$genomicRegions))



#test the enrichment of intergenic intronic and promoter regions in different modules
genomicRegionsNum.byMod=sweep(annotates$genomicRegions, 2, annotates$moduleSize, FUN="*")
mods.ubiAndMulti = rownames(ARE.modName.df)[grepl("Ubiquitious|Multi-tissue", ARE.modName.df[,"ARE.module.groups"])]
for(nm in c("prop. of AREs in intergenic", "prop. of AREs in promoters",  "prop. of AREs in introns"))
{
  print(nm)
  N.nm = sum(genomicRegionsNum.byMod[nm,])
  N.nm.UandM = sum(genomicRegionsNum.byMod[nm, mods.ubiAndMulti])
  N.UandM= sum(annotates$moduleSize[mods.ubiAndMulti])
  N.all= sum(annotates$moduleSize)
  
  alt="greater"
  if(N.nm.UandM/N.UandM < N.nm/N.all) alt="less"
  p.fisher=fisher.test(matrix(c(N.nm.UandM, 
                                N.nm-N.nm.UandM, 
                                N.UandM-N.nm.UandM,
                                N.all-N.nm-N.UandM+N.nm.UandM
                                ),
                              ncol=2), alternative = alt)$p.val
  
  print(paste(N.nm.UandM/N.UandM, N.nm/N.all, p.fisher))
}


# 
# #ARE length
# annotates$ARE.length = sapply(cluster2ARE[colnames(ARES.modules.ordered)],
# FUN=function(AREs)
# {
#   lengths = sapply(AREs,
#   FUN=function(x)
#   {
#     buf=unlist(strsplit(x, ":|-", perl = T))
#     as.numeric(buf[3])-as.numeric(buf[2])
#   })
#   
#   return(median(lengths))
# })

save(annotates,
     file=file.ou.RData)

#tissue ARE proportion
load(file.in.mergedPk2HMPk.RData)


annotates$HMPeakPropFC = sapply(cluster2ARE[colnames(ARES.modules.ordered)],
FUN=function(AREs)
{
  AREs.bed=gsub(":|-", "\t", AREs)   
  sapply(hm.mergedPeak2peaks,
  FUN=function(mp2pk)
  {
    sum(sapply(mp2pk[AREs.bed], length)!=0)/length(AREs.bed)/(sum(sapply(mp2pk, length)!=0)/length(mp2pk))
  })
})
rownames(annotates$tissPeakPropFC) = paste0("prop. of AREs in ", rownames(annotates$tissPeakPropFC))


#size of each ARE module
annotates$moduleSize = sapply(cluster2ARE[colnames(ARES.modules.ordered)], length)

# #tissue shared or not proportion FC
# annotates$tissSharingPropFC = sapply(cluster2ARE[colnames(ARES.modules.ordered)],
# FUN=function(AREs)
# {
#   AREs.bed=gsub(":|-", "\t", AREs)   
#   sapply(tissSharingGrp2mAREs,
#   FUN=function(pks)
#   {
#     pks.ovlp = intersect(AREs.bed, pks)
#     length(pks.ovlp)/length(AREs.bed)/(length(pks)/length(mergedPeak2tissPeaks))
#   })
# })

#rownames(annotates$tissSharingPropFC) 

#   
# #Epimap chromatin states annotations
# annotates$CS.prop = lapply(names(TISSUE2MatchedEpiMapSamp),
# FUN=function(tiss)
# {
#   samps = TISSUE2MatchedEpiMapSamp[[tiss]]
#   cs.prop=sapply(cluster2ARE[colnames(ARES.modules.ordered)],
#   FUN=function(AREs)
#   {
#    AREs.bed=gsub(":|-", "\t", AREs) 
#    cs = apply(mARE2samp.CS[AREs.bed, samps], 1, paste, collapse="|")
#    is.prom = grepl(CS.GROUP$promoter, cs, perl=T)
#    is.enh = grepl(CS.GROUP$enhancer, cs, perl=T)
#    #is.rest = !(is.prom|is.enh)
#    
#    return(c(promoter.prop = sum(is.prom)/length(AREs.bed),
#             enhancer.prop = sum(is.enh)/length(AREs.bed)))
#             #rest.prop = sum(is.rest)/length(AREs.bed)))
#   })
#   rownames(cs.prop) = paste(tiss, rownames(cs.prop), sep=".")
#   #colnames(cs.prop) = 
#   return(cs.prop)
# })
# names(annotates$CS.prop)=names(TISSUE2MatchedEpiMapSamp)
#   
#   
#eGTEx activity
# load(file.in.eGTEx.ARE.Activity.RData)
# annotates$eGTEx.modActivty = sapply(cluster2ARE[colnames(ARES.modules.ordered)],
# FUN=function(AREs)
# {
#   AREs.bed=gsub(":|-", "\t", AREs)
#   apply(ind2ARE.allTiss.log10.quant.scaled[,AREs.bed], 1, mean)
# })
# meta = read.table(file.in.meta.txt, sep="\t", row.names = 1, header = T)
# tiss2samps = split(rownames(meta), meta$tissue)
# annotates$eGTEx.tiss.modActivty = t(sapply(tiss2samps[c("Brain", "Heart", "Muscle", "Lung")],
# FUN=function(samps)
# {
#  apply(annotates$eGTEx.modActivty[samps,], 2, mean)
# }))
# annotates$eGTEx.modActivty = annotates$eGTEx.modActivty[unlist(tiss2samps[c("Brain", "Heart", "Muscle", "Lung")]),]

# #human accelerated region
# file.tmp.ovlp.bed= paste0(dir.tmp, "haDHS.ovlp.bed")
# cmd = paste0("intersectBed -u -a ", file.in.bg.peak.bed, " -b ", file.in.haDHS, ">", file.tmp.ovlp.bed)
# system(cmd)
# 
# ARE.haDHS.ovlp.bed = read.table(file.tmp.ovlp.bed, header = F, sep=",", row.names = NULL)[,1]
# annotates$haDHS.prop =sapply(cluster2ARE[colnames(ARES.modules.ordered)],
# FUN=function(AREs)
# {
#   AREs.bed=gsub(":|-", "\t", AREs) 
#   length(intersect(AREs.bed, ARE.haDHS.ovlp.bed))/length(AREs.bed)
# })

save(annotates,
     file=file.ou.RData)


####################################################################################################################################
#visualization
#sample annotations




load(file.ou.RData)
samps.nm.sorted = rownames(ARES.modules.ordered)[Epimap.samp.hclu.res$order]
samps.nm.slcted=get_summary_terms(samps.nm.sorted, 3)
samps.nm.slcted.3 = samps.nm.slcted[samps.nm.slcted!=""]
samps.nm.slcted=get_summary_terms(samps.nm.sorted, 4)
samps.nm.slcted.4 = samps.nm.slcted[samps.nm.slcted!=""]
samps.nm.slcted=union(samps.nm.slcted.3, samps.nm.slcted.4)
samps.nm.slcted_at = sapply(samps.nm.slcted, FUN=function(x) {which(x==rownames(ARES.modules.ordered))})




# #visualization



pdf(file.ou.otherAnnot.pdf, width=20, height=12)
par(mar=c(3,3,3,3))
#ha.genomic =rowAnnotation(df=data.frame(t(annotates$genomicRegions)), gp=gpar(col=))

genomic.prop = annotates$genomicRegions[, colnames(ARES.modules.ordered)]
genomic.prop.norm = t(apply(genomic.prop, 1, FUN=function(x) (x-min(x))/(max(x)-min(x))))
ht.annot.1 = Heatmap(genomic.prop.norm, 
               col=COLs.GENOMIC,
               cluster_rows=F,
               cluster_columns=F,
               show_row_names=T,
               show_column_names=F,
               height = unit(3, "cm"))

# ht.annot.2 = Heatmap(annotates$tissSharingPropFC, 
#                col=COLs.GENOMIC,
#                cluster_rows=F,
#                cluster_columns=F,
#                show_row_names=T,
#                show_column_names=F,
#                height = unit(3, "cm"))
#
annot2colorFun = list()
#annot2colorFun$moduleSize =  colorRamp2(c(0, max(annotates$moduleSize)), c("white", "black"))
#
# for(tiss in names(TISS2COLOR))
# {
#   annot2colorFun[[paste0("prop. of AREs in ", tiss)]]=colorRamp2(c(0, max(annotates$tissPeakPropFC[paste0("prop. of AREs in ", tiss),])), c("white", TISS2COLOR[tiss]))
# }
#

# annot2colorFun[["fullyShared"]]=colorRamp2(c(0, max(annotates$tissSharingPropFC["fullyShared",])), c("white", "black"))
# for(tiss in names(TISS2COLOR))
# {
#   annot2colorFun[[paste0(tiss, ".shared")]]=colorRamp2(c(0, max(annotates$tissSharingPropFC[paste0(tiss, ".shared"),])), c("white", TISS2COLOR[tiss]))
#   annot2colorFun[[paste0(tiss, ".unique")]]=colorRamp2(c(0, max(annotates$tissSharingPropFC[paste0(tiss, ".unique"),])), c("white", TISS2COLOR[tiss]))  
# }
#
# for(cs in rownames(annotates$CS.prop))
# {
#   annot2colorFun[[cs]]=colorRamp2(c(min(annotates$CS.prop[cs,]), max(annotates$CS.prop[cs,])), c("white", CS2color[cs]))
# }
# for(tiss in names(TISSUE2MatchedEpiMapSamp))
# {
#   prom.min = min(annotates$CS.prop[[tiss]][paste0(tiss, ".promoter.prop"),])
#   prom.max = max(annotates$CS.prop[[tiss]][paste0(tiss, ".promoter.prop"),])
#   enh.min = min(annotates$CS.prop[[tiss]][paste0(tiss, ".enhancer.prop"),])
#   enh.max = max(annotates$CS.prop[[tiss]][paste0(tiss, ".enhancer.prop"),])
#   
#   annot2colorFun[[paste0(tiss,".promoter.prop")]] = colorRamp2(c(prom.min, prom.max), c("white", "red"))
#   annot2colorFun[[paste0(tiss,".enhancer.prop")]] = colorRamp2(c(enh.min, enh.max), c("white", "orange"))
#   #annot2colorFun[[paste0(tiss,".rest.prop")]] = colorRamp2(c(0,1), c("white", "black"))
# }


ha.AREprop = HeatmapAnnotation(module_size=anno_barplot(annotates$moduleSize))#,
                               # df=data.frame(t(annotates$tissSharingPropFC),
                               #               #t(annotates$tissPeakPropFC),
                               #               t(annotates$CS.prop$Brain), 
                               #               t(annotates$CS.prop$Heart),
                               #               t(annotates$CS.prop$Muscle),
                               #               t(annotates$CS.prop$Lung),
                               #               #t(annotates$eGTEx.modActivty),
                               #               t(annotates$eGTEx.tiss.modActivty),
                               #               check.names = F),
                               # # Brain = ["prop. of AREs in Brain",],
                               # # Heart = annotates$tissPeakPropFC["prop. of AREs in Heart",],
                               # # Muscle = annotates$tissPeakPropFC["prop. of AREs in Muscle",],
                               # # Lung = annotates$tissPeakPropFC["prop. of AREs in Lung",],
                               # 
                               # col=annot2colorFun)

#word of bag and marker labels
ha.sampLabel = rowAnnotation(bar = samp2AREModGrp[rownames(ARES.modules.ordered)],
                              samp.nm=anno_mark(at = samps.nm.slcted_at, 
                                                side="right",
                                                labels = gsub("BSS\\d+_", "",samps.nm.slcted, perl=T), 
                                                labels_gp=gpar(fontsize=8)),
                              col = list(bar = ARE.ModGrp2COLOR))

#ha.genomic =HeatmapAnnotation(test=runif(125))
ht.EpiMap=Heatmap(ARES.modules.ordered, 
          col=COLs.MATRIX,
          cluster_rows=Epimap.samp.hclu.res,
          cluster_columns=F,
          #column_split = covar.value[colnames(pks.scaled.filt)],
          show_column_names = T,
          column_names_gp = gpar(fontsize = 6, col="black"),
          show_row_names=F,
          #column_names_gp = gpar(fontsize = 5, col="grey"),
          #left_annotation =ha.genomic,
          top_annotation=ha.AREprop,
          right_annotation=ha.sampLabel)

# 
# ht.eGTEx=Heatmap(annotates$eGTEx.modActivty, 
#            col=COLs.MATRIX,
#            cluster_rows=F, #Epimap.samp.hclu.res,
#            cluster_columns=F,
#            #column_split = covar.value[colnames(pks.scaled.filt)],
#            show_column_names = T,
#            column_names_gp = gpar(fontsize = 6, col="black"),
#            show_row_names=F)
           #column_names_gp = gpar(fontsize = 5, col="grey"),
           #left_annotation =ha.genomic,
           #top_annotation=ha.AREprop,
           #right_annotation=ha.sampLabel)

# lgd_list.line = list(Legend(labels = rownames(annotates$tissPeakPropFC),
#                        title = "mARE.FC", type = "points",
#                        pch = c(0,2,3,4),
#                        legend_gp = gpar(col = 2:5)))
# lgd_list.stackedbar = list(Legend(labels = rownames(annotates$tissPeakPropFC),
#                        title = "mARE.FC", type = "points",
#                        pch = 15,
#                        legend_gp = gpar(col = 2:5)))


ht_list = ht.annot.1 %v%  ht.EpiMap #%v% ht.eGTEx     #ht.annot.2 %v%

draw(ht_list)#, annotation_legend_list=lgd_list.stackedbar)


dev.off()


#

# 
file.ou.annotHeatmap.tsv = paste0(dir.ou, "modules.annotHeatmap.tsv")
cmb.mat=cbind(t(annotates$genomicRegions[, colnames(ARES.modules.ordered)]),
                 moduleSize= annotates$moduleSize[colnames(ARES.modules.ordered)],
                  t(ARES.modules.ordered))
write.table(cmb.mat, file=file.ou.annotHeatmap.tsv, sep="\t", quote=F)

