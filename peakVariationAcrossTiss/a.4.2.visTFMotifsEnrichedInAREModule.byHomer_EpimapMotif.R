#

if(!grepl("R version 4.0", R.version$version.string))
{
  stop("use R-4.0")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

library(ComplexHeatmap) #library(pheatmap)
library(circlize)


# colfunc <- colorRampPalette(c("white", "magenta3"))
# colors= colfunc(10)

QVAL.CUTOFF=0.05
# FC.CUTOFF=1.2
# FC.SATURATE=4
OR.CUTOFF=1.2

#
file.in.AREModule.RData="a_4_modules_annotations_mergedPk3actHMs_Epimap/mARE2ModuleGrpName.RData" #"a_3_AREModulesInEpimap_V1.4_mergedPk4Tiss/4TissMergedPeaks.inEpimap.H3K27acSignal.RData"

dir.in.motif = "a_4_TFMotifFor3actHMsAREModule.Epimap_byHomer_V1.1_EpimapMotif_0kb_randomBG/" #"a_4_TFMotifFor4TissAREModule.Epimap_byHomer_V1.1_EpimapMotif/"
files.in = system(paste0("find ", dir.in.motif, "*/*/knownResults.txt"), intern = T)
files.in.cond = sapply(files.in, 
FUN=function(f.i)
{
  buf= unlist(strsplit(dirname(f.i), split="/", fixed = T))
  #paste(buf[length(buf)-1], buf[length(buf)], sep="-")
  buf[length(buf)-1]
})
files.in.motifScore.txt = split(files.in, files.in.cond)

file.in.motif2archetype = "~/lhou.compbio/data/Epimap/motifs/collated_motifs_archetype.txt"

dir.ou=dir.in.motif

file.ou.RData = paste0(dir.ou, "motifEnrichInMod.qval", QVAL.CUTOFF, ".OR", OR.CUTOFF, ".RData")
file.ou.pdf = paste0(dir.ou, "motifEnrichInMod.qval", QVAL.CUTOFF, ".OR", OR.CUTOFF, ".pdf")


#order of module groups
load(file.in.AREModule.RData)
clu.sorted = rownames(ARE.modName.df)

#motif information
buf=read.table(file.in.motif2archetype, sep="\t", header = T, row.names = NULL, comment.char = "#", stringsAsFactors = F)
archetype2motif = split(buf[,"Motif"], buf[,"Cluster_ID"])
motif2archetype = buf[,"Cluster_ID"]
names(motif2archetype) = buf[,"Motif"]

#
motifs.all = c()
for(f.i in files.in.motifScore.txt[[1]])
{
  buf=read.table(f.i, sep="\t", head=T, row.names = 1, comment.char = "", quote="", check.names = F) 
  motifs.all=c(motifs.all, rownames(buf))
}

motif2archetype[setdiff(motifs.all, names(motif2archetype))]="unknown"


motif2TF = sapply(motifs.all,
FUN=function(m)
{
  toupper(unlist(strsplit(m, split="_"))[1])
}) 

#motif enrichment score
motif2pval.list <- motif2qval.list <- motif2OR.list  <- motif2FractDiff.list <- list()

for(type in c("promoter", "enhancer"))
{
  motif2pval.list[[type]] <- motif2OR.list[[type]] <- matrix(1, nrow=length(motifs.all), ncol=length(clu.sorted))
  motif2FractDiff.list[[type]] = matrix(0, nrow=length(motifs.all), ncol=length(clu.sorted))
  colnames(motif2pval.list[[type]]) <- colnames(motif2OR.list[[type]]) <- colnames(motif2FractDiff.list[[type]]) <- clu.sorted
  rownames(motif2pval.list[[type]]) <- rownames(motif2OR.list[[type]]) <- rownames(motif2FractDiff.list[[type]]) <- motifs.all

  for(clu in clu.sorted) 
  {
    
    cond = paste(clu, type, sep="_")
    
    for(f.i in files.in.motifScore.txt[[cond]])
    {
      buf=read.table(f.i, sep="\t", head=T, row.names = 1, comment.char = "", quote="", check.names = F)  
      motifs=rownames(buf)
      motif2pval.list[[type]][motifs, clu] = buf[motifs, 2]
      prop.fg = as.numeric(gsub("%", "", buf[motifs,6]))/100
      prop.bg = as.numeric(gsub("%", "", buf[motifs,8]))/100
      motif2OR.list[[type]][motifs, clu] = prop.fg/(1-prop.fg)/(prop.bg/(1-prop.bg))
      motif2FractDiff.list[[type]][motifs, clu] = prop.fg-prop.bg
    }
    
  }
  motif2qval.list[[type]] = motif2pval.list[[type]]*length(archetype2motif) #apply(motif2pval.list[[type]], 2, p.adjust, method="BH")
}



#filter
motif2pval.filt.list= motif2pval.list
motif2qval.filt.list= motif2qval.list
motif2OR.filt.list= motif2OR.list
motif2FractDiff.filt.list = motif2FractDiff.list
for(type in c("promoter", "enhancer"))
{
  #motif2pval.filt.list[[tiss]]
  
  
  is.filt=(motif2qval.list[[type]]<=QVAL.CUTOFF & motif2OR.list[[type]] > OR.CUTOFF)
  motif2pval.filt.list[[type]][!is.filt]=1
  motif2OR.filt.list[[type]][!is.filt]=1
  motif2FractDiff.filt.list[[type]][!is.filt]= 0
  
  motifs.slct = motifs.all[apply(is.filt, 1, any)]
  motif2pval.filt.list[[type]]= motif2pval.filt.list[[type]][motifs.slct,]
  motif2OR.filt.list[[type]]= motif2OR.filt.list[[type]][motifs.slct,]
  motif2FractDiff.filt.list[[type]]= motif2FractDiff.filt.list[[type]][motifs.slct,]
  # motifs.slct.order= order(apply(motif2pval.mat[motifs.slct, ], 1, which.max), decreasing = F)
  # #
  
}



save(motifs.all,
     motif2pval.list,
     motif2qval.list,
     motif2OR.list,
     motif2FractDiff.list,
     motif2TF,
     motif2pval.filt.list,
     motif2OR.filt.list,
     motif2FractDiff.filt.list,
     file=file.ou.RData)

pdf(file.ou.pdf, width=15, height=15)
par(mar=c(5,5,6,9))
# 
# for(tiss in TISSUEs)
# {
#   smoothScatter(-log10(motif2pval.list[[tiss]]), 
#                 log2(motif2OddsRatio.list[[tiss]])) 
# }

for(type in c("promoter", "enhancer"))
{
  mat  = log2(motif2OR.filt.list[[type]])
  mat.max = max(mat[mat!=Inf])
  mat[mat==Inf]=mat.max
  row.order = order(apply(mat, 1, which.max), decreasing = F)
  
  col_fun = colorRamp2(c(0, mat.max), c("grey80", "magenta3"))
  
  
  # 
  #######
  archType2motif= split(rownames(mat[row.order, ]), motif2archetype[rownames(mat[row.order, ])])
  #archType2motif = archType2motif[!is.na(names(archType2motif))]
  archType2OR.mat = t(sapply(archType2motif,
  FUN=function(motifs)
  {
    apply(rbind(mat[motifs,]), 2, mean)
  }))
  archType2OR.mat.rownames = sapply(archType2motif,
  FUN=function(motifs)
  {
    motifs.mat = rbind(mat[motifs,])
    column.top = which.max(apply(motifs.mat, 2, mean))[1]
    motifs[which.max(motifs.mat[,column.top])[1]]
  })
  
  archType2OR.rowOrder = order(apply(archType2OR.mat, 1, which.max), decreasing = F)
  
  x=Heatmap(archType2OR.mat[archType2OR.rowOrder, ], 
            name = type, 
            cluster_rows=F,
            cluster_columns=F,
            #row_split = 2,
            #row_split = factor(motif2archetype[rownames(mat[row.order, ])]),
            row_labels = motif2TF[archType2OR.mat.rownames[archType2OR.rowOrder]],
            #row_labels = pk2gene[pks],
            #top_annotation = column_ha,
            #right_annotation = mark_ha,
            col=col_fun,
            show_column_names=T) 
  print(x)
  
  
  #
  x=Heatmap(mat[row.order, ], 
            name = type, 
            cluster_rows=F,
            cluster_columns=F,
            #row_split = 2,
            row_split = factor(motif2archetype[rownames(mat[row.order, ])], levels=rownames(archType2OR.mat)[archType2OR.rowOrder]),
            row_labels = motif2TF[rownames(mat[row.order, ])],
            #row_labels = pk2gene[pks],
            #top_annotation = column_ha,
            #right_annotation = mark_ha,
            col=col_fun,
            show_column_names=T) 
  print(x)
}
dev.off()
# 
# for(type in c("promoter", "enhancer"))
# {
#   mat  = motif2FractDiff.filt.list[[type]]
#   row.order = order(apply(mat, 1, which.max), decreasing = F)
#   
#   
#   pheatmap(mat[row.order, ], 
#            cluster_cols=F,
#            cluster_rows=F,
#            labels_row = motif2TF[mat[row.order, ]],
#            fontsize_row = 2,
#            main=type)
#   
#   
# }



# for(type in c("promoter", "enhancer"))
# {
#   pheatmap(log10(motif2pval.list[[type]]),
#            cluster_cols=F,
#            cluster_rows=F,
#            labels_row = motif2TF[rownames(motif2pval.list[[type]])],
#            fontsize_row = 2,
#            main=type)
# }




