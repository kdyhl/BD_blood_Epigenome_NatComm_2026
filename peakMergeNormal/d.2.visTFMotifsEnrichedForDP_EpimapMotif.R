
if(!grepl("R version 4.0", R.version$version.string))
{
  stop("use R-4.0")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

library(ComplexHeatmap)
library(circlize)

#PVAL.CUTOFF = 0.01
QVAL.CUTOFF=0.1
# FC.CUTOFF=1.2
# FC.SATURATE=4
OR.CUTOFF=1.2

BG= "randomBG" #"givenBG" 

DP.V="DP.V1.4.1"
TOP.ATCH.N=15

#TRAITS.RM =c("isBad.byMismatch", "MHSRC", "MHNRTHEUR")
#TRAITS.LEFT.ORDER = c("")


#TISSUEs=c("Brain", "Heart", "Muscle", "Lung")

# tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
#                   Lung = "Lung",
#                   Muscle = "Muscle_Skeletal",
#                   Heart = "Heart_Left_Ventricle")


dir.in.motif =  paste0("d_TFMotifForPeak_DP.V1.4.1_V1.1_EpimapMotif_0kb_", BG, "/") #"b_2_TFMotifForPeak_SexAgeDis_V1.1_EpimapMotif_0kb_givenBG/" #"b_2_TFMotifForPeak_SexAgeDis_V1.1_EpimapMotif/"
files.in.motifScore.txt = system(paste0("find ", dir.in.motif, "/*/*/knownResults.txt"), intern = T)
files.in.motifScore.cond= sapply(files.in.motifScore.txt, 
FUN=function(f.i)
{
  buf= unlist(strsplit(dirname(f.i), split="/", fixed = T))
  #paste(buf[length(buf)-1], buf[length(buf)], sep="-")
  buf[length(buf)-1]
})
cond2files.in.motifScore = split(files.in.motifScore.txt, files.in.motifScore.cond)
# 
# files.in.gene=paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss, ".v8.egenes.txt.gz")
# names(files.in.gene) = names(tiss2GTExTiss)

#file.in.motif2gene = "~/lhou.compbio/data/motif/Homer.motif2TF.2020.txt" 

file.in.motif2archetype = "~/lhou.compbio/data/Epimap/motifs/collated_motifs_archetype.txt"

dir.ou= dir.in.motif
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData = paste0(dir.ou, DP.V, ".motifEnrich.q", QVAL.CUTOFF, ".OR", OR.CUTOFF, ".RData")
file.ou.pdf = paste0(dir.ou, DP.V, ".motifEnrich.q", QVAL.CUTOFF, ".OR", OR.CUTOFF, ".pdf")
file.ou.txt = paste0(dir.ou, DP.V, ".motifEnrich.q", QVAL.CUTOFF, ".OR", OR.CUTOFF, ".txt")

#

buf=read.table(file.in.motif2archetype, sep="\t", header = T, row.names = NULL, comment.char = "#")
archetype2motif = split(buf[,"Motif"], buf[,"Cluster_ID"])
motif2archetype = buf[,"Cluster_ID"]
names(motif2archetype) = buf[,"Motif"]

motifs.all = c()
for(f.i in cond2files.in.motifScore[[1]])
{
  buf=read.table(f.i, sep="\t", head=T, row.names = 1, comment.char = "", quote="", check.names = F) 
  motifs.all=c(motifs.all, rownames(buf))
}
motifs.left=setdiff(motifs.all, names(motif2archetype))
motif2archetype[motifs.left]="unknown"
archetype2motif[["unknown"]] = motifs.left


# #genes
# tiss.genes = lapply(files.in.gene,
# FUN=function(f.i.g)
# {
#   read.table(f.i.g, sep="\t", header = T, row.names = 1, stringsAsFactors = F)$gene_name
# })

#need to refine
motif2TF = sapply(motifs.all,
FUN=function(m)
{
  toupper(unlist(strsplit(m, split="_"))[1])
}) 

#
motif2pval.mat =matrix(1, ncol=length(motifs.all), nrow=length(cond2files.in.motifScore))
motif2OR.mat =matrix(1, ncol=length(motifs.all), nrow=length(cond2files.in.motifScore))
colnames(motif2pval.mat) <- colnames(motif2OR.mat) <- motifs.all
rownames(motif2pval.mat) <- rownames(motif2OR.mat) <- names(cond2files.in.motifScore)
for(cond in names(cond2files.in.motifScore))
{
  for(f.i in cond2files.in.motifScore[[cond]])
  {
    buf=read.table(f.i, sep="\t", head=T, row.names = 1, comment.char = "", quote="", check.names = F)  
    motifs=rownames(buf)
    motif2pval.mat[cond, motifs] = buf[motifs, 2]
    #motif2OR.mat[cond, motifs] = as.numeric(gsub("%", "", buf[motifs,6]))/as.numeric(gsub("%", "", buf[motifs,8
    prop.fg = as.numeric(gsub("%", "", buf[motifs,6]))/100
    prop.bg = as.numeric(gsub("%", "", buf[motifs,8]))/100
    motif2OR.mat[cond, motifs] = prop.fg/(1-prop.fg)/(prop.bg/(1-prop.bg))
    
  }
}

motif2qval.mat = motif2pval.mat*length(archetype2motif)
#motif2qval.mat = t(apply(motif2pval.mat, 1, p.adjust, method="BH"))

#motifs.slct = apply(motif2pval.mat[, !grepl("Ubiquitious|Multi-tissue",colnames(motif2pval.mat))]!=1, 1, any)
#motifs.slct = motifs[motif2TF[motifs, 1] %in% tiss.genes[[tiss]]]

is.kept= (motif2qval.mat<=QVAL.CUTOFF  & motif2OR.mat >= OR.CUTOFF)
motifs.slct = motifs.all[apply(is.kept, 2, any)]
motifArch.slct = levels(factor(motif2archetype[motifs.slct]))
cond.slct = rownames(motif2qval.mat)[apply(is.kept, 1, any)]
#motifs.slct.order= order(apply(motif2OR.mat[motifs.slct, ], 1, which.max), decreasing = F)

motif2pval.mat.filt= motif2pval.mat
motif2pval.mat.filt[!is.kept]=1
motif2pval.mat.filt=motif2pval.mat.filt[, motifs.slct]
motif2pval.mat.filt1=motif2pval.mat.filt[cond.slct,]
#
motif2OR.mat.filt= motif2OR.mat
motif2OR.mat.filt[!is.kept]=1
motif2OR.mat.filt=motif2OR.mat.filt[, motifs.slct]
motif2OR.mat.filt1=motif2OR.mat.filt[cond.slct,]


arch2motif= split(motifs.slct, motif2archetype[motifs.slct])
archType.mat.colnames = sapply(arch2motif,
FUN=function(motifs)
{
  motifs.mat = cbind(motif2OR.mat.filt[,motifs])
  row.top = which.max(apply(motifs.mat, 1, mean))[1]
  motifs[which.max(motifs.mat[row.top,])[1]]
})

# motifArch2pval.mat.filt = sapply(arch2motif,
# FUN=function(motifs)
# {
#   apply(cbind(motif2pval.mat.filt[, motifs]), 1, min)
# })

motifArch2OR.mat.filt = sapply(arch2motif,
FUN=function(motifs)
{
  #motifs = intersect(motifs.slct, archetype2motif[[arch]])  
  apply(cbind(motif2OR.mat.filt[, motifs]), 1, max)
})

motifArch2OR.mat.filt1 = sapply(arch2motif,
FUN=function(motifs)
{
 #motifs = intersect(motifs.slct, archetype2motif[[arch]])  
  apply(cbind(motif2OR.mat.filt1[, motifs]), 1, max)
})

topNArch=apply(motifArch2OR.mat.filt1, 1, 
FUN=function(x)
{
  list(colnames(motifArch2OR.mat.filt1)[ x>=OR.CUTOFF & rank(1/x)<=TOP.ATCH.N ] )
})
topNArch=levels(factor(unlist(topNArch)))

# cond2hm = t(sapply(rownames(motif2pval.mat.filt),
# FUN=function(cond)
# {
#   buf = unlist(strsplit(cond, split=".", fixed = T))  
#   hm = buf[1]
#   direct=buf[length(buf)]
#   #trait = paste(buf[c(-1, -length(buf))], collapse = ".")
#   return(c(hm=hm, direct=direct))
# }))

save(motif2pval.mat,
     motif2OR.mat,
     motif2qval.mat,
     motifs.slct,
     motif2TF,
     motif2pval.mat.filt,
     motif2pval.mat.filt1,
     motif2OR.mat.filt,
     motif2OR.mat.filt1,
     archType.mat.colnames,
     #motifArch2pval.mat.filt,
     motifArch2OR.mat.filt,
     motifArch2OR.mat.filt1,
     #cond2hm,
     file=file.ou.RData)



pdf(file.ou.pdf, width=6, height=12)

# ht1 = Heatmap(-log10(motifArch2pval.mat.filt), 
#               #cluster_columns=T,
#               #rect_gp = gpar(col = "white", lwd = 2),
#               row_names_gp = gpar(fontsize = 6) #,
#               # cell_fun = function(j, i, x, y, width, height, fill) 
#               # {
#               #   if(motif2qval.mat[, motifs.slct][i, j]<=QVAL.CUTOFF)
#               #     grid.text("*", x, y, gp = gpar(fontsize = 6))
#               # }
# )
# 
# print(ht1)

mat.max = max(motifArch2OR.mat.filt1[motifArch2OR.mat.filt1!=Inf])
motifArch2OR.mat.filt1[motifArch2OR.mat.filt1==Inf]=mat.max

#arch
ht = Heatmap(log2(t(motifArch2OR.mat.filt1)), 
              col= colorRamp2(c(0, mat.max), c("white", "red")),
              border =T,
              #cluster_rows=F,
              #row_split= factor(trait2annot[cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "trait"], "name"], levels = trait2annot[,"name"]),
              #cluster_row_slices = FALSE, 
              #row_labels = paste(cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "tiss"], cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "direct"], sep="_"),
             # row_names_gp = gpar(fontsize = 8),
              #row_title_rot = 0,
              row_labels = motif2TF[archType.mat.colnames]
              #column_names_gp = gpar(fontsize = 8),
              #row_gap = unit(3, "mm"),
              )

print(ht)



motifArch2log2OR.mat.filt1.ordered=log2(t(motifArch2OR.mat.filt1))
rownames(motifArch2log2OR.mat.filt1.ordered)=motif2TF[archType.mat.colnames]
motifArch2log2OR.mat.filt1.ordered=motifArch2log2OR.mat.filt1.ordered[row_order(ht),column_order(ht)]
write.table(motifArch2log2OR.mat.filt1.ordered, file=file.ou.txt, sep="\t", quote=F)


#top arch
ht = Heatmap(log2(t(motifArch2OR.mat.filt1[, topNArch])), 
             col= colorRamp2(c(0, mat.max), c("white", "red")),
             border =T,
             #cluster_rows=F,
             #row_split= factor(trait2annot[cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "trait"], "name"], levels = trait2annot[,"name"]),
             #cluster_row_slices = FALSE, 
             #row_labels = paste(cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "tiss"], cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "direct"], sep="_"),
             # row_names_gp = gpar(fontsize = 8),
             #row_title_rot = 0,
             row_labels = motif2TF[archType.mat.colnames[topNArch]]
             #column_names_gp = gpar(fontsize = 8),
             #row_gap = unit(3, "mm"),
)

print(ht)



#motif
mat.max = max(motif2OR.mat.filt1[motif2OR.mat.filt1!=Inf])
motif2OR.mat.filt1[motif2OR.mat.filt1==Inf]=mat.max
ht=Heatmap(log2(t(motif2OR.mat.filt1)), 
            col= colorRamp2(c(0, mat.max), c("white", "red")),
            border =T,
            #cluster_rows=F,
            #row_split= factor(trait2annot[cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "trait"], "name"], levels = trait2annot[,"name"]),
            #cluster_row_slices = FALSE, 
            #row_labels = paste(cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "tiss"], cond2traitAndTiss[rownames(motifArch2OR.mat.filt1), "direct"], sep="_"),
            row_names_gp = gpar(fontsize = 6),
           # row_title_rot = 0,
            row_labels = motif2TF[colnames(motif2OR.mat.filt1)]
            #column_names_gp = gpar(fontsize = 8),
            #row_gap = unit(3, "mm")
            ) 
print(ht)


dev.off()





# 
# 
# 
# ht2 = Heatmap(log2(motif2OR.mat.filt), 
#               #cluster_columns=T,
#               #rect_gp = gpar(col = "white", lwd = 2),
#               row_names_gp = gpar(fontsize = 6),
#               cell_fun = function(j, i, x, y, width, height, fill) 
#               {
#                 if(motif2qval.mat[, motifs.slct][i, j]<=QVAL.CUTOFF)
#                   grid.text("*", x, y, gp = gpar(fontsize = 6))
#               }
# )
# 
# print(ht2)
# 
# ht2 = Heatmap(enrich.logFC.matrix.filt1[, names(ARE.ModGrp2COLOR)], 
#               col= colorRamp2(c(-2, 0, 2), c("deepskyblue", "white", "red")),
#               cluster_columns=F,
#               row_split= factor(trait2annot[cond2traitAndTiss[rownames(enrich.logFC.matrix.filt1), "trait"], "name"], levels = trait2annot[,"name"]),
#               cluster_row_slices = FALSE, 
#               border = TRUE,
#               row_gap = unit(3, "mm"),
#               row_title_rot = 0,
#               row_labels = paste(cond2traitAndTiss[rownames(enrich.logFC.matrix.filt1), "tiss"], cond2traitAndTiss[rownames(enrich.logFC.matrix.filt1), "direct"], sep="_"),
#               #rect_gp = gpar(col = "white", lwd = 0),
#               row_names_gp = gpar(fontsize = 8)
#               # cell_fun = function(j, i, x, y, width, height, fill) 
#               #           {
#               #             if(enrich.qval.matrix.filt[i, j]<=QVAL.CUTOFF)
#               #             grid.text("*", x, y, gp = gpar(fontsize = 6))
#               #           }
# )
# 
# print(ht2)
# 
# 
# 
# 
# dev.off()
# 
# 
