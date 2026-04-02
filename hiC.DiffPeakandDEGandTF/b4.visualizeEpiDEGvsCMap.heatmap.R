#visualize the gseapy results



if(!grepl("R version 4.0", R.version$version.string))
{
  stop("use R-4.0")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

library(ComplexHeatmap)
library(circlize)
library(data.table)
#library(pheatmap)
library(RNOmni)

#pids.TFs=c(RORA="BRD-K20401833", FLI1="BRD-A62182663")

TFs.inferred=c("FLI1", "IRF8", "JDP2", "SPI1", "RFX4", "E2F4", "CTCF", "ATF1", "MBD2", "MECP2", "NRF1", "FEV", "EGR3", "CEBPB", "RFX2", "RORA", "ZBT14", "MYBL1", "NFIA")



dir.in="b_comp2CompoundSig_cMap_V1.1/"
file.in.slct.exp=paste0(dir.in, "b3_slctExp/exp.slct.p0.001.tsv")

file.in.gmt=paste0(dir.in, "EpiDEGs.gmt")
file.in.potentialTargets=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugTargets.txt")
file.in.drivers="../haQTLvseQTL/c_2_targetGenes_comb_V1.2_bulkAndcellsortedEQTL/215drivergenes.cellsorted.coloc.txt"
file.in.NES.RData=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.RData")

file.ou.RData=paste0(dir.in, "b3_slctExp/exp.p0.001.heatmap.RData")
file.ou.pdf=paste0(dir.in, "b3_slctExp/exp.p0.001.slct.DEGs.heatmap.pdf")
file.test.pdf=paste0(dir.in, "b3_slctExp/exp.p0.001.test.pdf")
file.ou.slct.DEGs=paste0(dir.in, "b3_slctExp/exp.p0.001.slct.DEGs.txt")
file.ou.slct.DEGs.id2name=paste0(dir.in, "b3_slctExp/exp.p0.001.slct.DEGs.id2name.txt") #from DAVID
file.ou.slct.DEGs.cluster=paste0(dir.in, "b3_slctExp/exp.p0.001.slct.DEGs.cluster") #from DAVID
#
load(file.in.NES.RData)
rownames(deg.gsea.res.df)=deg.gsea.res.df$experiment
#
exp.df= fread(file.in.slct.exp, sep="\t", data.table=F, header=T)
rownames(exp.df)=gsub("^\"\"|\"\"$", "", exp.df$rid, perl=T)
exp.mat=as.matrix(exp.df[,-1])

exp.mat.INT = apply(exp.mat,2, rankNorm)
#
buf= read.table(file.in.gmt, sep=",", row.names=NULL, header=F, stringsAsFactors=F)
gsets.list=lapply(buf[,1],
FUN=function(x)
{
    gs=unlist(strsplit(x, split="\t"))
    gs[-2]
})
names(gsets.list)=sapply(gsets.list, FUN=function(x)x[1])
gsets.list=lapply(gsets.list, FUN=function(x) x[-1])
gsets.ovlp.list=lapply(gsets.list, FUN=function(x) intersect(x, rownames(exp.df)))
gene2label= unlist(lapply(names(gsets.ovlp.list),
FUN=function(nm)
{
    x=rep(nm, length(gsets.ovlp.list[[nm]]))
    names(x)=gsets.ovlp.list[[nm]]
    return(x)
}))


#
pdf(file.test.pdf, width=20, height=12)
layout(matrix(1:15, ncol=3))
par(mar=c(3,3,3,3))
for(i in sample(1:ncol(exp.mat.INT), 100))
{
    density.all <- density(exp.mat.INT[,i])
    # density.dn <- density(exp.mat[gsets.ovlp.list$dnDEG,i])
    # density.up <- density(exp.mat[gsets.ovlp.list$upDEG,i])

    # # Plot the first density
    # plot(density.all, col = "red", lwd = 2, xlim = c(-10, 10), ylim = c(0, 0.5), 
    #       main = "", xlab = "Values", ylab = "Density")

    # # Add the second and third densities
    # lines(density.dn, col = "green", lwd = 2)
    # lines(density.up, col = "blue", lwd = 2)

    # hist(exp.mat[,i], breaks = 40, col = rgb(0.5, 0.5, 0.5, 0.1), freq = FALSE, xlim = c(-10, 10), ylim = c(0, 0.4),
    #  main = "Normalized Overlaid Histograms", xlab = "Values", ylab = "Density")

    # Add the second histogram
    hist(exp.mat.INT[gsets.ovlp.list$dnDEG,i], breaks = 40, col = rgb(0, 0, 1, 0.5), freq = FALSE,  
        xlim = c(-10, 10), ylim = c(0, 1),
        main="",  xlab = "Values", ylab = "Density")

    # Add the third histogram
    hist(exp.mat.INT[gsets.ovlp.list$upDEG,i], breaks = 40, col = rgb(1, 0, 0, 0.5), freq = FALSE, add = TRUE)

    lines(density.all, col = rgb(0, 0, 0, 1), lwd = 2)


    # Add a legend
    legend("topright", legend = c("all", "dn", "up"), 
           fill = c(rgb(.5, 0.5, 0.5, 0.5), rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))

}
dev.off()




# exp.mat[exp.mat>5]=5
# exp.mat[exp.mat< -5]=-5

# exp.z.mat= apply(exp.mat, 2,
# FUN=function(x)
# {
#     (x-mean(x))/sd(x)
# })
# exp.z.mat[exp.z.mat>5]=5
# exp.z.mat[exp.z.mat< -5]=-5
# #

#filter genes to find the robust differential genes
gs.ovlp=unlist(gsets.ovlp.list)
gs.ovlp.topDEG.list=apply(exp.mat.INT, 2, 
FUN=function(x)
{
    x.quant=quantile(x, prob=c(0.1, 0.9))
    top.upDEGs= names(x)[x>x.quant[2]]
    top.dnDEGs= names(x)[x<x.quant[1]]

    intersect(gs.ovlp, c(top.upDEGs, top.dnDEGs))
})
gs2topDEGCounts=summary(factor(unlist(gs.ovlp.topDEG.list)), maxsum=length(gs.ovlp))
gs.ovlp.topDEGs= names(gs2topDEGCounts)[rank(gs2topDEGCounts)>(length(gs2topDEGCounts)-499)]

write(gs.ovlp.topDEGs, file.ou.slct.DEGs)

#
buf=read.table(file.ou.slct.DEGs.id2name, sep="\t", header=T, row.names=1)
gid2name=buf[,1]
names(gid2name)=rownames(buf)
gname2id=rownames(buf)
names(gname2id)=buf[,1]

# #
# targets=read.table(file.in.potentialTargets)[,1]
# drivers=read.table(file.in.drivers, sep="\t", row.names=1, head=T)[,1]


# DEGs.ovlpTFDriver=intersect(gid2name, c(drivers, TFs.inferred))
# DEGs.ovlp.index=which(gs.ovlp.topDEGs %in% gname2id[DEGs.ovlpTFDriver])



#

row_dist <- dist(exp.mat.INT[gs.ovlp.topDEGs, ])
row_hclust <- hclust(row_dist)
# Cut the dendrogram into two clusters
row_clusters <- cutree(row_hclust, k = 2) 
row_split <- factor(row_clusters, levels = c(1, 2))
clu2genes=split(gid2name[names(row_split)], row_split)
for(clu in names(clu2genes))
{
    write(clu2genes[[clu]],file=paste0(file.ou.slct.DEGs.cluster, clu, ".txt"))
}
#
col_dist <- dist(t(exp.mat.INT[gs.ovlp.topDEGs, ]))
col_hclust <- hclust(col_dist)
# Cut the dendrogram into two clusters
col_clusters <- cutree(col_hclust, k = 2) 
col_split <- factor(col_clusters, levels = c(1, 2))
#
samp2type=rep("rest", ncol(exp.mat.INT))
names(samp2type)=colnames(exp.mat.INT)
samp2type[deg.gsea.res.df[colnames(exp.mat.INT), "upDEG.NES"]<0 & deg.gsea.res.df[colnames(exp.mat.INT), "dnDEG.NES"]>0]="reverseBDSignatures"

#
pdf(file.ou.pdf, width=6, height=6)

#pheatmap(exp.mat.INT[gs.ovlp,])
row_anno <- rowAnnotation(
  BD.DEGs = gene2label[gs.ovlp.topDEGs],
  col = list(BD.DEGs = c("dnDEG"="purple", "upDEG"="orange"))
  # mark = anno_mark( # Add the anno_mark annotation
  #   at = DEGs.ovlp.index,
  #   labels = gid2name[gs.ovlp.topDEGs[DEGs.ovlp.index]],
  #   labels_gp = gpar(fontsize = 5, col = "red"), # Custom font size and color
  #   link_width = unit(1, "mm") # Thickness of the arrow lines
  #)
)
col_anno <- columnAnnotation(
  SampleType = samp2type,  # Replace 'sample.types' with your vector
  col = list(SampleType = c("reverseBDSignatures" = "purple", "rest" = "orange")),
  annotation_name_side = "left"
)

ht = Heatmap(exp.mat.INT[gs.ovlp.topDEGs,], 
              col= colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              border =T,
              #cluster_rows=F,
              row_split= row_split,
              column_split=col_split,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              #cluster_row_slices = FALSE, 
              show_row_names=F,
              show_column_names=F,
              #row_labels = gid2name[gs.ovlp.topDEGs],
              #row_names_gp = gpar(fontsize = 3),
              #row_title_rot = 0,
              right_annotation = row_anno,
              top_annotation = col_anno
              #column_names_gp = gpar(fontsize = 8),
              #row_gap = unit(3, "mm"),
              )

print(ht)

dev.off()


save(exp.mat.INT,
    gene2label,
    gs.ovlp.topDEGs,
    row_hclust,
    row_split,
    col_hclust,
    col_split,
    samp2type,
    file=file.ou.RData)

#
f.o.tsv=paste0(dir.in, "b3_slctExp/exp.p0.001.slct.DEGs.heatmap.tsv")
write.table(exp.mat.INT[gs.ovlp.topDEGs,], f.o.tsv, sep="\t", quote=F)
