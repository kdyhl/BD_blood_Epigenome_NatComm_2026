#V1.1 
#b3.clusterPatientsBasedOnMOFA.factor.V1.1_allSamps.R
#all sample cluster
#try hcluster on the correlation matrix too


#use  PCs to cluster patient samples to understand their heterogeneity
#first find K nearest neighbor and then cluster 
#and compare the membership of patients across different marks
#use normalized mutual informaiton to 


R.LIBS.PATH = .libPaths()
.libPaths(c("~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/", R.LIBS.PATH))


# library(ggplot2)
# library(ggrepel)
library(MOFA)

library(igraph)
library(pheatmap)
library(ggplot2)
library(Rtsne)
library(entropy)
library(RColorBrewer)
library(ComplexHeatmap)
library(mclust)

set.seed(666)

COR.CUTOFF=0.4
NEAREST_K=10

colors = brewer.pal(10, "Set3")

#args = commandArgs(trailingOnly=TRUE)

#TYPE="HiC.DPs"

#file.in.factor.RData = paste0("b_factorAnalysis_across5HMs_MOFA/HiC.DPs/ind.MOFA.topPeak10000.varProp0.01.RData")
file.in.factor.RData =  paste0("b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.RData")#

file.in.meta.csv = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"



#dir.ou=paste0("b_factorAnalysis_acrossHMs_MOFA/", TYPE, "/")
#dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.RData = gsub(".RData", ".sampClu_V1.1.RData", file.in.factor.RData)
file.ou.pdf = gsub(".RData", ".sampClu_V1.1.pdf", file.in.factor.RData)
file.ou.diag.RDat = gsub(".RData", ".sampClu_V1.1.diag.RData", file.in.factor.RData)
file.ou.diag.pdf = gsub(".RData", ".sampClu_V1.1.diag.pdf", file.in.factor.RData)
file.ou.barplot.pdf = gsub(".RData", ".sampClu.barplot_V1.1.pdf", file.in.factor.RData)
file.ou.tsv = gsub(".RData", ".sampClu_V1.1.tsv", file.in.factor.RData)

# 
# file.ou.RData = paste0(dir.ou, "patientclu.RData")
# file.ou.clus.pdf = paste0(dir.ou, "patientHeter.pdf")
# file.ou.cluCmp.pdf = paste0(dir.ou, "patientHeter.cmpAcrossHMs.pdf")
# #
#

meta = read.table(file.in.meta.csv, sep=",", header=T, row.names=1)
rownames(meta)=paste("id", rownames(meta), sep="_")
samps.case = rownames(meta)[meta$group=="case"]





#cluster
#from correlation to cluster
#
getClusterFromCors_kNN_Louvain=function(cors, cor.cutoff, nearest_K, weighted=T)
{
  adj.matrix = cors
  adj.matrix[cors<=cor.cutoff] =0
  diag(adj.matrix)=0
  print(paste0("average #nearst neighbor: ", sum(adj.matrix!=0)/nrow(adj.matrix)))
  adj.matrix.isKNN = apply(adj.matrix, 2, 
  FUN=function(x)
  {
    w=rep(1, length(x))
    x.k = sort(x, decreasing = T)[nearest_K]
    w[x<x.k]=0
    
    return(w)
  })
  
  adj.matrix= adj.matrix *adj.matrix.isKNN
  adj.matrix=pmax(adj.matrix, t(adj.matrix))#make it symmetric
  
  net = graph_from_adjacency_matrix(adj.matrix, mode="undirected", weighted=T)
  res=cluster_louvain(net)
  
  membership(res)
}



#

load(file.in.factor.RData)  

mofa.fators = MOFAobject@Expectations$Z
samps.ovlp = intersect(rownames(meta), rownames(mofa.fators))
#samps.case.ovlp = intersect(samps.case, rownames(mofa.fators))
#samps.case.pcs = PCA.results$topCVPks$x[samps.case.ovlp, 1:PCs.TOP.N]
#
mofa.fators.impute= apply(mofa.fators,2,
FUN=function(x)
{
  x[is.na(x)] = mean(x, na.rm=T)
  return(x)
})
mofa.fators.impute.scaled = scale(mofa.fators.impute)
#

  
#cluster
samps.topPCs.cors = cor(t(mofa.fators.impute.scaled[samps.ovlp,]), use = "pairwise.complete.obs")
# samp2clu.louvain=getClusterFromCors_kNN_Louvain(samps.topPCs.cors[samps.ovlp, samps.ovlp],
#                                cor.cutoff=COR.CUTOFF,
#                                nearest_K=NEAREST_K)

cors.hclust = hclust(dist(samps.topPCs.cors))
samp2clu.corHclust.3 <- cutree(cors.hclust, k = 3)
samp2clu.corHclust.4 <- cutree(cors.hclust, k = 4)
samp2clu.corHclust.5 <- cutree(cors.hclust, k = 5)
samp2clu.corHclust.6 <- cutree(cors.hclust, k = 6)
samp2clu.corHclust.7 <- cutree(cors.hclust, k = 7)

# samp2clu.annot = rep("control", ncol(samps.topPCs.cors))
# names(samp2clu.annot) = colnames(samps.topPCs.cors)
# samp2clu.annot[samps.ovlp] = samp2clu
  
samp2grp = data.frame(group= factor(meta[samps.ovlp, c("group")]),
                      #clu.louvain = as.character(samp2clu.louvain[samps.ovlp]),
                      clu.corHclust3 = as.character(samp2clu.corHclust.3[samps.ovlp]),
                      clu.corHclust4 = as.character(samp2clu.corHclust.4[samps.ovlp]),
                      clu.corHclust5 = as.character(samp2clu.corHclust.5[samps.ovlp]),
                      clu.corHclust6 = as.character(samp2clu.corHclust.6[samps.ovlp]),
                      clu.corHclust7 = as.character(samp2clu.corHclust.7[samps.ovlp]))
rownames(samp2grp) = samps.ovlp






tsne = Rtsne(mofa.fators.impute.scaled[samps.ovlp,], perplexity=10)


grp2caseProp=lapply(colnames(samp2grp)[-1],
FUN=function(clu.nm)
{
  grp2ind = split(rownames(samp2grp), samp2grp[[clu.nm]])
  
  sapply(grp2ind,
  FUN=function(inds)
  {
    case.N= sum(samp2grp[inds, "group"]=="case")
    ctrl.N= sum(samp2grp[inds, "group"]=="control")
    
    return(c(prop=case.N/(case.N+ctrl.N), 
             N=case.N+ctrl.N) )
  })
}) 
names(grp2caseProp) =colnames(samp2grp)[-1]


save(mofa.fators,
     samps.ovlp,
     mofa.fators.impute.scaled,
     samps.topPCs.cors,
     # samp2clu.louvain,
     # samp2clu.corHclust,
     samp2grp,
     grp2caseProp,
     tsne,
     file=file.ou.RData)


write.table(samp2grp[,c("group","clu.corHclust5")], file=file.ou.tsv, sep="\t", quote=F)



#diagnosis about the robustness
#by subsampling

ARIs=t(sapply(1:50,
FUN=function(i)
{
  print(i)
  samps.ds=sample(samps.ovlp, round(0.9*length(samps.ovlp)), replace=F)

  samp2clu.ds=list()
  cors.hclust.ds = hclust(dist(samps.topPCs.cors[samps.ds, samps.ds]))
  for(k in 3:7)
  {
    samp2clu.ds[[paste0("clu.corHclust", k)]] <- cutree(cors.hclust.ds, k = k)  
  }
  
  ari=sapply(names(samp2clu.ds),
  FUN=function(nm)
  {
    adjustedRandIndex(samp2grp[samps.ds, nm], samp2clu.ds[[nm]][samps.ds])
  })
  names(ari)=names(samp2clu.ds)

  return(ari)
}))

pdf(file.ou.diag.pdf, height=4, width=6)
par(mar=c(8,4,2,2))
boxplot(ARIs, 
        las=2,
        ylab="adjusted rand index")

dev.off()

f.o.ARI.tsv=gsub(".RData", ".sampClu_V1.1.ARI.tsv", file.in.factor.RData)
write.table(ARIs, f.o.ARI.tsv, sep="\t", quote=F)


# if(!grepl("R version 4.0", R.version$version.string))
# {
#   stop("use R-4.0")
# }
# 
# R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
# R.LIBS.PATH = .libPaths()
# .libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

# load(file.ou.RData)


colors.louvain =colors[1:length(levels(factor(samp2grp$clu.louvain)))]
names(colors.louvain) = as.character(1:length(colors.louvain)) # paste0("clu", 1:length(colors.louvain))

colors.corHclust4 =colors[1:length(levels(factor(samp2grp$clu.corHclust4)))]
names(colors.corHclust4) =  as.character(1:length(colors.corHclust4)) #paste0("clu", 1:length(colors.corHclust4))#

colors.corHclust5 =colors[1:length(levels(factor(samp2grp$clu.corHclust5)))]
names(colors.corHclust5) =  as.character(1:length(colors.corHclust5)) #paste0("clu", 1:length(colors.corHclust4))#

colors.corHclust6 =colors[1:length(levels(factor(samp2grp$clu.corHclust6)))]
names(colors.corHclust6) = as.character(1:length(colors.corHclust6)) # paste0("clu", 1:length(colors.corHclust6))


col.list = list(clu.louvain=colors.louvain,
                clu.corHclust4=colors.corHclust4,
                clu.corHclust5=colors.corHclust5,
                clu.corHclust6=colors.corHclust6,
                group=c("case"="darkgoldenrod3", "control"="grey"))



pdf(file.ou.barplot.pdf, height=12, width=14)

layout(matrix(1:12, ncol=3))
for(clu.nm in names(grp2caseProp))
{
  #
  prop.mat= grp2caseProp[[clu.nm]]
  #chi test
  BDVSgrp.mat=rbind(round(prop.mat["N",] * prop.mat["prop",]),
                       round(prop.mat["N",] * (1-prop.mat["prop",]))) 
  rownames(BDVSgrp.mat)=c("case", "control")
  colnames(BDVSgrp.mat)=paste0("group", 1:ncol(BDVSgrp.mat))
  p.chi=chisq.test(BDVSgrp.mat)$p.value

  #prop.test
  #i.minProp=which.min(prop.mat["prop",])[1]
  avg.prop=sum(BDVSgrp.mat["case", ])/sum(BDVSgrp.mat)

  p.propTest=sapply(1:ncol(prop.mat),
  FUN=function(i)
  {
    #p=prop.test(BDVSgrp.mat["case", c(i, i.minProp)], prop.mat["N", c(i, i.minProp)], alternative="greater")$p.value
    p=prop.test(x=BDVSgrp.mat["case", i], n=prop.mat["N", i], p=avg.prop, alternative="two.sided")$p.value

  })
  p.propTest.sigLvl=rep("", length(p.propTest))
  p.propTest.sigLvl[p.propTest<0.1]="*"
  p.propTest.sigLvl[p.propTest<0.01]="**"
  p.propTest.sigLvl[p.propTest<0.001]="***"

  bp.x=barplot(prop.mat["prop",]*100, ylab="case % in each subgroup",
               names.arg=paste0("group", colnames(prop.mat)),
               col=col.list[[clu.nm]],
               main=paste0(clu.nm, " P=", p.chi))
  abline(h=avg.prop*100, lty=2)
  text(x=bp.x,
       y=prop.mat["prop",]/2*100,
       labels=prop.mat["N",])
  text(x=bp.x,
       y=prop.mat["prop",]*100,
       labels=p.propTest.sigLvl)


}


dev.off()

#
clu.nm="clu.corHclust5"
prop.mat= grp2caseProp[[clu.nm]]
f.o.tsv=gsub(".RData", paste0(".sampClu_V1.1.", clu.nm, ".barplot.tsv"), file.in.factor.RData)
write.table(prop.mat, f.o.tsv, sep="\t", quote=F)



pdf(file.ou.pdf, height=12, width=14)

pheatmap(samps.topPCs.cors,
         cluster_rows = cors.hclust,
         cluster_cols = cors.hclust,
         annotation_row=samp2grp,
         annotation_col=samp2grp,
         annotation_colors=col.list,
         fontsize_row=4,
         fontsize_col=4)

# #
# layout(matrix(1:12, ncol=3))
# for(clu.nm in names(grp2caseProp))
# {
#   #
#   prop.mat= grp2caseProp[[clu.nm]]
#   #chi test
#   BDVSgrp.mat=rbind(round(prop.mat["N",] * prop.mat["prop",]),
#                        round(prop.mat["N",] * (1-prop.mat["prop",]))) 
#   rownames(BDVSgrp.mat)=c("case", "control")
#   colnames(BDVSgrp.mat)=paste0("group", 1:ncol(BDVSgrp.mat))
#   p.chi=chisq.test(BDVSgrp.mat)$p.value

#   #prop.test
#   i.minProp=which.min(prop.mat["prop",])[1]
#   p.propTest=sapply(1:ncol(prop.mat),
#   FUN=function(i)
#   {
#     p=prop.test(BDVSgrp.mat["case", c(i, i.minProp)], prop.mat["N", c(i, i.minProp)], alternative="greater")$p.value

#   })
#   p.propTest.sigLvl=rep("", length(p.propTest))
#   p.propTest.sigLvl[p.propTest<0.05]="*"
#   p.propTest.sigLvl[p.propTest<0.01]="**"
#   p.propTest.sigLvl[p.propTest<0.001]="***"

#   bp.x=barplot(prop.mat["prop",]*100, ylab="case % in each subgroup",
#                names.arg=paste0("group", colnames(prop.mat)),
#                col=col.list[[clu.nm]],
#                main=paste0(clu.nm, " P=", p.chi))

#   text(x=bp.x,
#        y=prop.mat["prop",]/2*100,
#        labels=prop.mat["N",])
#   text(x=bp.x,
#        y=prop.mat["prop",]*100,
#        labels=p.propTest.sigLvl)
# }

# 
# ha=rowAnnotation(samp2grp,
#                      col=col.list)

Heatmap(mofa.fators.impute.scaled[samps.ovlp,], 
        split = paste0("grp", samp2grp$clu.louvain),
        name="clu.louvain")

Heatmap(mofa.fators.impute.scaled[samps.ovlp,], 
        split = paste0("grp", samp2grp$clu.corHclust4),
        name="clu.corHclust4")

Heatmap(mofa.fators.impute.scaled[samps.ovlp,], 
        split = paste0("grp", samp2grp$clu.corHclust5),
        name="clu.corHclust5")

Heatmap(mofa.fators.impute.scaled[samps.ovlp,], 
        split = paste0("grp", samp2grp$clu.corHclust6),
        name="clu.corHclust6")





#
samp.clu.df = data.frame(tsne$Y, 
                         clu=factor(paste0("grp", samp2grp[samps.ovlp, "clu.louvain"])),
                         group=samp2grp[samps.ovlp, "group"],
                         mofa.fators.impute.scaled[samps.ovlp,])
colnames(samp.clu.df)[1:2]=c("tSNE1", "tSNE2")
p=ggplot(samp.clu.df, aes(x=tSNE1, y=tSNE2)) + 
  geom_point(aes(col=clu, pch=group)) +
  scale_shape_manual(values=c(case=1, control=4))+
  theme_bw()+
  ggtitle("clu.louvain")
print(p)

#
samp.clu.df = data.frame(tsne$Y, 
                         clu=factor(paste0("grp", samp2grp[samps.ovlp, "clu.corHclust4"])),
                         group=samp2grp[samps.ovlp, "group"],
                         mofa.fators.impute.scaled[samps.ovlp,])
colnames(samp.clu.df)[1:2]=c("tSNE1", "tSNE2")
p=ggplot(samp.clu.df, aes(x=tSNE1, y=tSNE2)) + 
  geom_point(aes(col=clu, pch=group)) +
  scale_shape_manual(values=c(case=1, control=4))+
  theme_bw()+
  ggtitle("clu.corHclust4")
print(p)

samp.clu.df = data.frame(tsne$Y, 
                         clu=factor(paste0("grp", samp2grp[samps.ovlp, "clu.corHclust5"])),
                         group=samp2grp[samps.ovlp, "group"],
                         mofa.fators.impute.scaled[samps.ovlp,])
colnames(samp.clu.df)[1:2]=c("tSNE1", "tSNE2")
p=ggplot(samp.clu.df, aes(x=tSNE1, y=tSNE2)) + 
  geom_point(aes(col=clu, pch=group)) +
  scale_shape_manual(values=c(case=1, control=4))+
  theme_bw() +
  ggtitle("clu.corHclust5")
print(p)

#
samp.clu.df = data.frame(tsne$Y, 
                         clu=factor(paste0("grp", samp2grp[samps.ovlp, "clu.corHclust6"])),
                         group=samp2grp[samps.ovlp, "group"],
                         mofa.fators.impute.scaled[samps.ovlp,])
colnames(samp.clu.df)[1:2]=c("tSNE1", "tSNE2")
p=ggplot(samp.clu.df, aes(x=tSNE1, y=tSNE2)) + 
  geom_point(aes(col=clu, pch=group)) +
  scale_shape_manual(values=c(case=1, control=4))+
  theme_bw() +
  ggtitle("clu.corHclust6")
print(p)


# for(lf in colnames(mofa.fators.impute.scaled))
# {
#   p=ggplot(samp.clu.df, aes(x=tSNE1, y=tSNE2)) + 
#     geom_point(aes(col=lf)) +
#     theme_bw()
#   print(p)
# }



dev.off()


#
f.o.heatmap.tsv = gsub(".RData", ".sampClu_V1.1.heatmap.tsv", file.in.factor.RData)
write.table(samps.topPCs.cors, f.o.heatmap.tsv, sep="\t", quote=F)


