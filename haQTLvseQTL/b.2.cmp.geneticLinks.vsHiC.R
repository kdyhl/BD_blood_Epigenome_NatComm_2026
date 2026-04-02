#V1.1
#b.2.cmp.geneticLinks.vsHiC_V1.1_gARE.R

#only focusing on gAREs
#compare links with Hi-C
#also add permutation within distance bins

args = commandArgs(trailingOnly=TRUE)
#TISS = args[1]
HM = args[1]
linkType= "gene2Neighb" #=args[2] #"gene2hicReg"

# tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
#                   Lung = "Lung",
#                   Muscle = "Muscle_Skeletal",
#                   Heart = "Heart_Left_Ventricle",
#                   Blood="Whole_Blood")
# COND2HiCTISS = c(Brain = "Dorsolateral_Prefrontal_Cortex", #Hippocampus|
#                   Lung = "Lung",
#                   #Muscle = "Muscle_Skeletal",
#                   Heart = "Left_Ventricle")

#GWAS.PVAL.CUTOFF = 1e-6
# coloc.H4.ABF.CUTOFF=0.5
# MR.log10PVAL.CUTOFF = -log10(0.05)
# PGS.pcor.log10PVAL.CUTOFF=-log10(0.05)

set.seed(666)
PROMT.WIND = 2000
SCORE2COLOR=c(coloc.H4.ABF="red2",
              coloc.H4vsH3.ratio="orangered3",
              MR.qval="blueviolet",
              FMeQTL2ARE.DistanceInv="turquoise3"#,
              #gene2ARE.DistanceInv="grey"
              )
perm.N=50
DIS.BIN.NUM=50


#library(Sushi)
#library(igraph)
library(PRROC)


file.in.TSS = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.TSSUp0Down1bp.geneName.ENSG.bed"

file.in.geneticLinks.RData = paste0("b_summary_Links.coloc.MR.FMeQTL_V1.3_gARE/Blood_",  HM, ".linkFromColoc.MR.RData")
                                    
file.in.HiCLink.hg19.RData = "../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.link.RData"
file.in.HiCRegAnnot.RData = "../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.annot.RData"

file.in.pk.hg38to19.RDS= paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".pk.hg38tohg19.RDS")

dir.ou="b_2_cmp_geneticLinks.withHiC_V1.1_gARE/"
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.cmp.RData = paste(dir.ou, HM, ".cmpWithHiC_", linkType, ".RData", sep="")
file.ou.cmp.pdf = paste(dir.ou, HM, ".cmpWithHiC_", linkType, ".noPromoter.pdf", sep="")
#file.ou.haQTLDiffPks2GWAS.RData = paste(file.out.pref, ".haQTLDiffPks2GWAS.RData", sep="")
#file.ou.geneNeighb = paste(file.out.pref, ".txt", sep="")

#
load(file.in.geneticLinks.RData)
load(file.in.HiCLink.hg19.RData)
load(file.in.HiCRegAnnot.RData)
pk.hg38to19=readRDS(file.in.pk.hg38to19.RDS)

#
buf=read.table(file.in.TSS, sep="\t", header = F, row.names=NULL, stringsAsFactors = F)
gNameAndENSG= cbind(gName=buf[,5], ENSG=gsub("_\\d+", "", buf[,6], perl = T))
gName2ENSG=split(gNameAndENSG[,"ENSG"], gNameAndENSG[,"gName"])
gNames.multiHits= names(gName2ENSG)[sapply(gName2ENSG, FUN=function(x) length(levels(factor(x))))>1]
gNameAndENSG.filt=gNameAndENSG[!(gNameAndENSG[,"gName"] %in% gNames.multiHits), ]
ENSG2gName=sapply(split(gNameAndENSG.filt[,"gName"], gNameAndENSG.filt[,"ENSG"]), FUN=function(x) levels(factor(x))) 

ENSG2gName=ENSG2gName[sapply(ENSG2gName, length)==1]

#
if(linkType == "gene2hicReg")
{
  g2regs = gene2hicReg #1 step neigibhor
}else
{
  if(linkType == "gene2Neighb")
  {
    g2regs = gene2Neighb #2 step neighbor
  }
}
gene2pk.hg19.byHiC = lapply(g2regs,
FUN=function(hic.regs)
{
  pks.hg19=sapply(hic.regs,
  FUN=function(hic.reg)
  {
    hic.reg2HM.Pks[[HM]][[hic.reg]]
  })
  levels(factor(unlist(pks.hg19)))
})

save(gene2pk.hg19.byHiC,
     file=file.ou.cmp.RData)
######

scores.all.df=data.frame(scores.all)
buf=sapply(rownames(scores.all),
FUN=function(x)
{
  unlist(strsplit(x, split=";"))
})
scores.all.df$gene=buf[1,]
scores.all.df$pk.hg38=buf[2,]

scores.all.df$hic= sapply(1:nrow(scores.all.df),
FUN=function(i)
{
  if(i %% 1000==0) print(i)
  
  g=scores.all.df$gene[i]
  pk.hg38=scores.all.df$pk.hg38[i]
  gName=ENSG2gName[g]
  pk.hg19=pk.hg38to19[pk.hg38]
  
  return (gName %in% names(gene2pk.hg19.byHiC) &&
          pk.hg19 %in%  gene2pk.hg19.byHiC[[gName]] )
})


save(gene2pk.hg19.byHiC,
     ENSG2gName,
     scores.all.df,
     file=file.ou.cmp.RData)

#
f.o.tsv=paste0(dir.ou, HM, ".cmpWithHiC_", linkType, ".tsv")
write.table(scores.all.df, file=f.o.tsv, sep="\t", quote=F, row.names=F)


#permutation
pair2disRank=rank(1/scores.all.df$gene2ARE.DistanceIn)
pair2disBin=floor(pair2disRank/(length(pair2disRank)/DIS.BIN.NUM))
pair2disBin[pair2disBin==DIS.BIN.NUM]=DIS.BIN.NUM-1
disBin2pairs=split(rownames(scores.all.df), pair2disBin)
bin.pairs=do.call(c, disBin2pairs)
  
scores.all.perm.list=lapply(names(SCORE2COLOR),
FUN=function(nm)
{
  print(nm)
  sapply(1:perm.N,
  FUN=function(i)
  {
    score.perm.list=lapply(disBin2pairs,
    FUN=function(pairs)
    {
      s=scores.all.df[pairs, nm]
      sample(s, length(s), replace = F)
      
    })
    score.perm=do.call(c, score.perm.list)
    names(score.perm) = bin.pairs
    score.perm[rownames(scores.all.df)]
  })
})
names(scores.all.perm.list)=names(SCORE2COLOR)
######
#

classLab=rep("bg", nrow(scores.all.df))
names(classLab)=rownames(scores.all.df)
classLab[scores.all.df$hic] = "fg"



PRC.res.vsHiC = lapply(names(SCORE2COLOR),
FUN=function(nm)
{
  print(nm)
  samps.ovlp = rownames(scores.all.df)[!is.na(scores.all.df[, nm]) &
                scores.all.df$gene2ARE.DistanceInv <1/PROMT.WIND]
  print(length(samps.ovlp))
  
  scoreAndPerm=cbind(scores.all.df[samps.ovlp,nm],
                     scores.all.perm.list[[nm]][samps.ovlp,])
  apply(scoreAndPerm,2,
  FUN=function(scores)
  {
    prc = pr.curve(scores.class0 = scores[classLab[samps.ovlp]=="fg"],
                   scores.class1 = scores[classLab[samps.ovlp]=="bg"], curve=T)
    fg.prop = sum(classLab[samps.ovlp]=="fg")/length(samps.ovlp)  
    
    return(list(prc= prc,
                fg.prop = fg.prop))
    
  })

})
names(PRC.res.vsHiC) = names(SCORE2COLOR)

save(gene2pk.hg19.byHiC,
     ENSG2gName,
     scores.all.df,
     PRC.res.vsHiC,
     file=file.ou.cmp.RData)

#
pdf(file.ou.cmp.pdf, height=8, width = 8)
layout(matrix(1:4, ncol=2, byrow = T))

for(nm in names(PRC.res.vsHiC))
{
  plot(c(0,1), c(0,.5), 
       xlab="recall", 
       ylab= "precision",
       type="n",
       main=paste0("promoter captured Hi-C as Gold Standard Positive\n", nm))
  
  
  for(i in 2:length(PRC.res.vsHiC[[nm]]))
  {
    lines(PRC.res.vsHiC[[nm]][[i]]$prc$curve[,1:2], col="gray", lwd=.5, lty=1)
  }
  lines(PRC.res.vsHiC[[nm]][[1]]$prc$curve[,1:2], col=SCORE2COLOR[nm], lwd=2, lty=1)
  
  segments(x0=0,
           y0=PRC.res.vsHiC[[nm]][[1]]$fg.prop,
           x1=1,
           y1=PRC.res.vsHiC[[nm]][[1]]$fg.prop, col=SCORE2COLOR[nm], lwd=2, lty=2)
  
  aucs=sapply(PRC.res.vsHiC[[nm]], FUN=function(x) x$prc$auc.integral)
  bg= PRC.res.vsHiC[[nm]][[1]]$fg.prop
  hist(aucs[-1],col="gray", border="gray", 
       xlim=c(min(c(aucs,bg)), max(c(aucs,bg))),
       xlab="AUC")
  abline(v=aucs[1], col=SCORE2COLOR[nm], lwd=2, lty=1)
  abline(v=bg, col=SCORE2COLOR[nm], lwd=2, lty=2)
  # legend("topright", 
  #        legend=c(paste(1:length(PRC.res.vsHiC[[nm]]), "AUPRC:", signif(sapply(PRC.res.vsHiC[[nm]], FUN=function(x) x$prc$auc.integral),3)),
  #                 paste0("background: ", signif(PRC.res.vsHiC[[nm]][[1]]$fg.prop,3))),
  #        lwd=2,
  #        lty = c(1, rep(2, perm.N), 2),#rep(1:2, each= length(PRC.res.vsHiC)),
  #        col= c(SCORE2COLOR[nm], rep("gray", perm.N), SCORE2COLOR[nm])
  #       )
}

dev.off()



for(nm in names(PRC.res.vsHiC))
{
  f.o.tsv=paste(dir.ou, HM, ".", nm, ".cmpWithHiC_", linkType, ".noPromoter.PRC.tsv", sep="")  

  for(i in 1:length(PRC.res.vsHiC[[nm]]))
  {
    print(i)
    res=PRC.res.vsHiC[[nm]][[i]]$prc$curve
    if(i==1)
    {
      colnames(res)=paste0("observed_", c("recall", "precision", "cutoff"))
    }else
    {
      colnames(res)=paste0("perm_", i-1, "_", c("recall", "precision", "cutoff"))
    }

    res=t(res)

    write.table(res, f.o.tsv, sep="\t", quote=F, append=T)  
  }
  
  
}



