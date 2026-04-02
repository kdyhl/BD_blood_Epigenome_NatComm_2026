
#based on XUshen's dataframe
#estimate pvalue of directionality consistency with null hypothesis that the sign is 50%, 50% binomial distribution
library(reshape2)
library(pheatmap)

file.in.df.RDS="b.QTLSharingAcrossHMs_byXushen/QTL.comparison.between.histone.rds"

file.ou.directConsist.RData="b.QTLSharingAcrossHMs_byXushen/QTL.comparison.between.histone.directConsist_byLeiH.RData"
file.ou.directConsist.pdf="b.QTLSharingAcrossHMs_byXushen/QTL.comparison.between.histone.directConsist_byLeiH.pdf"


HMs=c("H3K27ac", "H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")

QTL.acrossHM.df=readRDS(file.in.df.RDS)
HM.discAndValid=paste(QTL.acrossHM.df$mark.disc, QTL.acrossHM.df$mark.valid, sep=";") 
HMComp2QTLEffect.list=split(QTL.acrossHM.df, HM.discAndValid)


directConsist.propAndPval=sapply(HMComp2QTLEffect.list,
FUN=function(df)
{
  N=nrow(df)
  N.sameSign= sum(sign(df$beta.disc)==sign(df$beta.valid))
  N.oppoSign= sum(sign(df$beta.disc)!=sign(df$beta.valid))
  
  p.biasedToSame = 1-pbinom(N.sameSign, size=N, prob=0.5)
  p.biasedToOppo = 1-pbinom(N-N.sameSign, size=N, prob=0.5)
  sharedProp.sameSign=(2*N.sameSign-N)/N
  sharedProp.OppoSign=(2*N.oppoSign-N)/N
  
  return(c(N=N,
           N.sameSign=N.sameSign,
           N.oppoSign=N.oppoSign,
           p.biasedToSame=p.biasedToSame,
           p.biasedToOppo=p.biasedToOppo,
           sharedProp.sameSign=sharedProp.sameSign,
           sharedProp.OppoSign=sharedProp.OppoSign
           ))
})

comp.discHM=sub("(.*);.*", "\\1", colnames(directConsist.propAndPval))
comp.validHM=sub(".*;(.*)", "\\1", colnames(directConsist.propAndPval))



directConsist.pval.df= data.frame(discHM=factor(comp.discHM, levels=HMs),
                                  validHM=factor(comp.validHM, levels=HMs),
                                  pval=directConsist.propAndPval["p.biasedToSame",]
                                  #stringsAsFactors = F
                                  )
#H3K27me3 expected to be opposite with each other
directConsist.pval.df$pval[grepl("H3K27me3", colnames(directConsist.propAndPval))]=1-directConsist.pval.df$pval[grepl("H3K27me3", colnames(directConsist.propAndPval))]
directConsist.pval.mat=acast(directConsist.pval.df, validHM~discHM)
diag(directConsist.pval.mat)=0

directConsist.prop.df= data.frame(discHM=factor(comp.discHM, levels=HMs),
                                  validHM=factor(comp.validHM, levels=HMs),
                                  sharedProp=directConsist.propAndPval["sharedProp.sameSign",]
                                  #stringsAsFactors = F
)
#H3K27me3 expected to be opposite with each other
directConsist.prop.df$sharedProp[grepl("H3K27me3", colnames(directConsist.propAndPval))]=-directConsist.prop.df$sharedProp[grepl("H3K27me3", colnames(directConsist.propAndPval))]
directConsist.prop.mat=acast(directConsist.prop.df, validHM~discHM)
directConsist.prop.mat[directConsist.pval.mat>=0.1]=0
diag(directConsist.prop.mat)=1

save(directConsist.prop.mat,
     #directConsist.prop.mat,
     file=file.ou.directConsist.RData
     )

pdf(file.ou.directConsist.pdf)
pheatmap(directConsist.prop.mat,
         display_numbers=signif(directConsist.prop.mat,2))
dev.off()


#
f.o.tsv="b.QTLSharingAcrossHMs_byXushen/QTL.comparison.between.histone.directConsist_byLeiH.heatmap.tsv"
write.table(directConsist.prop.mat, file=f.o.tsv, sep="\t", quote=F)


