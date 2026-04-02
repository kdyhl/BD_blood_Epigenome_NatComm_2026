#compare H3K27ac haQTL with other HMs
#a1 and a2 are defined in vcf files which are reference allele and alternative allele when mapping QTLs

source("~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R")
library(data.table)
library(ggplot2)
options(scipen = 999)

args = commandArgs(trailingOnly=TRUE)
HM = args[1]

#CHRs = "chr19"

PEAK.WIN = 100000
NOM.P.CUTOFF=0.001

#
files.in.haQTL.BP.a = system(paste0("find ./a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, "_nominal/chr*sorted.bed.gz"), intern=T)
names(files.in.haQTL.BP.a)= sapply(files.in.haQTL.BP.a , 
FUN=function(f.i)
{
  unlist(strsplit(basename(f.i), split=".", fixed=T))[1]
})

files.in.haQTL.BP.b = system(paste0("find ./a_4_haQTL_FixedFactorNum_V1.2.3.1_DP.BG_100k_genotypePCs/", HM, "_nominal/chr*sorted.bed.gz"), intern=T)
names(files.in.haQTL.BP.b)= sapply(files.in.haQTL.BP.b , 
FUN=function(f.i)
{
  unlist(strsplit(basename(f.i), split=".", fixed=T))[1]
})


dir.ou= "c_comphQTL_amongVersions/"
dir.create(dir.ou, showWarnings = F, recursive = T)  

file.ou.pdf = paste(dir.ou, HM, "V1.2.3.vs.V1.2.3.1_", PEAK.WIN/1000, "k_filtbyP", NOM.P.CUTOFF, ".pdf", sep="")
file.ou.RData = paste(dir.ou, HM, "V1.2.3.vs.V1.2.3.1_", PEAK.WIN/1000, "k_filtbyP", NOM.P.CUTOFF, ".RData", sep="")
#


#
topHaQTL.btw.AandB.list=lapply(names(files.in.haQTL.BP.a),
FUN=function(chr)
{
  print(chr)
  df.a=fread(files.in.haQTL.BP.a[chr], sep="\t", head=T, data.table=F)
  df.b=fread(files.in.haQTL.BP.b[chr], sep="\t", head=T, data.table=F)
  colnames(df.b)= c("#chr", "start", "end", "rs", "a1", "a2","peak", "distance", "pval_b", "beta_b")
  #
  df.a.filt=df.a[df.a$pval<=NOM.P.CUTOFF,]
  topQTLs.list=lapply(split(df.a.filt, df.a.filt$peak),
  FUN=function(df)
  {
    df[which.min(df$pval)[1],]
  })
  topQTLs.df=do.call(rbind, topQTLs.list)

  #pks.a=pks.a[!duplicated(pks.a)]

  topQTLs.merge.df=merge(topQTLs.df, df.b[,c("peak", "rs", "pval_b", "beta_b")], 
                        by=c("peak", "rs"),
                        all.x=T)
})

topHaQTL.btw.AandB.df=do.call(rbind, topHaQTL.btw.AandB.list)
topHaQTL.btw.AandB.df$z=sign(topHaQTL.btw.AandB.df$beta) * qnorm(1 - topHaQTL.btw.AandB.df$pval/2)
topHaQTL.btw.AandB.df$z_b=sign(topHaQTL.btw.AandB.df$beta_b) * qnorm(1 - topHaQTL.btw.AandB.df$pval_b/2)

save(topHaQTL.btw.AandB.df, 
     file=file.ou.RData)
#load(file.ou.RData)

  
pdf(file.ou.pdf, height=6, width=6)

beta.cor=signif(cor(topHaQTL.btw.AandB.df$beta, topHaQTL.btw.AandB.df$beta_b),5)
ggplot(topHaQTL.btw.AandB.df, aes(beta, beta_b)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  geom_point(aes(color=log10(pval)), size=0.2)+
  ggtitle(paste0("pearson corr. ", beta.cor)) +
  labs(x = "effect size wo Top 5 genotype PCs", y = "effect size with Top 5 genotype PCs")+
  theme_bw()

topHaQTL.btw.AandB.filt.df=topHaQTL.btw.AandB.df[!is.infinite(topHaQTL.btw.AandB.df$z) & 
                                                  !is.infinite(topHaQTL.btw.AandB.df$z_b), ]
z.cor=signif(cor(topHaQTL.btw.AandB.filt.df$z, topHaQTL.btw.AandB.filt.df$z_b, use="pairwise.complete.obs" ),5)
ggplot(topHaQTL.btw.AandB.df, aes(z, z_b)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  geom_point(size=0.2)+
  ggtitle(paste0("pearson corr. ", z.cor)) +
  labs(x = "z score wo Top 5 genotype PCs", y = "z score with Top 5 genotype PCs")+
  theme_bw()

dev.off()
  





