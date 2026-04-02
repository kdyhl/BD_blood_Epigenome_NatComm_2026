#V1.4
#c.getDiff_byLimma_V1.4_sva.R

#use sva to consider unwanted variables

args = commandArgs(trailingOnly=TRUE)
HM = args[1]

DRUG.N.CUTOFF=5
COUNT.CUTOFF=5
PRES.PROP = 0.2
TOP.NUM=1000

library(limma)
library(data.table)
#library("SmartSVA")
library(sva)
library(pheatmap)
# library(clusterProfiler)
# library(AnnotationHub)
# library(org.Hs.eg.db)
#

TYPE2COL=c("grey","red")
names(TYPE2COL)=c("control", "case") #c("IBDRA_Control", "CD", "UC", "RA")

file.in.meta.RData = paste0("a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.goodSamples.meta.RData")
file.in.norm.RData = paste0("a_2_peakHeightsNormalized_V1.4/", HM, ".MayoWithRef.heightsMean.depthGCNormed.RData")
file.in.ref2pkOvlprefSlct.combBP= paste0("a_2_peakHeightsNormalized_V1.4/", HM, ".ref2pkOvlpWithTest.refSlct.combBP.RDS")
file.in.covarite.csv = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
file.in.sampInfo.RData = "../sampleMeta/a_sampleInfo_lithiumAndTobacco_V1.1_researchMeta/samp.info.RData"



dir.ou=paste0("./c_bulkDiffPeak_limma_V1.4_sva_drugs/", HM, "/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.diffPeak.RData = paste(dir.ou, HM, ".diffPeak.RData", sep="")
file.ou.diffPeak.pdf = paste(dir.ou, HM, ".diffPeak.pdf", sep="")
file.ou.diffPeak.bed.pref = paste(dir.ou, HM, ".diffPeak.", sep="")
file.ou.diffPeak.BG.bed = paste(dir.ou, HM, ".diffPeak.BG.bed", sep="")



load(file.in.norm.RData)
#ref2pkOvlpWithTest.refSlct = get('ref2pkOvlpWithTest.refSlct', e1)
ref2pkOvlpWithTest.refSlct.combBP = readRDS(file.in.ref2pkOvlprefSlct.combBP)
#load(file.in.deconv.RData)
load(file.in.meta.RData)
load(file.in.sampInfo.RData)


covars = read.table(file.in.covarite.csv, sep=",", header=T, row.names=1)
rownames(covars)=paste("id", rownames(covars), sep="_")
#
pks.slct = rownames(ref2pkOvlpWithTest.refSlct.combBP)[apply(ref2pkOvlpWithTest.refSlct.combBP,1,sum)>=0.5]
samps.counts.disease = round(samp2HeightsMean$depth.GCnormTeng[pks.slct,unlist(dis2samps[c("control", "case")])])
pks.filt = pks.slct[apply(samps.counts.disease>=COUNT.CUTOFF, 1, sum)/ncol(samps.counts.disease)>=PRES.PROP]

meta.df = data.frame(disease=c(rep("control", length(dis2samps[["control"]])), 
                               rep("disease", length(dis2samps[["case"]]))),
                     covars[colnames(samps.counts.disease), c("Gender",  "age_atsamp")],
                     samp2info.all[colnames(samps.counts.disease), ] #c("samp.LiMed", "samp.vpa", "samp.tobacco", "samp.alcohol")
                     #"Ethnicity_Name", and race not informative
                     
)
# age = meta.goodsamples[colnames(samps.counts.disease), "Age"],
# sex = meta.goodsamples[colnames(samps.counts.disease), "Gender"])
current_drug_terms= colnames(meta.df)[grepl("_current$", colnames(meta.df))]
current_drug_terms.filt= current_drug_terms[apply(meta.df[, current_drug_terms],2,sum)>=DRUG.N.CUTOFF]

meta.filt.df = meta.df[,c("disease", "Gender", "age_atsamp", current_drug_terms.filt)] #"samp.tobacco", "samp.alcoholAbuse", 
  
design.matrix.mod0=model.matrix(~age_atsamp+Gender, meta.filt.df)
design.matrix.mod1=model.matrix(~., meta.filt.df)

counts.filt=samps.counts.disease[pks.filt, ]
# # #smartSVA to estimate the surrogate variates
# exp.r <- t(resid(lm(t(exp) ~., data=meta.df)))
# # ## Add one extra dimension to compensate potential loss of 1 degree of freedom
# # ## in confounded scenarios (very important)
# n.sv <- EstDimRMT(exp.r, FALSE)$dim + 1
# smartSVA.obj <- smartsva.cpp(as.matrix(exp), mod=design.matrix.mod, mod0=design.matrix.mod0, n.sv=n.sv)

pdf(file.ou.diffPeak.pdf, height=9, width=12)

#voom + sva 
lib.size = rep(1e6, nrow(design.matrix.mod1))
v <- voom(counts.filt, lib.size=lib.size, design =design.matrix.mod1, plot=TRUE)
n.sv = num.sv(v$E, design.matrix.mod1, method="leek")

if(n.sv !=0)
{
  sv.obj = sva(v$E, design.matrix.mod1, design.matrix.mod0, n.sv=n.sv)
  colnames(sv.obj$sv)=paste0("sv", 1:n.sv)
#DEGs
  PCC.svVSotherVar = apply(sv.obj$sv, 2,
  FUN=function(sv)
  {
    apply(samp2info.all[colnames(samps.counts.disease), c("samp.tobacco", "samp.alcoholAbuse")], 2,
    FUN=function(y)      
    {
      mod=lm(sv~y)
      summary(mod)$r.squared
    })
  })
  barplot(PCC.svVSotherVar, 
          ylab="R squared of variance of sv explained", 
          col=c("red", "blue"),
          beside=T)
  legend("topright", 
         legend=rownames(PCC.svVSotherVar),
         fill=c("red", "blue"))
  
  meta.filt.df= cbind(meta.filt.df,
               sv.obj$sv)
}
design.matrix.mod2=model.matrix(~., meta.filt.df)
v = voom(counts.filt, lib.size=lib.size, design =design.matrix.mod2, plot=TRUE)
fit <- lmFit(v, design.matrix.mod2)
fit <- eBayes(fit)

#

res=list()
for(nm in colnames(design.matrix.mod2)[-1])
{
  res[[nm]]=topTable(fit, coef=nm, adjust="BH", number=nrow(counts.filt))
  
}


#pvalue distribution
layout(matrix(1:9, ncol=3))
for(nm in names(res))
{
  hist(res[[nm]]$P.Value, xlab="p-value of DP", main=nm)  
}

#MA plot
layout(matrix(1:9, ncol=3))
for(nm in names(res))
{
  status=rep("peak", nrow(fit))
  status[res[[nm]][rownames(fit), "adj.P.Val"]<=0.2 & res[[nm]][rownames(fit), "t"]>0] = "DP.UP"
  status[res[[nm]][rownames(fit), "adj.P.Val"]<=0.2 & res[[nm]][rownames(fit), "t"]<0] = "DP.DW"
  
  values <- c("DP.UP","DP.DW")
  col <- c("red","blue")
  attr(status,"values") <- values
  attr(status,"col") <- col
  attr(status,"cex") <- 0.5
  plotMA(fit, coef=nm, 
         main=paste(HM, nm),
         status=status)
  
}

#disease vs drug
layout(matrix(1:9, ncol=3))

for(nm in setdiff(names(res), "diseasedisease"))
{
  smoothScatter(res$diseasedisease[pks.slct, "t"], 
                res[[nm]][pks.slct, "t"], 
                xlab="disease t", ylab=paste0(nm, " t"),
                main=paste0("PCC: ", signif(cor(res$diseasedisease[pks.slct, "t"], res[[nm]][pks.slct, "t"]), 2)))
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}

dev.off()


save(pks.filt, 
     samps.counts.disease, 
     design.matrix.mod2, 
     fit, 
     res,
     file=file.ou.diffPeak.RData)



pks.filt.bed= gsub(":|-", "\t", pks.filt)
write(pks.filt.bed, file=file.ou.diffPeak.BG.bed)

#FDR
for(cutoff in c(0.01, 0.05, 0.1, 0.2))
{
  for(type in names(res))
  {
    diffPk.up = rownames(res[[type]])[res[[type]]$adj.P.Val<=cutoff & res[[type]]$t >0]
    diffPk.up= diffPk.up[!is.na(diffPk.up)]
    write(gsub(":|-", "\t", diffPk.up), file=paste(file.ou.diffPeak.bed.pref, type, ".Q", cutoff, ".up.bed", sep=""))
    
    diffPk.dw = rownames(res[[type]])[res[[type]]$adj.P.Val<=cutoff & res[[type]]$t <0]
    diffPk.dw= diffPk.dw[!is.na(diffPk.dw)]
    write(gsub(":|-", "\t", diffPk.dw), file=paste(file.ou.diffPeak.bed.pref, type, ".Q", cutoff, ".dw.bed", sep=""))
  }
}

#top peaks
for(type in names(res))
{
  diffPk.up2adj.p = res[[type]]$adj.P.Val[res[[type]]$t >0]
  names(diffPk.up2adj.p) = rownames(res[[type]])[res[[type]]$t >0]
  diffPk.up= names(diffPk.up2adj.p)[order(diffPk.up2adj.p, decreasing = F)][1:TOP.NUM]
  write(gsub(":|-", "\t", diffPk.up), file=paste(file.ou.diffPeak.bed.pref, type, ".top", TOP.NUM, ".up.bed", sep=""))
  
  diffPk.dw2adj.p = res[[type]]$adj.P.Val[res[[type]]$t <0]
  names(diffPk.dw2adj.p) = rownames(res[[type]])[res[[type]]$t <0]
  diffPk.dw= names(diffPk.dw2adj.p)[order(diffPk.dw2adj.p, decreasing = F)][1:TOP.NUM]
  write(gsub(":|-", "\t", diffPk.dw), file=paste(file.ou.diffPeak.bed.pref, type, ".top", TOP.NUM, ".dw.bed", sep=""))
  
}
