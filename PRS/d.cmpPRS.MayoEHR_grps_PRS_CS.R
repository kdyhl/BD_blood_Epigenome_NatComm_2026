#V1.2
#d.cmpPRS.MayoEHR_grps_V1.2_PRS_CS.R

R.LIBS.PATH = .libPaths()
.libPaths(c("~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/", R.LIBS.PATH))


options(scipen=999)


#compare PRS with EHR data and BD subgroups
library(dplyr)
library(ggplot2)
library(pheatmap)
library(PRROC)
library(RColorBrewer)
library(MOFA)


# CLUMP.LD.R2.CUTOFF=0.4
# CLUMP.DIS.CUTOFF=100000
PEAK.WIND=2000

#SNP.PVAL.CUTOFF="pval_0.05"
set.seed(666)

subGRP2COL = brewer.pal(10, "Set3")[1:5]
names(subGRP2COL)=paste0("subgrp", 1:5)


dir.in.prs=paste0("c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/PRS_wind", PEAK.WIND/1000, "kb/")
file.in.PRS.RData = paste0(dir.in.prs, "plink_PRS_CS_genomicRegions.RData")
file.in.snp.all="c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/BD_Mayo_pst_eff_a1_b0.5_phiauto_hg19.bed"
files.in.snp=dir(dir.in.prs, "snp.peaks.window.hg19.bed", full.names=T)
names(files.in.snp)=gsub(".snp.peaks.window.hg19.bed", "", basename(files.in.snp))
#file.in.WGSid2phenotype = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"

file.in.MOFA.RData="../sampleManifold/b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.RData"
file.in.subgrp.RData="../sampleManifold/b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.sampClu_V1.1.RData"
file.in.EHR.RData="../sampleManifold/c_patientClu2EHR_5HMs_byMOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.patientCluvsEHR.clu.corHclust5_V1.1.2.RData"
file.in.icd = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/icd_codes.csv"
file.in.clinic = "/broad/compbio/lhou/data/Mayo_Bipolar/210415_new_broad_data_deidentified.csv"
file.in.meta = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"

file.in.cov.RData="../sampleManifold/b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.FactorvsCovar.RData"

dir.ou=paste0("d_vis_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/PRS.wind", PEAK.WIND/1000, "kb/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.snpOvlp.pdf=paste0(dir.ou, "snpOvlp.pdf")
file.ou.PRS.boxplot.pdf = paste0(dir.ou, "plink.PRS.caseVsCtrl.pdf")
file.ou.PRS.corr.pdf = paste0(dir.ou, "plink.PRS.corr.pdf")
#file.ou.PRS.EHR.pdf = paste0(dir.ou, "plink.PRS.EHR.snp", SNP.PVAL.CUTOFF, ".pdf")
file.ou.prsVSEHR.pdf = paste0(dir.ou, "plink.PRSvsEHR.snp.pdf")
file.ou.PRS.vsSubtype.pdf = paste0(dir.ou, "plink.PRSvssubtype.pdf")
file.ou.PRS.vsMOFA.pdf = paste0(dir.ou, "plink.PRSvsFactorbyMOFA.pdf")
#file.ou.PRS.RData = paste0(dir.ou, "plink.PRSvsphenotype.RData")


#snp overlap
snps.list=lapply(files.in.snp,
FUN=function(f.i.snp)
{
  read.table(f.i.snp, stringsAsFactor=F)[,1]
})
snps.list$all=read.table(file.in.snp.all, sep="\t", header=F, row.names=NULL, stringsAsFactor=F)[,4]
snps.set=setdiff(names(snps.list), "all")


jaccard.mat=matrix(0, ncol=length(snps.set), nrow=length(snps.set))
or.mat=matrix(1, ncol=length(snps.set), nrow=length(snps.set))
p.mat=matrix(1, ncol=length(snps.set), nrow=length(snps.set))
rownames(jaccard.mat) <- rownames(or.mat) <- rownames(p.mat) <- snps.set
colnames(jaccard.mat) <- colnames(or.mat) <- colnames(p.mat) <- snps.set

for(s1 in snps.set)
{
  for(s2 in snps.set)
  {
    snps.ovlp=intersect(snps.list[[s1]], snps.list[[s2]])
    snps.s1.specf=setdiff(snps.list[[s1]], snps.list[[s2]])
    snps.s2.specf=setdiff(snps.list[[s2]], snps.list[[s1]])
    snps.union=union(snps.list[[s1]], snps.list[[s2]])
    snps.rest=setdiff(snps.list$all, snps.union)

    res=fisher.test(matrix(c(length(snps.ovlp), 
                            length(snps.s1.specf), 
                            length(snps.s2.specf), 
                            length(snps.rest)),
                            ncol=2),
                    alternative="greater")
    
    jaccard.mat[s1, s2]=signif(length(snps.ovlp)/length(snps.union),2)
    or.mat[s1, s2]=res$estimate
    p.mat[s1, s2]=res$p.val


  }

}

pdf(file.ou.snpOvlp.pdf, width=10, height=10)

or.mat[or.mat==Inf]=max(or.mat[or.mat!=Inf])
rownames(or.mat) <- colnames(or.mat) <- rownames(jaccard.mat) <- colnames(jaccard.mat) <- gsub("snps.filtBy.", "", colnames(or.mat))
pheatmap(log10(or.mat),
          cluster_row=F,
          cluster_col=F,
          display_numbers=round(jaccard.mat,1))


dev.off()


#meta
load(file.in.subgrp.RData)

load(file.in.PRS.RData)


samp2PRS.df=samp2PRS.df[, grepl("PRS.z", colnames(samp2PRS.df))]
samp2PRS.df=cbind(samp2PRS.df, samp2grp[rownames(samp2PRS.df), c("group", "clu.corHclust5")])
samp2PRS.df=samp2PRS.df[!is.na(samp2PRS.df$group),]

samp2PRS.df$clu.corHclust5=paste0("group",samp2PRS.df$clu.corHclust5)


#PRS case vs control #################################
prs.z.nm=names(samp2PRS.df)[grepl("PRS.z", names(samp2PRS.df))]
res=lapply(prs.z.nm,
FUN=function(nm)
{
  data.frame(prs.z=samp2PRS.df[, nm],
            prs.nm=nm,
            group=samp2PRS.df$group,
            subgroup=samp2PRS.df$clu.corHclust5)
})
samp2PRS.forBoxplot.df = do.call(rbind, res)
 
df.list=split(samp2PRS.forBoxplot.df, samp2PRS.forBoxplot.df$prs.nm)
t.ps=sapply(prs.z.nm,
FUN=function(nm)
{
  df.nm=df.list[[nm]]
  p=t.test(df.nm$prs.z[df.nm$group=="case"], df.nm$prs.z[df.nm$group=="control"])$p.val
})
PRS.caseVSctrl.p.df=data.frame(t.p.caseVSctrol=t.ps,
                #x=(1:length(prs.z.nm))*2-.5,
                y=max(samp2PRS.forBoxplot.df$prs.z)+0.1,
                prs.nm=prs.z.nm,
                stringsAsFactors=F)
PRS.caseVSctrl.p.df$label=""
PRS.caseVSctrl.p.df$label[PRS.caseVSctrl.p.df$t.p.caseVSctrol<=0.001]="*"
PRS.caseVSctrl.p.df$label[PRS.caseVSctrl.p.df$t.p.caseVSctrol<=0.0005]="**"
PRS.caseVSctrl.p.df$label[PRS.caseVSctrl.p.df$t.p.caseVSctrol<=0.0001]="***"

PRS.nms.caseVSctrl.p.filt=PRS.caseVSctrl.p.df[PRS.caseVSctrl.p.df$t.p.caseVSctrol<=0.001 ,"prs.nm"]

PRS.nms.slcted=union(PRS.nms.caseVSctrl.p.filt, 
                c("filtBy.Brain_neuron_eGTEx.PRS.z",
                "filtBy.blood_gAREs_BDColoc.PRS.z", 
                "filtBy.blood_gAREs_BD.MR.PRS.z"))






pdf(file.ou.PRS.boxplot.pdf, width=20, height=6)

# ggplot(df, aes_string(x="all.PRS.z", y=nm)) +
#       geom_point(aes(shape=group, color=isGrp1or2))+ #aes(col=group), 
ggplot(data = samp2PRS.forBoxplot.df, aes(x = prs.nm, y = prs.z)) +
  geom_boxplot(aes(fill = group)) +
  #geom_jitter(width = 0.2, position = position_jitter(width=0.2)) +
  scale_fill_manual(values = c(control="deepskyblue1", case="red1"))+
  theme_minimal()+
  ggtitle("PRS distinguish BD and control") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_text(data = PRS.caseVSctrl.p.df, aes(x = prs.nm, y = y, label = label, color="red", size=10), vjust = -0.5)+
  ylim(NA, max(samp2PRS.forBoxplot.df$prs.z)+0.5)
  #+
  #facet_wrap(~cutoff, nrow=2)


df=samp2PRS.forBoxplot.df

prs.nm2df.list=split(df,df$prs.nm)
test.p.list=lapply(names(prs.nm2df.list),
FUN=function(prs.nm)
{
  df.nm=prs.nm2df.list[[prs.nm]]
  subgrp2df.nm.list=split(df.nm, df.nm$subgroup)
  t.ps=sapply(subgrp2df.nm.list,
  FUN=function(df.nm.grp)
  {
    p=t.test(df.nm.grp$prs.z[df.nm.grp$group=="case"], df.nm.grp$prs.z[df.nm.grp$group=="control"])$p.val
  })
  p.df=data.frame(t.p.caseVSctrol=t.ps,
                #x=(1:length(prs.z.nm))*2-.5,
                y=max(df$prs.z)+0.1,
                prs.nm=prs.nm,
                subgroup=names(t.ps)
                )
})
test.p.df=do.call(rbind, test.p.list)
test.p.df$label=""
test.p.df$label[test.p.df$t.p.caseVSctrol<=0.01]="*"
test.p.df$label[test.p.df$t.p.caseVSctrol<=0.005]="**"
test.p.df$label[test.p.df$t.p.caseVSctrol<=0.001]="***"

p=ggplot(data = df, aes(x = subgroup, y = prs.z)) +
    geom_boxplot(aes(fill = group)) +
    #geom_jitter(width = 0.2, position = position_jitter(width=0.2)) +
    scale_fill_manual(values = c(control="deepskyblue1", case="red1"))+
    theme_minimal()+
    #ggtitle(cutoff) +
    #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_text(data = test.p.df, aes(x = subgroup, y = y, label = label, color="red", size=10))+
    ylim(NA, max(df$prs.z)+0.5)+
    facet_wrap(~prs.nm, ncol=6)
print(p)


dev.off()




#PRS All vs PRS_filtered #################################

pdf(file.ou.PRS.corr.pdf, width=15, height=15)


df=samp2PRS.df

prs.z.nm=names(df)[grepl("PRS.z", names(df))]
res=lapply(setdiff(prs.z.nm, "all.PRS.z"),
FUN=function(nm)
{
  data.frame(prs.z.all=df$all.PRS.z,
            prs.z.epi=df[[nm]],
            prs.nm=nm,
            #cutoff=prs.df$cutoff,
            group=df$group,
            subgroup=df$clu.corHclust5)
})
res.df=do.call(rbind, res)

#r2
r2 <- res.df %>%
group_by(prs.nm) %>%
do({
  model <- lm(prs.z.epi ~ prs.z.all, data = .)
  data.frame(r.squared = summary(model)$r.squared)
})
r2$label <- paste0("R2=", round(r2$r.squared, 2))
r2$x <- -2
r2$y <- 2.5

p=ggplot(res.df, aes(x=prs.z.all, y=prs.z.epi)) +
  geom_point(aes(color=group))+ #aes(col=group), shape=, 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_text(data = r2, aes(x = x, y = y, label = label), inherit.aes = FALSE) + # Add this line
  theme_minimal()+
  #ggtitle(cutoff)+
  facet_wrap(~prs.nm, ncol=6)
print(p)


#
cors=cor(df[,prs.z.nm])
pheatmap(cors, main="across all ind.")

#
cors.case=cor(df[df$group=="case",prs.z.nm])
pheatmap(cors.case, main="across patients")




dev.off()



#
samp.meta=read.table(file.in.meta, sep=",", header=T, row.names = 1, stringsAsFactors = F)
rownames(samp.meta) =paste0("id_", rownames(samp.meta))


#correlation with MOFA factors
load(file.in.MOFA.RData)

mofa.factors = MOFAobject@Expectations$Z

samps.ovlp= intersect(rownames(mofa.factors), rownames(samp2PRS.df))
samps.ovlp.BD=samps.ovlp[samp2grp[samps.ovlp, "group"]=="case"]
samps.ovlp.CTL=samps.ovlp[samp2grp[samps.ovlp, "group"]=="control"]
prs.z.nms= names(samp2PRS.df)[grepl(".PRS.z", names(samp2PRS.df))]
corrs.factorVsPRS= cor(mofa.factors[samps.ovlp,], 
                      as.matrix(samp2PRS.df[samps.ovlp, prs.z.nms]), 
                      use="pairwise.complete.obs")

corrs.factorVsPRS.BD= cor(mofa.factors[samps.ovlp.BD,], 
                      as.matrix(samp2PRS.df[samps.ovlp.BD, prs.z.nms]), 
                      use="pairwise.complete.obs")

corrs.factorVsPRS.CTL= cor(mofa.factors[samps.ovlp.CTL,], 
                      as.matrix(samp2PRS.df[samps.ovlp.CTL, prs.z.nms]), 
                      use="pairwise.complete.obs")

#
PRS.vs.epiFactor.t=matrix(NA, nrow=ncol(mofa.factors), ncol=length(prs.z.nms))
PRS.vs.epiFactor.mark =matrix("", nrow=ncol(mofa.factors), ncol=length(prs.z.nms))
rownames(PRS.vs.epiFactor.mark) <- rownames(PRS.vs.epiFactor.t) <- colnames(mofa.factors)
colnames(PRS.vs.epiFactor.mark) <- colnames(PRS.vs.epiFactor.t) <- prs.z.nms
PRS.vs.epiFactor.p=PRS.vs.epiFactor.t


for(prs.z.nm in prs.z.nms)
{
  x=samp2PRS.df[samps.ovlp, prs.z.nm]
  for(lf in colnames(mofa.factors))
  {
    y=mofa.factors[samps.ovlp, lf]
    df=data.frame(y=y,
                  x=x,
                  disease=samp2grp[samps.ovlp, "group"],
                  sex=samp.meta[samps.ovlp, "Gender"],
                  age=samp.meta[samps.ovlp, "age_atsamp"]
                  )

    mod=lm(y~., df)
    PRS.vs.epiFactor.t[lf, prs.z.nm]= summary(mod)$coef[2, "t value"]
    PRS.vs.epiFactor.p[lf, prs.z.nm]= summary(mod)$coef[2, "Pr(>|t|)"]
  }
}
PRS.vs.epiFactor.mark[PRS.vs.epiFactor.p<=0.05]="*"
PRS.vs.epiFactor.mark[PRS.vs.epiFactor.p<=0.01]="**"
PRS.vs.epiFactor.mark[PRS.vs.epiFactor.p<=0.005]="***"
#
#BD patients only
PRS.vs.epiFactor.BD.t=matrix(NA, nrow=ncol(mofa.factors), ncol=length(prs.z.nms))
PRS.vs.epiFactor.BD.mark =matrix("", nrow=ncol(mofa.factors), ncol=length(prs.z.nms))
rownames(PRS.vs.epiFactor.BD.mark) <- rownames(PRS.vs.epiFactor.BD.t) <- colnames(mofa.factors)
colnames(PRS.vs.epiFactor.BD.mark) <- colnames(PRS.vs.epiFactor.BD.t) <- prs.z.nms
PRS.vs.epiFactor.BD.p=PRS.vs.epiFactor.BD.t
for(prs.z.nm in prs.z.nms)
{
  x=samp2PRS.df[samps.ovlp.BD, prs.z.nm]
  for(lf in colnames(mofa.factors))
  {
    y=mofa.factors[samps.ovlp.BD, lf]
    df=data.frame(y=y,
                  x=x,
                  #disease=samp2grp[samps.ovlp, "group"],
                  sex=samp.meta[samps.ovlp.BD, "Gender"],
                  age=samp.meta[samps.ovlp.BD, "age_atsamp"]
                  )

    mod=lm(y~., df)
    PRS.vs.epiFactor.BD.t[lf, prs.z.nm]= summary(mod)$coef[2, "t value"]
    PRS.vs.epiFactor.BD.p[lf, prs.z.nm]= summary(mod)$coef[2, "Pr(>|t|)"]
  }
}
PRS.vs.epiFactor.BD.mark[PRS.vs.epiFactor.BD.p<=0.05]="*"
PRS.vs.epiFactor.BD.mark[PRS.vs.epiFactor.BD.p<=0.01]="**"
PRS.vs.epiFactor.BD.mark[PRS.vs.epiFactor.BD.p<=0.005]="***"
#
pdf(file.ou.PRS.vsMOFA.pdf)
pheatmap(corrs.factorVsPRS, cluster_cols=F, main="corrs. PRS vs LF ")
pheatmap(corrs.factorVsPRS.BD, cluster_cols=F, main="corrs. PRS vs LF in BD samples")

pheatmap(corrs.factorVsPRS.CTL, cluster_cols=F, main="corrs. PRS vs LF in control samples")


pheatmap(PRS.vs.epiFactor.t, 
        display_numbers=PRS.vs.epiFactor.mark,
        cluster_cols=F, 
        main="LF ~ PRS + age +sex +disease")

pheatmap(PRS.vs.epiFactor.BD.t, 
        display_numbers=PRS.vs.epiFactor.BD.mark,
        cluster_cols=F, 
        main="LF ~ PRS + age +sex in BD patients")

pheatmap(t(PRS.vs.epiFactor.BD.t[, PRS.nms.slcted]), 
        display_numbers=t(PRS.vs.epiFactor.BD.mark[, PRS.nms.slcted]),
        cluster_cols=F, 
        main="LF ~ PRS + age +sex in BD patients")
dev.off()




#predict BD with PRS and EHR
# samp.meta=read.table(file.in.meta, sep=",", header=T, row.names = 1, stringsAsFactors = F)
# rownames(samp.meta) =paste0("id_", rownames(samp.meta))

buf= read.table(file.in.icd, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
icd2samp= split(paste0("id_", buf$newsubjectid), buf$Dx_Desc)
terms.icd=names(icd2samp)
ICD.GROUPS=list(diabetes=c(include="diabetes", exclude="screen"),
                hypertension=c(include="hypertension", exclude="screen|TRANSIENT"),
                pain=c(include="Pain|PAIN", exclude=NA),
                THROMBOSIS=c(include="EMBOLISM|THROMBO", exclude=NA),
                #infection=c(include="infection|viral", exclude=NA),
                #inflamm=c(include="chronic|nflamma|itis ", exclude="noninflamma|ALCOHOLIC|Immunization|ACUTE|INFECT"),
                infectionInflam=c(include="infection|nflamma|itis ", exclude="noninflamma|ALCOHOLIC|Immunization|ACUTE"),
                #inflammation=c(include="inflamma|acute", exclude="noninflamma|ALCOHOLIC|"),
                obesity=c(include="obesity", exclude=NA),
                glucose=c(include="glucose", exclude=NA),
                abuse=c(include="abuse", exlucde="PARTNER|PHYSICAL|emotional"))
aggGrp2id=lapply(ICD.GROUPS,
FUN=function(grp)
{
  terms.isInclude = grepl(grp["include"],
                        terms.icd,
                        ignore.case = T)
  if(!is.na(grp["exclude"]))
  {
    
    terms.isExclude = grepl(grp["exclude"],
                            terms.icd,
                            ignore.case = T)
    terms.isInclude=terms.isInclude & (!terms.isExclude)
  }
  print(terms.icd[terms.isInclude])
  print("\n")
  ids=icd2samp[terms.icd[terms.isInclude]]
  levels(factor(unlist(ids)))
})
samps.icd = levels(factor(unlist(icd2samp)))
samps.ovlp=intersect(samps.icd, rownames(samp2PRS.df))
icdGrp.mat=matrix(0, ncol=length(aggGrp2id), nrow=length(samps.ovlp))
rownames(icdGrp.mat)=samps.ovlp
colnames(icdGrp.mat)=names(aggGrp2id)
for(icdGrp in names(aggGrp2id))
{
  icdGrp.mat[intersect(samps.ovlp, aggGrp2id[[icdGrp]]), icdGrp]=1
}



prs.slcted.df=samp2PRS.df[samps.ovlp, PRS.nms.slcted]
#colnames(prs.slcted.df)=names(PRS.nms.slcted)

prsAndicdGrp.df= data.frame(y=factor(samp2PRS.df[samps.ovlp, "group"],levels=c("control", "case")),
                            prs.slcted.df,
                            sex=samp.meta[samps.ovlp, "Gender"],
                            age=samp.meta[samps.ovlp, "age_atsamp"],
                            icdGrp.mat,
                            stringsAsFactors=F)




#check corr with PRS and EHR
load(file.in.cov.RData)

stat.prsVScovar.patient = list()
pval.prsVScovar.patient=list()

samps.case=rownames(covars.all$meta)[covars.all$meta$group=="case"]
# prs.nms=names(samp2PRS.df)[grepl("PRS.z", names(samp2PRS.df))]
prs.slcted.df=samp2PRS.df[, PRS.nms.slcted]

for(categ in c("H3K27ac_CellFract", "H3K36me3_CellFract", "H3K4me1_CellFract", "med2.filt", "icd", "clinic"))
{
  ehr= covars.all[[categ]]
  statAndP=sapply(colnames(ehr), 
  FUN=function(ehr.nm)
  {
    print(ehr.nm)
   
    samps.case.ovlp = intersect(samps.case, intersect(rownames(prs.slcted.df), rownames(ehr)))

    df=data.frame(y=ehr[samps.case.ovlp, ehr.nm],
                  covars.all$meta[samps.case.ovlp, c("Gender", "age_atsamp")],
                  prs.slcted.df[samps.case.ovlp,])

    
    df$y[df$y==""]=NA

    y.class=class(df$y)
    y.levels=length(levels(factor(df$y[!is.na(df$y)])))
    if(y.levels<=1 || (y.class=="character" && y.levels>2))
    {
      return(c(rep(NA, ncol(prs.slcted.df)),  #R2
               rep(NA, ncol(prs.slcted.df)))) #pval
    }
    
    statAndp= sapply(colnames(prs.slcted.df),
    FUN=function(prs.nm)
    {
      fml=as.formula(paste0("y~",prs.nm, "+Gender+age_atsamp"))
      if(y.levels==2)
      {
        df$y=factor(df$y)
        if(all(summary(factor(df$y))>=5))
        {
          mod.sum=summary(glm(fml, data=df, family = binomial))
          return(mod.sum$coef[2, c("z value", "Pr(>|z|)")])   #&& length(mod.sum$coef)==4*4
        }
        return(c(NA,NA))
        
      }
      
      mod.sum = tryCatch(
        summary(lm(fml, data=df)),
        error = function(e) 
        {
          print(e)
          return(c(NA,NA))
        }
      )
      if((!is.na(mod.sum)) && nrow(mod.sum$coef)==4)
      {
        return(mod.sum$coef[2, c("t value", "Pr(>|t|)")])
      }else
      {
        return(c(NA,NA))
      }
    })
    statAndp=c(statAndp[1,], statAndp[2,])
  })
  
  stat=statAndP[1:(nrow(statAndP)/2), ]
  pval=statAndP[(nrow(statAndP)/2+1):nrow(statAndP),]
  
  stat.prsVScovar.patient[[categ]] = stat[, !is.na(stat[1,])]
  pval.prsVScovar.patient[[categ]]=pval[, !is.na(stat[1,])]
}
for(categ in names(stat.prsVScovar.patient))
{
  print(categ)
  colnames(stat.prsVScovar.patient[[categ]])=paste(colnames(stat.prsVScovar.patient[[categ]]), categ, sep="_")
  colnames(pval.prsVScovar.patient[[categ]])=paste(colnames(pval.prsVScovar.patient[[categ]]), categ, sep="_")
}
stat.prsVScovar.patient.matrix= do.call(cbind, stat.prsVScovar.patient)
pval.prsVScovar.patient.matrix= do.call(cbind, pval.prsVScovar.patient)
qval.prsVScovar.patient.matrix= apply(pval.prsVScovar.patient.matrix, 2, p.adjust, method="BH")



pdf(file.ou.prsVSEHR.pdf, height=10, width=15)  
QVAL.CUTOFF=0.1

stat.prsVScovar.patient.matrix.filt= t(stat.prsVScovar.patient.matrix[, apply(qval.prsVScovar.patient.matrix, 2, min, na.rm=T)<=QVAL.CUTOFF])
stat.prsVScovar.patient.isSig.matrix.filt=t(qval.prsVScovar.patient.matrix[, rownames(stat.prsVScovar.patient.matrix.filt)])
display_numbers=matrix("",
                       nrow=nrow(stat.prsVScovar.patient.isSig.matrix.filt),
                       ncol=ncol(stat.prsVScovar.patient.isSig.matrix.filt))
display_numbers[stat.prsVScovar.patient.isSig.matrix.filt<=QVAL.CUTOFF]="*"
rownames(display_numbers)=rownames(stat.prsVScovar.patient.isSig.matrix.filt)
colnames(display_numbers)=colnames(stat.prsVScovar.patient.isSig.matrix.filt)

prsVScovar.p=pheatmap(t(stat.prsVScovar.patient.matrix.filt),
         display_numbers=t(display_numbers),
         cluster_rows=T,
         cluster_cols=T,
         main=paste0("patients only, q ", QVAL.CUTOFF))

# stat.prsVScovar.patient.matrix.filt.rowOrder= order(apply(stat.prsVScovar.patient.matrix.filt, 1, which.max), decreasing = F)
# stat.prsVScovar.patient.matrix.filt.colRowSorted=stat.prsVScovar.patient.matrix.filt[stat.prsVScovar.patient.matrix.filt.rowOrder,]
# display_numbers.colRowSorted=display_numbers[stat.prsVScovar.patient.matrix.filt.rowOrder,]
# pheatmap(t(stat.prsVScovar.patient.matrix.filt.colRowSorted),
#          display_numbers=t(display_numbers.colRowSorted),
#          cluster_rows=F,
#          cluster_cols=F,
#          main=paste0("sorted, patients only, q ", QVAL.CUTOFF))

dev.off()

prs.nm.sorted=prsVScovar.p$tree_row$labels[prsVScovar.p$tree_row$order]
#





#
#predict BD subtypes with PRS and EHR
clinic= read.table(file.in.clinic, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
rownames(clinic)=paste0("id_", clinic$newsubjectid)


patients.subGrp_ehr_prs.df=prsAndicdGrp.df[prsAndicdGrp.df$y=="case",]
pats.id=rownames(patients.subGrp_ehr_prs.df)
patients.subGrp_ehr_prs.df=data.frame(patients.subGrp_ehr_prs.df, 
                      subgroup=paste0("subgrp", samp2grp[pats.id, "clu.corHclust5"]),
                      clinic[pats.id,],
                      stringsAsFactors=F)
# #inflammation markers
# cors=cor(patients.subGrp_ehr_prs.df$PRS.immune, 
#         patients.subGrp_ehr_prs.df[, c("toxoigg_z", "cmvigg_z", "cmvigm_z", "crp_z")],
#         use="pairwise.complete")

# model.addt <- lm(crp_z~PRS.immune+sex+age+inflamm+infection+pain+abuse+cmvigg_z, data = patients.subGrp_ehr_prs.df)
# summary(model.addt)$coef
model.addt <- glm(infectionInflam~PRS.immune+PRS.neuron+PRS.muscle+PRS.multiTiss+PRS.imm.gAREs.coloc+sex+age, data = patients.subGrp_ehr_prs.df,family = binomial)
summary(model.addt)$coef

model.addt <- glm(infectionInflam~PRS.immune+PRS.neuron+PRS.muscle+PRS.multiTiss+PRS.imm.gAREs.MR+sex+age, data = patients.subGrp_ehr_prs.df,family = binomial)
summary(model.addt)$coef

subgrps=levels(factor(patients.subGrp_ehr_prs.df$subgroup))
subgrp2prop=summary(factor(patients.subGrp_ehr_prs.df$subgroup))/nrow(patients.subGrp_ehr_prs.df)

#models
prs.mods=prs.nm.sorted #PRS.nms.slcted
#prs.mods=names(PRS.nms.slcted)
cov=c("infectionInflam", "glucose", "hypertension", "pain", "abuse", "age", "sex")
#cov=c( "age", "sex")#"infectionInflam", "glucose", "hypertension", "pain", "abuse",
#riskFacts.withPRS=c("PRS.immune", "infectionInflam", "glucose", "hypertension", "pain", "abuse")
#fml=as.formula(paste0("y~", paste(riskFacts, collapse="+"), "+age+sex"))
# riskFacts=c("infectionInflam", "glucose", "hypertension", "pain", "abuse")
# fml=as.formula(paste0("y~", paste(riskFacts, collapse="+"), "+age+sex"))



subGrp.riskFact.coef.z.mat=matrix(NA, ncol=length(prs.mods), nrow=length(subgrps))
rownames(subGrp.riskFact.coef.z.mat)=subgrps
colnames(subGrp.riskFact.coef.z.mat)=prs.mods
subGrp.riskFact.mark.mat <- subGrp.riskFact.coef.p.mat <-subGrp.riskFact.coef.z.mat
subGrp.riskFact.modAUPRC.mark.mat <- subGrp.riskFact.modAUPRC.p.mat <- subGrp.riskFact.coef.z.mat


#pdf(file.ou.PRS.EHR.pdf, width=5, height=5)
pdf(file.ou.PRS.vsSubtype.pdf, width=5, height=5)
layout(matrix(1:4, ncol=2))


for(subgrp in subgrps)
{
  df=patients.subGrp_ehr_prs.df
  df$y=rep("rest", nrow(df))
  df$y[df$subgroup==subgrp]="subgrp"
  df$y=factor(df$y, levels=c("rest", "subgrp"))

  for(mod in prs.mods)
  {
    fml=as.formula(paste0("y~", paste(c(mod, cov), collapse="+")))
    model <- glm(fml, data = df, family = binomial)
    subGrp.riskFact.coef.z.mat[subgrp, mod]=summary(model)$coef[mod, "z value"]
    subGrp.riskFact.coef.p.mat[subgrp, mod]=summary(model)$coef[mod, "Pr(>|z|)"]  

    plot(c(0,1), c(0,1), type="n", xlab="recall", ylab="precision", main=paste0(subgrp, " by ", mod))
    auprcs=sapply(1:50,
    FUN=function(i)
    {
      set.seed(i)
      
      ids.train <- sample(1:nrow(df), round(nrow(df)*0.7))
      train_data <- df[ids.train, ]
      test_data <- df[-ids.train, ]
      
      # Fit logistic regression model
      model.addt <- glm(fml, data = train_data, family = binomial)
      
      # Predict probabilities on the test set
      predicted_scores.addt <- stats::predict(model.addt, newdata = test_data, type = "response")
      #predicted_scores.intr <- predict(model.intr, newdata = test_data, type = "response")

      # Calculate Precision-Recall curve metrics
      pr_curve.addt <- pr.curve(scores.class0 = predicted_scores.addt[test_data$y=="subgrp"], 
                            scores.class1 = predicted_scores.addt[test_data$y=="rest"],
                            curve=T)
      
      # pr_curve.intr <- pr.curve(scores.class0 = predicted_scores.intr[test_data$y=="case"], 
      #                       scores.class1 = predicted_scores.intr[test_data$y=="control"])
      
      lines(pr_curve.addt$curve[, 1:2], col="red", lwd=0.2)    
      return(pr_curve.addt$auc.integral)#, intr=pr_curve.intr$auc.integral)
    })
    abline(h=subgrp2prop[subgrp], lty=2)
    ttest.res=t.test(auprcs-subgrp2prop[subgrp], alternative="greater")
    subGrp.riskFact.modAUPRC.p.mat[subgrp, mod]=ttest.res$p.val
  }
}


subGrp.riskFact.mark.mat[subGrp.riskFact.coef.p.mat<=0.1]="*"
subGrp.riskFact.mark.mat[subGrp.riskFact.coef.p.mat<=0.05]="**"
subGrp.riskFact.mark.mat[is.na(subGrp.riskFact.mark.mat)]=""

subGrp.riskFact.modAUPRC.mark.mat[subGrp.riskFact.modAUPRC.p.mat<=0.1]="*"
subGrp.riskFact.modAUPRC.mark.mat[subGrp.riskFact.modAUPRC.p.mat<=0.05]="**"
subGrp.riskFact.modAUPRC.mark.mat[is.na(subGrp.riskFact.modAUPRC.mark.mat)]=""
#layout(matrix(1, ncol=1))
# ys=barplot(t(subGrp.riskFact.coef.z.mat[, prs.mods]),
#         beside=T,
#         horiz=T,
#         las=2,
#         cex.axis=.6)

#pdf(file.ou.PRS.EHR.pdf, width=5, height=5)
pheatmap(t(subGrp.riskFact.coef.z.mat),
        cluster_rows=F,
        cluster_cols=F,
        display_numbers=t(subGrp.riskFact.mark.mat),
        main="RPS z to predict subgroup, * pvalue")

dev.off()

