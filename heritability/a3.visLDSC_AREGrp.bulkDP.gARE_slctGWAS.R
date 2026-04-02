
#combine previous result 
#remove clusters that are empty
#reorder the traits that shows up


#for the baseline they get different number of SNPs as I got.
#need to use the overlapping snps, ignore possible change to .M, .M_5_50
#should guarantee all the new annotations share the same set of SNPs




library(pheatmap)
library("RColorBrewer")
library(corrplot)
#

dir.in.LDSC = "a_2_heritPartition_byLDSC_AREGrp.bulkDP.gCRE_slctGWAS/"
files.in.GWAS=c()

files.in.GWAS["T1D_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1222_self_reported_type_1_diabetes.bed"
files.in.GWAS["COVID19.HG"] = "~/lhou.compbio/data/GWAS/COVID19_hg/COVID19_HGI_ANA_C2_V2.sumstats.gz"
files.in.GWAS["IBD_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1461_self_reported_inflammatory_bowel_disease.bed"
files.in.GWAS["RA_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1464_self_reported_rheumatoid_arthritis.bed"
files.in.GWAS["MS_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1261_self_reported_multiple_sclerosis.bed"
files.in.GWAS["Neurot_UKB"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20127_Neuroticism_score.bed"
files.in.GWAS["AD_NG2022"]="~/lhou.compbio/data/GWAS/AD/AD_Bellenguez_NG2022/GCST90027158_buildGRCh38.sumstats.gz"
files.in.GWAS["PD_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1262_self_reported_parkinsons_disease.bed"
files.in.GWAS["PGC_SCZ"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_pgc.scz2.bed"
files.in.GWAS["PGC_BP"]="~/lhou.compbio/data/GWAS/pgc_bipolar/daner_PGC_BIP32b_mds7a_0416a.hg38.sumstats.gz"
files.in.GWAS["MDD_PGC"]= "~/lhou.compbio/data/GWAS/PGC/PGC_UKB_MDD_2019_forLDSC.sumstats.gz"
files.in.GWAS["crossDisorder_PGC"]= "~/lhou.compbio/data/GWAS/PGC/pgc_cdg2_meta_no23andMe_oct2019_v2.daner.forLDSC.sumstats.gz"
files.in.GWAS["HighBloodPressure"]="~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/highBloodPressure_6150_4_logistic.EUR.sumstats.gz"
files.in.GWAS["SmokeInit"]="~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/SmokingInitiation.txt.hg38.sumstats.gz"
files.in.GWAS["DrinkPerWeek"]="~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/DrinksPerWeek.txt.hg38.sumstats.gz"
files.in.GWAS["T2D_UKBS"]= "~/lhou.compbio/data/GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_UKB_20002_1223_self_reported_type_2_diabetes.bed"
files.in.GWAS["height_GIANT"]= "~/lhou.compbio/data/GWAS/GTEx.v8.GWAS/GTEx.v8.imputedGWAS.withRSID/imputed_GIANT_HEIGHT.bed.gz"
files.in.GWAS["CRP"]= "~/lhou.compbio/data/GWAS/GWAS.Catalog/CRP_GCST90029070_buildGRCh37.pval_floored.forLDSC.sumstats.gz"
files.in.GWAS["neutrophil"]= "~/lhou.compbio/data/GWAS/GWAS.Catalog/neutropCounts_GCST90691823.h.forLDSC.sumstats.gz"




dirs.in.annot = c("a_annotations_bulkDP/limma_V1.4.1/", "a_annotations_gARE_emp.p.fdr.cutoff0.05/", "a_annotations_3actHM_AREGrp/")

annot.nms=unlist(lapply(dirs.in.annot,
FUN=function(d.i)
{
  levels(factor(sapply(dir(d.i, paste0("annot.gz"), full.names=F), #"haQTL.*\\.1.annot.gz"  "diffPeak.*\\.1.annot.gz""
  FUN=function(f.i)
  {
   return(gsub(".\\d+.annot.gz", "", f.i, perl=T))
  })))
}
)) #"haQTL.*\\.1.annot.gz"  "diffPeak.*\\.1.annot.g)


dir.ou = paste0("a_3_visHeritPartition_AREGrp.bulkDP.gARE_slctGWAS/")
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.RData = paste0(dir.ou, "enrichAndpval.RData")
file.ou.pval.pdf = paste0(dir.ou, "enrichAndpval.pdf") 
file.ou.enrich.pdf = paste0(dir.ou, "enrichment.pdf") 



#LDSC --h2
ldsc.pval = matrix(1, nrow=length(files.in.GWAS), ncol=length(annot.nms))
rownames(ldsc.pval) = names(files.in.GWAS)
colnames(ldsc.pval) = annot.nms
ldsc.enrich =matrix(0, nrow=length(files.in.GWAS), ncol=length(annot.nms))
rownames(ldsc.enrich) = names(files.in.GWAS)
colnames(ldsc.enrich) = annot.nms


for(annot.nm in annot.nms)
{
  for(pheno in names(files.in.GWAS))
  {
    f.i.res = paste0(dir.in.LDSC, pheno,"_", annot.nm, ".results")
    if(file.exists(f.i.res))  
    {
      d = as.matrix(read.table(f.i.res, sep="\t", row.names=1, header=T))
    
      ldsc.enrich[pheno, annot.nm] = d["L2_1", "Enrichment"]
      ldsc.pval[pheno, annot.nm] = d["L2_1", "Enrichment_p"]
    }
  }
}
# colnames(ldsc.enrich) = gsub("PASS_", "", colnames(ldsc.enrich))
# colnames(ldsc.pval) = gsub("PASS_", "", colnames(ldsc.pval))
ldsc.pval[ldsc.enrich<0]=1
ldsc.enrich[ldsc.enrich<0]=0
#gwas.kept = apply(ldsc.enrich, 2, FUN=function(x){any(x!=0)}) 

ldsc.qval = t(apply(ldsc.pval, 1, FUN=function(x){p.adjust(x, method = "BH")}))
  
save(#ARE.ModGrp2COLOR,
     ldsc.enrich,
     ldsc.pval,
     ldsc.qval,
     file=file.ou.RData)

# pdf(file.ou.enrich.pdf, width=20, height=9)
# pheatmap(t(log10(ldsc.enrich+1)),
#            col=  colorRampPalette(brewer.pal(n = 7, name = "OrRd" ))(200),#colorRampPalette(c('#FFFFFF', "#FF0000"))(n = 200),#
#            cluster_cols = T, show_rownames = T, show_colnames = T, cluster_rows=T, 
#            main = paste0("log10 enrichment of GWAS signal at peaks")
#   )
# 
# 
# 
# dev.off()



# row.dists = as.dist(1-cor(t(ldsc.enrich[, gwas.kept]), method = "spearman"))
# row.hclust = hclust(row.dists)

#row.ordered = c(setdiff(rownames(ldsc.enrich), ARES.Modules.ordered), ARES.Modules.ordered)

pdf(file.ou.pval.pdf, width=20, height=15)  

# ldsc.enrich.filt= ldsc.enrich
# ldsc.enrich.filt[ldsc.pval>0.05]=1
# row.order= order(apply(ldsc.enrich.filt, 1, which.max), decreasing = F)
# corrplot(log10(ldsc.enrich+1)[row.order, ], 
#          p.mat=ldsc.pval[row.order, ],
#          #addCoef.col=1,
#          is.corr=F, 
#          insig = "blank", 
#          sig.level = .05,
#          method = "square"
#          #, order = "hclust"
#          )



#show enrichment
ldsc.enrich.log.filt= log10(ldsc.enrich+1)
ldsc.enrich.log.filt[ldsc.qval>0.25]=0 #ldsc.pval<=0.05
ldsc.enrich.log.filt= ldsc.enrich.log.filt[apply(ldsc.enrich.log.filt>0,1,any), ] #apply(t(ldsc.pval)<=0.05,2,any)
ldsc.enrich.log.filt = ldsc.enrich.log.filt#[, names(ARE.ModGrp2COLOR)]
row.order= order(apply(ldsc.enrich.log.filt, 1, which.max), decreasing = F)
ldsc.enrich.log.filt= ldsc.enrich.log.filt[row.order,]
corrplot(ldsc.enrich.log.filt, 
        # p.mat=t(ldsc.pval)[rownames(ldsc.enrich.filt), colnames(ldsc.enrich.filt)],
        is.corr=F, 
        #insig = "blank", sig.level = .05,
        method = "square"
         #, order = "hclust"
         )


# #show log p
# ldsc.p.log.filt = -log10(ldsc.pval)
# #ldsc.p.log.filt[t(ldsc.enrich)<2]=0
# ldsc.p.log.filt[ldsc.qval>0.2]= 0 #ldsc.p.log.filt< -log10(0.05)
# ldsc.p.log.filt= ldsc.p.log.filt[apply(ldsc.p.log.filt>0,1,any), ] #apply(t(ldsc.pval)<=0.05,2,any)
# ldsc.p.log.filt = ldsc.p.log.filt#[, names(ARE.ModGrp2COLOR)]
# row.order= order(apply(ldsc.p.log.filt, 1, which.max), decreasing = F)
# ldsc.p.log.filt= ldsc.p.log.filt[row.order,]
# corrplot(ldsc.p.log.filt, 
#         # p.mat=t(ldsc.pval)[rownames(ldsc.enrich.filt), colnames(ldsc.enrich.filt)],
#         is.corr=F, 
#         #insig = "blank", sig.level = .05,
#         method = "square"
#          #, order = "hclust"
#          )


#
dev.off()

#
f.o.tsv=paste0(dir.ou, "enrichAndpval.tsv")
write.table(ldsc.enrich.log.filt, file=f.o.tsv, sep="\t", quote=F)
