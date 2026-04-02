#using EHR to anntoate clusters of patients
#and cell fractions
#and BD, sex, age


R.LIBS.PATH = .libPaths()
.libPaths(c("~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/", R.LIBS.PATH))
TYPE="HiC.DPs" #"DPs" #

# library(ggplot2)
# library(ggrepel)
library(MOFA)
library(pheatmap)
# args = commandArgs(trailingOnly=TRUE)
# DIR.IN = args[1] 
HMs= c("H3K27ac", "H3K36me3", "H3K4me3",  "H3K4me1",  "H3K27me3")
#file.in.factor.RData =  paste0("b_factorAnalysis_across5HMs_MOFA/", TYPE, "/ind.MOFA.topPeak10000.varProp0.01.RData")# paste0("b_factorAnalysis_acrossHMs_MOFA/HiC.DPs/ind.MOFA.topPeak5000.varProp0.01.RData")
#paste0("b_clusterInds_acrossHMs_MOFA/DPs/ind.MOFA.topPeak5000.RData")
file.in.factor.RData =  paste0("b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/", TYPE, "/ind.MOFA.DP-pval0.01.varProp0.005.RData")#

# file.in.meta = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
file.in.meta2 = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
#file.in.lab="~/lhou.compbio/data/Mayo_Bipolar/labs_w_dT_DawnOfTimesToSample_forSharing.csv"
#file.in.med1 = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/Med_source1.csv"
file.in.meta.csv = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"

file.in.med2 = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/Med_source2.csv"
file.in.icd = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/icd_codes.csv"
file.in.clinic = "/broad/compbio/lhou/data/Mayo_Bipolar/210415_new_broad_data_deidentified.csv"

# files.in.dp.RData= paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
# names(files.in.dp.RData) = HMs

files.in.deconv.RData =  system(paste0("find ../deconv/b_real_estFract_V1.4.4.1_nnPoisR_noIntcpt_peakNorm_6cellType_autoSig/*/estFraction.results.nnPoisR.FiltQ*.RData"), intern=T)
names(files.in.deconv.RData)= sub(".*/(H3K.*)", "\\1", dirname(files.in.deconv.RData))
#file.in.lab = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/labs.csv"
#file.in.quest = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/ppi.csv"

#gsub("b_clusterPatients_basedOnPCs", "c_patientClu2EHR", DIR.IN)

file.ou.RData = gsub(".RData", ".FactorvsCovar.RData", file.in.factor.RData)
file.ou.pdf = gsub(".RData", ".FactorvsCovar.pdf", file.in.factor.RData)


TERM.SampN.CUTOFF = 5

#
meta = read.table(file.in.meta.csv, sep=",", header=T, row.names=1)
rownames(meta)=paste("id", rownames(meta), sep="_")
samps.case = rownames(meta)[meta$group=="case"]
meta2=read.table(file.in.meta2, sep=",", header=T, row.names=1)
rownames(meta2)=paste0("id_", rownames(meta2))
samp.meta=cbind(group=meta$group, meta2[rownames(meta), c("Gender", "age_atsamp")])


#
# hm.coVar = lapply(files.in.dp.RData,
# FUN=function(f.i.dp.RData)
# {
#   load(f.i.dp.RData)    
#   design.matrix.mod2[, -1]
# })
# names(hm.coVar)=paste0(names(files.in.dp.RData), "_covar")



hm.cellFractEst=lapply(names(files.in.deconv.RData),
FUN=function(hm)
{
  f.i.RData = files.in.deconv.RData[hm]
  load(f.i.RData)
  Grp.slct.props.unNorm.mean = sapply(colnames(Grp.slct.props.sampling.unNorm[[1]]),
  FUN=function(ref)
  {
    ref.props= sapply(Grp.slct.props.sampling.unNorm, FUN=function(x) x[,ref])
    return(apply(ref.props, 1, mean))
  })

  colnames(Grp.slct.props.unNorm.mean)=paste0(hm, "_", colnames(Grp.slct.props.unNorm.mean))

  return(Grp.slct.props.unNorm.mean)
})
names(hm.cellFractEst)=paste0(names(files.in.deconv.RData), "_CellFract")
#
load(file.in.factor.RData)
samps.all = rownames(MOFAobject@Expectations$Z)

EHR =list()
buf= read.table(file.in.med2, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
buf = buf[buf$MED_GENERIC!="",]
term2samp= split(paste0("id_", buf$newsubjectid), buf$MED_GENERIC)
EHR$med2 = sapply(term2samp, 
function(x)
{
  buf = rep(0, length(samps.all))
  names(buf) = samps.all
  
  buf[intersect(samps.all, levels(factor(x)))]=1
  
  return(buf)
})
EHR$med2 = EHR$med2[,apply(EHR$med2, 2, sum) >=TERM.SampN.CUTOFF]

#for those medicine at least taken twice
samp2term= split(buf$MED_GENERIC,paste0("id_", buf$newsubjectid))
sampAndterm.filt = c()
for(samp in names(samp2term))
{
  term2freq= summary(factor(samp2term[[samp]]),maxsum=length(samp2term[[samp]]))
  
  if(sum(term2freq>=2)>=1)
  {
    sampAndterm.filt = rbind(sampAndterm.filt,
                           cbind(samp=samp, term=names(term2freq)[term2freq>=2]))
  }
}
term2samp=split(sampAndterm.filt[,"samp"], sampAndterm.filt[,"term"])
EHR$med2.filt = sapply(term2samp, 
function(x)
{
  buf = rep(0, length(samps.all))
  names(buf) = samps.all
  
  buf[intersect(samps.all, levels(factor(x)))]=1
  
  return(buf)
})
EHR$med2.filt = EHR$med2.filt[,apply(EHR$med2.filt, 2, sum) >=TERM.SampN.CUTOFF]


#
buf= read.table(file.in.icd, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
term2samp= split(paste0("id_", buf$newsubjectid), buf$Dx_Desc)
#EHR$icd = lapply(term2samp, function(x){levels(factor(x))})
EHR$icd = sapply(term2samp, 
function(x)
{
  buf = rep(0, length(samps.all))
  names(buf) = samps.all
  
  buf[intersect(samps.all, levels(factor(x)))]=1
  
  return(buf)
})
EHR$icd = EHR$icd[,apply(EHR$icd, 2, sum) >=TERM.SampN.CUTOFF]

#
buf= read.table(file.in.clinic, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
rownames(buf)=paste0("id_", buf$newsubjectid)
buf$newsubjectid=NULL
#buf=buf[intersect(rownames(buf),samps.all),]

EHR$clinic=buf[intersect(rownames(buf),samps.all),]


# #
# buf= read.table(file.in.lab, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
# term2samp= split(paste0("id_", buf$newsubjectid), buf$Dx_Desc)

covars.all=c(list(meta=samp.meta), hm.cellFractEst, EHR)

#remove clinic which is only for patients
R2.factorVScovar = lapply(covars.all[c("meta", "H3K27ac_CellFract", "H3K36me3_CellFract", "H3K4me1_CellFract", "med2.filt", "icd")],
FUN=function(ehr)
{
  R2=sapply(colnames(ehr), 
  FUN=function(ehr.nm)
  {
    print(ehr.nm)
    x=ehr[, ehr.nm]
    names(x) = rownames(ehr)
    if(length(levels(factor(x[x!="" & (!is.na(x))])))<=1)
    {
      return(rep(NA, ncol(MOFAobject@Expectations$Z)))
    }
    samps.ovlp = intersect(names(x), samps.all)
    sapply(colnames(MOFAobject@Expectations$Z[,]),
    FUN=function(factor.nm)
    {
      #print(factor.nm)
      y=MOFAobject@Expectations$Z[samps.ovlp, factor.nm]
      mod = lm(y~x[samps.ovlp])
      summary(mod)$r.squared
    })
    
  })
  
  
  R2 = R2[, !is.na(R2[1,])]
})

for(type in names(R2.factorVScovar))
{
  colnames(R2.factorVScovar[[type]])=paste(colnames(R2.factorVScovar[[type]]), type, sep="_")
}
R2.factorVScovar.matrix= do.call(cbind, R2.factorVScovar)


#
R2.factorVScovar.patient = list()
pval.factorVScovar.patient=list()

for(categ in c("meta", "H3K27ac_CellFract", "H3K36me3_CellFract", "H3K4me1_CellFract", "med2.filt", "icd", "clinic"))
{
  ehr= covars.all[[categ]]
  R2andP=sapply(colnames(ehr), 
  FUN=function(ehr.nm)
  {
    print(ehr.nm)
    x=ehr[, ehr.nm]
    names(x) = rownames(ehr)
    
    samps.ovlp = intersect(samps.case, intersect(names(x), samps.all))

    x.ovlp=x[samps.ovlp]
    x.ovlp[x.ovlp==""]=NA
    if(length(levels(factor(x.ovlp[!is.na(x.ovlp)])))<=1)
    {
      return(c(rep(NA, ncol(MOFAobject@Expectations$Z)),  #R2
               rep(NA, ncol(MOFAobject@Expectations$Z)))) #pval
    }
    # if(length(levels(factor(x.ovlp[x.ovlp!="" & (!is.na(x.ovlp))])))<=1)
    # {
    #   return(rep(NA, ncol(MOFAobject@Expectations$Z)))
    # }

    r2Andp= sapply(colnames(MOFAobject@Expectations$Z),
    FUN=function(factor.nm)
    {
      #print(factor.nm)
      y.ovlp=MOFAobject@Expectations$Z[samps.ovlp, factor.nm]
     # mod = lm(y~x.ovlp)
      #summary(mod)$r.squared
      
      mod.sum = tryCatch(
        summary(lm(y.ovlp~x.ovlp)),
        error = function(e) 
        {
          print(e)
          return(c(NA,NA))
        }
      )
      
      if((!is.na(mod.sum)) && nrow(mod.sum$coeff)==2)
      {
        return(c(mod.sum$r.squared,
                 mod.sum$coef[2,4]))
      }else
      {
        return(c(NA,NA))
      }
    })
    r2Andp=c(r2Andp[1,], r2Andp[2,])
  })
  
  R2=R2andP[1:(nrow(R2andP)/2), ]
  pval=R2andP[(nrow(R2andP)/2+1):nrow(R2andP),]
  
  R2.factorVScovar.patient[[categ]] = R2[, !is.na(R2[1,])]
  pval.factorVScovar.patient[[categ]]=pval[, !is.na(R2[1,])]
}
for(categ in names(R2.factorVScovar.patient))
{
  colnames(R2.factorVScovar.patient[[categ]])=paste(colnames(R2.factorVScovar.patient[[categ]]), categ, sep="_")
  colnames(pval.factorVScovar.patient[[categ]])=paste(colnames(pval.factorVScovar.patient[[categ]]), categ, sep="_")
}
R2.factorVScovar.patient.matrix= do.call(cbind, R2.factorVScovar.patient)
pval.factorVScovar.patient.matrix= do.call(cbind, pval.factorVScovar.patient)
qval.factorVScovar.patient.matrix= apply(pval.factorVScovar.patient.matrix, 2, p.adjust, method="BH")


# #
# pval.factorVScovar.patient = lapply(covars.all[c("meta", "H3K27ac_CellFract", "H3K36me3_CellFract", "H3K4me1_CellFract", "med2.filt", "icd", "clinic")],
# FUN=function(ehr)
# {
#   pval=sapply(colnames(ehr), 
#   FUN=function(ehr.nm)
#   {
#     print(ehr.nm)
#     x=ehr[, ehr.nm]
#     names(x) = rownames(ehr)
#     
#     samps.ovlp = intersect(samps.case, intersect(names(x), samps.all))
#     
#     x.ovlp[x.ovlp==""]=NA
#     if(length(levels(factor(x.ovlp[!is.na(x.ovlp)])))<=1)
#     {
#       return(rep(NA, ncol(MOFAobject@Expectations$Z)))
#     }
#     
#     sapply(colnames(MOFAobject@Expectations$Z[,]),
#     FUN=function(factor.nm)
#     {
#       #print(factor.nm)
#       #y=MOFAobject@Expectations$Z[samps.ovlp, factor.nm]
#       y.ovlp=MOFAobject@Expectations$Z[samps.ovlp, factor.nm]
#       # mod = lm(y~x.ovlp)
#       #summary(mod)$r.squared
#       
#       mod.sum = tryCatch(
#         summary(lm(y.ovlp~x.ovlp)),
#         error = function(e) 
#         {
#           print(e)
#           return(NA)
#         }
#       )
#      #mod.sum = summary(lm(y~x.ovlp))
#       
#       if((!is.na(mod.sum)) && nrow(mod.sum$coeff)==2)
#       {
#         return(mod.sum$coeff[2,4])
#       }else
#       {
#         return(NA)
#       }
#     })
#     
#   })
#   pval = pval[, !is.na(pval[1,])]
# })
# 
# for(type in names(pval.factorVScovar.patient))
# {
#   colnames(pval.factorVScovar.patient[[type]])=paste(colnames(pval.factorVScovar.patient[[type]]), type, sep="_")
# }
# pval.factorVScovar.patient.matrix= do.call(cbind, pval.factorVScovar.patient)

save(covars.all,
     R2.factorVScovar,
     R2.factorVScovar.patient,
     R2.factorVScovar.matrix,
     R2.factorVScovar.patient.matrix,
     pval.factorVScovar.patient.matrix,
     qval.factorVScovar.patient.matrix,
     file=file.ou.RData)

pdf(file.ou.pdf, height=15, width=10)  


R2.factorVScovar.matrix.filt= R2.factorVScovar.matrix[, apply(abs(R2.factorVScovar.matrix),2, max)>=0.2]
pheatmap(t(R2.factorVScovar.matrix.filt), main="with all samples")


# R2.factorVScovar.patient.matrix.filt= R2.factorVScovar.patient.matrix[, apply(abs(R2.factorVScovar.patient.matrix),2, max, na.rm=T)>=0.1]
# pheatmap(t(R2.factorVScovar.patient.matrix.filt), 
#          #display_numbers=t(display_numbers),
#          main="patients only")

# q<=0.05, R2>=0.1
R2.factorVScovar.patient.matrix.filt= t(R2.factorVScovar.patient.matrix[, apply(qval.factorVScovar.patient.matrix, 2, min, na.rm=T)<=0.05 &
                                                                        apply(abs(R2.factorVScovar.patient.matrix),2, max)>=0.1])
R2.factorVScovar.patien.isSig.matrix.filt=t(qval.factorVScovar.patient.matrix[, rownames(R2.factorVScovar.patient.matrix.filt)])
display_numbers=matrix("",
                       nrow=nrow(R2.factorVScovar.patien.isSig.matrix.filt),
                       ncol=ncol(R2.factorVScovar.patien.isSig.matrix.filt))
display_numbers[R2.factorVScovar.patien.isSig.matrix.filt<=0.05]="*"
rownames(display_numbers)=rownames(R2.factorVScovar.patien.isSig.matrix.filt)
colnames(display_numbers)=colnames(R2.factorVScovar.patien.isSig.matrix.filt)
#R2.factorVScovar.patient.matrix.filt.colSorted=R2.factorVScovar.patient.matrix.filt[, paste0("LF", 1:2)]
R2.factorVScovar.patient.matrix.filt.rowOrder= order(apply(R2.factorVScovar.patient.matrix.filt[, paste0("LF", 1:24)], 1, which.max), decreasing = F)
R2.factorVScovar.patient.matrix.filt.colRowSorted=R2.factorVScovar.patient.matrix.filt[R2.factorVScovar.patient.matrix.filt.rowOrder, paste0("LF", 1:24)]
display_numbers.colRowSorted=display_numbers[R2.factorVScovar.patient.matrix.filt.rowOrder, paste0("LF", 1:24)]
pheatmap(R2.factorVScovar.patient.matrix.filt.colRowSorted,
         display_numbers=display_numbers.colRowSorted,
         cluster_rows=F,
         cluster_cols=F,
         main="patients only, q 0.05, R^2 0.1")
dev.off()


#
f.ou.tsv = gsub(".RData", ".FactorvsCovar.heatmap.tsv", file.in.factor.RData)
write.table(R2.factorVScovar.patient.matrix.filt.colRowSorted, f.ou.tsv, sep="\t", quote=F)



#
# # q<=0.05, R2>=0.2
# R2.factorVScovar.patient.matrix.filt= R2.factorVScovar.patient.matrix[, apply(qval.factorVScovar.patient.matrix, 2, min, na.rm=T)<=0.05 &
#                                                                         apply(abs(R2.factorVScovar.patient.matrix),2, max)>=0.2]
# R2.factorVScovar.patien.isSig.matrix.filt=qval.factorVScovar.patient.matrix[, colnames(R2.factorVScovar.patient.matrix.filt)]
# display_numbers=matrix("",
#                        nrow=nrow(R2.factorVScovar.patien.isSig.matrix.filt),
#                        ncol=ncol(R2.factorVScovar.patien.isSig.matrix.filt))
# display_numbers[R2.factorVScovar.patien.isSig.matrix.filt<=0.05]="*"
# pheatmap(t(R2.factorVScovar.patient.matrix.filt),
#          display_numbers=t(display_numbers),
#          main="patients only, q 0.05, R^2 0.2")


# for(type in names(R2.factorVScovar))
# {
#   print(type)
#   R2 = R2.factorVScovar[[type]]
#   R2.filt= R2[, apply(abs(R2),2,max)>=0.05]#apply(abs(R2),2,max)>=0.1]
  
#   if(grepl("covar", type))
#   {
#     pheatmap(t(R2), main=type)
#   }else
#   {
#     pheatmap(t(R2.filt), main=type)
#   }
  
# }

# for(type in names(R2.factorVScovar.patient))
# {
#   print(type)
#   R2 = R2.factorVScovar.patient[[type]]
#   R2.filt= R2[, apply(abs(R2),2,max)>=0.05]#apply(abs(R2),2,max)>=0.1]
  
#   if(grepl("covar", type))
#   {
#     pheatmap(t(R2), main=paste0(type, " patient only"))
#   }else
#   {
#     pheatmap(t(R2.filt), main=paste0(type, " patient only"))
#   }
  
  
# }


#dev.off()


