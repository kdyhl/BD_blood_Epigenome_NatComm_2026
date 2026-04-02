#V1.1.2
#c.cmpAllSampClusterWithEHR_byMOFA_V1.1.2_addNumericFeature.R
#check all samples
#test other type of EHR data, including numeric, and multi-category ones

#V1.1.1
#fix bug at the  EHR fisher test 

#V1.1
#based on clusters from all samples
#add Yue's mixEHR, KS-test

#using EHR to annotate clusters of patients


R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/4.0/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))

library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)


# args = commandArgs(trailingOnly=TRUE)
# DIR.IN = args[1] 
# 
TYPE= "HiC.DPs"# "DPs" #"HiC.DPs"
CLU.RES = "clu.corHclust5" #"clu.louvain" #
TERM.SampN.CUTOFF=5 #codes with at least this  number would be considered
TERM.ENRICH.P.CUTOFF=0.02

#file.in.clu.RData=paste0("b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/", TYPE, "/ind.MOFA.DP-pval0.01.varProp0.005.patientclu.RData")
file.in.clu.RData=paste0("b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/", TYPE, "/ind.MOFA.DP-pval0.01.varProp0.005.sampClu_V1.1.RData")

# file.in.meta = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
# file.in.meta2 = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
#file.in.lab="~/lhou.compbio/data/Mayo_Bipolar/labs_w_dT_DawnOfTimesToSample_forSharing.csv"
#file.in.med1 = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/Med_source1.csv"
file.in.med2 = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/Med_source2.csv"
file.in.icd = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/icd_codes.csv"
#file.in.clinic = "/broad/compbio/lhou/data/Mayo_Bipolar/210415_new_broad_data_deidentified.csv"

# file.in.mixEHR = "~/lhou.compbio/data/Mayo_Bipolar/mayoclinic_BD_patient_topic_byYue_MixerEHR.txt"

file.in.LF.RData = "b_factorAnalysis_across5HMs_MOFA_V1.2_pvalCutoff/HiC.DPs/ind.MOFA.DP-pval0.01.varProp0.005.sampClu_V1.1.RData"


#file.in.lab = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/labs.csv"
#file.in.quest = "/broad/compbio/data/Mayo_BiPolar_ChIPseqAndEHR/outdata/ppi.csv"

#gsub("b_clusterPatients_basedOnPCs", "c_patientClu2EHR", DIR.IN)

file.ou.RData = gsub("b_factorAnalysis_across5HMs_MOFA", "c_sampClu2EHR_5HMs_byMOFA_V1.2_pvalCutoff", file.in.clu.RData)
file.ou.RData = gsub("sampClu_V1.1.RData", paste0("patientCluvsEHR.", CLU.RES, "_V1.1.2.RData"), file.ou.RData)
dir.create(dirname(file.ou.RData), showWarnings = F, recursive = T)
file.ou.EHR.enrich.pdf = gsub("RData", "enrich.pdf", file.ou.RData)
file.ou.hm.pdf = gsub("RData", "heatmap.pdf", file.ou.RData)


load(file.in.clu.RData)

#
EHR =list()
buf= read.table(file.in.med2, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
buf = buf[buf$MED_GENERIC!="",]
term2samp= split(paste0("id_", buf$newsubjectid), buf$MED_GENERIC) 
print(length(term2samp))#478 in total
EHR$med2 = lapply(term2samp, function(x){levels(factor(x))})
EHR$med2 = EHR$med2[sapply(EHR$med2, length) >=TERM.SampN.CUTOFF]
#
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
EHR$med2.filt = lapply(term2samp, function(x){levels(factor(x))})
EHR$med2.filt = EHR$med2.filt[sapply(EHR$med2.filt, length) >=TERM.SampN.CUTOFF]

buf= read.table(file.in.icd, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
term2samp= split(paste0("id_", buf$newsubjectid), buf$Dx_Desc)
print(length(term2samp)) #2680
EHR$icd = lapply(term2samp, function(x){levels(factor(x))})
EHR$icd = EHR$icd[sapply(EHR$icd, length) >=TERM.SampN.CUTOFF]
EHR$med2=NULL
#
# #
# buf= read.table(file.in.clinic, sep=",", header = T, row.names=NULL, stringsAsFactors = F, quote="\"", comment.char = "")
# rownames(buf)=paste0("id_", buf$newsubjectid)
# EHR$clinic=list()
# for(nm in colnames(buf))
# {
#   if(sum(is.na(buf[[nm]]))>=0.4*nrow(buf)) next
#   
#   if(sum(buf[[nm]]=="Yes", na.rm = T)>=TERM.SampN.CUTOFF &&  sum(buf[[nm]]=="No", na.rm = T)>=TERM.SampN.CUTOFF)
#   {
#     EHR$clinic[[nm]]=rownames(buf)[buf[[nm]]=="Yes"]
#     EHR$clinic[[nm]]=EHR$clinic[[nm]][!is.na(EHR$clinic[[nm]])]  
#   }
#   if(sum(buf[[nm]]==1, na.rm = T)>=TERM.SampN.CUTOFF &&  sum(buf[[nm]]==0, na.rm = T)>=TERM.SampN.CUTOFF)
#   {
#     EHR$clinic[[nm]]=rownames(buf)[buf[[nm]]==1]
#     EHR$clinic[[nm]]=EHR$clinic[[nm]][!is.na(EHR$clinic[[nm]])]
#   }
#   if(sum(buf[[nm]]=="Pos", na.rm = T)>=TERM.SampN.CUTOFF &&  sum(buf[[nm]]=="Neg", na.rm = T)>=TERM.SampN.CUTOFF)
#   {
#     EHR$clinic[[nm]]=rownames(buf)[buf[[nm]]=="Pos"]
#     EHR$clinic[[nm]]=EHR$clinic[[nm]][!is.na(EHR$clinic[[nm]])]
#   }
#   
# }
# #



#here use all samples
samps.mofa = rownames(mofa.fators.impute.scaled)
patient2clu= samp2grp[, CLU.RES]
names(patient2clu) = rownames(samp2grp)
# patient2clu.list$clu6 = samp2grp[samp2grp$group=="case", "clu.corHclust6"]
# names(patient2clu.list$clu6) = rownames(samp2grp)[samp2grp$group=="case"]
grp2patients= split(names(patient2clu), paste0("grp", patient2clu))


EHR.enrichment = list()
for(type in names(EHR))
{
  print(type)
  samps.ehr = levels(factor(unlist(EHR[[type]])))
  term2samp = EHR[[type]]
  
  EHR.enrichment[[type]]= data.frame()
  
  samps.ovlp = intersect(samps.ehr, samps.mofa)
  
  
  patients = names(patient2clu)
  patients.filt = intersect(samps.ovlp, patients) 
  
  
  for(grp in names(grp2patients))
  {
    patients.grp.filt = intersect(samps.ovlp, grp2patients[[grp]]) 
    
    for(term in names(term2samp))
    {
      #patient group vs rest patient
      samps.term.filt =  intersect(patients.filt, term2samp[[term]])
      samps.TermAndGrp = length(intersect(patients.grp.filt, samps.term.filt))
      samps.termSpec = length(setdiff(samps.term.filt, patients.grp.filt))
      samps.grpSpec = length(setdiff(patients.grp.filt, samps.term.filt))
      samps.bg = length(setdiff(patients.filt, union(samps.term.filt, patients.grp.filt)))
      
      fisher.res1 = fisher.test(matrix(c(samps.TermAndGrp, samps.termSpec, samps.grpSpec, samps.bg), ncol=2), alternative="greater")
      
      
      #
      #whole patients vs control
      samps.term.filt =  intersect(samps.ovlp, term2samp[[term]])
      samps.TermAndDis = length(intersect(samps.term.filt, patients.filt))
      samps.termSpec = length(setdiff(samps.term.filt, patients.filt))
      samps.disSpec = length(setdiff(patients.filt, samps.term.filt))
      samps.bg = length(setdiff(samps.ovlp, union(samps.term.filt, patients.filt)))
      
      fisher.res0 = fisher.test(matrix(c(samps.TermAndDis, samps.termSpec, samps.disSpec, samps.bg), ncol=2), alternative="greater")
      
      
      EHR.enrichment[[type]]=rbind(EHR.enrichment[[type]],
                                   data.frame(patientGroup = grp,
                                              EHR.term= term,
                                              grp.OR=fisher.res1$estimate,
                                              grp.OR.pval=fisher.res1$p.value,
                                              total.OR = fisher.res0$estimate,
                                              total.OR.pval = fisher.res0$p.val,
                                              stringsAsFactors = F)
      )
    }
  }
}

#clinic data 
# #patient only




save(EHR,
     EHR.enrichment,
     file=file.ou.RData)


#
load(file.in.LF.RData)

pdf(file.ou.EHR.enrich.pdf, height=5, width=20) 

EHR.terms.slct = list()
for(type in setdiff(names(EHR.enrichment), "mixEHR"))
{
  
  df = EHR.enrichment[[type]]
  df$total.OR[df$total.OR==Inf] = max(df$total.OR[df$total.OR!=Inf] )
  
  #pdf(paste0(file.ou.pdf.pref, type, ".pdf"), height=4, width=18)  
  terms.slct = df[df$grp.OR.pval<= TERM.ENRICH.P.CUTOFF, "EHR.term"]
  terms.slct.grp = df[df$grp.OR.pval<= TERM.ENRICH.P.CUTOFF, "patientGroup"]
  EHR.terms.slct[[type]]=data.frame(term=terms.slct,
                                    grp=terms.slct.grp,
                                    stringsAsFactors = F)
  
  p=ggplot(df, 
           aes(x=total.OR, y=grp.OR, label=EHR.term, color= (grp.OR.pval<=TERM.ENRICH.P.CUTOFF)))+
    #xlim(0, max(df$total.OR)*1.1) +
    geom_point(size=1) + #-log10(grp.OR.pval)
    facet_wrap(~patientGroup, nrow=1) +
    geom_hline(yintercept = 1, 
               color="black", 
               linetype="dashed", size=1)+
    theme_bw()+
    xlab("enrichment in patients vs healthy") +
    ylab("enrichment in patient group vs all patients") +
    geom_text_repel(
      #aes(color= (grp.OR.pval<=0.001)),
      data   = subset(df, EHR.term %in% terms.slct),
      #nudge_y       = 36 - subset(dat, mpg > 30)$mpg,
      size = 3,
      box.padding = 0.3
      #segment.size  = 0.2,
      #segment.color = "grey50",
      #direction     = "x"
    ) +
    ggtitle(paste(type, "(", paste(summary(factor(patient2clu)), collapse="," ),")"))
  print(p)
  #dev.off()
  
}


# #mixEHR
# df=EHR.enrichment[["mixEHR"]]
# terms.slct = levels(factor(df[df$grp.pval<= TERM.ENRICH.P.CUTOFF, "EHR.term"]))
# p=ggplot(df, 
#          aes(x=total.diff, y=grp.diff, label=EHR.term, color= (grp.pval<=TERM.ENRICH.P.CUTOFF)))+
#   #xlim(0, max(df$total.OR)*1.1) +
#   geom_point(size=1) + #-log10(grp.OR.pval)
#   facet_wrap(~patientGroup, nrow=1) +
#   # geom_hline(yintercept = 1, 
#   #            color="black", 
#   #            linetype="dashed", size=1)+
#   theme_bw()+
#   xlab("difference in patients vs healthy") +
#   ylab("difference in patient group vs others patients") +
#   geom_text_repel(
#     #aes(color= (grp.OR.pval<=0.001)),
#     data   = subset(df, EHR.term %in% terms.slct),
#     #nudge_y       = 36 - subset(dat, mpg > 30)$mpg,
#     size = 3,
#     box.padding = 0.3
#     #segment.size  = 0.2,
#     #segment.color = "grey50",
#     #direction     = "x"
#   ) +
#   ggtitle(paste("mixEHR Module", "(", paste(summary(factor(patient2clu)), collapse="," ),")"))
# print(p)


dev.off()

#visualize terms along with latent factor

EHR.terms.slct.matrix = c()
patients = names(patient2clu)
EHR.term.enrichedGrp=c()
for(type in names(EHR.terms.slct))
{
  for(term in EHR.terms.slct[[type]]$term)
  {
    pat2EHR = rep(0, length(patients))
    names(pat2EHR) = patients
    
    
    if(type !="clinic")
    {
      pat2EHR[intersect(patients, EHR[[type]][[term]])] =  1
      
    }else
    {
      if(class(clinic[patients, term])=="numeric")
      {
        x=clinic[patients, term]
        pat2EHR[patients] = (x-min(x, na.rm = T))/(max(x, na.rm=T)-min(x, na.rm=T))
      }else
      {
        pat2EHR[patients][clinic[patients, term] == "Yes" |
                            clinic[patients, term] == "Pos"|
                            clinic[patients, term] == 1]=1
      }
        
      
    }
    EHR.terms.slct.matrix = cbind(EHR.terms.slct.matrix, term=pat2EHR)
    colnames(EHR.terms.slct.matrix)[ncol(EHR.terms.slct.matrix)]=term
    
  }
 
  EHR.term.enrichedGrp=c(EHR.term.enrichedGrp, EHR.terms.slct[[type]]$grp)
  
}

EHR.terms.slct.matrix.filt=EHR.terms.slct.matrix[, ! (colnames(EHR.terms.slct.matrix) %in% c("crp_z", "crp"))] #redundant
EHR.term.enrichedGrp.filt= EHR.term.enrichedGrp[! (colnames(EHR.terms.slct.matrix) %in% c("crp_z", "crp"))]

library(circlize)
col_fun = colorRamp2(c(0, 1), c("grey90", "red"))
pdf(file.ou.hm.pdf, height=12, width=15) 

ht1=Heatmap(t(mofa.fators.impute.scaled), 
            column_split = paste0("grp", samp2grp[, "clu.corHclust5"]),
            cluster_column_slices = FALSE,
            name="clu.corHclust5")
ht2=Heatmap(t(EHR.terms.slct.matrix.filt),
            col=col_fun,
            row_split= factor(EHR.term.enrichedGrp.filt, levels=paste0("grp", 1:5)),
            cluster_row_slices =F,
            row_title_rot = 0
            )
ht_list = ht1 %v% ht2
draw(ht_list, 
     ht_gap = unit(1, "cm"),
     padding = unit(c(2, 2, 2, 60), "mm"))

dev.off()



save(EHR,
     EHR.enrichment,
     patient2clu,
     EHR.terms.slct.matrix.filt,
     EHR.term.enrichedGrp.filt,
     file=file.ou.RData)



