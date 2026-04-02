#visualize the gseapy results
library(ggrepel)

library(ggplot2)
library(data.table)


DRUGS.SLCT=c("aripiprazole", "trazdone", "topiramate", "quetiapine", "gabapentin", "valproate", 
        "haloperidol", "clozapine", "brexpiprazole")#BD

#pids.TFs=c(RORA="BRD-K20401833", FLI1="BRD-A62182663")


JOB.N=50
EXP.N=49550
EXP.N.PER.JOB= round(EXP.N/JOB.N)
# files.in.gsea=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood_", 
#                                 EXP.N.PER.JOB*(0:(JOB.N-1)), "-", EXP.N.PER.JOB*(1:JOB.N),
#                                 ".1000perm.tsv")
dir.in="b_comp2CompoundSig_cMap_V1.1/"
files.in.gsea=dir(dir.in, "1000perm.tsv", full.names=T)

file.in.meta= "~/lhou.compbio/data/expression/CMap/siginfo_beta.txt"
file.in.compound="~/lhou.compbio/data/expression/CMap/compoundinfo_beta.txt"
file.ou.pdf=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.pdf")
file.ou.RData=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.RData")
file.ou.txt=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugs.txt")
file.ou.2.txt=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugTargets.txt")
file.ou.3.txt=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.slct.sigID.txt")
#
compound.info=read.table(file.in.compound, quote="", comment="", sep="\t", header=T, row.names=NULL, stringsAsFactor=F)
pertID2cmap_name=lapply(split(compound.info$cmap_name, compound.info$pert_id), FUN=function(x) levels(factor(x)))
pertID2cmap_name.uniq=unlist(pertID2cmap_name[sapply(pertID2cmap_name, length)==1])


meta=fread(file.in.meta, sep="\t", header=T, stringsAsFactor=F, data.table=F)
rownames(meta)=meta$sig_id

deg.gsea.res.list=lapply(files.in.gsea,
FUN=function(f.i)
{
    d=read.table(f.i, sep="\t", header=T, row.names=1, stringsAsFactor=F)
})
deg.gsea.res.df=do.call(rbind, deg.gsea.res.list)


deg.gsea.res.df= data.frame(deg.gsea.res.df,
                        cell_iname=as.character(meta[deg.gsea.res.df$experiment, "cell_iname"]),
                        pert_id=as.character(meta[deg.gsea.res.df$experiment, "pert_id"]),
                        #sig_id=as.character(meta[deg.gsea.res.df$experiment, "sig_id"]),
                        stringsAsFactors=F)
deg.gsea.res.df$p.filt= #deg.gsea.res.df$upDEG.NES<0 & deg.gsea.res.df$dnDEG.NES>0 & 
                        (deg.gsea.res.df$upDEG.NOM_pval==0 & deg.gsea.res.df$dnDEG.NOM_pval==0)
deg.gsea.res.df$type="noSig"
deg.gsea.res.df$type[deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES<0]="opposite"
deg.gsea.res.df$type[deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES>0]="same"
deg.gsea.res.df$cmap_name=pertID2cmap_name.uniq[deg.gsea.res.df$pert_id]
deg.gsea.res.df$label= deg.gsea.res.df$cmap_name
deg.gsea.res.df$label[!(deg.gsea.res.df$label %in% DRUGS.SLCT)]=""
deg.gsea.res.df$label[!(deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES<0 & deg.gsea.res.df$dnDEG.NES>0)]=""
#deg.gsea.res.slct.df=deg.gsea.res.df[deg.gsea.res.df$slcted,]
#save(deg.gsea.res.df, 
    # deg.gsea.res.slct.df,
    # file=file.ou.RData)



#slect perturbation are quite robust
deg.gsea.res.filtbyp.df=deg.gsea.res.df[deg.gsea.res.df$p.filt, ]
exp.filt.NO=nrow(deg.gsea.res.filtbyp.df)
# #exp.passPvalCutoff=deg.gsea.res.df$upDEG.NOM_pval==0 & deg.gsea.res.df$dnDEG.NOM_pval==0
# exp.isConsist=deg.gsea.res.filtbyp.df$upDEG.NES<0 & deg.gsea.res.filtbyp.df$dnDEG.NES>0
# exp.isConsist.NO=sum(exp.isConsist)
# pertIDs=levels(factor(deg.gsea.res.df$pert_id))
# pertID2enrich=t(sapply(pertIDs,
# FUN=function(pid)
# {
#     exp.pid.NO=sum(deg.gsea.res.filtbyp.df$pert_id==pid)
#     exp.pid.consist.NO=sum(deg.gsea.res.filtbyp.df$pert_id==pid & exp.isConsist)
#     exp.rest.NO=nrow(deg.gsea.res.filtbyp.df)-sum(deg.gsea.res.filtbyp.df$pert_id==pid | exp.isConsist)
    
#     # exp.pid.slct.NO=sum(deg.gsea.res.df$pert_id==pid & deg.gsea.res.df$slcted)
#     # exp.rest.NO=nrow(deg.gsea.res.df)-sum(deg.gsea.res.df$pert_id==pid | deg.gsea.res.df$slcted)

#     fisher.res=fisher.test(matrix(c(exp.pid.consist.NO, 
#                                     exp.isConsist.NO-exp.pid.consist.NO, 
#                                     exp.pid.NO-exp.pid.consist.NO,
#                                     exp.rest.NO), 
#                                 ncol=2), alternative="greater")

#     return(c(p=fisher.res$p.value, or=fisher.res$estimate, prop=exp.pid.slct.NO/exp.pid.NO, size=exp.pid.NO))

# }))
# pertID2enrich=cbind(pertID2enrich, p.adj=p.adjust(pertID2enrich[,"p"], method="BH"))

save(#pertID2enrich,
    deg.gsea.res.df, 
    #deg.gsea.res.slct.df,
    file=file.ou.RData)


#
pert_ids.slct=deg.gsea.res.df$pert_id[deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES<0 & deg.gsea.res.df$dnDEG.NES>0]
compound.info.slct=compound.info[compound.info$pert_id %in% pert_ids.slct,]
write.table(compound.info.slct, file.ou.txt, sep="\t", quote=F, row.names=F)

targets=levels(factor(compound.info.slct$target))
targets=setdiff(targets, "\"\"")
write(targets, file=file.ou.2.txt)

sig.slct.id=deg.gsea.res.df$experiment[deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES *deg.gsea.res.df$dnDEG.NES<0]
write(sig.slct.id, file=file.ou.3.txt)


#
f.o.tsv=paste0(dir.in, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.tsv")
write.table(deg.gsea.res.df, f.o.tsv, sep="\t", quote=F)

#


#
pdf(file.ou.pdf, width=8, height=8)

layout(matrix(1:2, nrow=2))
hist(deg.gsea.res.df$upDEG.NOM_pval, main="upDEG.NOM_pval")
hist(deg.gsea.res.df$dnDEG.NOM_pval, main="dnDEG.NOM_pval")



#
pvalcutoff.NO=sum(deg.gsea.res.df$p.filt)
consist.NO=sum(deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES <0)
slct.NO=sum(deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES <0)
consist.p=fisher.test(matrix(c(slct.NO, 
                                consist.NO-slct.NO, 
                                pvalcutoff.NO-slct.NO,
                                nrow(deg.gsea.res.df)-sum(deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES <0 | deg.gsea.res.df$p.filt)), 
                                ncol=2), alternative="greater")$p.value


ggplot(deg.gsea.res.df, aes(x=upDEG.NES, y=dnDEG.NES)) +
  geom_point(aes(color=type, size=p.filt, alpha=p.filt))+
  scale_color_manual(values=c("noSig" = "gray" , "antagnostic"= "red", "synergistic"="blue"))+
  scale_size_manual(values=c("FALSE"=0.2, "TRUE"=1))+
  scale_alpha_manual(values=c("FALSE"=0.2, "TRUE"=1)) +
  ggtitle(paste0("consist effect enriched p.value: ", signif(consist.p, 3))) +
  geom_text_repel(aes(label = label), 
                    #max.overlaps = 10,  # Adjust for fewer overlaps
                    size = 3,          # Adjust text size
                    box.padding = 0.3, # Padding around labels
                    point.padding = 0.3) +
  theme_bw()



dev.off()

