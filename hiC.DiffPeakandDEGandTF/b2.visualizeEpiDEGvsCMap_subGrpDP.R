#V1.2
#b2.visualizeEpiDEGvsCMap_subGrpDP_V1.2_cytoscape.R
#output files for cytoscape 

#V1.1
#ATC code from drugbank


#DP from subgrp
#add ATC classification for the drugs and targets

#visualize the gseapy results

library(ggrepel)
library(ggplot2)
library(data.table)
library(igraph)
# library(httr)
# library(jsonlite)
library(dplyr)


DRUGS.SLCT=c("aripiprazole", "trazdone", "topiramate", "quetiapine", "gabapentin", "valproate", 
        "haloperidol", "clozapine", "brexpiprazole")#BD


#pids.TFs=c(RORA="BRD-K20401833", FLI1="BRD-A62182663")

dir.in="b_comp2CompoundSig_cMap_inflamSubGrps_V1.1/a3_V1.2_subGrps/" #a3_V1.2.1_orderedDEGs/"
files.in.gsea=dir(dir.in, "1000perm.tsv", full.names=T)

file.in.meta= "~/lhou.compbio/data/expression/CMap/siginfo_beta.txt"
file.in.compound="~/lhou.compbio/data/expression/CMap/compoundinfo_beta.txt"

file.in.drugBank.RDS="~/lhou.compbio/data/drugs/drugbank/drugbank_ATC_syn_v5.1.12.RDS"

dir.ou=paste0(dir.in, "b2_vis_V1.2_cytoscape/")
dir.create(dir.ou, recursive=T, showWarnings=F)
file.ou.RData=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.RData")
files.ou.drugs.txt=c(inflamSubGrps=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugs.inflamSubGrps.txt"),
                    nonInflamSubGrps=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugs.nonInflamSubGrps.txt"))
files.ou.targets.txt=c(inflamSubGrps=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugTargets.inflamSubGrps.txt"),
                    nonInflamSubGrps=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugTargets.nonInflamSubGrps.txt"))
files.ou.targets.uniq.txt=c(inflamSubGrps=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugTargets.inflamSubGrps.uniq.txt"),
                    nonInflamSubGrps=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.potentialDrugTargets.nonInflamSubGrps.uniq.txt"))

file.ou.GSEA.pdf=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.pdf")
file.ou.drugAndTargetsNet.pref=paste0(dir.ou, "Drug2Target.network.forSubgrps.drugbank.")
file.ou.nodes=paste0(file.ou.drugAndTargetsNet.pref, "nodes.tsv")
# file.ou.drugAndTargetsNet.both.slct=paste0(dir.ou, "Drug2Target.network.forSubgrps.drugbank.both.slct")
# file.ou.drugAndTargetsNet.subgrps=paste0(dir.ou, "Drug2Target.network.forSubgrps.drugbank.")
file.ou.drugAndATC.tsv=paste0(dir.ou, "Drug2ATC.drugbank.tsv")



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

#
deg.gsea.res.df= data.frame(deg.gsea.res.df,
                        cell_iname=as.character(meta[deg.gsea.res.df$experiment, "cell_iname"]),
                        pert_id=as.character(meta[deg.gsea.res.df$experiment, "pert_id"]),
                        stringsAsFactors=F)
deg.gsea.res.df$cmap_name=pertID2cmap_name.uniq[deg.gsea.res.df$pert_id]
deg.gsea.res.df$inflamSubGrps.dn_up.DEG.NES=deg.gsea.res.df$inflamSubGrps.dnDEG.NES-deg.gsea.res.df$inflamSubGrps.upDEG.NES
deg.gsea.res.df$nonInflamSubGrps.dn_up.DEG.NES=deg.gsea.res.df$nonInflamSubGrps.dnDEG.NES-deg.gsea.res.df$nonInflamSubGrps.upDEG.NES


deg.gsea.res.df$inflamSubGrps.p.filt= (deg.gsea.res.df$inflamSubGrps.upDEG.NOM_pval==0 & deg.gsea.res.df$inflamSubGrps.dnDEG.NOM_pval==0)
deg.gsea.res.df$inflamSubGrps.type="noSig"
deg.gsea.res.df$inflamSubGrps.type[deg.gsea.res.df$inflamSubGrps.p.filt & sign(deg.gsea.res.df$inflamSubGrps.upDEG.NES) * sign(deg.gsea.res.df$inflamSubGrps.dnDEG.NES)==-1]="opposite"
deg.gsea.res.df$inflamSubGrps.type[deg.gsea.res.df$inflamSubGrps.p.filt & sign(deg.gsea.res.df$inflamSubGrps.upDEG.NES) * sign(deg.gsea.res.df$inflamSubGrps.dnDEG.NES)==1]="same"
#
deg.gsea.res.df$inflamSubGrps.label= deg.gsea.res.df$cmap_name
deg.gsea.res.df$inflamSubGrps.label[!(deg.gsea.res.df$inflamSubGrps.label %in% DRUGS.SLCT)]=""
deg.gsea.res.df$inflamSubGrps.label[!(deg.gsea.res.df$inflamSubGrps.p.filt & 
                                        deg.gsea.res.df$inflamSubGrps.upDEG.NES<0 & 
                                        deg.gsea.res.df$inflamSubGrps.dnDEG.NES>0)]=""
#
deg.gsea.res.df$nonInflamSubGrps.p.filt= (deg.gsea.res.df$nonInflamSubGrps.upDEG.NOM_pval==0 & deg.gsea.res.df$nonInflamSubGrps.dnDEG.NOM_pval==0)
deg.gsea.res.df$nonInflamSubGrps.type="noSig"
deg.gsea.res.df$nonInflamSubGrps.type[deg.gsea.res.df$nonInflamSubGrps.p.filt & sign(deg.gsea.res.df$nonInflamSubGrps.upDEG.NES) * sign(deg.gsea.res.df$nonInflamSubGrps.dnDEG.NES)==-1]="opposite"
deg.gsea.res.df$nonInflamSubGrps.type[deg.gsea.res.df$nonInflamSubGrps.p.filt & sign(deg.gsea.res.df$nonInflamSubGrps.upDEG.NES) * sign(deg.gsea.res.df$nonInflamSubGrps.dnDEG.NES)==1]="same"
#
deg.gsea.res.df$nonInflamSubGrps.label= deg.gsea.res.df$cmap_name
deg.gsea.res.df$nonInflamSubGrps.label[!(deg.gsea.res.df$nonInflamSubGrps.label %in% DRUGS.SLCT)]=""
deg.gsea.res.df$nonInflamSubGrps.label[!(deg.gsea.res.df$nonInflamSubGrps.p.filt & 
                                        deg.gsea.res.df$nonInflamSubGrps.upDEG.NES<0 & 
                                        deg.gsea.res.df$nonInflamSubGrps.dnDEG.NES>0)]=""
#
deg.gsea.res.df$twoSubgrps="nonCandidate"
is.inflam.drugs=deg.gsea.res.df$inflamSubGrps.p.filt & 
                deg.gsea.res.df$inflamSubGrps.upDEG.NES<0 & 
                deg.gsea.res.df$inflamSubGrps.dnDEG.NES>0
is.nonInflam.drugs=deg.gsea.res.df$nonInflamSubGrps.p.filt &
                deg.gsea.res.df$nonInflamSubGrps.upDEG.NES<0 &
                deg.gsea.res.df$nonInflamSubGrps.dnDEG.NES>0
deg.gsea.res.df$twoSubgrps[is.inflam.drugs]="is.inflam.drugs"
deg.gsea.res.df$twoSubgrps[is.nonInflam.drugs]="is.nonInflam.drugs"
deg.gsea.res.df$twoSubgrps[is.inflam.drugs & is.nonInflam.drugs]="is.both.drugs"


# deg.gsea.res.df$label[!(deg.gsea.res.df$label %in% DRUGS.SLCT)]=""
# deg.gsea.res.df$label[!(deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES<0 & deg.gsea.res.df$dnDEG.NES>0)]=""

#deg.gsea.res.slct.df=deg.gsea.res.df[deg.gsea.res.df$inflamSubGrps.p.filt,]
#save(deg.gsea.res.df, 
    # deg.gsea.res.slct.df,
    # file=file.ou.RData)



# #slect perturbation are quite robust
# deg.gsea.res.filtbyp.df=deg.gsea.res.df[deg.gsea.res.df$p.filt, ]
# exp.filt.NO=nrow(deg.gsea.res.filtbyp.df)
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


# drugs and targets
#inflam
pert_ids.slct=list()
pert_ids.slct$inflamSubGrps=deg.gsea.res.df$pert_id[deg.gsea.res.df$inflamSubGrps.p.filt & 
                                    deg.gsea.res.df$inflamSubGrps.upDEG.NES<0 & 
                                    deg.gsea.res.df$inflamSubGrps.dnDEG.NES>0]                                    
pert_ids.slct$inflamSubGrps=levels(factor(pert_ids.slct$inflamSubGrps))
pert_ids.slct$nonInflamSubGrps=deg.gsea.res.df$pert_id[deg.gsea.res.df$nonInflamSubGrps.p.filt & 
                                    deg.gsea.res.df$nonInflamSubGrps.upDEG.NES<0 & 
                                    deg.gsea.res.df$nonInflamSubGrps.dnDEG.NES>0]                                    
pert_ids.slct$nonInflamSubGrps=levels(factor(pert_ids.slct$nonInflamSubGrps))

targets.slct=lapply(names(pert_ids.slct),
FUN=function(subgrp.nm)
{
    compound.info.slct=compound.info[compound.info$pert_id %in% pert_ids.slct[[subgrp.nm]],]
    write.table(compound.info.slct, files.ou.drugs.txt[subgrp.nm], sep="\t", quote=F, row.names=F)
    targets=levels(factor(compound.info.slct$target))    
    targets.filt=setdiff(targets, "\"\"")
    write(targets.filt, file=files.ou.targets.txt[subgrp.nm])

    return(targets.filt)
})
names(targets.slct)=names(pert_ids.slct)

write(setdiff(targets.slct$inflamSubGrps, targets.slct$nonInflamSubGrps), file=files.ou.targets.uniq.txt["inflamSubGrps"])
write(setdiff(targets.slct$nonInflamSubGrps, targets.slct$inflamSubGrps), file=files.ou.targets.uniq.txt["nonInflamSubGrps"])

#drug to ATC code
#
drugInfo=readRDS(file.in.drugBank.RDS)

synonyms.df <- as.data.frame(drugInfo$drugs$synonyms)
atc_codes.df <- as.data.frame(drugInfo$drugs$atc_codes)

synonyms.df <- synonyms.df %>%
  mutate(synonym = tolower(synonym))  # Convert synonyms to lower case for matching
synonyms_atc_mapped.df <- synonyms.df %>%
  left_join(atc_codes.df, by = "drugbank_id")  # Join with ATC codes

drug.syn2ATC=split(synonyms_atc_mapped.df$atc_code, synonyms_atc_mapped.df$synonym)
drug.syn2ATC.filt=drug.syn2ATC[sapply(drug.syn2ATC, #keep those with one mapped ATC or mapped to multiple ATCs from the same level 1 category
FUN=function(x) 
{
    if(all(is.na(x))) return(F)
    if(length(x)>=2)
    {
        return(all(substr(x, start=1, stop=1)== substr(x[1], start=1, stop=1)))
    }
    return(T)

})]
drug.syn2ATC_L1.filt=sapply(drug.syn2ATC.filt,
FUN=function(atc)
{
    substr(atc[1], start=1, stop=1)
})
# View the result


perti_ids.slct.cmb=unlist(pert_ids.slct)
perti_ids.slct.cmb=perti_ids.slct.cmb[!duplicated(perti_ids.slct.cmb)]

compounds.info.slct=compound.info[compound.info$pert_id %in% perti_ids.slct.cmb & compound.info$target!="\"\"",]

cmap_names.slct=list()
cmap_names.slct$inflamSubGrps=compound.info[compound.info$pert_id %in% pert_ids.slct$inflamSubGrps & compound.info$target!="\"\"", "cmap_name"]
cmap_names.slct$inflamSubGrps=cmap_names.slct$inflamSubGrps[!duplicated(cmap_names.slct$inflamSubGrps)]
cmap_names.slct$nonInflamSubGrps=compound.info[compound.info$pert_id %in% pert_ids.slct$nonInflamSubGrps & compound.info$target!="\"\"", "cmap_name"]
cmap_names.slct$nonInflamSubGrps=cmap_names.slct$nonInflamSubGrps[!duplicated(cmap_names.slct$nonInflamSubGrps)]

cmap_names.slct.cmb=unlist(cmap_names.slct)
cmap_names.slct.cmb=cmap_names.slct.cmb[!duplicated(cmap_names.slct.cmb)]


# cmap_name2ATC.slct=lapply(names(cmap_name2inchi_key),
# FUN=function(cname)
# {
#     print(cname) #paste(cname, cmap_name2alias[[cname]]))
#     inchi_keys=cmap_name2inchi_key[[cname]]
#     #get compound id
#     k2info.df=do.call(rbind, lapply(inchi_keys, getCIDfromInchiKey))

#     if(! ("Title" %in% colnames(k2info.df))) return(data.frame(cname=cname, CID=NA, ATC=NA))
    
#     cid.slct=k2info.df[(!is.na(k2info.df[,"Title"])) & toupper(k2info.df[,"Title"])==toupper(cname), "CID"]
#     if(length(cid.slct)==0) return(data.frame(cname=cname, CID=NA, ATC=NA))
#     print(cid.slct)
    
#     res=getATCfromCID(cid.slct)
#     data.frame(cname=cname, 
#                 CID=res[1],
#                 ATC=res[2]

#             )

# })
# cmap_name2ATC.slct.df=do.call(rbind, cmap_name2ATC.slct)
# cmap_name2ATC.slct.df=cmap_name2ATC.slct.df[!is.na(cmap_name2ATC.slct.df$ATC),]
# rownames(cmap_name2ATC.slct.df)=cmap_name2ATC.slct.df$cname

save(deg.gsea.res.df, 
    pert_ids.slct,
    cmap_names.slct,
    cmap_names.slct.cmb,
    drug.syn2ATC.filt,
    drug.syn2ATC_L1.filt,
    #cmap_name2ATC.slct.df,
    #deg.gsea.res.slct.df,
    file=file.ou.RData)

#network
cmap_name2target=sapply(split(compounds.info.slct$target, compounds.info.slct$cmap_name),FUN=function(x) levels(factor(x)))

cnameAndtarget.df=do.call(rbind, lapply(names(cmap_name2target),
FUN=function(cnm)
{
    data.frame(cname=cnm, target=cmap_name2target[[cnm]], stringsAsFactors=F)
}))
targets=cnameAndtarget.df[,"target"]
targets=targets[!duplicated(targets)]
cname.info.df=rbind(
                data.frame(
                    node=cmap_names.slct.cmb,
                    type="drug",
                    is.inflamSubgrps=cmap_names.slct.cmb %in% cmap_names.slct$inflamSubGrps,
                    is.nonInflamSubgrps=cmap_names.slct.cmb %in% cmap_names.slct$nonInflamSubGrps,
                    ATC=drug.syn2ATC_L1.filt[cmap_names.slct.cmb],
                    stringsAsFactors=F
                    ),
                data.frame(
                    node=targets,
                    type="target",
                    is.inflamSubgrps=targets %in% targets.slct$inflamSubGrps,
                    is.nonInflamSubgrps=targets %in% targets.slct$nonInflamSubGrps,
                    ATC=NA,
                    stringsAsFactors=F))
rownames(cname.info.df)=cname.info.df$node
cname.info.df$inflam.fract= cname.info.df$is.inflamSubgrps/(cname.info.df$is.inflamSubgrps+cname.info.df$is.nonInflamSubgrps)
cname.info.df$nonInflam.fract= cname.info.df$is.nonInflamSubgrps/(cname.info.df$is.inflamSubgrps+cname.info.df$is.nonInflamSubgrps)
# cname.info.df$grp.color=NA
# cname.info.df$grp.color[cname.info.df$is.inflamSubgrps]="#FFACAC"#"lightpink1"
# cname.info.df$grp.color[cname.info.df$is.nonInflamSubgrps]="#5CACEE" #"steelblue2"
# cname.info.df$grp.color[cname.info.df$is.inflamSubgrps & cname.info.df$is.nonInflamSubgrps]="#FF67FF" #"orange"
cname.info.df$ATC.color="gray"
cname.info.df$ATC.color[grepl("^C", cname.info.df$ATC)]="#FF0000" # "red" #cardiovascular
cname.info.df$ATC.color[grepl("^L", cname.info.df$ATC)]= "#00FF08" # "green" # antineoplastic
cname.info.df$ATC.color[grepl("^N", cname.info.df$ATC)]="#0000FF" # "blue" nervous system
cname.info.df$label=gsub("-", "_", cname.info.df$node) # #for better recognized, replace "-" by "_"
#cname.info.df$ATC.color[grepl("^J", cname.info.df$ATC)]="#" # "purple"
write.table(cname.info.df, file.ou.nodes, sep="\t", quote=F, row.names=F)



drugTarget.net <- graph_from_data_frame(d=cnameAndtarget.df, vertices=cname.info.df, directed=F)
drugTarget.net=simplify(drugTarget.net, remove.multiple = TRUE, edge.attr.comb = "first")

#V(drugTarget.net)$name=gsub("-", "_", V(drugTarget.net)$name)

#
net.components=components(drugTarget.net)
net.comp2node=split(names(net.components$membership), paste0("comp", net.components$membership))
net.comp2subGrp=sapply(net.comp2node,
FUN=function(nodes)
{
    if(all(cname.info.df[nodes, "is.inflamSubgrps"]) && ! any(cname.info.df[nodes, "is.nonInflamSubgrps"])) return("inflamSubgrps")
    if(all(cname.info.df[nodes, "is.nonInflamSubgrps"]) && !any(cname.info.df[nodes, "is.inflamSubgrps"])) return("nonInflamSubgrps")
    return("bothSubgrps")
})
subGrp2nodes=sapply(split(names(net.comp2subGrp), net.comp2subGrp),
FUN=function(cmps)
{
    unlist(net.comp2node[cmps])
})
subGrp2ATCGrpCounts=sapply(subGrp2nodes,
FUN=function(nodes)
{
    atcs=cname.info.df[nodes, "ATC"]
    atcs=atcs[!is.na(atcs)]
    atcs.grp=substr(atcs, 1, 1)

    summary(factor(atcs.grp, levels=c("A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V")))
})
write.table(subGrp2ATCGrpCounts, file=file.ou.drugAndATC.tsv, sep="\t", quote=F)

#
nodes.slct=c("EGFR", "quetiapine","nifedipine", "bezafibrate")#to select typical subgroups
net.comps.slct = net.comp2node[sapply(net.comp2node, FUN=function(cmp.nodes) any(cmp.nodes %in% nodes.slct))]
net.comps.nodes.slct=unlist(net.comps.slct)
subGrp2nodes[["bothSubgrps_slctedModules"]]=net.comps.nodes.slct



for(subGrp in names(subGrp2nodes))
{
    f.o.net=paste0(file.ou.drugAndTargetsNet.pref, subGrp, ".edges.tsv")

    nodes=subGrp2nodes[[subGrp]]
    nets= cnameAndtarget.df[ cnameAndtarget.df[,1] %in% nodes & cnameAndtarget.df[,2] %in% nodes  , ]
    
    #write.table(cname.info.df[nodes,], f.o.nodes, sep="\t", quote=F, row.names=F)
    write.table(nets, f.o.net, sep="\t", quote=F, row.names=F)
}


targets.slct=V(sub.net)$name[V(sub.net)$type=="target"]
write(targets.slct, file=paste0(dir.ou, "Drug2Target.network.forSubgrps.drugbank.shared.targets.txt"))



#test whether the targets of L (ATC code) drugs enriched for wnt pathways
#check with enrichR
drugs.L = cname.info.df$node[cname.info.df$ATC.color=="#00FF08"] #green"
drugs.L.targets= unlist(ego(drugTarget.net, order = 1, nodes = drugs.L))#as.character(neighbors(drugTarget.net, v=drugs.L))
drugs.L.targets=setdiff(names(drugs.L.targets), drugs.L)
drugs.L.targets=drugs.L.targets[!duplicated(drugs.L.targets)]
file.ou.L.targets.txt=paste0(dir.ou, "L_drugs.targets.txt")
write(drugs.L.targets, file=file.ou.L.targets.txt)

save(deg.gsea.res.df, 
    pert_ids.slct,
    cmap_names.slct,
    cmap_names.slct.cmb,
    targets.slct,
    drug.syn2ATC.filt,
    drug.syn2ATC_L1.filt,
    cnameAndtarget.df,
    cname.info.df,
    drugTarget.net,
    net.comp2node,
    subGrp2nodes,
    #deg.gsea.res.slct.df,
    file=file.ou.RData)


#
pdf(file.ou.GSEA.pdf, width=8, height=8)

layout(matrix(1:4, nrow=2))
hist(deg.gsea.res.df$inflamSubGrps.upDEG.NOM_pval, main="inflamSubGrps.upDEG.NOM_pval")
hist(deg.gsea.res.df$inflamSubGrps.dnDEG.NOM_pval, main="inflamSubGrps.dnDEG.NOM_pval")

hist(deg.gsea.res.df$nonInflamSubGrps.upDEG.NOM_pval, main="nonInflamSubGrps.upDEG.NOM_pval")
hist(deg.gsea.res.df$nonInflamSubGrps.dnDEG.NOM_pval, main="nonInflamSubGrps.dnDEG.NOM_pval")


# #
# pvalcutoff.NO=sum(deg.gsea.res.df$p.filt)
# consist.NO=sum(deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES <0)
# slct.NO=sum(deg.gsea.res.df$p.filt & deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES <0)
# consist.p=fisher.test(matrix(c(slct.NO, 
#                                 consist.NO-slct.NO, 
#                                 pvalcutoff.NO-slct.NO,
#                                 nrow(deg.gsea.res.df)-sum(deg.gsea.res.df$upDEG.NES * deg.gsea.res.df$dnDEG.NES <0 | deg.gsea.res.df$p.filt)), 
#                                 ncol=2), alternative="greater")$p.value

ggplot(deg.gsea.res.df, aes(x=inflamSubGrps.upDEG.NES, y=inflamSubGrps.dnDEG.NES)) +
  geom_point(aes(color=inflamSubGrps.type, size=inflamSubGrps.p.filt, alpha=inflamSubGrps.p.filt))+
  scale_color_manual(values=c("noSig" = "gray" , "opposite"= "red", "same"="blue"))+
  scale_size_manual(values=c("FALSE"=0.2, "TRUE"=1))+
  scale_alpha_manual(values=c("FALSE"=0.2, "TRUE"=1)) +
  ggtitle(paste0("inflamSubGrps")) +
  geom_text_repel(aes(label = inflamSubGrps.label), 
                    #max.overlaps = 10,  # Adjust for fewer overlaps
                    size = 3,          # Adjust text size
                    box.padding = 0.3, # Padding around labels
                    point.padding = 0.3)+
  #ggtitle(paste0("consist effect enriched p.value: ", signif(consist.p, 3))) +
  theme_bw()




ggplot(deg.gsea.res.df, aes(x=nonInflamSubGrps.upDEG.NES, y=nonInflamSubGrps.dnDEG.NES)) +
  geom_point(aes(color=nonInflamSubGrps.type, size=nonInflamSubGrps.p.filt, alpha=nonInflamSubGrps.p.filt))+
  scale_color_manual(values=c("noSig" = "gray" , "opposite"= "red", "same"="blue"))+
  scale_size_manual(values=c("FALSE"=0.2, "TRUE"=1))+
  scale_alpha_manual(values=c("FALSE"=0.2, "TRUE"=1)) +
  ggtitle(paste0("nonInflamSubGrps")) +
  geom_text_repel(aes(label = nonInflamSubGrps.label), 
                    #max.overlaps = 10,  # Adjust for fewer overlaps
                    size = 3,          # Adjust text size
                    box.padding = 0.3, # Padding around labels
                    point.padding = 0.3)+
  theme_bw()



dev.off()


f.o.tsv=paste0(dir.ou, "GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood.tsv")
write.table(deg.gsea.res.df, f.o.tsv, sep="\t", quote=F)
