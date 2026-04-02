#V1.2
#d.cmpBulkDPacrossHMs_byLimma_V1.2_cmpNumber.R
#compare number of peak overlapped at mCRE level

#compare differential signal across marks
#at merge peak level
#grouped by sharing and modules of mARE


HMs = c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3",  "H3K27me3")

#library(ggplot2)
library(pheatmap)

DP.Q.CUTOFF=0.05
DP.V="V1.4.1"


dir.in.dp ="c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/" #"c_bulkDiffPeak_limma_V1.0_addLiMed/"  #c_bulkDiffPeak_limma_V1.4_sva_drugs/"

files.in.dp.bed  = system(paste0("find ",dir.in.dp, "*/*disease.Q", DP.Q.CUTOFF, "*.bed"), intern=T)
names(files.in.dp.bed)=sapply(files.in.dp.bed,
FUN=function(f.i)
{
  x=sub(paste0(".diffPeak.diseasedisease.Q", DP.Q.CUTOFF), "", basename(f.i))
  x=sub(".bed", "", x)
})

file.in.5HM.mARE.RData = "../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2HMPeak.hg19.RData"

#file.in.sharing.RData = "../peakVariationAcrossTiss/b_mARE_annotation/mARE_annotation.RData"
#file.in.mARE.module.RData = "../peakVariationAcrossTiss/a_3_AREModulesInEpimap/H3K27ac/5HMMergedPeaks.inEpimap.H3K27acSignal.cluster.RData"

dir.ou = "d_cmpBulkDPAcrossHMs_limma_V1.2_CRE.number.ovlp/"
dir.create(dir.ou, recursive = T, showWarnings = F)
#file.ou.ovlp.pref = paste0(dir.ou, "peaks.ovlp.hg19.") 
file.ou.ovlp.RData = paste0(dir.ou, "DP.", DP.V, ".Q", DP.Q.CUTOFF, ".sharedAcrossHMs.hg19.RData") 
file.ou.ovlp.pdf = paste0(dir.ou, "DP.", DP.V, ".Q", DP.Q.CUTOFF, ".sharedAcrossHMs.hg19.pdf") 


#DP signal
hm.dps = lapply(files.in.dp.bed,
FUN=function(f.i)
{
  d=read.table(f.i, header = F, sep="\t", row.names = NULL, stringsAsFactors = F)
  paste0(d[,1], ":", d[,2], "-", d[,3])
})



#merged AREs
load(file.in.5HM.mARE.RData)

hm.dps.mARE = lapply(names(hm.dps),
FUN=function(cond)
{
  hm=sub(".dw|.up", "", cond)
  print(hm)
  levels(factor(hm.peak2mergedPeak[[hm]][hm.dps[[cond]]]))  
  
})
names(hm.dps.mARE) = names(hm.dps)

#
hm.dps.Jaccard=sapply(hm.dps.mARE,
FUN=function(mARE.1)
{
  sapply(hm.dps.mARE,
  FUN=function(mARE.2)
  {
    length(intersect(mARE.1, mARE.2))/length(union(mARE.1, mARE.2))
  })
})
diag(hm.dps.Jaccard)=NA

#overlap number
hm.dps.ovlp=sapply(hm.dps.mARE,
FUN=function(mARE.1)
{
  sapply(hm.dps.mARE,
  FUN=function(mARE.2)
  {
   length(intersect(mARE.1, mARE.2))
  })
})

#
hm.dps.ovlp.p=matrix(1, ncol=length(hm.dps.mARE),nrow=length(hm.dps.mARE))
hm.dps.ovlp.OR=matrix(1, ncol=length(hm.dps.mARE),nrow=length(hm.dps.mARE))
rownames(hm.dps.ovlp.p) <- rownames(hm.dps.ovlp.OR) <- names(hm.dps.mARE)
colnames(hm.dps.ovlp.p) <- colnames(hm.dps.ovlp.OR) <- names(hm.dps.mARE)

for(nm1 in names(hm.dps.mARE))
{
  hm1=sub(".dw|.up", "", nm1)
  for(nm2 in names(hm.dps.mARE))
  {
    hm2=sub(".dw|.up", "", nm2)
    mAREs.ovlp=intersect(hm.peak2mergedPeak[[hm1]],hm.peak2mergedPeak[[hm2]])
    mAREs.ovlp=levels(factor(mAREs.ovlp))

    mAREs.DP.ovlp= length(intersect(hm.dps.mARE[[nm1]], hm.dps.mARE[[nm2]]))
    mAREs.DP.hm1 = length((intersect(setdiff(hm.dps.mARE[[nm1]], hm.dps.mARE[[nm2]]), mAREs.ovlp)))
    mAREs.DP.hm2 = length(intersect(setdiff(hm.dps.mARE[[nm2]], hm.dps.mARE[[nm1]]), mAREs.ovlp))
    mAREs.rest=length(setdiff(mAREs.ovlp, union(hm.dps.mARE[[nm1]], hm.dps.mARE[[nm2]])))

    fish.res=fisher.test(matrix(c(mAREs.DP.ovlp, mAREs.DP.hm1, mAREs.DP.hm2, mAREs.rest), ncol=2), 
                        alternative="greater")
    print(paste(nm1, nm2, fish.res$estimate))
    hm.dps.ovlp.OR[nm1,nm2]=fish.res$estimate
    hm.dps.ovlp.p[nm1,nm2]=fish.res$p.value
    
  }
}
diag(hm.dps.ovlp.OR)=NA
diag(hm.dps.ovlp.p)=NA
hm.dps.ovlp.adjP= apply(hm.dps.ovlp.p, 2, p.adjust, method="BH")

#matrix(p.adjust(as.vector(hm.dps.ovlp.p), method="BH"), ncol=length(hm.dps.mARE))

# colnames(hm.dps.Jaccard) <- rownames(hm.dps.Jaccard) <- 
save(hm.dps.Jaccard,
    hm.dps.ovlp,
    hm.dps.mARE,

    hm.dps.ovlp.p,
    hm.dps.ovlp.OR,
    file=file.ou.ovlp.RData)



pdf(file.ou.ovlp.pdf, height=8, width=8)
# pheatmap(hm.dps.Jaccard,
#          cluster_row=F,
#          cluster_col=F,
#          labels_row=paste0(names(hm.dps), "(N=", sapply(hm.dps, length),")"),
#          display_numbers=hm.dps.ovlp)

# adjP.level=matrix("", ncol=ncol(hm.dps.ovlp.adjP), nrow=nrow(hm.dps.ovlp.adjP))
# adjP.level[hm.dps.ovlp.adjP<0.05]="*"

hm.dps.ovlp.display=hm.dps.ovlp
hm.dps.ovlp.display[hm.dps.ovlp.adjP>0.05]=""
pheatmap(hm.dps.ovlp.OR,
         cluster_row=F,
         cluster_col=F,
         labels_row=paste0(names(hm.dps), "(N=", sapply(hm.dps, length),")"),
         display_numbers=hm.dps.ovlp.display
        )
dev.off()



f.o.tsv=paste0(dir.ou, "DP.", DP.V, ".Q", DP.Q.CUTOFF, ".sharedAcrossHMs.hg19.tsv") 
write.table(hm.dps.ovlp.OR, file=f.o.tsv, sep="\t", quote=F)

