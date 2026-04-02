
#V1.1
#c.gneighb.VisGWASANDBulkAnnot_linkFromGenetics_V1.1_3HM.R
# only use H3K4me1, H3K27ac, and H3K27me3 
#which regulate gene expression

#links are between AREs and genes based on genetic evidence
#

# visualize information from different layer gene by gene by signal tracks
#Bulk annotation

#GWAS signal
#haQTL-ARE from Marks
#ARE MR and coloc results
#DP signal 

#use genes with promoter region overlapping with either TssA and TssFlank
#consider H3K27ac peaks associated within gene neighborhood
#genomic region is in format as "chr2\t1\t223"
#integrate signals in a gene centered manner
#a) how genomic regions is associated with genes
#b) how the haQTL, differential peak signal is related to each peak, which is in turn related to each gene
#check haQTL effect size VS GWAS for peaks for each gene

#prmotoer regions are defined as up/downstream 2kb exntended from tss
#gene neighborhood are defined as hic regions overlap with promoter regions, and Hi-C one step neighbors of these regions, if no overlapping found,  promoter itself is defined as gene neighborhood


# args = commandArgs(trailingOnly=TRUE)
# HM.LINK = args[1] 


#
library(Sushi)
library(MendelianRandomization)
source("~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R")


# DISEASES = c("CD", "UC", "RA")
# DIS2COLOR = c("brown", "purple", "orange")
# names(DIS2COLOR) = DISEASES
PRMOTER.WINDOW = 2000
PK.WIND = 10000
GWAS.PVAL.CUTOFF = 1e-6
GWAS.Z.CUTOFF = abs(qnorm(GWAS.PVAL.CUTOFF/2))
DP.Q.CUTOFF=0.05
COLOC.ABF.H4.CUTOFF=0.5
MR.P.CUTOFF=0.01
hQTL.P.CUTOFF = "emp.p.fdr.cutoff0.2"
gLinks.GSP.cutoffs = c("coloc.H4.ABF" =0.5,
                "coloc.H4vsH3.ratio" =5, 
                "MR.qval" =-0.05, 
                "FMeQTL2ARE.DistanceInv" = 1/2000)


HMs= c("H3K27ac","H3K4me1", "H3K27me3") #"H3K36me3" "H3K4me3", 

#HAQTL.COLs = c(rgb(1, 0, 0, 0.5), rgb(1, 1, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 1, 1, 0.5), rgb(0, 0, 1, 0.5), rgb(1, 0, 1, 0.5))

files.in.gLink.RData = paste0("b_summary_Links.coloc.MR.FMeQTL_V1.3_gARE/Blood_", HMs, ".linkFromColoc.MR.RData")
names(files.in.gLink.RData)=HMs
#file.in.g2pk.RData = paste0("./a_gAREsNeareGene_emp.p.fdr.cutoff0.2/Blood_1000k/", HM.LINK, ".eGene2gARE.hg38.RData")

file.in.eGene.info = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/Whole_Blood.v8.egenes.txt.gz")

# files.in.TssA = c(E047= "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE047TssAFlank.bed",
#                   E048="~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE048TssAFlank.bed",
#                   E039 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE039TssAFlank.bed",
#                   E040 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE040TssAFlank.bed",
#                   E029 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE029TssAFlank.bed",
#                   BPNeutrop = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE030TssAFlank.bed",#
#                   E032 = "~/lhou.compbio/data/Roadmap/15CoreChromStates127epigen_annot/gencodeV26.Promoter_ovlpWithE032TssAFlank.bed"
#                   )
#file.in.TSS = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.TSSUp0Down1bp.geneName.ENSG.bed"

files.in.pk.RData= paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".MayoWithRef.heightsMean.depthGCNormed.RData")
names(files.in.pk.RData)= HMs

files.in.diffPeak.RData = paste0("../peakMergeNormal/c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/", HMs, "/", HMs, ".diffPeak.RData")
names(files.in.diffPeak.RData) = HMs

files.in.haQTLPeak.hg38.RData = paste("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData", sep="")
names(files.in.haQTLPeak.hg38.RData)= HMs

files.in.haQTLPeak.coloc.hg38.RData=paste0("../haQTLvsGWAS/a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/PGC.BP/", HMs, "_100k/", HMs, ".bulkhQTLvsGWAS.RData")
names(files.in.haQTLPeak.coloc.hg38.RData) = HMs

files.in.haQTLPeak.MR.hg38.RData=paste0("../haQTLvsGWAS/a_bulkhaQTLvsGWAS_byMR_Egger_V1.1_allGWAS/PGC.BP/", HMs, "_100k/", HMs, ".bulkHaQTLvsGWAS.RData")
names(files.in.haQTLPeak.MR.hg38.RData) = HMs

file.in.BP.gz = "~/lhou.compbio/data/GWAS/pgc_bipolar/pgc_bip_hg19.sorted.bed.gz"

file.in.geneStructure.RData = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.annotation.forSushi_exonAndUTR.RData"

dir.tmp = "~/hptmp/MayoBipolar/c_geneNeighb_geneticLink/"
dir.create(dir.tmp, recursive = T, showWarnings = F)

#file.tmp.cellTypePeak = paste(dir.tmp, REF, ".hg19.nodup.filt.narrowPeak.merged.bed", sep="")
#file.tmp.TssA.bed = paste(dir.tmp, REF, ".TSSA.bed", sep="")
#file.tmp.hicRegs.bed = paste(dir.tmp, "hicRegions.bed", sep="")
file.tmp.gRegion.bed= paste(dir.tmp, "gRegion.bed", sep="")
file.tmp.gFilt.bed= paste(dir.tmp, "gFilt.bed", sep="")
file.tmp.genesStructure.bed = paste(dir.tmp, "genesStructure.bed", sep="")
file.tmp.gStructure.bed = paste(dir.tmp, "gStructure.bed", sep="")
file.tmp.pk.bed = paste(dir.tmp, "pk.bed", sep="")
file.tmp.BP.hits.bed = paste0(dir.tmp, "BP.hits.hg19.bed")
file.tmp.pk.ovlp.BP.hits.bed= paste0(dir.tmp, "pkOVLP.BP.hits.hg19.bed")

  
dir.ou="c_geneNeighb_visGWASANDBulkAnnot_byGene_GeneticLink_V1.1_3HM/"
dir.create(dir.ou, showWarnings = F, recursive = T)

#file.out.pref = paste(dir.ou, "geneNB", sep="")
file.ou.gNeighbor.link.RData = paste(dir.ou, "3HM.gLinks.RData", sep="")
file.ou.gNeighbor.annot.RData = paste(dir.ou, "3HM.geneNB.peaks.annot_gLinks.RData", sep="")
#file.ou.pdf = paste(dir.ou, "geneNB.GWASp", GWAS.PVAL.CUTOFF, ".haQTLDiffPks.", HM.LINK, "_geneticlink.pdf", sep="")
#file.ou.haQTLDiffPks2GWAS.RData = paste(file.out.pref, ".haQTLDiffPks2GWAS.RData", sep="")
#file.ou.geneNeighb = paste(file.out.pref, ".txt", sep="")


################################################################################
load(file.in.geneStructure.RData)
write.table(transcripts.bed, file=file.tmp.genesStructure.bed, sep="\t", col.names = F, row.names = F, quote=F)

###############################################################################
# gene2hicReg = list() #gene 2 neighbor regions based on hiC and promoters
# 
# hicReg2Annot =  matrix()

peaks2haQTL = list() #peak 2 haQTL, snpid and zscore
peaks2diffZ = list() #peaks 2 differential Zscore for each disease
g.neighb2peaks = list() #neighbor mapped 2 peaks
gene2peaks = list() #gene 2 peaks associated with its neighbors

#######################################################
HM.hg38tohg19= lapply(HMs,
FUN=function(HM)
{
  load(files.in.pk.RData[HM])
  return(pks.hg382hg19.uniq)
})
names(HM.hg38tohg19)=HMs

#annotate###################################################

HM.DP.score.hg19=lapply(HMs,
FUN=function(HM)
{
  load(files.in.diffPeak.RData[HM])
  
  dp.signals=cbind(log2FoldChange= res$diseasedisease$logFC, pvalue=res$diseasedisease$P.Value, padj = res$diseasedisease$adj.P.Val)
  rownames(dp.signals) = rownames(res$diseasedisease)
  
  return(dp.signals)
})
names(HM.DP.score.hg19)=HMs

HM.PK.slct.hg19 = lapply(HMs,
FUN=function(HM)
{
  rownames(HM.DP.score.hg19[[HM]])
})
names(HM.PK.slct.hg19)=HMs

HM.DP.hg19= lapply(HM.DP.score.hg19,
FUN=function(dp.signals)
{
  dps=rownames(dp.signals)[dp.signals[,"padj"]<=DP.Q.CUTOFF]
  # dp2direction=rep(1, length(dps))
  # names(dp2direction)=dps
  # dp2direction[dp.signals[dps,"log2FoldChange"]<0]=-1

  return(dp.signals[dps,"log2FoldChange"])

})
names(HM.DP.hg19) = names(HM.DP.score.hg19)


HM.hQTLPK.hg19 = lapply(HMs,
FUN=function(HM)
{
  load(files.in.haQTLPeak.hg38.RData[HM])
  pks.hg38=hQTLPeaks.list[[hQTL.P.CUTOFF]]
  #buf=read.table(f.i.bed, sep="\t", row.names = NULL, header = F, stringsAsFactors = F)
  #= paste0(buf[,1], ":", buf[,2], "-", buf[,3])
  pks.hg19 = HM.hg38tohg19[[HM]][pks.hg38]
  pks.hg19.slct = intersect(HM.PK.slct.hg19[[HM]], pks.hg19)
  return(pks.hg19.slct)
})
names(HM.hQTLPK.hg19) = HMs

HM.hQTLPK.coloc.hg19 = lapply(HMs,
FUN=function(HM)
{
  f.i.RData=files.in.haQTLPeak.coloc.hg38.RData[HM]
  load(f.i.RData)
  
  pks.hg38=names(hQTLPeaks.QTLvsGWAS)[sapply(hQTLPeaks.QTLvsGWAS, FUN=function(x){if(is.na(x)) return(NA);x$coloc$summary["PP.H4.abf"]>=COLOC.ABF.H4.CUTOFF})]
  pks.hg38=pks.hg38[!is.na(pks.hg38)]
  pks.hg19 = HM.hg38tohg19[[HM]][pks.hg38]
  pks.hg19.slct = intersect(HM.PK.slct.hg19[[HM]], pks.hg19)
  return(pks.hg19.slct)
  
})
names(HM.hQTLPK.coloc.hg19) = HMs

HM.hQTLPK.MR.hg19 = lapply(HMs,
FUN=function(HM)
{
  f.i.RData=files.in.haQTLPeak.MR.hg38.RData[HM]
  load(f.i.RData)

  pk2MR.est.hg38=sapply(hQTLPeaks.QTLvsGWAS, 
  FUN=function(x)
  {
    if(is.na(x)) return(NA)
    if(is.na(x$MR@Causal.pval)) return(NA)
    if(x$MR@Causal.pval<=MR.P.CUTOFF) return(x$MR@Estimate)
    return(NA)
  })

  pk2MR.est.hg38=pk2MR.est.hg38[!is.na(pk2MR.est.hg38)]
  pks.hg38=names(pk2MR.est.hg38)

  pks.hg19 = HM.hg38tohg19[[HM]][pks.hg38]
  pks.hg19.slct = intersect(HM.PK.slct.hg19[[HM]], pks.hg19)

  pk2MR.est.hg19=pk2MR.est.hg38[HM.hg38tohg19[[HM]][pks.hg38] %in% pks.hg19.slct]
  names(pk2MR.est.hg19)= HM.hg38tohg19[[HM]][names(pk2MR.est.hg19)]
  return(pk2MR.est.hg19)
})
names(HM.hQTLPK.MR.hg19)=HMs

#
cmd = paste0("zcat ", file.in.BP.gz, "|awk '($10<", GWAS.PVAL.CUTOFF, "){print }'>", file.tmp.BP.hits.bed)
system(cmd)
HM.PK.GWAS.hg19 = lapply(HMs,
FUN=function(HM)
{
  print(HM)
  write(paste(gsub(":|-", "\t", HM.PK.slct.hg19[[HM]]), HM.PK.slct.hg19[[HM]], sep="\t"), file.tmp.pk.bed)
  cmd = paste0("bedtools window -w ", PK.WIND, " -a ", file.tmp.pk.bed, " -b ", file.tmp.BP.hits.bed, ">", file.tmp.pk.ovlp.BP.hits.bed)
  system(cmd)
  buf = read.table(file.tmp.pk.ovlp.BP.hits.bed, "\t", header=F, row.names = NULL)
  peaks.ovlpGWAS = levels(factor(buf[,4]))


})
names(HM.PK.GWAS.hg19) = HMs


HM.PK.hg19.annot = lapply(HMs,
FUN=function(HM)
{
  pk.annot.matrix = matrix(0, ncol=5, nrow=length(HM.PK.slct.hg19[[HM]]))
  colnames(pk.annot.matrix)=c("DP.log2FC", "hQTLPk", "coloc", "MR.est", "GWAS")
  rownames(pk.annot.matrix) = HM.PK.slct.hg19[[HM]]

  pk.annot.matrix[names(HM.DP.hg19[[HM]]), "DP.log2FC"]=HM.DP.hg19[[HM]] #pk.annot.matrix[HM.DP.hg19[[HM]], "DP"]=1
  pk.annot.matrix[HM.hQTLPK.hg19[[HM]], "hQTLPk"]=1
  pk.annot.matrix[HM.hQTLPK.coloc.hg19[[HM]], "coloc"]=1
  pk.annot.matrix[names(HM.hQTLPK.MR.hg19[[HM]]), "MR.est"]=HM.hQTLPK.MR.hg19[[HM]]
  pk.annot.matrix[HM.PK.GWAS.hg19[[HM]], "GWAS"]=1


  pk.annot.matrix.filt=pk.annot.matrix[apply(pk.annot.matrix, 1, FUN=function(x) any(x!=0)),]


  return(pk.annot.matrix.filt)

})
names(HM.PK.hg19.annot)=HMs



save(HM.PK.hg19.annot,
     #hic.reg2HM.PK.annot.matrix,
     file=file.ou.gNeighbor.annot.RData)
# save(HM.hg38tohg19,
#      HM.DP.hg19,
#      HM.DP.score.hg19,
#      HM.hQTLPK.hg19,
#      HM.hQTLPK.coloc.hg19,
#      HM.hQTLPK.MR.hg19,
#      HM.PK.hg19.annot,
#      file=file.ou.gNeighbor.RData
#      )
################################################################################
#
eGene.info = read.table(gzfile(file.in.eGene.info), header = T, sep="\t", row.names = 1)
gESGid2gName = eGene.info$gene_name
names(gESGid2gName) = rownames(eGene.info)

#

print("generating gene neighborhood") ################################################################

hm.gLinks.hg19=lapply(HMs,
FUN=function(hm)
{
  load(files.in.gLink.RData[hm])
  gAndPkAndlinkScores = cbind(data.frame(t(sapply(rownames(scores.all), FUN=function(x){unlist(strsplit(x, split=";"))}))),
                              scores.all)
  gAndPkAndlinkScores.filt= gAndPkAndlinkScores[gAndPkAndlinkScores$coloc.H4.ABF >=gLinks.GSP.cutoffs["coloc.H4.ABF"] |
                                                  gAndPkAndlinkScores$coloc.H4vsH3.ratio >=gLinks.GSP.cutoffs["coloc.H4vsH3.ratio"] |
                                                  gAndPkAndlinkScores$MR.qval >=gLinks.GSP.cutoffs["MR.qval"] |
                                                  gAndPkAndlinkScores$FMeQTL2ARE.DistanceInv >=gLinks.GSP.cutoffs["FMeQTL2ARE.DistanceInv"],]
  print(nrow(gAndPkAndlinkScores.filt))
  gene2Neighb.scores.filt.hg19 = lapply(split(gAndPkAndlinkScores.filt[,2], gAndPkAndlinkScores.filt[,1]),
  FUN=function(pks.hg38)
  {
    HM.hg38tohg19[[hm]][pks.hg38]
  })
 
})
names(hm.gLinks.hg19) = HMs

for(hm in names(hm.gLinks.hg19))
{
  geneAndCRE.df = do.call(rbind, lapply(names(hm.gLinks.hg19[[hm]]), 
  FUN=function(g)
  {
    data.frame(gid=g, 
                gene=gESGid2gName[g],
                CRE=hm.gLinks.hg19[[hm]][[g]],
                HM=hm, 
                stringsAsFactors=F)
  }))

  file.ou.geneCRELink.tsv = paste(dir.ou, hm, ".geneCRELink.gLinks.hg19.tsv", sep="")
  write.table(geneAndCRE.df, file=file.ou.geneCRELink.tsv, sep="\t", quote=F, row.names=F)
}

save(hm.gLinks.hg19,
     file=file.ou.gNeighbor.link.RData)

#######################################################
# #gene neighbor filter
# eGene.info = read.table(gzfile(file.in.eGene.info), header = T, sep="\t", row.names = 1)
# gESGid2gName = eGene.info$gene_name
# names(gESGid2gName) = rownames(eGene.info)

buf= split(transcripts.bed, transcripts.bed$gene)
gName2range.hg19=lapply(buf,
FUN=function(x)
{
  if(any(x$chrom!=x$chrom[1]))
  {
    return(NA)
  }
  return(data.frame(chr=x$chrom[1],
             start=min(x$start),
             end=max(x$stop),
             strand = x$strand[1]))
})
gName2range.hg19=do.call(rbind, gName2range.hg19[!is.na(gName2range.hg19)])

hm.gNB2annotSum.list = lapply(names(hm.gLinks.hg19),
FUN=function(hm)
{
  t(sapply(hm.gLinks.hg19[[hm]],
  FUN=function(pks)
  {
    pks.ovlp = intersect(pks, rownames(HM.PK.hg19.annot[[hm]]))
    if(length(pks.ovlp)==0)
    {
      return(rep(0, ncol(HM.PK.hg19.annot[[hm]])))
    }
    apply(rbind(HM.PK.hg19.annot[[hm]][pks.ovlp, ]), 2 ,sum)
  }))
})
names(hm.gNB2annotSum.list) = names(hm.gLinks.hg19)

genes.all= levels(factor(unlist(lapply(hm.gNB2annotSum.list, rownames))))
gNB2annotSum=t(sapply(genes.all,
FUN=function(g)
{
  g.hm.annotSum=lapply(names(hm.gLinks.hg19),
  FUN=function(hm)
  {
    if(g %in% rownames(hm.gNB2annotSum.list[[hm]]))
    {
      return(hm.gNB2annotSum.list[[hm]][g,])
    }else
    {
      return(rep(0, ncol(hm.gNB2annotSum.list[[hm]])))
    }
  })
  names(g.hm.annotSum) = names(hm.gLinks.hg19)
  
  unlist(g.hm.annotSum)
}))

gNB2annotSumHM = cbind(apply(gNB2annotSum[, grepl(".DP", colnames(gNB2annotSum))], 1, sum),
                       apply(gNB2annotSum[, grepl(".hQTLPk", colnames(gNB2annotSum))], 1, sum),
                       apply(gNB2annotSum[, grepl(".coloc", colnames(gNB2annotSum))], 1, sum),
                       apply(gNB2annotSum[, grepl(".MR", colnames(gNB2annotSum))], 1, sum),
                       apply(gNB2annotSum[, grepl(".GWAS", colnames(gNB2annotSum))], 1, sum)
)
colnames(gNB2annotSumHM)=c("DP", "hQTLPk", "coloc", "MR", "GWAS")




genes.filt = rownames(gNB2annotSumHM)[gNB2annotSumHM[,"DP"]>=1 &
                                        #gNB2annotSumHM[,"hQTLPk"]>=1 &
                                        gNB2annotSumHM[,"GWAS"]>=1 &
                                        (gNB2annotSumHM[,"coloc"]>=1 &gNB2annotSumHM[,"MR"]>=1 )]


#enrichment test
#test enrichment
gs.bg = rownames(gNB2annotSumHM)
genes.DiffReg.NO=sum(gNB2annotSumHM[gs.bg,"DP"]>=1)
genes.Genetics.NO= sum(gNB2annotSumHM[gs.bg,"hQTLPk"]>=1 &
                         gNB2annotSumHM[gs.bg,"GWAS"]>=1 &
                         (gNB2annotSumHM[gs.bg,"coloc"]>=1  & gNB2annotSumHM[gs.bg,"MR"]>=1 ))
genes.ovlp.NO=length(genes.filt)
fish.p=fisher.test(matrix(c(genes.ovlp.NO,
                            genes.DiffReg.NO-genes.ovlp.NO,
                            genes.Genetics.NO-genes.ovlp.NO,
                            nrow(gNB2annotSumHM)-genes.DiffReg.NO-genes.Genetics.NO+genes.ovlp.NO), ncol=2))$p.val
print(fish.p)


#genes.filt = rownames(gNB2annotSumHM)[apply(gNB2annotSumHM>=1, 1, sum)>=2]
#genes.filt = rownames(gNB2annotSum)[apply(gNB2annotSum>=1, 1, sum)>2]

gNB2range.hg19=sapply(genes.all,
FUN=function(g)
{
  g.range = gName2range.hg19[gESGid2gName[g],]
    
  pks = unlist(lapply(hm.gLinks.hg19,
  FUN=function(gLinks)
  {
    if(g %in% names(gLinks)) return(gLinks[[g]])
    return(NULL)
  }))
    
  buf= sapply(pks, FUN=function(x) unlist(strsplit(x, split=":|-")))
  buf.chr= buf[1,]
  buf.start= min(c(as.numeric(buf[2,]), g.range$start))
  buf.end= max(c(as.numeric(buf[3,]), g.range$end))

  if(all(buf.chr==buf.chr[1]))
  {
    return(paste0(buf.chr[1], ":", buf.start, "-", buf.end))
  }else
  {
    return(NA)
  }
})

genes.filt = genes.filt[!is.na(gNB2range.hg19[genes.filt])]
write.table(cbind(gsub(":|-", "\t", gNB2range.hg19[genes.filt]), genes.filt), file=file.tmp.gFilt.bed, sep="\t", row.names = F, col.names = F, quote=F)
genes.filt.sorted = sapply(system(paste0("sortBed -i ", file.tmp.gFilt.bed), intern = T),
FUN=function(x)
{
  unlist(strsplit(x, split="\t"))[4]
})

# gNB2GWAS.topPval =sapply(genes.filt,
# FUN=function(g)
# {
#   print(g)
#   range=gNB2range[[g]]
#   if(is.na(range))
#   {
#     return(NA)
#   }
#   gwas.tab = seqminer::tabix.read.table(file.in.BP.gz, range)
#   min(gwas.tab$V10)
# })
# names(gNB2GWAS.topPval) =genes.filt

# genes.filt = rownames(gNB2annotSumHM)[apply(gNB2annotSumHM>=1, 1, sum)>=2 & gNB2GWAS.topPval<=GWAS.PVAL.CUTOFF]
# genes.filt = genes.filt[!is.na(genes.filt)]

save(HM.PK.hg19.annot,
     gESGid2gName,
     gName2range.hg19,
     gNB2annotSum,
     gNB2annotSumHM,
     gNB2range.hg19,
     #gNB2GWAS.topPval,
     genes.filt.sorted,
     file=file.ou.gNeighbor.annot.RData)






print("tracks for each gene neighborhood")##########################################
load(file.ou.gNeighbor.link.RData)
# load(file.ou.gNeighbor.annot.RData)
# #FUNCTION##############################################################################
# getChrFromRegions = function(regions) #vector of regions with format "chr\tstart\tend"
# {
#   chrs = sapply(regions, FUN=function(x){unlist(strsplit(x, split="\t"))[1]})
#   
#   return(levels(factor(chrs)))
# }
# 
# # getHiCLinks = function(regions, hic.reg2neighb)#vector of regions with format "chr\tstart\tend"
# # {
# #   nb.links = c()
# #   for(nb in intersect(regions, names(hic.reg2neighb)))
# #   {
# #     nb.partners = hic.reg2neighb[[nb]]
# #     for(nb.p in intersect(nb.partners, nbs))
# #     {
# #       nb.links=c(nb.links, paste(nb, nb.p, NA, 1, ".", ".", 1, sep="\t"))
# #       if(nb > nb.p)
# #       {
# #         nb.links=c(nb.links, paste(nb.p, nb, NA, 1, ".", ".", 1, sep="\t"))
# #       }
# #     }
# #   }
# #   nb.links = levels(factor(nb.links))
# #   return(nb.links)
# # }
# 
# getGenetLinksBEDPE.hg19=function(g, pks)
# {
#   #g is ensemb gid
#   pk.loc  = rbind(t(sapply(pks, FUN=function(x)unlist(strsplit(x, split=":|-", perl=T)))))
#   
#   gID.range= gName2range[gESGid2gName[g],]
#   link.type = gsub("NA", "", paste(c("", "coloc")[(gene2Neighb.scores.filt.hg19[[g]][pks,"coloc.H4.ABF"]>=COLOC.ABF.H4.CUTOFF)+1],
#                                           c("", "MR")[(gene2Neighb.scores.filt.hg19[[g]][pks,"MR.log10pval"]>= -log10(MR.P.CUTOFF)) +1],
#                                           sep=";"))
#   if(gID.range$strand=="+")
#   {
#     start1=gID.range$start-1
#     end1 = gID.range$start
#   }else
#   {
#     start1=gID.range$end
#     end1 = gID.range$end+1
#   }
#   nb.links.bedpe = data.frame(chrom1 = gID.range$chr,
#                               start1 = start1,
#                               end1 = end1,
#                               chrom2=pk.loc[,1],
#                               start2=as.numeric(pk.loc[,2]),
#                               end2=as.numeric(pk.loc[,3]),
#                               name= link.type,
#                               score = 1,
#                               strand1 = ".",
#                               strand2 =".",
#                               samplenumber=1
#                               )
#   return(nb.links.bedpe)
# }
# 
# getNBGenesBED = function(gName, g.chr, chr.start, chr.end)
# {
#   
#   write(paste(g.chr, "\t", chr.start, "\t", chr.end, sep=""), file.tmp.gRegion.bed)
#   cmd =paste( "intersectBed -a ", file.tmp.genesStructure.bed, " -b ", file.tmp.gRegion.bed," -u>", file.tmp.gStructure.bed, sep="")
#   system(cmd)
#   buf = read.table(file.tmp.gStructure.bed, head=F, sep="\t", row.names = NULL, stringsAsFactors = F)
#   colnames(buf) = c("chrom", "start", "stop", "gene", "gene.type", "transcript", "score", "strand", "type")
#   buf = buf[buf$gene.type=="protein_coding" | buf$gene==gName ,]#only keep protein coding genes and the gene we are focusing on now
#   buf$strand[as.character(buf$strand)=="+"] = 1
#   buf$strand[as.character(buf$strand)=="-"] = -1
#   buf$strand = as.numeric(buf$strand)
#   
#   NB.genes.bed = list()
#   NB.genes.bed$g = buf[buf$gene==gName, ]
#   NB.genes.bed$others = buf[buf$gene!=gName, ]
#   
#   return(NB.genes.bed)
# }
# 
# # gethaQTLFiltBedgraph = function(peaks2haQTL, g.pks, g.chr) #sorted
# # {
# #   haQTL.posAndT.bed = data.frame()
# #   for(i in 1:length(g.pks))
# #   {
# #     for(posID in names(peaks2haQTL[[g.pks[i]]]))
# #     {
# #       t= abs(peaks2haQTL[[g.pks[i]]][posID])
# #       pos = as.numeric(unlist(strsplit(posID, split="_"))[2])
# #       chr = unlist(strsplit(posID, split="_"))[1]
# #       haQTL.posAndT.bed = rbind(haQTL.posAndT.bed, data.frame(chrom = chr, start = pos-1, end = pos, value = t,stringsAsFactors = F))  
# #     }
# #     
# #   }
# #   haQTL.posAndT.bed =  haQTL.posAndT.bed[ haQTL.posAndT.bed$chrom==g.chr, ]
# #   haQTL.posAndT.bed =  haQTL.posAndT.bed[order(haQTL.posAndT.bed$start, decreasing = F),]
# #   return(haQTL.posAndT.bed)
# # }
# 
# 
# getGWASFiltBedgraph = function(g.haQTL.gwas, dis, g.chr)
# {
#   haQTL.posAndGWAS.bed = data.frame(stringsAsFactors = F)
#   for(posID in names(g.haQTL.gwas[[dis]]))
#   {
#     
#     z= abs(g.haQTL.gwas[[dis]][posID])
#     pos = as.numeric(unlist(strsplit(posID, split="_"))[2])
#     chr = unlist(strsplit(posID, split="_"))[1]
#     haQTL.posAndGWAS.bed = rbind(haQTL.posAndGWAS.bed, data.frame(chrom = chr, start = pos-1, end = pos, value = z,stringsAsFactors = F))
#     
#   }
#   haQTL.posAndGWAS.bed = haQTL.posAndGWAS.bed[haQTL.posAndGWAS.bed$chrom==g.chr,]
#   if(nrow(haQTL.posAndGWAS.bed)>=1)
#   {
#     haQTL.posAndGWAS.bed =  haQTL.posAndGWAS.bed[order(haQTL.posAndGWAS.bed$start, decreasing = F),]
#   }
#   return(haQTL.posAndGWAS.bed)
# }
# 
# 
# 
# #
# pdf(file.ou.pdf, width =25, height=12)
# layout(matrix(1:5, ncol=1))
# par(mar=c(3,15,2,3))
# for(g in genes.filt.sorted)
# {
#   print(g)
#   g.range = gNB2range.hg19[[g]]
#   if(is.na(g.range)) next
#     
#   g.chr = getChrFromRegions(gsub(":|-", "\t", g.range))
#   if(length(g.chr)>1)
#   {
#     print(paste(g, "mapping to multiple chroms"))
#     next
#   }
#   
#   #filt by gwas
#   gwas.tab = seqminer::tabix.read.table(file.in.BP.gz, g.range) 
#   if(length(gwas.tab)==0) next
#   gwas.bedgraph = data.frame(chr= gwas.tab$V1,
#                              start = gwas.tab$V2,
#                              end = gwas.tab$V3,
#                              value = -log10(gwas.tab$V10) 
#                              )
#   if(all(gwas.bedgraph$value <= -log10(GWAS.PVAL.CUTOFF))) next
#   
#   
#   #range
#   pks = rownames(gene2Neighb.scores.filt.hg19[[g]])
# 
#   nb.links.bedpe = getGenetLinksBEDPE.hg19(g, pks)
#   if(nrow(nb.links.bedpe)<1)
#   {
#     print(paste(g, "no links around on the same chrom"))
#     next
#   }
#   print(g)
#   
#   chr.start = max(1,min(c(nb.links.bedpe$start1, nb.links.bedpe$start2))-5000)
#   chr.end  = max(c(nb.links.bedpe$end1, nb.links.bedpe$end2))+5000
#   
#   
#   
#   #axis(side=2,las=2,tcl=.2)
#   #mtext("",side=2,line=1.75,cex=.75,font=2)
# 
#   #neighbor regions
#   # nbs.bed = data.frame(t(sapply(nbs, FUN=function(x){unlist(strsplit(x, split="\t"))})), score = 1 , strand = ".",stringsAsFactors = F)
#   # colnames(nbs.bed)[1:3] = c("chrom", "chromstart", "chromend")
#   # nbs.bed$chromstart = as.numeric(nbs.bed$chromstart)
#   # nbs.bed$chromend = as.numeric(nbs.bed$chromend)
#   # plotBed(beddata=nbs.bed, chrom = g.chr, chromstart=chr.start, chromend=chr.end, row  = "auto", wiggle=0.001)
#   # labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
#   #  
#   #gene
#   NB.genes.bed = getNBGenesBED(gESGid2gName[g], g.chr, chr.start, chr.end)
#   
#   
#   if(nrow(NB.genes.bed$others)==0)
#   {
#     plot.new()
#   }else
#   {
#     pg = plotGenes(NB.genes.bed$others[,c("chrom", "start", "stop", "gene", "score", "strand")],
#                chrom = g.chr, chromstart = chr.start,chromend=chr.end,
#                #colorby = c(0,1)[(NB.transcripts.bed$gene==g) +1],
#                #colorbycol=SushiColors(2),
#                col="blue",
#                types="exon",#NB.transcripts.bed$type, 
#                maxrows=8, bheight=0.2,
#                plotgenetype="box",bentline=FALSE, wigglefactor = 0.1,
#                labeloffset=.4, fontsize=0.6,# arrowlength = 0.025,
#                labeltext=TRUE)
#   }
#   
#   pg = plotGenes(NB.genes.bed$g[,c("chrom", "start", "stop", "gene", "score", "strand")],
#                  chrom = g.chr, chromstart = chr.start,chromend=chr.end,
#                  #colorby = c(0,1)[(NB.transcripts.bed$gene==g) +1],
#                  #colorbycol=SushiColors(2),
#                  col="red",
#                  types="exon",#NB.transcripts.bed$type, 
#                  maxrows=1, bheight=0.2,
#                  plotgenetype="box",bentline=FALSE, wigglefactor = 0.1,
#                  labeloffset=.4, fontsize=1.2,# arrowlength = 0.025,
#                  labeltext=TRUE)
#   labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
#   
#   
#   #axis(side=2,las=2,tcl=.2)
#   #mtext("",side=2,line=1.75,cex=.75,font=2)
#   
#   #hiC
#   pbpe = plotBedpe(nb.links.bedpe, chrom = g.chr, 
#                    chromstart = chr.start,
#                    chromend  = chr.end,
#                    heights = nb.links.bedpe$score,
#                    plottype="loops",
#                    main=g)
#   labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
#   
#   #annotations bed
#   annot.bed.df= data.frame()
#   pks.ovlp = intersect(pks, rownames(HM.PK.hg19.annot[[HM.LINK]]))
#   pks.annot = rbind(HM.PK.hg19.annot[[HM.LINK]][pks.ovlp,])
#   for(annot.nm in colnames(HM.PK.hg19.annot[[HM.LINK]]))
#   {
#     pks.filt = pks.ovlp[pks.annot[pks.ovlp, annot.nm]!=0]
#     if(length(pks.filt)==0) next
#   
#     buf= cbind(sapply(pks.filt, FUN=function(x) unlist(strsplit(x, split=":|-", perl=T))))
#     buf.chr = buf[1,]
#     buf.start = as.numeric(buf[2,])
#     buf.end = as.numeric(buf[3,])
#     
#     df =data.frame(chrom=buf.chr,
#                start = buf.start,
#                end = buf.end,
#                name=annot.nm,
#                score=0,
#                strand = ".", stringsAsFactors = F)
#     annot.bed.df = rbind(annot.bed.df, df)
#   }
#   annot.bed.df$row = as.numeric(factor(annot.bed.df$name))
#   annot.bed.df$color = maptocolors(annot.bed.df$row, col=SushiColors(6))
#   annot.bed.df.sorted = annot.bed.df[order(annot.bed.df$row),]
#   plotBed(beddata = annot.bed.df.sorted,
#           chrom = g.chr,
#           chromstart = chr.start,
#           chromend =chr.end,
#           rownumber = annot.bed.df.sorted$row, 
#           type = "region",
#           color=annot.bed.df.sorted$color,
#           row="given",
#           plotbg="white",
#           rowlabels=unique(annot.bed.df.sorted$name),
#           rowlabelcol=unique(annot.bed.df.sorted$color), 
#           rowlabelcex=0.75)
#   labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
#   mtext(paste0(HM.LINK, " bulk annotate"), side=3, adj=-0.065,line=0.5,font=2)
#   
#   #GWAS
#   
#   plotBedgraph(gwas.bedgraph,
#                chrom = g.chr, 
#                chromstart = chr.start,
#                chromend  = chr.end,
#                colorbycol= SushiColors(5))
#   labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
#   mtext("BP GWAS", side=3, adj=-0.065,line=0.5,font=2)
#   axis(side=2,las=2, tcl=.2)
#   
#   
#   
#   #
# 
#   
# }
# dev.off()




#


  