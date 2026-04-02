#V1.1.1
#a.2.linkGenePeak_MR-Egger_V1.1.1_gARE_perGene.R

#V1.1
#include enhancer 1MB around gene
#QTL around 100kb around enhancers
#filter hQTL based on multiple test for each enhancer
#further filter based on R2
#only check enhancers 10kb near eQTL hits

##use MR-Egger regression to prediction link between enhancer and gene 
#work on all genes 
#data loading  part is the same from PRS code
  
if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



args = commandArgs(trailingOnly=TRUE)
TISS = args[1] #Brain Heart Muscle Lung
HM = args[2]
FILE.IN.G2PKs.RData = args[3] #Brain Heart Muscle Lung
FILE.TMP.LINK.RData = args[4]
hQTL.WIND = as.numeric(args[5])
JOB.I = as.numeric(args[6])
JOBs.N = as.numeric(args[7])
JOB.I.inFN = paste(c(rep(0, nchar(JOBs.N)-nchar(JOB.I)), JOB.I), collapse="")

options(scipen = 999)

# library(dplyr)
# library(tidyr)
# library(zqtl)
source('~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R')
library(data.table)
# library(pheatmap)
# library("RColorBrewer")
# library(grid)
# library(gridExtra)
# library(ggplot2)
library(MendelianRandomization)
#library(coloc)
#library(ggforce)
#library(segclust2d)
#library("RSpectra")



haQTL.QVAL.CUTOFF =0.2
SNP.PRUNE.R2.CUTOFF = 0.6
SNP.MIN.NUM = 3
Gene2QTL.WIND.SIZE= 1000000
# Enh2QTL.WIND.SIZE = 100000
#Gene2QTL.WIND.SIZE = 2000000

tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
                  Lung = "Lung",
                  Muscle = "Muscle_Skeletal",
                  Heart = "Heart_Left_Ventricle",
                  Blood = "Whole_Blood")


#file.in.gene.pos = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.egenes.txt.gz")#"/broad/compbio/data/GTEx/v8/references/gencode.v26.GRCh38.genes.gtf"
dir.in.genotype.plink = "/broad/compbio/data/GTEx/GTEx_restricted/v8_plink/plink-relaxed/"
# file.in.totalVCF = "/broad/compbio/data/GTEx/GTEx_restricted/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"
#file.in.eGTEx.mergedPeaks.RData = "../peakMergeNormal/c_mergePeaksFromTissues_hg19/4Tiss.mergedPeak2tissuePeak.hg19.RData"
#file.in.histon.RData = paste0("/broad/compbio/data/eGTEx/H3K27ac/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2/", TISS, ".heightsMean.depthGCNormed.RData")
# file.in.hg382hg19.RDS =paste0("/broad/compbio/lhou/codes/eGTEX-H3K27ac/peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", TISS, ".pks.hg382hg19Uniq.RDS") 
#"/broad/compbio/ypp/eGTEx/result/step2/Brain/777.bed.gz"
file.in.eQTL = paste0("/broad/compbio/data/gtex_eqtl_tabix/hg38/", tiss2GTExTiss[TISS], ".bed.gz")
  #paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/", tiss2GTExTiss[TISS], ".allpairs.txt.gz") #"/broad/compbio/ypp/eGTEx/gtex_eqtl_v8/Brain_Frontal_Cortex_BA9.bed.gz"

dir.in.hQTL = paste0("../haQTL/a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, "_nominal/")

dir.tmp =paste0("~/hptmp/MayoBipolar/GeneEnhLink_MR-egger_V1.1/", TISS, "_", HM, "_", hQTL.WIND/1000, "k/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
# file.tmp.g.bed = paste0(dir.tmp, GENE,".bed")
# file.tmp.pks.hg38.bed = paste0(dir.tmp, CHR, ".peaks.hg38.bed")
# file.tmp.chr.g.SNP = paste0(dir.tmp, TISS, ".SNPs.txt")
# file.tmp.chr.g.SNP.r2 = paste0(dir.tmp, TISS, ".g.SNPsR2")

# dir.ou=paste0("a_geneEnhLink_MR-Egger_V1.1/", TISS, "/")
# dir.create(dir.ou, showWarnings = F, recursive = T)
#file.tmp.link.RData  = paste0(dir.tmp, TISS, ".geneEnhLink_", JOB.I.inFN,".RData")



# 
load(FILE.IN.G2PKs.RData)


JOB.geneN = length(eGene2gAREs.hg38)/JOBs.N
gene.startID = round((JOB.I-1)*JOB.geneN)+1
gene.endID = round(JOB.I*JOB.geneN)
print(paste(gene.startID, gene.endID, sep=" "))

print(paste(gene.startID, "-", gene.endID, " included in this job"))
#genes.slct = names(eGene2gAREs.hg38)[gene.startID:gene.endID]
set.seed(666)
genes.slct = sample(names(eGene2gAREs.hg38), length(eGene2gAREs.hg38), replace = F)[gene.startID:gene.endID] #to random gene list and avoid putting too many MHC genes in on ej ob



prunedSortedSNPWithR2 = function(snps.sorted,  #sorted with increasing p-value
                                 corr.r2,
                                 r2.cutoff)
{
  snps.pruned = snps.sorted[1]
  
  for(snp.test in snps.sorted[-1])
  {
    if(all(!is.na(corr.r2[snp.test, snps.pruned])) && all(corr.r2[snp.test, snps.pruned]<= r2.cutoff))
    {
      snps.pruned = c(snps.pruned, snp.test)
    }
  }
  
  return(snps.pruned)
}



gene2enh.MR_Egger=list()
for(gene in genes.slct)
{
  print(gene)
  g.chr = eGene2TSS[gene, "gene.chr"]
  g.pks=eGene2gAREs.hg38[[gene]]
  
  file.in.haQTL = system(paste0("find ", dir.in.hQTL, g.chr, ".*sorted.bed.gz"), intern = T) 
  ################################################################
  ## Read eQTL statistics
  eQTL.lb = max(1, eGene2TSS[gene, "gene.tss"]-Gene2QTL.WIND.SIZE)
  eQTL.ub = eGene2TSS[gene, "gene.tss"]+Gene2QTL.WIND.SIZE
  
  query = g.chr %&&% ':' %&&% eQTL.lb %&&% '-' %&&% eQTL.ub
  
  eqtl.tab = seqminer::tabix.read.table(file.in.eQTL, query) 
  eqtl.tab.filt = eqtl.tab[eqtl.tab$gene==gene,]
  rownames(eqtl.tab.filt) <- eqtl.tab.filt$rs <- gsub("_b38", "", eqtl.tab.filt$rs)
  
  # hqtl.tab = seqminer::tabix.read.table(file.in.haQTL, query)
  # hqtl.tab.filt = hqtl.tab[(hqtl.tab$peak %in% g.pks),] # &  !is.na(hqtl.tab$rs)
  

  ################################################################
  ##r2 extracted
  file.tmp.chr.g.SNP = paste0(dir.tmp, TISS, ".", gene, ".SNPs.txt")
  file.tmp.chr.g.SNP.r2 = paste0(dir.tmp, TISS, ".", gene, ".SNPsR2")
  
  buf = read.table(paste0(dir.in.genotype.plink, g.chr, ".bim"), sep="\t", header=F, row.names=NULL, stringsAsFactors = F)
  snpIDs.plink = buf[,2]
  
  #snpID2rs.ovlp = levels(factor(intersect(eqtl.tab.filt$rs, hqtl.tab.filt$rs)))
  snpID2rs.ovlp = levels(factor(eqtl.tab.filt$rs))
  #names(snpID2rs.ovlp) = gsub("_b38", "", snpID2rs.ovlp)
  names(snpID2rs.ovlp) = gsub("_", ":", snpID2rs.ovlp)
  
  snpIDs.ovlp.filt = intersect(snpIDs.plink, names(snpID2rs.ovlp)) #snps are in the order of those in plinks
  write(snpIDs.ovlp.filt, file.tmp.chr.g.SNP)
  
  
  plink.cmd = paste0("plink -bfile ", dir.in.genotype.plink, g.chr,
                     " --r2 square",
                     " --extract ", file.tmp.chr.g.SNP,
                     #" --ld-snp-list ", 
                     " --out ", file.tmp.chr.g.SNP.r2)
  
  if((!file.exists(paste0(file.tmp.chr.g.SNP.r2, ".ld"))) ||  
     (!system(paste0("grep \"End time\" ", file.tmp.chr.g.SNP.r2, ".log|wc -l"), intern=T)=="1"))
  {
    system(plink.cmd)
  }
  rs.ovlp.r2 = as.matrix(fread(paste0(file.tmp.chr.g.SNP.r2, ".ld"), sep="\t", header=F, data.table=F, colClasses="numeric",showProgress=T))
  #rs.ovlp.r2 = as.matrix(read.table(paste0(file.tmp.chr.g.SNP.r2, ".ld"), sep="\t", header=F, row.names = NULL))
  colnames(rs.ovlp.r2) = snpID2rs.ovlp[snpIDs.ovlp.filt]
  rownames(rs.ovlp.r2) = snpID2rs.ovlp[snpIDs.ovlp.filt]
  rs.ovlp.filt = colnames(rs.ovlp.r2)[apply(rs.ovlp.r2, 2, FUN=function(x){sum(is.na(x))==0})]
  file.remove(paste0(file.tmp.chr.g.SNP.r2, ".ld"))

  ################################################################
  ##r2 extracted
  #need to filter snps whose r2 is not available
  #prunings snps step by step on those with strong correlation with lead haQTLs


  pk2MR = lapply(g.pks,
  FUN=function(g.pk)
  {
    print(g.pk)
    buf=unlist(strsplit(g.pk, ":|-", perl=T))
    pk.lb = max(1, as.numeric(buf[2])-hQTL.WIND)
    pk.ub = as.numeric(buf[3])+hQTL.WIND
    query = g.chr %&&% ':' %&&% pk.lb %&&% '-' %&&% pk.ub
    hqtl.tab = seqminer::tabix.read.table(file.in.haQTL, query)
    hqtl.tab.filt = hqtl.tab[(hqtl.tab$peak == g.pk),] 
    rownames(hqtl.tab.filt) = hqtl.tab.filt$rs
    hqtl.tab.filt$hqtl.qval = p.adjust(hqtl.tab.filt$pval, method="BH")
    
    hqtl.tab.flip= hqtl.tab.filt
    hqtl.tab.flip$rs = paste(hqtl.tab.filt$chr, hqtl.tab.filt$end, hqtl.tab.filt$a1, hqtl.tab.filt$a2, sep="_")
    hqtl.tab.flip$beta = -hqtl.tab.filt$beta
    rownames(hqtl.tab.flip) = hqtl.tab.flip$rs
    
    hqtl.tab.merged= rbind(hqtl.tab.filt,hqtl.tab.flip)
    #filter by eQTL
    rs.gAndPk = intersect(eqtl.tab.filt$rs[eqtl.tab.filt$beta!=0], hqtl.tab.merged$rs[hqtl.tab.merged$hqtl.qval<= haQTL.QVAL.CUTOFF])
    rs.gAndPk.filter = intersect(rs.gAndPk, rs.ovlp.filt)
    
    if(length(rs.gAndPk.filter)< SNP.MIN.NUM)
    {
      print(length(rs.gAndPk.filter))
      return(NA)
    }
    #prune haQTLs snps
    rs.gAndPk.filter.sorted = rs.gAndPk.filter[order(hqtl.tab.merged[rs.gAndPk.filter, "pval"], decreasing = F)]
    rs.gAndPk.filter.sorted.pruned = prunedSortedSNPWithR2(rs.gAndPk.filter.sorted, 
                                                           rs.ovlp.r2[rs.gAndPk.filter.sorted, rs.gAndPk.filter.sorted],
                                                           SNP.PRUNE.R2.CUTOFF)
    
    if(length(rs.gAndPk.filter.sorted.pruned) < SNP.MIN.NUM)
    {
      return(NA)
    }
    
    
    b_se = cbind(bx = hqtl.tab.merged[rs.gAndPk.filter.sorted.pruned, "beta"],
                      bxse = abs(hqtl.tab.merged[rs.gAndPk.filter.sorted.pruned, "beta"]/qnorm(hqtl.tab.merged[rs.gAndPk.filter.sorted.pruned, "pval"]/2)),
                      by = eqtl.tab.filt[rs.gAndPk.filter.sorted.pruned, "beta"],
                      byse = abs(eqtl.tab.filt[rs.gAndPk.filter.sorted.pruned, "beta"]/qnorm(eqtl.tab.filt[rs.gAndPk.filter.sorted.pruned, "pval"]/2)))
    rownames(b_se) = rs.gAndPk.filter.sorted.pruned
    b_se.filt = b_se[apply(b_se, 1, FUN=function(x){sum(is.na(x))==0}), ]
    rs.gAndPk.filter.sorted.pruned = rownames(b_se.filt)
    
    
    
    
    
    MRInputObject.cor <- mr_input(bx=b_se.filt[,"bx"],
                                  bxse=b_se.filt[,"bxse"],
                                  by=b_se.filt[,"by"],
                                  byse=b_se.filt[,"byse"],
                                  snps = rs.gAndPk.filter.sorted.pruned,
                                  corr = rs.ovlp.r2[rs.gAndPk.filter.sorted.pruned, rs.gAndPk.filter.sorted.pruned]
                               )
    print(max(abs(rs.ovlp.r2[rs.gAndPk.filter.sorted.pruned, rs.gAndPk.filter.sorted.pruned][upper.tri(rs.ovlp.r2[rs.gAndPk.filter.sorted.pruned, rs.gAndPk.filter.sorted.pruned])])))
    
    result = tryCatch(
    {
      mr_egger(MRInputObject.cor,
                        correl = TRUE,
                        distribution = "normal",
                        alpha = 0.05)  
    }, error = function(e) {
        print(e)
      return(NA)
    })
    
    #print(result@Causal.pval)
    return(result)
    
  })
  names(pk2MR) = g.pks

  gene2enh.MR_Egger[[gene]]=pk2MR
  save(gene2enh.MR_Egger,
       file=paste0(FILE.TMP.LINK.RData, ".tmp"))
  
}
#names(gene2enh.MR_Egger) = genes.slct

save(gene2enh.MR_Egger,
     file=FILE.TMP.LINK.RData)

# eqtl.tab.filt, 
#      hqtl.tab.filt,
#      rs.ovlp.r2,
#      g.pks,

