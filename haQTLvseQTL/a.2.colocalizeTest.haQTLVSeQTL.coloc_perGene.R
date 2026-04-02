#V1.3
#a.2.colocalizeTest.haQTLVSeQTL.coloc_V1.3.gARE_perGene.R


#V1.2
#include enhancer 1MB around gene
#QTL around 100kb around enhancers
#filter hQTL based on multiple test for each enhancer
#further filter based on R2
#only check enhancers 10kb near eQTL hits
#support re-run

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
FILE.IN.G2PKs.RData = args[3]
FILE.TMP.COLOC.RData = args[4]
hQTL.WIND = as.numeric(args[5])
JOB.I = as.numeric(args[6])
JOBs.N = as.numeric(args[7])
JOB.I.inFN = paste(c(rep(0, nchar(JOBs.N)-nchar(JOB.I)), JOB.I), collapse="")

#eCAVIAR = "~/lhou.compbio/software/caviar/CAVIAR-C++/eCAVIAR"

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
#library(MendelianRandomization)
library(coloc)
#library(ggforce)
#library(segclust2d)
#library("RSpectra")



#haQTL.PVAL.CUTOFF =1e-3
#SNP.PRUNE.R2.CUTOFF = 0.6
SNP.MIN.NUM = 3
Gene2Enh.WIND.SIZE = 1000000
Gene2QTL.WIND.SIZE= 1000000
#hQTL.WIND = 100000
#Gene2QTL.WIND.SIZE = 2000000


tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
                  Lung = "Lung",
                  Muscle = "Muscle_Skeletal",
                  Heart = "Heart_Left_Ventricle",
                  Blood = "Whole_Blood")
HM2haQTL.sampN = c(H3K27ac= 176,
                   H3K4me3 = 165,
                   H3K4me1 = 149,
                   H3K36me3 = 165,
                   H3K27me3 = 178)
HM2Cov.Num = c(H3K27ac = 10,
              H3K36me3 = 10,
              H3K4me1 = 10,
              H3K27me3 = 10,
              H3K4me3 = 10)

tiss2eQTL.sampN = c(Brain= 175,
                    Heart = 386,
                    Muscle= 706,
                    Lung= 515,
                    Blood=670)


#file.in.gene.pos = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.egenes.txt.gz")#"/broad/compbio/data/GTEx/v8/references/gencode.v26.GRCh38.genes.gtf"
#dir.in.genotype.plink = "/broad/compbio/data/GTEx/GTEx_restricted/v8_plink/plink-relaxed/"
# file.in.totalVCF = "/broad/compbio/data/GTEx/GTEx_restricted/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"
#file.in.eGTEx.mergedPeaks.RData = "../peakMergeNormal/c_mergePeaksFromTissues_hg19/4Tiss.mergedPeak2tissuePeak.hg19.RData"
#file.in.histon.RData = paste0("/broad/compbio/data/eGTEx/H3K27ac/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2/", TISS, ".heightsMean.depthGCNormed.RData")
# file.in.hg382hg19.RDS =paste0("/broad/compbio/lhou/codes/eGTEX-H3K27ac/peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", TISS, ".pks.hg382hg19Uniq.RDS") 
#"/broad/compbio/ypp/eGTEx/result/step2/Brain/777.bed.gz"
file.in.eQTL = paste0("/broad/compbio/data/gtex_eqtl_tabix/hg38/", tiss2GTExTiss[TISS], ".bed.gz")
dir.in.haQTL = paste0("../haQTL/a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_100k/", HM, "_nominal/")
  #paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/", tiss2GTExTiss[TISS], ".allpairs.txt.gz") #"/broad/compbio/ypp/eGTEx/gtex_eqtl_v8/Brain_Frontal_Cortex_BA9.bed.gz"

#dir.tmp =dirname(FILE.TMP.COLOC.RData)
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
set.seed(666)
genes.slct = sample(names(eGene2gAREs.hg38), length(eGene2gAREs.hg38), replace = F)[gene.startID:gene.endID] #to random gene list and avoid putting too many MHC genes in on ej ob


eQTLvshaQTL.coloc=lapply(genes.slct,
FUN=function(gene)
{
  print(paste0("##### ", gene))
  g.chr = eGene2TSS[gene, "gene.chr"]
  g.pks=eGene2gAREs.hg38[[gene]]
  
  file.in.haQTL = system(paste0("find ", dir.in.haQTL, g.chr, ".*sorted.bed.gz"), intern = T) 
  ################################################################
  ## Read eQTL statistics
  eQTL.lb = max(1, eGene2TSS[gene, "gene.tss"]-Gene2QTL.WIND.SIZE)
  eQTL.ub = eGene2TSS[gene, "gene.tss"]+Gene2QTL.WIND.SIZE
  
  query = g.chr %&&% ':' %&&% eQTL.lb %&&% '-' %&&% eQTL.ub
  
  eqtl.tab = seqminer::tabix.read.table(file.in.eQTL, query) 
  eqtl.tab.filt = eqtl.tab[eqtl.tab$gene==gene,]
  rownames(eqtl.tab.filt) <- eqtl.tab.filt$rs <- gsub("_b38", "", eqtl.tab.filt$rs)
  
  #

  pks.coloc = lapply(g.pks,
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
    #hqtl.qval = p.adjust(hqtl.tab.filt$pval, method="BH")
    hqtl.tab.flip= hqtl.tab.filt
    hqtl.tab.flip$rs = paste(hqtl.tab.filt$chr, hqtl.tab.filt$end, hqtl.tab.filt$a1, hqtl.tab.filt$a2, sep="_")
    hqtl.tab.flip$beta = -hqtl.tab.filt$beta
    rownames(hqtl.tab.flip) = hqtl.tab.flip$rs
    
    hqtl.tab.merged= rbind(hqtl.tab.filt,hqtl.tab.flip)
    
    
    # #
    # if(!any(hqtl.tab.merged$pval<= haQTL.PVAL.CUTOFF))
    # {
    #   return(NA)
    # }
    #filter by eQTL
    snp.loc.ovlp = intersect(eqtl.tab.filt$rs, hqtl.tab.merged$rs)
    
    if(length(snp.loc.ovlp)==0)
    {
      return(NA)
    }
    
    #
    #pval = 2*pt(-abs(buf$beta/buf$se),df=samp.N-2)
    
    dataset.eQTL = list(beta=eqtl.tab.filt[snp.loc.ovlp, "beta"], 
                        varbeta=eqtl.tab.filt[snp.loc.ovlp, "se"]^2,
                        N=tiss2eQTL.sampN[TISS], type="quant")
    dataset.haQTL = list(beta=hqtl.tab.merged[snp.loc.ovlp, "beta"], 
                        varbeta=(-abs(hqtl.tab.merged[snp.loc.ovlp, "beta"])/qt(hqtl.tab.merged[snp.loc.ovlp, "pval"]/2,
                                                                              df=HM2haQTL.sampN[HM]-2-HM2Cov.Num[HM]))^2,
                        N=HM2haQTL.sampN[HM], type="quant")
    
    my.res <- coloc.abf(dataset1=dataset.eQTL,
                    dataset2=dataset.haQTL,
                    MAF=eqtl.tab.filt[snp.loc.ovlp, "maf"])#,
      
    
    #snp2CLPP = as.matrix(read.table(paste0(f.tmp.eCAV.res, "_col"), sep="\t", head=T, row.names = 1))[,2]
    return(my.res$summary)
  })
  names(pks.coloc) = g.pks

  return(pks.coloc)
})
names(eQTLvshaQTL.coloc) = genes.slct

save(eQTLvshaQTL.coloc,
     file=FILE.TMP.COLOC.RData)

# eqtl.tab.filt, 
#      hqtl.tab.filt,
#      rs.ovlp.r2,
#      g.pks,

