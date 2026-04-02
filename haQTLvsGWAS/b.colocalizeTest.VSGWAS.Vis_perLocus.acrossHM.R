

args = commandArgs(trailingOnly=TRUE)
HM=args[1] #check the loci showing coloc with this histone mark

GWAS.NM ="PGC.BP"



if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))





#eCAVIAR = "~/lhou.compbio/software/caviar/CAVIAR-C++/eCAVIAR"

options(scipen = 999)


# library(dplyr)
# library(tidyr)
# library(zqtl)
source('~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R')
library(data.table)
#ibrary(pheatmap)
library("RColorBrewer")
# library(grid)
# library(gridExtra)
library(ggplot2)
library(Sushi)
library(gtable)
#library(MendelianRandomization)
#library(coloc)
#library(ggforce)
#library(segclust2d)
#library("RSpectra")
MY.PALETTE <- colorRampPalette(c("royalblue2", "green1", "orange2", "darkred")) #rev(brewer.pal(11, "Spectral"))
#SC_COLOR = scale_colour_gradientn(colours = MY.PALETTE(8), limits=c(0,1)) 


haQTL.WIND=100000
PK.WIND = haQTL.WIND #200000#as.numeric(args[2])
HAQTL.CUTOFF="emp.p.fdr.cutoff0.2"


COLOC.H4.ABF.cutoff=0.5
COLOC.H4.ABF.lower.cutoff=0.1


HMs = c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3")

TYPE2PCH=c("rest"=19,
           "leadGWAS"=17,
           "haQTL"=4,
           "leadGWAS+haQTL"=13)
TYPE2SIZE=c("rest"=0.5,
           "leadGWAS"=2,
           "haQTL"=1,
           "leadGWAS+haQTL"=2)

dir.in.Mayo.genotype.plink = "~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/"

file.in.haQTLvsGWAS.RData =  paste0("a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS_summary/bulkhQTLvsGWAS.coloc.RData")

dirs.in.hQTL = paste0("../haQTL/a_4_haQTL_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, "_nominal/")
names(dirs.in.hQTL) = HMs

files.in.sig.haQTL = paste0("../haQTL/a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", HMs, ".haQTLPeak.cutoffs.hg38.RData")
names(files.in.sig.haQTL) = HMs


file.in.mAREs.hg19.RData="../peakVariationAcrossTiss/a_mergePeaksFromHMs_hg19/5HMs.mergedPeak2HMPeak.hg19.RData"
#file.in.mARE.modules.RData="../peakVariationAcrossTiss/a_4_modules_annotations_mergedPk5HMs_Epimap_CSandeGTExActivity/mARE2ModuleGrpName.RData"
files.in.ARE.hg38toHg19.RDS=paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HMs, ".pk.hg38tohg19.RDS")
names(files.in.ARE.hg38toHg19.RDS) = HMs

files.in.GWAS.gz = c(PGC.BP = "~/lhou.compbio/data/GWAS/PGC_leaveMayoBDout/BDlooMAYO_PGC3_hg38.sorted.bed.gz", #"~/lhou.compbio/data/GWAS/pgc_bipolar/pgc_bip_hg38.sorted.bed.gz",
                     GWASAtlas_ukb2.highBloodPressure="~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/highBloodPressure_6150_4_logistic.EUR.hg38.sorted.bed.gz",
                     GWASAtlas_ukb2.diabetes="~/lhou.compbio/data/GWAS/GWASAtlas_ukb2_NG_2019/diabetes_f.2443.0.0_logistic.EUR.hg38.sorted.bed.gz",
                     DrinksPerWeek = "~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/DrinksPerWeek.txt.hg38.sorted.bed.gz",
                     SmokingInitiation = "~/lhou.compbio/data/GWAS/Alcohol_tobacco_NG_2019/SmokingInitiation.txt.hg38.sorted.bed.gz")

#file.in.eQTL2pk.RData=paste0("../testPromEnhLink_eQTL/a_AREsNeareQTL/", TISS, "_10k/", TISS, ".gene2peaks_byeQTL.hg38.RData")
file.in.geneStructure.RData = "~/lhou.compbio/data/gene/human/gencode.v26.GRCh38.genes_GTEx.v8.forSushi.RData"

dir.tmp = paste0("~/hptmp/MayoBipolar/b_haQTL.GWAS.coloc_vis_byLocus/", GWAS.NM, "_", HM, "_", PK.WIND/1000, "kb/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
#file.tmp.FMeQTL.bed = paste0(dir.tmp, COND, ".fmeQTL.bed")
#file.tmp.eQTL.ARE.bed = paste0(dir.tmp, COND, ".eQTL-ARE.bed")
#file.tmp.ARE.OVLP.FMeQTL = paste0(dir.tmp, COND, ".FMeQTL-ARE.bed")
file.tmp.genesStructure.bed = paste(dir.tmp, GWAS.NM, ".genesStructure.bed", sep="")
file.tmp.gRegion.bed= paste(dir.tmp, GWAS.NM, ".", "gRegion.bed", sep="")
file.tmp.gStructure.bed = paste(dir.tmp, GWAS.NM, ".", "gStructure.bed", sep="")



dir.ou = paste0("b.coloc.haQTLvsGWAS_Vis_byLocus_acrossHMs/", GWAS.NM, "_", HM, "_", PK.WIND/1000, "kb/")
dir.create(dir.ou, showWarnings = F, recursive = T)
#file.ou.locID2dbSNPid.RDS = paste0("b.coloc.haQTLvsGWAS_Vis_byLocus_acrossHMs/locID2dbSNPid.RDS")
#file.ou.overall.pdf = paste0(dir.ou, GWAS.NM, ".mARE.haQTLvsGWAS.pdf")
file.ou.examp.pdf = paste0(dir.ou, HM, ".haQTLvsGWAS.examples.pdf")
file.ou.forVis.RData = paste0(dir.ou, HM, ".haQTLvseQTLvsGWAS.df.RData")
file.ou.forVis.RDS = paste0(dir.ou, HM, ".haQTLvseQTLvsGWAS.df.RDS")
# 

#

# #data load
# if(file.exists(file.ou.locID2dbSNPid.RDS))
# {
#   locID2dbSNPid=readRDS(file.ou.locID2dbSNPid.RDS)
# }else
# {
#   buf=fread(file.in.locID2dbSNPid, sep="\t", header = T, data.table=F, 
#             colClasses=c("character", "numeric", "character", "character", "character", "numeric", "character"),
#             showProgress=T)
#   locID2dbSNPid = buf[,"rs_id_dbSNP150_GRCh38p7"]
#   names(locID2dbSNPid) = buf[,"variant_id"]
#   locID2dbSNPid[locID2dbSNPid=="."] = names(locID2dbSNPid)[locID2dbSNPid=="."]
#   saveRDS(locID2dbSNPid,
#           file=file.ou.locID2dbSNPid.RDS)
# }


# eGene.id2name=c()
# for(f.i in files.in.eGene.info)
# {
#   gs.info=read.table(f.i, sep="\t", header = T, row.names = 1)
#   #gid2gName = gs.info$gene_name
#   #names(gid2gName)=rownames(gs.info)
#   eGene.id2name[rownames(gs.info)] = gs.info$gene_name
# }


hm.pk.hg38tohg19 = lapply(files.in.ARE.hg38toHg19.RDS, readRDS)
hm.pk.hg19tohg38 = lapply(hm.pk.hg38tohg19,
FUN=function(pk.hg38tohg19)
{
  pk.hg19tohg38 = names(pk.hg38tohg19)
  names(pk.hg19tohg38) = pk.hg38tohg19
  
  return(pk.hg19tohg38)
})

load(file.in.mAREs.hg19.RData)
# load(file.in.mARE.modules.RData)
# names(mARE2ModGrpName)=gsub(":|-", "\t", names(mARE2ModGrpName))
#
load(file.in.haQTLvsGWAS.RData)
mARE.hm2coloc.H4 = mARE.hm_gwas.coloc.H4.filt.BP[, grepl(GWAS.NM, colnames(mARE.hm_gwas.coloc.H4.filt.BP))]
colnames(mARE.hm2coloc.H4)=sub(paste0(";", GWAS.NM), "", colnames(mARE.hm2coloc.H4))
mARE.hm2coloc.H4.filt=mARE.hm2coloc.H4[!is.na(mARE.hm2coloc.H4[, HM]) & mARE.hm2coloc.H4[, HM]>=COLOC.H4.ABF.cutoff,]

hm.sig.haQTLs = lapply(files.in.sig.haQTL,
FUN=function(f.i)
{
  load(f.i)
  hQTLpeak2signifQTL.bycutoffs[[HAQTL.CUTOFF]]
})


save(
  hm.pk.hg19tohg38,
  hm.mergedPeak2peaks,
  mergedPeak2hmPeaks,
  #mARE2ModGrpName,
  mARE.hm2coloc.H4.filt,
  hm2gwas2coloc.hg38.list, 
  hm.sig.haQTLs,
  file=file.ou.forVis.RData)

load(file.in.geneStructure.RData)
write.table(transcripts.bed, file=file.tmp.genesStructure.bed, sep="\t", col.names = F, row.names = F, quote=F)


#for sushi function
getNBGenesBED.targetGenes = function(chr, chr.start, chr.end, genes=c())
{
  
  write(paste(chr, "\t", chr.start, "\t", chr.end, sep=""), file.tmp.gRegion.bed)
  cmd =paste( "intersectBed -a ", file.tmp.genesStructure.bed, " -b ", file.tmp.gRegion.bed," -u>", file.tmp.gStructure.bed, sep="")
  system(cmd)
  
  buf = tryCatch(
  {
    read.table(file.tmp.gStructure.bed, head=F, sep="\t", row.names = NULL, stringsAsFactors = F)
  },
  error=function(cond)
  {
    print(cond)
    return(c())
  }
  )
  if(length(buf)==0)
  {
    NB.genes.bed = c()
  }else{
    colnames(buf) = c("chrom", "start", "stop", "gene.ens", "gene", "gene.type", "transcript", "score", "strand", "type")
    buf = buf[buf$gene.type=="protein_coding",]#only keep protein coding genes and the gene we are focusing on now
    buf$strand[as.character(buf$strand)=="+"] = 1
    buf$strand[as.character(buf$strand)=="-"] = -1
    buf$strand = as.numeric(buf$strand)
    
    NB.genes.bed = buf
    NB.genes.bed$is.target = buf$gene.ens %in% genes
    # 
  }
  return(NB.genes.bed)
}

plotNBGenesBed = function(NB.genes.bed, chr, chr.start, chr.end)
{
  if(length(NB.genes.bed)==0 || nrow(NB.genes.bed)==0)
  {
    plot.new()
  }else
  {
    pg = plotGenes(NB.genes.bed[,c("chrom", "start", "stop", "gene", "score", "strand")],
                   chrom = chr, 
                   chromstart = chr.start,
                   chromend=chr.end,
                   #colorby = c(0,1)[(NB.transcripts.bed$gene==g) +1],
                   #colorbycol=SushiColors(2),
                   col=c("blue", "red")[NB.genes.bed$is.target +1],
                   types="exon",#NB.transcripts.bed$type, 
                   maxrows=8, bheight=0.2,
                   plotgenetype="box",bentline=FALSE, wigglefactor = 0.1,
                   labeloffset=.4, fontsize=1.2,# arrowlength = 0.025,
                   labeltext=TRUE)
  }
  mtext("genes",side=2,line=4,cex=1,font=1.5)
  
}




#
haQTL.GWAS.df.list=list()
#mARE2leadGWASishaQTL.topTiss=c()
#iseQTLColoc.topTiss =c()
# gAREs.topColoc=c()

pdf(file.ou.examp.pdf, width=5, height=10)
#pdf("test.pdf", width=15, height=6)
par(mar=c(4,8,2,10), xpd=TRUE)


for(i in 1:nrow(mARE.hm2coloc.H4.filt))
{
  print(i)
  mARE.hg19 = rownames(mARE.hm2coloc.H4.filt)[i]
  
  #mARE.modGrp = mARE2ModGrpName[mARE.hg19]
  
  AREs.hg19 = mergedPeak2hmPeaks[[mARE.hg19]]
  hm2AREs.hg38 =list()
  for(ARE in AREs.hg19)
  {
    buf = unlist(strsplit(ARE, split="_", fixed=T))
    hm=buf[1]
    
    hm2AREs.hg38[[hm]] = c(hm2AREs.hg38[[hm]], hm.pk.hg19tohg38[[hm]][buf[2]])
  }
  hm2AREs.hg38.maxH4 = lapply(names(hm2AREs.hg38),
  FUN=function(hm)
  {
    x=hm2gwas2coloc.hg38.list[[hm]][[GWAS.NM]][hm2AREs.hg38[[hm]]]
    x=x[!is.na(x)]
    x[which.max(x)[1]]
  })
  names(hm2AREs.hg38.maxH4) = names(hm2AREs.hg38)
  hm2AREs.hg38.maxH4=hm2AREs.hg38.maxH4[!(sapply(hm2AREs.hg38.maxH4, is.na))]

  gARE2coloc= unlist(hm2AREs.hg38.maxH4)
  gARE.topColoc = names(gARE2coloc)[which.max(gARE2coloc)[1]]
  buf = unlist(strsplit(gARE.topColoc, split="\\.|:|-"))
  hm.topColoc=buf[1]
  reg.chr = buf[2]
  peak.start = as.numeric(buf[3])
  peak.end= as.numeric(buf[4])
  gARE.topColoc=paste0(reg.chr, ":", peak.start, "-", peak.end)
  gARE.topColoc.sigHaQTLs=hm.sig.haQTLs[[hm.topColoc]][[gARE.topColoc]]$V2
  reg.start = max(0, peak.start-PK.WIND)
  reg.end = peak.end+PK.WIND
  query = reg.chr %&&% ':' %&&% reg.start %&&% '-' %&&% reg.end
  
  
  
  #
  #abline(v=haQTL.top.pos, xpd=F)
  
  
  #gwas
  gwas.tab = seqminer::tabix.read.table(files.in.GWAS.gz[GWAS.NM], query) 
  #rownames(gwas.tab) =  paste(gwas.tab[,1], gwas.tab[,3], gwas.tab[,5], gwas.tab[,4], sep=":")
  # gwas.p = gwas.tab$pval
  # names(gwas.p) = paste(gwas.tab[,1], gwas.tab[,3], gwas.tab[,5], gwas.tab[,4], sep="_")
  gwas.p.df = data.frame(pval=c(gwas.tab$P, gwas.tab$P),
                         pos=c(gwas.tab$stop, gwas.tab$stop),
                         rsID=c(gwas.tab$SNP, gwas.tab$SNP))
  
  rownames(gwas.p.df) = c(paste(gwas.tab[,1], gwas.tab[,3], gwas.tab[,5], gwas.tab[,6], sep="_"),
                           paste(gwas.tab[,1], gwas.tab[,3], gwas.tab[,6], gwas.tab[,5], sep="_"))

  #
  hm.haQTL.list=list()
  for (hm in names(hm2AREs.hg38.maxH4))
  {
    print(hm)
    
    #haQTL
    f.i.haQTL=system(paste0("find ", dirs.in.hQTL[hm], reg.chr, ".*sorted.bed.gz"), intern = T)  
    hqtl.tab = seqminer::tabix.read.table(f.i.haQTL, query) 
  
    gARE = names(hm2AREs.hg38.maxH4[[hm]])
    
    #sig.haQTLs= hm.sig.haQTLs[[hm]][[gARE]][,2]
    
    hqtl.tab.slct = hqtl.tab[(hqtl.tab$peak == gARE),] 
    rownames(hqtl.tab.slct) = hqtl.tab.slct$rs
      
    hm.haQTL.list[[paste0(hm, ".haQTL")]] = hqtl.tab.slct
   
  }
  
  QTL.snps.all = levels(factor(unlist(lapply(hm.haQTL.list, rownames))))
  
  #LD
  snps.all = union(rownames(gwas.p.df), QTL.snps.all)
  #names(snps.all) = gsub("_", ":", snps.all)
  
  buf = read.table(paste0(dir.in.Mayo.genotype.plink, reg.chr, ".bim"), sep="\t", header=F, row.names=NULL, stringsAsFactors = F)
  snpIDs.plink = buf[,2]
  
  file.tmp.mARE.SNP = paste0(dir.tmp, GWAS.NM, ".", gsub("\t", "_", mARE.hg19), ".SNPs.txt")
  file.tmp.mARE.SNP.r2 = paste0(dir.tmp, GWAS.NM, ".", gsub("\t", "_", mARE.hg19), ".SNPsR2")
  snpIDs.ovlp.filt = intersect(snpIDs.plink, snps.all) #snps are in the order of those in plinks
  write(snpIDs.ovlp.filt, file.tmp.mARE.SNP)
  
  plink.cmd = paste0("plink -bfile ", dir.in.Mayo.genotype.plink, reg.chr,
                     " --r2 square",
                     " --extract ", file.tmp.mARE.SNP,
                     #" --ld-snp-list ",
                     " --out ", file.tmp.mARE.SNP.r2)
  # if((!file.exists(paste0(file.tmp.mARE.SNP.r2, ".ld"))) ||
  #    system(paste0("grep \"End time\" ", file.tmp.mARE.SNP.r2, ".log|wc -l"), intern=T)!="1")
  # {
  system(plink.cmd)
  #}
  rs.ovlp.r2 = as.matrix(fread(paste0(file.tmp.mARE.SNP.r2, ".ld"), sep="\t", header=F, data.table=F, colClasses="numeric",showProgress=T))
  colnames(rs.ovlp.r2) = snpIDs.ovlp.filt#[snpIDs.ovlp.filt]
  rownames(rs.ovlp.r2) = snpIDs.ovlp.filt#[snpIDs.ovlp.filt]
  
  rs.ovlp.filt = colnames(rs.ovlp.r2)[apply(rs.ovlp.r2, 2, FUN=function(x){sum(is.na(x))==0})]
  
  gwas.rs.ovlp.filt= intersect(rownames(gwas.p.df),rs.ovlp.filt)
  GWAS.lead = gwas.rs.ovlp.filt[which.min(gwas.p.df[gwas.rs.ovlp.filt, "pval"])[1]]
  
  
  df.filt=list()
  
  df.filt$gwas = data.frame(position=gwas.p.df[gwas.rs.ovlp.filt, "pos"],
                            log10p = -log10(gwas.p.df[gwas.rs.ovlp.filt, "pval"]),
                            R2.2leadGWAS = rs.ovlp.r2[gwas.rs.ovlp.filt, GWAS.lead], 
                            stringsAsFactors = F
                            #is.lead = (gwas.rs.ovlp.filt==GWAS.lead)
                            )
  df.filt$gwas$type="rest"
  df.filt$gwas$type[gwas.rs.ovlp.filt %in% gARE.topColoc.sigHaQTLs] = "haQTL"
  df.filt$gwas$type[gwas.rs.ovlp.filt == GWAS.lead] = "leadGWAS"
  df.filt$gwas$type[gwas.rs.ovlp.filt == GWAS.lead & gwas.rs.ovlp.filt %in% gARE.topColoc.sigHaQTLs] = "leadGWAS+haQTL"
  
  #mARE2leadGWASishaQTL.topTiss[i]= any(df.filt$gwas$type=="leadGWAS+haQTL")
  
  for(nm in names(hm.haQTL.list))
  {
    if(nrow(hm.haQTL.list[[nm]])==0) next
    QTL.rs.ovlp.filt = intersect(rownames(hm.haQTL.list[[nm]]), rs.ovlp.filt)
    df.filt[[nm]] = data.frame(position = hm.haQTL.list[[nm]][QTL.rs.ovlp.filt, "end"],
                                 log10p = -log10(hm.haQTL.list[[nm]][QTL.rs.ovlp.filt, "pval"]),
                                 R2.2leadGWAS = rs.ovlp.r2[QTL.rs.ovlp.filt, GWAS.lead],
                                  stringsAsFactors = F)#,
                                 #is.lead = (QTL.rs.ovlp.filt==GWAS.lead))
    df.filt[[nm]]$type="rest"
    df.filt[[nm]]$type[QTL.rs.ovlp.filt %in% gARE.topColoc.sigHaQTLs] = "haQTL"
    df.filt[[nm]]$type[QTL.rs.ovlp.filt == GWAS.lead] = "leadGWAS"
    df.filt[[nm]]$type[QTL.rs.ovlp.filt == GWAS.lead & QTL.rs.ovlp.filt %in% gARE.topColoc.sigHaQTLs] = "leadGWAS+haQTL"
  }
  
                      
  #visualization
  #Sushi
  layout(matrix(c(1:7), ncol = 1, byrow = TRUE))
  
  #gene Annotation
  NB.genes.bed = getNBGenesBED.targetGenes(reg.chr, reg.start, reg.end, genes=NULL)
  plotNBGenesBed(NB.genes.bed, reg.chr, reg.start, reg.end)
  
  for(nm in names(df.filt))
  {
    plot(c(reg.start,reg.end), 
         c(0, max(df.filt[[nm]]$log10p)+1),
         xlim=c(reg.start, reg.end),
         #ylim=c(0,1),
         type ='n',bty='n',xaxt='n',ylab="",xlab="",xaxs="i"
       )
    if(nm == "gwas")
    {
      title(paste0("lead SNP: ", gwas.p.df[GWAS.lead, "rsID"]))
    }else
    {
      hm = gsub(".haQTL", "", nm)
      pk=hm.haQTL.list[[nm]]$peak[1]

      coloc= signif(hm2AREs.hg38.maxH4[[hm]],2)
      #coloc =signif(max(mARE.hm2coloc.H4.filt[[i]][[hm]]),2)
      title(paste0(pk, "\ncoloc PP=", coloc), adj = 0)#mARE.modGrp, 
             #" H4vsH3 Ratio=", signif(coloc.res[[nm]]["PP.H4.abf"]/coloc.res[[nm]]["PP.H3.abf"],2)), 
    }
    
    R2.cut = cut(c(1, df.filt[[nm]]$R2.2leadGWAS),breaks = 10)
    R2.breaks = sapply(levels(R2.cut),
    FUN=function(x)
    {
      #round(as.numeric(gsub("\\((\\d+\\.\\d+),.*?]", "\\1", x)),1)
      x=unlist(strsplit(x, split=","))[2]
      round(as.numeric(gsub("]", "", x)),1)
    })
    R2.breaks[seq(2, 10, 2)] =""
    points.col <- MY.PALETTE(10)[as.numeric(R2.cut)[-1]]
    points(x=df.filt[[nm]]$position, 
           y=df.filt[[nm]]$log10p, 
           pch=TYPE2PCH[df.filt[[nm]]$type], #c(19,17)[df.filt[[nm]]$is.lead+1], 
           col=points.col,
           cex=TYPE2SIZE[df.filt[[nm]]$type]) #c(0.5,2)[df.filt[[nm]]$is.lead +1])
    
    #abline(v=gene.pos, lty="dotted", col="black", lwd=3)
    if(nm=="gwas")
    {
      mtext("GWAS\n-log10(p value)", side=2, line=4, cex=1, font=1)
      
      
      legend("topright", 
             inset=c(-0.5,0), 
             legend=c(names(TYPE2PCH), "gARE"), 
             pch=c(TYPE2PCH, 15), 
             col = c(rep("black", length(TYPE2PCH)), "brown"),
             
             #title="R2 with lead eQTL",
             #cex=1,
             #pt.cex=1,
             #y.intersp=0.5#
             bty="n"
      )
      
    }else
    {
      pk=hm.haQTL.list[[nm]]$peak[1]
      peak.pos = as.numeric(unlist(strsplit(pk, split=":|-"))[2:3])
      rect(xleft=peak.pos[1], ybottom=0, xright=peak.pos[2], ytop=1/5*max(df.filt[[nm]]$log10p), col="brown")
           
      mtext(paste0(nm, "\n-log10(p value)"),side=2,line=4,cex=1,font=1)
      
      legend("topright",
             inset=c(-0.5,0),
             fill=rev(MY.PALETTE(10)),
             legend=rev(R2.breaks),
             title="R2 with lead GWAS SNP",
             cex=1,
             pt.cex=1,
             y.intersp=0.5,#
             bty="n"
      )
    }
    
    labelgenome(chrom=reg.chr, chromstart= reg.start, chromend = reg.end, n=5, scale="Mb")
    
    
  }
  
  
  haQTL.GWAS.df.list[[i]]=df.filt
}
dev.off()

# names(mARE2leadGWASishaQTL.topTiss)=names(mARE.hm2coloc.H4.filt)
# mARE2leadGWASNothaQTL.topTiss.prop = sum(!mARE2leadGWASishaQTL.topTiss)/length(mARE2leadGWASishaQTL.topTiss)



saveRDS(haQTL.GWAS.df.list, file.ou.forVis.RDS)


#
i=88

for(nm in names(haQTL.GWAS.df.list[[i]]))
{
  f.o.tsv = paste0(dir.ou, HM, ".haQTLvseQTLvsGWAS.df.", i, ".", nm, ".txt")

  write.table(haQTL.GWAS.df.list[[i]][[nm]], file=f.o.tsv, sep="\t", row.names=F, quote=F)

}




# #label gARE type
# buf=read.table(file.in.mARE2gARE.hg38, sep="\t", header = T, row.names=NULL, stringsAsFactors = F)
# hm.ARE2mARE.hg38=lapply(split(buf, buf$tis),
# FUN=function(x)
# {
#   ARE=sub("_",":", x$ARE)
#   ARE=sub("_","-", ARE)
#   ARE2mARE=x$mARE
#   names(ARE2mARE)=ARE
  
#   return(ARE2mARE)
# })

# hm.m_gARE.type=lapply(files.in.hm.gAREType,
# FUN=function(f.i)
# {
#   buf=read.table(f.i, header = T, sep="\t", row.names = 1, stringsAsFactors = F)
#   #rownames(buf)=buf$
#   #buf$gARE.hg38=paste0(buf$chr, ":", buf$start, "-", buf$end)
#   return(buf)
  
# })
# #m_gARE.type.all=do.call(rbind, hm.m_gARE.type)

# #gsub(".", "_", rownames(gARE.type.all), fixed = T)
# mAREs.bed=rownames(mARE2hm.gARE.coloc.matrix.filt)

# mAREs.bed2gARE.type=sapply(mAREs.bed,
# FUN=function(mARE)
# {
#   hm.gARE.coloc=unlist(mARE2hm.gARE.coloc.filt[[mARE]])
#   hm.gARE.topColoc = names(hm.gARE.coloc)[which.max(hm.gARE.coloc)[1]]
  
#   buf=unlist(strsplit(hm.gARE.topColoc, split=".", fixed = T))
#   hm.topColoc=buf[1]
#   gARE.topColoc=buf[2]
#   mARE.hg38 = hm.ARE2mARE.hg38[[hm.topColoc]][gARE.topColoc]
  
#   hm.m_gARE.type[[hm.topColoc]][mARE.hg38, "gARE_Type"]
# })


#
# pdf(file.ou.overall.pdf, width =7, height=15)

# # mAREs=sub("\t", ":", mAREs.bed)
# # mAREs=sub("\t", "-", mAREs)loa()
# annot.df=data.frame(missedByeQTL= factor(!iseQTLColoc.topTiss),
#                     leadGWASnothaQTL = factor(!mARE2leadGWASishaQTL.topTiss),
#                     hmSpecif = factor(mARE2hm.gARE.coloc.matrix.filt.H4.2nd<COLOC.H4.ABF.lower.cutoff),
#                     mAREGrp=factor(mARE2ModGrpName[mAREs.bed])
#                     )
# pheatmap(mARE2hm.gARE.coloc.matrix.filt,
#          cluster_rows=F,
#          cluster_cols=F,
#          annotation_colors=list(leadGWASnothaQTL=c("TRUE"="black", "FALSE"="white"),
#                                 hmSpecif=c("TRUE"="black", "FALSE"="white"),
#                                 missedByeQTL=c("TRUE"="black", "FALSE"="white"),
#                                 mAREGrp=ARE.ModGrp2COLOR
                                
#                                 ),
#          annotation_row=annot.df,
#          main=paste0(GWAS.NM, 
#                      "\nhm specific coloc prop.", mARE2hm.gARE.coloc.matrix.filt.hmSpecifProp,
#                     "\nleadGWAS not haQTL prop.", mARE2leadGWASNothaQTL.topTiss.prop,
#                     "\nmissed by eQTL coloc prop.", sum(!iseQTLColoc.topTiss)/length(iseQTLColoc.topTiss)))

# dev.off()


# write.table(mARE2hm.gARE.coloc.matrix.filt,
#             file=paste0(dir.ou, "Fig.4b.", GWAS.NM, ".haQTL.coloc.PP4.acrossTiss.hg19.tsv"),
#             sep="\t",
#             quote=F)








# 
# 
# 
# 
# # eqtl.tab.filt, 
# #      hqtl.tab.filt,
# #      rs.ovlp.r2,
# #      g.pks,
# 
