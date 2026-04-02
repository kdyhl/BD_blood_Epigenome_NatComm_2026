#use blueprint eQTL
#LD from GTEx



args = commandArgs(trailingOnly=TRUE)
CELLTYPE = args[1]  #tcel neut mono
  
if(!grepl("R version 3.6", R.version$version.string))
{
  stop("use .r-3.6.0-bioconductor")
}

R.LIB.MYPATH = "~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/"
R.LIBS.PATH = .libPaths()
.libPaths(c(R.LIB.MYPATH, R.LIBS.PATH))



options(scipen = 999)

library(data.table)
#library(pheatmap)
#library("RColorBrewer")
# library(grid)
# library(gridExtra)
library(ggplot2)
library(RColorBrewer)
MY.PALETTE <- colorRampPalette(c("royalblue2", "green1", "orange2", "darkred")) #rev(brewer.pal(11, "Spectral"))
SC_COLOR = scale_colour_gradientn(colours = MY.PALETTE(8), limits=c(0,1)) 

library(coloc)
source("~/lhou.compbio/codes/mylib/RCode/Util_zqtl_byYongjin.R")



eQTLOvlpGWAS.WIND = 10000

GWAS.P.CUTOFF = 1e-6
GWAS.Z.CUTOFF = abs(qnorm(GWAS.P.CUTOFF/2))
BP.GWAS.N = 20352+31358 #https://www.nature.com/articles/s41588-019-0397-8
BP.GWAS.CASE.PROP=20352/31358
#
GENE.WIND = 1000000
EQTL.P.CUTOFF = 1e-5
cellType2sampN = c(mono=194,
                   neut=192,
                   tcel=169)

# 


#haQTL.PVAL.CUTOFF =0.001
file.in.gene.info = "~/lhou.compbio/data/gene/human/gencode.v26.GRCh38.genes_GTEx.v8.gtf"

file.in.BP.gz = "~/lhou.compbio/data/GWAS/pgc_bipolar/pgc_bip_hg38.sorted.bed.gz"

file.in.eQTL = paste0("/broad/compbio/lhou/data/QTL/Blueprint/eQTL.bed.hg38/", CELLTYPE, "_gene_nor_combat_peer_10_all_summary.sorted.bed.gz")
#file.in.signif.eQTL = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.signif_variant_gene_pairs.txt.gz")
#file.in.eGene = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.egenes.txt.gz")

dir.in.GTEx.genotype.plink = "/broad/compbio/data/GTEx/GTEx_restricted/v8_plink/plink-relaxed/"#"~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/"
file.in.GTEx.genotype.vcf= "/broad/compbio/data/GTEx/GTEx_restricted/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz" #"~/lhou.compbio/data/Mayo_Bipolar/WGS/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID.vcf.gz"
#dir.in.1KG.genotype.plink ="/broad/compbio/data/1KG_phase3/eur/"

dir.tmp=paste0("~/hptmp/MayoBipolar/a_Blueprint.eQTLvsGWAS.coloc/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.eQTL.bed = paste0(dir.tmp, CELLTYPE, ".eQTL.bed")
file.tmp.BP.hits.bed = paste0(dir.tmp, CELLTYPE, ".BP.hits.hg38.bed")
file.tmp.eQTLGene.ovlp.BP.hits.bed = paste0(dir.tmp, CELLTYPE, ".eQTLGene.ovlp.BP.hits.hg38.bed")

#file.tmp.eQTL.chr.RDS = paste0(dir.tmp, HM, ".", CHR, ".eQTL.RDS")
#file.tmp.haQTL.chr.RDS = paste0(dir.tmp, HM, ".", CHR, ".haQTL.RDS")

#file.tmp.peak.bed = paste0(dir.tmp, HM, ".", CHR, ".g.hg38.bed")
#file.tmp.gene.bed = paste0(dir.tmp, HM, ".", CHR,".eGene.hg38.bed")

# file.tmp.haQTL.hg38.bed = paste0(dir.tmp, HM, ".", CHR,".haQTL.hg38.bed")
# file.tmp.haQTL.unmap.bed = paste0(dir.tmp, HM, ".", CHR,".haQTL.unmap.bed")

dir.ou = paste0("a_Blueprint.eQTLvsGWAS_bycoloc/", CELLTYPE, "_", GENE.WIND/1000, "k/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.pdf =paste0(dir.ou, CELLTYPE, "GWAS.P", GWAS.P.CUTOFF, ".eQTLvsGWAS.pdf")
file.ou.RData =paste0(dir.ou, CELLTYPE, "GWAS.P", GWAS.P.CUTOFF, ".eQTLvsGWAS.RData")
#


#focus on those peak overlapped with strong signals in GWAS
#load(file.in.haQTL.RData)


#
buf = read.table(file.in.gene.info, comment.char = "#", sep="\t", row.names=NULL, header = F)
colnames(buf) = c("chr", "source", "type", "start", "end", "v6", "strand", "v8", "v9")

gene.info = buf[buf[,3]=="gene", c("chr", "start", "end", "strand", "v9")]
gene.info$gene_id = sapply(gene.info$v9,
FUN=function(x)
{
  gsub("gene_id (ENSG\\d+\\..*?);.*", "\\1", x)
})
gene.info$gene_id_noVersion   = gsub("\\..*", "", gene.info$gene_id)

gene.info$gene_name = sapply(gene.info$v9,
FUN=function(x)
{
  gsub(".*gene_name (.*?);.*", "\\1", x)
})
gene.info$TSS = gene.info$start
gene.info$TSS[gene.info$strand=="-"] =  gene.info$end[gene.info$strand=="-"]
#rownames(gene.info) = gsub("\\..*", "", gene.info$gene_id)


#eQTL 

cmd=paste0("zcat ", file.in.eQTL, "| awk '($7<=", EQTL.P.CUTOFF, "){print $1 \"\\t\" $2 \"\\t\" $3 \"\\t\" $6}' >", file.tmp.eQTL.bed )
system(cmd)
# eGenes.df = read.table(gzfile(file.in.eGene), sep="\t", header=T, row.names = 1)
# ensG2gSymbol = eGenes.df[,"gene_name"]
# names(ensG2gSymbol) =  rownames(eGenes.df)
# 
# eQTLs.signif = read.table(gzfile(file.in.signif.eQTL), sep="\t", header=T, row.names = NULL)
# eQTLs.signif.bed = cbind(sapply(eQTLs.signif[,1], FUN=function(snp){ buf=unlist(strsplit(snp, split="_")); paste(buf[1], as.numeric(buf[2])-1, buf[2], sep="\t")}),
#                          eQTLs.signif[,"gene_id"],
#                          ensG2gSymbol[eQTLs.signif[,"gene_id"]]
#                          )
# write.table(eQTLs.signif.bed, file= file.tmp.eQTL.bed, sep="\t", quote=F, row.names = F, col.names = F)  
  
cmd = paste0("zcat ", file.in.BP.gz, "|awk '($10<", GWAS.P.CUTOFF, "){print }'>", file.tmp.BP.hits.bed)
system(cmd)

cmd = paste0("bedtools window -w ", eQTLOvlpGWAS.WIND, " -a ", file.tmp.eQTL.bed, " -b ", file.tmp.BP.hits.bed, ">", file.tmp.eQTLGene.ovlp.BP.hits.bed)
system(cmd)

tryCatch(
  expr = 
  {
    buf = read.table(file.tmp.eQTLGene.ovlp.BP.hits.bed, "\t", header=F, row.names = NULL)
  },
  error=function(e)
  {
    print("no overlap with eQTL and GWAS signal found")
    quit("no")
  }
  
)
eGenes.ovlpGWAS = levels(factor(gsub("\\..*", "", buf[,4])))
eGenes.ovlpGWAS=intersect(eGenes.ovlpGWAS, gene.info$gene_id_noVersion)
print(paste0(length(eGenes.ovlpGWAS), " genes overlap with GWAS"))
      
pdf(file.ou.pdf, width=12, height=8)
#layout(matrix(1:9, ncol=3))
eGenes.eQTLvsGWAS=lapply(eGenes.ovlpGWAS,
FUN=function(gene)
{
  print(gene)
  if(sum(gene.info$gene_id_noVersion==gene)>1)
  {
    reutrn(NA)
  }
  
  
  f.tmp.ovlp.snps =paste0(dir.tmp, CELLTYPE, ".", gene, ".", GENE.WIND/1000, "k.ovlpSNP.txt")
  f.tmp.ld.pre = paste0(dir.tmp, CELLTYPE, ".", gene, ".ovlpSNP")
  
  f.tmp.eQTL.z = paste0(dir.tmp, CELLTYPE, ".", gene, ".eQTL.z.txt")
  f.tmp.gwas.z = paste0(dir.tmp, CELLTYPE, ".", gene, ".gwas.z.txt")
  
  #f.o.eCaviar = paste0(d.o.gwas, GWAS.nm, ".", HM, ".", gsub(":", "-", gene), ".txt")
  #
  g.chr= gene.info[gene.info$gene_id_noVersion==gene, "chr"]
  g.TSS = gene.info[gene.info$gene_id_noVersion==gene, "TSS"]
  query = paste0(g.chr, ":", max(1, g.TSS-GENE.WIND), "-", g.TSS+GENE.WIND)
  #
  #SNP for MAF from GTEx snpid chr_loc_REF_ALT
  snp.tab = seqminer::tabix.read.table(file.in.GTEx.genotype.vcf, query)
  snp2maf = sapply(snp.tab$INFO,
  FUN=function(x)
  {
    af =as.numeric(gsub(".*AF=(.*?);.*", "\\1", x, perl=T))
    
  })
  names(snp2maf) = gsub("_b38", "", snp.tab$ID)
  
  
  #SNPs from eQTL
  #buf = unlist(strsplit(gene, ":|-", perl=T))
  #g.end = eGenes.df[gene, "gene_end"]
  eQTL.tab = seqminer::tabix.read.table(file.in.eQTL, query) 
  eQTL.tab$phenotypeID =gsub("\\..*", "", eQTL.tab$phenotypeID)
  eQTL.tab.filt = eQTL.tab[eQTL.tab$phenotypeID==gene,]
  eQTL.tab.filt = eQTL.tab.filt[!duplicated(eQTL.tab.filt$chr.pos_ref_alt),]
  rownames(eQTL.tab.filt) = eQTL.tab.filt$chr.pos_ref_alt
  eQTL.tab.filt$varbeta=(eQTL.tab.filt$std.error_of_beta)^2
  #
  eQTL.tab.filt.flip = eQTL.tab.filt
  rownames(eQTL.tab.filt.flip)  = sapply(rownames(eQTL.tab.filt),
  FUN=function(x)
  {
    buf= unlist(strsplit(x, split="_"))
    paste(buf[c(1,2,4,3)], collapse = "_")
  })
  eQTL.tab.filt.all = rbind(eQTL.tab.filt, eQTL.tab.filt.flip)
  #(-abs(eQTL.tab.filt[, "beta"])/qt(eQTL.tab.filt[, "pval"]/2, df=HM2sampN[HM]-2-HM2COVNum[HM]))^2

  
  #SNPs from GWAS
  gwas.tab = seqminer::tabix.read.table(file.in.BP.gz, query) 
  gwas.tab = gwas.tab[!duplicated(gwas.tab$stop), ]
  rownames(gwas.tab) = paste(gwas.tab$CHR, gwas.tab$stop, gwas.tab$a1, gwas.tab$a2, sep="_")
  gwas.tab.flip = gwas.tab
  gwas.tab.flip$a1= gwas.tab$a2
  gwas.tab.flip$a2= gwas.tab$a1
  rownames(gwas.tab.flip)  = paste(gwas.tab.flip$CHR, gwas.tab.flip$stop, gwas.tab.flip$a1, gwas.tab.flip$a2, sep="_")
  gwas.tab.all =rbind(gwas.tab, gwas.tab.flip)
  
  snp.loc.ovlp = intersect(names(snp2maf), intersect(rownames(eQTL.tab.filt.all), rownames(gwas.tab.all)))
  if(length(snp.loc.ovlp)==0) return(NA)
  
  
   
  
  dataset.gwas= list(beta=gwas.tab.all[snp.loc.ovlp, "lodds"], 
                     varbeta=gwas.tab.all[snp.loc.ovlp, "se"]^2,
                     N=BP.GWAS.N, 
                     type="cc",
                     s=BP.GWAS.CASE.PROP)
  dataset.eQTL = list(beta=eQTL.tab.filt.all[snp.loc.ovlp, "beta"], 
                        varbeta=eQTL.tab.filt.all[snp.loc.ovlp, "varbeta"],
                        N=cellType2sampN[CELLTYPE], type="quant")
    
  my.res <- coloc.abf(dataset1=dataset.gwas,
                  dataset2=dataset.eQTL,
                  MAF=snp2maf[snp.loc.ovlp])#,
    
    
  
  #overlaped SNPs and flapping
  
  print(paste0(length(snp.loc.ovlp), " of ", nrow(eQTL.tab.filt), " eQTL is with GWAS zscore"))
  write(gsub("_", ":", snp.loc.ovlp), f.tmp.ovlp.snps)
  
  # plot(eQTL.z[snps.locID.ovlp],
  #     gwas.z[snps.locID.ovlp], 
  #     xlab="eQTL zscore",
  #     ylab="GWAS zscore",
  #     main=paste(HM, GWAS.nm, hQTLPk, sep=" ")
  #    )
  
  #ld of overlapping SNP from eGTEx
  cmd = paste0("plink --r2 square  --bfile ", dir.in.GTEx.genotype.plink, g.chr, " --extract ",
                 f.tmp.ovlp.snps,  " --write-snplist  --out ", f.tmp.ld.pre)
  system(cmd)

  snp.loc.ovlp.kept=gsub(":","_", rownames(read.table(paste0(f.tmp.ld.pre, ".snplist"), sep="\t", row.names = 1, header = F)))
  #ld.matrix=as.matrix(read.table(paste0(f.tmp.ld.pre, ".ld"), sep="\t", row.names = NULL, header = F))
  
  ld.matrix = as.matrix(fread(paste0(f.tmp.ld.pre, ".ld"), sep="\t", header=F, data.table=F, colClasses="numeric",showProgress=T))
  
  rownames(ld.matrix)=snp.loc.ovlp.kept
  colnames(ld.matrix)=snp.loc.ovlp.kept
  snp.loc.ovlp.loc =sapply(snp.loc.ovlp, FUN=function(x) as.numeric(unlist(strsplit(x, "_"))[2]))

  if(length(snp.loc.ovlp)==0) return(NA)
  #if(!any(abs(gwas.z[snp.loc.ovlp]) >= GWAS.Z.CUTOFF)) next
  
  snp.topGWAS = snp.loc.ovlp.kept[which.min(gwas.tab.all[snp.loc.ovlp.kept, "p"])[1]]
  
  df =data.frame(log10p= c(-log10(gwas.tab.all[snp.loc.ovlp.kept, "p"]),
                           -log10(eQTL.tab.filt.all[snp.loc.ovlp.kept, "p.value"])),
                 type= c(rep("GWAS", length(snp.loc.ovlp.kept)),
                         rep("eQTL", length(snp.loc.ovlp.kept))),
                 r2.withTopGWASHit = rep(ld.matrix[snp.topGWAS,], 2),
                 is.topGWAS = rep(snp.loc.ovlp.kept==snp.topGWAS, 2),
                 loc = rep(snp.loc.ovlp.loc[snp.loc.ovlp.kept], 2))
  p=ggplot(df, aes(x=loc, y=log10p)) +
    geom_point(aes(col=r2.withTopGWASHit, shape=is.topGWAS)) +
    SC_COLOR + #scale_colour_gradientn(colours = cols) +
    ggtitle(paste0(CELLTYPE, " bipolar disorder ", gene, " ", gene.info[gene.info$gene_id_noVersion==gene,"gene_name"], "\nH4 ABF: ", signif(my.res$summary["PP.H4.abf"],3))) +
    facet_grid(type ~ .,scales="free_y") 
  print(p)

  # if(length(snp.loc.ovlp)==0) next
  # write(paste(snp.loc.ovlp, eQTL.z[snp.loc.ovlp], sep="\t"), file=f.tmp.haQTL.z)
  # write(paste(snp.loc.ovlp, gwas.z[snp.loc.ovlp], sep="\t"), file=f.tmp.gwas.z)


  
  # system("gunzip ", paste0(f.tmp.ld, ".ld.gz") )
  # 
  # #run eCaviar
  # 
  # cmd = paste0(eCAVIAR, " -o ", f.o.eCaviar,
  #              " -l ", paste0(f.tmp.ld.pre, ".ld"),
  #              " -l ", paste0(f.tmp.ld.pre, ".ld"),
  #              " -z ", f.tmp.gwas.z,
  #              " -z ", f.tmp.haQTL.z,
  #              #" -r ", 0.01,
  #              " -c 5"  
  #              )
  # 
  # write(cmd, file.ou.sh, append = T)
  return(list(qtl=df, coloc=my.res))
})
names(eGenes.eQTLvsGWAS) = eGenes.ovlpGWAS
dev.off()


save(gene.info,
     eGenes.eQTLvsGWAS, 
     file = file.ou.RData)
