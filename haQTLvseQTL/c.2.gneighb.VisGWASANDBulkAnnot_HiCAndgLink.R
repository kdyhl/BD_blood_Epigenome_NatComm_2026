
#show both gLink scores and HiC links
#


#V1.1
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
options(scipen=999)
library(Sushi)
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

HMs= c("H3K27ac","H3K4me1", "H3K27me3") #"H3K36me3" "H3K4me3", 


file.in.gNeighb.HiC.RData="../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneNB.annot.RData"
file.in.link.HiC.RData="../hiC.DiffPeakandDEGandTF/a_geneNeighb_visGWASANDBulkAnnot_byGene_V1.1_bulkLink/geneCRELink.HiC.hg19.RData"

file.in.gNeighb.gLinks.RData="c_geneNeighb_visGWASANDBulkAnnot_byGene_GeneticLink_V1.1_3HM/3HM.geneNB.peaks.annot_gLinks.RData"
file.in.glink.RData="c_geneNeighb_visGWASANDBulkAnnot_byGene_GeneticLink_V1.1_3HM/3HM.gLinks.RData"


file.in.BP.gz = "~/lhou.compbio/data/GWAS/pgc_bipolar/pgc_bip_hg19.sorted.bed.gz"

file.in.geneStructure.RData = "~/lhou.compbio/data/gene/human/gencode.v26lift37.basic.annotation.forSushi_exonAndUTR.RData"

dir.tmp = "~/hptmp/MayoBipolar/c_2_geneNeighb_HiCAndgLink/"
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

  
dir.ou="c_2_geneNeighb_visGWASANDBulkAnnot_byGene_HiCAndgLinks/"
dir.create(dir.ou, showWarnings = F, recursive = T)

#file.out.pref = paste(dir.ou, "geneNB", sep="")
file.ou.gNeighbor.pdf = paste0(dir.ou, "candiGene.gNeighbor.pdf")
file.ou.gNeighbor.RData = paste0(dir.ou, "candiGene.gNeighbor.RData")


#######################################
genes.candid =list()

load(file.in.gNeighb.HiC.RData)
genes.candid$HiC = genes.filt.sorted

load(file.in.gNeighb.gLinks.RData)
genes.candid$gLinks = genes.filt.sorted

gName2gESGid = names(gESGid2gName)
names(gName2gESGid)=gESGid2gName
genes.candid$ovlp=intersect(gESGid2gName[genes.candid$gLinks], genes.candid$HiC)

load(file.in.link.HiC.RData) #$hm.gene2CRE.HiC.hg19
load(file.in.glink.RData) #hm.gLinks.hg19

geneCand2hm.CRE.linked=lapply(genes.candid$ovlp,
FUN=function(g)
{
  gid=gName2gESGid[g]
  
  g.links.hic.df=do.call(rbind, 
  lapply(names(hm.gene2CRE.HiC.hg19),
  FUN=function(hm)
  {
    if(! (g %in% names(hm.gene2CRE.HiC.hg19[[hm]]))) return(NULL)
    data.frame(HM=hm,
               CREs=hm.gene2CRE.HiC.hg19[[hm]][[g]],
               type="HiC",
               stringsAsFactors = F)
  }))
  
  g.links.gLinks.df=do.call(rbind, 
  lapply(names(hm.gLinks.hg19),
  FUN=function(hm)
  {
    if(! (gid %in% names(hm.gLinks.hg19[[hm]]))) return(NULL)
    data.frame(HM=hm,
               CREs=hm.gLinks.hg19[[hm]][[gid]],
               type="gLinks",
               stringsAsFactors = F)
  }))
  
  rbind(g.links.hic.df, g.links.gLinks.df)
})
names(geneCand2hm.CRE.linked)=genes.candid$ovlp

save(geneCand2hm.CRE.linked,
     file=file.ou.gNeighbor.RData)


print("tracks for each gene neighborhood")##########################################
# load(file.ou.gNeighbor.annot.RData)

################################################################################
load(file.in.geneStructure.RData)
write.table(transcripts.bed, file=file.tmp.genesStructure.bed, sep="\t", col.names = F, row.names = F, quote=F)

#FUNCTION##############################################################################
getChrFromRegions = function(regions) #vector of regions with format "chr\tstart\tend"
{
  chrs = sapply(regions, FUN=function(x){unlist(strsplit(x, split="\t"))[1]})
  
  return(levels(factor(chrs)))
}

# getHiCLinks = function(regions, hic.reg2neighb)#vector of regions with format "chr\tstart\tend"
# {
#   nb.links = c()
#   for(nb in intersect(regions, names(hic.reg2neighb)))
#   {
#     nb.partners = hic.reg2neighb[[nb]]
#     for(nb.p in intersect(nb.partners, nbs))
#     {
#       nb.links=c(nb.links, paste(nb, nb.p, NA, 1, ".", ".", 1, sep="\t"))
#       if(nb > nb.p)
#       {
#         nb.links=c(nb.links, paste(nb.p, nb, NA, 1, ".", ".", 1, sep="\t"))
#       }
#     }
#   }
#   nb.links = levels(factor(nb.links))
#   return(nb.links)
# }

getLinksBEDPE.hg19=function(g, hm.pks.linked.df)
{
  #g is ensemb gid
  # pks=hm.pks.linked.df$CREs
  # pks.ovlp= pks[duplicated(pks)]
  # hm.pks.linked.df$type[hm.pks.linked.df$CREs %in% pks.ovlp]="both"
  # hm.pks.linked.df=hm.pks.linked.df[!duplicated(pks), ]
  
  hm.pks.linked.df$color="#FF853F"
  hm.pks.linked.df$color[hm.pks.linked.df$type == "gLinks"]="#0096FF"
  #hm.pks.linked.df$color[hm.pks.linked.df$type == "both"]="blue"
  
  g.range= gName2range.hg19[g,]
  
  if(g.range$strand=="+")
  {
    start1=g.range$start-1
    end1 = g.range$start
  }else
  {
    start1=g.range$end
    end1 = g.range$end+1
  }
  nb.links.bedpe = data.frame(chrom1 = g.range$chr,
                              start1 = start1,
                              end1 = end1,
                              chrom2=hm.pks.linked.df$chr,
                              start2=hm.pks.linked.df$start,
                              end2=hm.pks.linked.df$end,
                              name= hm.pks.linked.df$type,
                              color=hm.pks.linked.df$color,
                              score = 1,
                              strand1 = ".",
                              strand2 =".",
                              samplenumber=1
                              )
  return(nb.links.bedpe)
}

getNBGenesBED = function(gName, g.chr, chr.start, chr.end)
{
  
  write(paste(g.chr, "\t", chr.start, "\t", chr.end, sep=""), file.tmp.gRegion.bed)
  cmd =paste( "intersectBed -a ", file.tmp.genesStructure.bed, " -b ", file.tmp.gRegion.bed," -u>", file.tmp.gStructure.bed, sep="")
  system(cmd)
  buf = read.table(file.tmp.gStructure.bed, head=F, sep="\t", row.names = NULL, stringsAsFactors = F)
  colnames(buf) = c("chrom", "start", "stop", "gene", "gene.type", "transcript", "score", "strand", "type")
  buf = buf[buf$gene.type=="protein_coding" | buf$gene==gName ,]#only keep protein coding genes and the gene we are focusing on now
  buf$strand[as.character(buf$strand)=="+"] = 1
  buf$strand[as.character(buf$strand)=="-"] = -1
  buf$strand = as.numeric(buf$strand)
  
  NB.genes.bed = list()
  NB.genes.bed$g = buf[buf$gene==gName, ]
  NB.genes.bed$others = buf[buf$gene!=gName, ]
  
  return(NB.genes.bed)
}

getAnnotBed=function(pks.df, hm.pk.hg19.annot)  
{
  annot.bed.list= lapply(1:nrow(pks.df),
  FUN=function(i)
  {
    #print(i)
    hm= pks.df$HM[i]
    pk= pks.df$CREs[i]
    pk.start= pks.df$start[i]
    pk.end= pks.df$end[i]
    type=pks.df$type[i]
    
    if(pk %in% rownames(hm.pk.hg19.annot[[hm]]))
    {
      annot.nms= colnames(hm.pk.hg19.annot[[hm]])[hm.pk.hg19.annot[[hm]][pk, ]!=0]
      return(data.frame(
                      chrom=g.chr,
                      start=pk.start,
                      end=pk.end,                        
                      CRE=pk,
                      name=paste0(type, "_", annot.nms, "_", hm),
                      score=0,
                      strand=".",
                      stringsAsFactors = F))
    }else
    {
      return(NULL)
    }
  })
  annot.bed.df=do.call(rbind, annot.bed.list)
  annot.bed.df=annot.bed.df[ grepl("_DP|_MR|_coloc", annot.bed.df$name), ]  

  annot.bed.df$name=gsub("_MR|_coloc", "_upstream", annot.bed.df$name) 
  annot.bed.df$name=gsub("_DP", "_downstream", annot.bed.df$name) 
  annot.bed.df$row = as.numeric(factor(annot.bed.df$name))
  annot.bed.df$color = maptocolors(annot.bed.df$row, col=SushiColors(6))
  annot.bed.df.sorted = annot.bed.df[order(annot.bed.df$row),]  

  return(annot.bed.df.sorted)
}
  


getGWASFiltBedgraph = function(g.haQTL.gwas, dis, g.chr)
{
  haQTL.posAndGWAS.bed = data.frame(stringsAsFactors = F)
  for(posID in names(g.haQTL.gwas[[dis]]))
  {
    
    z= abs(g.haQTL.gwas[[dis]][posID])
    pos = as.numeric(unlist(strsplit(posID, split="_"))[2])
    chr = unlist(strsplit(posID, split="_"))[1]
    haQTL.posAndGWAS.bed = rbind(haQTL.posAndGWAS.bed, data.frame(chrom = chr, start = pos-1, end = pos, value = z,stringsAsFactors = F))
    
  }
  haQTL.posAndGWAS.bed = haQTL.posAndGWAS.bed[haQTL.posAndGWAS.bed$chrom==g.chr,]
  if(nrow(haQTL.posAndGWAS.bed)>=1)
  {
    haQTL.posAndGWAS.bed =  haQTL.posAndGWAS.bed[order(haQTL.posAndGWAS.bed$start, decreasing = F),]
  }
  return(haQTL.posAndGWAS.bed)
}



#
pdf(file.ou.gNeighbor.pdf, width =25, height=12)
layout(matrix(1:5, ncol=1))
par(mar=c(3,15,2,3))
for(gene in genes.candid$ovlp)#[27:29])
{
  print(gene)
  f.o.link.tsv=paste0(dir.ou, gene, ".gNeighbor.linked.tsv")
  f.o.gwas.tsv=paste0(dir.ou, gene, ".gNeighbor.gwas.tsv")
  
  #CRE linked
  hm.pks.linked.df= geneCand2hm.CRE.linked[[gene]]
  hm.pks.linked.df=cbind(hm.pks.linked.df, 
  do.call(rbind, lapply(hm.pks.linked.df$CREs,
  FUN=function(cre)
  {
    buf=unlist(strsplit(cre, split=":|-", perl=T))
    data.frame(chr=buf[1],
               start=as.numeric(buf[2]),
               end=as.numeric(buf[3]),
               stringsAsFactors = F
               )
  })))
  
  g.chr= hm.pks.linked.df$chr[1]
  chr.start = max(1,min(hm.pks.linked.df$start)-5000)
  chr.end  = max(hm.pks.linked.df$end)+5000
  
  #filt by gwas
  gwas.tab = seqminer::tabix.read.table(file.in.BP.gz, paste0(g.chr, ":", chr.start, "-", chr.end)) 
  if(length(gwas.tab)==0) next
  gwas.bedgraph = data.frame(chr= gwas.tab$V1,
                             start = gwas.tab$V2,
                             end = gwas.tab$V3,
                             value = -log10(gwas.tab$V10) 
                             )
  if(all(gwas.bedgraph$value <= -log10(GWAS.PVAL.CUTOFF))) next
  
  
  #linking 
  #pks = rownames(gene2Neighb.scores.filt.hg19[[g]])

  links.bedpe = getLinksBEDPE.hg19(gene, hm.pks.linked.df)
  if(nrow(links.bedpe)<1)
  {
    print(paste(g, "no links around on the same chrom"))
    next
  }
  #print(g)
  
  #
  write.table(gwas.bedgraph, f.o.gwas.tsv, sep="\t", quote=F, row.names=F)
  write.table(links.bedpe, f.o.link.tsv, sep="\t", quote=F, row.names=F)
  
  
  #axis(side=2,las=2,tcl=.2)
  #mtext("",side=2,line=1.75,cex=.75,font=2)

  #neighbor regions
  # nbs.bed = data.frame(t(sapply(nbs, FUN=function(x){unlist(strsplit(x, split="\t"))})), score = 1 , strand = ".",stringsAsFactors = F)
  # colnames(nbs.bed)[1:3] = c("chrom", "chromstart", "chromend")
  # nbs.bed$chromstart = as.numeric(nbs.bed$chromstart)
  # nbs.bed$chromend = as.numeric(nbs.bed$chromend)
  # plotBed(beddata=nbs.bed, chrom = g.chr, chromstart=chr.start, chromend=chr.end, row  = "auto", wiggle=0.001)
  # labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
  #  
  #gene
  NB.genes.bed = getNBGenesBED(gene, g.chr, chr.start, chr.end)
  
  
  if(nrow(NB.genes.bed$others)==0)
  {
    plot.new()
  }else
  {
    pg = plotGenes(NB.genes.bed$others[,c("chrom", "start", "stop", "gene", "score", "strand")],
               chrom = g.chr, chromstart = chr.start,chromend=chr.end,
               #colorby = c(0,1)[(NB.transcripts.bed$gene==g) +1],
               #colorbycol=SushiColors(2),
               col="blue",
               types="exon",#NB.transcripts.bed$type, 
               maxrows=8, bheight=0.2,
               plotgenetype="box",bentline=FALSE, wigglefactor = 0.1,
               labeloffset=.4, fontsize=0.6,# arrowlength = 0.025,
               labeltext=TRUE)
  }
  
  pg = plotGenes(NB.genes.bed$g[,c("chrom", "start", "stop", "gene", "score", "strand")],
                 chrom = g.chr, chromstart = chr.start,chromend=chr.end,
                 #colorby = c(0,1)[(NB.transcripts.bed$gene==g) +1],
                 #colorbycol=SushiColors(2),
                 col="red",
                 types="exon",#NB.transcripts.bed$type, 
                 maxrows=1, bheight=0.2,
                 plotgenetype="box",bentline=FALSE, wigglefactor = 0.1,
                 labeloffset=.4, fontsize=1.2,# arrowlength = 0.025,
                 labeltext=TRUE)
  labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
  
  
  #axis(side=2,las=2,tcl=.2)
  #mtext("",side=2,line=1.75,cex=.75,font=2)
  
  #hiC
  pbpe = plotBedpe(links.bedpe, 
                   chrom = g.chr, 
                   chromstart = chr.start,
                   chromend  = chr.end,
                   color= links.bedpe$color,
                   heights = links.bedpe$score,
                   plottype="loops",
                   main=gene)
  labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
  
  #annotations bed

  annot.bed.df.sorted=getAnnotBed(hm.pks.linked.df, HM.PK.hg19.annot)
    
  plotBed(beddata = annot.bed.df.sorted,
          chrom = g.chr,
          chromstart = chr.start,
          chromend =chr.end,
          rownumber = annot.bed.df.sorted$row, 
          type = "region",
          color=annot.bed.df.sorted$color,
          row="given",
          plotbg="white",
          rowlabels=unique(annot.bed.df.sorted$name),
          rowlabelcol=unique(annot.bed.df.sorted$color), 
          rowlabelcex=0.75)
  labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
  mtext("bulk annotate", side=3, adj=-0.065,line=0.5,font=2)
  
  #GWAS
  # 
  # plotBedgraph(gwas.bedgraph,
  #              chrom = g.chr, 
  #              chromstart = chr.start,
  #              chromend  = chr.end,
  #              colorbycol= SushiColors(5))
  # 
  plot(gwas.bedgraph$end, 
       gwas.bedgraph$value,
       xlim=c(chr.start, chr.end),
       type ='p',
       bty='n',
       pch=19,
       col=c("black", "red")[(gwas.bedgraph$value >= -log10(GWAS.PVAL.CUTOFF)) +1],
       cex=c(0.5, 1)[(gwas.bedgraph$value >= -log10(GWAS.PVAL.CUTOFF)) +1],
       xaxt='n',ylab="",xlab="",xaxs="i")
  labelgenome(chrom=g.chr, chromstart= chr.start, chromend = chr.end, n=5, scale="Mb")
  mtext(" -log10 (BP GWAS P)", side=3, adj=-0.065,line=0.5,font=2)
  #axis(side=2,las=2, tcl=.2)
  
  
  
  #

  
}
dev.off()




#


