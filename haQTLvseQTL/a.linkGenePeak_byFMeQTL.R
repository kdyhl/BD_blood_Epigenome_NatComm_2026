#V1.1
#a.linkGenePeak_byFMeQTL_V1.1_gARE.R
#gARE
#by distance to FM-eQTL
#

args = commandArgs(trailingOnly=TRUE)
TISS = args[1] #Brain Heart Muscle Lung
HM = args[2]
FLANK.WINDOW = as.numeric(args[3])

TISS2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
                  Lung = "Lung",
                  Muscle = "Muscle_Skeletal",
                  Heart = "Heart_Left_Ventricle",
                  Blood = "Whole_Blood")



#file.in.eQTL2pk.RData=paste0("./a_AREsNeareQTL/", COND, "_10k/", HM, ".gene2peaks_byeQTL.hg38.RData")
file.in.g2pk.RData = paste0("./a_gAREsNeareGene_emp.p.fdr.cutoff0.2/", TISS, "_1000k/", HM, ".eGene2gARE.hg38.RData")

file.in.finemap.eQTL = "/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_finemapping/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz"


dir.tmp= paste0("~/hptmp/MayoBipolar/a_genePeakLink_byFMeQTLDis/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.pks = paste0(dir.tmp, TISS,  ".", HM, ".g.pks.bed")
file.tmp.FMeQTL = paste0(dir.tmp, TISS, ".", HM, ".FMeQTL.bed")
file.tmp.pks.OVLP.FMeQTL = paste0(dir.tmp, TISS, ".", HM, ".OVLP.FMeQTL.bed")

dir.ou = paste0("a_genePeakLink_byFMeQTL_V1.1_gARE/", TISS, "_", FLANK.WINDOW/1000, "k/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData = paste0(dir.ou, TISS, ".", HM, ".gen2peaks.byFMeQTL.RData")
#file.ou.pdf = paste0(dir.ou, "gen2peaks.byFMeQTL.pdf")


#
load(file.in.g2pk.RData)

gAndPk.list = lapply(names(eGene2gAREs.hg38), 
FUN=function(g)
{
  pks = eGene2gAREs.hg38[[g]]
  as.matrix(cbind(pks, g))
})
gAndPk.matrix = do.call(rbind, gAndPk.list)
gAndPk.bed = paste0(gsub(":|-", "\t", gAndPk.matrix[,1]), "\t", gAndPk.matrix[,1], "\t", gAndPk.matrix[,2])
write(gAndPk.bed, file.tmp.pks)


#
data = read.table(gzfile(file.in.finemap.eQTL), sep="\t", header =T, row.names=NULL, stringsAsFactors = F)
fmeQTL.tiss = data[data$TISSUE==TISS2GTExTiss[TISS],]
fmeQTL.tiss.bed= paste0("chr", fmeQTL.tiss$CHROM, "\t", fmeQTL.tiss$POS-1, "\t", fmeQTL.tiss$POS, "\t", fmeQTL.tiss$eQTL, "\t", fmeQTL.tiss$GENE)
write(fmeQTL.tiss.bed, file=file.tmp.FMeQTL)

#
cmd = paste0("bedtools window -w ", FLANK.WINDOW, " -a ", file.tmp.FMeQTL, " -b ", file.tmp.pks, ">", file.tmp.pks.OVLP.FMeQTL)
system(cmd)

buf = read.table(file.tmp.pks.OVLP.FMeQTL, sep="\t", head=F, row.names = NULL, stringsAsFactors = F)
buf.filt = buf[buf[,5] == buf[,10],]
dist.pk2FMeQTL= apply(cbind(abs(buf.filt[,7]-buf.filt[,3]), abs(buf.filt[,8]-buf.filt[,3])), 1, min)
minDist.pk2FMeQTL = sapply(split(dist.pk2FMeQTL, paste(buf.filt[,5], buf.filt[,9], sep=";")), min)

  
geneAndPeak.FMeQTL_DistInv.df= data.frame(gene=gAndPk.matrix[,"g"], 
                                       peak=gAndPk.matrix[,"pks"], 
                                       FMeQTL.distInv = 0, stringsAsFactors = F) #distance inverse
rownames(geneAndPeak.FMeQTL_DistInv.df) = paste(geneAndPeak.FMeQTL_DistInv.df[,1], geneAndPeak.FMeQTL_DistInv.df[,2], sep=";")
geneAndPeak.FMeQTL_DistInv.df[names(minDist.pk2FMeQTL), "FMeQTL.distInv"] = 1/minDist.pk2FMeQTL
 #
save(geneAndPeak.FMeQTL_DistInv.df, 
     file=file.ou.RData)



