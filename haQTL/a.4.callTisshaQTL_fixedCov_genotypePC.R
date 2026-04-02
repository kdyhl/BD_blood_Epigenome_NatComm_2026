#V1.2.3.1
#a.4.callTisshaQTL_fixedCov_V1.2.3.1_dp.bg_genotypePC.R
#add genotype PCs 

#V1.2.3
#using  background peaks for QTL calling

#V1.2.2
#rm sex and age

#V1.2.1
#gINT+peer

#V1.2
#use peer input


#V1.1
#use quantile normalization version


#with hidden factor identified from RUV-g together with genotype PCs
#call haQTL first with different number of factors to test how many factor we should use


#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) 
{
  stop("need one argument for HMs.\n", call.=FALSE)
} 
HM = args[1]
CHR = args[2] #chr1 - chrX



sh.heads = c()
sh.heads["nominal"]= paste("#!/bin/bash",
                #"#$ -binding linear:1",
                #"#$ -P compbio_lab",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=10g",
                "#$ -l h_rt=0:30:00",
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                sep="\n")
sh.heads["permutation"]= paste("#!/bin/bash",
                #"#$ -binding linear:1",
                #"#$ -P compbio_lab",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=10g",
                "#$ -l h_rt=5:00:00",
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                sep="\n")

    

options(scipen = 999)
#
HM2COVNum = c(H3K27ac = 10,
              H3K36me3 = 10,
              H3K4me1 = 10,
              H3K27me3 = 10,
              H3K4me3 = 10)

topFactor.num =  HM2COVNum[HM]


FASTQTL="~/lhou.compbio/software/FastQTL/FastQTL-2.184/bin/fastQTL.static"
MAF = 0.05

PVAL.REPORT.CUTOFF = 1
#PVAL.VIS.CUTOFF = c(1e-4)#, 1e-10, 1e-12, 1e-15)
#PVAL.VIS.CUTOFF.COL = c(rgb(1,0,0,0.2), rgb(1,0,0,0.4), rgb(1,0,0,0.6), rgb(1,0,0,0.8))
CIS.WIND.LEN = 100000 #2e6 #5e3#

MODs = c("nominal")#, "permutation")

#

#
file.in.covar = system(paste0("find ./a_haQTLCalling_covariate_V1.2.3_pkFiltByMedianCV_gINT_Peer/", HM, ".*PEER_covariates.txt"), intern = T)#a_haQTLCalling_covariate_V1.2_quant.peer
file.in.HM.bed = paste0("./a_haQTLCalling_covariate_V1.2.3_pkFiltByMedianCV_gINT_Peer/", HM, ".bed") #a_haQTLCalling_covariate_V1.2_quant.peer
file.in.meta.csv = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
file.in.wgs.sampleInfo="~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
file.in.genotypePC="~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID.QCed.PCA.eigenvec"

file.in.pks.hg38tohg19.RDS = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".pk.hg38tohg19.RDS")
file.in.allPks.bed = paste0("../peakMergeNormal/a_2_peakHeightsNormalized_V1.4/", HM, ".diffPeak.BG.bed")




#file.in.totalVCF = "/broad/compbio/data/MayoWGS/Kellis.WGS.subsetFromBiobank.vcf.gz"
file.in.totalVCF = "~/lhou.compbio/data/Mayo_Bipolar/WGS/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID.vcf.gz"

#file.tmp.inds = paste("/broad/hptmp/lhou/m6A/inds.txt", sep="")
dir.tmp="~/hptmp/MayoBipolar/a_4_haQTL_V1.2.3.1_dp_bg_gPC/" #"~/hptmp/eGTEx-H3K27ac/a_2_haQTL_V1.2_quant_peer/"
dir.create(dir.tmp, showWarnings = F, recursive = T)
file.tmp.HM.inds = paste(dir.tmp, "inds.", HM, ".", CHR, ".txt", sep="")
file.tmp.HM.vcf = paste(dir.tmp, HM, ".", CHR, ".vcf", sep="")

file.tmp.phentp = paste(dir.tmp, HM, ".", CHR, ".bed", sep="")
file.tmp.phentp.sort = paste(dir.tmp, HM, ".", CHR, ".sorted.bed", sep="")
file.tmp.phentp.sort.gz = paste(dir.tmp, HM, ".", CHR, ".sorted.bed.gz", sep="")
file.tmp.cov.pref = paste(dir.tmp, HM, ".", CHR, ".", sep="")


dir.ou=paste0("a_4_haQTL_FixedFactorNum_V1.2.3.1_DP.BG_", CIS.WIND.LEN/1000, "k_genotypePCs/")#"a_2_haQTL_testFactorNumFromRUVg_withGTExMinCov_V1.2_quant_peer
dir.create(dir.ou, showWarnings = F, recursive = T)

for(mod in MODs)
{
  dir.ou.HM=paste0(dir.ou, HM, "_", mod, "/")
  dir.create(dir.ou.HM, showWarnings = F, recursive = T)
}  

#
meta = read.table(file.in.meta.csv, sep=",", header=T, row.names=1)
rownames(meta)=paste("id", rownames(meta), sep="_")

samp.info=read.table(file.in.wgs.sampleInfo, sep=",", header = T, row.names = 1, stringsAsFactors = F)
id2WGSid=samp.info$HA.Sample.Name
names(id2WGSid)= paste0("id_", rownames(samp.info))
wgsID.meta = meta[names(id2WGSid),]
rownames(wgsID.meta) = id2WGSid

genotype.PCs.df=read.table(file.in.genotypePC, header=T, row.names=2) #genotype PC
rownames(genotype.PCs.df)=paste0("s_", rownames(genotype.PCs.df))

covar.hidden = as.matrix(read.table(file.in.covar, header=T, row.names = 1, sep="\t", check.names = F))
posAndHM = read.table(file.in.HM.bed, sep="\t", header=T, row.names=NULL, comment.char = "", check.names=F)


#filter peaks
pks.hg38to19=readRDS(file.in.pks.hg38tohg19.RDS)
buf = read.table(file.in.allPks.bed, sep="\t", header = F, row.names = NULL, stringsAsFactors = F)
pks.bg=paste0(buf[,1], ":", buf[,2], "-", buf[,3])
pks.bg.hg38= names(pks.hg38to19)[pks.hg38to19 %in% pks.bg]

#vcf prepare for fastQTL#######################################################
#
subjs.withGentp = system(paste("bcftools query -l ", file.in.totalVCF,sep=""), intern=T)

HM.samps = colnames(posAndHM)[-c(1:4)]
HM.samps.subWithGentp = intersect(HM.samps, subjs.withGentp)

write.table(HM.samps.subWithGentp, file.tmp.HM.inds, sep="", row.names = F, col.names = F, quote=F)


system(paste("bcftools view  -r ", CHR, " -S ", file.tmp.HM.inds, " -m2 -M2 -v snps -o ", file.tmp.HM.vcf, " --min-ac=", round(length(HM.samps.subWithGentp)*2*MAF), ":minor --no-update ", file.in.totalVCF, sep=""))
if(file.exists(paste0(file.tmp.HM.vcf, ".gz")))
{
  file.remove(paste0(file.tmp.HM.vcf, ".gz"))    
}
system(paste("bgzip ", file.tmp.HM.vcf, " && tabix -p vcf ", file.tmp.HM.vcf, ".gz", sep=""))


#molecula phenotypes#######################################################
#posAndH3K27ac = read.table(file.in.H3K27ac.txt, sep="\t", header=T, row.names=NULL, comment.char = "")

peakCenterPos = round((posAndHM[,"start"] +posAndHM[, "end"])/2)
posAndHM[,"start"] = peakCenterPos-1
posAndHM[,"end"] = peakCenterPos
posAndHM.filt = posAndHM[(posAndHM$ID %in% pks.bg.hg38) & (posAndHM[[1]]==CHR),
                         c("#chr", "start", "end", "ID", HM.samps.subWithGentp)]

write.table(posAndHM.filt, file=file.tmp.phentp, sep="\t", quote=F, row.names=F)#, col.names=F)
write(paste(colnames(posAndHM.filt),collapse="\t"), file.tmp.phentp.sort)
system(paste("sortBed -i ", file.tmp.phentp, ">> ", file.tmp.phentp.sort))
if(file.exists(file.tmp.phentp.sort.gz))
{
  file.remove(file.tmp.phentp.sort.gz)    
}
system(paste("bgzip ", file.tmp.phentp.sort," && tabix -p bed ", file.tmp.phentp.sort.gz, sep=""))


#covariates #######################################################
file.tmp.cov = paste0(file.tmp.cov.pref, "top", topFactor.num, "cov.txt")
write(paste(c("id", HM.samps.subWithGentp), collapse="\t"), file.tmp.cov)

#genotype top PCs and hiddeen variables
write.table(rbind(covar.hidden[1:topFactor.num, HM.samps.subWithGentp],
                    age=t(wgsID.meta)["age_atsamp", HM.samps.subWithGentp],
                    sex=t(wgsID.meta)["Gender", HM.samps.subWithGentp],
                    t(genotype.PCs.df[HM.samps.subWithGentp, paste0("PC", 1:5)])
                    ),
            file=file.tmp.cov, sep="\t", row.names=T, col.names=F, quote=F, append=T)

# write.table(rbind(cov.GTEx[1:5, HM.samps.subWithGentp],
#                   t(cov)[1:topFactor.num, HM.samps.subWithGentp]), 
#             file=file.tmp.cov, sep="\t", row.names=T, col.names=F, quote=F, append=T)
if(file.exists(paste(file.tmp.cov, ".gz", sep="")))
{
  file.remove(paste(file.tmp.cov, ".gz", sep=""))
}  
system(paste("gzip ", file.tmp.cov, sep=""))



#fastqtl############################
#

for(mod in MODs)
{

  #file.tmp.cov = paste0(file.tmp.cov.pref, "top", topFactor.num, "cov.txt")
  file.ou.sh = paste0("a_4_haQTL.", HM, ".", CHR, ".", mod, ".", "top", topFactor.num, "var.sh")
  file.ou.fastqtl = paste0(dir.ou, HM, "_", mod, "/", CHR, ".cisQTL.Pval", format(PVAL.REPORT.CUTOFF, scientific=T), ".", "top", topFactor.num, "var.txt")
  
  write(sh.heads[mod], file.ou.sh)
  write("date", file.ou.sh, append = T)
  
  
  if(mod =="nominal")
  {
    
    cmd = paste0(FASTQTL," --vcf ", file.tmp.HM.vcf, ".gz", 
              " --bed ", file.tmp.phentp.sort.gz, 
              " --cov ", file.tmp.cov, ".gz", 
              " --window ", CIS.WIND.LEN, 
              " --threshold ", PVAL.REPORT.CUTOFF, 
              " --region ", CHR, 
              " --out ", file.ou.fastqtl)
  }else
  {
    cmd = paste0(FASTQTL," --vcf ", file.tmp.HM.vcf, ".gz", 
              " --bed ", file.tmp.phentp.sort.gz, 
              " --cov ", file.tmp.cov, ".gz", 
              " --window ", CIS.WIND.LEN, 
              " --permute 1000 10000", 
              " --region ", CHR, 
              " --out ", file.ou.fastqtl)
  }
             # " --permute ", 1000, sep="")
  #system(cmd)
  
  write(cmd, file.ou.sh, append = T)
  write(paste0("gzip ", file.ou.fastqtl), file.ou.sh, append = T)
  write("echo work Done!", file.ou.sh, append = T)
  write("date", file.ou.sh, append = T)
  system(paste0("qsub ", file.ou.sh))

}




#QC and visual the results ####################################################




