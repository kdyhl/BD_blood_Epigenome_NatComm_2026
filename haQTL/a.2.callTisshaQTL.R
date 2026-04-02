#V1.2.3
#a.2.callTisshaQTL_V1.2.3_filtPeaksByMedianCV_gINT.Peer_rmSexAge.R
#rm sex and age

#V1.2.2
#rm batch and age

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
HM = args[1]#Brain 
CHR = args[2] 



sh.head = paste("#!/bin/bash",
                "#$ -binding linear:1",
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

                
# 
# R.LIBS.MYPATHs = c(V3.5="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.5/",
#                 V3.4="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.4/",
#                 V3.3="~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.3/")
# R.LIBS.PATH = .libPaths()
# 
# if(grepl("R version 3.5", R.version$version.string))
# {
#   .libPaths(c(R.LIBS.MYPATHs["V3.5"], R.LIBS.PATH))
# }
# if(grepl("R version 3.4", R.version$version.string))
# {
#   .libPaths(c(R.LIBS.MYPATHs["V3.4"], R.LIBS.PATH))
# }
# if(grepl("R version 3.3", R.version$version.string))
# {
#   .libPaths(c(R.LIBS.MYPATHs["V3.3"], R.LIBS.PATH))
# }


options(scipen = 999)
#

FASTQTL="~/lhou.compbio/software/FastQTL/FastQTL-2.184/bin/fastQTL.static"
MAF = 0.05

PVAL.REPORT.CUTOFF = 0.001
#PVAL.VIS.CUTOFF = c(1e-4)#, 1e-10, 1e-12, 1e-15)
#PVAL.VIS.CUTOFF.COL = c(rgb(1,0,0,0.2), rgb(1,0,0,0.4), rgb(1,0,0,0.6), rgb(1,0,0,0.8))
CIS.WIND.LEN = 10000 # 1e6
TOP.FACTOR.NUM = c(1, 2, 5, 10, 15, 20, 25, 30)

#


file.in.covar = system(paste0("find ./a_haQTLCalling_covariate_V1.2.3_pkFiltByMedianCV_gINT_Peer/", HM, ".*PEER_covariates.txt"), intern = T)#a_haQTLCalling_covariate_V1.2_quant.peer
file.in.HM.bed = paste0("./a_haQTLCalling_covariate_V1.2.3_pkFiltByMedianCV_gINT_Peer/", HM, ".bed") #a_haQTLCalling_covariate_V1.2_quant.peer
file.in.meta.csv = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
file.in.wgs.sampleInfo="~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"



#file.in.totalVCF = "/broad/compbio/data/MayoWGS/Kellis.WGS.subsetFromBiobank.vcf.gz"
file.in.totalVCF = "~/lhou.compbio/data/Mayo_Bipolar/WGS/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.vcf.gz"


#file.tmp.inds = paste("/broad/hptmp/lhou/m6A/inds.txt", sep="")
dir.tmp="~/hptmp/MayoBipolar/a_haQTLcov_V1.2.3_gINT.peer/" #"~/hptmp/eGTEx-HM/a_2_haQTL_V1.2_quant_peer/"
dir.create(dir.tmp, showWarnings = F, recursive = T)

file.tmp.HM.inds = paste(dir.tmp, "inds.", HM, ".", CHR, ".txt", sep="")
file.tmp.HM.vcf = paste(dir.tmp, HM, ".", CHR, ".vcf", sep="")

file.tmp.phentp = paste(dir.tmp, HM, ".", CHR, ".bed", sep="")
file.tmp.phentp.sort = paste(dir.tmp, HM, ".", CHR, ".sorted.bed", sep="")
file.tmp.phentp.sort.gz = paste(dir.tmp, HM, ".", CHR, ".sorted.bed.gz", sep="")
file.tmp.cov.pref = paste(dir.tmp, HM, ".", CHR, ".", sep="")


dir.ou="a_2_haQTL_testFactorNum_V1.2.3_gINT_peer_rmAgeSex/"#"a_2_haQTL_testFactorNumFromRUVg_withGTExMinCov_V1.2_quant_peer
dir.create(dir.ou, showWarnings = F, recursive = T)

dir.ou.HM=paste(dir.ou, HM, "/", sep="")
dir.create(dir.ou.HM, showWarnings = F, recursive = T)


file.ou.sh.pref = paste0("a_2_haQTL.", HM, ".", CHR, ".")
file.ou.fastqtl.pref = paste(dir.ou.HM, CHR, ".cisQTL.Pval", format(PVAL.REPORT.CUTOFF, scientific=T), ".", sep="")
#files.ou.fastqtl.pref$RNA = paste(dir.ou.HM, ".", CHR, ".RNA_exon.cisQTL.Pval", format(PVAL.REPORT.CUTOFF, scientific=T), ".", sep="")

#file.ou.qtl.pdf = paste(dir.ou.HM, CHR, ".fastqtl.QC.pdf", sep="")
#file.ou.qtl.rm15.pdf = paste(dir.ou.tiss, ".", CHR, "Pval", format(PVAL.VIS.CUTOFF, scientific=T), ".rmTop15PCs.fastqtl.QC.pdf", sep="")

#
#load(file.in.covar.RData)
meta = read.table(file.in.meta.csv, sep=",", header=T, row.names=1)
rownames(meta)=paste("id", rownames(meta), sep="_")

samp.info=read.table(file.in.wgs.sampleInfo, sep=",", header = T, row.names = 1, stringsAsFactors = F)
id2WGSid=samp.info$HA.Sample.Name
names(id2WGSid)= paste0("id_", rownames(samp.info))
wgsID.meta = meta[names(id2WGSid),]
rownames(wgsID.meta) = id2WGSid

covar.hidden = as.matrix(read.table(file.in.covar, header=T, row.names = 1, sep="\t", check.names = F))
posAndHM = read.table(file.in.HM.bed, sep="\t", header=T, row.names=NULL, comment.char = "", check.names=F)

#vcf prepare for fastQTL#######################################################
#
subjs.withGentp = system(paste("bcftools query -l ", file.in.totalVCF,sep=""), intern=T)

HM.samps = colnames(posAndHM)[-c(1:4)]
HM.samps.subWithGentp = intersect(HM.samps, subjs.withGentp)

write.table(HM.samps.subWithGentp, file.tmp.HM.inds, sep="", row.names = F, col.names = F, quote=F)



system(paste("bcftools view  -r ", CHR, " -S ", file.tmp.HM.inds, " -m2 -M2 -v snps -o ", file.tmp.HM.vcf, " --min-ac=", round(length(HM.samps.subWithGentp)*2*MAF), ":minor --no-update ", file.in.totalVCF, sep="")) #only keep biallelic SNPS
if(file.exists(paste0(file.tmp.HM.vcf, ".gz")))
{
  file.remove(paste0(file.tmp.HM.vcf, ".gz"))    
}
system(paste("bgzip ", file.tmp.HM.vcf, " && tabix -p vcf ", file.tmp.HM.vcf, ".gz", sep=""))


#molecula phenotypes#######################################################
#posAndHM = read.table(file.in.HM.bed, sep="\t", header=T, row.names=NULL, comment.char = "")

peakCenterPos = round((posAndHM[,"start"] +posAndHM[, "end"])/2)
posAndHM[,"start"] = peakCenterPos-1
posAndHM[,"end"] = peakCenterPos
posAndHM.filt = posAndHM[,c("#chr", "start", "end", "ID", HM.samps.subWithGentp)]

write.table(posAndHM.filt[posAndHM.filt[[1]]==CHR,], file=file.tmp.phentp, sep="\t", quote=F, row.names=F)#, col.names=F)
write(paste(colnames(posAndHM.filt),collapse="\t"), file.tmp.phentp.sort)
system(paste("sortBed -i ", file.tmp.phentp, ">> ", file.tmp.phentp.sort))
if(file.exists(file.tmp.phentp.sort.gz))
{
  file.remove(file.tmp.phentp.sort.gz)    
}
system(paste("bgzip ", file.tmp.phentp.sort," && tabix -p bed ", file.tmp.phentp.sort.gz, sep=""))

#covariates #######################################################
#GTEx.minCov = as.matrix(read.table(file.in.GTExMinCov, header=T, row.names = 1, sep="\t", check.names = F))


for(topFactor.num in TOP.FACTOR.NUM)
{
  file.tmp.cov = paste0(file.tmp.cov.pref, "top", topFactor.num, "cov.txt")
  write(paste(c("id", HM.samps.subWithGentp), collapse="\t"), file.tmp.cov)
  
  #genotype top PCs and hiddeen variables
  write.table(rbind(#GTEx.minCov[, HM.samps.subWithGentp],
                    covar.hidden[1:topFactor.num, HM.samps.subWithGentp],
                    age=t(wgsID.meta)["age_atsamp", HM.samps.subWithGentp],
                    sex=t(wgsID.meta)["Gender", HM.samps.subWithGentp]
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

}

#fastqtl############################
#
for(topFactor.num in TOP.FACTOR.NUM)
{
  file.tmp.cov = paste0(file.tmp.cov.pref, "top", topFactor.num, "cov.txt")
  file.ou.sh = paste0(file.ou.sh.pref, "top", topFactor.num, "var.sh")
  file.ou.fastqtl = paste0(file.ou.fastqtl.pref, "top", topFactor.num, "var.txt")
  
  write(sh.head, file.ou.sh)
  
  write("date", file.ou.sh, append = T)
  cmd = paste0(FASTQTL," --vcf ", file.tmp.HM.vcf, ".gz", 
              " --bed ", file.tmp.phentp.sort.gz, 
              " --cov ", file.tmp.cov, ".gz", 
              " --window ", CIS.WIND.LEN, 
              " --threshold ", PVAL.REPORT.CUTOFF, 
              " --region ", CHR, 
              " --out ", file.ou.fastqtl)
             # " --permute ", 1000, sep="")
  #system(cmd)
  
  write(cmd, file.ou.sh, append = T)
  write("echo work Done!", file.ou.sh, append = T)
  write("date", file.ou.sh, append = T)
  system(paste0("qsub ", file.ou.sh))
}



#QC and visual the results ####################################################




