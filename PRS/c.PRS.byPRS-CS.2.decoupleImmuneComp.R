#https://choishingwan.github.io/PRS-Tutorial/plink/
#calculate PRS plink 
#and decouple the immune components based on intersection with our peak sets related to BD 
#library(ggplot2)



# cutoffs= rbind(c(0.001, 0, 0.001),
#            c(0.05, 0, 0.05),
#            c(0.1, 0, 0.1),
#            c(0.2, 0, 0.2),
#            c(0.3, 0, 0.3),
#            c(0.4, 0, 0.4),
#            c(0.5, 0, 0.5))

options(scipen=999)
library(data.table)



#CLUMP.LD.R2.CUTOFF=0.2
#CLUMP.DIS.CUTOFF=10000
PEAK.WIND=2000

#file.in.gwas.txt = "./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.forPRS-CS.txt"#"./BP_GWAS_QC/daner_PGC_BIP32b_mds7a_0416a.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.bed.gz"
file.in.PRS_CS.postEff="c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/BD_Mayo_pst_eff_a1_b0.5_phiauto_chr*.txt"
file.in.target.QC = "./b_Mayo.BP_PGC.leaveOutMayo.GWAS_MAF0.05.locID_V1.1_noSwapping_hg38_QC/plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC"#"./Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID_V1.1_noSwapping/plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC"
file.in.WGSid2phenotype = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
#file.in.meta2 = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
#file.in.clinic = "/broad/compbio/lhou/data/Mayo_Bipolar/210415_new_broad_data_deidentified.csv"
files.in.pks.bed=dir("b_genomicRegions_forPRS/", "hg19.bed", full.names=T) 
names(files.in.pks.bed)=gsub(".hg19.bed", "", basename(files.in.pks.bed))

#
dir.tmp=paste0("~/hptmp/MayoBipolar/c_PRS_CS_genomicRegions/wind", PEAK.WIND/1000, "kb/")
dir.create(dir.tmp, showWarnings = F, recursive = T)

dir.ou=paste0("c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/PRS_wind", PEAK.WIND/1000, "kb/")
dir.create(dir.ou, showWarnings = F, recursive = T)
#file.ou.gwas.logOR = "./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.logOR.bed"
file.ou.gwas.PRS_CS.postEff.hg19.bed="c_PRS_byPRS-CS_PGC.leaveOutMayo.GWAS_genomicRegions/BD_Mayo_pst_eff_a1_b0.5_phiauto_hg19.bed"
file.ou.PRS.pref = paste0(dir.ou, "plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC.vs.BDlooMAYO_PGC3.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.PRS_CS")
#file.ou.PRS.pdf = paste0(dir.ou, "plink.PRSvsphenotype.pdf")
file.ou.PRS.RData = paste0(dir.ou, "plink_PRS_CS_genomicRegions.RData")



#combine SNPs
cmd = paste0("cat ", file.in.PRS_CS.postEff, 
            "|awk 'BEGIN{OFS=\"\\t\"} {print \"chr\"$1, $3-1, $3, $2, $4, $5, $6}' ",
            ">", file.ou.gwas.PRS_CS.postEff.hg19.bed)
system(cmd)


#PRS calculation, SNP and pvalue 
#all SNPs
cmd=paste0("plink ",
    " --bfile ", file.in.target.QC,
    " --score ", file.ou.gwas.PRS_CS.postEff.hg19.bed, " 4 5 7",#" 4 5 13 header",
    " --out ", paste0(file.ou.PRS.pref, ".all"))
system(cmd)

#filter based on peaks and calculate

for(nm in names(files.in.pks.bed))
{
  print(nm)
  #filter 
  f.t.ovlp.bed=paste0(dir.tmp, "all.ovlp.", nm, ".hg19.bed")
  f.o.pks.ovlp.bed=paste0(dir.ou, nm, ".peaks.window.snp.hg19.bed")
  f.o.snp.filt = paste0(dir.ou, nm, ".snp.peaks.window.hg19.bed")

  # cmd=paste0("bedtools window -u -w ", PEAK.WIND, " -a ", file.in.gwas.gz, " -b ", files.in.pks.bed[[nm]], ">", f.t.ovlp.bed)
  # system(cmd)

  cmd=paste0("bedtools window -w ", PEAK.WIND, " -a ", file.ou.gwas.PRS_CS.postEff.hg19.bed, " -b ", files.in.pks.bed[[nm]], ">", f.t.ovlp.bed)
  system(cmd)
  df=read.table(f.t.ovlp.bed, sep="\t", head=F, row.names=NULL, stringsAsFactors=F)
  snps.filt=df[,4]
  snps.filt=snps.filt[!duplicated(snps.filt)]
  #pks.ovlp=paste0(df[,15], ":", df[,16], "-", df[,17]) #hg19
  pks.ovlp=df[,11]
  pks.ovlp=pks.ovlp[!duplicated(pks.ovlp)]
  #snps.filt=intersect(snps.clump, snps.filt)
  print(length(snps.filt))
  write(snps.filt, file=f.o.snp.filt)
  write(paste(gsub(":|-", "\t", pks.ovlp), pks.ovlp, sep="\t"), f.o.pks.ovlp.bed)

  #
  cmd=paste0("plink ",
    " --bfile ", file.in.target.QC,
    " --score ", file.ou.gwas.PRS_CS.postEff.hg19.bed, " 4 5 7",
    " --extract ", f.o.snp.filt,
    " --out ", paste0(file.ou.PRS.pref, ".filtBy.", nm))

  system(cmd)


}
print("done")


#
samp2meta=read.table(file.in.WGSid2phenotype, sep=",", header = T, row.names = 5, stringsAsFactors = F)
rownames(samp2meta)=gsub("s_", "", rownames(samp2meta))
#samp2meta$group[samp2meta$group=="case"] = "disease"
wgsId2sampID = paste0("id_", samp2meta$newsubjectid) 
names(wgsId2sampID)=rownames(samp2meta)
 


samp2PRS.list=lapply(c("all", paste0("filtBy.", names(files.in.pks.bed))),
FUN=function(nm)
{
  f.i=paste0(file.ou.PRS.pref, ".", nm, ".profile")
  buf=read.table(f.i, sep="", header=T, row.names=2)
  rownames(buf)=wgsId2sampID[rownames(buf)]
  df=data.frame(#group=samp2grp[rownames(buf), "group"], 
             buf$SCORE, 
             (buf$SCORE-mean(buf$SCORE))/sd(buf$SCORE),
             #paste0("pval_", cutoff),
             #nm,
             stringsAsFactors=F)
  colnames(df)=c(paste0(nm, ".PRS"), paste0(nm, ".PRS.z"))#, "cutoff", "type")
  rownames(df) = rownames(buf)

  return(df)
})
samp2PRS.df=do.call(cbind, samp2PRS.list)



save(samp2PRS.df, file=file.ou.PRS.RData)




