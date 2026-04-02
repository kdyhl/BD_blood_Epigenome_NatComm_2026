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

cutoffs= rbind(
          #c(0.001, 0, 0.001),
           c(0.05, 0, 0.05),
           c(0.1, 0, 0.1)
           #c(0.2, 0, 0.2),
           #c(0.3, 0, 0.3),
           #c(0.4, 0, 0.4),
           #c(0.5, 0, 0.5)
           )

CLUMP.LD.R2.CUTOFF=0.4
CLUMP.DIS.CUTOFF=100000
PEAK.WIND=2000

file.in.gwas.gz = "./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.bed.gz"#"./BP_GWAS_QC/daner_PGC_BIP32b_mds7a_0416a.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.bed.gz"
file.in.target.QC = "./b_Mayo.BP_PGC.leaveOutMayo.GWAS_MAF0.05.locID_V1.1_noSwapping_hg38_QC/plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC"#"./Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID_V1.1_noSwapping/plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC"
file.in.WGSid2phenotype = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
#file.in.meta2 = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"
#file.in.clinic = "/broad/compbio/lhou/data/Mayo_Bipolar/210415_new_broad_data_deidentified.csv"
files.in.pks.bed=dir("b_genomicRegions_forPRS/", "hg38.bed", full.names=T) 
names(files.in.pks.bed)=gsub(".hg38.bed", "", basename(files.in.pks.bed))

#
dir.tmp=paste0("~/hptmp/MayoBipolar/c_PRS_genomicRegions/ld", CLUMP.LD.R2.CUTOFF, ".clumpDis", CLUMP.DIS.CUTOFF/1000, "kb.wind", PEAK.WIND/1000, "kb/")
dir.create(dir.tmp, showWarnings = F, recursive = T)

dir.ou=paste0("c_PRS_byPlink_PGC.leaveOutMayo.GWAS_genomicRegions/ld", CLUMP.LD.R2.CUTOFF, ".clumpDis", CLUMP.DIS.CUTOFF/1000, "kb.wind", PEAK.WIND/1000, "kb/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.gwas.logOR = "./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.logOR.bed"
file.ou.target.QC.clump.snp = paste0(dir.ou, "plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC.clump.valid.snps")
file.ou.target.QC.clump.GWAS.bed = paste0(dir.ou, "plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC.clump.valid.GWAS.bed")

file.ou.base.snp.p = paste0(dir.ou, "BDlooMAYO_PGC3.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.P")
file.ou.p.range= paste0(dir.ou, "P_range.list")
file.ou.PRS.pref = paste0(dir.ou, "plink.maf0.01.hwe1e-06.geno0.1.mind0.1.QC.vs.BDlooMAYO_PGC3.hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.")
#file.ou.PRS.pdf = paste0(dir.ou, "plink.PRSvsphenotype.pdf")
file.ou.PRS.RData = paste0(dir.ou, "plink_PRS_genomicRegions.RData")
file.ou.PRS.tsv = paste0(dir.ou, "plink_PRS_genomicRegions.tsv")
#Update Effect Size by log odds
dat <- fread(file.in.gwas.gz, header=T, sep="\t", data.table=F)  #beta = log(OR
if(sum(grepl("beta", names(dat),  ignore.case=T))==0)
{
  dat$BETA <- log(dat$OR)
}else
{
  names(dat)[grepl("beta", names(dat),  ignore.case=T)]="BETA"
}
#names(dat)[ names(dat) =="SNP"] = "rs" #recoganize rs as default
if(!file.exists(file.ou.gwas.logOR))
{
  write.table(dat, file.ou.gwas.logOR, quote=F, row.names=F)
}


#clumping
#may need LD estimated from elsewhere
cmd =paste0("plink  --bfile ", file.in.target.QC, 
            " --clump-p1 1",
            " --clump-r2 ", CLUMP.LD.R2.CUTOFF, 
            " --clump-kb ", CLUMP.DIS.CUTOFF/1000, 
            " --clump ", file.ou.gwas.logOR, #By default, variant IDs are expected to be in the 'rs' column. You can change this with the --clump-snp-field flag, which takes a space-delimited sequence of field names to search for. With multiple field names, earlier names take precedence over later ones.
#By default, p-values are expected to be in the 'P' column; change this with --clump-field. This has the same semantics as --clump-snp-field
            #" --clump-snp-field rs", 
            # --clump-field P", 
            " --out ", file.in.target.QC)
system(cmd)


#SNPs
cmd= paste0("awk 'NR!=1{print $3}' ", file.in.target.QC, ".clumped >  ", file.ou.target.QC.clump.snp)
system(cmd)

#PRS calculation, SNP and pvalue 
#system(paste0("awk '{print $4,$11}' ", file.ou.gwas.logOR, " > ", file.ou.base.snp.p))
system(paste0("awk '{print $4,$7}' ", file.ou.gwas.logOR, " > ", file.ou.base.snp.p))
write.table(cutoffs, file=file.ou.p.range, sep=" ", row.names = F, col.names = F, quote=F)



#all SNPs
cmd=paste0("plink ",
    " --bfile ", file.in.target.QC,
    " --score ", file.ou.gwas.logOR, " 4 5 8 header",#" 4 5 13 header",
    " --q-score-range ", file.ou.p.range, " ", file.ou.base.snp.p,
    " --extract ", file.ou.target.QC.clump.snp,
    " --out ", paste0(file.ou.PRS.pref, ".all"))

system(cmd)

#filter based on peaks and calculate
snps.clump=read.table(file.ou.target.QC.clump.snp, stringsAsFactors=F)[,1]
snps.clump.bed=dat[dat$SNP %in% snps.clump,]
write.table(snps.clump.bed, file=file.ou.target.QC.clump.GWAS.bed, row.names=F, sep="\t", quote=F)

for(nm in names(files.in.pks.bed))
{
  print(nm)
  #filter 
  f.t.ovlp.bed=paste0(dir.tmp, "all.ovlp.", nm, ".hg38.bed")
  file.ou.target.QC.clump.snp.filt=paste0(file.ou.target.QC.clump.snp, ".filtBy.", nm)
  file.ou.pks.ovlp.bed=paste0(dir.ou, "peaks.window.index.snp.", nm, ".hg38.bed")

  # cmd=paste0("bedtools window -u -w ", PEAK.WIND, " -a ", file.in.gwas.gz, " -b ", files.in.pks.bed[[nm]], ">", f.t.ovlp.bed)
  # system(cmd)

  cmd=paste0("bedtools window -w ", PEAK.WIND, " -a ", file.ou.target.QC.clump.GWAS.bed, " -b ", files.in.pks.bed[[nm]], ">", f.t.ovlp.bed)
  system(cmd)
  df=read.table(f.t.ovlp.bed, sep="\t", head=F, row.names=NULL, stringsAsFactors=F)
  snps.filt=df[,4]
  snps.filt=snps.filt[!duplicated(snps.filt)]
  pks.ovlp=paste0(df[,15], ":", df[,16], "-", df[,17]) #hg38
  pks.ovlp=pks.ovlp[!duplicated(pks.ovlp)]
  #snps.filt=intersect(snps.clump, snps.filt)
  print(length(snps.filt))
  write(snps.filt, file=file.ou.target.QC.clump.snp.filt)
  write(paste(gsub(":|-", "\t", pks.ovlp), pks.ovlp, sep="\t"), file.ou.pks.ovlp.bed)


  snps.filt=read.table(f.t.ovlp.bed, sep="\t", head=F, row.names=NULL,stringsAsFactors=F)[,4]
  snps.filt=intersect(snps.clump, snps.filt)
  print(length(snps.filt))
  write(snps.filt, file=file.ou.target.QC.clump.snp.filt)

  #
  cmd=paste0("plink ",
    " --bfile ", file.in.target.QC,
    " --score ", file.ou.gwas.logOR, " 4 5 8 header",#" 4 5 13 header",
    " --q-score-range ", file.ou.p.range, " ", file.ou.base.snp.p,
    " --extract ", file.ou.target.QC.clump.snp.filt,
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
 


samp2PRS.list=lapply(cutoffs[,3],
FUN=function(cutoff)
{

  res=lapply(c("all", paste0("filtBy.", names(files.in.pks.bed))),
  FUN=function(nm)
  {
    f.i=paste0(file.ou.PRS.pref, ".", nm, ".", cutoff, ".profile")
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
  #names(res)=c("all", paste0("filtBy.", names(files.in.pks.bed)))
  res.df=do.call(cbind, res)
  #res.df=cbind(res.df, samp2grp[rownames(res.df), c("group", "clu.corHclust5")])
  #res.df=res.df[!is.na(res.df$group),]
  res.df$cutoff=paste0("pval_", cutoff)

  return(res.df)
})
names(samp2PRS.list)=paste0("pval_", cutoffs[,3])

save(samp2PRS.list, file=file.ou.PRS.RData)



#

samp.prs.df=samp2PRS.list$pval_0.05
samp.prs.df=samp.prs.df[,grepl("PRS.z", colnames(samp.prs.df))]
colnames(samp.prs.df)=gsub("filtBy.", "", colnames(samp.prs.df))
colnames(samp.prs.df)=gsub(".PRS.z", "", colnames(samp.prs.df))
write.table(samp.prs.df, file.ou.PRS.tsv, sep="\t", quote=F)
