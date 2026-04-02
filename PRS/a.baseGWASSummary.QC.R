#https://choishingwan.github.io/PRS-Tutorial/base/


TOTAL_N_PROP.CUTOFF = 0.9#(20352+31358)*0.9
MAF.CUTOFF=0.01
INFO.CUTOFF = 0.8


file.in.gwas= "~/lhou.compbio/data/GWAS/PGC_leaveMayoBDout/BDlooMAYO_PGC3_hg38.sorted.bed.gz"#"~/lhou.compbio/data/GWAS/pgc_bipolar/daner_PGC_BIP32b_mds7a_0416a.hg38.sorted.bed.gz"

dir.ou="./a_BP_PGC_leaveMayoBDout_GWAS_QC/"
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.filt.gwas = paste0(dir.ou, "BDlooMAYO_PGC3_hg38.sorted.maf", MAF.CUTOFF, ".INFO", INFO.CUTOFF, ".NProp", TOTAL_N_PROP.CUTOFF, ".bed")
file.ou.rmDup.gwas = paste0(dir.ou, "BDlooMAYO_PGC3_hg38.sorted.maf", MAF.CUTOFF, ".INFO", INFO.CUTOFF, ".NProp", TOTAL_N_PROP.CUTOFF, ".noDup.bed")
file.ou.dup.snp = paste0(dir.ou, "duplicated.snp")
file.ou.rmAmbg.gwas = paste0(dir.ou, "BDlooMAYO_PGC3_hg38.sorted.maf", MAF.CUTOFF, ".INFO", INFO.CUTOFF, ".NProp", TOTAL_N_PROP.CUTOFF, ".noDup.noAmbg.bed")


# Standard GWAS QC
N=as.numeric(system(paste0("zcat ", file.in.gwas, "|awk 'BEGIN{a= 0}{if (a<$12+0) a=$12} END{print a}' "), intern=T))
#system(paste0("zcat ", file.in.gwas, "|awk '(NR==1 ||$7 >= ", MAF.CUTOFF, " && $8 >= ", INFO.CUTOFF, " && $12>=", TOTAL_N_PROP.CUTOFF,"){print}'>", file.ou.filt.gwas))
system(paste0("zcat ", file.in.gwas, "|awk '(NR==1 ||$11 >= ", MAF.CUTOFF, " && $11 <=", 1-MAF.CUTOFF,
              " && $10 >= ", INFO.CUTOFF, " && $12>=", N*TOTAL_N_PROP.CUTOFF,"){print}'>", file.ou.filt.gwas))


#remove duplicate SNP
system(paste0("awk '{ print $4}' ", file.ou.filt.gwas, "|sort |uniq -d > ", file.ou.dup.snp))
system(paste0("cat ", file.ou.filt.gwas, "|grep -vf ", file.ou.dup.snp, ">", file.ou.rmDup.gwas))


#ambiguous SNP
system(paste0("awk '!( ($5==\"A\" && $6==\"T\") || ($5==\"T\" && $6==\"A\") || ($5==\"C\" && $6==\"G\") || ($5==\"G\" && $6==\"C\") ){ print}' ", file.ou.rmDup.gwas, "> ", file.ou.rmAmbg.gwas))

system(paste0("gzip ", file.ou.rmAmbg.gwas))

#sample overlap


#