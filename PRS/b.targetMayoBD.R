#V1.1
#b.targetMayoBD.V1.1.noSwap.QC.R

#https://choishingwan.github.io/PRS-Tutorial/target/
#as discussed with the Author
#allele swapping is wrong and not necessary, since downstream software would take care of it by couning effect allele from GWAS
#only work on complementary and recoding+ complementary
#to make sure the SNP id is the same in bim 

MAF.CUTOFF=0.01
HWE.CUTOFF = 1e-6
GENO.CUTOFF =0.1
MIND.CUTOFF = 0.1

#dir.in.plink = "~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/"
file.in.vcf = "~/lhou.compbio/data/Mayo_Bipolar/WGS/Kellis.WGS.subsetFromBiobank.hg38.snp.MAF0.05.locID.vcf.gz"
#file.in.gwas.gz = "./BP_GWAS_QC/daner_PGC_BIP32b_mds7a_0416a.hg38.sorted.maf0.01.INFO0.8.N46539.noDup.noAmbg.bed.gz"
file.in.gwas.gz = "./a_BP_PGC_leaveMayoBDout_GWAS_QC/BDlooMAYO_PGC3_hg38.sorted.maf0.01.INFO0.8.NProp0.9.noDup.noAmbg.bed.gz"
file.in.WGSid2indID = "~/lhou.compbio/data/Mayo_Bipolar/Batch_Feb.2019/broad_wgs_metadata.csv"
file.in.id2sex = "~/lhou.compbio/data/Mayo_Bipolar/Bipolar_ptDemographics.csv"

dir.ou =paste0("./b_Mayo.BP_PGC.leaveOutMayo.GWAS_MAF0.05.locID_V1.1_noSwapping_hg38_QC/")
dir.create(dir.ou, showWarnings = F, recursive = T)
#file.ou.filt.pref = paste0(dir.ou, "plink.maf", MAF.CUTOFF, ".hwe", HWE.CUTOFF, ".geno", GENO.CUTOFF, ".mind", MIND.CUTOFF)
file.ou.tmp.pref = paste0(dir.ou, "plink.maf", MAF.CUTOFF, ".hwe", HWE.CUTOFF, ".geno", GENO.CUTOFF, ".mind", MIND.CUTOFF, ".tmp")
file.ou.QC.pref = paste0(dir.ou, "plink.maf", MAF.CUTOFF, ".hwe", HWE.CUTOFF, ".geno", GENO.CUTOFF, ".mind", MIND.CUTOFF, ".QC")
#filt by MAF, HWE, geno, mind, 
cmd = paste0("plink --vcf ", file.in.vcf, " --maf ", MAF.CUTOFF, " --hwe ", HWE.CUTOFF, " --geno ", GENO.CUTOFF, " --mind ", MIND.CUTOFF, " --make-bed --write-snplist --out ", file.ou.tmp.pref)
system(cmd)

#SNP prune
cmd=paste0("plink --bfile ", file.ou.tmp.pref,  " --indep-pairwise 200 50 0.25  --out ", file.ou.QC.pref) 
system(cmd)

#Heterozygosity rates 
cmd = paste0("plink --bfile ", file.ou.tmp.pref, " --extract ", file.ou.QC.pref, ".prune.in --het --out ", file.ou.QC.pref)
system(cmd)

dat <- read.table(paste0(file.ou.QC.pref, ".het"), header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], paste0(file.ou.tmp.pref,".valid.sample"), quote=F, row.names=F) # print FID and IID for valid samples



#mismatching  SNPs###################################
#1. Load the bim file, the summary statistic and the QC SNP list into R
bim <- read.table(paste0(file.ou.tmp.pref, ".bim"), stringsAsFactors = F)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")  #A1 usually alternative, A2 usually major
# Read in QCed SNPs
qc <- read.table(paste0(file.ou.tmp.pref, ".snplist"), header = F, stringsAsFactors = F)

# Read in the GWAS data
gwas = read.table(gzfile(file.in.gwas.gz),
            header = T,
            stringsAsFactors = F, 
            comment.char = "",
            sep="\t")
# Change all alleles to upper case for easy comparison
colnames(gwas)[4:6] = c("rs", "A1", "A2") #A1 is alternative, corresponding to the effect size calculation

gwas$A1 <- toupper(gwas$A1)
gwas$A2 <- toupper(gwas$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)

#
# Merge summary statistic with target
colnames(gwas)[1] ="CHR"
gwas$CHR = gsub("chr", "", gwas$CHR)
colnames(gwas)[3] ="BP"

#Identify SNPs that require strand flipping
info <- merge(bim, gwas, by = c("CHR", "BP")) #due to different SNP ID, one is rsID, the other is locus id
# Filter QCed SNPs
info <- info[info$SNP %in% qc$V1,]
# Function for finding the complementary allele
complement <- function(x) {
    switch (
        x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}

#SNP name converted to rsid based on GWAS
locID2gwasID=bim$SNP
names(locID2gwasID)=bim$SNP #those cannot be matched by default
locID2gwasID[info$SNP] = info$rs #those could be matched
bim$SNP = locID2gwasID[bim$SNP]

# Get SNPs that have the same alleles across base and target
info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
snps.match <- bim$SNP %in% info.match$rs


#2. Identify SNPs that are complementary between base and target
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
# Update the complementary alleles in the bim file
# This allow us to match the allele in subsequent analysis
snps.complement <- bim$SNP %in% info.complement$rs
bim[snps.complement,]$B.A1 <-
    sapply(bim[snps.complement,]$B.A1, complement)
bim[snps.complement,]$B.A2 <-
    sapply(bim[snps.complement,]$B.A2, complement)

#3. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)
# identify SNPs that need recoding
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
# Update the recode SNPs
snps.recode <- bim$SNP %in% info.recode$rs
# tmp <- bim[recode.snps,]$B.A1
# bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
# bim[recode.snps,]$B.A2 <- tmp

# identify SNPs that need recoding & complement
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update the recode + strand flip SNPs
snps.crecode<- bim$SNP %in% info.crecode$rs
# bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
# bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
bim[snps.crecode,]$B.A1 <- as.character(sapply(bim[snps.crecode,]$B.A1, complement)) #do not swap allele
bim[snps.crecode,]$B.A2 <- as.character(sapply(bim[snps.crecode,]$B.A2, complement))

# Output updated bim file
write.table(
    bim,
    paste0(file.ou.QC.pref, ".PGC.leaveOutMayo.adj.bim"),
    quote = F,
    row.names = F,
    col.names = F,
    sep="\t"
)

#4. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)
mismatch <-
    bim$SNP[!(snps.match |
                snps.complement | 
                snps.recode |
                snps.crecode)]
write.table(
    mismatch,
    paste0(file.ou.QC.pref, ".PGC.leaveOutMayo.mismatch"),
    quote = F,
    row.names = F,
    col.names = F
)

#5  Replace EUR.bim with EUR.QC.adj.bim:
system(paste0("mv ", file.ou.tmp.pref, ".bim ", file.ou.tmp.pref, ".bim.bk"))
system(paste0("mv ", file.ou.QC.pref, ".PGC.leaveOutMayo.adj.bim ", file.ou.tmp.pref, ".bim"))




## Sex chromosomes#########################
#A sex check can be performed in PLINK, in which individuals are called as females if their X chromosome homozygosity estimate (F statistic) is < 0.2 and as males if the estimate is > 0.8.
#Before performing a sex check, pruning should be performed
cmd =paste0("plink  --bfile ", file.ou.tmp.pref, 
            " --extract ", file.ou.QC.pref, ".prune.in", 
            " --keep ", file.ou.tmp.pref, ".valid.sample",
            " --check-sex",
            " --out ", file.ou.QC.pref)
system(cmd)

#
buf =read.table(file.in.WGSid2indID, sep=",", row.names = NULL, header = T, stringsAsFactors = F)
meta.wgsID2id = paste0("id_", buf[,"newsubjectid"])
names(meta.wgsID2id) = buf[,"HA.Sample.Name"]
#
buf =read.table(file.in.id2sex, sep=",", row.names = NULL, header = T, stringsAsFactors = F)
id2sex =buf$Gender
names(id2sex)=paste0("id_", buf[["X...newsubjectid"]])
wgsID2sex = id2sex[meta.wgsID2id]
names(wgsID2sex) = gsub("s_", "", names(meta.wgsID2id))
  
dat <- read.table(paste0(file.ou.QC.pref, ".sexcheck"), header=T, stringsAsFactors = F)
valid <- read.table(paste0(file.ou.tmp.pref, ".valid.sample"), header=T)
is.valid = (dat$F>=0.5 & wgsID2sex[dat$IID]=="M") | (dat$F<=0.2 & wgsID2sex[dat$IID]=="F") 
valid <- subset(dat, is.valid)
write.table(valid[,c("FID", "IID")], paste0(file.ou.QC.pref, ".valid"), row.names=F, col.names=F, sep="\t", quote=F) 


#relatedness
cmd =paste0("plink  --bfile ", file.ou.tmp.pref, 
            " --extract ", file.ou.QC.pref, ".prune.in", 
            " --keep ", file.ou.QC.pref, ".valid",
            " --rel-cutoff 0.125", 
            " --out ", file.ou.QC.pref)
system(cmd)


#final step

cmd =paste0("plink  --bfile ", file.ou.tmp.pref, 
            " --make-bed",
            " --keep ", file.ou.QC.pref, ".rel.id",
            #" --extract ", file.ou.tmp.pref, ".snplist",
            " --exclude ", file.ou.QC.pref, ".PGC.leaveOutMayo.mismatch",
            " --out ", file.ou.QC.pref)
system(cmd)

