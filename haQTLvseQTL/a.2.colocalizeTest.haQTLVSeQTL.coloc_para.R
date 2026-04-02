#V1.3
#a.2.colocalizeTest.haQTLVSeQTL.coloc_V1.3.gARE_para.R
#use gARE around genes

#check only enhancers near significant eQTL
#adding norm or rerun 


sh.heads = c()
sh.heads["norm"]=paste("#!/bin/bash",
                "#$ -binding linear:1",
                #"#$ -P compbio_lab",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=20g",
                "#$ -l h_rt=2:00:00",
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "use .r-3.6.0-bioconductor",
                sep="\n")
sh.heads["rerun"]=paste("#!/bin/bash",
                "#$ -binding linear:1",
                #"#$ -P compbio_lab",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=30g",
                "#$ -l h_rt=6:00:00",
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "use .r-3.6.0-bioconductor",
                sep="\n")


tiss2GTExTiss = c(Brain = "Brain_Frontal_Cortex_BA9",
                  Lung = "Lung",
                  Muscle = "Muscle_Skeletal",
                  Heart = "Heart_Left_Ventricle",
                  Blood = "Whole_Blood")

CHRs = paste0("chr", c(1:22))

options(scipen = 999)
eQTL2Enh.WIND.SIZE = 10000

SCRIPT = "a.2.colocalizeTest.haQTLVSeQTL.coloc_perGene.R"

argv = commandArgs(trailingOnly = TRUE)
#GENO.DIR    = argv[1]

MOD = argv[1]
HM = argv[2]
TISS = "Blood" #argv[3]
hQTL.WIND = 100000 #as.numeric(argv[4]) #the window for coloclaization test 
JOBs.N = as.numeric(argv[3])


#file.in.eGene.pos = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.egenes.txt.gz")
#file.in.sig.eQTL = paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_significant/", tiss2GTExTiss[TISS], ".v8.signif_variant_gene_pairs.txt.gz")

# file.in.hg382hg19.RDS =paste0("/broad/compbio/lhou/codes/eGTEX-H3K27ac/peakMergeNormal/a_2_peakHeightsNormalized_hg19.narrowPeak.RSC0.8.RdsTot1e7.logQ2_V1.3_addRefPeak_CSFilt/", TISS, ".pks.hg382hg19Uniq.RDS") 
file.in.g2pk.RData = paste0("./a_gAREsNeareGene_emp.p.fdr.cutoff0.2/", TISS, "_1000k/", HM, ".eGene2gARE.hg38.RData")


dir.tmp =paste0("~/hptmp/MayoBipolar/a.2_haQTLvseQTL_coloc_V1.3/", TISS, "_", HM, "/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
# file.tmp.sigeQTL.bed = paste0(dir.tmp, TISS, ".eQTL.hg38.bed")
# file.tmp.pks.hg38.bed = paste0(dir.tmp, TISS, ".peaks.hg38.bed")
# file.tmp.eQTLOVLPpks.hg38.bed = paste0(dir.tmp, TISS, ".eQTLOVLPpeaks.hg38.bed")



for(JOB.I in 1:JOBs.N)
{
  JOB.I.inFN = paste(c(rep(0, nchar(JOBs.N)-nchar(JOB.I)), JOB.I), collapse="")
  file.tmp.coloc.RData  = paste0(dir.tmp, TISS, ".", HM, ".coloc_bycoloc_", JOB.I.inFN,".RData")
  
  if(MOD == "rerun" && file.exists(file.tmp.coloc.RData))
  {
    next
  }
  
  f.o.sh = paste0("a.2.coloc.", TISS, ".", HM, ".bycoloc.", JOB.I.inFN, ".sh")
  write(sh.heads[MOD], f.o.sh)
    
  write(paste("Rscript", SCRIPT, TISS, HM, file.in.g2pk.RData, file.tmp.coloc.RData, hQTL.WIND, JOB.I, JOBs.N, sep=" "), f.o.sh, append = T) #a.linkGenePeak_perGene_PRS.R
  
    
  write("echo jobs done", f.o.sh, append = T)
  system(paste0("qsub ", f.o.sh))
  
}
