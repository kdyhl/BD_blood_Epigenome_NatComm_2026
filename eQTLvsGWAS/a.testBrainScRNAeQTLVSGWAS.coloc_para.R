#use Brain scRNA eQTL from yongjin
#LD from GTEx


sh.head=paste("#!/bin/bash",
                "#$ -binding linear:1",
                #"#$ -P compbio_lab",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=20g",
                "#$ -l h_rt=5:00:00",
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "use .r-3.6.0-bioconductor",
                sep="\n")

CELLTYPEs = c("Astro", "Ex-L4", "Ex-NRGN", "In-PV-Chandelier", "In-SST", "In-VIP", "OPC", "Ex-L2or3", "Ex-L5or6", "In-PV-Basket", "In-Rosehip", "In-SV2C", "Microglia", "Oligo")

SCRIPT = "a.testBrainScRNAeQTLVSGWAS.coloc_byCellType.R"

for(ct in CELLTYPEs)
{
  f.o.sh = paste0("a.", ct, ".coloc.sh")
  write(sh.head, f.o.sh)
    
  write(paste("Rscript", SCRIPT, ct, sep=" "), f.o.sh, append = T) #a.linkGenePeak_perGene_PRS.R
  
    
  write("echo jobs done", f.o.sh, append = T)
  system(paste0("qsub ", f.o.sh))
}