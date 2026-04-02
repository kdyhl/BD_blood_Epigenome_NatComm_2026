#V1.4.1.2
#c.getDiff_byLimma_V1.4.1.2_sva_NoEHR_cmpBatch_paral.R

#use sva to consider unwanted variables
#use permutated batches

args = commandArgs(trailingOnly=TRUE)
HM = args[1]


sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=20g",
                "#$ -l h_rt=00:30:00", 
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                sep="\n")
script= "c.getDiff_byLimma_NoEHR_cmpBatch.R" #"a.4.callTisshaQTL_fixedCov_V1.2.2_filtPk_gINT.Peer_rmSexAge.R" #"a.4.callTisshaQTL_fixedCov_gINT.R"

for(batch.i in 1:70)
{
  f.o.sh = paste0("c.getDiff.V1.4.1.2_", HM, "_", batch.i, ".sh")
  
  write(sh.head, f.o.sh)
  
  cmd = paste("Rscript" , script, HM, batch.i, sep=" ")
  write(cmd, f.o.sh, append = T)
  
  write("echo work done", f.o.sh, append = T)
  
  system(paste("qsub ", f.o.sh, sep=""))
  
  
}
