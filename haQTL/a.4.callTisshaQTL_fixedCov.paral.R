#V1.2.3
#a.4.callTisshaQTL_fixedCov_V1.2.3.paral.R

#V1.2
#new version processed enhQTL
#V1.2 script 

#V1.1
#qvalue multiple test correction
#consider LD too
#answer 1 queston
#from brain tissue, if it is a haQTL, it and its proximal tags are more likely to be a eQTL by qq plot across different tissues
#the statistics of it, how many haQTLs or its proxmial tags are also eQTL, fisher pvalue 


#overlap between enhQTL and eQTL
#
#for this version no LD block is considered

#MAF.cutoff=0.05
#P.filter =0.01
#P.cutoff=1e-5
# Q.cutoff=0.05
# TAG.LD.cutoff=0.8
# TAG.WIND.cutoff=1e6#K
# eQTL.P.cutoff=1e-3
# TISS= "Brain" #"Muscle" #"" #"Heart" "Lung"

#P.cutoff.stringe=1e-5


args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) 
{
  stop("need one argument for organs, any of Heart, Brain, Lung, Muscle.\n", call.=FALSE)
} 
HM = args[1]
 

CHRs = paste0("chr", c(1:22))#, "X"))
            
sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                "##$ -V",
                "#$ -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=40g",
                "#$ -l h_rt=2:00:00", 
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                sep="\n")
script= "a.4.callTisshaQTL_fixedCov.R" #"a.4.callTisshaQTL_fixedCov_V1.2.2_filtPk_gINT.Peer_rmSexAge.R" #"a.4.callTisshaQTL_fixedCov_gINT.R"




for(chr in CHRs)
{
  f.o.sh = paste0("a_4_haQTL.", HM, ".", chr, ".sh")
  
  write(sh.head, f.o.sh)
  
  cmd = paste("Rscript" , script, HM, chr, sep=" ")
  write(cmd, f.o.sh, append = T)
  
  write("echo work done", f.o.sh, append = T)
  
  system(paste("qsub ", f.o.sh, sep=""))
  
  
}
