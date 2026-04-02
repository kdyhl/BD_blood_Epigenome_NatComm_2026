#V1.1
#b.testEpiDEGvsCMapSignatures_subGrpDP_V1.1_paral.R

sh.head = paste("#!/bin/bash",
                "#$ -S /bin/bash",
                "#$ -V -cwd",
                "#$ -e ./error.$JOB_NAME.$JOB_ID",
                "#$ -o ./outpt.$JOB_NAME.$JOB_ID",
                "#$ -l h_vmem=60g",
                "#$ -l h_rt=3:00:00", 
                "##$ -pe smp 4",
                "source ~/.bashrc",
                "source ~/.my.bashrc",
                "export LD_LIBRARY_PATH=$LIBRARY_PATH:/broad/uge/8.4.3/lib/lx-amd64",
                "export PATH=\"/broad/compbio/lhou/software/anaconda3_2023/envs/cMap/bin/:$PATH\"",
                sep="\n")
script= "b.testEpiDEGvsCMapSignatures_subGrpDP.py" #"a.4.callTisshaQTL_fixedCov_V1.2.2_filtPk_gINT.Peer_rmSexAge.R" #"a.4.callTisshaQTL_fixedCov_gINT.R"
JOB.N=50
EXP.N=49550
EXP.N.PER.JOB= EXP.N/JOB.N

files.in.degs=c(inflamSubGrps.dnDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/inflamSubgrps.DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
            inflamSubGrps.upDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/inflamSubgrps.DEG.metaDP.DEGs.up.gid_byDAVID.txt",
            nonInflamSubGrps.dnDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/nonInflamSubgrps.DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
            nonInflamSubGrps.upDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.5conds_subGrpDP/nonInflamSubgrps.DEG.metaDP.DEGs.up.gid_byDAVID.txt")
#dir.ou = "b_comp2CompoundSig_cMap_inflamSubGrps/"

# files.in.degs=c(inflamSubGrps.dnDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.1_subGrpDP_orderedDEGs/inflamSubgrps.DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
#             inflamSubGrps.upDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.1_subGrpDP_orderedDEGs/inflamSubgrps.DEG.metaDP.DEGs.up.gid_byDAVID.txt",
#             nonInflamSubGrps.dnDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.1_subGrpDP_orderedDEGs/nonInflamSubgrps.DEG.metaDP.DEGs.dn.gid_byDAVID.txt", 
#             nonInflamSubGrps.upDEG="a_3_DEG_meta.geneBodyAndgeneRegDom_V1.2.1_subGrpDP_orderedDEGs/nonInflamSubgrps.DEG.metaDP.DEGs.up.gid_byDAVID.txt")
# dir.ou = "b_comp2CompoundSig_cMap_inflamSubGrps/a3_V1.2.1_orderedDEGs/"
dir.ou = "b_comp2CompoundSig_cMap_inflamSubGrps_V1.1/a3_V1.2_subGrps/"
dir.create(dir.ou, showWarnings=F, recursive=T)
file.ou.gmt=paste0(dir.ou, "EpiDEGs.gmt")

#read in DEG list

gsets=t(sapply(names(files.in.degs),
FUN=function(nm)
{
  df= read.table(files.in.degs[nm], sep="\t", header=T, row.names=NULL, quote="")
  gids=df[,"To"]
  gids=gids[!duplicated(gids)]
  gids.quoted=sapply(gids, FUN=function(gid)paste0("\"", gid, "\""))
  return(c(nm, 
           "na", 
           paste(gids.quoted, collapse="\t")))
}))
write.table(gsets, file=file.ou.gmt, sep="\t", quote=F, col.names=F, row.names=F)



for (i in 1:JOB.N)
{
  i.start= round((i-1)*EXP.N.PER.JOB) #in python style
  i.end=round(i*EXP.N.PER.JOB)  #in python style since i.end is not included
  
  f.o.tsv=paste0(dir.ou, "/GSEA_EpiDEGs.vs.cMap_level5_beta_trt_cp_blood_", i.start, '-', i.end, ".1000perm.tsv")
  if(file.exists(f.o.tsv))
  {
    next
  }
  f.o.sh = paste0("b.testEpiDEGvsCMap_subGrpDP_V1.1.", i, ".sh")
  write(sh.head, f.o.sh)
  print(paste(i, i.start, i.end))

  

  cmd = paste("python" , script, dir.ou, i.start, i.end, sep=" ")
  write(cmd, f.o.sh, append = T)
  
  write("echo work done", f.o.sh, append = T) 
  
  system(paste("qsub ", f.o.sh, sep=""))
  
} 
