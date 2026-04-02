#V1.1
#a.4.getTFMotifsEnrichedForAREModule.byHomer_V1.1_EpiMapMotif.R
#use job array to run task for epimap motifs

args = commandArgs(trailingOnly=TRUE)
MOD = args[1]

PK.WIN=0 #2000
BG= "randomBG"  #"localBG" #"#"givenBG" ##"
PROMOT.WIN=2000
CLU.SIZE.CUTOFF=50
NTASK=10

qsub.pref = c()
qsub.pref["norm"]= paste0("qsub -cwd -t 1-", NTASK, " -l h_vmem=20G -l h_rt=3:00:00 -j y -r y")
qsub.pref["rerun"]= paste0("qsub -cwd -t 1-", NTASK, " -l h_vmem=30G -l h_rt=8:00:00 -j y -r y")

# LIFTOVER = "~/lhou.compbio/software/ucscUtil/liftOver"
# FILE.CHAIN = "~/lhou.compbio/software/ucscUtil/hg38ToHg19.over.chain"

FILE.SH = "a.4.runTFMotifsEnrichForAREModlue_byHomer.sh"

file.in.ARE.Mod.RData="a_3_AREModulesInEpimap_V1.6.1_merge3activeMarks/H3K27ac/3actHMMergedPeaks.inEpimap.H3K27acSignal.RData" #"a_3_AREModulesInEpimap_V1.4_mergedPk4Tiss/4TissMergedPeaks.inEpimap.H3K27acSignal.RData"
file.in.gene.info = "~/lhou.compbio/data/gene/human/gencode.v26lift37.annotation.gtf.gz"# paste0("/broad/compbio/data/GTEx/v8/eqtl/GTEx_Analysis_v8_eQTL_expression_matrices/", tiss2GTExTiss[TISS], ".v8.normalized_expression.bed.gz")

files.in.motif= system("find ~/lhou.compbio/data/Epimap/motifs/motif_forHomer/collated_motifs_bgFreq_homerFormat.pval_1e-4.10*.txt", intern=T)
names(files.in.motif) = sapply(files.in.motif,
FUN=function(f.i)
{
  gsub(".*10_(\\d+)\\.txt", "\\1", f.i)
})

# files.in.motif= system("find ~/lhou.compbio/data/Epimap/motifs/motif_forHomer/collated_motifs_homerFormat.pval_1e-4.10*.txt", intern=T)
# names(files.in.motif) = sapply(files.in.motif,
# FUN=function(f.i)
# {
#   gsub(".*10_(\\d+)\\.txt", "\\1", f.i)
# })


dir.tmp=paste0("~/hptmp/MayoBipolar/a_4_motif_module_EpimapMotif_byHomer_V1.1_", PK.WIN/1000, "kb_", BG, "/")
dir.create(dir.tmp, showWarnings = F, recursive = T)
#file.tmp.TSS.hg38.bed= paste0(dir.tmp, "gene.TSS.hg38.bed")
file.tmp.TSS.hg19.bed= paste0(dir.tmp, "gene.TSS.hg19.bed")
#file.tmp.TSS.unmap.bed= paste0(dir.tmp, "gene.TSS.unmap.bed")
file.tmp.mARE.bg.hg19.bed= paste0(dir.tmp, "AREs.bg.hg19.bed")
#
dir.ou= paste0("./a_4_TFMotifFor3actHMsAREModule.Epimap_byHomer_V1.1_EpimapMotif_", PK.WIN/1000, "kb_", BG, "/") #"bgFreq/")
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.RData = paste0(dir.ou, "mAREs.byMod.promAndenh.RData")
  
#data preparations
load(file.in.ARE.Mod.RData)
#load(file.in.mARE.RData)

#promoter and enhancer
#
buf =read.table(gzfile(file.in.gene.info), header = F, sep="\t", row.names = NULL, stringsAsFactors = F)#, comment.char = "")
buf=buf[buf[,3]=="transcript", ]
tss.bed =paste(buf[,1], buf[,4]-1, buf[,4], sep="\t")
tss.bed[buf[,7]=="-"] =paste(buf[,1], buf[,5]-1, buf[,5], sep="\t")[buf[,7]=="-"]

write(tss.bed, file=file.tmp.TSS.hg19.bed)
# 
# cmd= paste(LIFTOVER, file.tmp.TSS.hg38.bed, FILE.CHAIN, file.tmp.TSS.hg19.bed, file.tmp.TSS.unmap.bed)
# system(cmd)


#
#AREs.all.bed = gsub(":|-", "\t", names(tiss.peak2mergedPeak[[TISS]]))
AREs.all=unlist(cluster2ARE)
AREs.all.bed = gsub(":|-", "\t", AREs.all)
write(paste(AREs.all.bed,  AREs.all, sep="\t") , file.tmp.mARE.bg.hg19.bed)

cmd =paste0("bedtools window -u -w ", PROMOT.WIN, " -a ", file.tmp.mARE.bg.hg19.bed, " -b ", file.tmp.TSS.hg19.bed)
AREs.promt.bed = gsub(":|-", "\t", sapply(system(cmd,  intern = T), FUN=function(x)
{
  buf=unlist(strsplit(x, split="\t", perl = T))[4]
}))

clu2ARE.bed = lapply(names(cluster2ARE),
FUN=function(clu.nm)
{
  pks = cluster2ARE[[clu.nm]]
  pks.bed = gsub(":|-", "\t", pks)
  res=list(promoter=intersect(pks.bed, AREs.promt.bed),
       enhancer=setdiff(pks.bed,AREs.promt.bed))
  
  
  
  return(res)
})
names(clu2ARE.bed) = names(cluster2ARE)


clu2ARE.bed= c(list(bg=list(promoter=AREs.promt.bed,
                                 enhancer=setdiff(AREs.all.bed,AREs.promt.bed))),
                         clu2ARE.bed)

save(clu2ARE.bed,
     file=file.ou.RData)

#peak file
files.tmp.bed = lapply(names(clu2ARE.bed),
FUN=function(clu.nm)
{
  print(clu.nm)
  files.tmp = c()
  
  # paste0(dir.tmp, clu.nm, ".", c("promoter", "enhancer"), ".txt")
  # names(files.tmp) =  c("promoter", "enhancer")
  
  for(type in c("promoter", "enhancer"))
  {
    #if(length())
    if(length(clu2ARE.bed[[clu.nm]][[type]])>=CLU.SIZE.CUTOFF)
    {
      d = t(sapply(clu2ARE.bed[[clu.nm]][[type]], FUN=function(x){unlist(strsplit(x, split="\t", perl=T))}))  
      files.tmp[type] =  paste0(dir.tmp, clu.nm, ".", type, ".txt")
      
      d.forHomer = cbind(d[,1], as.numeric(d[,2])-PK.WIN, as.numeric(d[,3])+PK.WIN, paste(d[,1], d[,2], d[,3], sep="_"), 1, "+")
      write.table(d.forHomer, sep="\t", quote=F, file=files.tmp[type], col.names = F, row.names = F)
      
    }
  
  }
  
  
    
  return(files.tmp)
  
})
names(files.tmp.bed) = names(clu2ARE.bed)


save(clu2ARE.bed,
     files.tmp.bed,
     file=file.ou.RData)
#
for(clu.nm in names(clu2ARE.bed)[-1])
{
  for(type in c("promoter", "enhancer"))
  {
    print(paste(clu.nm, type))
    d.o = paste(dir.ou, clu.nm, "_", type, "/", sep="")
    
    if(BG=="givenBG")
    {
    
      cmd = paste0(qsub.pref[MOD],
                  " -N ", paste0(clu.nm,"_", type),
                  " -o ", paste0(dir.ou, "out/"),
                  " -v MOD=", MOD, ",fg=", files.tmp.bed[[clu.nm]][type], ",genome=hg19,dir_ou=", d.o, ",bg=", files.tmp.bed[["bg"]][type],
                  " ", FILE.SH)
                   
      
    }else if(BG=="randomBG")
    {
      cmd = paste0(qsub.pref[MOD],
                  " -N ", paste0(clu.nm,"_", type),
                  " -o ", paste0(dir.ou, "out/"),
                  " -v MOD=", MOD, ",fg=", files.tmp.bed[[clu.nm]][type], ",genome=hg19,dir_ou=", d.o, ",bg=", BG,
                  " ", FILE.SH)
                   
    }
    
    print(cmd)
    system(cmd)
   
  }
}
# 
# #for test
# for(clu.nm in c("clu_56", "clu_87", "clu_99"))
# {
#   type = "enhancer"
# 
#   print(paste(clu.nm, type))
#   d.o = paste(dir.ou, clu.nm, "_", type, "/", sep="")
# 
#   cmd = paste0(qsub.pref[MOD],
#                " -N ", paste0(clu.nm,"_", type),
#                " -o ", paste0(dir.ou, "out/"),
#                " -v fg=", files.tmp.bed[[clu.nm]][type], ",genome=hg19,dir_ou=", d.o, ",bg=", files.tmp.bed[["bg"]][type],
#                " ", FILE.SH)
# 
#   print(cmd)
#   system(cmd)
# 
# 
# }
# 














