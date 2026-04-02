#V1.4
#separate cluster visualization based on up or down for each disease

#V1.3
#add background peaks when comparing fractions
#and cluster annotation together while visualization

#V1.2 add neutrophil

#v1.1 for nnPoisR results
#annotate by ChIPseeker




library(rGREAT)
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# library(clusterProfiler)
#library(corrplot)
library(pheatmap)
colfunc <- colorRampPalette(c("white", "chocolate4"))
colors= colfunc(10)

FDR.CUTOFF=0.05
Hyper_Region_Set_Coverage.CUTOFF=0.1
Hyper_Fold_Enrichment.CUTOFF=1.5

dir.in.bed = "a_3_AREModulesInEpimap_V1.6.1_merge3activeMarks/H3K27ac/" #paste0("a_3_AREModulesInEpimap_V1.5_onlyCovered/") #paste0("a_3_AREModulesInEpimap_V1.4_mergedPk4Tiss/")# paste0("a_3_AREModulesInRoadmap_V1.2_addRefPeaks/", COND, "/")

file.in.bg.peak.bed = paste0(dir.in.bed, "AREs.all.bed")
files.in.clu.bed=dir(dir.in.bed, "AREs.clu.*?bed", full.names = T)
names(files.in.clu.bed) = sapply(files.in.clu.bed,
FUN=function(f.i)
{
  gsub(paste0("AREs.|.bed"), "", basename(f.i))
})
file.in.RData= paste0(dir.in.bed, "3actHMMergedPeaks.inEpimap.H3K27acSignal.cluster.RData")

dir.ou = "a_4_modules_rGREAT_merge3activeMarks_Epimap/" #"a_4_modules_rGREAT_mergedPk4Tiss_Epimap_onlyCovered/" #"a_4_modules_rGREAT_mergedPk4Tiss_Epimap/"
dir.create(dir.ou, showWarnings = F, recursive = T)
file.ou.clu.txt = paste0(dir.ou, "modules.term.txt")
file.ou.RData = paste0(dir.ou, "modules.term.RData")
file.ou.pdf = paste0(dir.ou, "modules.term", ".q", FDR.CUTOFF, ".cov",Hyper_Region_Set_Coverage.CUTOFF, ".FE", Hyper_Fold_Enrichment.CUTOFF, ".pdf")

#

#


bg.bed = read.table(file.in.bg.peak.bed, sep="\t", row.names = NULL, header=F)[,1:3]
colnames(bg.bed) = c("chr", "start", "end")
  
enrich.res=list()
if(file.exists(file.ou.RData))
{
  load(file.ou.RData)
}
for(nm in names(files.in.clu.bed))
{
  if(nm %in% names(enrich.res)) next
  print(paste0(" ", nm))
  mod.bed =  read.table(files.in.clu.bed[nm], sep="\t", row.names = NULL, header=F)
  colnames(mod.bed) = c("chr", "start", "end")
  
  job = submitGreatJob(mod.bed, bg.bed, species = "hg19")
  tb = getEnrichmentTables(job, category=c("GO", "Phenotype", "Genes"))
  tb.filt = lapply(tb,
  FUN=function(r)
  {
    r[r$Hyper_Adjp_BH<=FDR.CUTOFF, c("name", "Hyper_Fold_Enrichment", "Hyper_Region_Set_Coverage", "Hyper_Raw_PValue")]
  })
  enrich.res[[nm]]=tb.filt
  save(enrich.res, file=file.ou.RData)
}




load(file.in.RData)

categories = c("GO Biological Process")#,  "GO Cellular Component" ,
               #"GO Molecular Function", "Human Phenotype" ,
               #"Mouse Phenotype", "Mouse Phenotype Single KO")
enrichRes.merge= lapply(categories,
FUN=function(categ)
{
  terms.all = levels(factor(unlist(lapply(enrich.res, FUN=function(res){res[[categ]]$name[res[[categ]]$Hyper_Fold_Enrichment>=Hyper_Fold_Enrichment.CUTOFF &
                                                                                      res[[categ]]$Hyper_Region_Set_Coverage>=Hyper_Region_Set_Coverage.CUTOFF]}))))
  enrichMatrix = matrix(0, ncol=ncol(ARES.modules.ordered), nrow=length(terms.all))
  rownames(enrichMatrix) = terms.all
  colnames(enrichMatrix)=colnames(ARES.modules.ordered)
  for(clu in colnames(ARES.modules.ordered))
  {
    rownames(enrich.res[[clu]][[categ]]) = enrich.res[[clu]][[categ]][,"name"]
    terms = enrich.res[[clu]][[categ]][,"name"]
    terms.ovlp = intersect(terms, terms.all)
    enrichMatrix[terms.ovlp, clu] = enrich.res[[clu]][[categ]][terms.ovlp,"Hyper_Fold_Enrichment"]
  }
  
  return(enrichMatrix)
  
})
names(enrichRes.merge)=categories
     
pdf(file.ou.pdf, width=20, height=12)

for(categ in names(enrichRes.merge))
{
  print(categ)
  matrix.log2enrich = log2(enrichRes.merge[[categ]]+1)

  
  term2clu.row.order = order(apply(matrix.log2enrich, 1, which.max), decreasing = F)
  matrix.log2enrich.ordered = matrix.log2enrich[term2clu.row.order,]

 
  pheatmap(matrix.log2enrich.ordered,
          col=colors,
         show_rownames=T,
         show_colnames=T,
         cluster_rows=F,
         cluster_cols=F,
         #fontsize_base=0.2,
         fontsize_row=7,
         fontsize_col=5,
         main=categ)


  f.o.tsv=paste0(dir.ou, "modules.", categ, ".tsv")
  write.table(matrix.log2enrich.ordered, file=f.o.tsv, sep="\t", quote=F)
}

dev.off()


save(enrich.res,
     enrichRes.merge,
     file=file.ou.RData)






