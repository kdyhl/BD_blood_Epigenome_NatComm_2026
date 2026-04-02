#compare differential signals associated with BD and drug treatment

HMs = c("H3K27ac", "H3K36me3", "H3K4me1", "H3K4me3", "H3K27me3")
BD.COLOC.CUTOFF=0.5

dirs.in = c(BD ="c_bulkDiffPeak_limma_V1.4.1_sva_noEHR/",
            drugs ="c_bulkDiffPeak_limma_V1.4_sva_drugs/")

files.in.coloc.RData=paste0("../haQTLvsGWAS/a_bulkhQTLvsGWAS_bycoloc_V1.1_allGWAS/PGC.BP/", HMs, "_100k/", HMs, ".bulkhQTLvsGWAS.RData")
names(files.in.coloc.RData)=HMs

files.in.pk.hg38toHg19.RDS=paste0("a_2_peakHeightsNormalized_V1.4/", HMs, ".pk.hg38tohg19.RDS")
names(files.in.pk.hg38toHg19.RDS)=HMs


dir.ou = "d_cmpBulkDP_BDvsDrug_limma/"
dir.create(dir.ou, recursive = T, showWarnings = F)



#
pk.hg38toHg19.list=lapply(files.in.pk.hg38toHg19.RDS, readRDS)

BD.coloc.pks.hg19.list=lapply(names(files.in.coloc.RData),
FUN=function(hm)
{
  f.i.RData=files.in.coloc.RData[hm]

  load(f.i.RData)
  H4s=sapply(hQTLPeaks.QTLvsGWAS,
  FUN=function(res)
  {
    res$coloc$summary["PP.H4.abf"]
  })
  names(H4s)=names(hQTLPeaks.QTLvsGWAS)

  pks.coloc=names(H4s)[H4s>=BD.COLOC.CUTOFF]

  pks.coloc.hg19=pk.hg38toHg19.list[[hm]][pks.coloc]

  return(pks.coloc.hg19)
})
names(BD.coloc.pks.hg19.list)=names(files.in.coloc.RData)


for(hm in HMs)
{
  print(hm)
  f.o.cmp.pdf = paste0(dir.ou, hm, ".DP.cmp.hg19.pdf") 


  #BD
  f.i.RData = paste0(dirs.in["BD"], hm, "/", hm, ".diffPeak.RData")
  load(f.i.RData)
  BD.scores = res$diseasedisease
  rownames(BD.scores) =rownames(res$diseasedisease)
  

  #Drugs
  f.i.RData = paste0(dirs.in["drugs"], hm, "/", hm, ".diffPeak.RData")
  load(f.i.RData)
  drugs=names(res)[grepl("_current", names(res))]
  drug.scores.list = lapply(drugs,
  FUN=function(d)
  {
    ds = res[[d]][, c("logFC", "adj.P.Val")]
    rownames(ds) =rownames(res[[d]])
    return(ds)
  })
  names(drug.scores.list)=drugs
  

  pdf(f.o.cmp.pdf, height=9, width=6)
  layout(matrix(1:6, ncol=2, byrow=T))

  DPs.BD=rownames(BD.scores)[BD.scores$adj.P.Val<=0.05]
  PKs.coloc=BD.coloc.pks.hg19.list[[hm]]
  for(d in drugs)
  {

    DPs.drug=rownames(drug.scores.list[[d]])[drug.scores.list[[d]][,"adj.P.Val"]<=0.2]
    #DPs.drugAndBD=intersect(DPs.drug, DPs.BD)

    if(length(DPs.drug)==0) next


    #pks.ovlp=intersect(names(BD.tscores), names(drug.tscores.list[[d]]))
    pcc=signif(cor(BD.scores[DPs.BD, "logFC"], drug.scores.list[[d]][DPs.BD, "logFC"], use="pairwise.complete.obs"),2)
    mod=lm(drug.scores.list[[d]][DPs.BD, "logFC"]~BD.scores[DPs.BD, "logFC"])
    pval=signif(summary(mod)$coef[2, "Pr(>|t|)"],2)
    smoothScatter(BD.scores[DPs.BD, "logFC"], 
                drug.scores.list[[d]][DPs.BD, "logFC"], 
                xlab="logFC for differential peaks [BD vs control]",
                ylab="logFC for differential peaks for drug response",
                main=paste0(hm, ",", d, "\ndCREs N=", length(DPs.BD), "\nPCC=",pcc, " P:", pval))
    abline(h=0, lty=2)
    abline(v=0, lty=2)
    abline(lm(mod), col="red", lty=2, lwd=2)


    #
    df=cbind(BD.logFC=BD.scores[DPs.BD, "logFC"], 
             drug.logFC=drug.scores.list[[d]][DPs.BD, "logFC"])
    rownames(df)=DPs.BD
    f.o.tsv=paste0(dir.ou, hm, ".DP.cmp.hg19.tsv") 
    write.table(df, f.o.tsv, sep="\t", quote=F)

    #
    # pcc=cor(BD.scores[DPs.drug, "logFC"], drug.scores.list[[d]][DPs.drug, "logFC"], use="pairwise.complete.obs")
    # mod=lm(drug.scores.list[[d]][DPs.drug, "logFC"]~BD.scores[DPs.drug, "logFC"])
    # plot(BD.scores[DPs.drug, "logFC"], 
    #             drug.scores.list[[d]][DPs.drug, "logFC"], 
    #             xlab="logFC for differential peaks [BD vs control]",
    #             ylab="logFC for differential peaks for drug response",
    #             main=paste0(hm, ",", d, "\n drug response CREs N=", length(DPs.drug), "\nPCC=", signif(pcc,2)))
    # abline(h=0, lty=2)
    # abline(v=0, lty=2)
    # abline(lm(mod), col="red", lty=2, lwd=2)

    #
    pcc=signif(cor(BD.scores[PKs.coloc, "logFC"], drug.scores.list[[d]][PKs.coloc, "logFC"], use="pairwise.complete.obs"), 2)
    mod=lm(drug.scores.list[[d]][PKs.coloc, "logFC"]~BD.scores[PKs.coloc, "logFC"])
    pval=signif(summary(mod)$coef[2, "Pr(>|t|)"],2)
    plot(BD.scores[PKs.coloc, "logFC"], 
                drug.scores.list[[d]][PKs.coloc, "logFC"], 
                xlab="logFC for differential peaks [BD vs control]",
                ylab="logFC for differential peaks for drug response",
                main=paste0(hm, ",", d, "\n BD-coloc CREs N=", length(PKs.coloc), "\nPCC=", pcc, " P:", pval))
    abline(h=0, lty=2)
    abline(v=0, lty=2)
    abline(lm(mod), col="red", lty=2, lwd=2)


    df=cbind(BD.logFC=BD.scores[PKs.coloc, "logFC"], 
             drug.logFC=drug.scores.list[[d]][PKs.coloc, "logFC"])
    rownames(df)=PKs.coloc
    f.o.tsv=paste0(dir.ou, hm, ".BD.coloc.gCREs.cmp.tsv") 
    write.table(df, f.o.tsv, sep="\t", quote=F)

    #
    # DPs.BD=rownames(BD.scores)[BD.scores$adj.P.Val<=0.05]
    # pcc=cor(BD.scores[DPs.BD, "t"], drug.tscores.list[[d]][DPs.BD])
    # mod=lm(drug.tscores.list[[d]][DPs.BD]~BD.scores[DPs.BD, "t"])
    # smoothScatter(BD.scores[DPs.BD, "t"], 
    #             drug.tscores.list[[d]][DPs.BD], 
    #             xlab="t score for differential peaks [BD vs control]",
    #             ylab="t score for differential peaks for drug response",
    #             main=paste0(hm, ",", d, ",N=", length(DPs.BD), ", PCC=", signif(pcc,2)))
    # abline(h=0, lty=2)
    # abline(v=0, lty=2)
    # abline(lm(mod), color="red", lty=2)
    
  }
  dev.off()
}


