

#
HM = c("H3K27ac", "H3K4me1", "H3K36me3", "H3K4me3", "H3K27me3")
TOP.FACTOR.NUM = c(1, 2, 5, 10, 15, 20, 25, 30)
PVAL.CUTOFFs  = c(1e-3, 1e-4, 1e-5)
DIS.CUTOFF =10000#0

dir.in = "a_2_haQTL_testFactorNum_V1.2.3_gINT_peer_rmAgeSex/" #  "a_2_haQTL_testFactorNumFromRUVg_geneINT_10k/" # "a_2_haQTL_testFactorNumFromRUVg_withGTExMinCov_V1.1.1_geneINT/" # "a_2_haQTL_testFactorNum_withGTExMinCov_V1.3_quant_norm_peer/"#"a_2_haQTL_testFactorNumFromRUVg_withGTExMinCov_V1.2_quant_peer/"#"a_2_haQTL_testFactorNumFromRUVg_withGTExMinCov_V1.1_quant/" #"a_2_haQTL_testFactorNumFromRUVg_withGTExMinCov/"
dir.ou = "a_3_haQTLPeakvsCovNum/"
dir.create(dir.ou, showWarnings = F, recursive = T)

file.ou.pdf = paste0(dir.ou, "a_2_haQTL_testFactorNum_V1.2.3_gINT_peer_rmAgeSex_", DIS.CUTOFF/1000, "kb.pdf")


pdf(file.ou.pdf, width=9, height=15)
layout(matrix(1:15, ncol=3, byrow = T))

for(hm in HM)
{
  for(p.cutoff in PVAL.CUTOFFs)
  {
    haQTLPeaksNum = sapply(TOP.FACTOR.NUM,
    FUN=function(x)
    {
      cmd =paste0("cat ", dir.in, hm, "/chr*top", x, "var.txt|awk '(sqrt($3 ^2)<= ", DIS.CUTOFF, "&& $4<=", p.cutoff, "){print $1}'|sort|uniq|wc -l" )
      return(as.numeric(system(cmd, intern=T)))
      
    })
    plot(TOP.FACTOR.NUM, haQTLPeaksNum, type="b", main=paste0(hm, "\nnominal pvalue <=", p.cutoff),
         xlab = "number of hidden covariates included", ylab="number of peaks with QTL")
  }
}

dev.off()



