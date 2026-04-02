# Manhattan plot; Supp Figure 6b


.libPaths(c("~/lhou.compbio/libs/R/x86_64-pc-linux-gnu-library/3.6/", .libPaths()))

library(CMplot)
library(stringr)


HM2COL=c("H3K27ac"="#E5803C", 
         "H3K36me3"="#5FAA37",
         "H3K4me1"="#E974FF", 
         "H3K27me3"="#4474EA", 
         "H3K4me3"="#DD4545")

files.in.permutation.p = paste0("a_5_haQTL_multTest_FixedFactorNum_V1.2.3_DP.BG_100k/", names(HM2COL), ".permutation.all.txt")
names(files.in.permutation.p) = names(HM2COL)


dir.ou="b_manhattanPlot/"
dir.create(dir.ou, recursive=T, showWarnings=F)
file.ou.RDS=paste0(dir.ou, "HMs.hQTL.empP.RDS")


hQTL.list=lapply(names(files.in.permutation.p),
FUN=function(hm)
{
  f.i=files.in.permutation.p[hm]
  d=read.table(f.i, sep=" ", head=F, row.names = 1, stringsAsFactors = F)
  snps=d[,5]
  snps.buf=str_split_fixed(snps, "_", 4)
  data.frame(SNP=paste(snps.buf[,1], snps.buf[,2], sep="_"),
                #chr= snps.buf[,1],
                #position = as.numeric(snps.buf[,2]),
                emp.p=d[,10],
                stringsAsFactors=F
                )
})
names(hQTL.list) = names(files.in.permutation.p)

hQTL.df=hQTL.list[[1]]
for(hm in names(hQTL.list)[-1])
{
  hQTL.df = merge(hQTL.df, hQTL.list[[hm]], all=T, by="SNP", suffixes=c("",paste0(".",hm)))
}
snps.buf=str_split_fixed(hQTL.df$SNP, "_", 2)
hQTL.df=data.frame(SNP=hQTL.df$SNP,
                        chr=as.numeric(gsub("chr", "", snps.buf[,1])),
                        pos=as.numeric(snps.buf[,2]),
                        hQTL.df[,-1],
                        stringsAsFactors=F)
colnames(hQTL.df)[4:8] = names(files.in.permutation.p)
hQTL.sorted.df=hQTL.df[order(hQTL.df$chr, decreasing=F),]


saveRDS(hQTL.sorted.df, file.ou.RDS)
f.o.tsv=paste0(dir.ou, "HMs.hQTL.empP.tsv")
write.table(hQTL.sorted.df, file=f.o.tsv, sep="\t", quote=F, row.names=F)

CMplot(hQTL.sorted.df,#[sort(sample(1:nrow(hQTL.sorted.df), 10000, replace=F)),], 
        plot.type="m",
        multracks=T,
        col=c("gray"),#("grey30","grey60"), #c
        threshold=c(0.005),
        threshold.lty=c(1),
        cex=.1,
        threshold.lwd=c(1), 
        threshold.col=c("black"), 
        amplify=TRUE,
        chr.den.col=NULL, 
        signal.col=HM2COL,
        signal.cex=c(.4),
        file="pdf",
        memo="HMs.merged.manhattan2",
        trait.legend.cex=.5,
        trait.legend.pos="middle",
        #dpi=300,
        file.output=TRUE,
        verbose=TRUE)






