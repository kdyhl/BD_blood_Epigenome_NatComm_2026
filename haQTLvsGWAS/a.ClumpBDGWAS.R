#clump the GWAS signal 
# in hg38
#

library(data.table)

options(scipen=999)

GWAS.P.CUTOFF=1e-6
CLUMP.LD.R2.CUTOFF=0.2
CLUMP.DIS.CUTOFF=250

file.in.GWAS= "~/lhou.compbio/data/GWAS/PGC_leaveMayoBDout/BDlooMAYO_PGC3_hg38.sorted.bed.gz" #"~/lhou.compbio/data/GWAS/pgc_bipolar/pgc_bip_hg38.sorted.bed.gz"#"~/lhou.compbio/data/GWAS/AD/AD_Bellenguez_NG2022/GCST90027158_buildGRCh38.tsv.gz"
files.in.BDMayo.genotype=dir("~/lhou.compbio/data/Mayo_Bipolar/WGS/plink_locID_hg38/", "bim", full.names=T) #hg19
files.in.BDMayo.genotype=gsub(".bim", "", files.in.BDMayo.genotype)
names(files.in.BDMayo.genotype)=basename(files.in.BDMayo.genotype)
# # files.in.GTEx.genotype=dir("/broad/compbio/data/GTEx/GTEx_restricted/v8_plink/plink/", "bed", full.names=T)
# # files.in.GTEx.genotype=gsub(".bed", "", files.in.GTEx.genotype)
# # names(files.in.GTEx.genotype)=basename(files.in.GTEx.genotype)
# file.in.GTEx.lookup= "/broad/compbio/data/GTEx/GTEx_restricted/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"

dir.ou="a_GWAS_clump/BD_PGC_bip/"
dir.create(dir.ou, showWarnings=F, recursive=T)
file.ou.filt.GWAS=paste0(dir.ou, "BDlooMAYO_PGC3_hg38.sorted.P", GWAS.P.CUTOFF, ".tsv")
file.ou.clump=paste0(dir.ou, "BDlooMAYO_PGC3_hg38.R2", CLUMP.LD.R2.CUTOFF, ".dis", CLUMP.DIS.CUTOFF, "kb")
file.ou.clump.bed=paste0(dir.ou, "BDlooMAYO_PGC3_hg38.P", GWAS.P.CUTOFF, ".R2", CLUMP.LD.R2.CUTOFF, ".dis", CLUMP.DIS.CUTOFF, "kb.hg38.bed")


#GWAS with GTEx variants id
gwas.df=fread(file.in.GWAS, data.table=F, sep="\t", head=T,
             colClasses=c("character", rep("numeric",2), 
                            rep("character", 3),
                            rep("numeric", 8)))
gwas.pFilt.df=gwas.df[gwas.df$P<=GWAS.P.CUTOFF,]
colnames(gwas.pFilt.df)[1]="chr"
gwas.pFilt.df$loc_id=paste(gwas.pFilt.df$chr, gwas.pFilt.df$stop, gwas.pFilt.df$a1, gwas.pFilt.df$a2, sep="_")
gwas.pFilt.flip.df=gwas.pFilt.df
gwas.pFilt.flip.df$loc_id=paste(gwas.pFilt.flip.df$chr, gwas.pFilt.flip.df$stop, gwas.pFilt.flip.df$a2, gwas.pFilt.flip.df$a1, sep="_")
gwas.pFilt.cmb.df=rbind(gwas.pFilt.df, gwas.pFilt.flip.df)
#variants in 1KG
snps.BDMayo=unlist(lapply(files.in.BDMayo.genotype,
FUN=function(f.i.pref)
{
    snps=read.table(paste0(f.i.pref, ".bim"), sep="\t", header=F, row.names=NULL, stringsAsFactors=F)[,2]
    # snps.dup=snps[duplicated(snps)]
    # print(length(snps.dup))
    # setdiff(snps, snps.dup)

}))
gwas.pFilt.cmb.filt.df=gwas.pFilt.cmb.df[gwas.pFilt.cmb.df$loc_id %in% snps.BDMayo, ]
#colnames(gwas.pFilt.cmb.filt.df)[10]="P"
colnames(gwas.pFilt.cmb.filt.df)[4]="rs"
colnames(gwas.pFilt.cmb.filt.df)[15]="SNP" #for clumping
rownames(gwas.pFilt.cmb.filt.df)=gwas.pFilt.cmb.filt.df$SNP
write.table(gwas.pFilt.cmb.filt.df, file=file.ou.filt.GWAS, sep="\t", quote=F, row.names=F)

#clump 
#clumping
#may need LD estimated from elsewhere
for(chr in paste0("chr", 1:22))
{
    print(chr)
    cmd =paste0("plink  --bfile ", files.in.BDMayo.genotype[chr], 
            " --clump-p1 1",
            " --clump-r2 ", CLUMP.LD.R2.CUTOFF, 
            " --clump-kb ", CLUMP.DIS.CUTOFF, #kb
            " --clump ", file.ou.filt.GWAS, #By default, variant IDs are expected to be in the 'SNP' column. You can change this with the --clump-snp-field flag, which takes a space-delimited sequence of field names to search for. With multiple field names, earlier names take precedence over later ones.
#By default, p-values are expected to be in the 'P' column; change this with --clump-field. This has the same semantics as --clump-snp-field
          #  " --clump-snp-field rs", 
          #  " --clump-field p", 
            " --out ", paste0(file.ou.clump, ".", chr))
    system(cmd)  
}


#summary
write("#chr\tstart\tend\tSNP.loc_ID\trsID\ttagged.leadSNP\tP\tbeta\tse\tA1AF", file.ou.clump.bed)
files.in.clump=system(paste0("find ", file.ou.clump, "*.clumped"), intern=T)
for(f.i in files.in.clump)
{
    print(f.i)
    d=read.table(f.i, header=T, row.names=NULL, stringsAsFactors=F)

    d.list=lapply(1:nrow(d),
    FUN=function(i)
    {
        snp.lead = d[i, "SNP"]
        snps= unlist(strsplit(d[i, "SP2"], split=","))
        snps=setdiff(snps, "NONE")
        snps=gsub("(1)", "", snps, fixed=T)

        snps.all=c(snp.lead, snps)

        df=data.frame(
            chr=gwas.pFilt.cmb.filt.df[snps.all, "chr"],
            start=gwas.pFilt.cmb.filt.df[snps.all, "start"],
            end=gwas.pFilt.cmb.filt.df[snps.all, "stop"],
            SNP.loc_ID=gwas.pFilt.cmb.filt.df[snps.all, "SNP"],
            rsID=gwas.pFilt.cmb.filt.df[snps.all, "rs"],
            tagged.leadSNP=snp.lead,
            P=format(gwas.pFilt.cmb.filt.df[snps.all, "P"], scientific = T),
            beta=gwas.pFilt.cmb.filt.df[snps.all, "BETA"],
            se=gwas.pFilt.cmb.filt.df[snps.all, "SE"],
            A1AF=gwas.pFilt.cmb.filt.df[snps.all, "A1AF"],
            stringsAsFactors=F
            )
    })
    d.df=do.call(rbind, d.list)

    write.table(d.df, file.ou.clump.bed, sep="\t", quote=F, row.names=F, col.names=F, append=T)

}
d=read.table(file.ou.clump.bed, sep="\t", head=T, row.names=NULL, stringsAsFactors=F, comment.char="")
