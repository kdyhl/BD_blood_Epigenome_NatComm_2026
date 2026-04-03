#updated on 4.23.2018
#add another version of peak heights calcualtion: getPeakHeightsMeanFromTagAlign

#updated on 3.28.2018 add calculate counts on defined regions from bw
#functions used for H3K27ac analysis



#
getPeakHeightsMeanFromTagAlign =function(file.in.tagAlign.gz, frag.length, file.in.sorted.ref.bed, file.out, ref.region.ColId, regions2Len)
{
  
  cmd= paste("zcat ", file.in.tagAlign.gz, "|awk '{if($6==\"+\") {print $1 \"\\t\" $2 \"\\t\" $2+", frag.length, "} if($6==\"-\") {start =1>$3-", frag.length, "?1:$3-", frag.length, ";print $1 \"\\t\" start \"\\t\" $3 }}'|sortBed -i stdin| intersectBed -b stdin -a ", file.in.sorted.ref.bed, " -wo -sorted >", file.out, sep="")
  
  system(cmd)
  
  
  
  buf=read.table(file=file.out, sep="\t", header=F, row.names=NULL)
  #buf=t(sapply(buf, FUN=function(x){unlist(strsplit(x, split="\t"))}))
  buf = aggregate(as.numeric(buf[,ncol(buf)]), by=list(buf[,ref.region.ColId]), FUN=sum)
  counts=buf[,2]
  names(counts)= buf[,1]
  
  
  pk2Heights = rep(0, length(ref.regions))
  names(pk2Heights) = names(regions2Len)
  pk2Heights[names(counts)]= counts
  pk2Heights = pk2Heights/regions2Len
    
  return(pk2Heights)
}


getOvlpCounts_WeightedByOvlpFract = function(f.sorted.test, test.column.no, f.sorted.ref, ref.regions, f.o) #need sorted for both bed
{
  f.o.gz=paste(f.o, ".gz", sep="")
  if(file.exists(f.o.gz))
  {
    file.remove(f.o.gz)
  }
  #overlap & index
  cmd = paste("intersectBed -a ", f.sorted.test, " -b ", f.sorted.ref, " -wa -wb -sorted > ", f.o, sep="")
  system(cmd)
  
  cmd = paste("bgzip ", f.o, " && tabix -p bed ", f.o.gz, sep="")
  system(cmd)
  
  
  #counts, processing by chromosome
  reg.all2counts = rep(0, length(ref.regions))
  
  names(reg.all2counts) = ref.regions
  
  for(chr in c(paste("chr", c(1:22,"X"), sep="")))
  {
    #print(chr)
    cmd = paste("tabix ", f.o.gz, " ", chr, " |awk '{test_len = $3-$2; ovlp_min=$2>$", test.column.no+2, "?$2:$", test.column.no+2, "; ovlp_max=$3>$", test.column.no+3, "?$", test.column.no+3, ":$3; test_ovlp_frac=(ovlp_max-ovlp_min)/test_len;print test_ovlp_frac \"\\t\" $", test.column.no+1, "\":\" $", test.column.no+2, "\"-\" $", test.column.no+3, " }' ", sep="")#, f.o
    buf = system(cmd, intern = T)
    buf = t(sapply(buf, FUN=function(x){unlist(strsplit(x, split="\t"))}))
    df = data.frame(frac= as.numeric(buf[,1]), ref.region= buf[,2])
    buf= aggregate(df$frac, by=list(df$ref.region), FUN=sum)
    reg2counts = buf[,2]
    names(reg2counts) = buf[,1]  
   
    reg.all2counts[names(reg2counts)] = reg2counts 
  }
  return(reg.all2counts)
}


#GC normalized as in https://www.nature.com/nature/journal/v464/n7289/full/nature08872.html
#https://static-content.springer.com/esm/art%3A10.1038%2Fnature08872/MediaObjects/41586_2010_BFnature08872_MOESM66_ESM.pdf
getGCnormalizedPeakHeight = function(peak2GC, rawCounts, GC.binNum)
{
  peak2GC.sorted = sort(peak2GC)
  bin.peakNum = round(length(peak2GC.sorted)/GC.binNum)
  bin2pks = list()
  for(i in 1: GC.binNum)
  {
    if(i!=GC.binNum)
    {
      bin2pks[[i]]=names(peak2GC.sorted)[(1+(i-1)*bin.peakNum): (i*bin.peakNum)]
    }else
    {
      bin2pks[[i]]=names(peak2GC.sorted)[(1+(i-1)*bin.peakNum): length(peak2GC.sorted)]
    }
  
  }
  bin2GC= sapply(bin2pks, FUN=function(pks){mean(peak2GC[pks])})
  
  bin2reads = sapply(1:ncol(rawCounts),
  FUN=function(i)
  {
    print(i)
    b2reads = sapply(bin2pks,
    FUN=function(pks)
    {
      sum(rawCounts[pks,i])
    })
  })
  colnames(bin2reads)=colnames(rawCounts)
  samp2readsSum = apply(bin2reads,2,sum)
  samp2readsProp = samp2readsSum/sum(samp2readsSum)
  
  bin2log2RelEnrich = t(apply(bin2reads,1,
  FUN=function(x)
  {
    log2(x/sum(x)/samp2readsProp)
  }))
  
  samp2spline = apply(bin2log2RelEnrich, 2,
  FUN=function(relEnrich)
  {
    smooth.spline(bin2GC, relEnrich, spar=1)
  })
  
  GCNorm = sapply(1:ncol(rawCounts),
  FUN=function(i)
  {
    pk2counts = rawCounts[,i]
    pred.relEn = predict(samp2spline[[i]], peak2GC[names(pk2counts)])$y
    pk2counts*2^(-pred.relEn)
  })
  colnames(GCNorm)=colnames(rawCounts)
  rownames(GCNorm)=rownames(rawCounts)
  return(GCNorm)
}  


###############################################################################
#adapted from mingxiang Teng's script from github
#https://genome.cshlp.org/content/27/11/1930
#a. for GC, remove those positions with a "N"
#b. for GC, allow those regions larger than median and give them equal weights for GC, gctype="uniform"
# 
# @title Adjust ChIP-seq Read Count Table
# 
# @description
# For a given set of sites with the same/comparable width, their
# read count table from multiple samples are adjusted based on
# potential GC effects. For each sample separately, GC effects are
# estimated based on their effective GC content and
# reads count using generalized linear mixture models. Then, count
# table is adjusted based on estimated GC effects.
# It it important that the given sites includes both foreground and
# background regions, see \code{sites} below.
# 
# @param counts A count matrix with each row corresponding to each element
# in \code{sites} and each column corresponding to one sample. Every value
# in the matrix indicates the read counts for one site in one sample. It is
# noted that since effective GC content is used in this function, it is
# important to extend either original reads or original \code{sites} to
# consider reads that 5' starting in \code{flank} regions, when counting
# sequencing reads.
# 
# @param sites A GRanges object with length equivalent to number of rows
# in \code{counts} matrix. It is preferable that every GRange have the same
# width; otherwise, the mixture model is modeling different things with
# wider GRanges certainly have more reads. However, it is OK if only a
# minority of GRanges have different width, since the model is pretty robust
# to outliers. Also, it is important that \code{sites} including both
# foreground and background regions in each sample, otherwise the mixture
# model will fail to fit two components. Fortunately, if you are inputing
# a large collection of samples, foreground sites in one sample may play
# the role as background in other samples. In this case, manually selecting
# real background is not necessary.
# 
# @param flank A non-negative integer specifying the flanking width of
# ChIP-seq binding. This parameter provides the flexibility that reads
# appear in flankings by decreased probabilities as increased distance
# from binding region. This paramter helps to define effective GC
# content calculation.
# 
# @param outputidx A logical vector with the length equivalent to number
# of rows in \code{counts}. This provides which subset of adjusted count
# matrix should be outputed. This would be extremely useful if you have
# manually collected background sites and want to only export the sites
# you care about.
# 
# @param gcrange A non-negative numeric vector with length 2. This vector
# sets the range of GC content to filter regions for GC effect estimation.
# For human, most regions have GC content between 0.3 and 0.8, which is
# set as the default. Other regions with GC content beyond this range
# will be ignored. This range is critical when very few foreground regions
# are selected for mixture model fitting, since outliers could drive the
# regression lines. Thus, if possible, first make a scatter plot between
# counts and GC content to decide  this parameter. Alternatively,
# select a narrower range, e.g. c(0.35,0.7), to aviod outlier effects from
# both high and low GC-content regions.
# 
# @param emtrace A logical vector which, when TRUE (default), allows to
# print the trace of log likelihood changes in EM iterations.
# 
# @param plot A logical vector which, when TRUE (default), returns miture
# fitting plot.
# 
# @param model A character specifying the distribution model to be used in
# generalized linear model fitting. The default is negative
# binomial(\code{nbinom}), while \code{poisson} is also supported currently.
# More details see \code{gcEffects}.
# 
# @param mu0 A non-negative numeric initiating read count signals for
# background sites. This is treated as the starting value of background mean
# for poisson/nbinom fitting.
# 
# @param mu1 A non-negative numeric initiating read count signals for
# foreground sites. This is treated as the starting value of foreground mean
# for poisson/nbinom fitting.
# 
# @param theta0 A non-negative numeric initiating the shape parameter of
# negative binomial model for background sites. For more detail, see
# theta in \code{\link[MASS]{glm.nb}} function.
# 
# @param theta1 A non-negative numeric initiating the shape parameter of
# negative binomial model for foreground sites. For more detail, see
# theta in \code{\link[MASS]{glm.nb}} function.
# 
# @param p A non-negative numeric specifying the proportion of foreground
# sites in all estimated sites. This is treated as a starting value for
# EM algorithm.
# 
# @param converge A non-negative numeric specifying the condition of EM
# algorithm termination. EM algorithm stops when the ratio of log likelihood
# increment to whole log likelihood is less or equivalent to
# \code{converge}.
# 
# @param genome A \link[BSgenome]{BSgenome} object containing the sequences
# of the reference genome that was used to align the reads, or the name of
# this reference genome specified in a way that is accepted by the
# \code{\link[BSgenome]{getBSgenome}} function defined in the \pkg{BSgenome}
# software package. In that case the corresponding BSgenome data package
# needs to be already installed (see \code{?\link[BSgenome]{getBSgenome}} in
# the \pkg{BSgenome} package for the details).
# 
# @param gctype A character vector specifying choice of method to calculate
# effective GC content. Default \code{ladder} is based on uniformed fragment
# distribution. A more smoother method based on tricube assumption is also
# allowed. However, tricube should be not used if \code{flank} is too large.
# 
# @return The count matrix after GC adjustment. The matrix values are not
# integer any more.
# 
# @import S4Vectors
# @import IRanges
# @import GenomicRanges
# @import Biostrings
# @import MASS
# @importFrom BSgenome getBSgenome
# @importFrom BSgenome getSeq
# @importFrom matrixStats colMedians
# @importFrom splines ns
# @importFrom grDevices rgb
# @importFrom graphics plot
# @importFrom graphics axis
# @importFrom graphics lines
# @importFrom stats dpois
# @importFrom stats dnbinom
# @importFrom stats glm
# @importFrom stats predict
# 
# @export

refineSites <- function(counts,sites,flank=250L,
                        outputidx=rep(TRUE,nrow(counts)),
                        gcrange=c(0.3,0.8),emtrace=TRUE,plot=TRUE,
                        model=c('nbinom','poisson'),
                        mu0=1,mu1=50,theta0=mu0,theta1=mu1,
                        p=0.2,converge=1e-4,
                        genome="hg19",gctype=c("ladder","tricube", "uniform")){
   
    ### input sanity check
    genome <- getBSgenome(genome)
    model <- match.arg(model)
    gctype <- match.arg(gctype)
    sitew <- median(width(sites))
    if(sum(gcrange<0)>0 || sum(gcrange>1)>0 || sum(is.na(gcrange))>0)
        stop("Parameter 'gcrange' error.\n")
    if(mu0<=0 || mu1<=0 || mu0>=mu1)
        stop("Parameter 'mu0' or 'mu1' error.\n")
    if(model=='nbinom' && (theta1<=0 || theta0<=0))
        stop("Parameter 'theta0' or 'theta1' error in nbinom model.\n")
    if(p<=0 || p>=1)
        stop("'p' must be in (0,1).\n")
    if(converge<=0 || converge>0.1)
        stop("'converge' must be in (0,0.1].\n")
    if(gctype=="ladder"){
        weight <- c(seq_len(flank),rep(flank+1,sitew),rev(seq_len(flank)))
        weight <- weight/sum(weight)
    }else if(gctype=="tricube"){
        w <- flank+floor(sitew/2)
        weight <- (1-abs(seq(-w,w)/w)^3)^3
        weight <- weight/sum(weight)
    }
    ### effective gc content
    cat("Start to estimate GC effects.\n")
    cat("...... Calculating GC content with flanking",flank,"\n")
    nr <- IRanges::shift(resize(sites,width(sites) + flank*2),-flank)
    seqs <- getSeq(genome,nr)
    if(gctype=="ladder"|gctype=="tricube")
    {
      gcpos <- startIndex(vmatchPattern("S", seqs, fixed="subject"))
      gc <- round(sapply(gcpos,function(x) sum(weight[x])),3)
    }else if(gctype=="uniform") #allow for uniformly weights, updated by Lei 3.30.2018
    {
      gcpos = startIndex(vmatchPattern("S", seqs, fixed="subject"))
      Npos = startIndex(vmatchPattern("N", seqs, fixed=T))
      seq.len = sapply(seqs,length)
      gc =round(sapply(gcpos,length)/(seq.len-sapply(Npos,length)),3)
    }
      
    rm(nr,seqs,gcpos)
    ### em algorithms
    cat("...... Estimating GC effects\n")
    fitmu0 <- fitmu1 <- fitz <- matrix(NA,sum(outputidx),ncol(counts))
    for(rep in seq_len(ncol(counts))){
        cat("......... Estimating sample",rep,"\n")
        rc <- counts[,rep]
        idx <- gc>=gcrange[1] & gc<=gcrange[2] & !is.na(rc) & !is.na(gc)
        dat <- data.frame(y=rc[idx],gc=gc[idx])
        theta1E <- theta1
        theta0E <- theta0
        if(model=='poisson'){
            logp1 <- dpois(dat$y, lambda = mu1, log = TRUE)
            logp0 <- dpois(dat$y, lambda = mu0, log = TRUE)
        }else{
            logp1 <- dnbinom(dat$y, size=theta1E, mu=mu1, log = TRUE)
            logp0 <- dnbinom(dat$y, size=theta0E, mu=mu0, log = TRUE)
        }
        z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
        llf <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
        llgap <- llf
        i <- 0
        while(abs(llgap) > (abs(llf) * converge) && i < 100){
            p <- (2+sum(z))/(2*2+length(z))
            dat1 <- dat[z>=0.5,]
            dat0 <- dat[z<0.5,]
            
            if(nrow(dat0)==0 | nrow(dat1)==0) break #if no dat1 or dat0 population left, estimate GC effect as NA
            if(model=='poisson'){
                lmns0 <- glm(y ~ ns(gc, df = 2), data=dat0, family="poisson")
                lmns1 <- glm(y ~ ns(gc, df = 2), data=dat1, family="poisson")
                predY0 <- predict(lmns0,data.frame(gc=dat$gc),type="response")
                predY1 <- predict(lmns1,data.frame(gc=dat$gc),type="response")
                logp1 <- dpois(dat$y, lambda = predY1, log = TRUE)
                logp0 <- dpois(dat$y, lambda = predY0, log = TRUE)
            }else{
                lmns0 <- glm.nb(y ~ ns(gc, df=2),data=dat0,init.theta=theta0E)
                lmns1 <- glm.nb(y ~ ns(gc, df=2),data=dat1,init.theta=theta1E)
                predY0 <- predict(lmns0,data.frame(gc=dat$gc),type="response")
                predY1 <- predict(lmns1,data.frame(gc=dat$gc),type="response")
                theta1E <- lmns1$theta
                theta0E <- lmns0$theta
                logp1 <- dnbinom(dat$y, size=theta1E, mu=predY1, log = TRUE)
                logp0 <- dnbinom(dat$y, size=theta0E, mu=predY0, log = TRUE)
            }
            z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
            if(sum(z>=0.5) < length(gc)*0.0005 | sum(z<0.5) < length(gc)*0.0005) #just like only fit one population
                break;
            lli <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
            llgap <- lli - llf
            llf <- lli
            i <- i + 1
            if(emtrace)
                cat("......... Iteration",i,'\tll',llf,'\tincrement',llgap,'\n')
        }
        if(plot){
            idx0 <- sample.int(nrow(dat),min(50000,nrow(dat)))
            plot(dat$gc[idx0][z[idx0]<=0.5],dat$y[idx0][z[idx0]<=0.5]+0.5,col=rgb(0,0,0,alpha=0.1),
                 xlim=gcrange,ylim=c(min(dat$y[idx0])+0.5, max(dat$y[idx0])),
                 pch=20,main=paste("rep",rep),log='y', cex=0.5, yaxt='n',
                 xlab='Effective GC content',ylab="Read counts")
            points(dat$gc[idx0][z[idx0]>0.5],dat$y[idx0][z[idx0]>0.5]+0.5,col=rgb(1,0,0,alpha=0.2), #updated by Lei 3.30.2018
                 xlim=gcrange,pch=20,log='y', cex=0.5)
            idx00 <- sample.int(nrow(dat),min(1000,nrow(dat)))
            idx00 <- idx00[order(dat$gc[idx00])]
            lines(dat$gc[idx00],predY1[idx00]+0.5,col='red',lwd=3)
            lines(dat$gc[idx00],predY0[idx00]+0.5,col='blue',lwd=3)
            axis(side=2, at=c(0,2^(0:10))+0.5, labels=c(0,2^(0:10)))
        }
        ### gc effects
        fitmu0[idx[outputidx],rep] <- predY0[outputidx[idx]]
        fitmu1[idx[outputidx],rep] <- predY1[outputidx[idx]]
        fitz[idx[outputidx],rep] <- z[outputidx[idx]]
        cat("......... Sample",rep,"finished\n")
    }
    med1 <- fitmu1
    med0 <- fitmu0
    med1[fitz<0.5] <- NA
    med0[fitz<0.5] <- NA #should be original version which put more weight on sites with signal
    #med0[fitz>0.5] <- NA #updated by Lei 4.1.2018
    gce <- log2(t(t(fitmu1)/colMedians(med1,na.rm=TRUE))) * fitz +
           log2(t(t(fitmu0)/colMedians(med0,na.rm=TRUE))) * (1-fitz)
    counts[outputidx,] / 2^gce
}



################################################################################

#map a vector of number to colors from a range
mapNum2col=function(xs, col.a="white", col.b="red", bin.num=200)
{
  cols=colorRampPalette(c(col.a, col.b))(n = bin.num)
  
  x.min=min(xs)
  x.max=max(xs)
  if(x.min==x.max)
  {
    return(rep(col.a, length(xs)))
  }
  x.bin = (x.max-x.min)/bin.num
  x2col=sapply(xs,
  FUN=function(x)
  {
    b=round((x-x.min)/x.bin)
    if(b>bin.num)
    {
      b=bin.num
    }
    return(cols[b])
  })
  return(unlist(x2col))
  
}


