
library(ape)
library(seqinr)
library(RColorBrewer)

######################################################################
###    ANALYSIS BY ALIGNMENT COLUMNS (positions)

nt.nms <- c("A","C","G","T")

###  Nucleotide frequence by position given an alignment
###    and a vector of counts
FreqMat.w <- function(seqs,w)
{ strm <- t(sapply(seqs,function(sq) strsplit(sq,split="")[[1]]))
  res <-  apply(strm,2,function(x) 
             tapply(w,factor(x,levels=nt.nms),sum))
  res[is.na(res)] <- 0
  res
}  

###  Consensus sequence given alignment and frequencies
ConsSeq <- function(seqs,w)
{ ntm <- FreqMat.w(seqs,w)
  get.nt <- function(x)
  { mx <- max(x)
    idx <- which(x==mx)
    if(length(idx)<2) return(idx)
    return( sample(idx,1) )
  }  
  imx <- apply(ntm,2,get.nt)
  paste(nt.nms[imx],collapse="")
}  

###  Matrix of mutations frequencies per position
###  (the most frequent is set to 0 leaving only the mutations)
MutsTbl <- function(seq.tbl)
{ j <- apply(seq.tbl,2,which.max)
  seq.tbl[cbind(j,1:ncol(seq.tbl))] <- 0
  seq.tbl
}

###  Non-WT relative frequency in polymorphic sites
PolyDist <- function(seqs,w)
{ seq.tbl <- FreqMat.w(seqs,w)
  nt <- sum(seq.tbl[,1]) 
  seq.tbl <- MutsTbl(seq.tbl)
  seq.tbl <- seq.tbl[,apply(seq.tbl,2,function(x) sum(x)>0),drop=FALSE]
  colSums(seq.tbl)/nt
}  
  
###  Mutations barplot with a rug of wt nucleotides
subst.barplot <- function(seqs,nr,i0,ymx=NULL)
{
  ###  Compute the nucleotide frequencies and the mutations table
  pos.tbl <- FreqMat.w(seqs,nr)
  mut.tbl <- MutsTbl(pos.tbl)
  m <- ncol(pos.tbl)
  ### Barplot of mutations along with the consensus sequence rug
  prop.mut <- mut.tbl/sum(pos.tbl[,1])*100
  if(is.null(ymx))
    ymx <- max(apply(prop.mut[1:4,],2,sum))
  x <- i0-1+(1:nchar(seqs[1]))
  cls <- brewer.pal(4,"Dark2")
  plot(x,prop.mut[1,],type="h",col=cls[1],ylim=c(-ymx/40,ymx),
       lwd=2,xlab="position",ylab="mutations proportion (%)",
       lend=2)
  h0 <- prop.mut[1,]
  segments(x0=x,y0=h0,y1=h0+prop.mut[2,],col=cls[2],lwd=2,lend=2)
  h0 <- h0+prop.mut[2,]
  segments(x0=x,y0=h0,y1=h0+prop.mut[3,],col=cls[3],lwd=2,lend=2)
  h0 <- h0+prop.mut[3,]
  segments(x0=x,y0=h0,y1=h0+prop.mut[4,],col=cls[4],lwd=2,lend=2)
  legend("topleft",col=cls[1:4],legend=rownames(prop.mut)[1:4],
         cex=0.6,lwd=2)
  segments(x0=x+.25,y0=rep(-ymx/25,m),y1=rep(0,m),lwd=2,lend=2,
         col=cls[apply(pos.tbl,2,which.max)])
  abline(h=0)
}

###  Information content at a site
###    given the nucleotides frequency
InfContent <- function(v)
{ v <- v/sum(v)
  lgv <- ifelse(v==0,0,log2(v))
  2+sum(v*lgv)
}

###  Information contents profile of an alignment
get.inf.prof <- function(seqs,nr)
{ fm <- FreqMat.w(seqs,nr)
  apply(fm,2,InfContent)
}

###  Hill numbers profile by position for q=0,1,2,3
Hill.posProf <- function(seqs,nr)
{ nr <- nr/sum(nr)
  fm <- FreqMat.w(seqs,nr)
  ###  Hill number of first order
  foq <- function(v)
  { v <- v/sum(v)
    lgv <- ifelse(v==0,0,log(v))
    exp(-sum(v*lgv))
  }
  res <- matrix(0,nrow=4,ncol=nchar(seqs[1]))
  rownames(res) <- paste("q",0:3,sep="")
  res[1,] <- apply(fm,2,function(x) sum(x>0))
  res[2,] <- apply(fm,2,function(x) foq(x))
  res[3,] <- apply(fm,2,function(x) 1/sum(x^2))
  res[4,] <- apply(fm,2,function(x) 1/sqrt(sum(x^3)))
  res
}

######################################################################
###    ANALYSIS BY ALIGNMENT ROWS (haplotypes)

###  Shannon entropy 
###    w  vector of observed counts or frequencies
Shannon.entropy <- function(w)
{ h <- length(w)
  if(h<2) return(0)
  p <- w/sum(w)
  lgp <- ifelse(w==0,0,log(p))
  S <- -sum(p*lgp)
  ###  Correct bias only if vector of counts given
  if(all(w>=1)) S <- S-(h-1)/(2*sum(w)) ##Was S+(h-1)/(2*sum(w))??????
  if(S>log(h)) S <- log(h)  
  return(S)
}  

###  Shannon entropy asymptotic variance
###    given a vector of counts
Shannon.entropy.var <- function(w) 
{ h <- length(w) 
  if(h<2) return(0)
  N <- sum(w)
  if(N<2) return(NULL)
  w <- w/N
  lgw <- ifelse(w==0,0,log(w))
  S <- -sum(w*lgw)
  ((sum(w*lgw^2)-S^2) + (h-1)/(2*N)) / N
}  

###  Normalized Shannon entropy given a vector of counts
###    w  vector of observed counts or frequencies
N.Shannon.entropy <- function(w) 
{ h <- length(w)
  if(h<2) return(0)
  S <- Shannon.entropy(w)
  S/log(h)
}  

###  Normalized Shannon entropy asymptotic variance
###    given a vector of counts
N.Shannon.entropy.var <- function(w) 
{ h <- length(w) 
  if(h<2) return(0)
  N <- sum(w)
  if(N<2) return(NULL)
  w <- w/N
  lgw <- ifelse(w==0,0,log(w))
  S <- -sum(w*lgw)
  ((sum(w*lgw^2)-S^2+(h-1)/(2*N)) / (log(h)^2*N) ) 
}  


#finge-trimmed corrected Shannon entropy
correct.Shannon.entropy <- function(w)
{ h <- length(w)
if(h<2) return(0)
ifelse()

p <- w/sum(w)
lgp <- ifelse(w==0,0,log(p))
S <- -sum(p*lgp)
###  Correct bias only if vector of counts given
if(all(w>=1)) S <- S-(h-1)/(2*sum(w)) ##Was S+(h-1)/(2*sum(w))??????
if(S>log(h)) S <- log(h)  
return(S)
} 

#########################################################
###  GINI-SIMPSON INDEX

###  Gini-Simpson unbiased estimator 
###    w  vector of observed counts
GiniSimpson <- function(w)
{ n <- sum(w)
  if(n<2) return(NULL)
  p <- w/n
  (1 - sum(p^2))*n/(n-1)
}

###  Gini-Simpson asymptotic variance (Nayak 1983,1985) 
###    w  vector of observed counts
GiniSimpson.var <- function(w)
{ p <- w/sum(w)
 4/n*(sum(p^3)-sum(p^2)^2)
}

###  MVUE Gini-Simpson 
###    w  vector of observed counts or frequencies
GiniSimpson.MVUE <- function(w)
{ p <- w/sum(w)
  pm1 <- (w-1)/(sum(w)-1)
  1 - sum(p*pm1)
}

#########################################################
###  HIGHER ORDER HILL NUMBERS

###  MVUE third order Hill number
###    w  vector of observed counts or frequencies
D3 <- function(w)
{ p <- w/sum(w)
  pm1 <- (w-1)/(sum(w)-1)
  pm2 <- (w-2)/(sum(w)-2)
  pm2[pm2<0] <- 0
  1/sqrt(sum(p*pm1*pm2))
}

#########################################################
###  HILL NUMBERS and profile

###    w  vector of observed counts or frequencies
###    q  exponent
qD <- function(w,q)
{ if(q==0) return(length(w))
  if(q==1) return( exp(Shannon.entropy(w)) )
  p <- w/sum(w)
  if(q==Inf) return(1/max(p))
  if(q==-Inf) return(1/min(p))
  sum(p^q)^(1/(1-q))
}

###    w  vector of observed counts or frequencies
qD.profile <- function(w,q=NULL)
{ if(is.null(q))
    q <- c(seq(0,0.9,0.1),seq(1,1.8,0.2),seq(2,3.75,0.25),
           seq(4,10,1),Inf)
  dv <- numeric(length(q))
  for(i in 1:length(q))
    dv[i] <- qD(w,q[i])
  data.frame(q=q,qD=dv)
}  

#########################################################
###  Rényi entropies and profile

###    w  vector of observed counts or frequencies
###    q  exponent
renyi <- function(w,q)
{ if(q==0) return(log(length(w)))
  if(q==1) return(Shannon.entropy(w))
  p <- w/sum(w)
  if(q==Inf) return( -log(max(p)) )
  log(sum(p^q))/(1-q)
}

###    w  vector of observed counts or frequencies
renyi.profile <- function(w,q=NULL)
{ if(is.null(q))
    q <- c(seq(0,0.9,0.1),seq(1,1.8,0.2),seq(2,3.75,0.25),
           seq(4,10,1),Inf)
  dv <- numeric(length(q))
  for(i in 1:length(q))
    dv[i] <- renyi(w,q[i])
  data.frame(q=q,renyi=dv)
}  

#########################################################
###  Havrda-Charvat or Tsallis indices and profile

###  Havrda-Charvat estimator
###    w  vector of observed counts or frequencies
###    q  exponent
HCq <- function(w,q)
{ if(q==0) return(length(w)-1)
  if(q==1) return(Shannon.entropy(w))
  if(q==Inf) return(0)
  p <- w/sum(w)
  (1-sum(p^q))/(q-1)
}

###  Havrda-Charvat  asymptotic variance (Nayak 1983,1985)
###    w  vector of observed counts
###    q  exponent
HCq.var <- function(w,q)
{ n <- sum(w)
  if(n<2) return(NULL)
  p <- w/n
  1/n*(q/(q-1))^2*(sum(p^(2*q-1))-sum(p^q)^2)
}

###    w  vector of observed counts or frequencies
HCq.profile <- function(w,q=NULL)
{ if(is.null(q))
    q <- c(seq(0,0.9,0.1),seq(1,1.8,0.2),seq(2,3.75,0.25),
           seq(4,10,1),Inf)
  dv <- numeric(length(q))
  for(i in 1:length(q))
    dv[i] <- HCq(w,q[i])
  data.frame(q=q,HC=dv)
}  

#########################################################
###  POLYMORPHIC SITES & MUTATION FREQUENCY

###  Polymorphic sites and number of different mutations
poly.sites <- function(seqs)
{ nsq <- length(seqs)
  frqmat <- consensusMatrix(DNAStringSet(seqs))[c("A","C","G","T"),]
  psites <- apply(frqmat,2,function(nts) sum(nts>0))
  poly.pos <- which(psites>1)
  psites <- psites[poly.pos]
  tbl <- table(psites)  
  nmuts <- sum(psites-1)
  psites <- length(psites)
  list(psites=psites,nmuts=as.integer(nmuts),poly.pos=poly.pos)
}  

###  Aa polymorphic sites and number of different mutations
Aa.poly.sites <- function(seqs)
{ nsq <- length(seqs)
  frqmat <- consensusMatrix(AAStringSet(seqs))
  psites <- apply(frqmat,2,function(nts) sum(nts>0))
  poly.pos <- which(psites>1)
  psites <- psites[poly.pos]
  tbl <- table(psites)  
  nmuts <- sum(psites-1)
  psites <- length(psites)
  list(psites=psites,nmuts=as.integer(nmuts),poly.pos=poly.pos)
}  

###  Mutation frequency
MutationFreq <- function(w,nm,len)
{ mf <- sum(nm*w/sum(w))/len
  names(mf) <- NULL
  mf
}

###  Mutation frequency standard deviation
###    given a vector of counts (cf. M. Salicru)
MutationFreq.var <- function(w,nm,len)
{ 
  N <- sum(w)
  if(N<2) return(0)
  p <- w/sum(w)
  v <- ((sum(p*nm^2)-sum(p*nm)^2)/len^2) / N
  names(v) <- NULL
  v
} 

#########################################################
###  GENETIC DISTANCE

SortByMutations <- function(bseqs,nr)
{
  master <- bseqs[which.max(nr)]
  ##  Determine differences regarding the consensus string
  bseqs <- DNAStringSet(bseqs)
  psa <- pairwiseAlignment(pattern=bseqs,subject=master)
  nm <- nmismatch(psa)
  tnm <- table(nm)
  ##  Sort by number of mutations
  o <- order(nm)
  bseqs <- bseqs[o]
  nr <- nr[o]
  nm <- nm[o]
  ##  Order number within each number of mutations
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  ##  Sort by descending frequency within each number of mutations
  for(i in as.integer(names(tnm)))
  { idx <- which(nm==i)
    o <- order(nr[idx],decreasing=TRUE)
    bseqs[idx] <- bseqs[idx[o]]
    nr[idx] <- nr[idx[o]]
  }
  ##  Calculate relative frequency
  frq <- round(nr/sum(nr)*100,2)
  ## Full name for each haplotype
  nms <- paste("Hpl",nm,sprintf("%04d",isq),sep="_")
  ##  Header with name, number of reads, and relative frequency
  names(bseqs) <- nms
  list(bseqs=as.character(bseqs),nr=nr,nm=nm)
}

DNA.dist <- function(seqs,model="raw",gamma=FALSE)
{
  strm <- t(sapply(seqs,function(sq) strsplit(sq,split="")[[1]]))
  dst <- dist.dna(as.DNAbin(ape::as.alignment(strm)),model=model,gamma=gamma)
  dst[is.na(dst)] <- 0
  dst
}

#########################################################
###  RAO's QUADRATIC ENTROPY ==  NUCLEOTIDIC DIVERSITY

###  Rao's entropy  == Nucleotidic diversity
###   w vector of counts
###   dst a distance object or an aquare matrix of distances
###  Computed as eq 10.5 in Nei (1987) pg.256
rao <- function(w,dst)
{ n <- sum(w)
  if(n<2) return(0)
  p <- w/n
  D <- as.matrix(dst)
  #  El factor n/(n-1) corregeix de dividir per n^2 a
  #  dividir pel nombre de parells de seqüències n*(n-1)
  (n/(n-1)) * (t(p) %*% D %*% p)
}

###  Variance of the Rao estimator
###    (cf. S. Pavoine thesis)
###   w vector of counts
###   dst a distance object or an aquare matrix of distances
rao.var <- function(w,dst)
{ n <- sum(w)
  if(n<2) return(0)
  p <- w/n
  D <- as.matrix(dst)
  S <- -(p%*%t(p))
  diag(S) <- p*(1-p)
  4*t(p)%*%D%*%S%*%D%*%p/n
}

#########################################################
###  Rao entropy of order q
 
###  Rao entropy of order q
###   dst a distance object or an square matrix of 
###     distances or disimilarities.
###  Vectorized in q.
rao.pow <- function(w,dst,q)
{ n <- sum(w)
  if(n<2) return(NULL)
  p <- w/n
  D <- as.matrix(dst)
  Q <- (t(p) %*% D %*% p)    #  <p|D|p>
  O <- p %*% t(p)            #  |p><p|
  m <- length(q)
  res <- sapply(1:m,function(i) { sum(D*O^q[i]) })
  list(qQ=res,Q=as.numeric(Q))
}

###  Functional Hill numbers profile
###   w  vector of observed counts or frequencies
###   dst a distance object or an square matrix of 
###     distances or disimilarities.
rao.pow.profile <- function(w,dst,q=NULL)
{ if(is.null(q))
    q <- seq(0,2,0.1)
  dv <- rao.pow(w,dst,q)
  list(qQP=data.frame(q=q,qQ=dv$qQ),Q=dv$Q)
}  

#########################################################
###  FUNCTIONAL DIVERSITY

###  Functional Hill number of order q
###   See Chiu & Chao (2014)
###   dst a distance object or an square matrix of 
###     distances or disimilarities.
###  w is a vector of counts.
###  Vectorized in q.
func.qD <- function(w,dst,q)
{ n <- sum(w)
  if(n<2) return(NULL)
  p <- w/n
  D <- as.matrix(dst)
  Q <- (t(p) %*% D %*% p)  #  <p|D|p>
  Q <- Q * (n/(n-1))       #  correcting for pairs of individuals
  O <- p %*% t(p)          #  |p><p|
  m <- length(q)
  res <- numeric(m)
  for(i in 1:m)
  { if(q[i]==0)
    { res[i] <- sqrt(sum(D)/Q)
	  next
	}
	if(q[i]==1)
    { res[i] <- exp(-sum(D*O*log(O))/(2*Q))
	  next
	}
    res[i] <- (sum(D*O^q[i])/Q)^(1/(2*(1-q[i])))
  }	  
  list(fqD=res,Q=as.numeric(Q))
}

###  Functional Hill numbers profile
###   w  vector of observed counts or frequencies
###   dst a distance object or an square matrix of 
###     distances or disimilarities.
func.qD.profile <- function(w,dst,q=NULL)
{ if(is.null(q))
    q <- c(seq(0,1.8,0.2),seq(2,3.75,0.25),seq(4,10,1))
  dv <- func.qD(w,dst,q)
  list(fqD=data.frame(q=q,qD=dv$fqD),Q=dv$Q)
}  

###  Functional Hill number of order q given Q
###   See Chiu & Chao (2014)
###   dst a distance object or an square matrix of 
###     distances or disimilarities.
###  w is a vector of counts.
###  Vectorized in q.
func.u.qD <- function(w,dst,Q,q)
{ n <- sum(w)
  if(n<2) return(NULL)
  p <- w/n
  D <- as.matrix(dst)
  O <- p %*% t(p)          #  |p><p|
  m <- length(q)
  res <- numeric(m)
  for(i in 1:m)
  { if(q[i]==0)
    { res[i] <- sqrt(sum(D)/Q)
	  next
	}
	if(q[i]==1)
    { res[i] <- exp(-sum(D*O*log(O))/(2*Q))
	  next
	}
    res[i] <- (sum(D*O^q[i])/Q)^(1/(2*(1-q[i])))
  }	  
  res
}

###  Functional Hill numbers profile given Q
###   w  vector of observed counts or frequencies
###   dst a distance object or an square matrix of 
###     distances or disimilarities.
func.u.qD.profile <- function(w,dst,Q,q=NULL)
{ if(is.null(q))
    q <- c(seq(0,1.8,0.2),seq(2,3.75,0.25),seq(4,10,1))
  dv <- func.u.qD(w,dst,q) #??????? needs to include Q?
  data.frame(q=q,qD=dv)
}  


#########################################################
###  FULL SET OF INDICES

###  Computes a set of diversity indices for a given population
###   on nt sequences
###   seqs:  vector of strings (sequences)
###   nr:    counts of each sequence in seqs (same order)
###   frame: -1 no coding, 1,2,3 index of first nt in full codon
sq.diversity2 <- function(seqs,nr,frame=-1)
{
  inc.div <- integer(5)
  names(inc.div) <- c("Len","Counts","Hpl","PolySites","nMuts")
  ab.div <- numeric(4)
  names(ab.div) <- c("Mstr","Mpct","Shannon","GiniS")
  m <- ifelse(frame %in% 1:3,7,4)
  gnms <- c("raw","K80","Trans","Tranv","Synon","NonSyn","dN2dS")
  Mf <- numeric(m)
  names(Mf) <- gnms[1:m]
  Pi <- numeric(m)
  names(Pi) <- gnms[1:m]
  Mf.e <- numeric(m)
  names(Mf.e) <- gnms[1:m]
  Pi.e <- numeric(m)
  names(Pi.e) <- gnms[1:m]
  Ratio <-  numeric(m)
  names(Ratio) <- gnms[1:m]
  FAD <- numeric(m)
  names(FAD) <- gnms[1:m]

  o <- order(nr,decreasing=TRUE)
  nr <- nr[o]
  seqs <- seqs[o]
  len <- nchar(seqs[1])

  inc.div["Len"] <- nchar(seqs[1])
  inc.div["Counts"] <- sum(nr)
  inc.div["Hpl"] <- length(nr)

  ab.div["Mstr"] <- nr[1]
  ab.div["Mpct"] <- round(nr[1]/sum(nr)*100,2)

  ###  If there is only one haplotype
  if( length(nr) < 2 )  
   return(list(inc.div=inc.div,ab.div=ab.div,FAD=FAD,Mf=Mf,Pi=Pi,Ratio=Ratio))
   
  lst <- poly.sites(seqs)  
  inc.div["PolySites"] <- lst$psites
  inc.div["nMuts"] <- lst$nmuts

  ab.div["Shannon"] <- Shannon.entropy(nr)
  ab.div["GiniS"] <- GiniSimpson(nr)
  
  nru <- rep(1,length(nr))
  strm <- t(sapply(seqs,function(sq) strsplit(sq,split="")[[1]]))
  dna.bin <- as.DNAbin(ape::as.alignment(strm))
  dst <- dist.dna(dna.bin,model="N")
  dst[is.na(dst)] <- 0
  dst <- dst/len
  FAD["raw"] <- sum(as.matrix(dst))
  Mf["raw"] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
  Pi["raw"] <- rao(nr,dst)
  Ratio["raw"] <- Pi["raw"]/Mf["raw"]
  Mf.e["raw"] <- sum(as.matrix(dst)[1,]*nru)/sum(nru)
  Pi.e["raw"] <- rao(nru,dst)

  dst <- dist.dna(dna.bin,model="K80")
  dst[is.na(dst)] <- 0
  FAD["K80"] <- sum(as.matrix(dst))
  Mf["K80"] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
  Pi["K80"] <- rao(nr,dst)
  Ratio["K80"] <- Pi["K80"]/Mf["K80"]
  Mf.e["K80"] <- sum(as.matrix(dst)[1,]*nru)/sum(nru)
  Pi.e["K80"] <- rao(nru,dst)

  dst <- dist.dna(dna.bin,model="TS")
  dst[is.na(dst)] <- 0
  dst <- dst/len
  FAD["Trans"] <- sum(as.matrix(dst))
  Mf["Trans"] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
  Pi["Trans"] <- rao(nr,dst)
  Ratio["Trans"] <- Pi["Trans"]/Mf["Trans"]
  Mf.e["Trans"] <- sum(as.matrix(dst)[1,]*nru)/sum(nru)
  Pi.e["Trans"] <- rao(nru,dst)

  dst <- dist.dna(dna.bin,model="TV")
  dst[is.na(dst)] <- 0
  dst <- dst/len
  FAD["Tranv"] <- sum(as.matrix(dst))
  Mf["Tranv"] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
  Pi["Tranv"] <- rao(nr,dst)
  Ratio["Tranv"] <- Pi["Tranv"]/Mf["Tranv"]
  Mf.e["Tranv"] <- sum(as.matrix(dst)[1,]*nru)/sum(nru)
  Pi.e["Tranv"] <- rao(nru,dst)

  if( frame %in% 1:3 )
  { ll <- len-frame+1
    ll <- ll-(ll%%3)+(frame-1)
	seqs <- substr(seqs,frame,ll)
	algn <- list(nb=length(seqs),nam=names(seqs),seq=seqs,
	             com=rep("",length(seqs)))
	attr(algn,"class") <- "alignment"
	lst <- kaks(algn)
	dst <- lst$ks
	dst[is.na(dst)] <- 0
    FAD["Synon"] <- sum(as.matrix(dst))
    Mf["Synon"] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
    Pi["Synon"] <- rao(nr,dst)
    Ratio["Synon"] <- Pi["Synon"]/Mf["Synon"]
    Mf.e["Synon"] <- sum(as.matrix(dst)[1,]*nru)/sum(nru)
    Pi.e["Synon"] <- rao(nru,dst)
	dst <- lst$ka
	dst[is.na(dst)] <- 0
    FAD["NonSyn"] <- sum(as.matrix(dst))
    Mf["NonSyn"] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
    Pi["NonSyn"] <- rao(nr,dst)
    Ratio["NonSyn"] <- Pi["NonSyn"]/Mf["NonSyn"]
    Mf.e["NonSyn"] <- sum(as.matrix(dst)[1,]*nru)/sum(nru)
    Pi.e["NonSyn"] <- rao(nru,dst)
    FAD["dN2dS"] <- FAD["NonSyn"]/FAD["Synon"]
    Mf["dN2dS"] <- Mf["NonSyn"]/Mf["Synon"]
    Pi["dN2dS"] <- Pi["NonSyn"]/Pi["Synon"]
    Mf.e["dN2dS"] <- Mf.e["NonSyn"]/Mf.e["Synon"]
    Pi.e["dN2dS"] <- Pi.e["NonSyn"]/Pi.e["Synon"]
    Ratio["dN2dS"] <- Pi["dN2dS"]/Mf["dN2dS"]
  }
  list(inc.div=inc.div,ab.div=ab.div,FAD=FAD,Mf=Mf,Pi=Pi,Ratio=Ratio,
       Mf.e=Mf.e,Pi.e=Pi.e)
}

###  Computes a set of diversity indices for a given population
###   on nt sequences, given an evolutionary model.
sq.diversity <- function(seqs,nr,nm,model="raw",gamma=FALSE)
{
  res <- data.frame(Counts=integer(1),Hpl=integer(1),PolySites=integer(1),
                    nMuts=integer(1),Mstr=numeric(1),
					Shannon=numeric(1),GiniS=numeric(1),
					Mf=numeric(1),Rao=numeric(1),Ratio=numeric(1))
  res$Counts[1] <- sum(nr)
  res$Hpl[1] <- length(nr)
  lst <- poly.sites(seqs)  
  res$PolySites[1] <- lst$psites
  res$nMuts[1] <- lst$nmuts
  res$Mstr[1] <- round(nr[1]/sum(nr)*100,1)
  res$Shannon[1] <- Shannon.entropy(nr)
  res$GiniS[1] <- GiniSimpson(nr)
  strm <- t(sapply(seqs,function(sq) strsplit(sq,split="")[[1]]))
  dst <- dist.dna(as.DNAbin(ape::as.alignment(strm)),model=model,gamma=gamma)
  dst[is.na(dst)] <- 0
  res$Mf[1] <- sum(as.matrix(dst)[1,]*nr)/sum(nr)
  res$Rao[1] <- rao(nr,dst)
  res$Ratio[1] <- res$Rao[1]/res$Mf[1]
  res
}

###  Computes a set of diversity indices for a given population
###   on Aa sequences
Aa.sq.diversity <- function(seqs,nr,nm)
{
  res <- data.frame(Counts=integer(1),Hpl=integer(1),PolySites=integer(1),
                    nMuts=integer(1),Mstr=numeric(1),
					Shannon=numeric(1),GiniS=numeric(1),
					Mf=numeric(1),Rao=numeric(1),Ratio=numeric(1))
  res$Counts[1] <- sum(nr)
  res$Hpl[1] <- length(nr)
  lst <- Aa.poly.sites(seqs)  
  res$PolySites[1] <- lst$psites
  res$nMuts[1] <- lst$nmuts
  res$Mstr[1] <- round(nr[1]/sum(nr)*100,1)
  res$Shannon[1] <- Shannon.entropy(nr)
  res$GiniS[1] <- GiniSimpson(nr)
  res$Mf[1] <- MutationFreq(nr,nm,nchar(seqs[1]))
  dst <- stringDist(seqs)/nchar(seqs[1])
  res$Rao[1] <- rao(nr,dst)
  res$Ratio[1] <- res$Rao[1]/res$Mf[1]
  res
}

#########################################################
###  BOOTSTRAP

### Sample n individuals from a population with counts
###  by class given by w.
###   indexes of haplotypes
boot.cloning <- function(all.ind,rnd,cum.dist,div.fn,...)
{ dots <- list(...)  
  idx <- sapply(rnd,function(x) which.min(abs(cum.dist-x))[1])
  idx <- ifelse(rnd>cum.dist[idx],idx+1,idx)
  tbl <- table(idx)
  ihpl <- as.integer(names(tbl))
  w <- as.vector(tbl)
  names(w) <- names(cum.dist)[ihpl]
  if(!is.null(dots$nm)) dots$nm <- dots$nm[ihpl]
  if(!is.null(dots$dst)) dots$dst <- dots$dst[ihpl,ihpl]
  largs <- c(list(w=w),dots)
  do.call(div.fn,args=largs)
}

###  SIMPLE EXAMPLE:
###
### frq <- c(200,50,80,120,10,5)
### vals <- runif(6,20,50)
### names(frq) <- names(vals) <- paste("V",1:length(frq),sep="")
### (mean.val <- sum(frq*vals)/sum(frq))
### sum(frq)
### mn.fun <- function(data,w,vals)
### { #print(w); print(sum(w)); 
###   sum(w*vals[names(w)])/sum(w) 
### }
### boot(1:sum(frq),boot.cloning,R=999,cum.dist=cumsum(frq),
###      div.fn=mn.fun,vals=vals) 
###
###  CURRENT USAGE:
###
### library(boot)
### boot(1:sum(w),boot.cloning,R=999,cum.dist=cumsum(w),
###      div.fn=Shannon.entropy)
### boot(1:sum(w),boot.cloning,R=999,cum.dist=cumsum(w),
###      div.fn=MutationFreq,nm=nm,len=len)
### boot(1:sum(w),boot.cloning,R=999,cum.dist=cumsum(w),
###      div.fn=rao,dst=dst.matrix)


#Corrected segregating sites S/a1 nucleotides
correct.seg.sites <- function(seqs, count){
  if(length(count)<2) return(0)
  S <- poly.sites(seqs)$psites
  h <- length(count) #number of haplotypes
  a1 <- numeric(h-1)
  for(i in 1: (h-1)){
    a1[i] <- 1/i
  }
  a1 <- sum(a1)
  cor.seg.sites <- S/a1
  cor.seg.sites
}

#Corrected segregating sites S/a1 amino acids
AA.correct.seg.sites <- function(seqs, count){
  if(length(count)<2) return(0)
  S <- Aa.poly.sites(seqs)$psites
  h <- length(count) #number of haplotypes
  a1 <- numeric(h-1)
  for(i in 1:(h-1)){
    a1[i] <- 1/i
  }
  a1 <- sum(a1)
  cor.seg.sites <- S/a1
  cor.seg.sites
}