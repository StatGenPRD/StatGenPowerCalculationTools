# Last Update: 12/23/2008 10:00AM

#-------------------------------------------------------------------------------
# title: Power.R
# $Source: Power.R,v $
# $Revision: 1.0 $
# $Date: 2008/12/23 10:00:00 $
# $Author: Silviu-Alin Bacanu & Judong Shen
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# README for Power.R
# Power.R includes all the functions that can be used to do both power 
#   estimation and sample size estimation.
#
# Input files
# Parameter file contains all parameters needed.
#
# Output file:
# *.pdf: Power plots or sample size estimation plot
# *.txt: a table that contains the power or estimated sample size results.
# *.log: log file
#-------------------------------------------------------------------------------

# Start of functions

#-------------------------------------------------------------------------------
# Function: parcheck
# This function checks the logic of parameter setting
#-------------------------------------------------------------------------------
parcheck <- function() {
  if (sum(is.null(hsqmat),is.null(deltamat))!=1) {
    write("\n(Only) one of hsqmat and deltamat should be 'NULL'.\n", file=logfile,
          append=T)
  }
  if (n < 0 | n == 0) {
    write("\nWarning: n should NOT be set as negative or 0.\n", 
            file=logfile, append=T) 
    stop("Parameter n (sample size) should NOT be set as negative or 0")   
  }
  if (!(method %in% c("A","S","M"))) {
    write("\nWarning: the parameter 'method' should be one of 'A','S' and 'M', please reset it.\n", 
        file=logfile, append=T)
    stop("The parameter 'method' should be one of 'A','S' and 'M'.")      
  } 
  if (!(outcome %in% c("power","sample"))) {
    write("\nWarning: the parameter 'outcome' should be set as one of 'power' and 'sample', please reset it.\n", 
        file=logfile, append=T)
    stop("The parameter 'outcome' should be one of 'power' and 'sample'.")  
  }  
  if (nmin < 5) {
    write("\nWarning: 'nmin' parameter should be set as greater than or equal to 5, please reset it.\n", 
        file=logfile, append=T)
    stop("The parameter 'nmin' is set as too small (<5).")      
  } 
  if (outcome == "sample" & (targetpow <= 0 | targetpow > 1 | 
      !exists("targetpow") | is.na(targetpow) | is.null(targetpow))) {
    write("\nWarning: the parameter 'targetpow' should be set b/w 0 and 1 when outcome is set as 'sample', please reset it.\n", 
          file=logfile, append=T)
    stop("Please reset 'targetpow' parameter.")
  } 
}

#-------------------------------------------------------------------------------
# Function: getmeans
# This function computes penetrances for a causal polymorphism with the above
#  parameters.
#   p: MAF
#   sdmat: the vector containing the sd for the 3 genotypes
#   hsq: heritabilty (R^2)
#   delta: difference in means between the 2 homozygots
#   frac: for use with general and it is the fraction of the diatance
#         between the 2 homozygots where the hetozygot sits
#   modi: mode of inheritance with values "a" (additive), "d" (dominant),
#         "r" (recessive) and "g" is general
# Please note effect size is given by either hsq and delta (note one must be null).
#-------------------------------------------------------------------------------
getmeans <- function(p=p,sdmat=rep(1,3),hsq=NULL,delta=NULL,frac=1/4, modi=modi){
	q <- 1-p
	mu0<-0
	mu1<-0
	mu2<-0
	if(modi == "r"){
		frac=0
	} else if(modi == "d") {
		frac=1
	} else if (modi == "a") {
		frac=1/2
	}
	if(!is.null(hsq)){
	 delta=sqrt(-(hsq*((-1+p)^2*sdmat[[1]]^2+p*(-2*(-1+p)*sdmat[[2]]^2+p*sdmat[[3]]^2)))/
     		((1-hsq)*(-1+p)*p*(2*frac^2+p-4*frac^2*p +(p-2*frac*p)^2)))
	}
	mu1<-frac*delta
	mu2<-delta
	mu<-c(mu0,mu1,mu2)   #   mu: means vector
	mu
}

#-------------------------------------------------------------------------------
# Function: getsampQT
# This function simulates random samples of QT
#   n: sample size
#   p: MAF (for disease SNP)
#   mu: means vector (mu0,mu1,mu2)
#   sdmat: standard deviations vector for sd's of distributions for the signal 
#          at the 3 genotypes
#-------------------------------------------------------------------------------
getsampQT <- function(n=n, p=p, mu=rep(0.0,3),sdmat=rep(1,3)){
	pheno <- NULL
	geno <- NULL
	if(n> 0){
		geno <- sample(c(0,1,2),size=n,prob=c((1-p)^2,2*p*(1-p),p^2),
	          replace=T)
		pheno<-rnorm(n,mean=mu[geno+1],sd=sdmat[geno+1])
	}
	data <- cbind(pheno,geno)
	data
}

#-------------------------------------------------------------------------------
# Function: getpvalssimQT
# This function gets the exact p-values from simulated samples
#   n: sample size
#   p: MAF (for disease SNP)
#   mu: means vector (mu0,mu1,mu2)
#   sdmat: standard deviations vector for sd's of distributions for the signal 
#          at the 3 genotypes
#   nmin: the number of sims at each level in the design
#-------------------------------------------------------------------------------
getpvalssimQT<-function(n=n,p=p,mu=rep(0,3),sdmat=rep(1,3),nmin=5){
	data<-getsampQT(n=n, p=p, mu=mu,sdmat=sdmat)
	mod0<-lm(data[,1]~1)
	if(sum(as.numeric(data[,2]>=1))>=nmin){
		mod1<-lm(data[,1]~data[,2]+as.numeric(data[,2]>=1))
	} else {
		mod1<-lm(data[,1]~data[,2])
	}
	pval<-anova(mod0,mod1)[2,"Pr(>F)"]
	pval
}

#-------------------------------------------------------------------------------
# Function: getnoncQT
# This function gets the asymptotic noncentrality of QT analysis.
#   n: sample size
#   p: MAF (for disease SNP)
#   mu: means vector (mu0,mu1,mu2)
#   sdmat: standard deviations vector for sd's of distributions for the signal 
#          at the 3 genotypes
#   minsd: minimum value of sd
#-------------------------------------------------------------------------------
getnoncQT <- function(n=n,p=p, mu=mu,sdmat=rep(1,3),minsd=10^-5){
	sdmat[sdmat<minsd]<-minsd
	pvec<-c((1-p)^2,2*p*(1-p),p^2)
	grandmu<-sum(pvec*mu)
	nonc<-(sum(pvec*n*(mu-grandmu)^2))/sum(pvec*sdmat^2)
	nonc
}

#-------------------------------------------------------------------------------
# Function: getpowasymptQT
# This function gets the asymptotic power for a QT analysis
#   n: sample size
#   p: MAF (for disease SNP)
#   mu: means vector (mu0,mu1,mu2)
#   sdmat: standard deviations vector for sd's of distributions for the signal 
#          at the 3 genotypes
#   alpha: type I error
#-------------------------------------------------------------------------------
getpowasymptQT <- function(n=n, p=p, mu=mu, sdmat=sdmat, alpha=10^-7){
  r <- 3   # r is the number of groups
  nonc <- getnoncQT(n=n,p=p,mu=mu,sdmat=sdmat)
  pow <- pchisq(qchisq(alpha,df=r-1,low=F),df=r-1,ncp=nonc,low=F)
  c(pow, nonc)
}

#-------------------------------------------------------------------------------
# Function: pownoncchisqeq
# This function calculates the ncp (non-centrality parameter) or lambda to get
#  the targetpow
#   alpha: type I error
#   df: degree of freedom
#   lambda: ncp - non-centrality parameter
#   targetpow: the target power
#-------------------------------------------------------------------------------
pownoncchisqeq <- function(alpha=alpha, df=2, lambda=lambda, targetpow=0.8){
	pchisq(qchisq(alpha,df=df,low=F),df=df,ncp=lambda,low=F)-targetpow
}


#-------------------------------------------------------------------------------
# Function: getpowsampleQT
# This is the master function for power for QT analysis of random samples.
#   n: sample size
#   p: MAF (for disease SNP)
#   hsqmat: hsq (heritabilty or R^2) vector
#   deltamat: delta vector (difference between means of the two homozygotes)
#   sdmat: the vector of sds for distributions for the 3 genotypes
#   modi: mode of inheritance with values "a" (additive), "d" (dominant) and
#         "r" (recessive) and "g" (general, see also frac to specify this model)
#   frac: for use with general and it is the fraction of the diatance
#       between the 2 homozygots where the hetozygot sits
#   alphamat: alpha (type I error) vector
#   method: "A" for asymptotic, "S" for simulations and "M" for mixt,
#           i.e. if the counts of genotypes having higher than minimum
# 	        penetrance are greater or equalt to nmin (set to 5 here)
#   outcome: either power (for a given sample size) or sample size for a given
#            power
#   nmin: min number of subjects with geonotype of interest that will trigger
#         asymptotic calculations for power under method = "M"
#   targetpow: the target power
#   nsim: the number of sims at each level in the design
#-------------------------------------------------------------------------------
getpowsampleQT <- function(n=n, p=p,hsqmat=NULL,deltamat=NULL,sdmat=rep(1,3),
                         modi=modi, frac=1/4, alphamat=alphamat, method="M",
                         outcome="power", nmin=5, targetpow=0.8, nsim=100){
  ##	neff is the number of levels for the effect
  neff <- max(c(length(hsqmat),length(deltamat)))
  q <- 1-p
  outfit <- NULL
  if(outcome == "power"){
    out <- NULL
    for(i in 1: neff){
      hsq <- hsqmat[[i]]
      delta <- deltamat[[i]]
      mu <- getmeans(p=p,sdmat=sdmat,hsq=hsq,delta=delta,frac=frac, modi=modi)
      ## mu: means vector (mu0,mu1,mu2)
      if(sum(is.na(mu)) == 0){
        ## nhtmp is the counts of genotypes having higher than
        ##	minimum penetrance is greater than nmin (set to 5 here)
        ## i.e. genotypes with 2 disease alleles for recessive and
        ## genotypes with one and 2 alleles for additive and dominant
        if(modi == "r"){
          nhtmp <- n*p^2
        } else {
          nhtmp <- n*(1-q^2)
        }
        outpow <- NULL
        if(method == "A" | (method == "M" & nhtmp >= nmin) ){
          for(alpha in alphamat){
	    pow<-getpowasymptQT(n=n, p=p, mu=mu, sdmat=sdmat, alpha=alpha)
            outpow <- rbind(outpow,c(0,alpha,pow))
          }
        } else {
          pvals <- replicate(nsim,getpvalssimQT(n=n,p=p,mu=mu,sdmat=sdmat,nmin=nmin))
          for(alpha in alphamat){
            outpow <- rbind(outpow,c(1,alpha,sum(pvals<alpha,na.rm=T)/nsim, NA))  
            # The last value is the non-centrality parameter placeholder
          }
        }
        vexp<-mu[[1]]^2*(-1+p)^2-2*mu[[2]]^2*(-1+p)*p+mu[[3]]^2*p^2-
	(mu[[1]]*(-1+p)^2+p*(-2*mu[[2]]*(-1+p)+mu[[3]]*p))^2
	verr<-q^2*sdmat[[1]]^2+2*p*q*sdmat[[2]]^2+p^2*sdmat[[2]]^2
	hsq<-vexp/(verr+vexp)
	delta<-mu[3]-mu[1]
        out <- rbind(out,cbind(p,hsq,delta,mu[[1]],mu[[2]],mu[[3]],sdmat[[1]],
                sdmat[[2]],sdmat[[3]],outpow))
      }
    }
    if(!is.null(out)){
      out <- data.frame(out)
      ##	is.sim is the indicator vector containing the effects that were simulated
      dimnames(out)[[2]] <- c("p","hsq","delta","mu0","mu1","mu2","sd0","sd1",
                               "sd2","is.sim","alpha","pow", "noncenp")
      out$powfit <- out$pow
      if(sum(out$is.sim) >= 1 & max(out$pow) >= 0.1 & max(table(out$alpha))>=3){
        ## if any levels simulated estimated power via	glm method using data from
        ##  all effects
        for(alpha in alphamat){
          out1 <- out[out$alpha==alpha,]
          outlong <- NULL
          for(j in 1:nrow(out1)){
            nind <- round(nsim*(1-out1$pow[[j]])) + round(nsim*out1$pow[[j]])
            outlong <- rbind(outlong,
                             cbind(matrix(out1[j,], nind, ncol(out1),
                                          byrow=T),
                                   rep(c(0,1),c(round(nsim*(1-out1$pow[[j]])),
                                                round(nsim*out1$pow[[j]])))))
          }
          outlong <- data.frame(outlong)
          dimnames(outlong)[[2]] <- c("p","hsq","delta","mu0","mu1","mu2","sd0",
                   "sd1","sd2","is.sim","alpha","pow","noncenp","powfit","powb")
          outlong$noncenp <- rep(NA, nrow(outlong))
          outlong$hsq <- as.numeric(outlong$hsq)
          outlong$powb <- as.numeric(outlong$powb)
          modglm <- glm(powb~log(1+hsq),data=outlong,
                        family <- binomial(link=logit))
          outlong$powfit <- modglm$fitted
          cat(modglm$fitted,"\n")
          outlong1<-unique(outlong)
          cat(out$is.sim,"\n")
          outlong1$powfit[outlong1$is.sim==0] <- outlong1$pow[outlong1$is.sim==0]
          outlong1 <- outlong1[,dimnames(outlong1)[[2]]!="powb"]
          outfit <- rbind(outfit,outlong1)
        }
      } else {
        outfit <- out
      }
      if(!is.null(outfit)){
        outfit$n <- n
      }
    }
  } else {
    ## compute sample size asymptotically
    ## get the noncentrality parameter, per subject used
    for(i in 1: neff){
      hsq <- hsqmat[[i]]
      delta <- deltamat[[i]]
      mu <- getmeans(p=p,sdmat=sdmat,hsq=hsq,delta=delta,frac=frac, modi=modi)
      nonc<-getnoncQT(n=1, p=p, mu=mu,sdmat=sdmat,minsd=10^-5)
      q <- 1-p
      if(sum(is.na(mu)) == 0){
        vexp<-mu[[1]]^2*(-1+p)^2-2*mu[[2]]^2*(-1+p)*p+mu[[3]]^2*p^2-
	(mu[[1]]*(-1+p)^2+p*(-2*mu[[2]]*(-1+p)+mu[[3]]*p))^2
	verr<-q^2*sdmat[[1]]^2+2*p*q*sdmat[[2]]^2+p^2*sdmat[[2]]^2
	hsq<-vexp/(verr+vexp)
	delta<-mu[3]-mu[1]
        for(alpha in alphamat){
          res <- uniroot(pownoncchisqeq, interval=c(1,100), alpha=alpha, df=2,
                         targetpow=targetpow)
          nall <- ceiling(res$root/nonc)
          outfit <- rbind(outfit,c(p,hsq,delta,mu,sdmat,alpha,targetpow,nall, nonc))
        }
      }
    }
    if(!is.null(outfit)){
      outfit <- data.frame(outfit)
      dimnames(outfit)[[2]] <- c("p", "hsq", "delta", "mu0", "mu1", "mu2","sd0",
                                 "sd1","sd2","alpha", "targetpow", "n","noncenp")
    }
  }
  outfit
}

#-------------------------------------------------------------------------------
# Function: getpowallQT
# This is the function for power for all levels of the design (i.e., many MAFs
#  and many modes of inheritance)
#   n: sample size
#   pmat: MAF vector (for disease SNP)
#   hsqmat: hsq (heritabilty or R^2) vector
#   deltapmat: deltap vector
#   sdmat: the vector of sds for distributions for the 3 genotypes 
#   modimat: disease mode vector of inheritance with values "a" (additive),
#            "d" (dominant) and "r" (recessive)
#   frac: the fraction of the distance between low risk homozygote and high 
#         risk homozygote where the heterozygote mean sits, frac =0.25 is 
#         halfway between additive and recessive
#   alphamat: alpha (type I error) vector
#   method: "A" for asymptotic, "S" for simulations and "M" for mixt,
#           i.e. if the counts of genotypes having higher than minimum
# 	        penetrance are greater or equalt to nmin (set to 5 here)
#   outcome: either power (for a given sample size) or sample size for a given
#            power
#   nmin: min number of subjects with geonotype of interest that will trigger
#         asymptotic calculations for power under method = "M"
#   targetpow: the target power
#   nsim: the number of sims at each level in the design
#-------------------------------------------------------------------------------
getpowallQT <- function(n=n, pmat=pmat,  hsqmat=NULL,
                      deltamat=NULL, sdmat=rep(1,3),modimat=modimat,frac=frac,
                      alphamat=alphamat, method="M", outcome="power", nmin=10,
                      targetpow=0.8, nsim=100) {
  outall <- NULL
  for(p in pmat){
    for(modi in modimat){
      out <- getpowsampleQT(n=n, p=p,hsqmat=hsqmat,deltamat=deltamat,sdmat=sdmat,
                         modi=modi, frac=frac, alphamat=alphamat, method=method,
                         outcome=outcome, nmin=nmin, targetpow=targetpow,nsim=nsim)

      if(!is.null(out)){
        out$modi <- rep(modi, nrow(out))
      }
      outall <- rbind(outall, out)
    }
  }
  outall
}

#-------------------------------------------------------------------------------
# Function: plotpowQT
# This is the function to plot for power for all levels of the design
#  (i.e., many MAFs and many modes of inheritance)
#   n: sample size
#   pmat: MAF vector (for disease SNP)
#   K: prevalence of disease
#   hsqmat: hsq (heritabilty or R^2) vector
#   deltapmat: deltap vector
#   sdmat: the vector of sds for distributions for the 3 genotypes 
#   modimat: disease mode vector of inheritance with values "a" (additive),
#            "d" (dominant) and "r" (recessive)
#   frac: the fraction of the distance between low risk homozygote and high 
#         risk homozygote where the heterozygote mean sits, frac =0.25 is 
#         halfway between additive and recessive
#   alphamat: alpha (type I error) vector
#   method: "A" for asymptotic, "S" for simulations and "M" for mixt,
#           i.e. if the counts of genotypes having higher than minimum
# 	        penetrance are greater or equalt to nmin (set to 5 here)
#   outcome: either power (for a given sample size) or sample size for a given
#            power
#   nmin: min number of subjects with geonotype of interest that will trigger
#         asymptotic calculations for power under method = "M"
#   targetpow: the target power
#   nsim: the number of sims at each level in the design
#   datafile: the file name where the data frame containing power is writen
#   plotfile: the file name where the PDF (figures) is writen
#   colors: the colors for each MAF power/sample size curve
#-------------------------------------------------------------------------------
plotpowQT <- function(n=n, pmat=pmat, K=K, hsqmat=NULL,
                    deltamat=NULL,sdmat=rep(1,3), modimat=modimat,frac=1/4,
                    alphamat=alphamat, method="M", outcome="power", nmin=10,
                    targetpow=0.8, nsim=100, datafile=datafile,
                    plotfile=plotfile, colors=colors){
  # Parameter check
  parcheck()
	require(lattice)
	outall <- getpowallQT(n=n, pmat=pmat,  hsqmat=hsqmat,sdmat=sdmat,
                      deltamat=deltamat, modimat=modimat,frac=frac,
                      alphamat=alphamat, method=method, outcome=outcome, nmin=nmin,
                      targetpow=targetpow, nsim=nsim)
  outall0 <- outall
 	if (!is.null(outall)) {
 	  if (is.list(outall)) {
      outall1 <- NULL
      for (i in 1:nrow(outall)) {
        tmp <- unlist(outall[i,])
        outall1 <- rbind(outall1,tmp)
      }
      outall <- outall1
    }
    write.table(outall, file=datafile, sep="\t", na="", col.names=T,
                quote=FALSE, row.names=F ,append=F)
    write("\nSuccessfully writing the result table file. \n", file=logfile,
          append=T)
    outall <- outall0
    rm(outall0); rm(outall1)

    # Now add some variable for better plotting
  	outall$pexp <- factor(paste("RAF=", 100*unlist(outall$p), sep=""),
                             levels=paste("RAF=", 100*pmat, sep=""),ordered=T)
    outall$modiexp <- rep("Additive",nrow(outall))
    outall$modiexp[outall$modi == "d"] <- "Dominant"
    outall$modiexp[outall$modi == "r"] <- "Recessive"
    outall$modiexp[outall$modi == "g"] <- paste("Gen fr=",frac,sep="")
    cat(dim(outall),"\n")
    cat(dimnames(outall)[[2]],"\n")
    pdf(file=plotfile, height=6, width=5)
    if(outcome == "power"){
      if(!is.null(hsqmat)){
			 for(alpha in alphamat){
			   title = paste("Power @ alpha=", alpha," & n=", n, sep="")
         myplot <- xyplot(powfit~100*hsq|modiexp, groups=pexp,
					           data = outall[outall$alpha==alpha,], type="l",
                     ylab = "Power", xlab = expression(paste(R^2,"(%)")),
                     panel = function(x,y,...) {
                            panel.abline(h = c(0.5, 0.8, 0.9), col = "gray", 
                                lwd = 1, lty=2)
                            panel.xyplot(x, y,...) },                        
              			 key = list(title="RAF-Risk Allele Frequecy (%)", lines=1,
             				 lty = c(1:length(levels(as.factor(outall$pexp)))),
            				 col = colors[1:length(levels(as.factor(outall$pexp)))],
              			 text = list(levels(as.factor(outall$pexp))),columns=2),
            				 col = colors[1:length(levels(as.factor(outall$pexp)))],
                     lwd = 1.25,
                     lty = c(1:length(levels(as.factor(outall$pexp)))),
             				 main = title, as.table = T,
                     layout = c(length(levels(factor(outall$modiexp))),1))
              			 print(myplot)
			 }
		  } else if(!is.null(deltamat)){
			 for(alpha in alphamat){
			   	title <- paste("Power @ alpha=", alpha, " & n=", n, sep="")
  				myplot <- xyplot(powfit~delta|modiexp,groups=pexp,
                				data=outall[outall$alpha==alpha,], type="l",
                        ylab="Power", xlab=expression(paste(Delta)),
                        panel = function(x,y,...) {
                                panel.abline(h = c(0.5, 0.8, 0.9), col = "gray", 
                                    lwd = 1, lty=2)
                                panel.xyplot(x, y,...) },                           
              					key=list(title="RAF-Risk Allele Frequecy (%)", lines=1,
              					lty=c(1:length(levels(as.factor(outall$pexp)))),
              					col=colors[1:length(levels(as.factor(outall$pexp)))],
              					text=list(levels(as.factor(outall$pexp))), columns=2),
              					col=colors[1:length(levels(as.factor(outall$pexp)))],
                        lwd=1.25,
                        lty=c(1:length(levels(as.factor(outall$pexp)))),
              					main=title, as.table=T,
                        layout=c(length(levels(factor(outall$modiexp))),1))
				  print(myplot)
			 }
	   }
     write("\nSuccessfully generating the power pdf plots. \n", file=logfile,
          append=T)
	  } else if(outcome=="sample"){
		 if(!is.null(hsqmat)){
			for(alpha in alphamat){
				title <- paste("# of subjects alpha=", alpha, " & power=",
                      targetpow, sep="")
				myplot <- xyplot(n~100*hsq|modiexp, groups=pexp,
					               data = outall[outall$alpha==alpha,], type = "l",
                         ylab = "# of cases (log10 scale)",
              				   xlab = expression(paste(R^2,"(%)")),
                         panel = function(x,y,...) {
                                  panel.abline(h = log10(n), col = "gray", 
                                         lwd = 1, lty=2)
                                  panel.xyplot(x, y,...) },                  				   
                         scales = list(y=list(log=T)),
                				 key = list(title="RAF-Risk Allele Frequecy (%)",
                         lines=1,
                				 lty=c(1:length(levels(as.factor(outall$pexp)))),
              					 col=colors[1:length(levels(as.factor(outall$pexp)))],
               					 text=list(levels(as.factor(outall$pexp))),columns=2),
              					 col = colors[1:length(levels(as.factor(outall$pexp)))],
                         lwd=1.25,
               					 lty=c(1:length(levels(as.factor(outall$pexp)))),
               					 main=title, as.table=T,
                         layout=c(length(levels(factor(outall$modiexp))),1))
				print(myplot)
			}
		 } else if(!is.null(deltamat)){
			for(alpha in alphamat){
				title <- paste("# of subjects alpha=", alpha, " & power=",
                      targetpow, sep="")
				myplot <- xyplot(n~delta|modiexp, groups = pexp,
              					data = outall[outall$alpha == alpha,], type = "l",
                        ylab = "# of cases (log10 scale)",
                        xlab=expression(paste(Delta)),
                        panel = function(x,y,...) {
                                  panel.abline(h = log10(n), col = "gray", 
                                         lwd = 1, lty=2)
                                  panel.xyplot(x, y,...) },                            
                        scales = list(y=list(log=T)),
					              key = list(title="RAF-Risk Allele Frequecy (%)",
                        				lines = 1,
              					lty = c(1:length(levels(as.factor(outall$pexp)))),
              					col = colors[1:length(levels(as.factor(outall$pexp)))],
              					text = list(levels(as.factor(outall$pexp))),columns=2),
              					col = colors[1:length(levels(as.factor(outall$pexp)))],lwd=1.25,
              					lty = c(1:length(levels(as.factor(outall$pexp)))),
              					main = title, as.table=T,
                        layout=c(length(levels(factor(outall$modiexp))),1))
				print(myplot)
			}
		 }
     write("\nSuccessfully generating the sample size estimation pdf plots. \n",
            file=logfile, append=T)
	 }
	} else {
    write("\nThe result table is empty, nothing will be written to the result file. \n",
       file=logfile, append=T)
  }
	dev.off()
}

# End of functions
