# Last Update: 12/01/2008 10:00AM

#-------------------------------------------------------------------------------
# title: Power_CC.R
# $Source: Power_CC.R,v $
# $Revision: 1.0 $
# $Date: 2008/11/20 10:00:00 $
# $Author: Silviu-Alin Bacanu & Judong Shen
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# README for Power_CC.R
# Power_CC.R includes all the functions that can be used to do both power 
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
  if (sum(is.null(hsqmat),is.null(ORmat))!=1) {
    write("\n(Only) one of hsqmat and ORmat should be 'NULL'.\n", file=logfile,
          append=T)
  }
  if (K < 0 | K > 1) {
    write("\n Warning: the K (disease prevalence) should be set b/w 0 and 1.\n", 
          file=logfile, append=T)    
    stop("K (disease prevalence) is not b/w 0 and 1")          
  }
  if (any(nca<0,nco<0,npc<0)) {
    write("\nWarning: any of 'nca', 'nco' and 'npc' should NOT be set as negative.\n", 
            file=logfile, append=T) 
    stop("Negative number found in 'nca', 'nco' and 'npc' parameters")   
  }
  if (sum(nca>0,nco>0,npc>0) < 2) {
    write("\nWarning: at least two of 'nca', 'nco' and 'npc' should be set as greater than 0.\n", 
            file=logfile, append=T)
    stop("At least two of 'nca', 'nco' and 'npc' parameters should be > 0 ")    
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
# Function: getsampbin
# This function simulates random samples of cases, cntls and pop cntls
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   p: MAF (for disease SNP)
#   f: penetrance vector (f0,f1,f2)
#   K: prevalence of disease
#-------------------------------------------------------------------------------
getsampbin <- function(nca=0, nco=0, npc=0, p=p, f=rep(0.01,3), K=K){
	f0 <- f[[1]]
	f1 <- f[[2]]
	f2 <- f[[3]]
	pheno <- NULL
	geno <- NULL
	if(nca > 0){
		ca <- sample(c(0,1,2),size=nca,prob=c((1-p)^2*f0/K,2*p*(1-p)*f1/K,p^2*f2/K),
          replace=T)
		geno <- c(geno,ca)
		pheno <- c(pheno,rep(2,nca))
	}
	if(nco > 0){
		co <- sample(c(0,1,2),size=nco,prob=c((1-p)^2*(1-f0)/(1-K),
          2*p*(1-p)*(1-f1)/(1-K),p^2*(1-f2)/(1-K)),replace=T)
		geno <- c(geno,co)
		pheno <- c(pheno,rep(1,nco))
	}
	if(npc > 0){
		pc <- sample(c(0,1,2),size=npc,prob=c((1-p)^2,2*p*(1-p),p^2),replace=T)
		geno <- c(geno,pc)
		pheno <- c(pheno,rep(0,npc))
	}
	data <- cbind(pheno,geno)
	data
}

#-------------------------------------------------------------------------------
# Function: getpens
# This function computes penetrances for a causal polymorphism with the above
#  parameters.
#   pd: MAF
#   K: prevalence of disease
#   hsq: heritabilty (R^2)
#   OR: odds ratios of the count table formed by disease status and the 2
#       homozygots
#   phom: NOT yet implemented
#   modi: mode of inheritance with values "a" (additive), "d" (dominant) and
#         "r" (recessive)
#   f: penetrance vector
# Please note effect size is given by either hsq and OR (note one must be null).
#-------------------------------------------------------------------------------
getpens <- function(pd=pd, K=K, hsq=NULL, OR=NULL, phom=c(NA,NA), modi=modi){
	qd <- 1-pd
	if(!is.null(hsq)){
		if(hsq == 0){
			delta1 <- 0; delta2 <- 0
		} else {
			if(modi == "r"){
				delta2 <- sqrt(hsq)/(pd*sqrt(1/((1-K)*K)-pd^2/((1-K)*K)))
				delta1 <- 0
			} else if(modi == "d") {
				delta1 <- sqrt(hsq)/sqrt((2*pd)/((1-K)*K)-(5*pd^2)/((1-K)*K)+
				          (4*pd^3)/((1-K)*K)-pd^4/((1-K)*K))
				delta2 <- 0
			} else if (modi == "a") {
				delta1 <- sqrt(hsq)/sqrt((2*pd)/((1-K)*K)-(2*pd^2)/((1-K)*K))
				delta2 <- delta1
			}
		}
	} else if(!is.null(OR)){
		if(OR == 1){
			delta1 <- 0; delta2 <- 0
		} else {
			if(modi == "r"){
				 delta2 <- (-1+K-K*OR+(-1+2*K)*(-1+OR)*pd^2+
     					     sqrt((1+K*(-1+OR))^2-2*(-1+OR)*(-1+K+K*OR)*pd^2+
					         (-1+OR)^2*pd^4))/(2.*(-1+OR)*pd^2*(-1+pd^2))
				delta1<-0
			} else if (modi == "d") {
				 delta1 <- (-1+(-1+OR)*(-2+pd)*pd-K*(-1+OR)*(1+2*(-2+pd)*pd)+
          				  sqrt(K^2*(-1+OR)^2+(-1+(-1+OR)*(-2+pd)*pd)^2+
				            2*K*(-1+OR)*(1+(1+OR)*(-2+pd)*pd)))/
          					(2.*(-1+OR)*(-2+pd)*(-1+pd)^2*pd)
				delta2<-0
			} else if (modi == "a") {
				 delta1 <- (-1+pd-OR*pd+K*(-1+OR)*(-1+2*pd)+
          					sqrt((1+K*(-1+OR))^2-2*(-1+OR)*(-1+K+K*OR)*pd+
					         (-1+OR)^2*pd^2))/(4.*(-1+OR)*(-1+pd)*pd)
				delta2<-delta1
			}
		}
	} else if (sum(!is.na(phom))==2){
		p1 <- (K*n2*p0*p2 + (-1 + K)*n1*(-(frac*p2) + p0*(-1 + frac + p2)))/
		      ((-1 + K)*n1*(-1 + frac*(p0 - p2) + p2)+ K*n2*(frac*(p0 - p2) + p2))
	}
	f0 <- K-(1-qd^2)*delta1-pd^2*delta2
	f1 <- f0+delta1
	f2 <- f1+delta2
	f <- c(f0,f1,f2)
	if(min(f)<0 | max(f)>1){
		rep(NA,3)
	} else {
		f
	}
}

#-------------------------------------------------------------------------------
# Function: getpvalssim
# This function gets the exact p-values from simulated samples
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   p: MAF (for disease SNP)
#   f: penetrance vector (f0,f1,f2)
#   K: prevalence of disease
#-------------------------------------------------------------------------------
getpvalssim <- function(nca=0, nco=0, npc=0, p=p, f=f, K=K){
	f0 <- f[[1]]
	f1 <- f[[2]]
	f2 <- f[[3]]
	if(min(f) < 0 | max(f) > 1){
		pvals<-NA
	} else {
		data <- getsampbin(nca=nca,nco=nco,npc=npc,p=p,f=f,K=K)
		if(sd(data[,2],na.rm=T) < 0.01) {
			pvals<-NA
		} else {
      #	compute the fisher exact test p-value
			pvals<-fisher.test(table(data[,1],data[,2]))$p.value
		}
	}
	pvals
}

#-------------------------------------------------------------------------------
# Function: getnonc
# This function gets the asymptotic noncentrality of case-controls or
#  case-population controls.
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   p: MAF (for disease SNP)
#   f: penetrance vector (f0,f1,f2)
#   K: prevalence of disease
# Please note if both clinical controls and population controls are available
#  then the two are combined
#-------------------------------------------------------------------------------
getnonc <- function(nca=0, nco=0, npc=0, p=p, f=f, K=K){
	f0 <- f[[1]]
	f1 <- f[[2]]
	f2 <- f[[3]]
	pvec <- c((1-p)^2,2*p*(1-p),p^2)
	tbca <- rep(0,3)
	tbco <- rep(0,3)
	tbpc <- rep(0,3)
	if(nca > 0){
		probs <- f*pvec/K
		tbca <- probs/sum(probs)
	}
	if(nco > 0){
		probs <- (1-f)*pvec/(1-K)
		tbco <- probs/sum(probs)
	}
	if(npc > 0){
		tbpc <- pvec
	}
	if(nca > 0){
		alltab <- rbind(tbca,(nco*tbco+npc*tbpc)/(nco+npc))
		Na <- nca
		Nu <- (nco+npc)
	} else {
		alltab <- rbind(tbpc,tbco)
		Na <- npc
		Nu <- nco
	}
	nonc <- Na*Nu*sum((alltab[1,]-alltab[2,])^2/(Na*alltab[1,]+Nu*alltab[2,]))
	nonc
}

#-------------------------------------------------------------------------------
# Function: getpowasympt
# This function gets the asymptotic power from the sample of case-controls or
#  case-population controls.
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   p: MAF (for disease SNP)
#   f: penetrance vector (f0,f1,f2)
#   K: prevalence of disease
#   alpha: type I error
# Please note if both clinical controls and population controls are available
#  then the two are combined
#-------------------------------------------------------------------------------
getpowasympt <- function(nca=0, nco=0, npc=0, p=p, f=f, K=K, alpha=10^-7){
  r <- 3   # r is the number of groups
  nonc <- getnonc(nca=nca,nco=nco,npc=npc,p=p,f=f,K=K)
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
# Function: getpowsample
# This is the master function for power for samples of case-controls or
#  case-population controls.
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   p: MAF (for disease SNP)
#   K: prevalence of disease
#   hsqmat: hsq (heritabilty or R^2) vector
#   ORmat: OR vector
#   deltapmat: deltap vector
#   modi: mode of inheritance with values "a" (additive), "d" (dominant) and
#         "r" (recessive)
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
# Please note if both clinical controls and population controls are available
#  then the two are combined
#-------------------------------------------------------------------------------
getpowsample <- function(nca=0, nco=0, npc=0, p=p, K=K, hsqmat=NULL, ORmat=NULL,
                         deltapmat=NULL, modi=modi, alphamat=alphamat, method="M",
                         outcome="power", nmin=10, targetpow=0.8, nsim=100){
  ##	neff is the number of levels for the effect
  neff <- max(c(length(hsqmat),length(ORmat), length(deltapmat)))
  q <- 1-p
  outfit <- NULL
  if(outcome == "power"){
    out <- NULL
    for(i in 1: neff){
      hsq <- hsqmat[[i]]
      OR <- ORmat[[i]]
      deltap <- deltapmat[[i]]
      f <- getpens(pd=p,K=K,hsq=hsq,OR=OR,phom=c(NA,NA),modi=modi)
      ## f: penetrance vector (f0,f1,f2)
      f0 <- f[[1]]
      f1 <- f[[2]]
      f2 <- f[[3]]
      if(sum(is.na(f)) == 0){
        ## nhtmp is the counts of genotypes having higher than
        ##	minimum penetrance is greater than nmin (set to 5 here)
        ## i.e. genotypes with 2 disease alleles for recessive and
        ## genotypes with one and 2 alleles for additive and dominant
        if(modi == "r"){
          nhtmp <- nca*p^2*f2/K+nco*p^2*(1-f2)/(1-K)+npc*p^2
        } else {
          nhtmp <- nca*(p^2*f2/K+2*p*q*f1/K)+nco*(p^2*(1-f2)/
                                                  (1-K)+2*p*q*(1-f1)/(1-K))+npc*(1-q^2)
        }
        outpow <- NULL
        if(method == "A" | (method == "M" & nhtmp >= nmin) ){
          for(alpha in alphamat){
            pow <- getpowasympt(nca=nca,nco=nco,npc=npc,p=p,f=f,K=K,alpha)
            outpow <- rbind(outpow,c(0,alpha,pow))
          }
        } else {
          pvals <- replicate(nsim,getpvalssim(nca=nca,nco=nco,npc=npc,p=p,f=f,K=K))
          for(alpha in alphamat){
            outpow <- rbind(outpow,c(1,alpha,sum(pvals<alpha,na.rm=T)/nsim, NA))  # The last value is the non-centrality parameter placeholder
          }
        }
        ORr <- (1-f0)*f2/f0/(1-f2)
        deltaf <- f2-f0
        vexp <- f0^2*(-1+p)^2-2*f1^2*(-1+p)*p+f2^2*p^2-
          (f0*(-1+p)^2+p*(-2*f1*(-1+p)+f2*p))^2
        hsq <- vexp/(K*(1-K))
        out <- rbind(out,cbind(p,hsq,ORr,deltaf,f0,f1,f2,outpow))
      }
    }
    if(!is.null(out)){
      out <- data.frame(out)
      ##	is.sim is the indicator vector containing the effects that were simulated
      dimnames(out)[[2]] <- c("p","hsq","OR","deltaf","f0","f1","f2","is.sim",
                              "alpha","pow", "noncenp")
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
          dimnames(outlong)[[2]] <- c("p","hsq","OR","deltaf","f0","f1","f2",
                                      "is.sim", "alpha","pow","noncenp","powfit","powb")
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
        outfit$nca <- nca
        outfit$nco <- nco
        outfit$npc <- npc
      }
    }
  } else {
    ## compute sample size asymptotically
    ## get the noncentrality parameter, per subject used
    rca <- nca/(nca+nco+npc)
    rco <- nco/(nca+nco+npc)
    rpc <- npc/(nca+nco+npc)
    for(i in 1: neff){
      hsq <- hsqmat[[i]]
      OR <- ORmat[[i]]
      deltap <- deltapmat[[i]]
      hsq <- hsqmat[[i]]
      OR <- ORmat[[i]]
      deltap <- deltapmat[[i]]
      f <- getpens(pd=p,K=K,hsq=hsq,OR=OR,phom=c(NA,NA),modi=modi)
      nonc <- getnonc(nca=rca,nco=rco,npc=rpc,p=p,f=f,K=K)
      f0 <- f[[1]]
      f1 <- f[[2]]
      f2 <- f[[3]]
      q <- 1-p
      if(sum(is.na(f)) == 0){
        ORr <- (1-f0)*f2/f0/(1-f2)
        deltaf <- f2-f0
        vexp <- f0^2*(-1+p)^2-2*f1^2*(-1+p)*p+f2^2*p^2-
          (f0*(-1+p)^2+p*(-2*f1*(-1+p)+f2*p))^2
        hsq <- vexp/(K*(1-K))
        for(alpha in alphamat){
          res <- uniroot(pownoncchisqeq, interval=c(1,100), alpha=alpha, df=2,
                         targetpow=targetpow)
          nall <- res$root/nonc
          nca <- ceiling(nall*rca)
          nco <- ceiling(nall*rco)
          npc <- ceiling(nall*rpc)
          outfit <- rbind(outfit,c(p,hsq,ORr,deltaf,f,alpha,targetpow,nca,nco,npc, nonc))
        }
      }
    }
    if(!is.null(outfit)){
      outfit <- data.frame(outfit)
      dimnames(outfit)[[2]] <- c("p", "hsq", "OR", "deltaf", "f0", "f1", "f2",
                                 "alpha", "targetpow", "nca", "nco", "npc",
                                 "noncenp")
    }
  }
  outfit
}

#-------------------------------------------------------------------------------
# Function: getpowall
# This is the function for power for all levels of the design (i.e., many MAFs
#  and many modes of inheritance)
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   pmat: MAF vector (for disease SNP)
#   K: prevalence of disease
#   hsqmat: hsq (heritabilty or R^2) vector
#   ORmat: OR vector
#   deltapmat: deltap vector
#   modimat: disease mode vector of inheritance with values "a" (additive),
#            "d" (dominant) and "r" (recessive)
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
# Please note if both clinical controls and population controls are available
#  then the two are combined
#-------------------------------------------------------------------------------
getpowall <- function(nca=0, nco=0, npc=0, pmat=pmat, K=K, hsqmat=NULL,
                      ORmat=NULL, deltapmat=NULL,modimat=modimat,
                      alphamat=alphamat, method="M", outcome="power", nmin=10,
                      targetpow=0.8, nsim=100) {
  outall <- NULL
  for(p in pmat){
    for(modi in modimat){
      out <- getpowsample(nca=nca, nco=nco, npc=npc, p=p, K=K, hsqmat=hsqmat,
                          ORmat=ORmat, deltapmat=deltapmat, modi=modi,
                          alphamat=alphamat, method=method, outcome=outcome,
                          nmin=nmin, targetpow=targetpow, nsim=nsim)
      if(!is.null(out)){
        out$modi <- rep(modi, nrow(out))
      }
      outall <- rbind(outall, out)
    }
  }
  outall
}

#-------------------------------------------------------------------------------
# Function: plotpow
# This is the function to plot for power for all levels of the design
#  (i.e., many MAFs and many modes of inheritance)
#   nca: number of cases
#   nco: number of clinical controls
#   npc: number of population controls
#   pmat: MAF vector (for disease SNP)
#   K: prevalence of disease
#   hsqmat: hsq (heritabilty or R^2) vector
#   ORmat: OR vector
#   deltapmat: deltap vector
#   modimat: disease mode vector of inheritance with values "a" (additive),
#            "d" (dominant) and "r" (recessive)
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
# Please note if both clinical controls and population controls are available
#  then the two are combined
#-------------------------------------------------------------------------------
plotpow <- function(nca=0, nco=0, npc=0, pmat=pmat, K=K, hsqmat=NULL,
                    ORmat=NULL, deltapmat=NULL, modimat=modimat,
                    alphamat=alphamat, method="M", outcome="power", nmin=10,
                    targetpow=0.8, nsim=100, datafile=datafile,
                    plotfile=plotfile, colors=colors){
  # Parameter check                    
  parcheck()       
	outall <- getpowall(nca=nca, nco=nco, npc=npc, pmat=pmat, K=K, hsqmat=hsqmat,
                      ORmat=ORmat, deltapmat=deltapmat, modimat=modimat,
                      alphamat=alphamat, method=method, outcome=outcome,
                      nmin=nmin, targetpow=targetpow, nsim=nsim)
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
    pdf(file=plotfile, height=6, width=5)
    if(outcome == "power"){
      if(!is.null(hsqmat)){
			 for(alpha in alphamat){
			   title = paste("Power @ alpha=", alpha," & nca=", nca,", nco=",
                      nco, ", npc=", npc, sep="")
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
		  } else if(!is.null(ORmat)){
			 for(alpha in alphamat){
			   	title <- paste("Power @ alpha=", alpha, " & nca=", nca, ", nco=", nco,
                      ", npc=", npc, sep="")
  				myplot <- xyplot(powfit~OR|modiexp,groups=pexp,
                				data=outall[outall$alpha==alpha,], type="l",
                        ylab="Power", xlab="OR",
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
				title <- paste("Number of cases for alpha=", alpha, " & power=",
                      targetpow, sep="")
				myplot <- xyplot(nca~100*hsq|modiexp, groups=pexp,
					               data = outall[outall$alpha==alpha,], type = "l",
                         ylab = "# of cases (log10 scale)",
              				   xlab = expression(paste(R^2,"(%)")), 
                         scales = list(y=list(log=T)),
		                      panel = function(x,y,...) {
                                  panel.abline(h = log10(nca), col = "gray", 
                                         lwd = 1, lty=2)
                                  panel.xyplot(x, y,...) },                            
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
		 } else if(!is.null(ORmat)){
			for(alpha in alphamat){
				title <- paste("Number of cases for alpha=",alpha," & power=",targetpow,sep="")
				myplot <- xyplot(nca~OR|modiexp, groups = pexp,
              					data = outall[outall$alpha == alpha,], type = "l",
                        ylab = "# of cases (log10 scale)", xlab="OR",
                        scales = list(y=list(log=T)),
                        panel = function(x,y,...) {
                                  panel.abline(h = log10(nca), col = "gray", 
                                         lwd = 1, lty=2)
                                  panel.xyplot(x, y,...) },                                        
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
