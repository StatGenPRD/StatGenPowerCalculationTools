## file: Power_contour.R
## author: Matthew R. Nelson
## date created: 08/22/2014


## ncp: Calculate NCP for a binary trait test
## n: total sample size (both groups)
## p: allele frequency if additive, otherwise genotype frequency
## eff: effect in standard deviations, OR, HR (per allele if lAdd = TRUE) for
##  continuous, binary, and survival traits, respectively
## f: proportion of the sample with target response (i.e. cases or events)
## trait: one of "cont", "bin", or "surv" for continuous, binary, or survival traits,
##  respectively
## lAdd: whether genetic model is additive, otherwise dom/rec
ncp = function(n, p, eff, f = 1, trait = c("cont", "bin", "surv")[2], lAdd = TRUE) {
  x = 1
  if(lAdd)
    x = 2
  
  if(trait == "cont")
    ncp = n * x * p * (1 - p) * eff^2
  else if(trait == "bin")
    ncp = n * x * p * (1 - p) * f * (1 - f) * log(eff)^2
  else if(trait == "surv")
    ncp = n * x * p * (1 - p) * f * log(eff)^2
  
  return(ncp)
}
  
## power.calc.ncp: Calculate the power of a test at a given testwise alpha level
## alpha: testwise significance threshold
## ncp: non-centrality parameter
power.calc.ncp = function(alpha, ncp) {
  critval = qchisq(alpha, df = 1, lower.tail = FALSE)
  pchisq(critval, df = 1, ncp = ncp, lower.tail = FALSE)
}

## power.calc: Calculate the power to test a given trait effect
## n: total sample size (both groups)
## p: allele frequency if additive, otherwise genotype frequency
## eff: effect in standard deviations, OR, HR (per allele if lAdd = TRUE) for
##  continuous, binary, and survival traits, respectively
## f: proportion of the sample with target response (i.e. cases or events)
## trait: one of "cont", "bin", or "surv" for continuous, binary, or survival traits,
##  respectively
## alpha: the testwise significance threshold
power.calc = function(n, p, eff, f = 1, trait = c("cont", "bin", "surv")[2], 
                          lAdd = TRUE, alpha = 0.05) {
  ncp = ncp(n, p, eff, f, trait, lAdd)
  power.calc.ncp(alpha, ncp)
}


## ncp.find: Find the NCP for a test to have the specified power at the specified 
##  significance
## alpha: testwise significance threshold
## pwr: desired power of the test
ncp.find = function(alpha, pwr) {
  critval = qchisq(alpha, df = 1, lower.tail = FALSE)
  fxn = function(ncp, q) pchisq(q, ncp = ncp, df = 1, lower.tail = FALSE) - pwr
  ncp = uniroot(fxn, c(0, 1000), q = critval)$root
  return(ncp)
}


## eff.size: Derive the effect size corresponding to a test with specified
##  parameters and NCP to yield the desired power (from ncp.find)
## n: total sample size (both groups)
## p: allele frequency if additive, otherwise genotype frequency
## f: proportion of the sample with target response (i.e. cases or event)
## ncp: non-centrality parameter
## trait: one of "cont", "bin", or "surv" for continuous, binary, or survival traits,
## lAdd = whether genetic model is additive, otherwise dom/rec
eff.size = function(n, p, f = 1, ncp, trait = c("cont", "bin", "surv")[2], 
                   lAdd = TRUE) {
  x = 1
  if(lAdd)
    x = 2
    
  if(trait == "cont")
    eff = sqrt(ncp/(n * x * p * (1 - p)))
  else if(trait == "bin")
    eff = exp(sqrt(ncp/(n * x * p * (1 - p) * f * (1 - f))))
  else if(trait == "surv")
    eff = exp(sqrt(ncp/(n * x * p * (1 - p) * f)))
              
  return(eff)
}


## Create a power contour plot of the variant effects that the study had the specified
##  power to identify assuming an additive model
## power = vector of desired power levels
## alpha = testwise significance threshold
## n = total sample size (both groups)
## f = proportion of the sample with target response (i.e. cases)
## lAdd = whether genetic model is additive, otherwise dom/rec
## ylim = custom limits on y-axis, otherwise drawn from observed data
## p = allele frequency if additive, otherwise genotype frequency
power.contour = function(power = c(0.5, 0.8, 0.95), alpha, n, f, 
                         trait = c("cont", "bin", "surv")[2], lAdd = TRUE,
                         ylim = NULL, p = NULL) {
  require(ggplot2)
  require(reshape2)
  
  if(is.null(p))
    p = seq(0.01, 0.5, by = 0.01)
  
  eff = matrix(NA, nrow = length(p), ncol = length(power),
               dimnames = list(NULL, paste("Power", power, sep = "_")))
  for(i in seq(along = power)) {
    ncp = ncp.find(alpha, power[i])
    eff[,i] = eff.size(n, p, f, ncp, trait, lAdd)
  }
  data = data.frame(Freq = p, eff)
  data.long = melt(data, id = "Freq", value.name = "Effect")
  data.long$Power = factor(rep(power, each = length(p)))
  
  
  if(lAdd)
    xlab = "Minor Allele Frequency"
  else
    xlab = "Genotype Frequency"
  
  if(trait == "bin")
    ylab = "Odds Ratio"
  else if(trait == "surv")
    ylab = "Hazard Ratio"
  else
    ylab = "Standard Deviations"
  
  if(trait %in% c("bin", "surv")) {
    if(is.null(ylim))
      ylim = c(1, max(data.long$Effect, na.rm = TRUE))
    
    g = ggplot(data.long, aes(y = Effect, x = Freq, group = Power, colour = Power)) + 
      geom_line(size = 1.25) + xlab(xlab) + ylab(ylab) + 
      theme_bw() + theme(legend.position = "top") + 
      scale_y_log10(limits = ylim, expand = c(0, 0.025)) +
      annotation_logticks()
  }
  else {
    if(is.null(ylim))
      ylim = c(0, max(data.long$Effect, na.rm = TRUE))
    
    g = ggplot(data.long, aes(y = Effect, x = Freq, group = Power, colour = Power)) + 
      geom_line(size = 1.25) + xlab(xlab) + ylab(ylab) + 
      ylim(ylim[1], ylim[2]) + 
      theme_bw() + theme(legend.position = "top")
  }
  
  print(g)
  
  invisible(data)
}


## laprash.example: Example use of power.contour.bin based on the TEACH rash 
##  association study
laprash.example = function(lPdf = TRUE){
  n = 72 + 651 # grade 3+4 rash versus none
  f = 72/n
  
  ntests = 157
  alpha = 0.05/ntests
  if(lPdf) {
    pdf("lapatinib_teach_rash_power.pdf", height = 4.5, width = 6.5)
    on.exit(dev.off())
  }
  power.contour(alpha = alpha, n = n, f = f, trait = "bin", lAdd = TRUE, 
                    p = seq(0.01, 0.2, by = 0.001))
  power.contour(alpha = alpha, n = n, f = f, trait = "bin", lAdd = FALSE, 
                    p = seq(0.01, 0.2, by = 0.001))
  
}