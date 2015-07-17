
## The statistical power of estimating genetic variance or genetic correlation for binary trait 
## using genome-wide SNPs
## N: general sample size
## n.study.sample: your total sample size
## var.A: Variance of the SNP-derived genetic relationships
## realh2: The power calculation requires the true SNP-heritability, so that the power is the 
## probability of estimating a SNP-heritability that is greater than zero. 
## K: disease risk in population
## fcase: proportion of the sample with target response (i.e. cases)

powerCCh2<-within(expand.grid(N = seq(500,10000,500),
      n.study.sample=3000,
      realh2=c(0.1,0.2,0.3,0.4),
      KEEP.OUT.ATTRS = FALSE),
      { alpha <- 0.05
        K <- 0.1
        fcase= 0.5
        var.A=0.00002
        height<-dnorm(qnorm(K,lower.tail=F))
        critval <- qchisq(alpha, df = 1, lower.tail = FALSE)
        var.sampling <- (2*(1-K)^4*(N*fcase+N*(1-fcase))^2)/((N*fcase)^2*(N*(1-fcase))^2*(height/K)^4*var.A)
        NCP <- realh2^2/var.sampling
        powerCCh2 <- pchisq(critval, df = 1, ncp = NCP, lower.tail = FALSE)})

library(ggplot2) #required library for plotting
pow <- powerCCh2
pow$alpha <- paste("alpha", pow$alpha, sep = "=")
pow$h2 <- as.factor(as.character(pow$realh2))

#main title
title <- paste("The probability of detecting (h^2 > 0) \n for binary trait. \n Vertical dashed Red line: Cases=", pow$n.study.sample[1]*pow$fcase," & Controls=",pow$n.study.sample[1]*(1-pow$fcase)," \n Disease in population: ", pow$K*100,"%",", Proportion of cases: ",pow$fcase[1]*100,"%",sep = "")
xlab <- "Total number of cases and controls"
#y label
ylab <- "Power"
#output file name
pdf.file <- "PowerPlot.pdf"
pdf(file=pdf.file, height = 4.5, width = 5)
zp1 <- ggplot(pow ,  # data for plotting
  aes( y = powerCCh2, x = N, colour = h2,linetype = h2)) + 
  facet_grid(. ~ alpha , scale = "free") +  #grid by alpha level 
  geom_line() + 
  ylab(ylab)+  #y label
  xlab(xlab) +  #x label 
  ggtitle(title) +  theme_bw() +
  geom_line(size = 1.25)+theme(plot.title=element_text(family="Times", face="bold", size=10))+
  geom_hline(yintercept = 0.8,color="lightskyblue",linetype=2)+geom_hline(yintercept = 0.9,color="lightpink",linetype=2)+
  geom_vline(aes(xintercept=n.study.sample[1]), colour="#BB0000", linetype="dashed")
print(zp1)
graphics.off()

