## The statistical power of estimating genetic variance or genetic correlation for continuous trait 
## using genome-wide SNPs
## N: general sample size
## n.study.sample: your total sample size
## var.A: Variance of the SNP-derived genetic relationships
## realh2: The power calculation requires the true SNP-heritability, so that the power is the 
## probability of estimating a SNP-heritability that is greater than zero. 

powerh2<-within(expand.grid(N = seq(500,10000,500),
      n.study.sample=4900,
      realh2=seq(0.1,0.4,0.02),
      KEEP.OUT.ATTRS = FALSE),
      { alpha <- 0.05
     	var.A <- 0.00002
        var.sampling <- 2/(N^2*var.A)
        NCP <- realh2^2/var.sampling
	    critval <- qchisq(alpha, df = 1, lower.tail = FALSE)
	    powerh2 <- pchisq(critval, df = 1, ncp = NCP, lower.tail = FALSE)})


library(ggplot2) #required library for plotting
pow <- powerh2
pow$alpha <- paste("alpha", pow$alpha, sep = "=")
pow$h2 <- as.factor(as.character(pow$realh2))

#main title
title <- paste("The probability of detecting (h^2 > 0) \n for continuous trait. \n Vertical dashed Red line: n=",pow$n.study.sample[1],sep = "")
#x label
xlab <- "Sample size"
#y label
ylab <- "Power"
#output file name
pdf.file <- "PowerPlot.pdf"
pdf(file=pdf.file, height = 4.5, width = 5)
zp1 <- ggplot(pow ,  # data for plotting
  aes( y = powerh2, x = N, colour = h2,linetype = h2)) + 
  facet_grid(. ~ alpha , scale = "free") +  #grid by alpha level 
  geom_line() + 
  ylab(ylab)+  #y label
  xlab(xlab) +  #x label 
  ggtitle(title) +  theme_bw() +
  geom_line(size = 1.25)+theme(plot.title=element_text(family="Times", face="bold", size=12))+
  geom_hline(yintercept = 0.8,color="lightskyblue",linetype=2)+geom_hline(yintercept = 0.9,color="lightpink",linetype=2)+
  geom_vline(aes(xintercept=n.study.sample[1]), colour="#BB0000", linetype="dashed")
print(zp1)
graphics.off()
