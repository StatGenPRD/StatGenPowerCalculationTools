# Please note: only one of hsqmat or deltamat should be specified and the other
  #  one should be set to NULL
  # hsqmat is the vector(matrix) of heritabilities if given else NA
  # heritability is the prop of variance explained by disease SNP (R^2 in LD lingo)
  # Other setting example: hsqmat = NULL
hsqmat = seq(0.005,0.1,0.005)
  # deltamat is the vector of deltas's, i.e. the difference between the means of 
  # high risk and low risk homozygots.
  # Other setting example: deltamat = c(0.5,1,1.5,3,5)
  #  or deltamat = seq(0.25,5,0.25) or deltamat = 0.1*2^c(1:6)
deltamat = NULL
  # sd is the common standard deviation of distributions for the 3 genotypes or 
  # the vector of the 3 different standard deviations of the distributions of 
  # the QT for the 3 genotypes sd of length 3 triggers simulation.
sd = 1
  # sdmat: the vector of sds for distributions for the 3 genotypes  
sdmat = rep(1,3)
  # pmat is the vector of disease MAF
pmat = c(0.02,0.05,0.1,0.2,0.3,0.5)
  # alphamat-vector of type I error
alphamat = c(10^-2,10^-4,10^-7)
  # modimat is the vector of genetic models we need to run
  # "g' stands for general model and it is available just for when deltamat is given
modimat = c("a","d","r","g")
  # frac is the fraction of the distance between low risk homozygote and high 
  # resk homozygote where the heterozygote mean sits
  # frac =0.25 is halfway between additive and recessive
frac=0.25
  # n sample size > 0
  # n should be positive.
n = 1000
  # method can be set as "A" for asymptotic, "M" for mixt and "S" for simulation
method = "S"
  # outcome denotes if power or sample size estimation is needed
  # outcome can be set as either "power" or "sample"
outcome = "sample"
  # min # of subjects with geonotype of interest that will trigger
  # asymptotic calculations for power under method="M"
nmin = 5
  # nsim is the number of sims at each level in the design,
  # 1000 is the recommended number, but you can use 100 for fast run.
nsim = 100
  # targetpow is the target power (only for sample size calculation,
  # i.e., while outcome="sample"), and the value should be b/w 0 and 1
  # when outcome="power", the targetpow value specified below will be ignored
  # whatever it is set.
targetpow = 0.8
  # colors are the colors for each MAF power/sample size curve
colors = c("black","dark red","orange","green","blue","magenta")
  # datafile is the file where to write the data
datafile = "Power_result.txt"
  # plotfile is the PDF file where to plot the power/sample size
plotfile = "Power_plot.pdf"
logfile = "Power.log"
