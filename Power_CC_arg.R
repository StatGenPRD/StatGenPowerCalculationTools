
  # Please note: only one of hsqmat and ORmat should be specified and the other 
  #  one should be set to NULL
  # hsqmat is the vector(matrix) of heritabilities if given else NA
  # heritability is the prop of variance explained by disease SNP (R^2 in LD lingo)
  # Other setting example: hsqmat = NULL
hsqmat = seq(0.005,0.1,0.005)
  # ORmat is the vector of OR's of the table of the two homozygotes
  # Other setting example: ORmat = c(1.5, 3, 4.5, 6, 9, 15) 
  #  or ORmat = seq(2,50,1) or ORmat = 2^c(1:6)
ORmat = NULL
  # pmat is the vector of disease MAF
pmat = c(0.02,0.05,0.1,0.2,0.3,0.5)
  # deltap is unused as of now-Mike Mosteller's diff of proportion in
  #	cases and controls within the 2 homozygots
deltap = NA
  # alphamat-vector of type I error
alphamat = c(10^-2,10^-4,10^-7)
  # modimat is the vector of genetic models we need to run
modimat = c("a","d","r")
  # K is disease prevalence: b/w 0 and 1
K = 0.05
  # nca-# of cases, nco-# of controls, npc-# of population controls
  # All should not be negative and at least two should be positive.
nca = 100
npc = 0
nco = 290
  # method can be set as "A" for asymptotic, "M" for mixt and "S" for simulation
method = "A"  
  # outcome denotes if power or sample size estimation is needed
  # outcome can be set as either "power" or "sample"
outcome = "sample"  
  # min # of subjects with geonotype of interest that will trigger
  # asymptotic calculations for power under method="M"
nmin = 5
  # nsim is the number of sims at each level in the design, 
  # 1000 is the recommended number, but you can use 100 for fast run.
nsim = 1000 
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
