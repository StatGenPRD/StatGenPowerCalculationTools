
# Title: Power_run.R
# Author: Silviu-Alin Bacanu & Judong Shen
# Date:  2008/12/01

# Setup working directory, user SHOULD change it to his/her working directory.
setwd("C:\\Power_CC")

library(lattice)

# Source R functions and R parameter file.
source("Power_CC.R")
source("Power_CC_arg.R")

# Start to write in log file
write("\nLog file for running of Power.R\n",file=logfile,append=F)
write(paste("R.version()=",R.Version(),sep=""),file=logfile,append=T)

# Run main program
plotpow(nca=nca, nco=nco, npc=npc, pmat=pmat, K=K, hsqmat=hsqmat, ORmat=ORmat,
        deltapmat=NULL, modimat=modimat, alphamat=alphamat, method=method,
        outcome=outcome, nmin=nmin, targetpow=targetpow, nsim=nsim, datafile=datafile,
        plotfile=plotfile, colors=colors)           
