
########################################
### Read in simulation results and format results

path.sim <- paste0(path,'Simulation Results/')
source(paste0(path,'GitHub_Functions_ComplianceMixture.R'))

###########################################################################################
###Fixed method sim results, examining estimation and performance when dose-levels a fixed to IND or REL
##Used to create Figure 1 in manuscript

resfix <- read.table(paste0(path.sim,'fixresults_091618.txt'),header=T, stringsAsFactors=F)
fit <- c(0,1000,1001,1111,1)
fix.plot(dat=resfix,fit=fit)

###########################################################################################
###Sim results included in manuscript and supplementary materials
res1 <- read.table(paste0(path.sim,'simresults_set1_091618.txt'),header=T, stringsAsFactors=F)

##Dose-level 3 figures included in the manuscript for all other levels correctly specified (nic_dose3, p3 results) and all other dose-levels misspecified to varying degrees (nic_dose3_var, p3_var)
#Figure 2 in manuscript
nic_dose3 <- c('4/3/4/1/0','4/3/3.5/1/0','4/3/3/1/0','4/3/2.5/1/0','4/3/2/1/0','4/3/1.5/1/0','4/3/1/1/0','4/3/0.5/1/0','4/3/0/1/0')
p3 <- sim.plot(dat=res1,nic_scen=nic_dose3,doses=seq(0,4,.5),dose.ch=3,out='mse.c95',mainlab='95th Percentile MSE',meanlog2_use=4,sdlog_use=0.668,fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95','rjmcmc99'))

#calculuate percent change in MSE compared to model with IND
(p3$dose3 - matrix(rep(p3$dose3['allcomp',],each=6),nrow=6)) / matrix(rep(p3$dose3['allcomp',],each=6),nrow=6) 

#Figure 3 in manuscript
nic_dose3_var <- c('4/3/4/1/0','4/3/3.5/1/0','4/3/3/1/0','4/3/2.5/1/0','4/3/2/1/0','4/3/1.5/1/0','4/3/1/1/0','4/3/0.5/1/0','4/3/0/1/0')
p3_var <- sim.plot(dat=res1,nic_scen=nic_dose3_var,doses=seq(0,4,.5),dose.ch=3,out='mse.c95',mainlab='95th Percentile MSE',meanlog2_use=4,sdlog_use=0.668,fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95','rjmcmc99'))

##Dose-level 1, 2, and 4 figure results included in the Supplementary Materials as Figures 1-3
nic_dose1 <- c('6/3/2/1/0','5.5/3/2/1/0','5/3/2/1/0','4.5/3/2/1/0','4/3/2/1/0','3.5/3/2/1/0','3/3/2/1/0','2.5/3/2/1/0','2/3/2/1/0','1.5/3/2/1/0','1/3/2/1/0','0.5/3/2/1/0','0/3/2/1/0')
p1 <- sim.plot(dat=res1,nic_scen=nic_dose1_ex,doses=seq(0,6,.5),dose.ch=1,out='mse.c95',mainlab='95th Percentile MSE',meanlog2_use=4,sdlog_use=0.668,fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95','rjmcmc99'))

nic_dose2 <- c('4/5/2/1/0','4/4.5/2/1/0','4/4/2/1/0','4/3.5/2/1/0','4/3/2/1/0','4/2.5/2/1/0','4/2/2/1/0','4/1.5/2/1/0','4/1/2/1/0')
p2 <- sim.plot(dat=res1,nic_scen=nic_dose2,doses=seq(1,5,.5),dose.ch=2,out='mse.c95',mainlab='95th Percentile MSE',meanlog2_use=4,sdlog_use=0.668,fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95','rjmcmc99'))

nic_dose4 <- c('4/3/2/3/0','4/3/2/2.5/0','4/3/2/2/0','4/3/2/1.5/0','4/3/2/1/0','4/3/2/0.5/0','4/3/2/0/0')
p4 <- sim.plot(dat=res1,nic_scen=nic_dose4,doses=seq(0,3,.5),dose.ch=4,out='mse.c95',mainlab='95th Percentile MSE',meanlog2_use=4,sdlog_use=0.668,fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95','rjmcmc99'))

###########################################################################################
###Additional simulation studies for SM
##Used to create Tables 1-6 in the Supplementary Materials
res2 <- read.table(paste0(path.sim,'simresults_set2_091618.txt'),header=T, stringsAsFactors=F)

#Take simulation results and format for presentation, use $resq for Supplementary Materials
s1 <- sim.tab_dist(dat=res2,sdlog_use=0.668,meanlog2_use=4,nic_true='4/3/2/1/0',fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95'))
s2 <- sim.tab_dist(dat=res2,sdlog_use=0.668,meanlog2_use=4,nic_true='4/3/2/0.5/0',fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95'))
s3 <- sim.tab_dist(dat=res2,sdlog_use=0.668,meanlog2_use=4,nic_true='4/3/1/0.5/0',fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95'))
s4 <- sim.tab_dist(dat=res2,sdlog_use=0.668,meanlog2_use=4,nic_true='4/1.5/1/0.5/0',fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95'))
s5 <- sim.tab_dist(dat=res2,sdlog_use=0.668,meanlog2_use=4,nic_true='2/1.5/1/0.5/0',fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95'))
s6 <- sim.tab_dist(dat=res2,sdlog_use=0.668,meanlog2_use=4,nic_true='6/4.5/3/1.5/0',fituse=c('allcomp','meanrel','rjmcmc05','rjmcmc50','rjmcmc95'))

