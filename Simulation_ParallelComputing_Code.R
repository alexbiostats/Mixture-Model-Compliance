###Code to run simulation scenarios with parallel computing via snowfall package
asd <- 5 #number of dose-levels for simulations

library(snowfall)
sfInit(cpus = 12, type = 'SOCK', parallel = T)

sfExportAll()

sfLibrary(gtools)
sfLibrary(matrixStats)

sfSource("./GitHub_Functions_ComplianceMixture.R")

##################################################################################################################################
##################################################################################################################################
###Fixed simulation scenarios (results produce Figure 1 in Manuscript)
write.table(matrix(c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('norm_mu_comp',1:asd),paste0('norm_mu_nonc',1:asd)), nrow=1), file = 'fixresults_091618.txt', col.name=FALSE, row.name=FALSE, append=T)
sfSapply(1:500, function(x) fix_sim(x,nic=c(6,4,3,1,0),es=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate='091618') )

##################################################################################################################################
##################################################################################################################################
###Set 1 simulation scenarios (results produce Figures 2 and 3 in Manuscript, and Figures 1-3 in Supplementary Materials)
write.table(matrix(c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1))), nrow=1), file = 'simresults_set1_091618.txt', col.name=FALSE, row.name=FALSE, append=T)
simdate <- 'set1_091618'

#Dose-level 3 varies, other dose-levels correctly specified for REL, results used for Figure 2 in manuscript
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,4,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,3.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,3,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,1.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,1,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,0.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,0,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )

#Dose-level 3 varies, other dose-levels misspecified to varying degrees (i.e., nic argument), results used for Figure 3 in manuscript
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,4,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,3.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,3,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,2.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,1.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,1,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,0.5,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(6,2,2,0.1,0),es=c(4,3,0,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )

#Dose-level 1 varies, over extended range, results used for Figure 1 in Supplementary Materials
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(6,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(5.5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4.5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(3.5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(3,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(2.5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(2,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(1.5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(1,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(0.5,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(0,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )

#Dose-level 2 varies, results used for Figure 2 in Supplementary Materials
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,5,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,4.5,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,4,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3.5,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,2.5,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,2,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,1.5,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,1,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )

#Dose-level 4 varies, results used for Figure 3 in Supplementary Materials
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,3,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,2.5,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,2,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,1.5,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,0.5,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )
sfSapply(1:500, function(x) mix_sim(x,nic=c(4,3,2,1,0),es=c(4,3,2,0,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate) )


##################################################################################################################################
##################################################################################################################################
###Set 2 simulation scenarios (Tables 1-6 in Supplementary Materials)
write.table(matrix(c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1))), nrow=1), file = 'simresults_set2_091618.txt', col.name=FALSE, row.name=FALSE, append=T)
simdate2 <- 'set2_091618'

sfSapply(1:1000, function(x) mix_sim(x,es=c(4,3,2,1,0),nic=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate2) )
sfSapply(1:1000, function(x) mix_sim(x,es=c(4,3,2,0.5,0),nic=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate2) )
sfSapply(1:1000, function(x) mix_sim(x,es=c(4,3,1,0.5,0),nic=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate2) )
sfSapply(1:1000, function(x) mix_sim(x,es=c(4,1.5,1,0.5,0),nic=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate2) )
sfSapply(1:1000, function(x) mix_sim(x,es=c(2,1.5,1,0.5,0),nic=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate2) )
sfSapply(1:1000, function(x) mix_sim(x,es=c(6,4.5,3,1.5,0),nic=c(4,3,2,1,0),meanlog2=4,sigma=.668,pcomp=c(rep(.3,4),1),nsamp=rep(100,5),nburn=1000,nchain=10000, simdate=simdate2) )


sfStop()