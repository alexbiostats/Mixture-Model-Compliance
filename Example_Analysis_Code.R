#########################################################################
### Example Analysis Code
### Simulate an example data set of data and demonstrate different functions
#########################################################################

library(xtable)

source('~/GitHub_Functions_ComplianceMixture.R')

#####################################################################
### Specify settings for calculations (e.g., chain length, burn-in period, initial values, prior values)

plot.ind <- F #set equal to T if plots wanted with output
export.ind <- F #set equal to T if you want to export results to wd
nburn <- 100 #burn-in period to discard
nchain <- 1000 #chain length
seed <- 515 #specify seed for functions for reproducibility
hp <- list( tau.a=.001, tau.b=.001, theta.prec=0.001, mu.prec=0.001, p.alpha=1, p.beta=1 ) #specify priors
inits <- list( theta=10, tau=1, p=0.5, mu_comp=10, mu_nonc=20, model=0) #specify initial values
prob <- rep(.95,5) #prior probability for RJ assuming level j should be estimated with IND vs. REL



#####################################################################
### Simulate data for example data set to use in examples of implementing functions
# Note: CENIC-p1s1 example assumed log-normally distributed biomarker, so posterior estimates have been provided with the translation back to the non-log scale

set.seed(515) 

nic <- c(5,10,15,25,50) #specify some sort of biomarker relationship, note here dose-levels 1 and 2 (i.e. 5 and 10) are correct, but 3 and 4 are misspecified

v1 <- log(sort(c(rnorm(n=10,mean=5,sd=1), rnorm(n=40,mean=50,sd=1))))
v2 <- log(sort(c(rnorm(n=10,mean=10,sd=1), rnorm(n=40,mean=50,sd=1))))
v3 <- log(sort(c(rnorm(n=10,mean=13,sd=1), rnorm(n=40,mean=50,sd=1))))
v4 <- log(sort(c(rnorm(n=10,mean=16.4,sd=1), rnorm(n=40,mean=50,sd=1))))
v5 <- log(sort(rnorm(n=50,mean=50,sd=1)))

vals <- c(v1,v2,v3,v4,v5)
grp <- c(rep(1,length(v1)),rep(2,length(v2)),rep(3,length(v3)),rep(4,length(v4)),rep(5,length(v5)))
grp.lab <- paste0('d',1:5) #give names to dose-levels/groups to fit mixtures to
grp.title <- paste0('d',1:5) #give names to dose-levels/groups to fit mixtures to



#########################################################################
###Fit simulated data with IND, REL, RJ, and "fixed" models without diagnostics
#This section fits the IND, REL, RJ, and fixed models without diagnostics on 1 chain

ind_model <- mcmc.gaussianmix(vals=vals,grp=grp,fit='ind',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,inits=inits,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind)	
rel_model <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rel',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,inits=inits,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind)
rj_model <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rjmcmc',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,inits=inits,prob=prob,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind)
fix_model <- mcmc.gaussianmix(vals=vals,grp=grp,fit='fixed',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,prob=prob,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind,inits=list( theta=10, tau=1, p=0.5, mu_comp=10, mu_nonc=20, model=c(1,0,0,1,1)))

#Extract tab.hpd object which returns parameter posterior estimates and 95% HPD intervals for each dose-level
ind_tab <- ind_model$tab.hpd
rel_tab <- rel_model$tab.hpd
rj_tab <- rj_model$tab.hpd
fix_tab <- fix_model$tab.hpd

#Combine posterior estimates from IND, REL, and RJ models by dose-level/group into one table
mat.res <- NULL
for( i in 1:5 ){
	mat.res <- rbind(mat.res, ind_tab[i,], rel_tab[i,], rj_tab[i,])
}

print( xtable(mat.res, digits=c(0,3,3,3,3,3,3,3,3)) ) #print as LaTeX table

#Plot compliance probability figure for IND and RJ
compliance_fig(tab1=ind_model$log.tab, tab2=rj_model$log.tab, export=F, xmax=80, bio.lab=c(5,10,15,25))



#########################################################################
###Diagnostics included in application to CENIC-p1 data
#This section fits the above IND, REL, and RJ models, but with diagnostics calculated and it uses 3 chains instead of 1

inits_diag <- list( list( tau=1, p=0.5, mu_comp=0, mu_nonc=1, model=c(0,0,0,0,1)),
	list( tau=1, p=0.5, mu_comp=-10, mu_nonc=10, model=c(1,1,1,1,1)),
	list( tau=1, p=0.5, mu_comp=10, mu_nonc=20, model=c(1,0,1,0,1)) )

ind_model_diag <- mcmc.gaussianmix_diag(vals=vals,grp=grp,fit='ind',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,inits=inits_diag,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind,diag=T)	
rel_model_diag <- mcmc.gaussianmix_diag(vals=vals,grp=grp,fit='rel',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,inits=inits_diag,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind,diag=T)
rj_model_diag <- mcmc.gaussianmix_diag(vals=vals,grp=grp,fit='rjmcmc',nburn=nburn,nchain=nchain,plot.ind=plot.ind,plot.grp=grp.title,hp=hp,inits=inits_diag,prob=prob,seed=seed,nic=nic,grp.lab=grp.lab,export=export.ind,diag=T)

#Gelman-Rubin diagnostic for IND and REL
ind_model_diag$gr_test
rel_model_diag$gr_test

ind_tab_diag <- ind_model_diag$tab.hpd
rel_tab_diag <- rel_model_diag$tab.hpd
rj_tab_diag <- rj_model_diag$tab.hpd

mat.res_diag <- NULL
for( i in 1:5 ){
	mat.res_diag <- rbind(mat.res_diag, ind_tab_diag[i,], rel_tab_diag[i,], rj_tab_diag[i,])
}

#Combine posterior estimates from IND, REL, and RJ models by dose-level/group into one table
print( xtable(mat.res_diag, digits=c(0,3,3,3,3,3,3,3,3)) )
print( xtable(mat.res_diag[,c(1:3,6:8)], digits=c(0,3,3,3,3,3,3)) ) #select specific parameters

#Plot compliance probability figure for IND and RJ
compliance_fig(tab1=ind_model_diag$log.tab, tab2=rj_model_diag$log.tab, export=F, xmax=80, bio.lab=c(5,10,15,25))



