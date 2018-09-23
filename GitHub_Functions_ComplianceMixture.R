library(matrixStats)

####################################################################
## Application for MCMC for Mixture Distribution with Diagnostics ##
####################################################################

mcmc.gaussianmix_diag <- function(vals,fit,nburn,nchain,plot.ind=F,plot.grp='',hp,hp.rj=NULL,inits,prob=NULL,nic=NULL,grp,notmix=T,seed,export=F,grp.lab='',diag=T){
###Function to calculate mixture components with diagnostics from 3 chains
#vals: vector of measurements to fit mixture model to
#fit: model to fit:
##ind: estimate all components in Gaussian mixture
##rel: mean values for each Gaussian component are given with relationship assumed amongst compliant means between groups
##rjmcmc: use RJMCMC between the two models
#nburn: burn-in period with default of 1000 to discard from chain
#nchain: chain length with default of 3000
#plot.ind: histogram of data with fitted model densities plotted if TRUE
#plot.grp: group to add to plot legend
#inits: list with 3 lists of different initial values to use for lambda, theta, tau, p [note: tau is precision]
#hp: list of hyperparameter values (for lambda, theta, tau, and p) 
#hp.rj: list of hyperparameter values for parameters added/removed in RJ step (for lambda and theta) 
#prob: vector with prior model probability for model being correct (prob for REL for each level/group)
#nic: the nicotine group having its mixture distribution estimated as a nicotine value: e.g., 0.4, 1.3, 2.4, 5.2, 15.8
#grp: a vector identifying which group each data belongs to 
#notmix: indicator if last group shouldn't estimate entire mixture since all are assumed compliant, default is TRUE
#seed: seed to set samplers
#export: write plots to PDF if plot=T
#grp.lab: label to include in PDF title to identify group
#diag: if TRUE, return diagnostic plots/meausures

	if(fit %in% c('ind','rel','rjmcmc')){}else{print('The value for fit is not valid.'); break} #break out of function if invalid 'fit' type entered
	j <- length(unique(grp))
	yj <- sapply( split(vals, grp), sort)

	#Fit sampler/algorithm
	if(fit=='ind'){
		res1 <- mix2normals(y=vals,grp=grp,M=nchain,nburn=0,hp=hp,inits=inits[[1]],seed=seed,notmix=notmix,type='ind')
		res2 <- mix2normals(y=vals,grp=grp,M=nchain,nburn=0,hp=hp,inits=inits[[2]],seed=seed,notmix=notmix,type='ind')
		res3 <- mix2normals(y=vals,grp=grp,M=nchain,nburn=0,hp=hp,inits=inits[[3]],seed=seed,notmix=notmix,type='ind')
	}
	if(fit=='rel'){
		res1 <- mix2normals(y=vals,grp=grp,M=nchain,nburn=0,hp=hp,inits=inits[[1]],seed=seed,type='rel',nic=nic,notmix=T)
		res2 <- mix2normals(y=vals,grp=grp,M=nchain,nburn=0,hp=hp,inits=inits[[2]],seed=seed,type='rel',nic=nic,notmix=T)
		res3 <- mix2normals(y=vals,grp=grp,M=nchain,nburn=0,hp=hp,inits=inits[[3]],seed=seed,type='rel',nic=nic,notmix=T)
	}
	if(fit=='rjmcmc'){
		res1 <- rjmcmc_2normals(y=vals,grp=grp,inits=inits[[1]],hp=hp,hp.rj=hp.rj,prob=prob,M=nchain,seed=seed,nburn=0,nic=nic)
		res2 <- rjmcmc_2normals(y=vals,grp=grp,inits=inits[[2]],hp=hp,hp.rj=hp.rj,prob=prob,M=nchain,seed=seed,nburn=0,nic=nic)
		res3 <- rjmcmc_2normals(y=vals,grp=grp,inits=inits[[3]],hp=hp,hp.rj=hp.rj,prob=prob,M=nchain,seed=seed,nburn=0,nic=nic)
	}

	#Diagnostics
	samp1 <- res1$samp; samp2 <- res2$samp; samp3 <- res3$samp
	gr_test <- NULL #intialize

	if((fit %in% c('ind','rel')) & diag==T){

		#Gelman-Rubin calculations/plots
		gr_test <- sapply(1:ncol(samp1), function(z) gelman.rubin( out = list(samp1[,z], samp2[,z], samp3[,z]) ))

		#autocorrelation calculations/plots
		auto_test <- lapply(1:ncol(samp1), function(z) autocor(out = list(samp1[,z], samp2[,z], samp3[,z]),h.max=50,nburn=nburn))
		names(auto_test) <- colnames(samp1)
		for(d in 1:ncol(samp1)){
			if(d %in% seq(1,ncol(samp1),by=9)){dev.new();par(mfrow=c(3,3))}
			plot(auto_test[[d]][[1]], ylim=c(-1,1), main=names(auto_test)[d], type='l', ylab='Autocorrelation')
			lines(auto_test[[d]][[2]], col='blue')
			lines(auto_test[[d]][[3]], col='orangered2')		
		}

		#trace plots
		for(d in 1:ncol(samp1)){
			if(d %in% seq(1,ncol(samp1),by=9)){dev.new();par(mfrow=c(3,3))}
			plot(samp1[,d], main=names(auto_test)[d], type='l', ylab=names(auto_test)[d], ylim=range(c(samp1[,d],samp2[,d],samp3[,d])) )
			lines(samp2[,d], col='blue', lty=2)
			lines(samp3[,d], col='orangered2',lty=3)		
			abline( v=nburn, col='gray65' )
		}
	}
	if((fit %in% 'rjmcmc') & diag==T){
		s1nb <- samp1[-(1:nburn),c('model','dev')]
		s2nb <- samp2[-(1:nburn),c('model','dev')]
		s3nb <- samp3[-(1:nburn),c('model','dev')]

		iter.seq <- seq(100,nrow(s1nb),by=10)
		bg_test <- t(sapply(iter.seq, function(x) brooks.giudici( out=list(s1nb[1:x,], s2nb[1:x,], s3nb[1:x,]) ) ))

		dev.new()
		par(mfrow=c(2,2))

		plot( y=bg_test[,'Vhat'], x=iter.seq, xlab='Iteration' , type='l', ylim=range(bg_test[,c('Vhat','Wc')]) )
		lines( x=iter.seq, y=bg_test[,'Wc'], lty=2)

		plot( y=bg_test[,'Wm'], x=iter.seq, xlab='Iteration' , type='l', ylim=range(bg_test[,c('Wm','WmWc')]) )
		lines( x=iter.seq, y=bg_test[,'WmWc'], lty=2)

		plot( y=bg_test[,'Bm'], x=iter.seq, xlab='Iteration' , type='l', ylim=range(bg_test[,c('Bm','BmWc')]) )
		lines( x=iter.seq, y=bg_test[,'BmWc'], lty=2)

	}

	#Extract mixture results, remove burn-in period
	samp.all <- rbind( res1$samp[-(1:nburn),], res2$samp[-(1:nburn),], res3$samp[-(1:nburn),] )
				
	#Extract mixture results; for log-normal transformed data, calculate estimates at each iteration to use
	ps <- samp.all[, paste0('p', 1:j) ]
	mucomps <- samp.all[, paste0('mu_comp', 1:j) ]
	munoncs <- samp.all[, paste0('mu_nonc', 1:j) ]
	sds <- samp.all[, paste0('sigma', 1:j) ]
	if(fit=='rjmcmc'){ model <- samp.all[,'model'] }
	if(fit=='rjmcmc'){ dev <- samp.all[,'dev'] }

	#For use in plotting compliance figure, estimates based on log(Y) ~ Normal
	ps.log <- 1-colMeans(ps)
	mc.log <- colMeans(mucomps)
	mn.log <- colMeans(munoncs)
	sd.log <- colMeans(sds)
	log.tab <- cbind(ps.log, mc.log, mn.log, sd.log); colnames(log.tab) <- c('p_comp','mean.mu_comp','mean.mu_nonc','sd')

	ln_mucomps <- exp( mucomps + 0.5*sds^2 ); colnames(ln_mucomps) <- paste0('mu_comp',1:j)
	ln_munoncs <- exp( munoncs + 0.5*sds^2 ); colnames(ln_munoncs) <- paste0('mu_nonc',1:j)
	ln_sd1 <- sqrt( exp(2*mucomps + sds^2) * (exp(sds^2)-1) ); colnames(ln_sd1) <- paste0('sigma1',1:j)
	ln_sd2 <- sqrt( exp(2*munoncs + sds^2) * (exp(sds^2)-1) ); colnames(ln_sd2) <- paste0('sigma2',1:j)
	ln_c80 <- qlnorm(meanlog=mucomps, sdlog=sds, 0.8); colnames(ln_c80) <- paste0('c80_',1:j)
	ln_c90 <- qlnorm(meanlog=mucomps, sdlog=sds, 0.9); colnames(ln_c90) <- paste0('c90_',1:j)
	ln_c95 <- qlnorm(meanlog=mucomps, sdlog=sds, 0.95); colnames(ln_c95) <- paste0('c95_',1:j)

	samp <- cbind( ps, ln_mucomps, ln_munoncs, ln_sd1, ln_sd2, ln_c80, ln_c90, ln_c95)
	if(fit=='rjmcmc'){ samp <- cbind( samp, model, dev) }
		
	#Extract column means and standard deviations
	samp.m <- colMeans(samp)
	samp.s <- colSds(samp); names(samp.s) <- names(samp.m) #add names to SD vector

	#Extract estimate proportion for models (note that we should have p1=1-p2 and p2=1-p1)
	p_comp <- 1-samp.m[ paste0('p',1:j) ] #proportion in complier group
	p_nonc <- samp.m[ paste0('p',1:j) ] #proportion in non-complier group
	hpd.p_comp <- sapply(1:j, function(x) boa.hpd(1-samp[ , paste0('p',x) ], alpha=0.05))

	#Compliant group mean
	mean.mu_comp <- samp.m[ paste0('mu_comp', 1:j) ] 
	sd.mu_comp <- samp.s[ paste0('mu_comp', 1:j) ] 
	prec.mu_comp <- 1/sd.mu_comp
	hpd.mu_comp <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('mu_comp',x) ], alpha=0.05))

	#Compliant group percentiles
	mean.c80 <- samp.m[ paste0('c80_', 1:j) ]
	mean.c90 <- samp.m[ paste0('c90_', 1:j) ]
	mean.c95 <- samp.m[ paste0('c95_', 1:j) ]
	hpd.c90 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('c90_',x) ], alpha=0.05))
	hpd.c95 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('c95_',x) ], alpha=0.05))

	#Non-compliant group mean
	mean.mu_nonc <- samp.m[ paste0('mu_nonc', 1:j) ] 
	sd.mu_nonc <- samp.s[ paste0('mu_nonc', 1:j) ] 
	prec.mu_nonc <- 1/sd.mu_nonc
	hpd.mu_nonc <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('mu_nonc',x) ], alpha=0.05))

	#Standard deviation for Gaussian components
	sd1 <- samp.m[ paste0('sigma1', 1:j) ]
	hpd.sd1 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('sigma1',x) ], alpha=0.05))
	sd2 <- samp.m[ paste0('sigma2', 1:j) ]
	hpd.sd2 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('sigma2',x) ], alpha=0.05))

	if(fit=='rjmcmc'){
		mods <- samp[,'model']
		mods_avg <- colMeans(t(sapply(1:length(mods), function(x) base10toK(mods[x], K=2, size=(j-1) ))))
		if( j == 2 ){ mods_avg <- mean(mods) }
		mods_tab <- table(mods)
	}else{
		mods_avg <- rep(NA, (j-1))
		mods_tab <- NULL				
	}

	tab <- cbind( p_comp, mean.mu_comp, mean.mu_nonc, sd1, sd2, mean.c80, mean.c90, mean.c95 )
	
	h1 <- sapply(1:j, function(x) paste0(round2(p_comp[x]),' (',round2(hpd.p_comp[1,x]),',',round2(hpd.p_comp[2,x]),')') )
	h2 <- sapply(1:j, function(x) paste0(round2(mean.mu_comp[x]),' (',round2(hpd.mu_comp[1,x]),',',round2(hpd.mu_comp[2,x]),')') )
	h3 <- sapply(1:j, function(x) paste0(round2(mean.mu_nonc[x]),' (',round2(hpd.mu_nonc[1,x]),',',round2(hpd.mu_nonc[2,x]),')') )
	h4 <- sapply(1:j, function(x) paste0(round2(sd1[x]),' (',round2(hpd.sd1[1,x]),',',round2(hpd.sd1[2,x]),')') )
	h5 <- sapply(1:j, function(x) paste0(round2(sd2[x]),' (',round2(hpd.sd2[1,x]),',',round2(hpd.sd2[2,x]),')') )
	h6 <- sapply(1:j, function(x) paste0(round2(mean.c90[x]),' (',round2(hpd.c90[1,x]),',',round2(hpd.c90[2,x]),')') )
	h7 <- sapply(1:j, function(x) paste0(round2(mean.c95[x]),' (',round2(hpd.c95[1,x]),',',round2(hpd.c95[2,x]),')') )

	tab.hpd <- cbind(h1,h2,h3,h4,h5,h6,h7,c(mods_avg,NA) )
	colnames(tab.hpd) <- c('p_c','mu_comp','mu_nonc','sd1','sd2','c90','c95','propREL')
	rownames(tab.hpd) <- rownames(tab)

	if(plot.ind==T){
		mdy <- format(Sys.Date(), format='%m%d%y')
		yj.range <- c(floor(min(yj)), ceiling(max(yj)))
		breaks <- seq(yj.range[1],yj.range[2],by=.1)

		###Create histograms overlaid with estimated compliant and non-compliant densities
		for(a in 1:j){
			if(export==T){ pdf(paste0('hist_',fit,'_',grp.lab[a],'_',mdy,'.pdf'))}else{ dev.new() }

			#Define plot title
			if(fit=='ind'){plot.title <- paste0('IND Model, ',plot.grp[a])}
			if(fit=='rel'){plot.title <- paste0('REL Model, ',plot.grp[a])}
			if(fit=='rjmcmc'){plot.title <- paste0('RJMCMC, ',plot.grp[a])}

			h<-hist(yj[,a], freq=F, xlab='log(biomarker)', breaks=breaks, xlim=yj.range,main=plot.title, cex.axis=1.4, cex.main=1.4, cex.lab=1.4)
			legend('topleft', lty=c(2,2,1), col=c('orangered2','blue','green'), legend=c('Complier Component','Non-complier Component','Mixture Density'), bty='n', cex=1.2, lwd=2)

			xfit<-seq(-2,6,length.out=1000)
			yfit1<- dnorm(xfit,mean=mc.log[a],sd=sd.log[a]) * p_comp[a]
			lines(xfit, yfit1, col='orangered2', lty=2, lwd=2)

			yfit2 <- dnorm(xfit, mean=mn.log[a], sd=sd.log[a]) * p_nonc[a]
			lines(xfit, yfit2, col='blue', lty=2, lwd=2)

			yfit.mix <- yfit1 + yfit2
			lines(xfit, yfit.mix, col='green', lty=1, lwd=2)
			if(export==T){ dev.off() }
		}
		
	}


	##################
	#Calculate compliance-related values to return
	z1.prop <- NA #lapply(1:j, function(x) res$zj.sum[[x]] / (nchain-nburn) )

	ret <- list(p_comp=p_comp,p_nonc=p_nonc,mu1=mean.mu_comp,mu2=mean.mu_nonc,sd1=sd1,sd2=sd2,mu1_80=mean.c80,mu1_90=mean.c90,mu1_95=mean.c95,tab=tab,mods_avg=mods_avg,mods_tab=mods_tab,z1.prop=z1.prop,yj=yj,tab.hpd=tab.hpd,gr_test=gr_test,log.tab=log.tab)	
	return(ret)

}

###################################
## MCMC for Mixture Distribution ##
###################################

mcmc.gaussianmix <- function(vals,fit,nburn,nchain,plot.ind=F,plot.grp='',hp,hp.rj=NULL,inits,prob=NULL,nic=NULL,grp,notmix=T,seed,export=F,grp.lab='',pair=F){
###Function to run self-programmed MCMC for Gaussian mixture
#vals: vector of measurements to fit mixture model to
#fit: model to fit:
##ind: estimate all components in Gaussian mixture
##rel: mean values for each Gaussian component are given with relationship assumed amongst compliant means between groups
##rjmcmc: use RJMCMC between the two models
##fixed: fit dose-levels by ind or rel as specified (can be different for dose-levels, but doesn't implement RJ)
#nburn: burn-in period with default of 1000 to discard from chain
#nchain: chain length with default of 3000
#plot.ind: histogram of data with fitted model densities plotted if TRUE
#plot.grp: group to add to plot legend
#inits: list of initial values to use for lambda, theta, tau, p [note: tau is precision]
#hp: list of hyperparameter values (for lambda, theta, tau, and p) 
#hp.rj: list of hyperparameter values for parameters added/removed in RJ step (for lambda and theta) 
#prob: vector with prior model probability for model being correct (prob for REL for each level/group)
#nic: the nicotine group having its mixture distribution estimated as a nicotine value: e.g., 0.4, 1.3, 2.4, 5.2, 15.8
#grp: a vector identifying which group each data belongs to (currently 'c' vs. 'ub')
#notmix: indicator if last group shouldn't estimate entire mixture since all are assumed compliant, default is TRUE
#seed: seed to set samplers
#export: write plots to PDF if plot=T
#grp.lab: label to include in PDF title to identify group
#pair: if TRUE, return pairwise matrices for concordance of compliance and non-compliance	

	if(fit %in% c('ind','rel','rjmcmc','fixed')){}else{print('The value for fit is not valid.'); break} #break out of function if invalid 'fit' type entered
	j <- length(unique(grp))
	yj <- sapply( split(vals, grp), sort)

	#Fit sampler/algorithm
	if(fit=='ind'){res <- mix2normals(y=vals,grp=grp,M=nchain,nburn=nburn,hp=hp,inits=inits,seed=seed,notmix=notmix,type='ind',pair=pair)}
	if(fit=='rel'){res <- mix2normals(y=vals,grp=grp,M=nchain,nburn=nburn,hp=hp,inits=inits,seed=seed,type='rel',nic=nic,notmix=T,pair=pair)}
	if(fit=='rjmcmc'){ res <- rjmcmc_2normals(y=vals,grp=grp,inits=inits,hp=hp,hp.rj=hp.rj,prob=prob,M=nchain,seed=seed,nburn=nburn,nic=nic,pair=pair)}
	if(fit=='fixed'){  res <- setmcmc_2normals(y=vals,grp=grp,inits=inits,hp=hp,prob=prob,M=nchain,seed=seed,nburn=nburn,nic=nic,pair=pair)}

	#Extract mixture results; for log-normal transformed data, calculate estimates at each iteration to use
	ps <- res$samp[, paste0('p', 1:j) ]
	mucomps <- res$samp[, paste0('mu_comp', 1:j) ]
	munoncs <- res$samp[, paste0('mu_nonc', 1:j) ]
	sds <- res$samp[, paste0('sigma', 1:j) ]
	if(fit %in% c('rjmcmc','fixed')){ model <- res$samp[,'model'] }
	if(fit %in% c('rjmcmc','fixed')){ dev <- res$samp[,'dev'] }

	#For use in plotting compliance figure, estimates based on log(Y) ~ Normal
	ps.log <- 1-colMeans(ps)
	mc.log <- colMeans(mucomps)
	mn.log <- colMeans(munoncs)
	sd.log <- colMeans(sds)
	log.tab <- cbind(ps.log, mc.log, mn.log, sd.log); colnames(log.tab) <- c('p_comp','mean.mu_comp','mean.mu_nonc','sd')

	#Estimates on original Y ~ LN scale
	ln_mucomps <- exp( mucomps + 0.5*sds^2 ); colnames(ln_mucomps) <- paste0('mu_comp',1:j)
	ln_munoncs <- exp( munoncs + 0.5*sds^2 ); colnames(ln_munoncs) <- paste0('mu_nonc',1:j)
	ln_sd1 <- sqrt( exp(2*mucomps + sds^2) * (exp(sds^2)-1) ); colnames(ln_sd1) <- paste0('sigma1',1:j)
	ln_sd2 <- sqrt( exp(2*munoncs + sds^2) * (exp(sds^2)-1) ); colnames(ln_sd2) <- paste0('sigma2',1:j)
	ln_c80 <- qlnorm(meanlog=mucomps, sdlog=sds, 0.8); colnames(ln_c80) <- paste0('c80_',1:j)
	ln_c90 <- qlnorm(meanlog=mucomps, sdlog=sds, 0.9); colnames(ln_c90) <- paste0('c90_',1:j)
	ln_c95 <- qlnorm(meanlog=mucomps, sdlog=sds, 0.95); colnames(ln_c95) <- paste0('c95_',1:j)

	samp <- cbind( ps, ln_mucomps, ln_munoncs, ln_sd1, ln_sd2, ln_c80, ln_c90, ln_c95)
	if(fit %in% c('rjmcmc','fixed')){ samp <- cbind( samp, model, dev) }
		
	#Extract column means and standard deviations
	samp.m <- colMeans(samp)
	samp.s <- colSds(samp); names(samp.s) <- names(samp.m) #add names to SD vector

	mucomps.m <- colMeans(mucomps)
	munoncs.m <- colMeans(munoncs)

	#Extract estimate proportion for models (note that we should have p1=1-p2 and p2=1-p1)
	p_comp <- 1-samp.m[ paste0('p',1:j) ] #proportion in complier group
	p_nonc <- samp.m[ paste0('p',1:j) ] #proportion in non-complier group
	hpd.p_comp <- sapply(1:j, function(x) boa.hpd(1-samp[ , paste0('p',x) ], alpha=0.05))

	#Compliant group mean
	mean.mu_comp <- samp.m[ paste0('mu_comp', 1:j) ] 
	sd.mu_comp <- samp.s[ paste0('mu_comp', 1:j) ] 
	prec.mu_comp <- 1/sd.mu_comp
	hpd.mu_comp <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('mu_comp',x) ], alpha=0.05))

	#Compliant group percentiles
	mean.c80 <- samp.m[ paste0('c80_', 1:j) ]
	mean.c90 <- samp.m[ paste0('c90_', 1:j) ]
	mean.c95 <- samp.m[ paste0('c95_', 1:j) ]
	hpd.c90 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('c90_',x) ], alpha=0.05))
	hpd.c95 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('c95_',x) ], alpha=0.05))

	#Non-compliant group mean
	mean.mu_nonc <- samp.m[ paste0('mu_nonc', 1:j) ] 
	sd.mu_nonc <- samp.s[ paste0('mu_nonc', 1:j) ] 
	prec.mu_nonc <- 1/sd.mu_nonc
	hpd.mu_nonc <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('mu_nonc',x) ], alpha=0.05))

	#Standard deviation for Gaussian components
	sd1 <- samp.m[ paste0('sigma1', 1:j) ]
	hpd.sd1 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('sigma1',x) ], alpha=0.05))
	sd2 <- samp.m[ paste0('sigma2', 1:j) ]
	hpd.sd2 <- sapply(1:j, function(x) boa.hpd(samp[ , paste0('sigma2',x) ], alpha=0.05))

	if(fit %in% c('rjmcmc','fixed') ){
		mods <- samp[,'model']
		mods_avg <- colMeans(t(sapply(1:length(mods), function(x) base10toK(mods[x], K=2, size=(j-1) ))))
		if( j == 2 ){ mods_avg <- mean(mods) }
		mods_tab <- table(mods)
	}else{
		mods_avg <- rep(NA, (j-1))
		mods_tab <- NULL				
	}

	tab <- cbind( p_comp, mean.mu_comp, mean.mu_nonc, sd1, sd2, mean.c80, mean.c90, mean.c95 )
	
	h1 <- sapply(1:j, function(x) paste0(round2(p_comp[x]),' (',round2(hpd.p_comp[1,x]),',',round2(hpd.p_comp[2,x]),')') )
	h2 <- sapply(1:j, function(x) paste0(round2(mean.mu_comp[x]),' (',round2(hpd.mu_comp[1,x]),',',round2(hpd.mu_comp[2,x]),')') )
	h3 <- sapply(1:j, function(x) paste0(round2(mean.mu_nonc[x]),' (',round2(hpd.mu_nonc[1,x]),',',round2(hpd.mu_nonc[2,x]),')') )
	h4 <- sapply(1:j, function(x) paste0(round2(sd1[x]),' (',round2(hpd.sd1[1,x]),',',round2(hpd.sd1[2,x]),')') )
	h5 <- sapply(1:j, function(x) paste0(round2(sd2[x]),' (',round2(hpd.sd2[1,x]),',',round2(hpd.sd2[2,x]),')') )
	h6 <- sapply(1:j, function(x) paste0(round2(mean.c90[x]),' (',round2(hpd.c90[1,x]),',',round2(hpd.c90[2,x]),')') )
	h7 <- sapply(1:j, function(x) paste0(round2(mean.c95[x]),' (',round2(hpd.c95[1,x]),',',round2(hpd.c95[2,x]),')') )

	tab.hpd <- cbind(h1,h2,h3,h4,h5,h6,h7,c(mods_avg,NA) )
	colnames(tab.hpd) <- c('p_c','mu_comp','mu_nonc','sd1','sd2','c90','c95','propREL')
	rownames(tab.hpd) <- rownames(tab)

	if(plot.ind==T){
		mdy <- format(Sys.Date(), format='%m%d%y')
		yj.range <- c(floor(min(yj)), ceiling(max(yj)))
		breaks <- seq(yj.range[1],yj.range[2],by=.1)

		###Create histograms overlaid with estimated compliant and non-compliant densities
		for(a in 1:j){
			if(export==T){ pdf(paste0('hist_',fit,'_',grp.lab[a],'_',mdy,'.pdf'))}else{ dev.new() }

			#Define plot title
			if(fit=='ind'){plot.title <- paste0('IND Model, ',plot.grp[a])}
			if(fit=='rel'){plot.title <- paste0('REL Model, ',plot.grp[a])}
			if(fit=='rjmcmc'){plot.title <- paste0('RJMCMC, ',plot.grp[a])}

			h<-hist(yj[,a], freq=F, xlab='log(biomarker)', breaks=breaks, xlim=yj.range,main=plot.title, cex.axis=1.4, cex.main=1.4, cex.lab=1.4)
			legend('topleft', lty=c(2,2,1), col=c('orangered2','blue','green'), legend=c('Complier Component','Non-complier Component','Mixture Density'), bty='n', cex=1.2, lwd=2)

			xfit<-seq(-2,6,length.out=1000)
			yfit1<- dnorm(xfit,mean=mc.log[a],sd=sd.log[a]) * p_comp[a]
			lines(xfit, yfit1, col='orangered2', lty=2, lwd=2)

			yfit2 <- dnorm(xfit, mean=mn.log[a], sd=sd.log[a]) * p_nonc[a]
			lines(xfit, yfit2, col='blue', lty=2, lwd=2)

			yfit.mix <- yfit1 + yfit2
			lines(xfit, yfit.mix, col='green', lty=1, lwd=2)
			if(export==T){ dev.off() }
		}
		
	}


	##################
	#Calculate compliance-related values to return
	z1.prop <- lapply(1:j, function(x) res$zj.sum[[x]] / (nchain-nburn) )

	if(pair==T){
		comp.mat <- res$comp.mat/(nchain-nburn) #lapply(1:j, function(x) res$comp.mat[[x]] / (nchain-nburn) )		
		nonc.mat <- res$nonc.mat/(nchain-nburn) #lapply(1:j, function(x) res$nonc.mat[[x]] / (nchain-nburn) )		
	}else{
		comp.mat <- nonc.mat <- NA
	}

	ret <- list(p_comp=p_comp,p_nonc=p_nonc,mu1=mean.mu_comp,mu2=mean.mu_nonc,sd1=sd1,sd2=sd2,mu1_80=mean.c80,mu1_90=mean.c90,mu1_95=mean.c95,tab=tab,mods_avg=mods_avg,mods_tab=mods_tab,z1.prop=z1.prop,comp.mat=comp.mat,nonc.mat=nonc.mat,yj=yj,tab.hpd=tab.hpd, mucomps.m=mucomps.m, munoncs.m=munoncs.m, log.tab=log.tab)	
	return(ret)

}


mix2normals = function(type,y,grp,nic=NULL,inits,hp,M,nburn,notmix=T,seed,pair=F){
###Gibbs samplers for mixture of 2 normal distributions assuming equal variance
#type: Gibbs sampler to use:
##'ind' estimates all parameters
##'rel' assumes relationship between UB/15.8 group and compliant group means for other nicotine levels
#y: data in a vector
#grp: vector identifying group membership for each y_i numerically (i.e., 1=VLNC,2=1.3 mg/g, etc.)
#nic: the nicotine group having its mixture distribution estimated as a nicotine value: e.g., 0.4, 1.3, 2.4, 5.2, 15.8
#inits: list of initial values to use for mu, theta, tau, p [note: tau is precision]
##if one value given, it is assumed to be the same across all groups; otherwise give vector with value for each group!!!
#hp: list of hyperparameter values (for mu, theta, tau, and p) 
#M: number of iterations for Gibbs sampler
#nburn: number of iterations to discard from beginning of chain
#notmix: indicator if last group shouldn't estimate entire mixture since all are assumed compliant, default is FALSE
#seed: set seed to replicate results
#pair: if TRUE, return pairwise matrices for concordance of compliance and non-compliance	

	set.seed(seed)

	yj <- lapply( split(y, grp), sort) #split data and sort in one step
	j <- length(yj) #number of groups
	nj <- sapply( yj, length) #calculate number of observations within each group

	draws=matrix(NA,M,7*j)
	colnames(draws) <- c(paste0('p',1:j),paste0('mu_comp',1:j),paste0('mu_nonc',1:j),paste0('sigma',1:j),paste0('c80_',1:j),paste0('c90_',1:j),paste0('c95_',1:j))

	#initial values (if only one given, assume it is initial value for each level/group)
	if( length(inits$mu_comp)==1 ){ mu_compj <- rep(inits$mu_comp, j) }else{ mu_compj <- inits$mu_comp }
	if( length(inits$mu_nonc)==1 ){ mu_noncj <- rep(inits$mu_nonc, j) }else{ mu_noncj <- inits$mu_nonc }
	if( length(inits$tau)==1 ){ tauj <- rep(inits$tau, j) }else{ tauj <- inits$tau }
	if( length(inits$p)==1 ){ pj <- rep(inits$p, j) }else{ pj <- inits$p }
	if( (length(mu_compj)!=j | length(mu_noncj)!=j | length(tauj)!=j | length(pj)!=j) == TRUE ){ stop('Please either give inits of length 1 or length j. At least one is not equal to 1 or j!') } #break from function if inits are not of correct length

	#Convert hyperparameter values to be for each group if separate values not provided for each group (if only one given, assume it is used for each level/group)
	if( length(hp$tau.a)==1 ){ hp$tau.a <- rep(hp$tau.a,j) }	
	if( length(hp$tau.b)==1 ){ hp$tau.b <- rep(hp$tau.b,j) }	
	if( length(hp$p.alpha)==1 ){ hp$p.alpha <- rep(hp$p.alpha,j) }	
	if( length(hp$p.beta)==1 ){ hp$p.beta <- rep(hp$p.beta,j) }	
	if( length(hp$mu.prec)==1 ){ hp$mu.prec <- rep(hp$mu.prec,j) }	
	if( length(hp$theta.prec)==1 ){ hp$theta.prec <- rep(hp$theta.prec,j) }	
	if( (length(hp$tau.a)!=j | length(hp$tau.b)!=j | length(hp$p.alpha)!=j | length(hp$p.beta)!=j | length(hp$mu.prec)!=j | length(hp$theta.prec)!=j) == TRUE ){ stop('Please give hp of length 1 or length j. At least one is not equal to 1 or j!')}

	#intialize latent indicator variable 
	zj <- lapply(1:j, function(x) c(0, rep(NA, nj[x]-2), 1) ) #initialize latent indicator for each group

	#initialize latent pairwise classification matrices
	zj.sum <- lapply(1:j, function(x) rep(0, nj[x]) ) #list keeping track of number of iterations after burn-in where y_ij is classified as k=1
	comp.mat <- matrix(0, nrow=sum(nj), ncol=sum(nj) ) 
	nonc.mat <- matrix(0, nrow=sum(nj), ncol=sum(nj) ) 
	#comp.mat <- lapply(1:j, function(x) matrix(0, nrow=nj[x], ncol=nj[x]) )
	#nonc.mat <- lapply(1:j, function(x) matrix(0, nrow=nj[x], ncol=nj[x]) )

	for (iter in 1:M){
		pvecj <- lapply(1:length(yj), function(x) pj[x]*dnorm(yj[[x]], mean=mu_noncj[x], sd=sqrt(1/tauj[x])) / ( (1-pj[x])*dnorm(yj[[x]], mean=mu_compj[x], sd=sqrt(1/tauj[x])) + pj[x]*dnorm(yj[[x]], mean=mu_noncj[x], sd=sqrt(1/tauj[x])) ) )
		if(notmix==T){ pvecj[[j]] <- rep(0, nj[j] ) } #if notmix=TRUE, set all values as 0 to indicate membership in compliant group

		if(type=='ind'){gibbs.step <- gibbs_normalmixture_ind(yj=yj,nj=nj,pvecj=pvecj,muj=mu_compj,thetaj=mu_noncj-mu_compj,tauj=tauj,hpj=hp,notmix=notmix,j=j)}
		if(type=='rel'){gibbs.step <- gibbs_normalmixture_rel(yj=yj,nj=nj,pvecj=pvecj,muj=mu_compj,thetaj=mu_noncj-mu_compj,tauj=tauj,hpj=hp,j=j,nic=nic)}

		mu_compj <- gibbs.step$mu_comp
		mu_noncj <- gibbs.step$mu_nonc
		tauj <- gibbs.step$tau
		sigmaj <- gibbs.step$sigma
		zj <- gibbs.step$z
		pj <- gibbs.step$p

		#Save pairwise comparisons
		if(iter > nburn){
			zj.sum <- lapply(1:j, function(x) zj[[x]] + zj.sum[[x]]) #add this iterations zj values to sum

			#Updating and storing pairwise matrices can be time intensive, only calculate if pair==T
			if(pair==T){
				zj.unlist <- unlist(zj)
				comp.mat <- pairwiseZmat(z=zj.unlist,zval=0) + comp.mat
				nonc.mat <- pairwiseZmat(z=zj.unlist,zval=1) + nonc.mat
			}
		}

		comp.perc <- sapply(1:j, function(x) qnorm( c(.8,.9,.95), mean=mu_compj[x], sd=sigmaj[x]))
		comp.percj <- as.vector( t(comp.perc) )

		draws[iter,] <- c(pj,mu_compj,mu_noncj,sigmaj,comp.percj)
	}

	if(nburn == 0){ ret <- list( samp=draws, zj.sum=zj.sum, comp.mat=comp.mat, nonc.mat=nonc.mat ) }
	if(nburn > 0){ ret <- list( samp=draws[-(1:nburn),], zj.sum=zj.sum, comp.mat=comp.mat, nonc.mat=nonc.mat ) }

	return(ret)

}



rjmcmc_2normals = function(y,inits,hp,hp.rj=NULL,prob,M,seed,nburn,grp,nic,pair=F){
###Gibbs samplers for mixture of 2 normal distributions assuming equal variance
#y: data in a vector in order from smallest to largest
#inits: list of initial values to use for the model, compliant mean, noncompliant mean, tau, p [note: tau is precision]
#hp: list of hyperparameter values (for lambda, theta, tau, and p) 
#hp.rj: list of hyperparameter values for parameters added/removed in RJ step (currently NULL since proposal distribution is conditional posterior of mu_j)
#prob: prior probability for level j using all components model (model 0)
#M: number of iterations for Gibbs sampler
#seed: seed value to set (in order to replicate results)
#nburn: number of iterations to discard from beginning of chain
#grp: a vector identifying which group each data belongs to (currently 'c' vs. 'ub')
#nic: nicotine level for group j
#pair: if TRUE, return pairwise matrices for concordance of compliance and non-compliance	

	set.seed(seed)

	yj <- lapply( split(y, grp), sort) #split data and sort in one step
	j <- length(yj) #number of groups
	nj <- sapply( yj, length) #calculate number of observations within each group

	draws=matrix(NA,M,(7*j+2))
	colnames(draws) <- c('model',paste0('p',1:j),paste0('mu_comp',1:j),paste0('mu_nonc',1:j),paste0('sigma',1:j),paste0('c80_',1:j),paste0('c90_',1:j),paste0('c95_',1:j),'dev')

	#initial values (if only one given, assume it is initial value for each level/group)
	if( length(inits$mu_comp)==1 ){ mu_compj <- rep(inits$mu_comp, j) }else{ mu_compj <- inits$mu_comp }
	if( length(inits$mu_nonc)==1 ){ mu_noncj <- rep(inits$mu_nonc, j) }else{ mu_noncj <- inits$mu_nonc }
	if( length(inits$tau)==1 ){ tauj <- rep(inits$tau, j) }else{ tauj <- inits$tau }
	if( length(inits$p)==1 ){ pj <- rep(inits$p, j) }else{ pj <- inits$p }
	if( length(inits$model)==1 ){model <- c(rep(inits$model, (j-1)),1)}else{ model <- c(inits$model[1:(j-1)],1) } #always assume Jth model is 1 for relationship
	if( (length(mu_compj)!=j | length(mu_noncj)!=j | length(tauj)!=j | length(pj)!=j | length(model)!=j ) == TRUE ){ stop('Please either give inits of length 1 or length j. At least one is not equal to 1 or j!') } #break from function if inits are not of correct length

	#Convert hyperparameter values to be for each group if separate values not provided for each group (if only one given, assume it is used for each level/group)
	if( length(hp$tau.a)==1 ){ hp$tau.a <- rep(hp$tau.a,j) }	
	if( length(hp$tau.b)==1 ){ hp$tau.b <- rep(hp$tau.b,j) }	
	if( length(hp$p.alpha)==1 ){ hp$p.alpha <- rep(hp$p.alpha,j) }	
	if( length(hp$p.beta)==1 ){ hp$p.beta <- rep(hp$p.beta,j) }	
	if( length(hp$mu.prec)==1 ){ hp$mu.prec <- rep(hp$mu.prec,j) }	
	if( length(hp$theta.prec)==1 ){ hp$theta.prec <- rep(hp$theta.prec,j) }	
	if( (length(hp$tau.a)!=j | length(hp$tau.b)!=j | length(hp$p.alpha)!=j | length(hp$p.beta)!=j | length(hp$mu.prec)!=j | length(hp$theta.prec)!=j) == TRUE ){ stop('Please give hp of length 1 or length j. At least one is not equal to 1 or j!')}

	#intialize latent indicator variable 
	zj <- lapply(1:length(yj), function(x) c(0, rep(NA, nj[x]-2), 1) ) #initialize latent indicator for each group

	#initialize latent pairwise classification matrices
	zj.sum <- lapply(1:j, function(x) rep(0, nj[x]) ) #list keeping track of number of iterations after burn-in where y_ij is classified as k=1
	comp.mat <- lapply(1:j, function(x) matrix(0, nrow=nj[x], ncol=nj[x]) )
	nonc.mat <- lapply(1:j, function(x) matrix(0, nrow=nj[x], ncol=nj[x]) )

	for (iter in 1:M){
		#Update MCMC according to model state currently in
		pvecj <- lapply(1:length(yj), function(x) pj[x]*dnorm(yj[[x]], mean=mu_noncj[x], sd=sqrt(1/tauj[x])) / ( (1-pj[x])*dnorm(yj[[x]], mean=mu_compj[x], sd=sqrt(1/tauj[x])) + pj[x]*dnorm(yj[[x]], mean=mu_noncj[x], sd=sqrt(1/tauj[x])) ) )
		pvecj[[j]] <- rep(0, nj[j] ) #replace group J with all compliant indicators

		#update step with Gibbs sampler for given model assumed for each level
		gibbs.step <- gibbs_normalmixture_rjmcmc(yj=yj,nj=nj,pvecj=pvecj,muj=mu_compj,thetaj=mu_noncj-mu_compj,tauj=tauj,hpj=hp,j=j,model=model,nic=nic)

		mu_compj <- gibbs.step$mu_comp
		mu_noncj <- gibbs.step$mu_nonc
		tauj <- gibbs.step$tau
		sigmaj <- gibbs.step$sigma
		zj <- gibbs.step$z
		pj <- gibbs.step$p

		#Save pairwise comparisons
		if(iter > nburn){
			zj.sum <- lapply(1:j, function(x) zj[[x]] + zj.sum[[x]]) #add this iterations zj values to sum

			#Updating and storing pairwise matrices can be time intensive, only calculate if pair==T
			if(pair==T){
				comp.mat <- lapply(1:j, function(x) pairwiseZmat(z=zj[[x]],zval=0) + comp.mat[[x]] ) #determine which pairwise zj values are both compliant
				nonc.mat <- lapply(1:j, function(x) pairwiseZmat(z=zj[[x]],zval=1) + nonc.mat[[x]] ) #determine which pairwise zj values are both non-compliant
			}
		}

		#Complete RJ step to determine if we switch models and update 
		rjstep <- updatemodel(yj=yj,pj=pj,j=j,nj=nj,tauj=tauj,zj=zj,muj=mu_compj,thetaj=mu_noncj-mu_compj,model=model,hpj=hp,prob=prob,nic=nic)
		model <- rjstep$model
		mu_compj <- rjstep$muj
		dev <- rjstep$dev

		mod10 <- baseKto10(x=model[1:(j-1)],K=2)

		comp.perc <- sapply(1:j, function(x) qnorm( c(.8,.9,.95), mean=mu_compj[x], sd=sigmaj[x]))
		comp.percj <- as.vector( t(comp.perc) )

		draws[iter,] <- c(mod10,pj,mu_compj,mu_noncj,sigmaj,comp.percj, dev)
	}

	if(nburn == 0){ ret <- list( samp=draws, zj.sum=zj.sum, comp.mat=comp.mat, nonc.mat=nonc.mat) }
	if(nburn > 0){ ret <- list( samp=draws[-(1:nburn),], zj.sum=zj.sum, comp.mat=comp.mat, nonc.mat=nonc.mat) }

	return(ret)

}


setmcmc_2normals = function(y,inits,hp,prob,M,seed,nburn,grp,nic,pair=F){
###Gibbs samplers for mixture of 2 normal distributions assuming equal variance where the modeling assumption for each dose-level is FIXED (as opposed to using the RJMCMC to consider exploring both relationships)
#y: data in a vector in order from smallest to largest
#inits: list of initial values to use for the model, compliant mean, noncompliant mean, tau, p [note: tau is precision]
##NOTE: model in inits is the fixed model indicators to use for the entire chain (no jumping)
#hp: list of hyperparameter values (for lambda, theta, tau, and p) 
#prob: prior probability for level j using all components model (model 0)
#M: number of iterations for Gibbs sampler
#seed: seed value to set (in order to replicate results)
#nburn: number of iterations to discard from beginning of chain
#grp: a vector identifying which group each data belongs to (currently 'c' vs. 'ub')
#nic: nicotine level for group j
#pair: if TRUE, return pairwise matrices for concordance of compliance and non-compliance	

	set.seed(seed)

	yj <- lapply( split(y, grp), sort) #split data and sort in one step
	j <- length(yj) #number of groups
	nj <- sapply( yj, length) #calculate number of observations within each group

	draws=matrix(NA,M,(7*j+2))
	colnames(draws) <- c('model',paste0('p',1:j),paste0('mu_comp',1:j),paste0('mu_nonc',1:j),paste0('sigma',1:j),paste0('c80_',1:j),paste0('c90_',1:j),paste0('c95_',1:j),'dev')

	#initial values (if only one given, assume it is initial value for each level/group)
	if( length(inits$mu_comp)==1 ){ mu_compj <- rep(inits$mu_comp, j) }else{ mu_compj <- inits$mu_comp }
	if( length(inits$mu_nonc)==1 ){ mu_noncj <- rep(inits$mu_nonc, j) }else{ mu_noncj <- inits$mu_nonc }
	if( length(inits$tau)==1 ){ tauj <- rep(inits$tau, j) }else{ tauj <- inits$tau }
	if( length(inits$p)==1 ){ pj <- rep(inits$p, j) }else{ pj <- inits$p }
	if( length(inits$model)==1 ){model <- c(rep(inits$model, (j-1)),1)}else{ model <- c(inits$model[1:(j-1)],1) } #always assume Jth model is 1 for relationship
	if( (length(mu_compj)!=j | length(mu_noncj)!=j | length(tauj)!=j | length(pj)!=j | length(model)!=j ) == TRUE ){ stop('Please either give inits of length 1 or length j. At least one is not equal to 1 or j!') } #break from function if inits are not of correct length

	#Convert hyperparameter values to be for each group if separate values not provided for each group (if only one given, assume it is used for each level/group)
	if( length(hp$tau.a)==1 ){ hp$tau.a <- rep(hp$tau.a,j) }	
	if( length(hp$tau.b)==1 ){ hp$tau.b <- rep(hp$tau.b,j) }	
	if( length(hp$p.alpha)==1 ){ hp$p.alpha <- rep(hp$p.alpha,j) }	
	if( length(hp$p.beta)==1 ){ hp$p.beta <- rep(hp$p.beta,j) }	
	if( length(hp$mu.prec)==1 ){ hp$mu.prec <- rep(hp$mu.prec,j) }	
	if( length(hp$theta.prec)==1 ){ hp$theta.prec <- rep(hp$theta.prec,j) }	
	if( (length(hp$tau.a)!=j | length(hp$tau.b)!=j | length(hp$p.alpha)!=j | length(hp$p.beta)!=j | length(hp$mu.prec)!=j | length(hp$theta.prec)!=j) == TRUE ){ stop('Please give hp of length 1 or length j. At least one is not equal to 1 or j!')}

	#intialize latent indicator variable 
	zj <- lapply(1:length(yj), function(x) c(0, rep(NA, nj[x]-2), 1) ) #initialize latent indicator for each group

	#initialize latent pairwise classification matrices
	zj.sum <- lapply(1:j, function(x) rep(0, nj[x]) ) #list keeping track of number of iterations after burn-in where y_ij is classified as k=1
	comp.mat <- lapply(1:j, function(x) matrix(0, nrow=nj[x], ncol=nj[x]) )
	nonc.mat <- lapply(1:j, function(x) matrix(0, nrow=nj[x], ncol=nj[x]) )

	for (iter in 1:M){
		#Update MCMC according to model state currently in
		pvecj <- lapply(1:length(yj), function(x) pj[x]*dnorm(yj[[x]], mean=mu_noncj[x], sd=sqrt(1/tauj[x])) / ( (1-pj[x])*dnorm(yj[[x]], mean=mu_compj[x], sd=sqrt(1/tauj[x])) + pj[x]*dnorm(yj[[x]], mean=mu_noncj[x], sd=sqrt(1/tauj[x])) ) )
		pvecj[[j]] <- rep(0, nj[j] ) #replace group J with all compliant indicators

		#update step with Gibbs sampler for given model assumed for each level
		gibbs.step <- gibbs_normalmixture_rjmcmc(yj=yj,nj=nj,pvecj=pvecj,muj=mu_compj,thetaj=mu_noncj-mu_compj,tauj=tauj,hpj=hp,j=j,model=model,nic=nic)

		mu_compj <- gibbs.step$mu_comp
		mu_noncj <- gibbs.step$mu_nonc
		tauj <- gibbs.step$tau
		sigmaj <- gibbs.step$sigma
		zj <- gibbs.step$z
		pj <- gibbs.step$p

		#Save pairwise comparisons
		if(iter > nburn){
			zj.sum <- lapply(1:j, function(x) zj[[x]] + zj.sum[[x]]) #add this iterations zj values to sum

			#Updating and storing pairwise matrices can be time intensive, only calculate if pair==T
			if(pair==T){
				comp.mat <- lapply(1:j, function(x) pairwiseZmat(z=zj[[x]],zval=0) + comp.mat[[x]] ) #determine which pairwise zj values are both compliant
				nonc.mat <- lapply(1:j, function(x) pairwiseZmat(z=zj[[x]],zval=1) + nonc.mat[[x]] ) #determine which pairwise zj values are both non-compliant
			}
		}

		dev <- NA #used for RJ step
		mod10 <- baseKto10(x=model[1:(j-1)],K=2)

		comp.perc <- sapply(1:j, function(x) qnorm( c(.8,.9,.95), mean=mu_compj[x], sd=sigmaj[x]))
		comp.percj <- as.vector( t(comp.perc) )

		draws[iter,] <- c(mod10,pj,mu_compj,mu_noncj,sigmaj,comp.percj, dev)
	}

	if(nburn == 0){ ret <- list( samp=draws, zj.sum=zj.sum, comp.mat=comp.mat, nonc.mat=nonc.mat) }
	if(nburn > 0){ ret <- list( samp=draws[-(1:nburn),], zj.sum=zj.sum, comp.mat=comp.mat, nonc.mat=nonc.mat) }

	return(ret)

}


gibbs_normalmixture_ind <- function(yj,nj,j,pvecj,muj,thetaj,tauj,hpj,notmix=T){
###This function provides the code for updating the Gibbs sampler for a mixture of two normal distributions where all components are updated across j-levels
#yj: data in a list where each item (vector) includes all observations (y_i) for group j (should be ordered smallest to largest within each group)
#nj: sample size of each group in vector by group
#j: number of groups
#pvecj: corresponding list to y for probability of each person belonging to non-compliant group within group j
#muj: current estimate for mu parameters for compliant group components
#thetaj: current estimate for theta parameter for non-compliant group components
#tauj: current estimate for shared estimate of precision in mixture distribution within each group j
#hpj: list of hyperparameter values (for mu, theta, tau, and p) 
#notmix: indicator if last group shouldn't estimate entire mixture since all are assumed compliant, default is TRUE

	zj <- lapply(1:j, function(x) c(0,rbinom(n=(nj[x]-2), size=1, pvecj[[x]][2:(nj[x]-1)]),1) ) #generate new latent indicators from Bernoulli with smallest observation fixed as compliant and largest fixed as non-compliant
	if(notmix==T){ zj[[j]] <- rep(0, length(zj[[j]])) } #if highest group/level is only compliant, convert all indicators to 0
	
	ncj <- sapply(1:j, function(x) sum(zj[[x]]==0) )
	nnj <- sapply(1:j, function(x) sum(zj[[x]]==1) )

	ycj <- lapply(1:j, function(x) yj[[x]][which(zj[[x]]==0)] )
	ynj <- lapply(1:j, function(x) yj[[x]][which(zj[[x]]==1)] )

	###calculate component means and variances for each group
	#compliant components
	ycm <- sapply(1:j, function(x) mean(ycj[[x]]) )
	ycv <- sapply(1:j, function(x) if( length(ycj[[x]])==1 ){ 0 }else{ var(ycj[[x]]) }) #place variance of 0 if only one observation to avoid errors
	#non-compliant components
	ynm <- sapply(1:j, function(x) if( length(ynj[[x]])==0 ){0}else{ mean(ynj[[x]]) })
	ynv <- sapply(1:j, function(x) if( length(ynj[[x]])%in%c(0,1) ){ 0 }else{ var(ynj[[x]]) })  #place variance of 0 if only one observation to avoid errors

	# sampling p: probability of membership in group 2 (non-compliers)
	pj <- rbeta(j, nnj + hpj$p.alpha, ncj + hpj$p.beta)

	# sampling tau
	shape.val <- (nj/2) + hpj$tau.a
	rate.val <- 0.5 * ((ncj-1)*ycv + (nnj-1)*ynv + ncj*(ycm-muj)^2 + nnj*(ynm-(muj+thetaj))^2 + 2*hpj$tau.b )
	tauj <- rgamma(j,shape=shape.val,rate=rate.val)
	sigmaj <- 1/sqrt(tauj)

	# sampling mu
	m.mean <- -( tauj * (nnj*(thetaj - ynm) - ncj*ycm) ) / ( nj*tauj + hpj$mu.prec )
	m.prec <- ( nj*tauj + hpj$mu.prec )
	muj <- rnorm(j, mean=m.mean, sd=sqrt(1/m.prec))

	# sampling theta
	t.mean <- -( nnj*tauj*(muj - ynm) ) / (nnj*tauj + hpj$theta.prec)
	t.prec <- (nnj*tauj + hpj$theta.prec)
	thetaj <- abs( rnorm(j, mean=t.mean, sd=sqrt(1/t.prec)))
	if(notmix==T){ thetaj[j] <- 0 } #replace largest group theta with 0 if notmix is TRUE

	ret <- list(mu_comp=muj, mu_nonc=muj+thetaj, mu=muj, theta=thetaj, tau=tauj, sigma=sigmaj, z=zj, p=pj)

}


gibbs_normalmixture_rel <- function(yj,nj,j,pvecj,muj,thetaj,tauj,hpj,nic){
###This function provides the code for updating the Gibbs sampler for a mixture of two normal distributions where the compliant group means are specified in relationship to the UB/15.8 group
#yj: data in a list where each item (vector) includes all observations (y_i) for group j (should be ordered smallest to largest within each group)
#nj: sample size of each group in vector by group
#j: number of groups (the final Jth group should be UB/15.8 or whatever reference group is defined)
#pvecj: corresponding list to y for probability of each person belonging to non-compliant group within group j
#muj: current estimate for mu_J parameter for compliant group
#thetaj: current estimate for theta parameter for non-compliant group components
#tauj: current estimate for shared estimate of precision in mixture distribution within each group j
#hpj: list of hyperparameter values (for mu, theta, tau, and p) 
#nic: nicotine level for group j

	dj <- log(nic/nic[length(nic)])

	zj <- lapply(1:j, function(x) c(0,rbinom(n=(nj[x]-2), size=1, pvecj[[x]][2:(nj[x]-1)]),1) ) #generate new latent indicators from Bernoulli with smallest observation fixed as compliant and largest fixed as non-compliant
	zj[[j]] <- rep(0, length(zj[[j]])) #highest group/level is only compliant, convert all indicators to 0
	
	ncj <- sapply(1:j, function(x) sum(zj[[x]]==0) )
	nnj <- sapply(1:j, function(x) sum(zj[[x]]==1) )

	ycj <- lapply(1:j, function(x) yj[[x]][which(zj[[x]]==0)] )
	ynj <- lapply(1:j, function(x) yj[[x]][which(zj[[x]]==1)] )

	###calculate component means and variances for each group
	#compliant components
	ycm <- sapply(1:j, function(x) mean(ycj[[x]]) )
	ycv <- sapply(1:j, function(x) if( length(ycj[[x]])==1 ){ 0 }else{ var(ycj[[x]]) }) #place variance of 0 if only one observation to avoid errors
	#non-compliant components
	ynm <- sapply(1:j, function(x) if( length(ynj[[x]])==0 ){0}else{ mean(ynj[[x]]) })
	ynv <- sapply(1:j, function(x) if( length(ynj[[x]])%in%c(0,1) ){ 0 }else{ var(ynj[[x]]) })  #place variance of 0 if only one observation to avoid errors

	# sampling p: probability of membership in group 2 (non-compliers)
	pj <- rbeta(j, nnj + hpj$p.alpha, ncj + hpj$p.beta)

	# sampling tau
	shape.val <- (nj/2) + hpj$tau.a
	rate.val <- 0.5 * ((ncj-1)*ycv + (nnj-1)*ynv + ncj*(ycm - (dj + muj[j]))^2 + nnj*(ynm-(dj+muj[j]+thetaj))^2 + 2*hpj$tau.b )
	tauj <- rgamma(j,shape=shape.val,rate=rate.val)
	sigmaj <- 1/sqrt(tauj)

	# sampling mu_J
	g <- tauj*(nj*dj + nnj*thetaj - ncj*ycm - nnj*ynm)
	H1 <- sum(tauj*nj) + hpj$mu.prec[j]
	H2 <- sum(g[1:(j-1)]) - tauj[j]*nj[j]*ycm[j]
	m.mean <- -H2/H1
	m.prec <- H1
	muJ <- rnorm(1, mean=m.mean, sd=sqrt(1/m.prec))
	muj <- muJ+dj

	# sampling theta
	t.mean <- -( nnj*tauj*(dj + muj[j] - ynm) ) / (nnj*tauj + hpj$theta.prec)
	t.prec <- (nnj*tauj + hpj$theta.prec)
	thetaj <- abs( rnorm(j, mean=t.mean, sd=sqrt(1/t.prec)))
	thetaj[j] <- 0 #replace largest group theta with 0 since everyone is compliant

	ret <- list(mu_comp=muj, mu_nonc=muj+thetaj, mu=muj, theta=thetaj, tau=tauj, sigma=sigmaj, z=zj, p=pj)

}


gibbs_normalmixture_rjmcmc <- function(yj,nj,j,pvecj,muj,thetaj,tauj,hpj,model,nic){
###This function provides the code for updating the Gibbs sampler for a mixture of two normal distributions where all components are updated across j-levels with chosen approach
#yj: data in a list where each item (vector) includes all observations (y_i) for group j (should be ordered smallest to largest within each group)
#nj: sample size of each group in vector by group
#j: number of groups
#pvecj: corresponding list to y for probability of each person belonging to non-compliant group within group j
#muj: current estimate for mu parameters for compliant group components
#thetaj: current estimate for theta parameter for non-compliant group components
#tauj: current estimate for shared estimate of precision in mixture distribution within each group j
#hpj: list of hyperparameter values (for mu, theta, tau, and p) 
#model: a vector of indicators of the current model we are in for level j (0=all components, 1=relationship between UB/15.8 and level j)
#nic: nicotine level for group j

	j.all <- which( model==0 )
	j.rel <- which( model==1 )
	dj <- log(nic/nic[length(nic)])

	zj <- lapply(1:j, function(x) c(0,rbinom(n=(nj[x]-2), size=1, pvecj[[x]][2:(nj[x]-1)]),1) ) #generate new latent indicators from Bernoulli with smallest observation fixed as compliant and largest fixed as non-compliant
	zj[[j]] <- rep(0, length(zj[[j]])) #if highest group/level is only compliant, convert all indicators to 0
	
	ncj <- sapply(1:j, function(x) sum(zj[[x]]==0) )
	nnj <- sapply(1:j, function(x) sum(zj[[x]]==1) )

	ycj <- lapply(1:j, function(x) yj[[x]][which(zj[[x]]==0)] )
	ynj <- lapply(1:j, function(x) yj[[x]][which(zj[[x]]==1)] )

	###calculate component means and variances for each group
	#compliant components
	ycm <- sapply(1:j, function(x) mean(ycj[[x]]) )
	ycv <- sapply(1:j, function(x) if( length(ycj[[x]])==1 ){ 0 }else{ var(ycj[[x]]) }) #place variance of 0 if only one observation to avoid errors
	#non-compliant components
	ynm <- sapply(1:j, function(x) if( length(ynj[[x]])==0 ){0}else{ mean(ynj[[x]]) })
	ynv <- sapply(1:j, function(x) if( length(ynj[[x]])%in%c(0,1) ){ 0 }else{ var(ynj[[x]]) })  #place variance of 0 if only one observation to avoid errors

	# sampling p: probability of membership in group 2 (non-compliers)
	pj <- rbeta(j, nnj + hpj$p.alpha, ncj + hpj$p.beta)

	# sampling tau
	shape.val <- (nj/2) + hpj$tau.a
	rate.val <- 0.5 * ((ncj-1)*ycv + (nnj-1)*ynv + ncj*(ycm-muj)^2 + nnj*(ynm-(muj+thetaj))^2 + 2*hpj$tau.b )
	tauj <- rgamma(j,shape=shape.val,rate=rate.val)
	sigmaj <- 1/sqrt(tauj)

	### sampling mu
	#levels estimating compliant mean with own information only
	m.mean_ind <- -( tauj * (nnj*(thetaj - ynm) - ncj*ycm) ) / ( nj*tauj + hpj$mu.prec )
	m.prec_ind <- ( nj*tauj + hpj$mu.prec )
	muj_ind <- rnorm( length(j.all) , mean=m.mean_ind[j.all], sd=sqrt(1/m.prec_ind[j.all]))

	#groups calculating compliant mean based on relationship to mu_J
	g <- tauj[j.rel]*(nj[j.rel]*dj[j.rel] + nnj[j.rel]*thetaj[j.rel] - ncj[j.rel]*ycm[j.rel] - nnj[j.rel]*ynm[j.rel])
	H1 <- sum(tauj[j.rel]*nj[j.rel]) + hpj$mu.prec[j]
	H2 <- sum(g[0:(length(j.rel)-1)]) - tauj[j]*nj[j]*ycm[j]
	m.mean_rel <- -H2/H1
	m.prec_rel <- H1
	muJ <- rnorm(1, mean=m.mean_rel, sd=sqrt(1/m.prec_rel))
	muj_rel <- muJ+dj[j.rel]

	#combine muj estimates from both approaches
	muj[j.all] <- muj_ind
	muj[j.rel] <- muj_rel

	# sampling theta
	t.mean <- -( nnj*tauj*(muj - ynm) ) / (nnj*tauj + hpj$theta.prec)
	t.prec <- (nnj*tauj + hpj$theta.prec)
	thetaj <- abs( rnorm(j, mean=t.mean, sd=sqrt(1/t.prec)))
	thetaj[j] <- 0 #replace largest group theta with 0 if notmix is TRUE

	ret <- list(mu_comp=muj, mu_nonc=muj+thetaj, mu=muj, theta=thetaj, tau=tauj, sigma=sigmaj, z=zj, p=pj)

}


updatemodel <- function(yj,nj,j,pj,tauj,zj,muj,thetaj,model,hpj,hp.rj=NULL,prob,nic){
###Function to complete the "reversible jump" portion of the RJMCMC to move between our two proposed models
#yj: data in a list where each item (vector) includes all observations (y_i) for group j (should be ordered smallest to largest within each group)
#nj: sample size of each group in vector by group
#j: number of groups
#pj: current estimate for proportion in non-compliant group for each level j
#tauj: current estimate for shared estimate of precision in mixture distribution within each group j
#zj: latent variable classifying y_i subject in compliant or non-compliant group
#muj: current estimate for mu parameters for compliant group components
#thetaj: current estimate for theta parameter for non-compliant group components
#model: a vector of indicators of the current model we are in for level j (0=all components, 1=relationship between UB/15.8 and level j)
#hpj: list of hyperparameter values (for mu, theta, tau, and p) 
#hp.rj: list of hyperparameter values for RJ step (for each muj) (currently NULL since proposal distribution is conditional posterior with same priors assumed)
#prob: prior probability for level j using all components model (model 0)
#nic: the nicotine group having its mixture distribution estimated as a nicotine value: e.g., 0.4, 1.3, 2.4, 5.2, 15.8

	ncj <- sapply(1:j, function(x) sum(zj[[x]]==0) )
	nnj <- sapply(1:j, function(x) sum(zj[[x]]==1) )

	ycj <- sapply(1:j, function(x) yj[[x]][which(zj[[x]]==0)] )
	ynj <- sapply(1:j, function(x) yj[[x]][which(zj[[x]]==1)] )

	###calculate component means and variances for each group
	#compliant components
	ycm <- sapply(1:j, function(x) mean(ycj[[x]]) )
	ycv <- sapply(1:j, function(x) if( length(ycj[[x]])==1 ){ 0 }else{ var(ycj[[x]]) }) #place variance of 0 if only one observation to avoid errors
	#non-compliant components
	ynm <- sapply(1:j, function(x) if( length(ynj[[x]])==0 ){0}else{ mean(ynj[[x]]) })
	ynv <- sapply(1:j, function(x) if( length(ynj[[x]])%in%c(0,1) ){ 0 }else{ var(ynj[[x]]) })  #place variance of 0 if only one observation to avoid errors

	###Identify current model state and propose new variables for other model
	#Determine which level to propose different model for
	r <- sample(1:(j-1),1)
	r.state <- model[r]

	muj_prop <- muj #create proposal muj to update with value below depending on current r.state
	model_prop <- model; if(r.state==0){model_prop[r] <- 1}else{model_prop[r] <- 0}
	m.mean_r <- (-( tauj * (nnj*(thetaj - ynm) - ncj*ycm) ) / ( nj*tauj + hpj$mu.prec ))[r]
	m.prec_r <- ( nj*tauj + hpj$mu.prec )[r]

	if(r.state==0){ #model 0 is all components, so propose values under relationship assumption
		muj_prop[r] <- muj[j] + log(nic[r]/nic[length(nic)])		
	}else{ #model 1 is relationship, so propose values under all components model
		muj_prop[r] <- rnorm(1, mean=m.mean_r, sd=sqrt(1/m.prec_r))
	}

	modprob.num <- prod( ( prob^abs(model_prop-1) * (1-prob)^model_prop)[1:(j-1)] )
	modprob.den <- prod( ( prob^abs(model-1) * (1-prob)^model)[1:(j-1)] )

	###Calculate acceptance probability
	#Calculate log-likelihoods
	loglikj <- sapply(1:j, function(x) sum( log( (1-pj[x])*dnorm(yj[[x]], mean=muj[x], sd=sqrt(1/tauj[x])) + pj[x]*dnorm(yj[[x]], mean=muj[x]+thetaj[x], sd=sqrt(1/tauj[x])) ) ) )
	loglik <- sum(loglikj)

	newloglikj <- sapply(1:j, function(x) sum( log( (1-pj[x])*dnorm(yj[[x]], mean=muj_prop[x], sd=sqrt(1/tauj[x])) + pj[x]*dnorm(yj[[x]], mean=muj_prop[x]+thetaj[x], sd=sqrt(1/tauj[x])) ) ) )
	newloglik <- sum(newloglikj)

	#if currently in model 0 (est all), look at num being proposed move for level r to model 1 (relationship assumed)
	if(r.state == 0){ 
		num <- newloglik + log(dnorm(muj_prop[r], mean=m.mean_r, sd=sqrt(1/m.prec_r))) + log(modprob.num)
		den <- loglik + log(dnorm(muj[r], mean=0, sd=sqrt(1/hpj$mu.prec[r]))) + log(modprob.den)
	}else{ #if currently in model 1 (relationship), look at num being proposed move for level r to model 0 (estimate level with level's data)
		num <- newloglik + log(dnorm(muj_prop[r], mean=0, sd=sqrt(1/hpj$mu.prec[r]))) + log(modprob.num)
		den <- loglik + log(dnorm(muj[r], mean=m.mean_r, sd=sqrt(1/m.prec_r))) + log(modprob.den)
	}

	###Acceptance probability
	A <- min(1, exp(num-den))
	u <- runif(1)
	if(u <= A){
		if(r.state == 0){model[r] <- 1}else{model[r] <- 0}
		muj <- muj_prop
		dev <- -2*newloglik #deviance for diagnostics
	}else{
		dev <- -2*loglik
	}

	ret <- list(model=model, muj=muj, mu1j=muj, dev=dev)

}



#####################################################################
## Simulation Functions to Implement MCMC for Mixture Distribution ##
#####################################################################

mix_sim <- function(x, meanlog2, sigma, pcomp, nsamp, es, nic, nburn, nchain, simdate ){
###Function to fit Gaussian mixture distributions for simulated data
#x: seed to set
#meanlog2: non-compliant component meanlog
#sigma: value to use for sdlog (lognormal)
#pcomp: proportion compliant (used to determine how many to simulate from each group)
#nsamp: sample size to simulate
#es: effect size for each level to simulate compliant group
#nic: effect size to match to biomarker levels of each group for simulation (i.e., if es=nic, then we would expect the relationship approach to be the most appropriate model)
#simdate: date to include in result output file

	set.seed(515+x)

	#Calculate meanlog and sdlog parameters for each dose-level given meanlog2, es, sigma, and nic values
	mu1 <- meanlog2 - es*sigma
	mu1.nic <- meanlog2 - nic*sigma
	biom <- 100 * exp(mu1.nic - meanlog2)
	asd <- length(mu1)

	meanlog <- mu1
	sdlog <- rep(sigma,asd)

	#Determine the true mean and sd for the Gaussian distribution
	mean.truth <- exp(meanlog + 0.5*sdlog^2)
	sd.truth <- sqrt( (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2) )

	#Simulate data
	vals <- log(as.vector( sapply(1:asd, function(x) c(rlnorm(n=nsamp[x]*pcomp[x], meanlog=meanlog[x], sdlog=sdlog[x]), rlnorm(n=nsamp[x]*(1-pcomp[x]), meanlog=meanlog[asd], sdlog=sdlog[asd])) )))
	grp <- as.vector( sapply(1:asd, function(x) rep( (1:length(meanlog))[x], nsamp[x]) ))

	#Define hyperparameters and initial starting values
	hp_normal <- list( tau.a=.001, tau.b=.001, theta.prec=0.001, mu.prec=0.001, p.alpha=1, p.beta=1 )
	inits_normal <- list( theta=1, tau=.001, p=0.5, mu_comp=0, mu_nonc=5, model=0)

	#Convert information for simulation to store in output for reference
	esx <- paste(es,collapse="/")
	nicx <- paste(nic,collapse="/")
	px <- paste(pcomp,collapse="/")

	#Mixture model assuming IND (i.e., 'ind' estimated independently)
	samp_noa <- mcmc.gaussianmix(vals=vals,grp=grp,fit='ind',nburn=nburn,nchain=nchain,hp=hp_normal,inits=inits_normal,seed=(515+x), nic=biom)
	samp_noa$mods_avg <- rep(NA, length(mu1)-1)
	res_noa <- matrix( c((515+x), 'ind',esx,nicx,meanlog2,sigma,px,samp_noa$p_comp,samp_noa$mu1,samp_noa$mu2,samp_noa$sd1,samp_noa$sd2,samp_noa$mu1_80,samp_noa$mu1_90,samp_noa$mu1_95,samp_noa$mods_avg), nrow=1)
	colnames(res_noa) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1)))
	write.table(res_noa, file = paste0('simresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

	#Mixture model assuming REL
	samp_nor <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rel',nburn=nburn,nchain=nchain,hp=hp_normal,inits=inits_normal,seed=(515+x),nic=biom)
	samp_nor$mods_avg <- rep(NA, length(mu1)-1)
	res_nor <- matrix( c((515+x), 'rel',esx,nicx,meanlog2,sigma,px,samp_nor$p_comp,samp_nor$mu1,samp_nor$mu2,samp_nor$sd1,samp_nor$sd2,samp_nor$mu1_80,samp_nor$mu1_90,samp_nor$mu1_95,samp_nor$mods_avg), nrow=1)
	colnames(res_nor) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1)))
	write.table(res_nor, file = paste0('simresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

	#Mixture model with RJMCMC and 0.5 prior on approach
	samp_no <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rjmcmc',nburn=nburn,nchain=nchain,hp=hp_normal,inits=inits_normal,prob=rep(.5,5),seed=(515+x),nic=biom)
	res_no <- matrix( c((515+x), 'rjmcmc50',esx,nicx,meanlog2,sigma,px,samp_no$p_comp,samp_no$mu1,samp_no$mu2,samp_no$sd1,samp_no$sd2,samp_no$mu1_80,samp_no$mu1_90,samp_no$mu1_95,samp_no$mods_avg), nrow=1)
	colnames(res_no) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1)))
	write.table(res_no, file = paste0('simresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

	#Mixture model with RJMCMC and 0.99 prior on approach (favors REL)
	samp_no99 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rjmcmc',nburn=nburn,nchain=nchain,hp=hp_normal,inits=inits_normal,prob=rep(.99,5),seed=(515+x),nic=biom)
	res_no99 <- matrix( c((515+x), 'rjmcmc99',esx,nicx,meanlog2,sigma,px,samp_no99$p_comp,samp_no99$mu1,samp_no99$mu2,samp_no99$sd1,samp_no99$sd2,samp_no99$mu1_80,samp_no99$mu1_90,samp_no99$mu1_95,samp_no99$mods_avg), nrow=1)
	colnames(res_no99) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1)))
	write.table(res_no99, file = paste0('simresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

	#Mixture model with RJMCMC and 0.95 prior on approach (favors REL)
	samp_no95 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rjmcmc',nburn=nburn,nchain=nchain,hp=hp_normal,inits=inits_normal,prob=rep(.95,5),seed=(515+x),nic=biom)
	res_no95 <- matrix( c((515+x), 'rjmcmc95',esx,nicx,meanlog2,sigma,px,samp_no95$p_comp,samp_no95$mu1,samp_no95$mu2,samp_no95$sd1,samp_no95$sd2,samp_no95$mu1_80,samp_no95$mu1_90,samp_no95$mu1_95,samp_no95$mods_avg), nrow=1)
	colnames(res_no95) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1)))
	write.table(res_no95, file = paste0('simresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

	#Mixture model with RJMCMC and 0.05 prior on approach (favors IND)
	samp_no05 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='rjmcmc',nburn=nburn,nchain=nchain,hp=hp_normal,inits=inits_normal,prob=rep(.05,5),seed=(515+x),nic=biom)
	res_no05 <- matrix( c((515+x), 'rjmcmc05',esx,nicx,meanlog2,sigma,px,samp_no05$p_comp,samp_no05$mu1,samp_no05$mu2,samp_no05$sd1,samp_no05$sd2,samp_no05$mu1_80,samp_no05$mu1_90,samp_no05$mu1_95,samp_no05$mods_avg), nrow=1)
	colnames(res_no05) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('sd_comp',1:asd),paste0('sd_nonc',1:asd),paste0('c80_',1:asd),paste0('c90_',1:asd),paste0('c95_',1:asd),paste0('modavg',1:(asd-1)))
	write.table(res_no05, file = paste0('simresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

}

fix_sim <- function(x, meanlog2, sigma, pcomp, nsamp, es, nic, nburn, nchain, simdate ){
###Function to fit models where dose-levels are fixed at REL or IND
#x: seed to set
#meanlog2: non-compliant component meanlog
#sigma: value to use for sdlog (lognormal)
#pcomp: proportion compliant (used to determine how many to simulate from each group)
#nsamp: sample size to simulate
#es: effect size for each level to simulate compliant group
#nic: effect size to match to biomarker levels of each group for simulation (i.e., if es=nic, then we would expect the relationship approach to be the most appropriate model)
#simdate: date to include in result output file

	set.seed(515+x)

	#Calculate meanlog and sdlog parameters for each dose-level given meanlog2, es, sigma, and nic values
	mu1 <- meanlog2 - es*sigma
	mu1.nic <- meanlog2 - nic*sigma
	biom <- 100 * exp(mu1.nic - meanlog2)
	asd <- length(mu1)

	meanlog <- mu1
	sdlog <- rep(sigma,asd)

	#Determine the true mean and sd for the Gaussian distribution
	mean.truth <- exp(meanlog + 0.5*sdlog^2)
	sd.truth <- sqrt( (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2) )

	#Simulate data
	vals <- log(as.vector( sapply(1:asd, function(x) c(rlnorm(n=nsamp[x]*pcomp[x], meanlog=meanlog[x], sdlog=sdlog[x]), rlnorm(n=nsamp[x]*(1-pcomp[x]), meanlog=meanlog[asd], sdlog=sdlog[asd])) )))
	grp <- as.vector( sapply(1:asd, function(x) rep( (1:length(meanlog))[x], nsamp[x]) ))

	#Define hyperparameters
	hp_normal <- list( tau.a=.001, tau.b=.001, theta.prec=0.001, mu.prec=0.001, p.alpha=1, p.beta=1 )

	#Convert information for simulation to store in output for reference
	esx <- paste(es,collapse="/")
	nicx <- paste(nic,collapse="/")
	px <- paste(pcomp,collapse="/")

	#Mixture model where all dose-levels are estimated with IND
 	samp_0000 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='fixed',nburn=nburn,nchain=nchain,hp=hp_normal,seed=(515+x), nic=biom, inits=list( theta=1, tau=.001, p=0.5, mu_comp=0, mu_nonc=5, model=c(0,0,0,0,1)))
	res_0000 <- matrix( c((515+x), '0000', esx,nicx,meanlog2,sigma,px,samp_0000$p_comp,samp_0000$mu1,samp_0000$mu2,samp_0000$mucomps.m,samp_0000$munoncs.m), nrow=1)
	colnames(res_0000) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('norm_mu_comp',1:asd),paste0('norm_mu_nonc',1:asd))
	write.table(res_0000, file = paste0('fixresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

 	#Mixture model where all dose-levels are estimated with REL
	samp_1111 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='fixed',nburn=nburn,nchain=nchain,hp=hp_normal,seed=(515+x), nic=biom, inits=list( theta=1, tau=.001, p=0.5, mu_comp=0, mu_nonc=5, model=c(1,1,1,1,1)))
	res_1111 <- matrix( c((515+x), '1111', esx,nicx,meanlog2,sigma,px,samp_1111$p_comp,samp_1111$mu1,samp_1111$mu2,samp_1111$mucomps.m,samp_1111$munoncs.m), nrow=1)
	colnames(res_1111) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('norm_mu_comp',1:asd),paste0('norm_mu_nonc',1:asd))
	write.table(res_1111, file = paste0('fixresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

 	#Mixture model where dose-levels are estimated by REL, IND, IND, and REL for dose-levels 1, 2, 3, and 4, respectively (dose-level 5 serves as "control")
	samp_1001 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='fixed',nburn=nburn,nchain=nchain,hp=hp_normal,seed=(515+x), nic=biom, inits=list( theta=1, tau=.001, p=0.5, mu_comp=0, mu_nonc=5, model=c(1,0,0,1,1)))
	res_1001 <- matrix( c((515+x), '1001', esx,nicx,meanlog2,sigma,px,samp_1001$p_comp,samp_1001$mu1,samp_1001$mu2,samp_1001$mucomps.m,samp_1001$munoncs.m), nrow=1)
	colnames(res_1001) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('norm_mu_comp',1:asd),paste0('norm_mu_nonc',1:asd))
	write.table(res_1001, file = paste0('fixresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

 	#Mixture model where dose-levels are estimated by REL, IND, IND, and IND for dose-levels 1, 2, 3, and 4, respectively (dose-level 5 serves as "control")
	samp_1000 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='fixed',nburn=nburn,nchain=nchain,hp=hp_normal,seed=(515+x), nic=biom, inits=list( theta=1, tau=.001, p=0.5, mu_comp=0, mu_nonc=5, model=c(1,0,0,0,1)))
	res_1000 <- matrix( c((515+x), '1000', esx,nicx,meanlog2,sigma,px,samp_1000$p_comp,samp_1000$mu1,samp_1000$mu2,samp_1000$mucomps.m,samp_1000$munoncs.m), nrow=1)
	colnames(res_1000) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('norm_mu_comp',1:asd),paste0('norm_mu_nonc',1:asd))
	write.table(res_1000, file = paste0('fixresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

 	#Mixture model where dose-levels are estimated by IND, IND, IND, and REL for dose-levels 1, 2, 3, and 4, respectively (dose-level 5 serves as "control")
	samp_0001 <- mcmc.gaussianmix(vals=vals,grp=grp,fit='fixed',nburn=nburn,nchain=nchain,hp=hp_normal,seed=(515+x), nic=biom, inits=list( theta=1, tau=.001, p=0.5, mu_comp=0, mu_nonc=5, model=c(0,0,0,1,1)))
	res_0001 <- matrix( c((515+x), '0001', esx,nicx,meanlog2,sigma,px,samp_0001$p_comp,samp_0001$mu1,samp_0001$mu2,samp_0001$mucomps.m,samp_0001$munoncs.m), nrow=1)
	colnames(res_0001) <- c('seed','fit','nic_true','nic_rel','meanlog2','sdlog','pcomp',paste0('p',1:asd),paste0('mu_comp',1:asd),paste0('mu_nonc',1:asd),paste0('norm_mu_comp',1:asd),paste0('norm_mu_nonc',1:asd))
	write.table(res_0001, file = paste0('fixresults_',simdate,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

}


########################################
## Table to format simulation results ##
########################################

sim.tab_dist <- function(dat,fituse,meanlog2_use,nic_true,sdlog_use){
###Function to create table from MCMC simulation results for different models
#dat: data frame to generate table from
#fituse: fits to include in the table (NOTE: first fit listed is what all other fits are compared to for relative SD)
#meanlog2_use: value for non-compliant mean to use in calculating bias
#sdlog_use: standard deviation value to pull results for
#nic_true: true nicotine dose-relationship level to pull results for

	result.tab <- result.tab_quantiles <- est.tab <- NULL

	for(fit in fituse){

		samp.dat <- dat[which(dat$fit==fit & dat$nic_true==nic_true & dat$meanlog2==meanlog2_use & dat$sdlog==sdlog_use),]
		j <- length(strsplit(samp.dat[['nic_true']],"/")[[1]])
		ptruth <- as.numeric(strsplit(samp.dat[['pcomp']][1],"/")[[1]])
		nictruth <- as.numeric(strsplit(samp.dat[['nic_true']][1],"/")[[1]])
		sdtruth <- as.numeric(samp.dat[['sdlog']][1])

		meanlog <- meanlog2_use - nictruth*sdtruth
		sdlog <- rep(sdtruth,j)

		c80truth <- qlnorm(.80, meanlog=meanlog, sdlog=sdlog)
		c90truth <- qlnorm(.90, meanlog=meanlog, sdlog=sdlog)
		c95truth <- qlnorm(.95, meanlog=meanlog, sdlog=sdlog)
		mean_comp.truth <- coluse.true <- exp(meanlog + 0.5*sdlog^2)
		mean_nonc.truth <- exp(meanlog + 0.5*sdlog^2)[j]
		sd_comp.truth <- sqrt( (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2) )		
		sd_nonc.truth <- sqrt( (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2) )[j]

		samp <- suppressWarnings(as.matrix( samp.dat[,!is.na(as.numeric(samp.dat[1,]))] )) #extract only numeric columns, supresses warnings created with non-numeric columns

		#Extract column means and standard deviations
		samp.m <- colMeans(samp)
		samp.s <- colSds(samp); names(samp.s) <- names(samp.m) #add names to SD vector

		if(fit==fituse[1]){
			sd.ref_p <- samp.s[ paste0('p',1:j) ]
			sd.ref_mu_comp <- samp.s[ paste0('mu_comp', 1:j) ] 
			sd.ref_mu_nonc <- samp.s[ paste0('mu_nonc', 1:j) ] 
			sd.ref_sd_comp <- samp.s[ paste0('sd_comp', 1:j) ]
			sd.ref_sd_nonc <- samp.s[ paste0('sd_nonc', 1:j) ]
			sd.ref_c80 <- samp.s[ paste0('c80_', 1:j) ]
			sd.ref_c90 <- samp.s[ paste0('c90_', 1:j) ]
			sd.ref_c95 <- samp.s[ paste0('c95_', 1:j) ]
		}

		#Extract estimate proportion for models (note that we should have p1=1-p2 and p2=1-p1)
		mean.p_comp <- samp.m[ paste0('p',1:j) ] #proportion in complier group
		bias.p_comp <- mean.p_comp - ptruth
		sd.p_comp <- samp.s[ paste0('p',1:j) ]
		mse.p_comp <- sd.p_comp^2 + bias.p_comp^2
		if(fit!=fituse[1]){ sd.p_comp <- samp.s[ paste0('p',1:j) ] / sd.ref_p }
		p_comp_tab  <- round( cbind( bias.p_comp,sd.p_comp,mse.p_comp), 2)

		#Compliant group mean
		mean.mu_comp <- samp.m[ paste0('mu_comp', 1:j) ] 
		bias.mu_comp <- mean.mu_comp - mean_comp.truth
		sd.mu_comp <- samp.s[ paste0('mu_comp', 1:j) ] 
		mse.mu_comp <- sd.mu_comp^2 + bias.mu_comp^2
		if(fit!=fituse[1]){ sd.mu_comp <- samp.s[ paste0('mu_comp', 1:j) ] / sd.ref_mu_comp } 
		mu_comp_tab  <- round(cbind( bias.mu_comp, sd.mu_comp, mse.mu_comp ),2)

		#Non-compliant group mean
		mean.mu_nonc <- samp.m[ paste0('mu_nonc', 1:j) ] 
		bias.mu_nonc <- mean.mu_nonc - mean_nonc.truth
		sd.mu_nonc <- samp.s[ paste0('mu_nonc', 1:j) ] 
		mse.mu_nonc <- sd.mu_nonc^2 + bias.mu_nonc^2
		if(fit!=fituse[1]){ sd.mu_nonc <- samp.s[ paste0('mu_nonc', 1:j) ] / sd.ref_mu_nonc } 
		mu_nonc_tab  <- round( cbind( bias.mu_nonc, sd.mu_nonc, mse.mu_nonc),2)

		#Standard deviation for compliant component
		mean.sd_comp <- samp.m[ paste0('sd_comp', 1:j) ]
		bias.sd_comp <- mean.sd_comp - sd_comp.truth
		sd.sd_comp <- samp.s[ paste0('sd_comp', 1:j) ]
		mse.sd_comp <- bias.sd_comp^2 + sd.sd_comp^2
		if(fit!=fituse[1]){ sd.sd_comp <- samp.s[ paste0('sd_comp', 1:j) ] / sd.ref_sd_comp } 
		sd_comp_tab <- round( cbind(bias.sd_comp, sd.sd_comp, mse.sd_comp), 2)

		#Standard deviation for noncompliant component
		mean.sd_nonc <- samp.m[ paste0('sd_nonc', 1:j) ]
		bias.sd_nonc <- mean.sd_nonc - sd_nonc.truth
		sd.sd_nonc <- samp.s[ paste0('sd_nonc', 1:j) ]
		mse.sd_nonc <- bias.sd_nonc^2 + sd.sd_nonc^2
		if(fit!=fituse[1]){ sd.sd_nonc <- samp.s[ paste0('sd_nonc', 1:j) ] / sd.ref_sd_nonc } 
		sd_nonc_tab <- round( cbind(bias.sd_nonc, sd.sd_nonc, mse.sd_nonc), 2)

		#80th percentile of compliant group mean
		mean.c80 <- samp.m[ paste0('c80_', 1:j) ]
		bias.c80 <- mean.c80 - c80truth
		sd.c80 <- samp.s[ paste0('c80_', 1:j) ]
		mse.c80 <- bias.c80^2 + sd.c80^2
		if(fit!=fituse[1]){ sd.c80 <- samp.s[ paste0('c80_', 1:j) ] / sd.ref_c80 } 
		c80_tab <- round( cbind(bias.c80, sd.c80, mse.c80), 2)

		#90th percentile of compliant group mean
		mean.c90 <- samp.m[ paste0('c90_', 1:j) ]
		bias.c90 <- mean.c90 - c90truth
		sd.c90 <- samp.s[ paste0('c90_', 1:j) ]
		mse.c90 <- bias.c90^2 + sd.c90^2
		if(fit!=fituse[1]){ sd.c90 <- samp.s[ paste0('c90_', 1:j) ] / sd.ref_c90 } 
		c90_tab <- round( cbind(bias.c90, sd.c90, mse.c90), 2)

		#95th percentile of compliant group mean
		mean.c95 <- samp.m[ paste0('c95_', 1:j) ]
		bias.c95 <- mean.c95 - c95truth
		sd.c95 <- samp.s[ paste0('c95_', 1:j) ]
		mse.c95 <- bias.c95^2 + sd.c95^2
		if(fit!=fituse[1]){ sd.c95 <- samp.s[ paste0('c95_', 1:j) ] / sd.ref_c95 } 
		c95_tab <- round( cbind(bias.c95, sd.c95, mse.c95), 2)

		#RJMCMC results
		if(fit %in% c('rjmcmc05','rjmcmc50','rjmcmc95')){
			mean.rj <- samp.m[ paste0('modavg', 1:(j-1)) ]
			sd.rj <- samp.s[ paste0('modavg', 1:(j-1)) ]
			rj_tab <- c( paste0( round(mean.rj,2), ' (', round(sd.rj,2),')'), '-')
		}else{
			rj_tab <- c(rep('-',5))
		}

		roundtruth <- round(cbind(mean_comp.truth, mean.p_comp, mean.mu_comp, mean.mu_nonc, mean.sd_comp, mean.sd_nonc, mean.c80, mean.c90, mean.c95),3)
		est.tab <- rbind(est.tab, cbind(fit, roundtruth))

		result.tab <- rbind(result.tab, cbind( fit, nictruth, rbind(p_comp_tab[1:(j-1),],rep('-',3)), mu_comp_tab, rbind(mu_nonc_tab[1:(j-1),],rep('-',3)), sd_comp_tab, rbind(sd_nonc_tab[1:(j-1),],rep('-',3)), rj_tab) )
		result.tab_quantiles <- rbind(result.tab_quantiles, cbind( fit, nictruth, rbind(c80_tab[1:(j-1),],rep('-',3)), rbind(c90_tab[1:(j-1),],rep('-',3)), rbind(c95_tab[1:(j-1),],rep('-',3)), rj_tab)[-j,] )

	}

	colnames(result.tab) <- c('fit','ES',
		paste0('pc.', c('bias','sd','mse')),
		paste0('mu_comp.', c('bias','sd','mse')),
		paste0('mu_nonc.', c('bias','sd','mse')),
		paste0('sd_comp.', c('bias','sd','mse')),
		paste0('sd_nonc.', c('bias','sd','mse')),
		'prop_in_rel' )

	colnames(result.tab_quantiles) <- c('fit','ES',
		paste0('c80.', c('bias','sd','mse')),
		paste0('c90.', c('bias','sd','mse')),
		paste0('c95.', c('bias','sd','mse')),
		'prop_in_rel' )

	colnames(est.tab) <- c('fit','truemu_comp','p.est','mu_comp.est','mu_nonc.est','sd_comp.est','sd_nonc.est','c80.est','c90.est','c95.est')

	restab <- xtable(result.tab)
	restab_q <- xtable(result.tab_quantiles)

	ret <- list(res=restab, resq=restab_q, mat=result.tab, matq=result.tab_quantiles, est.tab=est.tab)

	return(ret)

}


sim.plot <- function(dat,fituse,meanlog2_use,nic_scen,sdlog_use,doses,dose.ch,out,mainlab,singleplot=F){
###Function to create figures from MCMC simulation results for different models
#dat: data frame to generate figure from
#fituse: fits to include in the table (NOTE: first fit listed is what all other fits are compared to for relative SD)
#meanlog2_use: value for non-compliant mean to use in calculating bias
#sdlog_use: standard deviation value to pull results for
#nic_scen: scenario to pull nicotine results for
#doses: the value for the ES for the different scenarios used (helps to create matrix and plots)
#dose.ch: identify which dose-level is changing [NOTE: 1=lowest dose level, J=highest dose level]
#out: the name of the object you want to plot
#mainlab: title for the figure(s)
#singleplot: if TRUE provide a single plot for each dose level, otherwise only 2x2 figure is produced

	mse.tab_1 <- mse.tab_2 <- mse.tab_3 <- mse.tab_4 <- mse.tab_5 <- matrix( nrow=length(fituse), ncol=length(doses), dimnames=list(fituse,doses)) #create matrices to store results for each dose level

	for(nic_true in nic_scen){
	for(fit in fituse){

		samp.dat <- dat[which(dat$fit==fit & dat$nic_true==nic_true & dat$meanlog2==meanlog2_use & dat$sdlog==sdlog_use),]
		j <- length(strsplit(samp.dat[['nic_true']],"/")[[1]])
		ptruth <- as.numeric(strsplit(samp.dat[['pcomp']][1],"/")[[1]])
		nictruth <- as.numeric(strsplit(samp.dat[['nic_true']][1],"/")[[1]])
		sdtruth <- as.numeric(samp.dat[['sdlog']][1])

		meanlog <- meanlog2_use - nictruth*sdtruth
		sdlog <- rep(sdtruth,j)

		c80truth <- qlnorm(.80, meanlog=meanlog, sdlog=sdlog)
		c90truth <- qlnorm(.90, meanlog=meanlog, sdlog=sdlog)
		c95truth <- qlnorm(.95, meanlog=meanlog, sdlog=sdlog)
		mean_comp.truth <- coluse.true <- exp(meanlog + 0.5*sdlog^2)
		mean_nonc.truth <- exp(meanlog + 0.5*sdlog^2)[j]
		sd_comp.truth <- sqrt( (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2) )		
		sd_nonc.truth <- sqrt( (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2) )[j]

		samp <- suppressWarnings(as.matrix( samp.dat[,!is.na(as.numeric(samp.dat[1,]))] )) #extract only numeric columns, supresses warnings created with non-numeric columns

		#Extract column means and standard deviations
		samp.m <- colMeans(samp)
		samp.s <- colSds(samp); names(samp.s) <- names(samp.m) #add names to SD vector

		#Extract estimate proportion for models (note that we should have p1=1-p2 and p2=1-p1)
		mean.p_comp <- samp.m[ paste0('p',1:j) ] #proportion in complier group
		bias.p_comp <- mean.p_comp - ptruth
		sd.p_comp <- samp.s[ paste0('p',1:j) ]
		mse.p_comp <- sd.p_comp^2 + bias.p_comp^2
		p_comp_tab  <- round( cbind( bias.p_comp,sd.p_comp,mse.p_comp), 2)

		#Compliant group mean
		mean.mu_comp <- samp.m[ paste0('mu_comp', 1:j) ] 
		bias.mu_comp <- mean.mu_comp - mean_comp.truth
		sd.mu_comp <- samp.s[ paste0('mu_comp', 1:j) ] 
		mse.mu_comp <- sd.mu_comp^2 + bias.mu_comp^2
		mu_comp_tab  <- round(cbind( bias.mu_comp, sd.mu_comp, mse.mu_comp ),2)

		#Non-compliant group mean
		mean.mu_nonc <- samp.m[ paste0('mu_nonc', 1:j) ] 
		bias.mu_nonc <- mean.mu_nonc - mean_nonc.truth
		sd.mu_nonc <- samp.s[ paste0('mu_nonc', 1:j) ] 
		mse.mu_nonc <- sd.mu_nonc^2 + bias.mu_nonc^2
		mu_nonc_tab  <- round( cbind( bias.mu_nonc, sd.mu_nonc, mse.mu_nonc),2)

		#Standard deviation for compliant component
		mean.sd_comp <- samp.m[ paste0('sd_comp', 1:j) ]
		bias.sd_comp <- mean.sd_comp - sd_comp.truth
		sd.sd_comp <- samp.s[ paste0('sd_comp', 1:j) ]
		mse.sd_comp <- bias.sd_comp^2 + sd.sd_comp^2
		sd_comp_tab <- round( cbind(bias.sd_comp, sd.sd_comp, mse.sd_comp), 2)

		#Standard deviation for noncompliant component
		mean.sd_nonc <- samp.m[ paste0('sd_nonc', 1:j) ]
		bias.sd_nonc <- mean.sd_nonc - sd_nonc.truth
		sd.sd_nonc <- samp.s[ paste0('sd_nonc', 1:j) ]
		mse.sd_nonc <- bias.sd_nonc^2 + sd.sd_nonc^2
		sd_nonc_tab <- round( cbind(bias.sd_nonc, sd.sd_nonc, mse.sd_nonc), 2)

		#80th percentile of compliant group mean
		mean.c80 <- samp.m[ paste0('c80_', 1:j) ]
		bias.c80 <- mean.c80 - c80truth
		sd.c80 <- samp.s[ paste0('c80_', 1:j) ]
		mse.c80 <- bias.c80^2 + sd.c80^2
		c80_tab <- round( cbind(bias.c80, sd.c80, mse.c80), 2)

		#90th percentile of compliant group mean
		mean.c90 <- samp.m[ paste0('c90_', 1:j) ]
		bias.c90 <- mean.c90 - c90truth
		sd.c90 <- samp.s[ paste0('c90_', 1:j) ]
		mse.c90 <- bias.c90^2 + sd.c90^2
		c90_tab <- round( cbind(bias.c90, sd.c90, mse.c90), 2)

		#95th percentile of compliant group mean
		mean.c95 <- samp.m[ paste0('c95_', 1:j) ]
		bias.c95 <- mean.c95 - c95truth
		sd.c95 <- samp.s[ paste0('c95_', 1:j) ]
		mse.c95 <- bias.c95^2 + sd.c95^2
		c95_tab <- round( cbind(bias.c95, sd.c95, mse.c95), 2)

		mse.tab_1[fit,paste0(nictruth[dose.ch])] <- get(out)[1]
		mse.tab_2[fit,paste0(nictruth[dose.ch])] <- get(out)[2]
		mse.tab_3[fit,paste0(nictruth[dose.ch])] <- get(out)[3]
		mse.tab_4[fit,paste0(nictruth[dose.ch])] <- get(out)[4]
		mse.tab_5[fit,paste0(nictruth[dose.ch])] <- get(out)[5]

	}
	}

	p1.mat <- cbind(t(mse.tab_1) / mse.tab_1[1,])
	p2.mat <- cbind(t(mse.tab_2) / mse.tab_2[1,])
	p3.mat <- cbind(t(mse.tab_3) / mse.tab_3[1,])
	p4.mat <- cbind(t(mse.tab_4) / mse.tab_4[1,])
	p5.mat <- cbind(t(mse.tab_5) / mse.tab_5[1,])

	###Plot one 2x2 figure comparing dose levels 1, 2, 3, and 4
	par(mfrow=c(2,2),oma=c(4,1,1,1))
	matplot( p1.mat, main='Dose Level 1', type=c('b'), pch=1, lty=1:6, col=1:6,ylim=c(0,5), ylab='MSE Relative to the IND Model', xlab=paste0('Dose-Level ',dose.ch,' ES'), xaxt='n'); axis(1, at=1:length(doses), doses)
	matplot( p2.mat, main='Dose Level 2', type=c('b'), pch=1, lty=1:6, col=1:6,ylim=c(0,5), ylab='MSE Relative to the IND Model', xlab=paste0('Dose-Level ',dose.ch,' ES'), xaxt='n'); axis(1, at=1:length(doses), doses)
	matplot( p3.mat, main='Dose Level 3', type=c('b'), pch=1, lty=1:6, col=1:6,ylim=c(0,5), ylab='MSE Relative to the IND Model', xlab=paste0('Dose-Level ',dose.ch,' ES'), xaxt='n'); axis(1, at=1:length(doses), doses)
	matplot( p4.mat, main='Dose Level 4', type=c('b'), pch=1, lty=1:6, col=1:6,ylim=c(0,5), ylab='MSE Relative to the IND Model', xlab=paste0('Dose-Level ',dose.ch,' ES'), xaxt='n'); axis(1, at=1:length(doses), doses)
	#matplot( p5.mat, main='Dose Level 5', type=c('b'), pch=1, col=1:6,ylim=c(0,2), ylab='Measure Relative to the IND Model', xlab=paste0('Dose-Level ',dose.ch,' ES'), xaxt='n'); axis(1, at=1:length(doses), doses)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend('bottom',xpd=T,horiz=T,inset = c(0,0),lty=1:6,legend=c('IND','REL',expression('RJ'['05']),expression('RJ'['50']),expression('RJ'['95']),expression('RJ'['99'])), col=1:6, pch=1, lwd=1, bty='n', cex=1.4)
	#legend('bottom',xpd=T,horiz=T,inset = c(0,0),legend=c('IND','REL',expression('RJ'['05']),expression('RJ'['50']),expression('RJ'['95'])), lty=1:5, col=1:5, pch=1, lwd=1, bty='n', cex=1.4)
	title(main=mainlab, line=-1.25)

	subtitle <- nictruth ;subtitle[dose.ch] <- 'X'
	mtext(paste0('ES: ',paste0(subtitle,collapse='/')), line=-2.75)

	###Plot figures for each dose level
	ylim.list <- sapply(1:length(nictruth), function(x) NULL)
	ylim.list[[dose.ch]] <- c(0,5)

	if(singleplot==T){
		for(m in 1:5){
			dev.new()
			par(oma=c(2,0,0,0))
			matplot( get(paste0('p',m,'.mat')), ylim=ylim.list[[m]], main=paste0('Dose Level ',m), type=c('b'), pch=1, col=1:6, ylab='Measure Relative to the IND Model', xlab=paste0('Dose-Level ',dose.ch,' ES'), xaxt='n'); axis(1, at=1:length(doses), doses)
			mtext(mainlab, line=0.25)
			par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
			plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
			legend('bottom',xpd=T,horiz=T,inset = c(0,0),lty=1:6, legend=c('IND','REL',expression('RJ'['05']),expression('RJ'['50']),expression('RJ'['95']),expression('RJ'['99'])), col=1:6, pch=1, lwd=1, bty='n', cex=1)
		}
	}

	ret <- list(dose1=mse.tab_1, dose2=mse.tab_2, dose3=mse.tab_3, dose4=mse.tab_4, dose5=mse.tab_5)

	return(ret)

}


fix.plot <- function(dat,fit){
###Function to create table from MCMC simulation results for different models from fix_sim for 5 fits and the simulated truth/proposed relationship
#dat: data frame to generate table from
#fit: the assumed independence/relative relationship for each level to include in the table/plots

	plot.list <- rep( list(NA), length(fit) )
	names(plot.list) <- fit

	#Extract scenario parameters to determine 
	meanlog_noncomp <- as.numeric(dat[['meanlog2']][1])
	nictruth <- as.numeric(strsplit(dat[['nic_true']][1],"/")[[1]])
	nicrel <- as.numeric(strsplit(dat[['nic_rel']][1],"/")[[1]])
	sdtruth <- as.numeric(dat[['sdlog']][1])

	#Determine the true meanlog for the simulation and the resulting meanlog for assuming the relationship is true
	meanlog <- meanlog_noncomp - nictruth*sdtruth
	rellog <- meanlog_noncomp - nicrel*sdtruth

	for(j in fit){
		samp <- dat[which(dat$fit==j),c(paste0('norm_mu_comp',1:5),paste0('norm_mu_nonc',1:4))]

		#Extract column means and standard deviations
		samp.m <- colMeans(samp)
		samp.hpd <- sapply(1:ncol(samp), function(x) boa.hpd(samp[,x], alpha=0.05) )
		samp.mat <- t( rbind( samp.m, samp.hpd ) )
		rownames(samp.mat) <- c('d1c','d2c','d3c','d4c','d5','d1n','d2n','d3n','d4n')
		colnames(samp.mat) <- c('est','lhpd','uhpd')

		plot.list[[paste0(j)]] <- samp.mat
	}


	ci.plot <- function( sm, xs=0, pch.use=1 ){
	###Plot to add points and HPD interval
	#sm is samp.mat data, xs is how much to shift x-axis point by
		points( sm[1:5,'est'], x=(1:5)+xs, pch=pch.use )
		lines( x=t(matrix((1:5)+xs,nrow=5,ncol=3)), y=t(cbind(sm[1:5,c('lhpd','uhpd')],NA)) )

		points( sm[6:9,'est'], x=(1:4)+xs, pch=pch.use )
		lines( x=t(matrix((1:4)+xs,nrow=4,ncol=3)), y=t(cbind(sm[6:9,c('lhpd','uhpd')],NA)) )
	}

	#Add points to figure
	plot(xlim=c(0.5,5.5), ylim=c(0,5), xlab='Randomized Group', ylab='Biomarker Level',x=-10,y=-10)
	points(x=(1:5)-0.3, y=meanlog, pch=3, lwd=2)
	points(x=(1:4)-0.3, y=rep(4,4), pch=3, lwd=2)
	points(x=(1:4)-0.3, y=rellog[-5], pch=4, lwd=2)
	for(i in 1:length(plot.list)){ ci.plot(sm=plot.list[[i]], xs=c(-.15,-.05,.05,.15,.25)[i], pch.use=c(0,1,2,5,6)[i] ) }
	legtext <- c('True Mean','Predicted Mean',expression(paste(xi,"=(0,0,0,0)")),expression(paste(xi,"=(1,0,0,0)")),expression(paste(xi,"=(1,0,0,1)")),expression(paste(xi,"=(1,1,1,1)")),expression(paste(xi,"=(0,0,0,1)")),'95% HPD Interval')
	legend('bottomright', bty='n', pch=c(3,4,0,1,2,5,6,NA), lty=c(rep(NA,7),1), legend=legtext )

	ret <- list(plot.list=plot.list)

	return(ret)

}


base10toK <- function(x,K=2,size=NA){
###Function to convert base 10 number to vector representing base K value
#x: base 10 value
#K: base to convert to
#size: if provided, provide base K representation with leading 0's (e.g., base 2: 0010 vs. 10)

	if( is.na(size)==FALSE ){
		num <- x
		remainder <- NA
		numK <- rep(NA, size)
		for( i in size:1 ){
			remainder <- num %% K
			num <- num %/% K
			numK[ i ] <- remainder
		}
		if( num != 0 ){print("Size not long enough to represent base K number. Answer doesn't include leading 0's."); size <- NA}
	}

	if( is.na(size)==TRUE ){
		num <- x
		remainder <- NA
		numK <- NULL
		while( num != 0 ){
			remainder <- num %% K
			num <- num %/% K
			numK <- c(remainder, numK)
		}
	}
	
	return( numK )
}


baseKto10 <- function(x,K=2){
###Function to convert base K vector of values into base 10 number
#x: vector with values in base K (e.g., base 2 has 0/1)
#K: base to convert from, default is base 2
	sum(x * K^((length(x)-1):0))
}


pairwiseZmat <- function(z,zval){
###Function to take vector of compliance and convert to matrices with pairwise comparisons
#z: vector of compliance for each y
#zval: returns pairwise matrix for vector for values of z equal to zval

	asd <- as.matrix(abs(dist(z,diag=T)-1)) #determine which pairwise z values are equal (i.e., both compliant or both non-compliant)
	asd[,which(z!=zval)] <- 0 #change all columns that aren't zval to 0
	return(asd)
}


round2 <- function(x,digits=2){
###Function to round a given value to have 'digits' decimal places (default is 2)
	format(round(x, digits), nsmall = digits)
}


boa.hpd <- function(x, alpha){
###Function to calculate HPD interval
    n <- length(x)
    m <- max(1, ceiling(alpha * n))
    y <- sort(x)
    a <- y[1:m]
    b <- y[(n - m + 1):n]
    i <- order(b - a)[1]
    structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}


compliance_fig <- function(tab1, tab2, plot.title='Estimated Compliance as a Function of TNE', export=F, bio.lab=c(0.4,1.3,2.4,5.2), gray=F, xmax){
###Function to create compliance figure on non-log scale based on log.tab output
#tab1: data from all components (IND) 'tab' output from mcmc.gaussianmix function
#tab2: data from RJMCMC 'tab' output from mcmc.gaussianmix function
#plot.title: title of plot
#export: if TRUE, export to PDF file
#bio.lab: biomarker levels to use for labels
#gray: use gray-scale or colors
#xmax: value to use for upper limit of xlim in plot

	###Estimated probability of compliance
	if(export==T){ pdf(paste0('compliance_allVrj_',datestr(),'.pdf')) }else{ dev.new() }
	tne.grid <- seq(0.1,xmax,by=.2) #generate grid of TNE values to plot
	j <- length(bio.lab)

	if(gray==T & j==4){col1 <- col2 <- c('black','gray25','gray50','gray75')}
	if(gray==F & j==4){col1 <- col2 <- c('orangered2','darkgreen','blue','darkmagenta')}
	if(gray==T & j==3){col1 <- col2 <- c('black','gray25','gray75')}
	if(gray==F & j==3){col1 <- col2 <- c('orangered2','darkgreen','darkmagenta')}

	plot(xlim=c(0,xmax),ylim=c(0,1),xlab='TNE (nmol/mL)',ylab='Pr(C=1)',x=-100,y=-100,main=plot.title, cex.axis=1.4, cex.main=1.4, cex.lab=1.4)	
	for(aa in 1:j){
		pc.grid <- ( tab1[aa,'p_comp']*dnorm(log(tne.grid), mean=tab1[aa,'mean.mu_comp'], sd=tab1[aa,'sd']) ) / ( tab1[aa,'p_comp']*dnorm(log(tne.grid), mean=tab1[aa,'mean.mu_comp'], sd=tab1[aa,'sd']) + (1-tab1[aa,'p_comp'])*dnorm(log(tne.grid), mean=tab1[aa,'mean.mu_nonc'], sd=tab1[aa,'sd']) )
		lines(x=tne.grid, y=pc.grid,lty=1,col=col1[aa], lwd=2)

		pc.grid2 <- ( tab2[aa,'p_comp']*dnorm(log(tne.grid), mean=tab2[aa,'mean.mu_comp'], sd=tab2[aa,'sd']) ) / ( tab2[aa,'p_comp']*dnorm(log(tne.grid), mean=tab2[aa,'mean.mu_comp'], sd=tab2[aa,'sd']) + (1-tab2[aa,'p_comp'])*dnorm(log(tne.grid), mean=tab2[aa,'mean.mu_nonc'], sd=tab2[aa,'sd']) )
		lines(x=tne.grid, y=pc.grid2,lty=4,col=col2[aa], lwd=2)

		if(bio.lab[aa]==0.4){
			text(x=tne.grid[50],y=pc.grid[50],'0.4', pos=4, col=col1[aa], cex=.8)
			text(x=tne.grid[56],y=pc.grid[55],'IND', pos=4, col=col1[aa], cex=.8)
			text(x=tne.grid[20],y=pc.grid2[17],'0.4', pos=2, col=col2[aa], cex=.8)
			text(x=tne.grid[23],y=pc.grid2[18],expression('RJ'[95]), pos=2, col=col2[aa], cex=.8)
		}
		if(bio.lab[aa]==1.3){
			text(x=tne.grid[46],y=pc.grid[48],'1.3', pos=4, col=col1[aa],cex=.8)
			text(x=tne.grid[50],y=pc.grid[52],'IND', pos=4, col=col1[aa],cex=.8)
			text(x=tne.grid[38],y=pc.grid2[40],expression('1.3 RJ'[95]), pos=4, col=col2[aa],cex=.8)
		}
		if(bio.lab[aa]==2.4){
			text(x=tne.grid[16],y=pc.grid[15],'2.4 IND', pos=2, col=col1[aa],cex=.8)
			text(x=tne.grid[120],y=pc.grid2[120],expression('2.4 RJ'[95]), pos=3, col=col2[aa],cex=.8)
		}
		if(bio.lab[aa]==5.2){
			text(x=tne.grid[101],y=pc.grid[101],'5.2 IND', pos=2, col=col1[aa],cex=.8)
			text(x=tne.grid[101],y=pc.grid2[101],expression('5.2 RJ'[95]), pos=4, col=col2[aa],cex=.8)}
	}

	legend(x=xmax-xmax*0.37, y=1.02, bty='n', legend=bio.lab[1:j], lty=1, lwd=2, col=col1, title='IND', cex=1.2)
	legend(x=xmax-xmax*0.17, y=1.02, bty='n', legend=bio.lab[1:j], lty=4, lwd=2, col=col2, title=expression('RJ'[95]), cex=1.2)

	if(export==T){ dev.off() }

}


gelman.rubin <- function(out){
###Function to calculate the Gelman-Rubin (1992) potential scale reduction factor (PSRF) across multiple chains
#out: list of MCMC output to analyze from M chains

	M <- length(out)
	n <- length(out[[1]])
	theta_m_star <- sapply(1:M, function(x) mean(out[[x]]) )
	var_m <- sapply(1:M, function(x) var(out[[x]]) )
	theta_star_star <- sum(theta_m_star) / M

	B <- (n / (M-1)) * sum( (theta_m_star - theta_star_star)^2 )  #between-chain variance
	W <- (1 / M) * sum( var_m )  #within-chain variance
	Vhat <- ((n-1)/n)*W + (1/n)*B
	
	psrf <- sqrt( Vhat / W )
	return(psrf)
}


autocor <- function(out,h.max=50,nburn){
###Function to calculate autocorrelation within chains
#out: list of MCMC output to analyze from M chains
#h.max: largest lag value to go out to (default of 50)
#nburn: burn-in period length

	M <- length(out)
	n <- length(out[[1]])
	rho <- lapply(1:M, function(m) sapply(1:h.max, function(h) gammah(h=h, theta=out[[m]][nburn:n] )))

	return(rho)
}	


gammah <- function(h,theta){ 
###Function for use within autocor to calculate gamma function
#h: lag length
#theta: output from chain to calculate function

	n <- length(theta)
	t.m <- mean(theta)
	gamh <- (1/(n-h)) * sum( (theta[(h+1):n] - t.m)*(theta[1:(n-h)] - t.m))
	gam0 <- (1/n) * sum( (theta[1:n] - t.m)*(theta[1:n] - t.m))
	rhoh <- (gamh / gam0)
	return(rhoh)
}


brooks.giudici <- function(out){
###Function to calculate the extension to the Gelman-Rubin diagnostic proposed by Brooks and Giudici (2000)
#out: list of MCMC output to analyze from M chains excluding the burn-in period
	
	I <- length(out) #total number of chains
	for(chain in 1:I){
		out[[chain]][,'model'] <- out[[chain]][,'model']+1 #augment model counts to be from 1:x instead of 0:(x-1)
	}

	T <- nrow(out[[1]]) #total iterations/time
	mods <- 1:16 #all possible model values
	M <- length(mods) #total number of models
	
	Kim <- lapply(1:I, function(x) sapply(mods, function(y) sum(out[[x]][,'model'] == y) ))
	Ki <- sapply(1:I, function(x) sum(Kim[[x]]))
	Km <- sapply(1:M, function(x) sum(do.call(rbind, Kim)[,x]) )

	outc <- do.call(rbind, out) #combine all chains together in one data set for some calculations
	theta_star.m_star <- sapply(1:M, function(x) (1 / Km[x]) * sum(outc[which(outc[,'model']==x),'dev']) )
	theta_i.star_star <- sapply(1:I, function(x) (1 / Ki[x]) * sum(out[[x]][,'dev']) )
	theta_i.m_star <- lapply(1:I, function(x) sapply(1:M, function(y) (1/Kim[[x]][y])*sum(out[[x]][which(out[[x]][,'model']==y),'dev']) ) )
	theta_star.star_star <- (1/(I*T)) * sum(outc[,'dev'])
	theta_i_star <- (1/T) * sapply(1:I, function(x) sum( out[[x]][,'dev'] ) )
	theta_star_star <- (1/I) * sum(theta_i_star)

	Vhat <- (1 / (I*T - 1)) * sum( (outc[,'dev'] - theta_star_star)^2 )
	Bm <- sum((theta_star.m_star - theta_star.star_star)^2 / (M-1), na.rm=T)
	Wm <- (1/M) *  sum(na.rm=T,unlist(lapply(1:I, function(x) sapply(1:M, function(y) ((out[[x]][which(out[[x]][,'model']==y),'dev'] - theta_star.m_star[y])^2)/(I*Km[y] - 1) )) ))
	Wc <- (1/I) *  sum(na.rm=T,unlist(lapply(1:I, function(x) sapply(1:M, function(y) ((out[[x]][which(out[[x]][,'model']==y),'dev'] - theta_i.star_star[x])^2)/(M*Ki[x] - 1) )) ))
	BmWc <- sum(unlist(lapply(1:I, function(x) sapply(1:M, function(y) ((theta_i.m_star[[x]][y] - theta_i.star_star[x])^2) / (I*(M-1)) ))), na.rm=T)
	WmWc <- (1/(I*M)) * sum(na.rm=T,unlist(lapply(1:I, function(x) sapply(1:M, function(y) ((out[[x]][which(out[[x]][,'model']==y),'dev'] - theta_i.m_star[[x]][y])^2)/(Kim[[x]][y] - 1) )) ))

	ret <- c(Vhat,Wc,Wm,WmWc,Bm,BmWc)
	names(ret) <- c('Vhat','Wc','Wm','WmWc','Bm','BmWc')
	return(ret)
}

