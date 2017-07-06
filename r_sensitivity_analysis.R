# require(deSolve)
# require(bbmle)
# source("data.R")
# hostk2fit1 <- 0.1330609  # from fig 2.R
# hostk2fit2 <- 6.568825e-02  # from fig 2.R
# #************************DATA**************************************
# tsetse <- dat[,c(1:4,11)] 
# names(tsetse) <- c("year","month","tsetse_over_time","lusulu","time")
# #***********Likelihood function********************
# nll.pois <- function(parms=tsetse_params(),dat=tsetse){ 
#   if(parms$pf > 1 | parms$pf <= 0.001 | parms$mu.b > 0.1 ) {
#     ll <- -100000
#   }
#   if(parms$pf <= 1 & parms$pf > 0.001 & parms$mu.b <= 0.1 ) {
#     initial <- c(H=as.numeric(parms$hosts.zero), P=as.numeric(parms$pupae.zero), A=as.numeric(parms$adults.zero)) # initial conditions
#     simulate <- simPop(parms=parms) # model output 
#     time <- dat$time
#     num.tsetse <- round(dat$tsetse_over_time,0) # data 
#     est <- round(simulate[simulate$time %in% time,4],0) # model output at timepoints have data for
#     ll <- sum(dpois(x=num.tsetse+1,lambda=est+1,log=T))
#   }
#   if (ll == "NaN") {
#     ll <- -100000
#   }
#   if (ll == Inf ){
#     ll <- -100000
#   }
#   return(-ll)
# }
# #********enable multiple parameters to be fit and substituted back in parameter list****
# subsParms <- function(fit.params, fixed.params=tsetse_params())
#   within(fixed.params, {
#     loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
#     unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]        
#     for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
#     for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
#     rm(nm, loggedParms, unloggedParms)
#   }) ## Make likelihood a function of fixed and fitted parameters.
# objFXN <- function(fit.params ## paramters to fit
#                    , fixed.params =tsetse_params() ## fixed paramters
#                    , dat=tsetse) {
#   parms <- subsParms(fit.params, fixed.params)
#   nll.pois(parms, dat = dat) ## then call likelihood
# }
# #********************PARAMETERS AND CONDITIONS*************************************
# tseqDay <- tsetse$time                     ## times to solve at
# tsetse_params <- function(  area = 541     ## total area 
#                             , larv = 1/11  ## rate larvae are produced
#                             , pup = 1/45   ## rate pupae emerge as adults
#                             , mu.b = 0.015 ## female adult background death rate 
#                             , mu.p = 0.006 ## pupal initial mortality rate
#                             , mud.p = 0 
#                             , pf = 1 ## average probability of not feeding on any given day during feeding cycle given one host present in 1km2
#                             , S = 6 ## number of days a fly can go without feeding before it starve
#                             , h.r = 0.007 / 30 # host growth rate
#                             , hk_2 = as.numeric(hostk2fit2 / 30) #rate is fit to months and need to approx. for days
#                             , hk_1 = as.numeric(hostk2fit1 / 30)
#                             , hosts.zero = 2100
#                             , adults.zero = 250
#                             , pupae.zero = 250
# )
# return(as.list(environment()))
# initial <- c(H=tsetse_params()$hosts.zero, P=tsetse_params()$pupae.zero, A=tsetse_params()$adults.zero) # initial conditions
# #********************POPULATION DYNAMICS FUNCTION***********************************
# tsetse_mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
#   if(tt < 1125){
#     lambda <- H / area  
#   }
#   if(tt >= 1125 & tt < 1605){
#     lambda <- (H+90) / area   # accounting for introduction of 90 head of cattle during this time
#   }
#   if(tt >= 1605){
#     lambda <- H / area
#   }
#   pnf.S.lambda <- exp(-H*(pf)*(S)/area)
#   mu.S <- -log(1 - pnf.S.lambda) / S   # mortality rate due to starvation
#   mu.a <- mu.S + mu.b                  # female adult mortality
#   # ODEs
#   deriv <- rep(NA,3)
#   if(tt <= 285){
#     deriv[1] <- H*0 # no change in host population before time 285
#   }
#   if(tt > 285 & tt < 345) { hk = hk_1 }
#   if(tt >=  345) {hk = hk_2  }
#   if(tt > 285){
#     deriv[1] <-  H*(h.r- (h.r+hk)) # change in host population after time 285
#   }
#   deriv[2] <-  larv*A*0.5 - mu.p*P - mud.p*P*P - pup*P
#   deriv[3] <-  P*pup - mu.a*A # change in adults over time
#   return(list(deriv))
# })
# #**************SIMULATE***************************************
# simPop <- function(init=initial, tseq = tseqDay, modFunction=tsetse_mod, parms = tsetse_params()) {
#   simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
#   return(simDat)
# }
# 
# 
# 
# #*************Sensitivitiy analyses******************************
# #*********************Other parameters***************************************
# larv <- c(1/8,1/10,1/12)
# pup <- c(1/30,1/45,1/50)
# mu.p <- c(0.0025,0.006,0.01)
# S <- c(4,6,8)
# mud.p <- c(0,0.00001,0.0001,0.001,0.01)
# pupae.zero <- c(50,250,500)
# 
# mat <- expand.grid(larv,pup,mu.p,S,mud.p,pupae.zero)
# names(mat) <- c("larv","pup","mu.p","S","mud.p","pupae.zero")
# 
# 
# #**************Initial parameter(s)***********************************
# init.pars <- c(log_pf=log(0.1)
#                ,log_mu.b=log(0.01))
# #**************Optimise*******************************************
# optimise.func <- function(vlarv,vpup,vmu.p,vS,vmud.p,vpupae.zero){
# trace <- 3
# optim.vals <- optim(par = init.pars
#                     , objFXN
#                     , fixed.params = tsetse_params(larv=vlarv,pup=vpup,mu.p=vmu.p,S=vS,mud.p=vmud.p,pupae.zero=vpupae.zero)
#                     , dat = tsetse
#                     , control = list(trace = trace, maxit = 150)
#                     , method = "SANN")
# exp(optim.vals$par) # 
# optim.vals <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = tsetse_params(larv=vlarv,pup=vpup,mu.p=vmu.p,S=vS,mud.p=vmud.p,pupae.zero=vpupae.zero)
#                     , dat = tsetse
#                     , control = list(trace = trace, maxit = 150)
#                     , method = "SANN")
# exp(optim.vals$par) # 
# # want to then feed output from SANN to Nelder-Mead
# optim.vals <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = tsetse_params(larv=vlarv,pup=vpup,mu.p=vmu.p,S=vS,mud.p=vmud.p,pupae.zero=vpupae.zero)
#                     , dat = tsetse
#                     , control = list(trace = trace, maxit = 1000, reltol = 10^-7)
#                     , method = "Nelder-Mead" # instead of Nelder-Mead L-BFSG-B allows box constraints could be used below
#                     , hessian = T)
# optim.vals # convergence 0 means algorithm converged
# MLEfits <- optim.vals$par
# exp(MLEfits)
# return(optim.vals)
# }
# 
# 
# fits <- mapply(optimise.func,vlarv=mat[,1],vpup=mat[,2],vmu.p=mat[,3],vS=mat[,4],vmud.p=mat[,5],vpupae.zero=mat[,6])
# 
# loglik <- as.numeric(fits[2,])
# pars <- fits[1,]
# 
# 
# 
# p1 <- as.numeric(lapply(pars, "[[", 1))
# p1 <- exp(p1)
# p2 <- as.numeric(lapply(pars,"[[",2))
# p2 <- exp(p2)
# 
# summ.fits <- cbind.data.frame(loglik,p1,p2,mat)
# sub.summ.fits <- summ.fits[summ.fits$loglik < 395,]
# write.csv(summ.fits,"summ.fits.csv")
