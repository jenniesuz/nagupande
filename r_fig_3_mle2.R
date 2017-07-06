require(deSolve)
require(bbmle)
source("r_data.R")
hostk2fit1 <- 0.1330609  # from r_fig_2_3.R
hostk2fit2 <- 6.568825e-02  # from r_fig_2_3.R
#************************DATA**************************************
tsetse <- dat[,c(1:4,11)] 
names(tsetse) <- c("year","month","tsetse_over_time","lusulu","time")
#***********Likelihood function********************
#******fitting both feeding-related and background mortality
nll.pois <- function(log_pf,log_mu.b,parms=tsetse_params(),dat=tsetse){ 
  parms$pf <- exp(log_pf)
  parms$mu.b <- exp(log_mu.b)
  if(parms$pf >= 1 | parms$pf <= 0.00001 | parms$mu.b > 0.1 ) { 
    ll <- -100000
  }
  if(parms$pf < 1 & parms$pf > 0.00001 & parms$mu.b <= 0.1 ) {
  initial <- c(H=as.numeric(parms$hosts.zero), P=as.numeric(parms$pupae.zero), A=as.numeric(parms$adults.zero)) # initial conditions
  simulate <- simPop(parms=parms) # model output 
  time <- dat$time
  num.tsetse <- round(dat$tsetse_over_time,0) # data 
  est <- round(simulate[simulate$time %in% time,4],0) # model output at timepoints have data for
  ll <- sum(dpois(x=num.tsetse+1,lambda=est+1,log=T)) # adding one so don't get zero flies.
  }
  if (ll == "NaN") {
    ll <- -100000
  }
  if (ll == Inf ){
    ll <- -100000
  }
  return(-ll)
}

# fitting just background mortality
nll.pois.mub <- function(log_mu.b,parms=tsetse_params(),dat=tsetse){ 
  parms$mu.b <- exp(log_mu.b)
  if(parms$mu.b > 0.1 ) {  
    ll <- -100000
  }
  if( parms$mu.b <= 0.1 ) {
    initial <- c(H=as.numeric(parms$hosts.zero), P=as.numeric(parms$pupae.zero), A=as.numeric(parms$adults.zero)) # initial conditions
    simulate <- simPop(parms=parms) # model output 
    time <- dat$time
    num.tsetse <- round(dat$tsetse_over_time,0) # data 
    est <- round(simulate[simulate$time %in% time,4],0) # model output at timepoints have data for
    ll <- sum(dpois(x=num.tsetse+1,lambda=est+1,log=T)) # adding one so don't get zero flies.
  }
  if (ll == "NaN") { 
    ll <- -100000
  }
  if (ll == Inf ){
    ll <- -100000
  }
  return(-ll)
}
#********************Parameters and conditions*************************************
tseqDay <- tsetse$time                        ## times to solve at
tsetse_params <- function(  area = 541        ## total area 
                            , larv = 1/11     ## rate larvae are produced
                            , pup = 1/45      ## rate pupae emerge as adults
                            , mu.b = 0.015    ## female adult background death rate 
                            , mu.p = 0.006    ## pupal initial mortality rate
                            , mud.p = 0.000   ## pupal density-dependent mortality coefficient
                            , pf = 1          ## average probability of not feeding on any given day during feeding cycle given one host present in 1km2
                            , S = 6           ## number of days a fly can go without feeding before it starve
                            , h.r = 0.007 / 30 # host growth rate
                            , hk_2 = as.numeric(hostk2fit2 / 30) #rate is fit to months and need to approx. for days
                            , hk_1 = as.numeric(hostk2fit1 / 30)
                            , hosts.zero = 2100
                            , adults.zero = 250
                            , pupae.zero = 250
)
return(as.list(environment()))
initial <- c(H=tsetse_params()$hosts.zero, P=tsetse_params()$pupae.zero, A=tsetse_params()$adults.zero) # initial conditions
#********************Population dynamics model***********************************
tsetse_mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  pnf.S.lambda <- exp(-H*(pf)*(S)/area)
  mu.S <- -log(1 - pnf.S.lambda) / S   # mortality rate due to starvation
  mu.a <- mu.S + mu.b                  # female adult mortality
  # ODEs
  deriv <- rep(NA,3)
  if(tt <= 285){
    deriv[1] <- H*0 # no change in host population before time 285
  }
  if(tt > 285 & tt < 345) { hk = hk_1 }
  if(tt >=  345) {hk = hk_2  }
  if(tt > 285){
    deriv[1] <-  H*(h.r- (h.r+hk)) # change in host population after time 285
  }
  deriv[2] <-  larv*A*0.5 - mu.p*P - mud.p*P*P - pup*P
  deriv[3] <-  P*pup - mu.a*A # change in adults over time
  return(list(deriv))
})
#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = tseqDay, modFunction=tsetse_mod, parms = tsetse_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#**************Optimise*******************************************
fitmod <- mle2(function(par1,par2){nll.pois(par1,par2)},
               start=list(par1=log(0.124),par2=log(0.018))
               ,data=list(tsetse)) # run with with method = "SANN" first

fitmodmub <- mle2(function(par1){nll.pois.mub(par1)},
                  start=list(par1=-3.299187 )
                  ,data=list(tsetse)) # run with with method = "SANN" first


fit <- simPop(init=initial,tseq=0:2055,parms=tsetse_params(
                                  pf = exp(coef(fitmod)[1])
                                  ,mu.b=exp(coef(fitmod)[2])  
))

fit.mub <- simPop(init=initial,tseq=0:2055,parms=tsetse_params(
 mu.b=exp(coef(fitmodmub)[1])                                 
))
# 
# tiff("fig_4_odefit.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
# par(cex=0.5,mar=c(4,4,1,1))
# par(new = T)
# #plot(tsetse$time, tsetse$tsetse_over_time,bty="n",pch=20,log="y",xlim=c(0,2055),ylim=c(1,1000),
# #     cex.lab="1.2",col="#525252",xlab="Day of experiment",ylab="Mean monthly numbers of tsetse caught")
# par(new = T)
# plot(0:2055, fit[,4], type="l",col="#000000", axes=F, xlab=NA, ylab=NA,xlim=c(0,2055),ylim=c(1,2000),log="y")
# par(new = T)
# plot(0:2055, fit.mub[,4], type="l"
#      ,col="#000000", axes=F, xlab=NA
#      , ylab=NA,xlim=c(0,2055),ylim=c(1,1000),log="y",lty=2)
# 
# legend("bottomleft",legend=c("Observed", 
#                              "Model fit A"
#                              ,"Model fit B"
#                              ),bty="n",
#        lty=c(NA,1,2),pch=c(20,NA,NA),lwd=1,col=c("#525252","#000000","#000000"))
# dev.off()
# 
#***************************confidence intervals*****************************
prof <- profile(fitmod)
ci <- confint(prof)

ci
sl <- slice(fitmod,dim=2)
#*****************************
profmub <- profile(fitmodmub)
ci <- confint(profmub)

sl <- slice(fitmodmub,dim=1)

