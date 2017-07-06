require(deSolve)
require(bbmle)
source("r_data.R")
hostk2fit1 <- 0.1330609  
hostk2fit2 <- 6.568825e-02  
#************************DATA**************************************
tsetse <- dat[,c(1:4,11)] 
names(tsetse) <- c("year","month","tsetse_over_time","lusulu","time")

#********************Parameters and conditions*************************************
tseqDay <- tsetse$time                     ## times to solve at
tsetse_params <- function(  area = 541     ## total area 
                            , larv = 1/11  ## rate larvae are produced
                            , pup = 1/45   ## rate pupae emerge as adults
                            , mu.b = 0.01 ## female adult background death rate 
                            , mu.p = 0.006 ## pupal initial mortality rate
                            , mud.p = 0.0001
                            , pf = 0.17 ## average probability of  feeding on any given day during feeding cycle given one host present in 1km2
                            , S = 6 ## number of days a fly can go without feeding before it starve
                            , h.r = 0.007 / 30 # host growth rate
                            , hk_2 = as.numeric(hostk2fit2 / 30) #rate is fit to months and need to approx. for days
                            , hk_1 = as.numeric(hostk2fit1 / 30)
                            , hosts.zero = 2100
                            , adults.zero = 250#250
                            , pupae.zero = 250#250
)
return(as.list(environment()))
initial <- c(H=tsetse_params()$hosts.zero, P=tsetse_params()$pupae.zero
             , A=tsetse_params()$adults.zero) # initial conditions
#********************Population dynamics model***********************************
tsetse_mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  pnf.S.lambda <- exp(-2100*(pf)*(S)/area)
  mu.S <- -log(1 - pnf.S.lambda) / S   # mortality rate due to starvation
  mu.a <- mu.S + mu.b                  # adult mortality
  
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
simPop <- function(init=initial, tseq = tseqDay, modFunction=tsetse_mod
                   , parms = tsetse_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}

test1 <- simPop(#tseq=seq(1,5000,30))
)
plot(test1[,1],test1[,4])
plot(test1[,1],test1[,3])
test1[69,]


