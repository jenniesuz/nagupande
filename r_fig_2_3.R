require(deSolve)
require(bbmle)
# Estimating host numbers in Nagupande during game elimination experiment
source("r_data.R") # containing host and tsetse data for Lusulu and Nagupande
#*************************Mammal data***************************************
dat$tot_kill <- rowSums(dat[,5:10]) # add all species to give totals shot per month

tot.kill <- dat$tot_kill[!dat$tot_kill %in% NA] # all months
#tot.kill.1 <- tot.kill[-c(39:61)]
#tot.kill.2 <- tot.kill[-c(1:2,39:61)] # without first two months
tot.kill.1 <- tot.kill[]
tot.kill.2 <- tot.kill[-c(1:2)] # without first two months
#**********************Fitting functions*********************************
#**************************Model parameters**************************
host_params <- function( k_1 = 2000 # = -KN(0) / k_2 estimated (K - carrying capacity, N0 - number of animals at start)
                        ,k_2 = 0.07 # = r - k        estimated (r - growth rate)
                        ,number_months = c(1:length(tot.kill.1))
)
  return(as.list(environment()))
#*********************Model of numbers killed per month***********
month_kill <- function(parms) with(c(parms), {
  N_k <- numeric(length(number_months)) # empty vector to put results in
  for (i in 1:length(number_months)) {  # loop through each month
    N_k[i] <- k_1 * ( exp(-k_2*(number_months[i] - 1)) - exp(-k_2 * number_months[i])) # equation to estimate the number killed each month 
  }
  return(N_k)
})
#***********Likelihood function********************
nll.pois <- function(log_k1,log_k2,parms=host_params(),data=tot.kill.2){ 
  parms$k_1 <- exp(log_k1)
  parms$k_2 <- exp(log_k2)
  ll <- sum(dpois(x=data,lambda=month_kill(parms=parms),log=T))
  return(-ll)
}
#**************Optimise with first two months*******************************************
fit.mod1 <- mle2(function(par1,par2){nll.pois(par1,par2,data=tot.kill.1)},
               start=list(par1=log(1000),par2=log(0.007)),data=list(tot.kill.1))
#***************** fit***********************************
fit.1 <- month_kill(host_params(k_1=exp(coef(fit.mod1)[1]),k_2=exp(coef(fit.mod1)[2])))
#****************confidence intervals**********************
prof1 <- profile(fit.mod1)
ci1 <- confint(prof1)
ci1 <- exp(ci1)

#**************Optimise without first two months*******************************************
fit.mod2 <- mle2(function(par1,par2){nll.pois(par1,par2,parms=host_params(number_months = c(1:length(tot.kill.2))))},
                 start=list(par1=log(1875),par2=log(0.07)),data=list(tot.kill.2))
fit.mod2
#***************** fit***********************************
fit.2 <- month_kill(host_params(k_1=exp(coef(fit.mod2)[1])
                                ,k_2=exp(coef(fit.mod2)[2])
                                ,number_months=c(1:length(tot.kill.2))))
#****************confidence intervals**********************
prof2 <- profile(fit.mod2)
ci2 <- confint(prof2)
ci2 <- exp(ci2)

#*********************PLOT**********************************
# plot of numbers of mammals shot by month of the study
tiff("fig_2_hosts_shot.tiff", height = 4, width = 4, units = 'in', compression="lzw", res=400)
par(cex=0.7,mar=c(4,4,1,1))
plot(dat$time[10:length(dat$time)], 
     dat$tot_kill[10:length(dat$time)],
     pch=20,col="#969696",bty="n", 
     xlab="Day of experiment",ylab="Numbers of hosts shot"
)
lines(dat$time[10:length(dat$time)],fit.1,col="#000000",lwd=2,lty=2)
lines(dat$time[12:length(dat$time)],fit.2,col="#000000",lwd=2,lty=1)
legend("topright", bty="n",legend=c("Observed","Model fit (including first two months)","Model fit (excluding first two months)"),
       cex=0.8,col=c("#969696","#000000","#000000"),pch=c(pch=19,NA,NA),lty=c(NA,2,1),lwd=c(NA,2,2))
dev.off()
#**********Estimating number of hosts at month 12 of experiment (day 345, December 1962)************
# parameters
growth.rate <- seq(0,0.015,0.001)  # range of potential mammal growth rates (months)
num.start <- function(k_1,k_2,r){
  N0 <- -k_1 * -k_2 / (r - -k_2)
  return(N0)
}
number.hosts.1 <- num.start(exp(coef(fit.mod1)[1]),exp(coef(fit.mod1)[2]),growth.rate) # including first two months depending on growth rate
number.hosts.2 <- num.start(exp(coef(fit.mod2)[1]),exp(coef(fit.mod2)[2]),growth.rate) # excluding first two months depending on growth rate
#*******************Estimating number of hosts over time***********************************
host_params <- function( N0 = number.hosts.2[8] # = number of animals at start
                         ,r = growth.rate[8] # = chosen growth rate
                         ,k_2 = exp(coef(fit.mod2)[2]) # = 
                         ,t = 1 # = time
)
  return(as.list(environment()))

N_t <- function(parms) with(c(parms), {
  K <- r + k_2 # kill rate
  num.anim <- N0*exp((r-K)*t)
  return(num.anim)
})
est.mammals.r0 <- numeric(length(dat$time[12:length(dat$time)])) # estimating only numbers of mammals from month 12 onwards
for (i in 1:length(est.mammals.r0)){
  est.mammals.r0[i] <- N_t(host_params(t=i,N0=number.hosts.2[1],r=growth.rate[1]))
}
est.mammals.r0.015 <- numeric(length(dat$time[12:length(dat$time)]))
for (i in 1:length(est.mammals.r0)){
  est.mammals.r0.015[i] <- N_t(host_params(t=i,N0=number.hosts.2[length(number.hosts.2)],r=growth.rate[length(growth.rate)]))
}
est.mammals.r0.007 <- numeric(length(dat$time[12:length(dat$time)]))
for (i in 1:length(est.mammals.r0.007)){
  est.mammals.r0.007[i] <- N_t(host_params(t=i,N0=number.hosts.2[8],r=growth.rate[8]))
}
#**************************estimating numbers of hosts in first two months***********************
# if at start of month 12 (before lower kill rate) there were an estimated 1763 hosts assuming growth rate of 0.007 
dat$tot_kill[10:11] # numbers killed in first two months were 300 and 361
n11 <- est.mammals.r0.007[1]/(1+0.007) + 112 
# estimating numbers of hosts at start of month 11 accounting for chosen growth rate and numbers killed
n10 <- n11/(1+0.007) + 361 # and for month 10
n9 <- n10/(1+0.007) + 300

time.points <- c(255,285,315,345)
hosts.3mths <- c(n9,n10,n11,est.mammals.r0.007[1])
plot(time.points,hosts.3mths,bty="n",pch=19) 
#***********Likelihood function********************
nll.pois <- function(log_k2,parms=host_params(),dat=hosts.3mths[2:4]){ 
  parms$k_2 <- exp(log_k2)
  ll <- sum(dpois(x=round(dat,0),lambda=round(N_t(parms=parms),0),log=T))
  return(-ll)
}
#**************Optimise without first two months*******************************************
fit.mod3 <- mle2(function(par1){nll.pois(par1,parms=host_params(t=1:3,N0=n9))},
                 start=list(par1=log(0.08)),data=list(hosts.3mths[2:4]))
fit.mod3
#****************confidence intervals**********************
ci3 <- confint(fit.mod3)
ci3 <- exp(ci3)
#*********************plot figures****************************************
tiff("S2_fig.tiff", height = 4, width = 4, units = 'in', compression="lzw", res=400)
par(cex=0.7,mar=c(4,4,1,1))
plot(dat$time[12:length(dat$time)],est.mammals.r0,type="l",bty="n",ylim=c(0,2000),lwd=2,
     xlab="Day of experiment",ylab="Estimated number of hosts",col="#000000")
par(new = T)
plot(dat$time[12:length(dat$time)], est.mammals.r0.015, type="l",
     col="#000000",lty=2, axes=F, xlab=NA, ylab=NA,lwd=2,ylim=c(0,2000))
legend("topright",
       legend=c("Assuming monthly growth rate of 0","Assuming monthly growth rate of 0.015"),
       lwd=2, lty=c(1,2), col="#000000",bty="n",cex=0.8)
dev.off()

#***************fitmodel**********
source("r_fig_3_mle2.R")
#*******************************

#*******************PLOT***********************************************
tiff("fig_3_ode_fit.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
par(cex=0.6,mar=c(4,4,1,4))

plot(dat$time, dat$nagupande,pch=20,bty="n", cex.main=0.8,log="y",
     xlab="Day of experiment",ylab="Mean monthly number of tsetse caught",ylim=c(1,2000),xlim=c(0,2055))
points(dat$time,dat$lusulu,pch=1)
abline(v=285,lty=2,col="grey")

par(new = T)
plot(0:2055, fit[,4], type="l",col="#000000", axes=F, xlab=NA, ylab=NA,xlim=c(0,2055),ylim=c(1,2000),log="y")
par(new = T)
plot(0:2055, fit.mub[,4], type="l"
     ,col="#000000", axes=F, xlab=NA
     , ylab=NA,xlim=c(0,2055),ylim=c(1,2000),log="y",lty=2)


legend(x=10,y=5,
       legend=c("Total numbers of hosts at Nagupande"
                ,"Start of host elimination"
                , "Tsetse at control site"
                ,"Tsetse at Nagupande"
                ,"Model fit A"
                ,"Model fit B"
       ),
       lwd=c(1,1,NA,NA,1,1), pch=c(NA,NA,1,20,NA,NA),
       lty=c(1,2,NA,NA,1,2),cex=0.6,bg="white",col=c("grey","grey","black","black","black","black")
       ,box.lty=0)

par(new = T)
plot(dat$time, c(rep(NA,9),n10,n11,est.mammals.r0.007), 
     type="l", axes=F, xlab=NA, ylab=NA,lwd=1,ylim=c(20,2500),xlim=c(0,2055),col="grey")
axis(side = 4)
mtext(side = 4, line = 3, 'Estimated number of hosts',cex=0.6)
dev.off()



