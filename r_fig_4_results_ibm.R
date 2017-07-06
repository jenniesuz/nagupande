library("plyr")
source("r_data.R")
tsetse <- dat[,c(1:4,11)] 
names(tsetse) <- c("year","month","tsetse_over_time","lusulu","time")



param.values.2 <- expand.grid("model-host-decline" = "true"
                              ,"closed-system" = c("true","false")
                              #,"distance-move" = c(0.25,0.5,0.75,1,1.25,1.5)
                              ,"distance-move" = c(0.5)
                              ,"prob-feed" = 0.135
                              ,"days-to-starve" = 6
                              ,"dd-coef" = 0.00001
                              
)

sim.results.1 <- read.csv("ibm.results1.csv",header=T)
sim.results.2 <- read.csv("ibm.results2.csv",header=T)
sim.results.3 <- read.csv("ibm.results3.csv",header=T)
sim.results.4 <- read.csv("ibm.results4.csv",header=T)
sim.results.5 <- read.csv("ibm.results5.csv",header=T)
sim.results.6 <- read.csv("ibm.results6.csv",header=T)
sim.results.7 <- read.csv("ibm.results7.csv",header=T)
sim.results.8 <- read.csv("ibm.results8.csv",header=T)
sim.results.9 <- read.csv("ibm.results9.csv",header=T)
sim.results.10 <- read.csv("ibm.results10.csv",header=T)

summarise.func <- function(data=sim.results.1 ,ncolumns=8){
sim.results.1a <- matrix(data[-1],ncol=ncolumns)
sim.results.1a <- as.data.frame(sim.results.1a)
times <- dat$time # bin counts by month time
months <- c(1:69)
days <- c(1:2055)
month <- numeric(2055)
month[1:30] <- months[1]

for (i in 1:length(months)){
  days[i]
  month[days <= times[i+1] & days > times[i]] <- months[i + 1]
}

sim.results.1a$month <- as.factor(month)
for (i in 1:ncolumns){
sim.results.1a[,i] <- as.numeric(sim.results.1a[,i])
names(sim.results.1a)[i] <- paste0("R",i)
}
sim.1.means <- ddply(sim.results.1a,.(month)
                     ,summarise,mean=mean(R1,na.rm=T))


observed <- tsetse$tsetse_over_time
resids <- numeric(ncolumns)
sim.1.means <- matrix(ncol=ncolumns,nrow=69)

for (i in 1:ncolumns){
  temp <- cbind.data.frame(sim.results.1a[,i],sim.results.1a$month)
  names(temp) <- c("count","month")
  temp2 <- ddply(temp,.(month),summarise,mean=mean(count,na.rm=T))
  #relcount <- temp[,2]/max(temp[,2])*100
  sim.1.means[,i] <- temp2[,2]
  resids[i] <- sum((log(observed+1) - log(temp2[,2]+1))^2)
}

resids.params <- cbind.data.frame(resids,param.values.2[2:6])

names(resids.params) <- c("resids","system","distance","probfeed",
                          "daystostarve","maxpupae")

open <- resids.params[resids.params$system %in% "false",]
closed <- resids.params[resids.params$system %in% "true",]
  
closed.sims <- sim.results.1a[,resids.params$system %in% "true"]
open.sims <- sim.results.1a[,resids.params$system %in% "false"]
return(list(open,closed,closed.sims,open.sims))
}

summary.1 <- summarise.func(sim.results.1,ncolumns=12)
summary.2 <- summarise.func(sim.results.2,ncolumns=12)
summary.3 <- summarise.func(sim.results.3,ncolumns=12)
summary.4 <- summarise.func(sim.results.4,ncolumns=12)
summary.5 <- summarise.func(sim.results.5,ncolumns=12)
summary.6 <- summarise.func(sim.results.6,ncolumns=12)
summary.7 <- summarise.func(sim.results.7,ncolumns=12)
summary.8 <- summarise.func(sim.results.8,ncolumns=12)
summary.9 <- summarise.func(sim.results.9,ncolumns=12)
summary.10 <- summarise.func(sim.results.10,ncolumns=12)

mean.resids.open <- rowMeans(data.frame(summary.1[[1]][1]
                                        ,summary.2[[1]][1]
                                        ,summary.3[[1]][1]
                                        ,summary.4[[1]][1]
                                        ,summary.5[[1]][1]
                                        ,summary.6[[1]][1]
                                        ,summary.7[[1]][1]
                                        ,summary.8[[1]][1]
                                        ,summary.9[[1]][1]
                                        ,summary.10[[1]][1]
))


open.1 <- summary.1[[4]][2]
open.2 <- summary.2[[4]][2]
open.3 <- summary.3[[4]][2]
open.4 <- summary.4[[4]][2]
open.5 <- summary.5[[4]][2]
open.6 <- summary.6[[4]][2]
open.7 <- summary.7[[4]][2]
open.8 <- summary.8[[4]][2]
open.9 <- summary.9[[4]][2]
open.10 <- summary.10[[4]][2]

all.open <- cbind.data.frame(open.1,open.2,open.3,open.4,open.5,open.6,open.7)
openmeans <- rowMeans(all.open)

mean.resids.closed <- rowMeans(data.frame(summary.1[[2]][1]
                                        ,summary.2[[2]][1]
                                        ,summary.3[[2]][1]
                                        ,summary.4[[2]][1]
                                        ,summary.5[[2]][1]
                                        ,summary.6[[2]][1]
                                        ,summary.7[[2]][1]
))


closed.1 <- summary.1[[3]][2]
closed.2 <- summary.2[[3]][2]
closed.3 <- summary.3[[3]][2]
closed.4 <- summary.4[[3]][2]
closed.5 <- summary.5[[3]][2]
closed.6 <- summary.6[[3]][2]
closed.7 <- summary.7[[3]][2]

all.closed <- cbind.data.frame(closed.1,closed.2,closed.3,closed.4,closed.5,closed.6,closed.7)
closedmeans <- rowMeans(all.closed)

#test <- summary.1[[4]]/max(summary.1[[4]])*100+1


tiff("fig_5_IBM.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
par(cex=0.5,mar=c(4,4,1,1))
par(new = T)
plot(1:2055,
     openmeans+1
     ,bty="n",log="y",type="l",axes=F
     ,ylim=c(1,1000),col="grey57",lwd=1, xlab=NA, ylab=NA,
     lty=3)
par(new = T)

plot(1:2055,closedmeans
     ,bty="n",log="y",type="l",ylim=c(1,1000),col="darkgrey",lwd=1
     ,axes=F, xlab=NA, ylab=NA)
par(new = T)
plot(tsetse$time
     ,tsetse$tsetse_over_time
     ,bty="n",log="y"
     ,pch=20,ylim=c(1,1000)
     ,ylab="Mean monthly numbers of tsetse caught"
     ,xlab="Day of experiment")
#axis(2, at = yticks+1, labels = yticks, col.axis="black", las=2)

legend(y=5, x=100,
       legend=c("Observed","Model A","Model B"),bty="n",
       lwd=c(NA,1,1), pch=c(20,NA,NA),
       lty=c(NA,3,1),cex=0.8,col=c("black","grey57","darkgrey"))
dev.off()
