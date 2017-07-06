source("r_fig_5_simulate_ibm.R")
library("rJava")
library("RNetLogo")

p.names <- c("model-host-decline","closed-system","distance-move","prob-feed","days-to-starve"
             ,"dd-coef"
           )

reps <- 1

param.values.2 <- expand.grid("model-host-decline" = "true"
                              ,"closed-system" = c("true","false")
                              #,"distance-move" = c(0.25,0.5,0.75,1,1.25,1.5)
                              ,"distance-move" = c(0.5)
                              ,"prob-feed" = 0.135
                              ,"days-to-starve" = 6
                              ,"dd-coef" = 0.00001

)

gui <- T
nl.path <- "C:/Program Files/NetLogo 6.0.1/app"
# put in model path here:
#model.path <- ".../nagupande.nlogo"

cl <- makeCluster(processors)
clusterEvalQ(cl, library(RNetLogo))
invisible(parLapply(cl, 1:processors, prepro, gui=gui,
                    nl.path=nl.path, model.path=model.path))


sim.results <- parApply(cl,param.values.2,1,sim.nagupande,parameter.names=p.names,
                                 no.repeated.sim=reps,nl.obj="nagupande",trace.progress=F)

invisible(parLapply(cl, 1:processors, postpro))
stopCluster(cl)

write.csv(as.data.frame(sim.results),"ibm.resultsa.csv")

