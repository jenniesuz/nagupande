library("parallel")
processors <- detectCores()

prepro <- function(dummy, gui
                   ,nl.path="C:/Program Files/NetLogo 6.0.1"   # make sure path correct
                   ,model.path="C:/Users/Jennifer.Lord/Documents/Github/tsetse2016/Manuscript/figures_r_code/nagupande.nlogo"){
  NLStart(nl.path,nl.obj="nagupande",nl.jarname = 'netlogo-6.0.1.jar')
  NLLoadModel(model.path,nl.obj="nagupande")
}

postpro <- function(x){
  NLQuit("nagupande")
}

# a function to handle a simulation
# gets a set of parameters
# returns results of evaluation criteria
# NOTE: It runs repeated simulations for stochastic models.
#       To control stochasticity it runs replicated simulations for current parameter combination
#       and calculates the mean simulation output.
#       If your model is deterministic, just set no.repeated.sim to 1.

#  no.repeated.sim - number of times to repeat simulation with given set of params
#  trace progress if = T prints progress
sim.nagupande <- function(param.set, parameter.names, no.repeated.sim, nl.obj, 
                     trace.progress=F, iter.length, function.name,
                     reporter="cnt-flies") {
  
  if (length(param.set) != length(parameter.names)) # checking parameter set and names are same length
  { stop("Wrong length of param.set!") }
  if (no.repeated.sim <= 0) # simulations must be repeated
  { stop("Number of repetitions must be > 0!") }
  if (length(parameter.names) <= 0) # must give the function some parameter names
  { stop("Length of parameter.names must be > 0!") }
  
  eval.values <- NULL # an empty list to save the simulation results
  
  for (i in 1:no.repeated.sim) # repeated simulations (to control stochasticity)
  {
    NLCommand("random-seed",runif(1,-2147483648,2147483647), nl.obj=nl.obj) # create random seed
    
    # set NetLogo parameters to current parameter values
    lapply(seq(1:length(parameter.names)), function(x) {NLCommand("set ",parameter.names[x], param.set[x], nl.obj=nl.obj)})
    NLCommand("setup", nl.obj=nl.obj) 
    NLCommand("go",nl.obj=nl.obj) # run simulation
    
    output <- NLReport(reporter, nl.obj=nl.obj)
    
    # append to former results
    eval.values <- rbind.data.frame(eval.values,output)
  }
  
  # print the progress if requested
  if (trace.progress == TRUE)
  {
    already.processed <- get("already.processed",env=globalenv()) + 1
    assign("already.processed", already.processed, env=globalenv())
    print(paste("processed (",function.name,"): ", already.processed / iter.length * 100, "%", sep = ""))
  }
  
  # return the mean of the repeated simulation results
  if (no.repeated.sim > 1)
  {
    #return(list(param.set,colMeans(eval.values)))
    return(colMeans(eval.values))
  }
  else {
    return(eval.values)
  }
  
}