### bmeta - A package for Bayesian Meta-Analysis and Meta-Regression in R
###
bmeta<- function(data,outcome=c("bin","ctns","count"),model=c("std.norm","std.dt","reg.norm","reg.dt",
                                                              "std.ta","std.mv","reg.ta","reg.mv",
                                                              "std","std.unif","std.hc","reg","reg.unif","reg.hc"),
                 type=c("fix","ran"),n.iter=10000,n.burnin=5000,n.samples=1000,n.chains=2,model.file="model.txt")  UseMethod("bmeta")


### Call to the actual function ###
bmeta.default<- function(data,outcome=c("bin","ctns","count"),model=c("std.norm","std.dt","reg.norm","reg.dt",
                                                                      "std.ta","std.mv","reg.ta","reg.mv","std","std.unif","std.hc","reg","reg.unif","reg.hc"),
                         type=c("fix","ran"),n.iter=10000,n.burnin=5000,n.samples=1000,n.chains=2,model.file="model.txt"){
  
  ### Defines the "quiet" function (which prevents from showing the messages)
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  
  #############################################
  ## 0. Set up the required path and libraries
  requireNamespace("R2jags",quietly=TRUE)
  working.dir<- getwd()
  
  
  
  ############################################
  # 1. Define length of vector and data list 
  
  if (model=="std.mv"| model=="reg.mv"){
    S<-length(data$y)                    
  } else {
    S<-length(data$y0)                       #sample size for a study
  }
  
  ##########################################
  ## 2. Checks that if there are covariates
  chk<- is.null(data$X) 
  #if (chk==TRUE) {
  #    J=1
  # }
  if (chk==FALSE) {
    J <-dim(data$X)[2]
    
    # then checks if covariates are centered
    if(!is.null(J)) {
      threshold<- 1e-12
      Z0<-data$X
      for (j in 1:J) {
        if (mean(data$X[,j])>threshold) {
          Z0[,j]<-scale(data$X[,j],scale=FALSE)
        }
      }      
    }
  }
  
  # Define the data list
  
  if (outcome=="bin" & (model=="std.norm"|model=="std.dt")){
    dataJags<-list(S=S,y0=data$y0,n0=data$n0,y1=data$y1,n1=data$n1)
  } 
  
  if (outcome=="bin" & (model=="reg.norm"|model=="reg.dt")){
    dataJags<-list(S=S,y0=data$y0,n0=data$n0,y1=data$y1,n1=data$n1,m.beta0=rep(0,J),
                   tau.beta0=.0001*diag(J),J=J,Z0=Z0)
  }
  
  if (outcome=="ctns" & model=="std.ta"){
    dataJags<-list(S=S,y0=data$y0,se0=data$se0,y1=data$y1,se1=data$se1)
  }
  
  if (outcome=="ctns" & model=="std.mv"){
    dataJags<-list(S=S,y=data$y,prec=data$prec)
  }
  
  if (outcome=="ctns" & model=="reg.ta"){
    dataJags<-list(S=S,y0=data$y0,se0=data$se0,y1=data$y1,se1=data$se1,m.beta0=rep(0,J),
                   tau.beta0=.0001*diag(J),J=J,Z0=Z0)
  }
  
  if (outcome=="ctns" & model=="reg.mv"){
    dataJags<-list(S=S,y=data$y,prec=data$prec,m.beta0=rep(0,J),
                   tau.beta0=.0001*diag(J),J=J,Z0=Z0)
  }
  
  if (outcome=="count" & (model=="std"|model=="std.unif"|model=="std.hc")){
    dataJags<-list(S=S,y0=data$y0,pt0=data$pt0,y1=data$y1,pt1=data$pt1)
  }
  
  if (outcome=="count" & (model=="reg"|model=="reg.unif"|model=="reg.hc")){
    dataJags<-list(S=S,y0=data$y0,pt0=data$pt0,y1=data$y1,pt1=data$pt1,J=J,Z0=Z0,m.beta0=rep(0,J),
                   tau.beta0=.0001*diag(J))
  } 
  
  
  #######################################################################################################
  # 3. writes the model code to a file using the function 'writemodel'. This selects the required modules
  # and then assign the name of the file to the variable 'filein'
  
  writeModel(outcome=outcome, model=model, type=type, model.file=model.file, data=dataJags)
  filein<-model.file
  
  
  ####################################
  # 4. Defines the parameters vector
  
  if (outcome=="bin" & type=="fix"){
    params <- c("rho","alpha","delta")
  } 
  
  if (outcome=="bin" & type=="ran"){
    params <- c("rho","alpha","delta","tau.delta","gamma","mu")
  }
  
  if (outcome=="ctns" & type=="fix" & (model=="std.ta" | model=="reg.ta")){
    params <- c("alpha0","alpha1","delta")
  }
  
  if (outcome=="ctns" & type=="ran" & (model=="std.ta" | model=="reg.ta")){
    params <- c("alpha0","alpha1","delta","mu","tau.delta")
  }
  
  if (outcome=="ctns" & type=="fix" & model=="std.mv" ){
    params <- c("mu")
  }
  
  if (outcome=="ctns" & type=="ran" & model=="std.mv"){
    params <- c("mu","delta","tau.delta")
  }  
  
  if (outcome=="ctns" & type=="fix" & model=="reg.mv" ){
    params <- c("alpha")
  }    
  
  if (outcome=="ctns" & type=="ran" & model=="reg.mv"){
    params <- c("alpha","tau.alpha","mu")
  }  
  
  if (outcome=="count" & type=="fix"){
    params <- c("lambda0","lambda1","delta","IRR")
  }
  
  if (outcome=="count" & type=="ran"){
    params <- c("lambda0","lambda1","mu","delta","tau.delta","IRR","gamma")
  }
  
  
  
  
  ##################################################################
  # 5. Defines the initial values for the random nodes in the model
  
  ### Binary outcome###
  
  if (outcome=="bin" & type=="fix" & (model=="std.norm"|model=="std.dt")){  
    list.temp <- list(alpha=rnorm(S,0,1),delta=rnorm(1))
  }
  
  if (outcome=="bin" & type=="ran" & (model=="std.norm"|model=="std.dt")){
    list.temp <- list(alpha=rnorm(S,0,1),delta=rnorm(S,0,1),mu=rnorm(1),sigma.delta=runif(1))  
  }
  
  if (outcome=="bin" & type=="fix" & (model=="reg.norm"|model=="reg.dt")){
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),alpha=rnorm(S,0,1),delta=rnorm(1))  
  }  
  
  if (outcome=="bin" & type=="ran" & (model=="reg.norm"|model=="reg.dt")){
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),alpha=rnorm(S,0,1),delta=rnorm(S,0,1),
                      mu=rnorm(1),sigma.delta=runif(1))  
  }
  
  ### Continuous outcome ###
  
  if (outcome=="ctns" & type=="fix" & model=="std.ta"){  
    list.temp <- list(alpha0=rnorm(S,0,1),delta=rnorm(1))
  }
  
  if (outcome=="ctns" & type=="ran" & model=="std.ta"){  
    list.temp <- list(alpha0=rnorm(S,0,1),mu=rnorm(1),sigma=runif(1))
  }
  
  if (outcome=="ctns" & type=="fix" & model=="std.mv"){  
    list.temp <- list(mu=rnorm(1))
  }
  
  if (outcome=="ctns" & type=="ran" & model=="std.mv"){  
    list.temp <- list(mu=rnorm(1),sigma=runif(1))
  }
  
  if (outcome=="ctns" & type=="fix" & model=="reg.ta"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),gamma0=rnorm(S,0,1),delta=rnorm(1))
  }
  
  if (outcome=="ctns" & type=="ran" & model=="reg.ta"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),gamma0=rnorm(S,0,1),mu=rnorm(1),sigma=runif(1))
  }
  
  if (outcome=="ctns" & type=="fix" & model=="reg.mv"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),alpha=rnorm(1))
  }
  
  if (outcome=="ctns" & type=="ran" & model=="reg.mv"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),mu=rnorm(1),sigma=runif(1))
  }
  
  ### Count outcome ###
  
  if (outcome=="count" & model=="std"){  
    list.temp <- list(lograte0=runif(S,0,1),delta=rnorm(1))
  } 
  
  if (outcome=="count" & model=="std.unif"){  
    list.temp <- list(lograte0=runif(S,0,1),mu=rnorm(1),sigma=runif(1))
  } 
  
  if (outcome=="count" & model=="std.hc"){  
    list.temp <- list(lograte0=runif(S,0,1),mu=rnorm(1),z.sigma=rnorm(1),epsilon.sigma=rgamma(1,0.5,0.5),B.sigma=runif(1,0,0.5))
  } 
  
  if (outcome=="count" & model=="reg"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),logbrate=runif(S,0,1),delta=rnorm(1))
  } 
  
  if (outcome=="count" & model=="reg.unif"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),logbrate=runif(S,0,1),mu=rnorm(1),sigma=runif(1))
  } 
  
  if (outcome=="count" & model=="reg.hc"){  
    list.temp <- list(beta0=rnorm(dataJags$J,0,1),logbrate=runif(S,0,1),mu=rnorm(1),z.sigma=rnorm(1),epsilon.sigma=rgamma(1,0.5,0.5),B.sigma=runif(1,0,0.5))
  } 
  
  inits <- function(){
    list.temp
  }
  
  
  
  
  #######################################################
  # 6. Runs JAGS to produce the posterior distributions
  n.iter <- n.iter                           # number of iterations
  n.burnin <- n.burnin                       # number of burn-in
  if (n.burnin>n.iter) {n.burnin <- n.iter/2}
  n.samples <- n.samples 
  n.thin <- floor((n.iter-n.burnin)/n.samples)     # number of thinning so that 1000 iterations are stored
  if(n.thin < 1) {n.thin <- 1}
  n.chains <- n.chains                       # number of Markov chains
  
  mod<- R2jags::jags(dataJags,inits,params,model.file=filein,n.chains=2,
                     n.iter,n.burnin,n.thin,DIC=TRUE,working.directory=working.dir, progress.bar="text")
  
  
  #######################################################
  # 7. Runs JAGS for the "null" independence model
  ##writeModel(outcome=outcome, model="indep", type=type, model.file=model.file, data=dataJags)
  file0 <- "file0.txt"
  if (outcome=="bin") {
    cat("model {
        for (s in 1:S) {
        y0[s] ~ dbin(pi0[s],n0[s])
        y1[s] ~ dbin(pi1[s],n1[s])
        logit(pi0[s]) <- alpha0[s]
        logit(pi1[s]) <- alpha1[s]
        alpha0[s] ~ dnorm(0,1.45)                       ### assume OR>10 is very unlikely
        alpha1[s] ~ dnorm(0,1.45)
        or[s] <- exp(alpha1[s]-alpha0[s])
        lor[s] <- log(or[s])
        }
  }",file=file0)
    dataJags0<-list(S=S,y0=data$y0,n0=data$n0,y1=data$y1,n1=data$n1)
    param0 <- c("or","lor")
    inits0 <- function(){
      list(alpha0=rnorm(S,0,1),alpha1=rnorm(S,0,1))
    }
  }
  
  
  if (outcome=="ctns" & (model=="std.ta" | model=="reg.ta")) {
    cat("model {
        for (s in 1:S) {
        y0[s] ~ dnorm(mu0[s],prec0[s])
        y1[s] ~ dnorm(mu1[s],prec1[s])
        mu0[s] ~ dnorm(0,0.0001)
        mu1[s] ~ dnorm(0,0.0001)
        prec0[s] <- pow(se0[s],-2)
        prec1[s] <- pow(se1[s],-2)
        diff[s] <- mu1[s]-mu0[s]
        }
  }",file=file0)
    dataJags0<-list(S=S,y0=data$y0,se0=data$se0,y1=data$y1,se1=data$se1)
    param0 <- c("diff")
    inits0 <- function(){
      list(mu0=rnorm(S,0,1),mu1=rnorm(S,0,1))
    }
  }
  
  if (outcome=="ctns" & (model=="std.mv" | model=="reg.mv")) {
    cat("model {
        for (s in 1:S) {
        y[s] ~ dnorm(mu[s],prec[s])
        mu[s] ~ dnorm(0,0.0001)
        }
  }",file=file0)
    dataJags0<-list(S=S,y=data$y,prec=data$prec)
    param0 <- c("mu")
    inits0 <- function(){
      list(mu=rnorm(S,0,1))
    }
  }
  
  if (outcome=="count") {
    cat("model {
        for (s in 1:S) {
        y0[s] ~ dpois(lambda0[s])
        y1[s] ~ dpois(lambda1[s])
        log(lambda0[s]) <- lograte0[s]+log(pt0[s])
        log(lambda1[s]) <- lograte1[s]+log(pt1[s])
        lograte0[s] ~ dunif(-5,5)
        lograte1[s] ~ dunif(-5,5)
        irr[s] <- exp(lograte1[s]-lograte0[s])
        lirr[s] <- lograte1[s]-lograte0[s]
        }
  }",file=file0)
    dataJags0<-list(S=S,y0=data$y0,pt0=data$pt0,y1=data$y1,pt1=data$pt1)
    param0 <- c("irr","lirr")
    inits0 <- function(){
      list(lograte0=runif(S,0,1),lograte1=runif(S,0,1))
    }
  }
  
  # Now also run the "null" model
  quiet(mod0 <- R2jags::jags(dataJags0,inits0,param0,model.file="file0.txt",n.chains=2,
                             n.iter,n.burnin,n.thin,DIC=TRUE,working.directory=working.dir, progress.bar="text"))
  
  unlink(file0,force=TRUE)  # Now deletes the file "file0.txt" from current directory
  
  #########################################
  # 8. Defines the output of the function
  out <- list(mod=mod, params=params, data=dataJags, inits=inits, outcome=outcome, type=type, 
              model=model,mod0=mod0)
  class(out) <- "bmeta"
  out
}



########################################
print.bmeta<-function(x,...){
  print(x$mod,interval=c(0.025,0.975),digit=3)
}


##### POSTERIOR PLOT ##### 
posterior.plot<-function(x,xlim=NULL,xlab="",main="Posterior distribution Plot",scale="log"){
  requireNamespace("R2jags",quietly=TRUE)
  
  param.sims <- x$mod$BUGSoutput$sims.matrix
  if (x$outcome=="bin" & x$type=="fix"){ 
    if(scale=="log") {param <- param.sims[,"delta"]} 
    if(scale=="exp") {param <- param.sims[,"rho"]}
    
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)}
  }
  
  
  if (x$outcome=="bin" & x$type=="ran"){ 
    if(scale=="log") {param <- param.sims[,"mu"]} 
    if(scale=="exp") {param <- param.sims[,"rho"]}
    
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)} 
  }
  
  
  
  if (x$outcome=="ctns" & x$type=="fix" & (x$model=="std.ta"|x$model=="reg.ta")){
    param <- param.sims[,"delta"]	
    if(is.null(xlim)){  
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)
  }
  
  
  if (x$outcome=="ctns" & x$type=="ran" & (x$model=="std.ta"|x$model=="reg.ta")){
    param <- param.sims[,"mu"]
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
      
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)
  }
  
  
  if (x$outcome=="ctns" & x$model=="std.mv"){ 
    param <- param.sims[,"mu"] 
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)  
  }
  
  
  if (x$outcome=="ctns" & x$model=="reg.mv" & x$type=="fix"){  
    param <- param.sims[,"alpha"]
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)  
  }
  
  
  if (x$outcome=="ctns" & x$model=="reg.mv" & x$type=="ran"){  
    param <- param.sims[,"mu"]
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)  
  }
  
  
  if (x$outcome=="count" & x$type=="fix"){
    if(scale=="log") {param <- param.sims[,"delta"]} 
    if(scale=="exp") {param <- param.sims[,"IRR"]}
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)}
  }
  
  
  
  if (x$outcome=="count" & x$type=="ran"){
    if(scale=="log") {param <- param.sims[,"mu"]} 
    if(scale=="exp") {param <- param.sims[,"IRR"]}
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)}
  }
  
  
  
}



##### FOREST PLOT #####
forest.plot <- function(x,title=NULL,xlab=NULL,log=FALSE,study.label=NULL,clip=c(-3,3),
                        lines="black",box="blue",summary="orange",box.symb="box",label.cex=.8,
                        xlab.cex=1,ticks.cex=.8,...) {
  
  # x = a bmeta object including the MCMC simulations for the results
  # title = possibly a string including the title for the forest plot
  # log = indicates whether to report the analysis on the log (TRUE) or the natural (FALSE, default) scale
  # study.label = a vector of study labels. If NULL then forest.plot will create one 
  # clip = graphical parameter to determine the range of values showed for the x-axis
  # lines = selects the colour for the lines of the intervals
  # box = selects the colour for the mean estimate
  # summary = selects the colour for the pooled estimate
  # box.symb = selects the symbol used to plot the mean. Options are "box" (default) or "circle"
  # label.cex = defines the size of the text for the label. Defaults at .8 of normal size
  # xlab.cex = defines the size of the text for the x-label. Defaults at 1 of the normal size
  # ticks.cex = defines the size of the text for the x-axis ticks. Defaults at .8 of the normal size
  # ... = additional arguments, including
  #  - add.null = TRUE/FALSE. If set to true, adds a plot of the null (no-pooling model)
  #  - line.margin = the distance between lines in case multiple graphs are shown on the same plot
  #  - box.size = the size of the summary box
  #  - new.page = TRUE/FALSE. If set to true, then a new graph overwrite the existing one
  #  - zero (x-axis coordinate for zero line. If you provide a vector of length 2 it will print 
  #         a rectangle instead of just a line. Default at 0 or 1 depending on log scale)
  #  - legend = a legend for the multi-graph plot
  
  exArgs <- list(...)  
  
  tab0 <- x$mod0$BUGSoutput$summary 
  tab <- x$mod$BUGSoutput$summary 
  
  if (is.null(study.label)) {
    study.label <- c(paste0("Study",1:x$data$S),"Summary")
  }
  
  is.sum=c(rep(FALSE,length(study.label)-1),TRUE)
  
  # Defines the default values for the original forestplot options
  if(exists("line.margin",where=exArgs)){line.margin=exArgs$line.margin} else {line.margin=0.1}
  if(exists("boxsize",where=exArgs)){boxsize=exArgs$boxsize} else {boxsize=0.4}
  if(exists("new_page",where=exArgs)){new_page=exArgs$new_page} else {new_page=TRUE}
  
  
  # If no xlabel has been specified, selects a suitable one
  if(is.null(xlab)) {
    cond <- c((x$outcome=="bin" & log==FALSE),(x$outcome=="bin" & log==TRUE),
              x$outcome=="ctns",(x$outcome=="count" & log==FALSE),
              (x$outcome=="count" & log==TRUE))
    labs <- c("OR","log(OR)","Mean","IRR","log(IRR)")
    xlab <- labs[which(cond==TRUE)]
  }
  
  
  ### different color/box symbol option for simple plots and plots showing results from two different models 
  # Checks whether the null (no-pooling) model should be also plotted & sets the relevant parameters for either cases
  if(exists("add.null",where=exArgs)) {add.null=exArgs$add.null} else {add.null=FALSE}
  if (add.null==FALSE){
    col <- fpColors(lines=lines,box=box,summary=summary)
    fn.ci_norm <- fpDrawNormalCI
    if(box.symb=="circle") {fn.ci_norm <- fpDrawCircleCI}
  } else {
    if(length(lines)!=2){lines <- c("grey","black")}
    if(length(box)!=2){box <-c("blue","darkred")}
    col <- fpColors(lines=c(lines[1],lines[2]),box=c(box[1],box[2]),summary=summary)
    fn.ci_norm <- c(fpDrawNormalCI, fpDrawCircleCI)     
  }
  
  txt_gp=fpTxtGp(xlab=grid::gpar(cex=xlab.cex),ticks=grid::gpar(cex=ticks.cex),label=grid::gpar(cex=label.cex))
  
  
  ### forest plot for fixed-effects model ###
  if (x$type=="fix") {
    
    if(x$outcome=="bin"){
      ind <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      ind0 <- which((substr(rownames(tab0),1,2)=="or")==TRUE)
      ind1 <- which((substr(rownames(tab),1,3)=="rho")==TRUE)
      ind2 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      
      if(log==FALSE){
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab0[ind0,1],tab[ind1,1]),lower=c(tab0[ind0,3],tab[ind1,3]),
                               upper=c(tab0[ind0,7],tab[ind1,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,xlog=xlog,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
      } else {
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab0[ind,1],tab[ind2,1]),lower=c(tab0[ind,3],tab[ind2,3]),
                               upper=c(tab0[ind,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,clip=clip,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
        
      }
      
    } 
    
    if(x$outcome=="ctns") {
      log==FALSE
      
      if(x$model=="std.ta" | x$model=="reg.ta"){
        ind0 <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
        ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)    
      }
      
      
      if(x$model=="std.mv"){
        ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
        ind1 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
      }
      
      if(x$model=="reg.mv"){
        ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
        ind1 <- which((substr(rownames(tab),1,5)=="alpha")==TRUE)
      }
      
      if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
      
      forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab0[ind0,1],tab[ind1,1]),lower=c(tab0[ind0,3],tab[ind1,3]),
                             upper=c(tab0[ind0,7],tab[ind1,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                             txt_gp=txt_gp,xlab=xlab,col=col,clip=clip,zero=zero,fn.ci_norm=fn.ci_norm)  
      
      
    }
    
    if(x$outcome=="count"){
      ind <- which((substr(rownames(tab0),1,4)=="lirr")==TRUE)
      ind0 <- which((substr(rownames(tab0),1,3)=="irr")==TRUE)
      ind1 <- which((substr(rownames(tab),1,3)=="IRR")==TRUE)
      ind2 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      
      if(log==FALSE){
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab0[ind0,1],tab[ind1,1]),lower=c(tab0[ind0,3],tab[ind1,3]),
                               upper=c(tab0[ind0,7],tab[ind1,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,xlog=xlog,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
      } else {
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab0[ind,1],tab[ind2,1]),lower=c(tab0[ind,3],tab[ind2,3]),
                               upper=c(tab0[ind,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,clip=clip,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
        
      }
      
      
    }
    
    
  }
  
  
  ### forest plot for random-effects models ###   
  if(x$type=="ran"){
    if(x$outcome=="bin" & log==TRUE){
      ind0 <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
    } 
    
    if(x$outcome=="bin" & log==FALSE){
      ind0 <- which((substr(rownames(tab0),1,2)=="or")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="gamma")==TRUE)
      ind2 <- which((substr(rownames(tab),1,3)=="rho")==TRUE)
    }
    
    if(x$outcome=="ctns" & (x$model=="std.ta"|x$model=="reg.ta")){
      log=TRUE
      ind0 <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="ctns" & x$model=="std.mv"){
      log=TRUE
      ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }  
    
    if(x$outcome=="ctns" & x$model=="reg.mv"){
      log==TRUE
      ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="alpha")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="count" & log==TRUE){
      ind0 <- which((substr(rownames(tab0),1,4)=="lirr")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
    } 
    
    if(x$outcome=="count" & log==FALSE){
      ind0 <- which((substr(rownames(tab0),1,3)=="irr")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="gamma")==TRUE)
      ind2 <- which((substr(rownames(tab),1,3)=="IRR")==TRUE)
    }
    
    
    if (add.null==FALSE){
      
      if((x$outcome=="bin" & log==FALSE)|(x$outcome=="count" & log==FALSE)){
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab[ind1,1],tab[ind2,1]),lower=c(tab[ind1,3],tab[ind2,3]),
                               upper=c(tab[ind1,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,xlog=xlog,col=col,zero=zero,fn.ci_norm=fn.ci_norm)
      } else {
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=c(tab[ind1,1],tab[ind2,1]),lower=c(tab[ind1,3],tab[ind2,3]),
                               upper=c(tab[ind1,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,clip=clip,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
        
      }
      
      
    } else {
      
      ### adding results from model (study-specific estimates and summary estimate) together
      comb1<-c(tab[ind1,1],tab[ind2,1]) 
      comb2<-c(tab[ind1,3],tab[ind2,3])
      comb3<-c(tab[ind1,7],tab[ind2,7])
      
      ### results from null model but how to avoid the summary estimate appear twice???   
      comb4<-c(tab0[ind0,1],500) 
      comb5<-c(tab0[ind0,3],500)
      comb6<-c(tab0[ind0,7],500)
      
      mean0<-as.matrix(comb4)
      lower0<-as.matrix(comb5)
      upper0<-as.matrix(comb6)
      
      mean1<-as.matrix(comb1)
      lower1<-as.matrix(comb2)
      upper1<-as.matrix(comb3)
      
      if((x$outcome=="bin" & log==FALSE)|(x$outcome=="count" & log==FALSE)){
        
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        if(exists("legend",where=exArgs)){
          legend=exArgs$legend
        } else {
          legend=c("Estimates from random-effects model","Estimates from no-pooling effects model")
        }
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=cbind(mean1[,1],mean0[,1]),
                               lower=cbind(lower1[,1],lower0[,1]),upper=cbind(upper1[,1],upper0[,1]),
                               is.summary=is.sum,new_page=new_page,boxsize=boxsize,line.margin=line.margin,xlog=xlog,
                               fn.ci_norm=fn.ci_norm,col=col,title=title,legend=legend,txt_gp=txt_gp,
                               xlab=xlab,clip=clip,zero=zero)  
      } else{
        
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        if(exists("legend",where=exArgs)){
          legend=exArgs$legend
        } else {
          legend=c("Estimates from random-effects model","Estimates from no-pooling effects model")
        }
        
        forestplot::forestplot(labeltext=c(study.label,""),mean=cbind(mean1[,1],mean0[,1]),
                               lower=cbind(lower1[,1],lower0[,1]),upper=cbind(upper1[,1],upper0[,1]),
                               is.summary=is.sum,new_page=new_page,boxsize=boxsize,line.margin=line.margin,
                               fn.ci_norm=fn.ci_norm,col=col,title=title,legend=legend,txt_gp=txt_gp,
                               xlab=xlab,clip=clip,zero=zero)  
        
      }
      
      
      
      
    }
    
  }
  
}   





###### Funnel Plot #####
funnel.plot<-function(x,xlab=NULL,ylab=NULL,title=NULL){
  
  tab0 <- x$mod0$BUGSoutput$summary
  tab <- x$mod$BUGSoutput$summary 
  
  if (is.null(xlab)){
    xlab="effect"
  }
  
  if (is.null(ylab)){
    ylab="size"
  }
  
  if (x$type=="fix") {
    
    if(x$outcome=="bin"){
      ind <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      plot(tab0[ind,1],1/tab0[ind,2],xlab=xlab,ylab=ylab,main=title)
      abline(v=0,lty=2)
    } 
    
    if(x$outcome=="ctns") {
      
      if(x$model=="std.ta" | x$model=="reg.ta"){
        ind <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
      }
      
      if(x$model=="std.mv" | x$model=="reg.mv"){
        ind <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)   
      }     
      
      plot(tab0[ind,1],1/tab0[ind,2],xlab=xlab,ylab=ylab,main=title)
      abline(v=0,lty=2)
    }
    
    if(x$outcome=="count"){
      
      ind <- which((substr(rownames(tab0),1,4)=="lirr")==TRUE)      
      plot(tab0[ind,1],1/tab0[ind,2],xlab=xlab,ylab=ylab,main=title)
      abline(v=0,lty=2)
      
    } 
    
  }
  
  if(x$type=="ran"){
    if(x$outcome=="bin"){
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
    } 
    
    
    if(x$outcome=="ctns" & (x$model=="std.ta"|x$model=="std.mv"|x$model=="reg.ta")){
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="ctns" & x$model=="reg.mv"){
      ind1 <- which((substr(rownames(tab),1,5)=="alpha")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="count"){
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
    } 
    
    plot(tab[ind1,1],1/tab[ind1,2],xlab=xlab,ylab=ylab,main=title)
    abline(v=0,lty=2)
    
  } 
  
}





##### DIAG PLOT #####
diag.plot <- function(x,diag="Rhat") {
  tab <- x$mod$BUGSoutput$summary
  
  if(diag=="Rhat") {
    plot(tab[,"Rhat"],xlab="Parameters",ylab="Gelman-Rubin statistic (Rhat)",col="white",
         main="Convergence diagnostics")
    text(1:dim(tab)[1],tab[,"Rhat"],rownames(tab),cex=.6)
    abline(h=1.1,lwd=2,lty=2)
  }
  if (diag=="n.eff") {
    plot(tab[,"n.eff"],xlab="Parameters",ylab="Effective sample size",col="white",
         main="Convergence diagnostics")
    text(1:dim(tab)[1],tab[,"n.eff"],rownames(tab),cex=.6)
    abline(h=x$mod$BUGSoutput$n.sims,lwd=2,lty=2)    
  }
}


##### TRACEPLOT #####
traceplot.bmeta <- function(x,node,title="",lab=""){
  requireNamespace("R2jags",quietly=TRUE)
  ## node is a string with the name of the node to be plotted
  ## x is the name of the bmeta object containing the MCMC simulations
  xlab <- "Iteration"
  ## this way works with R2jags as well as with R2WinBUGS
  ###  cmd <- ifelse(class(model)=="rjags",mdl <- model$BUGSoutput,mdl <- model)
  mdl <- x$mod$BUGSoutput
  col <- colors()
  if (mdl$n.chains==1) {
    plot(mdl$sims.array[,1,node],t="l",col=col[which(col=="blue")],xlab=xlab,
         ylab=lab,main=title)
  }
  if (mdl$n.chains==2) {
    cols <- c("blue","red")
    plot(mdl$sims.array[,1,node],t="l",col=col[which(col==cols[1])],xlab=xlab,
         ylab=lab,main=title,ylim=range(mdl$sims.array[,1:2,node]))
    points(mdl$sims.array[,2,node],t="l",col=col[which(col==cols[2])])
  }
  if (mdl$n.chains>2) {
    cols <- c("blue","red","green","magenta","orange","brown","azure")
    plot(mdl$sims.array[,1,node],t="l",col=col[which(col==cols[1])],xlab=xlab,
         ylab=lab,main=title,
         ylim=range(mdl$sims.array[,1:mdl$n.chains,node]))
    for (i in 2:mdl$n.chains) {
      points(mdl$sims.array[,i,node],t="l",col=col[which(col==cols[i])])
    }
  }
}


##### ACF PLOT #### 
acf.plot <- function(x,node,title="Autocorrelation function") {
  requireNamespace("R2jags",quietly=TRUE)
  ##cmd <- ifelse(class(m)=="rjags",mdl <- m$BUGSoutput,mdl <- m)
  mdl <- x$mod$BUGSoutput
  ind2 <- 1:dim(mdl$summary)[1]
  ind <- which(colnames(mdl$sims.matrix)==node)
  if (var(mdl$sims.matrix[,ind])==0) {return(invisible())}
  acf(mdl$sims.matrix[,ind],ylab="Autocorrelation",main=title)
}




#### WRITE MODEL #####
writeModel<-function(outcome,model,type,model.file,data){
  
  ## Model selection
  ## Selects modules in the model code, according to the distributional assumptions and whether it is a standard meta-analysis or meta-regression
  
  
  
  
  ##############################################
  ## (i) standard meta-analysis for binary data
  
  
  sel.mod.bin.std.norm.fix<-" 
  model {   \n## a. binary standard fixed-effects meta-analysis with normal prior
  for (s in 1:S){
  y0[s]~dbin(pi0[s],n0[s])
  y1[s]~dbin(pi1[s],n1[s])
  logit(pi0[s])<-alpha[s]
  logit(pi1[s])<-alpha[s]+delta
  alpha[s]~dnorm(0,0.0001)
  }
  ###prior###
  delta~dnorm(0,0.0001)
  rho<-exp(delta)
  } "
  
  
  
  
  
  sel.mod.bin.std.dt.fix<-" 
  model {  \n## b. binary standard fixed-effects meta-analysis with t-distribution prior
  for (s in 1:S){
  y0[s]~dbin(pi0[s],n0[s])
  y1[s]~dbin(pi1[s],n1[s])
  logit(pi0[s])<-alpha[s]
  logit(pi1[s])<-alpha[s]+delta
  alpha[s]~dnorm(0,0.0001)
  }
  ###prior###
  delta~dt(0,0.5,v)
  v~dunif(0,8)
  rho<-exp(delta)
  } "
  
  
  
  
  
  
  sel.mod.bin.std.norm.ran<-"
  model{  \n## c. binary standard random-effects meta-analysis with normal prior 
  for (s in 1:S){
  y0[s]~dbin(pi0[s],n0[s])
  y1[s]~dbin(pi1[s],n1[s])
  logit(pi0[s])<-alpha[s]
  logit(pi1[s])<-alpha[s]+delta[s]
  
  alpha[s]~dnorm(0,0.0001)
  delta[s]~dnorm(mu,tau.delta)
  gamma[s]<-exp(delta[s])
  }
  ###prior###
  mu~dnorm(0,0.0001)
  sigma.delta~dunif(0,5)
  tau.delta<-pow(sigma.delta,-2)
  rho<-exp(mu)
  } "
  
  
  
  
  
  
  sel.mod.bin.std.dt.ran<-"
  model{  \n## d. binary standard random-effects meta-analysis with t-distribution prior
  for (s in 1:S){
  y0[s]~dbin(pi0[s],n0[s])
  y1[s]~dbin(pi1[s],n1[s])
  logit(pi0[s])<-alpha[s]
  logit(pi1[s])<-alpha[s]+delta[s]
  
  delta[s]~dnorm(mu,tau.delta)
  alpha[s]~dnorm(0,0.0001)
  gamma[s]<-exp(delta[s])
  }
  ###Prior###
  sigma.delta~dunif(0,5)
  tau.delta<-pow(sigma.delta,-2)
  mu~dt(0,0.5,v)
  v~dunif(0,8)
  rho<-exp(mu)
  } " 
  
  
  
  
  #########################################
  ## (ii) meta-regression for binary data
  
  if(!is.null(data$J)) {
    
    sel.mod.bin.reg.norm.fix<-if(data$J==1) { " 
      model {  \n## a. binary fixed-effects meta-regression with normal prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
      alpha[s]~dnorm(0,0.0001)
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      
      delta~dnorm(0,0.0001)
      rho<-exp(delta)
      } "
      
    } else {
      " 
      model {  \n## a. binary fixed-effects meta-regression with normal prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
      alpha[s]~dnorm(0,0.0001)
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      
      delta~dnorm(0,0.0001)
      rho<-exp(delta)
      } "   
      
    }   
    
    
    
    
    
    
    
    
    sel.mod.bin.reg.dt.fix<-if(data$J==1) {
      " 
      model {  \n## b. binary fixed-effects meta-regression with t-distribution prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
      alpha[s]~dnorm(0,0.0001)
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      
      delta~dt(0,0.5,v)
      v~dunif(0,8)
      rho<-exp(delta)
      } "
    } else {
      "model {  \n## b. binary fixed-effects meta-regression with t-distribution prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
      alpha[s]~dnorm(0,0.0001)
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      
      delta~dt(0,0.5,v)
      v~dunif(0,8)
      rho<-exp(delta)
    } "
    }
    
    
    
    sel.mod.bin.reg.norm.ran<-if(data$J==1) {
      "
      model {  \n## c. binary random-effects meta-regression with normal prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      alpha[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau.delta)
      gamma[s]<-exp(delta[s])
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      
      mu~dnorm(0,0.0001)
      sigma.delta~dunif(0,5)
      tau.delta<-pow(sigma.delta,-2)
      rho<-exp(mu)
      } "
    } else {
      "
      model {  \n## c. binary random-effects meta-regression with normal prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      alpha[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau.delta)
      gamma[s]<-exp(delta[s])
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      
      mu~dnorm(0,0.0001)
      sigma.delta~dunif(0,5)
      tau.delta<-pow(sigma.delta,-2)
      rho<-exp(mu)
      } "
    }
    
    
    
    
    sel.mod.bin.reg.dt.ran<-if(data$J==1) {
      " 
      model {  \n## d. binary random-effects meta-regression with t-distribution prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      delta[s]~dnorm(mu,tau.delta)
      alpha[s]~dnorm(0,0.0001)
      gamma[s]<-exp(delta[s])
      }
      ###Prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      
      sigma.delta~dunif(0,5)
      tau.delta<-pow(sigma.delta,-2)
      mu~dt(0,0.5,v)
      v~dunif(0,8)
      rho<-exp(mu)
      } "
    } else {
      "
      model {  \n## d. binary random-effects meta-regression with t-distribution prior
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      delta[s]~dnorm(mu,tau.delta)
      alpha[s]~dnorm(0,0.0001)
      gamma[s]<-exp(delta[s])
      }
      ###Prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      
      sigma.delta~dunif(0,5)
      tau.delta<-pow(sigma.delta,-2)
      mu~dt(0,0.5,v)
      v~dunif(0,8)
      rho<-exp(mu)
      } "
    }
    
  }
  
  
  
  #################################################### 
  ## (iii) standard meta-analysis for continuous data
  
  sel.mod.ctns.std.ta.fix<-"
  model{   \n## a. continuous fixed-effects meta-analysis with data available for two arms separately
  for(s in 1:S) { 
  y0[s]~dnorm(alpha0[s],prec0[s])       
  y1[s]~dnorm(alpha1[s],prec1[s])       
  
  alpha1[s]<-alpha0[s]+delta
  prec0[s]<-pow(se0[s],-2)
  prec1[s]<-pow(se1[s],-2)
  
  alpha0[s]~dnorm(0,0.0001)
  }
  
  ###prior###
  delta~dnorm(0,0.0001)
  
  }"
  
  
  
  
  sel.mod.ctns.std.ta.ran<-"
  model{  \n## b. continuous random-effects meta-analysis with data available for two arms separately
  for(s in 1:S) { 
  y0[s]~dnorm(alpha0[s],prec0[s])       
  y1[s]~dnorm(alpha1[s],prec1[s])       
  
  alpha1[s]<-alpha0[s]+delta[s]
  prec0[s]<-pow(se0[s],-2)
  prec1[s]<-pow(se1[s],-2)
  
  alpha0[s]~dnorm(0,0.0001)
  delta[s]~dnorm(mu,tau.delta)
  }
  
  ###prior###
  mu~dnorm(0,0.0001)
  tau.delta<-pow(sigma,-2)
  sigma~dunif(0,10)
  }"
  
  
  
  
  sel.mod.ctns.std.mv.fix<-"
  model{  \n## c. continuous fixed-effects meta-analysis for studies reporting mean difference and pooled variance
  for(s in 1:S){
  y[s]~dnorm(mu,prec[s])
  
  }
  ###prior###
  mu~dnorm(0,0.0001)
  }"
  
  
  
  
  
  
  sel.mod.ctns.std.mv.ran<-"
  model{  \n## d. continuous random-effects meta-analysis for studies reporting mean difference and pooled variance
  for(s in 1:S){
  y[s]~dnorm(delta[s],prec[s])        
  delta[s]~dnorm(mu,tau.delta) 
  
  }
  ###prior###
  mu~dnorm(0,0.0001)
  tau.delta<-pow(sigma,-2)
  sigma~dunif(0,10)
  }" 
  
  
  
  
  ############################################
  ## (iv) meta-regression for continuous data
  
  if(!is.null(data$J)) {
    sel.mod.ctns.reg.ta.fix<-if(data$J==1){ 
      "
      model{ \n## a. continuous fix-effects meta-regression with data available for two arms separately
      for(s in 1:S) { 
      y0[s]~dnorm(alpha0[s],prec0[s])       
      y1[s]~dnorm(alpha1[s],prec1[s])       
      
      alpha1[s]<-alpha0[s]+delta
      prec0[s]<-pow(se0[s],-2)
      prec1[s]<-pow(se1[s],-2)
      
      
      alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
      gamma0[s]~dnorm(0,0.0001)
      
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])      
      delta~dnorm(0,0.0001)
      }"
    } else { " 
      model{ \n## a. continuous fix-effects meta-regression with data available for two arms separately
      for(s in 1:S) { 
      y0[s]~dnorm(alpha0[s],prec0[s])       
      y1[s]~dnorm(alpha1[s],prec1[s])       
      
      alpha1[s]<-alpha0[s]+delta
      prec0[s]<-pow(se0[s],-2)
      prec1[s]<-pow(se1[s],-2)
      
      
      alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
      gamma0[s]~dnorm(0,0.0001)
      
      }
      
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])     
      delta~dnorm(0,0.0001)
      }"
    }
    
    
    
    sel.mod.ctns.reg.ta.ran<-if(data$J==1){
      "
      model{  \n## b. continuous random-effects meta-regression with data available for two arms separately
      for(s in 1:S) { 
      y0[s]~dnorm(alpha0[s],prec0[s])       
      y1[s]~dnorm(alpha1[s],prec1[s])       
      
      alpha1[s]<-alpha0[s]+delta[s]
      prec0[s]<-pow(se0[s],-2)
      prec1[s]<-pow(se1[s],-2)
      
      alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
      gamma0[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau.delta)
      
      }
      
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])     
      mu~dnorm(0,0.0001)
      tau.delta<-pow(sigma,-2)
      sigma~dunif(0,10)
      
      }" 
    } else {
      "model{ \n## b. continuous random-effects meta-regression with data available for two arms separately
      for(s in 1:S) { 
      y0[s]~dnorm(alpha0[s],prec0[s])       
      y1[s]~dnorm(alpha1[s],prec1[s])       
      
      alpha1[s]<-alpha0[s]+delta[s]
      prec0[s]<-pow(se0[s],-2)
      prec1[s]<-pow(se1[s],-2)
      
      alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
      gamma0[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau.delta)
      
      }
      
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])      
      mu~dnorm(0,0.0001)
      tau.delta<-pow(sigma,-2)
      sigma~dunif(0,10)
      
    }"   
    }  
    
    
    
    
    sel.mod.ctns.reg.mv.fix<-if(data$J==1){
      "
      model{ \n## d. continuous fixed-effects meta-regression for studies reporting mean difference and pooled variance 
      for(s in 1:S){
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha+Z0[s,]%*%beta0
      
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])       
      alpha~dnorm(0,0.0001)
      }"
    } else {
      "model{ \n## d. continuous fixed-effects meta-regression for studies reporting mean difference and pooled variance
      for(s in 1:S){ 
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha+Z0[s,]%*%beta0
      
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])       
      alpha~dnorm(0,0.0001)
      
    }"
    }
    
    
    
    
    sel.mod.ctns.reg.mv.ran<-if(data$J==1){
      "
      model{ \n## d. continuous random-effects meta-regression for studies reporting mean difference and pooled variance
      for(s in 1:S){
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha[s]+Z0[s,]%*%beta0
      alpha[s]~dnorm(mu,tau.alpha)
      
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])       
      mu~dnorm(0,0.0001)
      tau.alpha<-pow(sigma,-2)
      sigma~dunif(0,10)
      }"
    } else {
      "model{  \n## d. continuous random-effects meta-regression for studies reporting mean difference and pooled variance
      for(s in 1:S){
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha[s]+Z0[s,]%*%beta0
      alpha[s]~dnorm(mu,tau.alpha)
      
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])       
      mu~dnorm(0,0.0001)
      tau.alpha<-pow(sigma,-2)
      sigma~dunif(0,10)
      
    }" 
    }
  }
  
  
  
  
  ##############################################
  ## (v) standard meta-analysis for count data
  
  sel.mod.count.std.fix<-"
  model{ \n## a. count fixed-effects meta-analysis 
  for(s in 1:S){
  y0[s]~dpois(lambda0[s])
  y1[s]~dpois(lambda1[s])
  
  log(lambda0[s])<-lograte0[s]+log(pt0[s])
  log(lambda1[s])<-lograte1[s]+log(pt1[s])
  
  lograte0[s]~dunif(-5,5)
  lograte1[s]<-lograte0[s]+delta      ### delta is the log rate-ratio
  
  }
  ### Prior ###
  delta~dnorm(0,0.0001)
  IRR<-exp(delta)
  }"
  
  
  
  
  
  sel.mod.count.std.unif.ran<-"
  model{ \n## b. count random-effects meta-analysis with uniform prior
  for(s in 1:S){
  y0[s]~dpois(lambda0[s])
  y1[s]~dpois(lambda1[s])
  
  log(lambda0[s])<-lograte0[s]+log(pt0[s])
  log(lambda1[s])<-lograte1[s]+log(pt1[s])
  
  lograte0[s]~dunif(-5,5)
  lograte1[s]<-lograte0[s]+delta[s]
  delta[s]~dnorm(mu,tau.delta)
  gamma[s]<-exp(delta[s])
  }
  ### Prior ###
  mu~dnorm(0,0.0001)
  tau.delta<-pow(sigma,-2)
  sigma~dunif(0,10)
  IRR<-exp(mu)
  }"
  
  
  
  
  sel.mod.count.std.hc.ran<-"
  model{ \n## c. count random-effects meta-analysis with halfcauchy prior
  for(s in 1:S){
  y0[s]~dpois(lambda0[s])
  y1[s]~dpois(lambda1[s])
  
  log(lambda0[s])<-lograte0[s]+log(pt0[s])
  log(lambda1[s])<-lograte1[s]+log(pt1[s])
  
  lograte0[s]~dunif(-5,5)
  lograte1[s]<-lograte0[s]+delta[s]
  delta[s]~dnorm(mu,tau.delta)
  gamma[s]<-exp(delta[s])
  }
  ### Prior ###
  mu~dnorm(0,0.0001)
  tau.delta<-pow(sigma,-2)
  
  sigma<-abs(z.sigma)/pow(epsilon.sigma,0.5)
  z.sigma~dnorm(0,tau.z.sigma)
  epsilon.sigma~dgamma(0.5,0.5)
  tau.z.sigma<-pow(B.sigma,-2)
  B.sigma~dunif(0,0.5)
  
  IRR<-exp(mu) 
  }"
  
  
  
  
  
  
  ##########################################
  ## (vi) meta-regression for count data 
  
  if(!is.null(data$J)) {
    sel.mod.count.reg.fix<-if(data$J==1){
      "
      model{ \n## a. count fixed-effects meta-regression
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-lograte0[s]+log(pt0[s])
      log(lambda1[s])<-lograte1[s]+log(pt1[s])
      
      lograte0[s]<-logbrate[s]+Z0[s,]%*%beta0   
      logbrate[s]~dunif(-5,5)
      lograte1[s]<-lograte0[s]+delta
      
      }
      ### Prior ###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      delta~dnorm(0,0.0001)
      
      IRR<-exp(delta)
      }"
    } else {"
      model{ \n## a. count fixed-effects meta-regression
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-lograte0[s]+log(pt0[s])
      log(lambda1[s])<-lograte1[s]+log(pt1[s])
      
      lograte0[s]<-logbrate[s]+Z0[s,]%*%beta0   
      logbrate[s]~dunif(-5,5)
      lograte1[s]<-lograte0[s]+delta
      
      }
      ### Prior ###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      delta~dnorm(0,0.0001)
      IRR<-exp(delta)
      
      }" 
    }
    
    
    
    
    
    sel.mod.count.reg.unif.ran<-if(data$J==1){
      "
      model{ \n## b. count random-effects meta-regression with uniform prior
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-lograte0[s]+log(pt0[s])
      log(lambda1[s])<-lograte1[s]+log(pt1[s])
      
      lograte0[s]<-logbrate[s]+Z0[s,]%*%beta0   
      logbrate[s]~dunif(-5,5)
      lograte1[s]<-lograte0[s]+delta[s]
      delta[s]~dnorm(mu,tau.delta)
      gamma[s]<-exp(delta[s])
      }
      ### Prior ###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,]) 
      mu~dnorm(0,0.0001)
      tau.delta<-pow(sigma,-2)
      sigma~dunif(0,10)
      IRR<-exp(mu)
      }"
    } else {
      "model{ \n## b. count random-effects meta-regression with uniform prior
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-lograte0[s]+log(pt0[s])
      log(lambda1[s])<-lograte1[s]+log(pt1[s])
      
      lograte0[s]<-logbrate[s]+Z0[s,]%*%beta0   
      logbrate[s]~dunif(-5,5)
      lograte1[s]<-lograte0[s]+delta[s]
      delta[s]~dnorm(mu,tau.delta)
      gamma[s]<-exp(delta[s])
      }
      ### Prior ###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      mu~dnorm(0,0.0001)
      tau.delta<-pow(sigma,-2)
      sigma~dunif(0,10)
      IRR<-exp(mu)
    }" 
    }
    
    
    
    
    
    sel.mod.count.reg.hc.ran<-if(data$J==1){"
      model{ \n## c. count random-effects meta-regression with halfcauchy prior
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-lograte0[s]+log(pt0[s])
      log(lambda1[s])<-lograte1[s]+log(pt1[s])
      
      lograte0[s]<-logbrate[s]+Z0[s,]%*%beta0   
      logbrate[s]~dunif(-5,5)
      lograte1[s]<-lograte0[s]+delta[s]
      delta[s]~dnorm(mu,tau.delta)
      gamma[s]<-exp(delta[s])                                                                                       
      }
      ### Prior ###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,]) 
      mu~dnorm(0,0.0001)
      tau.delta<-pow(sigma,-2)
      
      sigma<-abs(z.sigma)/pow(epsilon.sigma,0.5)
      z.sigma~dnorm(0,tau.z.sigma)
      epsilon.sigma~dgamma(0.5,0.5)
      tau.z.sigma<-pow(B.sigma,-2)
      B.sigma~dunif(0,0.5)
      IRR<-exp(mu)
      
      }"
    } else {"
      model{ \n## c. count random-effects meta-regression with halfcauchy prior
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-lograte0[s]+log(pt0[s])
      log(lambda1[s])<-lograte1[s]+log(pt1[s])
      
      lograte0[s]<-logbrate[s]+Z0[s,]%*%beta0   
      logbrate[s]~dunif(-5,5)
      lograte1[s]<-lograte0[s]+delta[s]
      delta[s]~dnorm(mu,tau.delta)
      gamma[s]<-exp(delta[s])
      }
      ### Prior ###
      
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      mu~dnorm(0,0.0001)
      tau.delta<-pow(sigma,-2)
      
      sigma<-abs(z.sigma)/pow(epsilon.sigma,0.5)
      z.sigma~dnorm(0,tau.z.sigma)
      epsilon.sigma~dgamma(0.5,0.5)
      tau.z.sigma<-pow(B.sigma,-2)
      B.sigma~dunif(0,0.5)
      
      IRR<-exp(mu)         
      
      }"  
    }
    
  }
  
  
  
  
  
  
  ###############################################################
  ## Combines the required modules in a singel model code ##
  mod.outcome<-c("bin","ctns","count")
  mod.sel<-c("std.norm","std.dt","reg.norm","reg.dt",
             "std.ta","std.mv","reg.ta","reg.mv",
             "std","std.unif","std.hc","reg","reg.unif","reg.hc")
  mod.type<-c("fix","ran")
  
  lab.outcome<-c("Model for binary data ", "Model for continuous data ", "Model for count data ")
  lab.sel<-c("Meta-analysis with normal prior ","Meta-analysis with t-distribution prior ",
             "Meta-regression with normal prior ", "Meta-regression with t-distribution prior ",
             "Meta-analysis with data avalaible for two arms separately ", 
             "Meta-analysis for studies reporting mean difference and pooled variance ",
             "Meta-regression with data avalaible for two arms separately ", 
             "Meta-regression for studies reporting mean difference and pooled variance ",
             "Meta-analysis ","Meta-analysis with uniform prior ","Meta-analysis with halfcauchy prior ",
             "Meta-regression ","Meta-regression with uniform prior ","Meta-regression with halfcauchy prior ")
  lab.type<-c("Fixed-effects ", "Random-effects ")
  time <- Sys.time()
  
  for (mo in 1:length(mod.outcome)){
    for (ms in 1:length(mod.sel)){
      for (mt in 1:length(mod.type)){
        if (outcome==mod.outcome[mo] & model==mod.sel[ms] & type==mod.type[mt]){
          summary.text<- paste0("#",lab.outcome[mo],lab.type[mt],lab.sel[ms])
          print.time <- paste0("#",time)
          txt<- paste0("summary.text, print.time, sel.mod.", mod.outcome[mo], ".", mod.sel[ms], ".", mod.type[mt])  
          cmd<- paste0("model<- paste(",txt,")")
          eval(parse(text=cmd))
        }
      }
    }
  }
  
  filein <- model.file
  writeLines(model, con=filein)
  
}
