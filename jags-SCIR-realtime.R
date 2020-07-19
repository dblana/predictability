# Sat Jul 18 17:57:49 2020 ------------------------------
# Mario Castro

require(rjags)

# Some auxiliary functions
reverse <- function(var) var[seq(length(var),1,-1)] # reverse an array (Auxilirary function)
colMed <- function(X) apply(X, 2, median) # Compute the median of a matrix, by column (similar to the standard colMeans)
colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Compute de quantile "q" of a matrix, by column

# Calculate the analytical solution of SCIR model
exact.SCIR <- function(param,data){
  I0 <- data$I0
  Iq <- data$Iq
  t0 <- data$t0
  tq <- data$tq
  tf <- data$tf
  beta <- param$beta
  rmu <- param$rmu
  p <- param$p
  q <- param$q
  
  t <- (t0+1):(tq-1) # From second day to quarantine 
  out <- c(I0,I0+(beta-rmu)*(t-t0))
  tout <- t
  t <- (tq+1):tf # From day after quarantine to the end
  out <- c(out,Iq+ ((beta*q)/(p+q)^2*(1-exp(-(p+q)*(t-tq)))+(beta-rmu-beta*q/(q+p))*(t-tq)))
  tout <- c(tout,tq,t)-1
  return(data.frame(tout,out))
}

# Plot prediction
plot.SCIR.output <- function(output,data) {
  x11("",8,7) # Create new window for plotting
  par(mar=c(5.1,5.1,4.1,2.1)) # Change margins
  
  mat.output <- as.matrix(output) # Convert MCMC samples to matrix
  t0 <- data$t0
  tq <- data$tq
  tmax <- data$tmax
  tf <- data$tf
  ty <- (t0):tf # Vector of times (Days since first confirmed case)
  m <- mat.output[,-(1:6)] # Remove first 6 columns (corresponding to model parameters)
  y <- m[,1:length(ty)] # Posterior predictive for active cases: y[1]..y[tf]
  
  # Create a plot with the original data (converted to log10 scale)
  plot(ty[t0:tmax],data$I[t0:tmax]/log(10),ylim=c(1,8),xlim=c(1,tf),
       xlab='Days since first confirmed case',
       ylab=expression(paste('Active cases (',log[10],')')),
       cex.axis=1.75,cex.lab=2,col='black',las=1); 
  
  # Auxiliary variables to create a polygon for the 95% posterior density
  tt <- c(ty,reverse(ty)) # Patch times and reversed times
  yM <- colq(y,0.975)/log(10) # 97.5% quantile (and convert to log10 scale)
  ym <- colq(y,0.0225)/log(10)# 2.5% quantile
  yy <- c(yM,reverse(ym)) # patch 97.5 and (reversed) 2.5 quantiles
  polygon(tt,yy,col=rgb(0,0,0,.2),border='gray')
  
  param <- data.frame(t(colMed(mat.output[,1:4]))) # Median parameters from posteriors ("beta","p","q","rmu")
  # Create a vector with the exact solution of the SCIR model for the median parameters
  exact<- exact.SCIR(param,data)
  t.median <- exact$tout
  y.median <- exact$out
  lines(t.median,y.median/log(10),col='darkorange',lwd=4)
  
  arrows(tq,5,tq,3,angle=10,lwd=2) # Annotate the plot
  text(tq,5,"Confiment begins",pos=3)
}

plot.posteriors<- function(output) {
  x11("",8,7) # Create new window for plotting
  
  mat <- as.matrix(output)
  ############ PLOTTING ###########
  par(mfrow=c(3,2))
  par(mar=c(5.1,5.1,4.1,2.1))
  hist(mat[,"beta"],xlab=expression(beta),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"rmu"],20,xlab=expression(r+mu),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"q"],20,xlab=expression(q),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"p"],xlab=expression(p),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(1/sqrt(mat[,"tauI"]),xlab=expression(sigma[I]),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(1/sqrt(mat[,"tauX"]),xlab=expression(sigma[X]),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  
  gd <- gelman.diag(output[,1:5])
  print(gd$mpsrf)
  
  title(sprintf("Gelman diagnostic: %.2f",gd$mpsrf), line = -52, outer = TRUE,cex.main=1.5)
  return(gd$mpsrf)
}

# tmax: Last data-point used for estimation
# plot.flag: To plot or not to plot
# save.plot: To export figure or not
scir.bayesian <- function(tmax=33,plot.flag=TRUE,save.plot=TRUE) {
  # Read datasets
  Day <- read.csv('datasets/confirmed-march31.csv')$Day
  confirmed <- read.csv('datasets/confirmed-march31.csv')$Spain
  recovered <- read.csv('datasets/recovered-march31.csv')$Spain
  deaths <- read.csv('datasets/deaths-march31.csv')$Spain
  
  I <- log(confirmed-recovered-deaths) # Active cases in logarithmic scale
  X <- c(0,log(diff(deaths+recovered)))  # New deaths+recovered (daily derivative) in logarithmic scale 
  I[is.infinite(I) | is.na(I)] <- 0 # Remove infinites and NaNs
  X[is.infinite(X) | is.na(X)] <- 0 # Remove infinites and NaNs
  
  t <- 1:length(I)-1 # Vector of times (0 to length of data-1)
  tf <- 90 # Last day to project the data
  t0 <- which(I>0)[1] # First day with more than 1 confirmed case
  tX0 <- which(X>0)[1] # First day with more than 1 new death or recovered
  
  tq <- max(11,t0+2) # Begin quarantine (9 March): schools closures are announced 
  
  if(plot.flag==TRUE) {
    x11(width=16,height=7)
    par(mar=c(7.1,5.1,4.1,2.1))
    barplot(height=I,names=Day,las=2,col='skyblue',border = F,
            main='Total Number of Confirmed cases (Spain)')
  }
  
  # Build a list for JAGS
  # I0: First datum in the Active cases series
  # Iq: First datum after quarantine in the Active cases series
  data <- list(I=I,X=X,tq=tq,tf=tf,t0=t0,tX0=tX0,tmax=tmax,I0=I[t0],Iq=I[tq]) 
  
  modelstring="
  # JAGS model
  model {
  # Regression before global confinement (populations are in log scale)
  for(t in (t0+1):(tq)) {
  I[t] ~ dnorm(I0+(beta-rmu)*(t-t0),tauI) # Active cases
  y[t] ~ dnorm(I0+(beta-rmu)*(t-t0),tauI) # Posterior predictive   
  }
  # Regression for active cases post-confinement (populations are in log scale)
  for(t in (tq+1):tmax) {
  I[t] ~ dnorm(Iq+ ((beta*q)/(p+q)^2*(1-exp(-(p+q)*(t-tq)))+ (beta-rmu-beta*q/(q+p))*(t-tq)),tauI)
  }
  # Posterior predictive for active cases post-confinement
  for(t in (tq+1):tf) {
  y[t] ~ dnorm(y[tq]+ ((beta*q)/(p+q)^2*(1-exp(-(p+q)*(t-tq)))+(beta-rmu-beta*q/(q+p))*(t-tq)),tauI)
  }
  # Regression for new death+recovered cases
  for(t in tX0:tmax) {
  X[t] ~ dnorm(log(rmu)+I[t],tauX)
  }
  for(t in (tX0):tf) {
  z[t] ~ dnorm(log(rmu)+y[t],tauX) # New Death + Recovered
  }
  # Priors for parameters
  p ~ dunif(0,5)
  q ~ dunif(0,5)
  beta ~ dunif(0,1)
  rmu ~ dunif(0,1)
  # Priors for precissions (inverse of variance)
  tauI ~ dgamma(.01,.01)
  tauX ~ dgamma(0.01,0.01)
  y[t0] <- I0
  }
  "
  
  model=jags.model(textConnection(modelstring), data=data,n.chains = 3,n.adapt = 5000) # burn-in=5000
  update(model,n.iter=10000) # Update MCMC 
  output.scir=coda.samples(model=model,variable.names=c("p","q","z","y","rmu","tauI","tauX","beta"), 
                           n.iter=10000, thin=1) # Extract MCMC data in matrix form 
  
  
  
  if(plot.flag==TRUE) {
    plot.SCIR.output(output.scir,data)
    plot.posteriors(output.scir)
    if(save.plot==TRUE)  {
      dev.copy2pdf(file='output/bayesian-SCIR-fit.pdf')
      dev.copy2pdf(file='output/posteriors-SCIR.pdf')
      }
  }
  
  return(list(output.scir=output.scir,data=data)) # Return MCMC samples and data in a list
}
simulation <- scir.bayesian() # Run code

#saveRDS(simulation,'output/simulation.rda') # Uncomment to save all into a binary file
