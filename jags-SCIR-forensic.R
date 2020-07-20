# Sat Jul 18 17:57:49 2020 ------------------------------
# Mario Castro

require(rjags)

# Some auxiliary functions
reverse <- function(var) var[seq(length(var),1,-1)] # reverse an array (Auxilirary function)
colMed <- function(X) apply(X, 2, median) # Compute the median of a matrix, by column (similar to the standard colMeans)
colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Compute de quantile "q" of a matrix, by column

# Calculate the analytical solution of SCIR model
exact.SCIR.forensic <- function(param,data){
  I0 <- data$I0
  Iq <- data$Iq
  Iq2 <- data$Iq2
  t0 <- data$t0
  tq <- data$tq
  tq2 <- data$tq2
  tf <- data$tf
  beta <- param$beta
  rmu <- param$rmu
  p <- param$p
  q <- param$q
  p2 <- param$p2
  q2 <- param$q2
  
  t <- (t0+1):(tq-1) # From second day to first quarantine 
  out <- c(I0,I0+(beta-rmu)*(t-t0))
  tout <- t
  t <- (tq+1):tq2 # From day after first quarantine to second quarantine
  out <- c(out,Iq+ ((beta*q)/(p+q)^2*(1-exp(-(p+q)*(t-tq)))+(beta-rmu-beta*q/(q+p))*(t-tq)))
  tout <- c(tout,tq,t)
  t <- (tq2+1):tf # From day after first quarantine to second quarantine
  out <- c(out,Iq2+ ((beta*q2)/(p2+q2)^2*(1-exp(-(p2+q2)*(t-tq2)))+(beta-rmu-beta*q2/(q2+p2))*(t-tq2)))
  tout <- c(tout,t)-1
  return(data.frame(tout,out))
}

# Plot prediction
plot.SCIR.forensic.output <- function(output,data) {
  x11("",8,7) # Create new window for plotting
  par(mar=c(5.1,5.1,4.1,2.1)) # Change margins
  
  mat.output <- as.matrix(output) # Convert MCMC samples to matrix
  t0 <- data$t0
  tq <- data$tq
  tq2 <- data$tq2
  tmax <- data$tmax
  tf <- data$tf
  ty <- (t0):tf # Vector of times (Days since first confirmed case)
  m <- mat.output[,-(1:8)] # Remove first 8 columns (corresponding to model parameters)
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
  
  param <- data.frame(t(colMed(mat.output[,1:6]))) # Median parameters from posteriors ("beta","p","q","rmu")
  # Create a vector with the exact solution of the SCIR model for the median parameters
  exact<- exact.SCIR.forensic(param,data)
  t.median <- exact$tout
  y.median <- exact$out
  lines(t.median,y.median/log(10),col='darkorange',lwd=4)
  
  arrows(tq,5,tq,3.3,angle=10,lwd=2) # Annotate the plot
  text(tq,5,"1st Confiment begins",pos=3)
  arrows(tq2,7,tq2,5,angle=10,lwd=2) # Annotate the plot
  text(tq2,7,"2nd Confiment begins",pos=3)
}

plot.posteriors.forensic <- function(output) {
  x11("",8,7) # Create new window for plotting
  
  mat <- as.matrix(output)
  ############ PLOTTING ###########
  par(mfrow=c(3,3))
  par(mar=c(5.1,5.1,4.1,2.1))
  hist(mat[,"beta"],xlab=expression(beta),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"rmu"],20,xlab=expression(r+mu),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"q"],20,xlab=expression(q),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"p"],xlab=expression(p),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"q2"],20,xlab=expression(q[2]),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(mat[,"p2"],xlab=expression(p[2]),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(1/sqrt(mat[,"tauI"]),xlab=expression(sigma[I]),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  hist(1/sqrt(mat[,"tauX"]),xlab=expression(sigma[X]),main="",col="skyblue",probability=T,cex.axis=1.5,cex.lab=2,border=F)
  
  gd <- gelman.diag(output[,1:5])
  print(gd$mpsrf)
  
  title(sprintf("Gelman diagnostic: %.2f",gd$mpsrf), line = -2, outer = TRUE,cex.main=1.5)
  return(gd$mpsrf)
}

# tmax: Last data-point used for estimation
# plot.flag: To plot or not to plot
# save.plot: To export figure or not
scir.bayesian.forensic <- function(tmax=34,plot.flag=TRUE,save.plot=TRUE) {
  # Read datasets
  covid.es <- read.csv('datasets/covid-19-es.csv') # Load dataset
  confirmed <- covid.es$Confirmed # Number of confirmed cases
  recovered <- covid.es$Recovered
  deaths <- covid.es$Death
  day <- covid.es$Day # Day in calendar
  
  idx <- confirmed>10 # Comment if you want to include them all.
  recovered <- recovered[idx]  # Comment if you want to include them all.
  deaths <- deaths[idx]  # Comment if you want to include them all.
  confirmed <- confirmed[idx]  # Comment if you want to include them all.
  tf <- which(is.na(recovered))[1]-1 # After the 18th of May the Spanish government stopped publishing data about recovered
  day <- day[idx][1:tf]

  I <- log(confirmed-recovered-deaths)[1:tf] # Active cases in logarithmic scale
  X <- c(0,log(diff(deaths+recovered)))[1:tf]  # New deaths+recovered (daily derivative) in logarithmic scale 
  I[is.infinite(I) | is.na(I)] <- 0 # Remove infinites and NaNs
  X[is.infinite(X) | is.na(X)] <- 0 # Remove infinites and NaNs
  
  t0 <- which(I>0)[1] # First day with more than 1 confirmed case
  tX0 <- which(X>0)[1] # First day with more than 1 new death or recovered
  
  tq <- max(11,t0+2) # Begin quarantine (9 March): schools closures are announced 
  tq2 <- 33 # Begin second quarantine (30 March)
  
  if(plot.flag==TRUE) {
    x11(width=16,height=7)
    par(mar=c(7.1,5.1,4.1,2.1))
    barplot(height=confirmed,names=covid.es$Day[idx],las=2,col='skyblue',border = F,
            main='Total Number of Confirmed cases (Spain)')
  }
  
  # Build a list for JAGS
  # I0: First datum in the Active cases series
  # Iq: First datum after quarantine in the Active cases series
  data <- list(I=I,X=X,tq=tq,tq2=tq2,tf=tf,t0=t0,tX0=tX0,tmax=tmax,I0=I[t0],Iq=I[tq],Iq2=I[tq2]) 
  
  modelstring="
  # JAGS model
  model {
  # Regression before first global confinement (populations are in log scale)
  for(t in (t0+1):(tq)) {
  I[t] ~ dnorm(I0+(beta-rmu)*(t-t0),tauI) # Active cases
  y[t] ~ dnorm(I0+(beta-rmu)*(t-t0),tauI) # Posterior predictive   
  }
  # Regression for active cases after the first confinement (populations are in log scale)
  for(t in (tq+1):tq2) {
  I[t] ~ dnorm(Iq+ ((beta*q)/(p+q)^2*(1-exp(-(p+q)*(t-tq)))+ (beta-rmu-beta*q/(q+p))*(t-tq)),tauI)
  y[t] ~ dnorm(y[tq]+ ((beta*q)/(p+q)^2*(1-exp(-(p+q)*(t-tq)))+ (beta-rmu-beta*q/(q+p))*(t-tq)),tauI) # posterior predictive
  }
  # Regression for active cases after the second confinement (populations are in log scale)
  for(t in (tq2+1):tmax) {
  I[t] ~ dnorm(Iq2+ ((beta*q2)/(p2+q2)^2*(1-exp(-(p2+q2)*(t-tq2)))+ (beta-rmu-beta*q2/(q2+p2))*(t-tq2)),tauI)
  }
  # Posterior predictive for active cases after the second confinement (extended until tf)
  for(t in (tq2+1):tf) {
  y[t] ~ dnorm(y[tq2]+ ((beta*q2)/(p2+q2)^2*(1-exp(-(p2+q2)*(t-tq2)))+ (beta-rmu-beta*q2/(q2+p2))*(t-tq2)),tauI)
  }
  # Regression for new deaths+recovered cases
  for(t in tX0:tmax) {
  X[t] ~ dnorm(log(rmu)+I[t],tauX)
  }
  # Posterior predictive for new deaths+recovered (extended until tf)
  for(t in (tX0):tf) {
  z[t] ~ dnorm(log(rmu)+y[t],tauX) # New Death + Recovered
  }
  # Priors for parameters
  p ~ dunif(0,5) 
  q ~ dunif(0,5) 
  p2 ~ dunif(0,5)
  q2 ~ dunif(0,5)
  beta ~ dunif(0,1) # Doubling time is less than 1 per day
  rmu ~ dunif(0,1) # rmu is lower than beta (so R0>1)
  # Priors for precision (inverse of variance)
  tauI ~ dgamma(0.01,0.01) # Non-informative prior
  tauX ~ dgamma(0.01,0.01) # Non-informative prior
  y[t0] <- I0
  }
  "
  
  model=jags.model(textConnection(modelstring), data=data,n.chains = 3,n.adapt = 5000) # burn-in=5000
  update(model,n.iter=10000) # Update MCMC 
  output.scir.forensic=coda.samples(model=model,variable.names=c("p","q","p2","q2","z","y","rmu","tauI","tauX","beta"), 
                           n.iter=10000, thin=1) # Extract MCMC data in matrix form 
  
  
  
  if(plot.flag==TRUE) {
    plot.SCIR.forensic.output(output.scir.forensic,data)
    if(save.plot==TRUE)  
      dev.copy2pdf(file='output/forensic-bayesian-SCIR-fit.pdf')
    plot.posteriors.forensic(output.scir.forensic)
    if(save.plot==TRUE)  
      dev.copy2pdf(file='output/forensic-posteriors-SCIR.pdf')
      
  }
  
  return(list(output.scir.forensic=output.scir.forensic,data=data)) # Return MCMC samples and data in a list
}
simulation.forensic <- scir.bayesian.forensic(tmax = 33 ) # Run code 

#saveRDS(simulation,'output/simulation.rda') # Uncomment to save all into a binary file
