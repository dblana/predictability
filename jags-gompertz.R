# Sat Jul 18 17:57:49 2020 ------------------------------
# Mario Castro

require(rjags)

# Some auxuliary functions
reverse <- function(var) var[seq(length(var),1,-1)] # reverse an array (Auxilirary function)
colMed <- function(X) apply(X, 2, median) # Compute the median of a matrix, by column (similar to the standard colMeans)
colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Compute de quantile "q" of a matrix, by column

# delta.peak: distance to the empirical peak used to fit the model
# plot.flag: Boolean, to plot or not to plot
# save.plot: Boolean, save or not the second plot
# prior: log(2) as in any country the doubling time has been larger than 1 day (non-inf would be, for instance 10, 1000, ...
fit.gompertz <- function(delta.tpeak=-7,plot.flag=TRUE,save.plot=TRUE,prior=0.6931472)  {
  
covid.es <- read.csv('datasets/covid-19-es.csv') # Load dataset
y <- covid.es$Confirmed # Number of confirmed cases
day <- covid.es$Day # Day in calendar

# As the initial numbers in the series are too small, we remove initial reports with less than 5 cases per day.
# These are 26 data points full of repeated 1's and 2's and removing them do not affect the fits.
idx <- y>10 # Comment if you want to include them all.
y <- y[idx]  # Comment if you want to include them all.
day <- day[idx]

N <- length(y) # Number of obsrevations

tpeak <- 47 # This corresponds to the 13th of April.
tmax <- tpeak + delta.tpeak # Last point used to fit the data (-7: a week before, 14: two weaks later, ...) 
tf <- N-2 # Last day to project the data
time <- (1:tf) # Vector of times

# Barplot with the data
if(plot.flag==TRUE) {
  x11(width=16,height=7) # Create new window for plotting
  par(mar=c(7.1,5.1,4.1,2.1)) # Change margins
  barplot(height=y,names=day,las=2,col='skyblue',border = F,main='Total Number of Confirmed cases (Spain)')
}

data <- list(y=y,tmax=tmax,tf=tf,prior=prior) # Build a list for JAGS

modelstring="
model {
  for(t in 1:tmax) {
    y[t] ~ dnorm(K*exp(-b/c*exp(-c*(t-1))),tau) # (t-1) because jags arrays start at [1]
  }  
  for(t in 1:tf) {
    ypred[t] ~ dnorm(K*exp(-b/c*exp(-c*(t-1))),tau) # Posterior predictive
  }   

  b ~ dunif(0,prior)
  c ~ dunif(0,10) # Non-informative prior
  m ~ dunif(3,8) # Order of magnitud of steady state: From 1000 to 10^8 total infected (Spanish population: 45*10^6)
  K <- 10^m # Gompertz steady-state parameter
  tau ~ dgamma(.01,.01) # Precision of the gaussian: tau = 1/variance
}
"
model=jags.model(textConnection(modelstring), data=data,n.chains = 3,n.adapt=5000) # burn-in=5000
update(model,n.iter=15000) # Update MCMC
output=coda.samples(model=model,variable.names=c("b","K","c","tau",'ypred'), n.iter=15000, thin=1) #Sample MCMC

mat <- as.matrix(output) # Extract MCMC data in matrix form 


# Plot the results
if(plot.flag==TRUE) {
  x11(width=9,height=7) # Create new window for plotting
  par(mar=c(5.1,5.1,4.1,2.1))
  yp <- mat[,-(1:4)] # Posterior predictive (ypred[1] ...ypred[90])
  plot((data$y),xlim=c(0,tf),ylim=c(0,325000),pch=19,cex=.5,
       ylab='Total number of confirmed cases',xlab='Days since 27th of february',lwd=2,
       cex.lab=1.5,cex.axis=1.5)
  
  yM <- colq(yp,.975) # 97.5% quantile
  ym <- colq(yp,0.0225) # 2.5% quantile
  tt <- c(time,reverse(time))
  yy <- c(yM,reverse(ym))
  polygon(tt,yy,col=rgb(0,0,0,.2),border=rgb(0,0,0,.2),lwd=3)
  points((data$y[1:tmax]),col=rgb(.33,.55,.67,1),cex=1.5,pch=19) # Plot only used points in bayesian regression
  lines((colMed(yp)),col=rgb(1,.55,0,1),lwd=4) # Median of posterior predictive
  # Estimate theoretical point of inflection for the Gompertz curve
  median.b <- median(mat[,"b"]) # Extract median value of parameter "b"
  median.c <- median(mat[,"c"]) # Extract median value of parameter "c"
  print(gompertz.peak <- -log(median.c/median.b)/median.c) # Theoretical inflection point
  
  # Add some annotations to the plot
  arrows(75,190000,70,y[70],angle=10,lwd=3)
  text(70,180000,'After the 6th of may',pos=4,cex=1.5)
  text(70,160000,'there is a new',pos=4,cex=1.5)
  text(70,140000,'inflection point',pos=4,cex=1.5)
  abline(v=tpeak,lty=2,col=rgb(1,.55,0,1))
  text(tpeak,50000,"Real inflection ",pos=4,col=rgb(1,.55,0,1),cex=1.5)
  text(tpeak,30000,"point (epidemic  'peak')",pos=4,col=rgb(1,.55,0,1),cex=1.5)
  abline(v=gompertz.peak,lwd=2,lty=3,col=rgb(.33,.55,.67,1),cex=1.5)
  text(gompertz.peak,150000,'Gompertz',pos=2,col=rgb(.33,.55,.67,1),cex=1.5)
  text(gompertz.peak,130000,'inflection point',pos=2,col=rgb(.33,.55,.67,1),cex=1.5)
  # Export figure?
  if(save.plot==TRUE) dev.copy2pdf(file=sprintf("output/gompertz-fit_%+d_prior_%f.pdf",delta.tpeak,prior))
  }
if(plot.flag==TRUE) { # Plot posterior distributions of the model parameters
  x11(width=9,height=7) # Create new window for plotting
  par(mfrow=c(2,2)) # A 2x2 plot 
  par(mar=c(5.1,5.1,4.1,2.1)) # Change margins

  # Plot posteriors for model parameters
  hist(mat[,"b"],xlab="b",main="",col="skyblue",probability = T,border=FALSE)
  hist(mat[,"c"],xlab="c",main="",col="skyblue",probability = T,border=FALSE)
  hist(mat[,"K"],xlab="K",main="",col="skyblue",probability = T,border=FALSE)
  hist(1/sqrt(mat[,"tau"]),xlab=expression(sigma),main="",col="skyblue",probability = T,border=FALSE)
  print(gelman.diag(output[,1:4])) # Check if the MCMC chains converged
  }
}

# Run example with default values
fit.gompertz()
# Run more examples
fit.gompertz(prior = 10)
fit.gompertz(delta.tpeak = 7)
fit.gompertz(delta.tpeak = 7,prior=10)
fit.gompertz(delta.tpeak = 21)
fit.gompertz(delta.tpeak = 21,prior=10)
