# Sat Jul 18 17:57:49 2020 ------------------------------
# Mario Castro

require(rjags)

# Some auxuliary functions
reverse <- function(var) var[seq(length(var),1,-1)] # reverse an array (Auxilirary function)
colMed <- function(X) apply(X, 2, median) # Compute the median of a matrix, by column (similar to the standard colMeans)
colq <- function(X,q) apply(X, 2, function(X) as.numeric(quantile(X,q))) # Compute de quantile "q" of a matrix, by column

scir.bayesian <- function(tmax=33,plot.flag=TRUE) {
# Read datasets
confirmed <- read.csv('datasets/confirmed-march31.csv')
recovered <- read.csv('datasets/recovered-march31.csv')
deaths <- read.csv('datasets/deaths-march31.csv')

print(names <- colnames(read.csv('datasets/names.csv'))) # Print regions (21 for all Spain aggregated)
region <- 21 
name <- names[region]

I <- log(confirmed[,region]-recovered[,region]-deaths[,region]) # Active cases in logarithmic scale
X <- c(0,log(diff(deaths[,region]+recovered[,region])))  # New deaths+recovered (daily derivative) in logarithmic scale 
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
  barplot(height=I,names=confirmed$Day,las=2,col='skyblue',border = F,main='Total Number of Confirmed cases (Spain)')
}

data <- list(I=I,X=X,tq=tq,tf=tf,t0=t0,tX0=tX0,tmax=tmax,I0=I[t0],Iq=I[tq]) # Build a list for JAGS

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
  # Priors for precissiones (inverse of variance)
  tauI ~ dgamma(.01,.01)
  tauX ~ dgamma(0.01,0.01)
  y[t0] <- I0
}
"

model=jags.model(textConnection(modelstring), data=data,n.chains = 3,n.adapt = 5000) # burn-in=5000
update(model,n.iter=10000) # Update MCMC 
output.scir=coda.samples(model=model,variable.names=c("p","q","z","y","rmu","tauI","tauX","beta"), 
                          n.iter=10000, thin=1) # Extract MCMC data in matrix form 



mat.scir <- as.matrix(output.scir)
if(plot.flag==TRUE) {
print(colnames(mat.scir))  
}
return(list(output.scir=output.scir,data=data)) # Return MCMC samples and data in a list
}
simulation <- scir.bayesian() # Run code
#saveRDS(simulation,'output/simulation.rda') # Save to binary file


pinta(output.scir,data)
pinta(restored$output,restored$data)
pintaLinear(restored$output,restored$data)
gelm <- ploteaPanel(restored$output)
gelm <- ploteaPanel(output.video)
# par(mfrow=c(1,2))
# pinta(restored$output)
# pinta(output.video)
pico <- pintaPeak(output.video)

# dev.copy2pdf(file=sprintf('post-%d.pdf',tmax))
