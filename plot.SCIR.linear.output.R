plot.SCIR.linear.output <- function(output,data) {
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
  plot(ty[t0:tmax],exp(data$I[t0:tmax]),ylim=c(0,2.6e5),xlim=c(1,tf),
       xlab='Days since first confirmed case',
       ylab='Active cases',
       cex.axis=1.75,cex.lab=2,col='black'); 
  
  
  # Auxiliary variables to create a polygon for the 95% posterior density
  tt <- c(ty,reverse(ty)) # Patch times and reversed times
  yM <- colq(y,0.975) # 97.5% quantile (and convert to log10 scale)
  ym <- colq(y,0.0225)# 2.5% quantile
  yy <- exp(c(yM,reverse(ym))) # patch 97.5 and (reversed) 2.5 quantiles
  ix <- seq(1,180,6)
  polygon(tt[ix],yy[ix],col=rgb(0,0,0,.2),border='gray')
  
  param <- data.frame(t(colMed(mat.output[,1:4]))) # Median parameters from posteriors ("beta","p","q","rmu")
  # Create a vector with the exact solution of the SCIR model for the median parameters
  exact<- exact.SCIR(param,data)
  t.median <- exact$tout
  y.median <- exp(exact$out)
  lines(t.median,y.median,col='darkorange',lwd=4)
  points(ty[tmax:tf],exp(data$I[tmax:tf]),pch=19,cex=.8,col=2)
  
  # Annotate the plot
  arrows(tq,50000,tq,20,angle=10,lwd=2) 
  text(tq,62000,"1st Confiment",pos=3,cex=1.5)
  text(tq,50000,"begins",pos=3,cex=1.5)
  arrows(33,150000,33,80000,angle=10,lwd=2) 
  text(33,150000,"2nd Confiment begins",pos=3,cex=1.5)
  arrows(47,50000,47,95000,angle=10,lwd=2) 
  text(47,50000,"Epidemic peak",pos=1,cex=1.5)
}

for(i in 1:4) {
  temp <- simulations[[i]]
plot.SCIR.output(temp$output.scir,temp$data)
  dev.copy2pdf(file=sprintf('output/bayesian-SCIR-fit-%d.pdf',temp$data$tmax))
}

for(i in 1:4) {
  temp <- simulations[[i]]
  plot.SCIR.linear.output(temp$output.scir,temp$data)
  dev.copy2pdf(file=sprintf('output/bayesian-SCIR-linear-fit-%d.pdf',temp$data$tmax))
}
