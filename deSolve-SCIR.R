# Sun Jul 19 15:22:52 2020 ------------------------------
# Mario Castro

require(deSolve)

# Desolve function to integrate SCIR model
SCIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS  <- -beta/N * I * S - q*S + p*C
    dC <- q*S - p*C
    dI  <- beta/N * I * S -(rmu) * I
    dX <- rmu*I
    list(c(dS,dC,dI,dX))
  })
}
I0 <- 1
N <- 45000000
init <- c(S = N-I0,  C=0, I = I0,X=0)
params <- c(beta=0.43, rmu = 0.021, p=0.007*0, q=0.062*0 ,I0 = I0)

t <- 1:11 # time in days before confinement
fit0 <- data.frame(ode(y = init, times = t, func = SCIR, parms = params))
I0 <- fit0$I[length(fit0$I)]
t <- 11:33
init <- c(S = N-I0,  C=0, I = I0,X=0)
params <- c(beta=0.43, rmu = 0.021 , p=0.007, q=0.062 ,I0 = 18)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params))

