library(adaptMCMC)
library(deSolve)
library(arm)


#compartments
sir <- function(time, state, parameters) { ## this is the ODE system
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

#initial vals and actual parameters
init <- c(S = 100000-10, I = 10, R = 0.0)
true.beta <- 0.53
true.gamma <- 1/6
parameters <- c(beta = true.beta, gamma = true.gamma)

#time point seq
times <- seq(0, 77, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
out$time <- NULL #out is the array of values from running model, time not needed so discarded

true.sigma <- .006 #the SD
y_obs <- rnorm(n = nrow(out), mean = out$I, sd = true.sigma) # data

matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovered", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(40, 0.7, c("Susceptibles", "Infected", "Recovered"), pch = 1, col = 2:4) #SIr trajectories

#########
### BEGIN MACHINERY
# par is an array containing (I[0], R[0], S[0], beta, gamma, sigma)
# In unconstrained space, we have transf(par) = c(log(beta), log(gamma), logit(I[0]), logit(R[0]), logit(S[0]), log(sigma))
# 'sigma' here is the standard deviation of the normal likelihood
# Note: here, I chose to include the initial conditions as part of the parameters to be estimated
getTransformedParameter <- function(par){
  ## this function takes the parameter on an unconstrained scale and returns parameters on their natural scale
  transfpar <- rep(NA, length(par))
  unnorm <- arm::invlogit(par[1:3])
  transfpar[1:3] <- unnorm/sum(unnorm)  ## transform initial conditions
  transfpar[4:6] <- exp(par[4:6])
  names(transfpar) <- c("S", "I", "R", "beta", "gamma", "sigma")
  return(transfpar)
}
#
getSolution <- function(par, times){
  ## this function takes the parameters and times and returns the solution to the ODEs
  allParameters <- getTransformedParameter(par)
  sol <- as.data.frame(ode(y = allParameters[1:3], times = times, func = sir, parms = allParameters[4:5]))
  return(sol$I)
}
#
Likelihood <- function(par, times, data){
  sol <- getSolution(par, times)
  return( 
    sum(dnorm(data, mean = sol, sd = exp(par[6]), log = TRUE)) 
  ) 
}
#
Prior <- function(par){
  tpars <- getTransformedParameter(par)
  lpr <- dbeta(tpars[1], 1, 1, log = TRUE) ## pi(S0)
  lpr <- lpr + dbeta(tpars[2], 1, 1, log = TRUE) ## pi(I0)
  lpr <- lpr + dbeta(tpars[3], 1, 1, log = TRUE) ## pi(R0)
  lpr <- lpr + dgamma(tpars[4], 1, 1, log = TRUE) ## pi(beta)
  lpr <- lpr + dgamma(tpars[5], 1, 1, log = TRUE) ## pi(gamma)
  lpr <- lpr + dgamma(tpars[6], .1, .1, log = TRUE) ## pi(sigma)
  return(lpr)
}
#
Target <- function(pars, times, data){
  return(
    Likelihood(par = pars, times = times, data = data) + Prior(pars)
  )
}
### END MACHINERY
#############

### Running the algorithm

init.par <- rnorm(6) ## initial guess for the parameters
chain <- MCMC(p = Target, init = init.par,  adapt = TRUE, acc.rate = .234,  n = 5E4, times = times, data = y_obs)

### Annotating results

Samples <- chain$samples
for(i in 1:nrow(Samples)){
  Samples[i, ] <- getTransformedParameter(chain$samples[i, ])
}
burnin <- .2
Samples.bnin <- Samples[round(burnin * nrow(Samples)):nrow(Samples), ]

hist(Samples.bnin[, 4], xlab = expression(beta), main = "Posterior of infection rate")
abline(v = true.beta, lwd = 2, lty = 2)

hist(Samples.bnin[, 5], xlab = expression(gamma), main = "Posterior of recovery rate")
abline(v = true.gamma, lwd = 2, lty = 2)

hist(Samples.bnin[, 6], xlab = expression(sigma), main = "Posterior of likelihood standard deviation")
abline(v = true.sigma, lwd = 2, lty = 2)


