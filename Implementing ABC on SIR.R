#Basic ABC algorithm
#1.Simulate θ from prior
#2.Use θ to simulate values from the likelihood
#3.Accept θ if it is close enough to a value from our dataset
#Note: As we have large data, we can use the mean and sd of the data to decide on whether we accept or reject θ

#In this example we consider b=0.01 recovery=0.005 to be the real parameters 
#We simulate 10000 data from these parameters

library('deSolve')

epi203 <- function(pars){
  
  times <- seq(from = 0, to = 600, by = 1) #time steps
  yinit <- c(Susc = 0.9, Infected = 0.1, Recovered = 0) #initial conditions
  
  SIR_model <- function(times, yinit, pars){
    with(as.list(c(yinit,pars)), { #with needs to take a list, data frame or integer as input
      dSusc <- - beta*Infected*Susc
      dInfected <- beta*Infected*Susc - recovery*Infected
      dRecovered <- recovery*Infected
      
      return(list(c(dSusc, dInfected, dRecovered))) #ode function needs func argument to return a list
    })
  }
  
  results <-  ode(func = SIR_model, y = yinit, times = times, parms = pars) #returns a matrix
  results <-  as.data.frame(results) #turn the matrix to a data frame
}

pars <- c(beta = 0.1, recovery = 0.005)
results <- epi203(pars)

matplot(results[,1], results[,2:4], type = 'l', lty = 1, lwd = 4, xlab = 'Time', ylab = 'Percentage')
legend("topright", legend = c('S', 'I', 'R'), lwd = 1, col = 1:3)









































