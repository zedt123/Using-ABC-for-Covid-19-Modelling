library(deSolve)
SIR <- function(t,x, parameters){
  ##Inputs:
  #t: time
  #x: state variables
  #parameters: parameter vector
  S <- x[1]
  I <- x[2]
  R <- x[3]
  N <- S+I+R
  
  #The with statement means we can use the parameter names as in parameters
  with(as.list(parameters),{
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    res <- c(dS, dI, dR)
    list(res)
  })
}

run_model <- function(S0, I0, R0, beta, gamma, t_total){
  out <- lsoda(func = SIR,
               y = c(S = S0, I = I0, R = R0),
               parms = c(beta = beta, gamma = gamma),
               times = seq(1, t_total, by = 1))
  return(data.frame(out[,c('S','I','R')]))
}

calc_distance <- function(D, D_star){
  return(sqrt(sum((D - D_star)**2)))
}

