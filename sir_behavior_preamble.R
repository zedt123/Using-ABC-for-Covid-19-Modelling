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
    dS <- -(alpha_1*exp(-alpha_2*I/N)) * S * I / N
    dI <-  (alpha_1*exp(-alpha_2*I/N)) * S * I / N - gamma * I
    dR <- gamma * I
    res <- c(dS, dI, dR)
    list(res)
  })
}

run_model <- function(S0, I0, R0, alpha_1, alpha_2, gamma, t_total){
  out <- lsoda(func = SIR,
               y = c(S = S0, I = I0, R = R0),
               parms = c(alpha_1=alpha_1, alpha_2=alpha_2, gamma = gamma),
               times = seq(1, t_total, by = 1))
  return(data.frame(out[,c('S','I','R')]))
}

calc_distance <- function(D, D_star){
  return(sqrt(sum((D - D_star)**2)))
}