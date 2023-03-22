library(deSolve)
SEIIR <- function(t,x, parameters){
  ##Inputs:
  #t: time
  #x: state variables
  #parameters: parameter vector
  S <- x[1]
  E <- x[2]
  Is <- x[3]
  Ia <- x[4]
  R <- x[5]
  N <- S+E+Is+Ia+R
  
  #The with statement means we can use the parameter names as in parameters
  with(as.list(parameters),{
    dS <- - S/N * (beta_e * E + beta_s * Is + beta_a * Ia)
    dE <- S/N * (beta_e * E + beta_s * Is + beta_a * Ia) - delta * E
    dIs <- f * delta * E - gamma_s * Is
    dIa <- (1 - f) * delta * E - gamma_a * Ia
    dR <- gamma_s * Is + gamma_a * Ia
    res <- c(dS, dE, dIs, dIa, dR)
    list(res)
  })
}

run_model <- function(S0, E0, Is0, Ia0, R0, beta_e, beta_s, beta_a, delta, gamma_s, gamma_a, f, t_total){
  out <- lsoda(func = SEIIR,
               y = c(S = S0, E = E0, Is = Is0, Ia = Ia0, R = R0),
               parms = c(beta_e = beta_e, beta_s = beta_s, beta_a = beta_a, delta = delta, gamma_s = gamma_s, gamma_a = gamma_a, f = f),
               times = seq(1, t_total, by = 1))
  return(data.frame(out[,c('S','E','Is','Ia','R')]))
}

calc_distance <- function(D, D_star){
  return(sqrt(sum((D - D_star)**2)))
}