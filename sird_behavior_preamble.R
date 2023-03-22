library(deSolve)
SIRD <- function(t,x, parameters){
  ##Inputs:
  #t: time
  #x: state variables
  #parameters: parameter vector
  S <- x[1]
  I <- x[2]
  R <- x[3]
  Dw <- x[4]
  D <- x[5]
  N <- S+I+R+Dw+D
  

  
  #The with statement means we can use the parameter names as in parameters
  with(as.list(parameters),{
    
    dS <- -alpha_1 * exp(-alpha_2 * Dw / N) * S * I / N
    dI <- alpha_1 * exp(-alpha_2 * Dw / N) * S * I / N - gamma * I - zeta_1 * I
    dR <- gamma * I
    dDw <- zeta_1 * I - zeta_2 * Dw
    dD <- zeta_2 * Dw
    
    res <- c(dS, dI, dR, dDw, dD)
    list(res)
  })
}

run_model <- function(S0, I0, R0, Dw0, D0, alpha_1, alpha_2, gamma, zeta_1, zeta_2, t_total){
  out <- lsoda(func = SIRD,
               y = c(S = S0, I = I0, R = R0, Dw = Dw0, D = D0),
               parms = c(alpha_1 = alpha_1, alpha_2 = alpha_2, gamma = gamma, zeta_1 = zeta_1, zeta_2 = zeta_2),
               times = seq(1, t_total, by = 1))
  return(data.frame(out[,c('S','I','R', 'Dw','D')]))
}



calc_distance <- function(D, D_star){
  return(sqrt(sum((D - D_star)**2)))
}

