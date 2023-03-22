source('seir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIR_data.csv')

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','E','I','R')]

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 100 # Number of accepted particles
epsilon <- 575396.5 # Epsilon value
n_par <- 4 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('beta_1','beta_2','delta', 'gamma', 'distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  beta_1_0 <- runif(1, 0, 10)
  beta_2_0 <- runif(1, 0, 10)
  delta_0 <- runif(1, 0, 10)
  gamma_0 <- runif(1, 0, 10)
  D_star <- run_model(data$S[1], data$E[1], data$I[1], data$R[1], beta_1_0, beta_2_0, delta_0, gamma_0, t_final)
  distance <- calc_distance(data$I, D_star$I)
  if (distance<epsilon){
    break
  }
}

res[1,] <- c(0.016110018, 1.0942616, 0.3299727, 0.9662622, 575396.5)
#res[1,] <- c(beta_1_0, beta_2_0, delta_0, gamma_0, distance)

while(i < N){
  # Generate a proposal
  beta_1_star <- rnorm(1,res[i,1],0.1)
  beta_2_star <- rnorm(1,res[i,2],0.1)
  delta_star <- rnorm(1,res[i,3],0.1)
  gamma_star <- rnorm(1,res[i,4],0.1)
  
  if (beta_1_star < 0){
    next
  }
  
  if (beta_2_star < 0){
    next
  }
  
  if (delta_star < 0){
    next
  } 
  
  if (gamma_star < 0){
    next
  }  
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$E[1], data$I[1], data$R[1], beta_1_star, beta_2_star, delta_star, gamma_star, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)
  
  # Accept or reject the parameters
  if(distance <= epsilon){
    res[i+1,] <- c(beta_1_star, beta_2_star, delta_star, gamma_star, distance)
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 10000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)
beta_1 <- 0.2697847640
beta_2 <- 0.3677300
delta <- 0.4870130
gamma <- 0.6515517

D_star <- run_model(data$S[1], data$E[1], data$I[1], data$R[1], beta_1, beta_2, delta, gamma, t_final)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,300000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

write.csv(res, 'England_seir_mcmc_results_analysis_2.csv')

# SIR Model to compare results

source('case_1/case_1_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

beta <- 0.5086747
gamma <- 0.4683127
D_star_2 <- run_model(data$S[1], data$I[1], data$R[1], beta, gamma, t_final)

points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     seir = D_star$I,
                     sir = D_star_2$I
)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = seir, color = "SEIR Model", shape = "SEIR Model"), size = 3) +
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3) +
  scale_color_manual(values=c("blue","red","green")) +
  scale_shape_manual(values=c(16, 15, 17)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 


