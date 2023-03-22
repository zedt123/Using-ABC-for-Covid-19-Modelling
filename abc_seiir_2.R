source('seiir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIIR_data.csv')

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','E','Ia','Is','R')]

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 100 # Number of accepted particles
epsilon <- 604081.3 # Epsilon value
n_par <- 6 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('beta_e','beta_s','beta_a','delta','gamma_s','gamma_a','distance')
f <- 0.4 # https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2787098

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  beta_e_0 <- runif(1, 0, 10)
  beta_s_0 <- runif(1, 0, 1)
  beta_a_0 <- runif(1, 0, 10)
  delta_0 <- runif(1, 0, 10)
  gamma_s_0 <- runif(1, 0, 10)
  gamma_a_0 <- runif(1, 0, 10)
  
  D_star <- run_model(data$S[1], data$E[1], data$Is[1], data$Ia[1], data$R[1], beta_e_0, beta_s_0, beta_a_0, delta_0, gamma_s_0, gamma_a_0, f, t_final)
  distance <- calc_distance(data$Is, D_star$Is) + calc_distance(data$Ia,D_star$Ia) # I is the sum of Is and Ia
  if (distance<epsilon){
    break
  }
}

res[1,] <- c(2.369735, 0.3078401, 0.35549869, 9.525765, 0.2657443, 0.5857766, 604081.3)
#res[1,] <- c(beta_e_0, beta_s_0, beta_a_0, delta_0, gamma_s_0, gamma_a_0, distance)

while(i < N){
  # Generate a proposal
  beta_e_star <- rnorm(1,res[i,1],0.01)
  beta_s_star <- rnorm(1,res[i,2],0.01)
  beta_a_star <- rnorm(1,res[i,3],0.01)
  delta_star <- rnorm(1,res[i,4],0.01)
  gamma_s_star <- rnorm(1,res[i,5],0.01)
  gamma_a_star <- rnorm(1,res[i,6],0.01)
  
  if(beta_e_star<0){
    next
  }
  
  if(beta_s_star<0){
    next
  }
  
  if(beta_a_star<0){
    next
  }
  
  if(delta_star<0){
    next
  }  

  if(gamma_s_star<0){
    next
  } 
  
  if(gamma_a_star<0){
    next
  }   
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$E[1], data$Is[1], data$Ia[1],  data$R[1], beta_e_star, beta_s_star, beta_a_star, delta_star, gamma_s_star, gamma_a_star, f, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$Is, D_star$Is) + calc_distance(data$Ia,D_star$Ia)
  
  # Accept or reject the parameters
  if(distance <= epsilon){
    res[i+1,] <- c(beta_e_star, beta_s_star, beta_a_star, delta_star, gamma_s_star, gamma_a_star, distance)
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 10000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)

beta_e <- 0.0004211976
beta_s <- 12.51536
beta_a <- 5.705052
delta <- 0.03770838
gamma_s <- 2.329986
gamma_a <- 4.849722

beta_e <- 2.364649
beta_s <- 0.3017949
beta_a <- 0.3557738
delta <- 9.516286
gamma_s <- 0.2618110
gamma_a <- 0.5801952

D_star <- run_model(data$S[1], data$E[1], data$Is[1], data$Ia[1], data$R[1], beta_e, beta_s, beta_a, delta, gamma_s, gamma_a, f, t_final)

calc_distance(data$Is+data$Ia, D_star$Is+D_star$Ia)

plot(seq(1,t_final,by=1), data$Is+data$Ia, pch=16, ylim = c(0,300000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$Is+D_star$Ia, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

write.csv(res, 'England_seiir_mcmc_results_4.csv')

# SIR, SEIR Models to compare results

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

source('seir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIR_data.csv')


data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','E','I','R')]

beta_1 <- 0.2697847640
beta_2 <- 0.3677300
delta <- 0.4870130
gamma <- 0.6515517

D_star_3 <- run_model(data$S[1], data$E[1], data$I[1], data$R[1], beta_1, beta_2, delta, gamma, t_final)


points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     seiir = D_star$Is + D_star$Ia,
                     seir = D_star_3$I,
                     sir = D_star_2$I
)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = seiir, color = "SEIIR Model", shape = "SEIIR Model"), size = 4) +
  geom_point(aes(y = seir, color = "SEIR Model", shape = "SEIR Model"), size = 3) +
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3) +
  scale_color_manual(values=c("blue",'purple',"red","green")) +
  scale_shape_manual(values=c(16, 18, 15, 17)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 

