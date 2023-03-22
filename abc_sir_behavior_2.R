source('sir_behavior_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 300 # Number of accepted particles
epsilon <- 834392.3 # Epsilon value
n_par <- 3 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('alpha_1', 'alpha_2', 'gamma', 'distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  alpha_1_0 <- runif(1, 0, 10)
  alpha_2_0 <- runif(1, 0, 1000)
  gamma_0 <- runif(1, 0, 10)
  
  if(alpha_1_0<0){
    next
  }
  
  D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_1_0, alpha_2_0, gamma_0, t_final)
  distance <- calc_distance(data$I, D_star$I)
  j <- j+1
  if (distance<epsilon){
    res[1,] <- c(alpha_1_0, alpha_2_0, gamma_0, distance)
    epsilon <- distance
    print(j)
    print(res[1,])
  }
}


res[1,] <- c(alpha_1_0, alpha_2_0, gamma_0, distance)
res[1,] <- c(0.3470168, 119.4436, 0.2681696, 834392.3)

while(i < N){
  # Generate a proposal
  alpha_1_star <- rnorm(1,res[i,1],0.01)
  alpha_2_star <- rnorm(1,res[i,2],1)
  gamma_star <- rnorm(1,res[i,3],0.01)
  
  #alpha_1_star <- res[[i,1]]
  #alpha_2_star <- res[[i,2]]
  #gamma_star <- res[[i,4]]
  
  if(gamma_star<0){
    next
  }
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_1_star, alpha_2_star, gamma_star, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)
  
  # Accept or reject the parameters
  if(distance < epsilon){
    res[i+1,] <- c(alpha_1_star, alpha_2_star, gamma_star, distance)
    epsilon <- distance
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 5000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)

alpha_1 <- 0.6841056
alpha_2 <- -20.16637
gamma <- 0.6513138

D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, gamma, t_final)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,300000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

ts <- seq(1,t_final,by=1)
plot(ts,alpha_1 * exp(-alpha_2 * D_star$I/(D_star$S[1]+D_star$I[1]+D_star$R[1])),type = 'o', xlab = 'Day',ylab = 'Transmission Rate')


# Compare to constant beta models

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


source('smeiir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIIR_data.csv')

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','E','Ia','Is','R')]

beta_e <- 0.04573545
beta_s <- 1.327204
beta_a <- 0.5656386
delta <- 0.3133744
gamma_s <- 0.8620612
gamma_a <- 1.950651

beta_Me <- beta_e/2
beta_Ms <- beta_s/2
beta_Ma <- beta_a/2

f <- 0.4

D_star_5 <- run_model(0.85*data$S[1], 0.15*data$S[1], data$E[1], data$Is[1], data$Ia[1], data$R[1], beta_e, beta_s, beta_a, beta_Me, beta_Ms, beta_Ma, delta, gamma_s, gamma_a, f, t_final)


points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$Is + data$Ia,
                     sir_cos = D_star$I,
                     smeiir = D_star_5$Is + D_star_5$Ia,
                     sir = D_star_2$I)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = sir_cos, color = "SIR Cosine Model", shape = "SIR Cosine Model"), size = 3) +
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3) +
  geom_point(aes(y = smeiir, color = "SMEIIR Model", shape = "SMEIIR Model"), size = 3) +
  scale_color_manual(values=c("blue","purple",'green', "orange")) +
  scale_shape_manual(values=c(16, 18, 17, 19)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 


