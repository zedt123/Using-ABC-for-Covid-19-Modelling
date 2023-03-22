source('sir_behavior_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]
# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 300 # Number of accepted particles
epsilon <- 5363.447 # Epsilon value
n_par <- 3 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('alpha_1', 'alpha_2', 'gamma', 'distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  alpha_1_0 <- runif(1, 4, 8)
  alpha_2_0 <- runif(1, 0, 1000)
  gamma_0 <- runif(1, 4, 8)
  if(alpha_2_0<0){
    next
  }
  
  
  D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_1_0, alpha_2_0, gamma_0, t_final)
  distance <- calc_distance(data$I, D_star$I)
  if (distance<epsilon){
    break
  }
}


res[1,] <- c(alpha_1_0, alpha_2_0, gamma_0, distance)
res[1,] <- c(2.586383, 860.6431, 2.248640, 58187.72)
#res[1,] <- c(6.059977, 3.142311, 5.979825, 81420.05)

while(i < N){
  # Generate a proposal
  alpha_1_star <- rnorm(1,res[i,1],0.1)
  alpha_2_star <- rnorm(1,res[i,2],1)
  gamma_star <- rnorm(1,res[i,3],0.1)
  
  #alpha_2_star <- res[[i,2]]
  
  if(alpha_1_star<0){
    next
  }
  
  if(alpha_2_star<0){
    next
  }
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_1_star, alpha_2_star, gamma_star, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)
  
  # Accept or reject the parameters
  if(distance < epsilon){
    res[i+1,] <- c(alpha_1_star, alpha_2_star, gamma_star, distance)
    i <- i + 1
    epsilon <- distance
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 5000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)
alpha_1 <- 4.306061
alpha_2 <- 88.75176
gamma <- 4.094363
D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, gamma, t_final)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,30000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

ts <- seq(1,t_final,by=1)
plot(ts,alpha_1 * exp(-alpha_2 * D_star$I/(D_star$S[1]+D_star$I[1]+D_star$R[1])),type = 'o', xlab = 'Day',ylab = 'Transmission Rate')


# Compare to constant beta models

source('case_1/case_1_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]
# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

beta <- 5.784274
gamma <- 5.622316
D_star_2 <- run_model(data$S[1], data$I[1], data$R[1], beta, gamma, t_final)


source('case_1/sir_sin_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]
# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

alpha_1 <- 1.788685
alpha_2 <- 0.1359589
alpha_3 <- 0.04717872
gamma <- 1.718247
D_star_4 <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, alpha_3, gamma, t_final)

source('sir_exp_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]
# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

alpha_1 <- 3.703882
alpha_2 <- 5.666301
alpha_3 <- 0.8975742   
gamma <- 3.613029

D_star_3 <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, alpha_3, gamma, t_final)

points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     sir_sin = D_star_4$I,
                     sir = D_star_2$I,
                     sir_exp = D_star_3$I,
                     sir_behavior = D_star$I)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = sir_behavior, color = "SIR Behavior Model", shape = "SIR Behavior Model"), size = 3) +
  geom_point(aes(y = sir_exp, color = "SIR Exponential Model", shape = "SIR Exponential Model"), size = 4) +
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3) +
  scale_color_manual(values=c("blue",'red',"yellow",'green')) +
  scale_shape_manual(values=c(16, 15, 20, 17)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 

