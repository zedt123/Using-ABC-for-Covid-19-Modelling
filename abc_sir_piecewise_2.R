source('sir_piecewise_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

t_final_1 <- 74
t_final_2 <- 101
t_final_3 <- 134
t_final_4 <- length(data$S)

#### ABC Setup ####
N <- 500 # Number of accepted particles
epsilon <- 147393.4 # Epsilon value
n_par <- 6 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('alpha_11', 'alpha_12','alpha_13','alpha_14','alpha_2', 'gamma', 'distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  alpha_11_0 <- runif(1, 0, 10)
  alpha_12_0 <- runif(1, 0, 10)
  alpha_13_0 <- runif(1, 0, 10)
  alpha_14_0 <- runif(1, 0, 10)
  alpha_2_0 <- runif(1, 0, 1000)
  gamma_0 <- runif(1, 0, 10)
  
  D_star_1 <- run_model(data$S[1], data$I[1], data$R[1], alpha_11_0, alpha_2_0, gamma_0, t_final_1)
  D_star_2 <- run_model(D_star_1$S[t_final_1], D_star_1$I[t_final_1], D_star_1$R[t_final_1], alpha_12_0, alpha_2_0, gamma_0, t_final_2-t_final_1+1)
  D_star_3 <- run_model(D_star_2$S[t_final_2-t_final_1+1], D_star_2$I[t_final_2-t_final_1+1], D_star_2$R[t_final_2-t_final_1+1], alpha_13_0, alpha_2_0, gamma_0, t_final_3-t_final_2+1)
  D_star_4 <- run_model(D_star_3$S[t_final_3-t_final_2+1], D_star_3$I[t_final_3-t_final_2+1], D_star_3$R[t_final_3-t_final_2+1], alpha_14_0, alpha_2_0, gamma_0, t_final_4-t_final_3+1)
  D_star <- data.frame(S = c(D_star_1$S, D_star_2$S[2:(t_final_2-t_final_1+1)], D_star_3$S[2:(t_final_3-t_final_2+1)], D_star_4$S[2:(t_final_4-t_final_3+1)]),
                       I = c(D_star_1$I, D_star_2$I[2:(t_final_2-t_final_1+1)], D_star_3$I[2:(t_final_3-t_final_2+1)], D_star_4$I[2:(t_final_4-t_final_3+1)]),
                       R = c(D_star_1$R, D_star_2$R[2:(t_final_2-t_final_1+1)], D_star_3$R[2:(t_final_3-t_final_2+1)], D_star_4$R[2:(t_final_4-t_final_3+1)]))
  
  distance <- calc_distance(data$I, D_star$I)
  
  if(is.na(distance)){
    next
  }
  
  j <- j+1
  if (distance<epsilon){
    res[1,] <- c(alpha_11_0, alpha_12_0, alpha_13_0, alpha_14_0,alpha_2_0, gamma_0, distance)
    epsilon <- distance
    print(j)
    print(res[1,])
  }
}

#res[1,] <- c(alpha_11_0, alpha_12_0, alpha_13_0, alpha_2_0, gamma_0, distance)
res[1,] <- c(0.09411848, 0.006326168, 0.1292559, 0.002343643, 121.0366, 0.04108736, 147393.4)

while(i < N){
  # Generate a proposal
  alpha_11_star <- rnorm(1,res[i,1],0.0001)
  alpha_12_star <- rnorm(1,res[i,2],0.0001)
  alpha_13_star <- rnorm(1,res[i,3],0.0001)
  alpha_14_star <- rnorm(1,res[i,4],0.0001)
  alpha_2_star <- rnorm(1,res[i,5],0.001)
  gamma_star <- rnorm(1,res[i,6],0.0001)
  
  #alpha_1_star <- res[[i,1]]
  #alpha_2_star <- res[[i,2]]
  #gamma_star <- res[[i,4]]
  
  if(alpha_11_star<0){
    next
  }
  
  if(alpha_12_star<0){
    next
  }
  
  if(alpha_13_star<0){
    next
  }
  
  if(gamma_star<0){
    next
  }
  
  # Simulate dataset from the model
  D_star_1 <- run_model(data$S[1], data$I[1], data$R[1], alpha_11_star, alpha_2_star, gamma_star, t_final_1)
  D_star_2 <- run_model(D_star_1$S[t_final_1], D_star_1$I[t_final_1], D_star_1$R[t_final_1], alpha_12_star, alpha_2_star, gamma_star, t_final_2-t_final_1+1)
  D_star_3 <- run_model(D_star_2$S[t_final_2-t_final_1+1], D_star_2$I[t_final_2-t_final_1+1], D_star_2$R[t_final_2-t_final_1+1], alpha_13_star, alpha_2_star, gamma_star, t_final_3-t_final_2+1)
  D_star_4 <- run_model(D_star_3$S[t_final_3-t_final_2+1], D_star_3$I[t_final_3-t_final_2+1], D_star_3$R[t_final_3-t_final_2+1], alpha_14_star, alpha_2_star, gamma_star, t_final_4-t_final_3+1)
  D_star <- data.frame(S = c(D_star_1$S, D_star_2$S[2:(t_final_2-t_final_1+1)], D_star_3$S[2:(t_final_3-t_final_2+1)], D_star_4$S[2:(t_final_4-t_final_3+1)]),
                       I = c(D_star_1$I, D_star_2$I[2:(t_final_2-t_final_1+1)], D_star_3$I[2:(t_final_3-t_final_2+1)], D_star_4$I[2:(t_final_4-t_final_3+1)]),
                       R = c(D_star_1$R, D_star_2$R[2:(t_final_2-t_final_1+1)], D_star_3$R[2:(t_final_3-t_final_2+1)], D_star_4$R[2:(t_final_4-t_final_3+1)]))
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)
  
  # Accept or reject the parameters
  if(distance < epsilon){
    res[i+1,] <- c(alpha_11_star, alpha_12_star, alpha_13_star, alpha_14_star, alpha_2_star, gamma_star, distance)
    epsilon <- distance
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 5000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}


alpha_11 <- 0.09411848
alpha_12 <- 0.006326168
alpha_13 <- 0.1292559
alpha_14 <- 0.002343643
alpha_2 <- 121.0366
gamma <- 0.04108736

D_star_1 <- run_model(data$S[1], data$I[1], data$R[1], alpha_11, alpha_2, gamma, t_final_1)
D_star_2 <- run_model(D_star_1$S[t_final_1], D_star_1$I[t_final_1], D_star_1$R[t_final_1], alpha_12, alpha_2, gamma, t_final_2-t_final_1+1)
D_star_3 <- run_model(D_star_2$S[t_final_2-t_final_1+1], D_star_2$I[t_final_2-t_final_1+1], D_star_2$R[t_final_2-t_final_1+1], alpha_13, alpha_2, gamma, t_final_3-t_final_2+1)
D_star_4 <- run_model(D_star_3$S[t_final_3-t_final_2+1], D_star_3$I[t_final_3-t_final_2+1], D_star_3$R[t_final_3-t_final_2+1], alpha_14, alpha_2, gamma, t_final_4-t_final_3+1)
D_star <- data.frame(S = c(D_star_1$S, D_star_2$S[2:(t_final_2-t_final_1+1)], D_star_3$S[2:(t_final_3-t_final_2+1)], D_star_4$S[2:(t_final_4-t_final_3+1)]),
                     I = c(D_star_1$I, D_star_2$I[2:(t_final_2-t_final_1+1)], D_star_3$I[2:(t_final_3-t_final_2+1)], D_star_4$I[2:(t_final_4-t_final_3+1)]),
                     R = c(D_star_1$R, D_star_2$R[2:(t_final_2-t_final_1+1)], D_star_3$R[2:(t_final_3-t_final_2+1)], D_star_4$R[2:(t_final_4-t_final_3+1)]))

plot(seq(1,t_final_4,by=1), data$I, pch=16, ylim = c(0,300000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final_4,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)


ts <- seq(1,t_final_4,by=1)
beta <- c(alpha_11 * exp(-alpha_2 * D_star$I[1:t_final_1]/(D_star$S[1]+D_star$I[1]+D_star$R[1])),
          alpha_12 * exp(-alpha_2 * D_star$I[(t_final_1+1):t_final_2]/(D_star$S[1]+D_star$I[1]+D_star$R[1])),
          alpha_13 * exp(-alpha_2 * D_star$I[(t_final_2+1):t_final_3]/(D_star$S[1]+D_star$I[1]+D_star$R[1])),
          alpha_14 * exp(-alpha_2 * D_star$I[(t_final_3+1):t_final_4]/(D_star$S[1]+D_star$I[1]+D_star$R[1]))
)
plot(ts,beta,type = 'o', xlab = 'Day',ylab = 'Transmission Rate')


# Comparisons
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

source('case_1/sir_sin_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

alpha_1 <- 0.6336706
alpha_2 <- 0.03869089
alpha_3 <- 0.05035212  
gamma <- 0.5795639

D_star_3 <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, alpha_3, gamma, t_final)


source('sir_behavior_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

alpha_1 <- 0.6841056
alpha_2 <- -20.16637
gamma <- 0.6513138

D_star_4 <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, gamma, t_final)

points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     sir_piece = D_star$I,
                     sir = D_star_2$I,
                     sir_cos = D_star_3$I,
                     sir_behave = D_star_4$I)


library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +  
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 2.5) +
  geom_point(aes(y = sir_behave, color = "SIR Behavioral  Model", shape = "SIR Behavioral  Model"), size = 3) +
  geom_point(aes(y = sir_cos, color = "SIR Cosine  Model", shape = "SIR Cosine  Model"), size = 3) +
  geom_point(aes(y = sir_piece, color = "SIR Piecewise  Model", shape = "SIR Piecewise  Model"), size = 4) +
  scale_color_manual(values=c("blue",'grey','purple',"green","black")) +
  scale_shape_manual(values=c(16, 19 ,18, 17, 20)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 

