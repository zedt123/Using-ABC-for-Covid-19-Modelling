source('sir_piecewise_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

t_final <- 500

alpha_11 <- 0.09411848
alpha_2 <- 121.0366
gamma <- 0.04108736
D_star <- run_model(data$S[1], data$I[1], data$R[1], alpha_11, alpha_2, gamma, t_final)

plot(seq(1,t_final,by=1), D_star$I, pch=15, ylim = c(0,400000), xlab = 'Day', ylab = 'Number of Cases', col = 'blue')
points(seq(1,length(data$I),by=1), data$I, pch=16, col = 'black')
legend("topright", legend=c("No Lockdown Data", "Real Data"),
       col=c("blue", "black"), pch=c(15,16), cex=0.8)

ts <- seq(1,t_final,by=1)
plot(ts,alpha_11 * exp(-alpha_2 * D_star$I/(D_star$S[1]+D_star$I[1]+D_star$R[1])),type = 'o', xlab = 'Day',ylab = 'Transmission Rate')
