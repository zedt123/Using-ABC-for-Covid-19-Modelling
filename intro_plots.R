library(ggplot2)
library(tidyquant)
library(dplyr)
library(scales)

data <- read.csv('../Data/england_SIR_data.csv')
data <- data[c('date', 'daily_I')]
data <- data[data$date <= '2022-02-25',]
data$date <- as.Date(data$date)

tests <- read.csv('../Data/england_tests.csv')
tests <- tests[tests$date >= '2020-01-31' & tests$date <= '2022-02-25',]
tests$date <- as.Date(tests$date)



event <- c('2020-01-31', '2020-03-05', '2020-03-26', '2020-05-10',
           '2020-06-01', '2020-06-15', '2020-07-03', '2020-07-17',
           '2020-07-26', '2020-07-31', '2020-09-14', '2020-09-17',
           '2020-09-21', '2020-11-05', '2020-12-02', '2021-01-04',
           '2021-03-08', '2021-04-12', '2021-07-19', '2021-12-08')
event <- as.Date(event)

lockdowns <- c('2020-03-26', '2020-11-05', '2021-01-04')
lockdowns <- as.Date(lockdowns)

data_ma <- data %>%
  mutate(I_ma = zoo::rollmean(daily_I, k = 7, fill = 0, align = "right"))

data_ma$daily_tests <- tests$daily_tests
data_ma$daily_tests_scaled <- data_ma$daily_tests / 10

bar_fill_color <- "#6B8E23"  # Olive green for histogram bars
ma_line_color <- "#8B4513"  # Saddle brown for moving average line
point_color <- "#483D8B"   # Dark slate blue for dot points
lockdown_color <- '#0C2D48'
ma_line_size <- 1.5
point_size <- 2
lockdown_size <- 4



p <- ggplot(data_ma, aes(x=date)) +
  geom_bar(aes(y=daily_tests_scaled, color='Daily Tests'), stat = 'identity', alpha=0.3) +
  geom_bar(aes(y=daily_I, color='Daily Cases'), stat='identity', fill=bar_fill_color) +
  geom_line(aes(y=I_ma, color='Moving Average'), lty='dashed', lwd=ma_line_size) +
  geom_point(data=data[data$date %in% event,],aes(x=date, y=daily_I, color='Events'), size=point_size, shape=16) +
  geom_point(data=data[data$date %in% lockdowns,],aes(x=date, y=daily_I, color='Lockdowns'), size=lockdown_size, shape=17) +
  scale_color_manual(name = "", values = c("Daily Cases" = bar_fill_color, "Moving Average" = ma_line_color, 'Events'=point_color, 'Lockdowns'= lockdown_color,'Daily Tests' = 'grey')) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1), legend.background = element_blank()) +
  xlab("") +
  scale_y_continuous(labels = scales::scientific,name = "Daily Cases",sec.axis = sec_axis(trans=~.*10, name="Daily Tests", labels = scales::scientific))

p

### Plot for 1st analysis ###

data <- read.csv('../Data/england_SIR_data.csv')
data <- data[c('date', 'daily_I')]
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
data$date <- as.Date(data$date)

tests <- read.csv('../Data/england_tests.csv')
tests <- tests[tests$date >= '2020-03-06' & tests$date <= '2020-06-30',]
tests$date <- as.Date(tests$date)



event <- c('2020-01-31', '2020-03-05', '2020-03-26', '2020-05-10',
           '2020-06-01', '2020-06-15', '2020-07-03', '2020-07-17',
           '2020-07-26', '2020-07-31', '2020-09-14', '2020-09-17',
           '2020-09-21', '2020-11-05', '2020-12-02', '2021-01-04',
           '2021-03-08', '2021-04-12', '2021-07-19', '2021-12-08')
event <- as.Date(event)

lockdowns <- c('2020-03-26', '2020-11-05', '2021-01-04')
lockdowns <- as.Date(lockdowns)

data_ma <- data %>%
  mutate(I_ma = zoo::rollmean(daily_I, k = 7, fill = 0, align = "right"))

data_ma$daily_tests <- tests$daily_tests
data_ma$daily_tests_scaled <- data_ma$daily_tests / 10

bar_fill_color <- "#6B8E23"  # Olive green for histogram bars
ma_line_color <- "#8B4513"  # Saddle brown for moving average line
point_color <- "#483D8B"   # Dark slate blue for dot points
lockdown_color <- '#0C2D48'
ma_line_size <- 1.5
point_size <- 4
lockdown_size <- 5



p <- ggplot(data_ma, aes(x=date)) +
  geom_bar(aes(y=daily_tests_scaled, color='Daily Tests'), stat = 'identity', alpha=0.3) +
  geom_bar(aes(y=daily_I, color='Daily Cases'), stat='identity', fill=bar_fill_color) +
  geom_line(aes(y=I_ma, color='Moving Average'), lty='dashed', lwd=ma_line_size) +
  geom_point(data=data[data$date %in% event,],aes(x=date, y=daily_I, color='Events'), size=point_size, shape=16) +
  geom_point(data=data[data$date %in% lockdowns,],aes(x=date, y=daily_I, color='Lockdowns'), size=lockdown_size, shape=17) +
  scale_color_manual(name = "", values = c("Daily Cases" = bar_fill_color, "Moving Average" = ma_line_color, 'Events'=point_color, 'Lockdowns'= lockdown_color,'Daily Tests' = 'grey')) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1), legend.background = element_blank()) +
  xlab("") +
  scale_y_continuous(labels = scales::scientific,name = "Daily Cases",sec.axis = sec_axis(trans=~.*10, name="Daily Tests", labels = scales::scientific))

p


### Plot for 2st analysis ###

library(ramify)
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[c('date', 'daily_I')]
data <- data[data$date >= '2020-06-29' & data$date <= '2021-04-11',]
data$date <- as.Date(data$date)

tests <- read.csv('../Data/england_tests.csv')
tests <- tests[tests$date >= '2020-06-29' & tests$date <= '2021-04-11',]
tests$date <- as.Date(tests$date)



event <- c('2020-01-31', '2020-03-05', '2020-03-26', '2020-05-10',
           '2020-06-01', '2020-06-15', '2020-07-03', '2020-07-17',
           '2020-07-26', '2020-07-31', '2020-09-14', '2020-09-17',
           '2020-09-21', '2020-11-05', '2020-12-02', '2021-01-04',
           '2021-03-08', '2021-04-12', '2021-07-19', '2021-12-08')
event <- as.Date(event)

lockdowns <- c('2020-03-26', '2020-11-05', '2021-01-04')
lockdowns <- as.Date(lockdowns)

data_ma <- data %>%
  mutate(I_ma = zoo::rollmean(daily_I, k = 7, fill = 0, align = "right"))

data_ma$daily_tests <- tests$daily_tests
data_ma$daily_tests_scaled <- data_ma$daily_tests / 10

bar_fill_color <- "#6B8E23"  # Olive green for histogram bars
ma_line_color <- "#8B4513"  # Saddle brown for moving average line
point_color <- "#483D8B"   # Dark slate blue for dot points
lockdown_color <- '#0C2D48'
ma_line_size <- 1.5
point_size <- 4
lockdown_size <- 5

data_ma$daily_tests_scaled <- clip(data_ma$daily_tests_scaled,0,1.0e+05)

p <- ggplot(data_ma, aes(x=date)) +
  geom_bar(aes(y=daily_tests_scaled, color='Daily Tests'), stat = 'identity', alpha=0.3) +
  geom_bar(aes(y=daily_I, color='Daily Cases'), stat='identity', fill=bar_fill_color) +
  geom_line(aes(y=I_ma, color='Moving Average'), lty='dashed', lwd=ma_line_size) +
  geom_point(data=data[data$date %in% event,],aes(x=date, y=daily_I, color='Events'), size=point_size, shape=16) +
  geom_point(data=data[data$date %in% lockdowns,],aes(x=date, y=daily_I, color='Lockdowns'), size=lockdown_size, shape=17) +
  scale_color_manual(name = "", values = c("Daily Cases" = bar_fill_color, "Moving Average" = ma_line_color, 'Events'=point_color, 'Lockdowns'= lockdown_color,'Daily Tests' = 'grey')) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1), legend.background = element_blank()) +
  xlab("") +
  scale_y_continuous(labels = scales::scientific,name = "Daily Cases",sec.axis = sec_axis(trans=~.*10, name="Daily Tests", labels = scales::scientific), limits = c(0,1.0e+05))
p

