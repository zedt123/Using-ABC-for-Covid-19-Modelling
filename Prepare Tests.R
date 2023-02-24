# Prepare the data with the number of COVID-19 tests per day

tests <- read.csv('../Data/nation_2022-11-03.csv')
tests$date <- as.Date(tests$date)
tests <- tests[order(as.Date(tests$date, format="%Y/%m/%d")),]
names(tests)[names(tests) == "newVirusTestsBySpecimenDate"] <- "daily_tests"
tests <- tests[c('areaName','date', 'daily_tests')]

write.csv(tests,"../Data/england_tests.csv")
