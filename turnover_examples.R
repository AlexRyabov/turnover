#' Alexey Ryabov 2020
#' examples of using turnover and turnover_s

source("turnover.R")

##read and normalize data
##
data <- read.csv( file = "Species.csv")
# Date,X,Y,Species1,Species2,...

#number of rows and columns
M <- nrow(data)
N <- ncol(data)

#Define columns with species abundance data
SpecColumns = 4:N;


#turnover index richness based 
SERr = turnover_s(data[, SpecColumns]) #default parameters
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
#note that our richness based turnover index is equivalent to the binary distance in R
#so we can get the same result using 
SERr2 = as.matrix(dist(data[, SpecColumns][, ], method = "binary"))
#both functions give the same result
plot(SERr, SERr2)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")

#plot SERr vs SERs
plot(SERr, SERa)

#richness based turnover characteristics (list of turnover index + other metrics)
lSERr = turnover(data[, SpecColumns], method = "SERr") 

#abundance based turnover characteristics (list of turnover index + other metrics) 
lSERa = turnover(data[, SpecColumns], method = "SERa")

#plot effective number of extinct species vs turnover index
plot(lSERa$SER, lSERa$S_ext)

#plot effective number of extinct species vs common number of species
plot(lSERa$S_common, lSERa$S_ext)


#include observation dates
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #check if you use input data have class date
lSERa = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates)
#plot turnover as a function time intervals
plot(lSERa$TimeIntv/365, lSERa$SER)

#plot turnover index as a function of distance between observations
XY = data[, 2:3]; 
lSERa = turnover(data[, SpecColumns], method = "SERa", locations =  XY)
plot(lSERa$Dist, lSERa$SER)

