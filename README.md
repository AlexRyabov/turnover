# Species turnover index
R functions for calculating the species turnover, based on Hillebrand et al. (2018). *Biodiversity change is uncoupled from species richness trends: Consequences for conservation and monitoring.* J Appl Ecol, 55, 169–184. https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12959

Input is an (M x N) table or matrix with M observations of N species and (optionally) observation dates and locations.

There are two main function  `turnover` and `turnover_s`

`SER = turnover_s(X, method)` is a short version which returns an upper triangle matrix with turnover indexes between different rows. SER[i, j] is turnover between X[i, ] and X[j, ]. 
method = “SERr” for richness based turnover and “SERa” for abundance based turnover

`lSER = turnover(X, method)` returns a list of upper triangle matrices with the turnover indexes and  effective numbers of common/extinct/immigrating/total species. Again element i,j of each matrix describes changes in species composition from observation X[i, ] to X[j, ].

`lSER = turnover(X, method, dates, locations)` additionally returns upper triangle matrices of temporal and spatial distances between observations i and j.



```R
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
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #convert the input dates from string into class "date"
lSERa = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates)
#plot turnover as a function time intervals
plot(lSERa$TimeIntv/365, lSERa$SER)

#plot turnover index as a function of distance between observations
XY = data[, 2:3]; 
lSERa = turnover(data[, SpecColumns], method = "SERa", locations =  XY)
plot(lSERa$Dist, lSERa$SER)

```
