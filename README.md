# Species turnover index
R functions for calculating the species turnover, based on Hillebrand et al. (2018). *Biodiversity change is uncoupled from species richness trends: Consequences for conservation and monitoring.* J Appl Ecol, 55, 169–184. https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12959

Input is an (M x N) table or matrix with M observations of N species and (optionally) observation dates and locations.

There are three main function  `turnover`, `turnover_s` and `turnover_g`

`SER = turnover_s(X, method)` is a short version which returns an upper triangle matrix with turnover indexes between different rows. SER[i, j] is turnover between X[i, ] and X[j, ]. 
method = “SERr” for richness based turnover and method =“SERa” for abundance based turnover.

`lSER = turnover(X, method)` returns a list of upper triangle matrices with the turnover indexes and  effective numbers of common/extinct/immigrating/total species. Again element i,j of each matrix describes changes in species composition from observation X[i, ] to X[j, ].

`lSER = turnover(X, method, dates, locations)` additionally returns upper triangle matrices of temporal and spatial distances between observations i and j.

 The function result contains:
 
 `lSER$SER` is an upper triangle M x M matrix with turnover index. `lSER$SER[i,j]` is turnover between row i and j in X
 
 `lSER$S_common, lSER$S_ext, lSER$S_imm, lSER$S_total` are upper triangle MxM matrices  
 with the effective number of common/extinct/immigrating/total species,
 e.g., `SER$S_total[i,j]` is the total effective number of species in observation i and j
 and  `lSER$S_ext[i,j]` is the effective number of species present in X[i, ] but not in X[j, ]
 
 Note that the effective species numbers are richness based when  `method="SERr"` and 
 abundance (Simpson index) based when `method="SERa"`, 
 see details in Hillebrand, H. et al. J Appl Ecol 55, 169–184 (2018).
 
 `lSER$TimeIntv[i,j]` time intervals between observation i and j
 
 `lSER$Dist[i,j]`  spatial distances between observation i and j
 

`lSER_gr = turnover_g(X, method = "SERa", dates = NULL, locations = NULL, groupby)`
This function calculates turnover within groups defined by `groupby` dataframe. 
`groupby` must contain the same number of rows as `X`. The function splits `X` 
into groups and returns a list with two fields: `lSER_gr$groupnames` contains 
a list unique group names, `lSER_gr$turnover` contains a list with information 
on species turnover within each group. 
 



```R
#' Alexey Ryabov 2020
#' examples of using turnover, turnover_s and turnover_g

source("turnover.R")

##read data
##
data = read.csv( file = "Species.csv")
# Date,X,Y,Species1,Species2,...

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 4:N;


#richness based turnover index  
SERr = turnover_s(data[, SpecColumns]) #default parameters
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
#note that our richness based turnover index is equivalent to the binary distance in R
#so you can get the same result using 
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


#Group by some factors 
#define stations and areas 
nr = nrow(data)
Area = sample(c("Area 1", "Area 2"), size = nr, replace = TRUE)
Area = sort(Area)
Station = sample(c("A", "B", "C"), size = nr, replace = TRUE)
#
StatArea = data.frame(Area, Station);
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #convert the input dates from string into class "date"
lSERa_gr = turnover_g(data[, SpecColumns], method = "SERa", dates =  SampleDates, groupby = StatArea)

#lSERa_gr is a list with the first element containing the names of groups and turnover information 
#for this group
#plot time interval SER for the second group
GrID = 1;
lSERa_gr$groupnames[GrID]
plot(lSERa_gr$turnover[[GrID]]$TimeIntv, lSERa_gr$turnover[[GrID]]$SER)

#group by area only
lSERa_gr = turnover_g(data[, SpecColumns], method = "SERa", dates =  SampleDates, groupby = StatArea[, "Area"])
GrID = 1;
lSERa_gr$groupnames[GrID]
plot(lSERa_gr$turnover[[GrID]]$TimeIntv, lSERa_gr$turnover[[GrID]]$SER)



```
