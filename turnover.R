#' Alexey Ryabov 2020
#' Calculate species turnover index
#' see details in 
#' Hillebrand, H. et al. J Appl Ecol 55, 169-184 (2018).

#' 
#' @description
#' 
#' 
#' lSER = turnover(X, method = "SERr",  combinations = "i<j", 
#'               dates = NULL, locations = NULL, measure = "euclidean", 
#'               ext_inv = FALSE, groupby = NULL)
#' This function returns a data frame lSER of the turnover between rows i and j 
#' of matrix X. lSER$From = i, lSER$To = j, lSER$SER is the turnover between these rows. 
#' 
#' The function result contains:
#' lSER$From is the row number of 'from' observation in X
#' lSER$To is the row number of 'to' observation
#' lSER$SER is the turnover between these observations
#' if ext_inv = TRUE, then the result includes additionally 
#' the effective number of common/extinct/immigrating/total species for 
#' observation i and j, lSER$S_common, lSER$S_ext, lSER$S_imm, lSER$S_total 
#' e.g., lSER$S_total is the total number of species in observation i and j
#' and  lSER$S_ext the number of species observed in i but not in j
#' Note that the effective species number is richness based when  method="SERr" and 
#' abundance (Simpson index) based when method="SERa", 
#' 
#' if dates are specified, lSER$TimeIntv is the time interval between observation i and j, 
#' 
#' if locations are specified, lSER$Dist is the spatial distances between observation i and j, 

#' lSER = turnover(X, method = "SERa", dates = NULL, locations = NULL, measure = "euclidean", groupby)
#' calculates turnover within groups defined by groupby dataframe. 
#' groupby must contain the same number of rows as X. Function turnover_g(X, ...) splits X 
#' into groups and returns a data frame with the same structure as returned by  function turnover(X, ...), 
#' plus an additional field "groupname"   
#' 
#' 
#' 
#' 
#' @param X  is  (M x N) data frame or matrix with M observations of N species. 
#'   X can contain species abundances, frequencies or, only for SERr, presence/absence data 
#'   if X contains abundance then for calculating SERa the abundances will be normalized, 
#'   so that the sum in each row equal 1. 
#' @param method  must be one of "SERr" and  "SERa". 
#' If method = "SERr", then richness based turnover index is calculated 
#' SER_ij  = (S_immigrant + S_extinct)/S_total, 
#' and the effective numbers of species are richness based
#' 
#' If method = "SERa", then  abundance based turnover index is  calculated
#'  SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2  sum_k( pik  pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )
#' and the effective numbers of species are based on Simpson index
#'  
#' @param combinations	indicates for which combinations of 
#' rows i and j of matrix X the turnover must be calculated. The available 
#' values are "i<j", "i!=j" and "i,j". 
#' For combinations="i<j", the turnover is calculated between each row i and all subsequent 
#' rows j of matrix X. This type is useful for calculating temporal turnover, 
#' because in this case (if the data are sorted by date) we will calculate 
#' turnover from past to future, but not vice versa.  
#' For combinations="i!=j", turnover is calculated for all combinations of i and j, 
#' except for i=j. 
#' For combinations="i,j", the turnover is computed for all combinations 
#' of i and j, in which case the turnover for i=j is zero. This option is useful 
#' if the results should be converted later into a square matrix T, such that 
#' element T[i,j] shows the turnover between rows X[i, :] and X[j, :].
#'  
#' @param dates observation dates (datetime) or any numbers e.g. observation years, 
#' 
#' @param locations spatial coordinates, typically should be Mx1, for 1D gradients, 
#' Mx2 matrix if 2 coordinates are know, or Mx3 matrix for a 3D habitats
#' or  longitude and latitude (in this order!) in this case measure = "lonlat"
#' 
#' @param measure measure for measuring distance can be any of measures used 
#' by dist(X, ..) function, or can be "lonlat" in this case the geodesic distance 
#' in meters is measured  using 'geodist' package. 
#' This package should be installed before. the 1st column of locations should 
#' contain longitude and the second column latitude
#' 
#' @param ext_inv should the information about invaded, extinct, common and total 
#' species number be included?
#'  
#' @param groupby  is a column or data frame with M rows  which are used to  
#'  group the data  (see split(x, ...) function). In this case turnover is 
#'  calculated within groups defined by "groupby" and resutl contains column 
#'  groupname.
#'  
#'  

#' @examples
#' SER = turnover_s(X, method = "SERa")
#' 
#' SER = turnover_s(X, method = "SERr")
#' 
#' lSER = turnover_s(X, method = "SERa", dates)
#' spatial turnover, when coordinates are longitude and lattitude 
#' lSER = turnover_s(X, method = "SERa", location = LonLat, measure = "lonlat")
#' 
#' Calculating temporal species turnover separately for multiple 
#' stations 
#' lSER_gr = turnover_g(X, method = "SERa", dates = ObservationDates, groupby = StationList)
#' where StationList is a data frame with station names for each observation. 
#' Output: lSER_gr$groupname contains station names, 
#' lSER_gr$SER is turnover information calculated within each station, but not between them. 
#' 


turnover_s <- function (X, method = "SERr", combinations = "i<j") {
  #short version, returns only turnover index
  #between rows of the X 
  
  #X is  (M x N) table or matrix with M observations of abundances of N species, 
  #M observation -- rows
  #N species -- columns

  #method = 'SERa' output is based on relative species frequencies
  # SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2 * sum_k( pik * pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )

  #method = 'SERr' output is based on presence/absence data
  #SER_ij  = (S_immigrant + S_extinct)/S_total
  
  M = nrow(X); 
  Res = turnover(X, method,  combinations)
  return(Res)
}




turnover <- function (X, method = "SERr",  combinations = "i<j", 
                      dates = NULL, locations = NULL, measure = "euclidean", 
                      ext_inv = FALSE, groupby = NULL, RowID = NULL) {
  M = nrow(X); 
  if (!is.null(groupby))  {#if groups, then split and call turnover for each group
    RowID = c(1:M); #assign row ID
    RowIDSplit = split(RowID, groupby);
    Xsplit = split(X, groupby);
    if (!is.null(dates))     { DatesSplit = split(dates,     groupby)}
    DatesPiece = NULL;
    if (!is.null(locations)) { LocSplit   = split(locations, groupby)}
    LocPiece = NULL;  
    Res = NULL;
    GroupNames = names(Xsplit)
    for(i in 1:length(Xsplit)){
      if (!is.null(dates))      { DatesPiece = DatesSplit[[i]]}
      if (!is.null(locations))  { LocPiece =   LocSplit[[i]]}
      #get turnover
      TO =  turnover(Xsplit[[i]], method, combinations, DatesPiece, LocPiece, measure, ext_inv, NULL, RowIDSplit[[i]]);
      #group name
      TO$groupname = rep(GroupNames[[i]], nrow(TO));
      #rbind results
      if (is.null(Res)) {
        Res = TO;
      }else {
        Res <- rbind(Res, TO)
      }
    }
    return(Res)  
  }else {#if no groups, then calculate turnover
    switch(method, 
           SERa={
             # abundance based turnover...
             Res = turnover_a(X, method);
           },
           SERr={
             # richness based turnover...
             Res = turnover_r(X, method);
           },
           {
             stop('turnover method is undefined')
           }
    )  
    
    if(is.null(RowID)) { #assign rowID is not assigned before
      RowID = c(1:M);
    }
    RowIDs = rep.int(RowID, M);
    i_ind = matrix(RowIDs, nrow=M,byrow = FALSE);
    j_ind = matrix(RowIDs, nrow=M,byrow = TRUE);

  
        
    ord_ind = mat_index(Res$SER, combinations);
    
    From    = i_ind[ord_ind];
    To      = j_ind[ord_ind];
    SER     = Res$SER[ord_ind];
    if (ext_inv){
      S_total = Res$S_total[ord_ind];
      S_common= Res$S_common[ord_ind];
      S_imm   = Res$S_imm[ord_ind];
      S_ext   = Res$S_ext[ord_ind];
      Result = data.frame(From, To, SER, S_total, S_common, S_ext);
    }
    else {
      Result = data.frame(From, To, SER);
    }

   if (!is.null(dates)){
      #get number of records 
      M <- nrow(X);
      M2 = length(dates)
      if (M != M2){
        stop("The length of the dates array must be equal to the number of rows in X")
      }
      
      #find differences in dates (if not null) or interval between columns
      D = rep(dates, M);
      D1 = matrix(D, M, M, byrow = FALSE);
      D2 = matrix(D, M, M, byrow = TRUE);
      
      #D1 = as.Date(D1, origin ="1970-01-01")  unwrap D1 and D2
      #D2 = as.Date(D2, origin ="1970-01-01")
      #TimeIntv = difftime(D2, D1, units);
      #need to pack to back into a matrix
      #TimeIntv = matrix(TimeIntv, M, M, byrow = FALSE)
      
      TimeIntv = D2 - D1; #interval in days or integer numbers
      #set elements below the main diagonal to zero
      Result$TimeIntv = TimeIntv[ord_ind];
    }
    
    if (!is.null(locations)){
      M <- nrow(X);
      M2 = nrow(locations)
      if (M != M2){
        stop("The number of rows in the location matrix must be equal to the number of rows in X.")
      }    
      if (tolower(measure)  == "lonlat") {
        library("geodist");
        Dist = geodist(locations, paired = TRUE, measure = "geodesic");
        Result$Dist = Dist[ord_ind];
      }else{
        Dist = as.matrix(dist(locations, method = measure, diag = TRUE, upper = TRUE, p = 2))
        Result$Dist = Dist[ord_ind];
      }
    }  
    return(Result);
  }
}


  

turnover_a <- function(X, method) {
  X = as.matrix(X)
  
  dims = dim(X)
  M = dims[1];
  N = dims[2];

    
  #to calculate abundance based turnover index, the data should be normalized
  #each row should contain species frequencies, for row j, \sum_i(p_ij)=1
  s = rowSums(X);
  Sums = matrix(rep(s, N), M, N);
  X = X/Sums;
  
  
  
  SimpsInd = as.matrix(rowSums(X^2));  #sum of squares, Simpson index
  P_sq =rep(SimpsInd, M);
  P1_sq =matrix(P_sq, M, M, byrow = FALSE);
  S1 = 1/P1_sq; #in each row the number of species in row i
  P2_sq =matrix(P_sq, M, M, byrow = TRUE);
  S2 = 1/P2_sq; #in each column the number of species in row i
  
  P12  = X %*% t(X);  #overlap 

  #find invaded, common and etc species 
  S_common  = (P12)/(P1_sq * P2_sq) ;
  S_ext   = S1 - S_common;
  S_imm   = S2 - S_common;
  S_total = S1 + S2 - S_common;

  SER = (S_imm + S_ext)/S_total;
  #SER2 = (P1_sq + P2_sq - 2 * P12)/(P1_sq + P2_sq - P12);
  
  Result = list(SER=SER, S_total = S_total, S_common = S_common, S_imm = S_imm, S_ext = S_ext);
 
  return(Result)
}



turnover_r <- function(X, method) {
  dims = dim(X)
  M = dims[1];
  N = dims[2];
  
  P1 = as.matrix(X > 0)*1;
  P2 = t(P1);

  S_common = P1 %*% P2;

  S_sample = as.matrix(rowSums(X>0));  #number of species in each sample
  S_sample =rep(S_sample, M);
  S1 =matrix(S_sample, M, M, byrow = FALSE); #put number of species columnwise  
  S2 =matrix(S_sample, M, M, byrow = TRUE);  #put number of species rowwise 
  
  S_total = S1 + S2 - S_common;
  S_imm = ((!P1) + 0) %*% (P2);
  S_ext =  P1 %*% ((!P2) + 0);
  SER = (S_imm + S_ext)/(S_total);
  Result = list(SER=SER, S_total = S_total, S_common = S_common, S_imm = S_imm, S_ext = S_ext);
  return(Result)
}

# Mat2DataFrame <- function(Res, diag, upper, xnames){
#   #Set to NULL all elements below the main diagonal, because this elements correspond  
#   #to 'back-in-time' turnover rates
#   M <- nrow(Res$SER);
#   
#   
#   ij = rep.int(c(1:M), M);
#   i_ind = matrix(ij, nrow=M,byrow = FALSE);
#   j_ind = matrix(ij, nrow=M,byrow = TRUE);
#   
#   ind = mat_index(Res$SER, diag, upper);
# 
#   From    = i_ind[ind];
#   To      = j_ind[ind];
#   SER     = Res$SER[ind];
#   S_total = Res$S_total[ind];
#   S_common= Res$S_common[ind];
#   S_imm   = Res$S_imm[ind];
#   S_ext   = Res$S_ext[ind];
#   
#   return(data.frame(From, To, SER, S_total, S_common, S_ext));
# }

mat_index <- function(SampleMatrix, combinations){
  if (!is.na(pmatch(combinations, "i<j"))) 
    combinations <- "i<j"
  COMBs <- c("i<j", "i!=j", "i,j")
  combID <- pmatch(combinations, COMBs)
  if (is.na(combID)) 
    stop("invalid i,j combinations")
  if (combID == -1) 
    stop("ambiguous i,j combinations")
  
  indU = upper.tri(SampleMatrix, diag=FALSE);
  if (combID == 1){ #upper part
    ind = indU;
  } else if (combID == 2) {  #upper and lower part
    indL = lower.tri(SampleMatrix, diag=FALSE);
    ind = indL | indU;
  } else { #everything
    indL = lower.tri(SampleMatrix, diag=TRUE);
    ind = indL | indU;
  }
          
return(ind);
}