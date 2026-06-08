#' Alexey Ryabov 2026
#' Calculate community turnover (overlap) indices
#' see details in 
#' Hillebrand, H. et al. J Appl Ecol 55, 169-184 (2018).
#' Gotelli, N.J., Chao, A. & Colwell, R.K. (2024). Measuring and Estimating Species Richness, Taxonomic and Phylogenetic Diversity, and Related Biotic (Dis)similarity Indices From Sampling Data, 314–339.


#' 
#' @description
#' 
#' 
#' turnover_df = turnover(X, method = "SERr",  combinations = "i<j", 
#'               dates = NULL, locations = NULL, measure = "euclidean", 
#'               ext_inv = FALSE, groupby = NULL)
#' This function returns a data frame turnover_df of the community turnover between rows i and j 
#' of matrix X. turnover_df$From = i, turnover_df$To = j, turnover_df$SER is the turnover value between these rows. 
#' 
#' The function result contains:
#' turnover_df$From is the row number of 'from' observation in X
#' turnover_df$To is the row number of 'to' observation
#' turnover_df$SER is the turnover value between these observations
#' if ext_inv = TRUE, then the result includes additionally 
#' the effective number of common/extinct/immigrating/total species for 
#' observation i and j, turnover_df$S_common, turnover_df$S_nonshared,
#' turnover_df$S_ext, turnover_df$S_imm, turnover_df$S_total
#' e.g., turnover_df$S_total is the total number of species in observation i and j
#' and turnover_df$S_ext is the number of species observed in i but not in j.
#' turnover_df$S_nonshared is calculated as turnover_df$S_total - turnover_df$S_common.
#' S_ext and S_imm are calculated only for method="SERr" and method="SERa".
#' For method="Horn", S_ext and S_imm are returned as NA because Horn
#' dissimilarity is symmetric and does not naturally decompose into directional
#' "extinct" and "immigrant" components.
#' Note that the effective species number is richness based when method="SERr",
#' abundance (Simpson index) based when method="SERa", and Shannon/Horn based
#' when method="Horn". For method="U2", these effective species components are
#' returned as NA because U2 is a symmetric turnover (overlap) index rather than a
#' directional turnover decomposition.
#' 
#' if dates are specified, turnover_df$dateFrom and turnover_df$dateTo are the dates or time
#' values for observations i and j, and turnover_df$TimeIntv is the time interval
#' between them.
#' 
#' if locations are specified, turnover_df$Dist is added to the output. It contains
#' the spatial distance between the same pair of observations i and j. The
#' spatial distance is not used to calculate turnover_df$SER; it is returned so that
#' species turnover can be related to the spatial separation of communities.

#' turnover_df = turnover(X, method = "SERa", dates = NULL, locations = NULL, measure = "euclidean", groupby)
#' calculates turnover within groups defined by the groupby data frame or vector. 
#' groupby must contain the same number of rows as X. In this case turnover() 
#' splits X into groups and returns a data frame with the same structure as the 
#' ungrouped result, plus an additional field "groupname".
#' 
#' 
#' 
#' 
#' @param X  is  (M x N) data frame or matrix with M observations of N species. 
#'   X can contain species abundances, frequencies or, only for SERr, presence/absence data 
#'   if X contains abundance then for calculating SERa the abundances will be normalized, 
#'   so that the sum in each row equal 1. 
#' @param method  must be one of "SERr", "U0", "Horn", "U1", "SERa", and "U2". 
#' If method = "SERr", then richness based turnover values are calculated 
#' SER_ij  = (S_immigrant + S_extinct)/S_total, 
#' and the effective numbers of species are richness based. method = "U0" is an
#' alias for method = "SERr". It corresponds to the Hill-number order q = 0
#' Jaccard-type overlap family, returned here as the complement of Jaccard
#' similarity, i.e. as Jaccard dissimilarity.
#' If method = "Horn", then Horn dissimilarity is calculated from relative
#' abundances using Shannon entropy. S_ext and S_imm are returned as NA when
#' ext_inv = TRUE because Horn dissimilarity is symmetric and does not naturally
#' decompose into directional "extinct" and "immigrant" components. method =
#' "U1" is an alias for method = "Horn". It corresponds to the Hill-number
#' order q = 1 Shannon-entropy-based Horn measure, returned here as Horn
#' dissimilarity.
#'
#' If method = "SERa", then abundance based turnover values are calculated
#'  SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2  sum_k( pik  pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )
#' and the effective numbers of species are based on Simpson index. This is the
#' complement of Wishart similarity ratio, a Simpson-based abundance-weighted
#' analogue of Jaccard similarity as used by Hillebrand et al. (2018).
#'
#' If method = "U2", then the Hill-number order q = 2 turnover (overlap) index is
#' calculated from relative abundances:
#' U2_ij = 4 * sum_k(pik * pjk) /
#'         (sum_k(pik^2) + sum_k(pjk^2) + 2 * sum_k(pik * pjk)).
#' This is the q = 2 Jaccard-type overlap measure, also known as the regional
#' species-overlap measure. It is the q = 2 member of the Chao Uq family, not
#' the Wishart index. U2 is returned as a turnover (overlap) value, so identical
#' communities have U2 = 1. S_total, S_common, S_nonshared, S_ext, and S_imm are
#' returned as NA when ext_inv = TRUE because these turnover components are not
#' defined for U2 in this implementation.
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
#' @param dates Optional vector with M observation dates or numeric time values,
#' one value per row of X. If provided, the output includes dateFrom, dateTo,
#' and TimeIntv. dateFrom and dateTo are the original date or time values for
#' the two observations in each pair. TimeIntv is calculated as dates[j] -
#' dates[i]. For Date values the interval is in days; for numeric values it is
#' in the same units as the input, for example years.
#' 
#' @param locations Optional M x K matrix or data frame with spatial coordinates,
#' one row per observation in X. Typical inputs are M x 1 for a one-dimensional
#' gradient, M x 2 for two-dimensional coordinates, or M x 3 for three-dimensional
#' habitats. If measure = "lonlat", locations must have longitude in the first
#' column and latitude in the second column. If locations is provided, the output
#' includes Dist, the distance between the two communities in each returned pair.
#' 
#' @param measure Distance measure used for locations. It can be any method
#' accepted by stats::dist(), for example "euclidean", "manhattan", "maximum",
#' or "minkowski". If measure = "lonlat", geodesic distances in meters are
#' calculated with geodist::geodist(); in this case the geodist package must be
#' installed and locations must contain longitude and latitude in the first two
#' columns.
#' 
#' @param ext_inv Logical. If TRUE, the output includes effective numbers of
#' total, common, non-shared, extinct, and immigrant species: S_total, S_common,
#' S_nonshared, S_ext, and S_imm. S_nonshared is calculated as S_total -
#' S_common. S_ext and S_imm are calculated only for method = "SERr" and method
#' = "SERa". For method = "SERr" these values are richness-based counts; for
#' method = "SERa" they are abundance/effective-species components. For method =
#' "Horn", S_total, S_common, and S_nonshared are Shannon/Horn based, while
#' S_ext and S_imm are NA because Horn dissimilarity is symmetric and does not
#' naturally decompose into directional "extinct" and "immigrant" components.
#' For method = "U2", all five components are returned as NA because U2 is a
#' symmetric turnover (overlap) index and these effective turnover components are not
#' defined in this implementation.
#'  
#' @param groupby Optional vector or data frame with M rows used to split the
#' observations into groups, for example station, region, treatment, or any
#' combination of grouping variables. If provided, turnover is calculated only
#' within each group, not between groups. The output includes groupname, and
#' From and To keep the original row numbers from X.
#'  
#'  

#' @examples
#' 
#' richness based turnover characteristics (list of turnover values + other metrics)
#' turnover_df_r = turnover(data[, SpecColumns], method = "SERr", ext_inv = TRUE) 
  #abundance based turnover characteristics (list of turnover values + other metrics) 
#'  turnover_df_a = turnover(data[, SpecColumns], method = "SERa", ext_inv = TRUE)
#' 
#' turnover_s() is a convenience wrapper for species data only. It accepts X,
#' method, and combinations, and does not add dates, locations, or grouping.
#' TO_SERa = turnover_s(X, method = "SERa")
#' 
#' TO_SERr = turnover_s(X, method = "SERr")
#'
#' TO_Horn = turnover_s(X, method = "Horn")
#'
#' TO_U2 = turnover_s(X, method = "U2")
#' 
#' turnover_df = turnover(X, method = "SERa", dates = dates)
#' spatial turnover, when coordinates are longitude and lattitude 
#' turnover_df = turnover(X, method = "SERa", locations = LonLat, measure = "lonlat")
#' 
#' Calculating temporal species turnover separately for multiple 
#' stations 
#' turnover_df_gr = turnover(X, method = "SERa", dates = ObservationDates, groupby = StationList)
#' where StationList is a vector or data frame with station names for each observation. 
#' Output: turnover_df_gr$groupname contains station names, 
#' turnover_df_gr$SER is turnover information calculated within each station, but not between them. 
#' 


turnover_s <- function (X, method = "SERr", combinations = "i<j") {
  #a short wrapper function for turnover between rows of the X, 
  #when you have only species data, but no dates, coordinates etc. 

  #X is  (M x N) table or matrix with M observations of abundances of N species, 
  #M observation -- rows
  #N species -- columns
  
  #method = 'SERa' output is based on relative species frequencies
  # SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2 * sum_k( pik * pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )
  
  #method = 'SERr' or 'U0' output is based on presence/absence data
  #SER_ij  = (S_immigrant + S_extinct)/S_total
  #method = 'Horn' or 'U1' output is Horn dissimilarity based on relative abundances
  #method = 'U2' output is Hill-number order q = 2 overlap based on relative abundances
  # see a more detailed description of the input paramters in function turnover
  
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
    Xsplit = split(as.data.frame(X), groupby);
    if (!is.null(dates))     { DatesSplit = split(dates,     groupby)}
    DatesPiece = NULL;
    if (!is.null(locations)) { LocSplit   = split(as.data.frame(locations), groupby)}
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
    method_key <- method
    if (method == "U0") {
      method_key <- "SERr"
    } else if (method == "U1") {
      method_key <- "Horn"
    }
    
    switch(method_key, 
           SERa={
             # abundance based turnover...
             Res = turnover_a(X, method);
           },
           SERr={
             # richness based turnover...
             Res = turnover_r(X, method);
           },
           Horn={
             # Horn dissimilarity...
             Res = turnover_h(X, method);
           },
           U2={
             # Hill-number order q = 2 overlap...
             Res = turnover_u2(X, method);
           },
           {
             stop('turnover method is undefined. Use one of "SERr", "U0", "SERa", "Horn", "U1", or "U2".')
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
      S_nonshared = S_total - S_common;
      S_imm   = Res$S_imm[ord_ind];
      S_ext   = Res$S_ext[ord_ind];
      Result = data.frame(From, To, SER, S_total, S_common, S_nonshared, S_ext, S_imm);
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
      
      date_i = matrix(rep(seq_len(M), M), M, M, byrow = FALSE);
      date_j = matrix(rep(seq_len(M), M), M, M, byrow = TRUE);
      Result$dateFrom = dates[date_i[ord_ind]];
      Result$dateTo = dates[date_j[ord_ind]];
      Result$TimeIntv = Result$dateTo - Result$dateFrom;
    }
    
    if (!is.null(locations)){
      M <- nrow(X);
      M2 = nrow(locations)
      if (M != M2){
        stop("The number of rows in the location matrix must be equal to the number of rows in X.")
      }    
      if (tolower(measure)  == "lonlat") {
        if (!requireNamespace("geodist", quietly = TRUE)) {
          stop("Package 'geodist' is required for measure='lonlat'. Install it via install.packages('geodist').")
        }
        Dist <- geodist::geodist(locations, sequential = FALSE, paired = FALSE, measure = "geodesic");
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
  
  
  #to calculate abundance based turnover values, the data should be normalized
  #each row should contain species frequencies, for row j, \sum_i(p_ij)=1
  s <- rowSums(X, na.rm = TRUE);
  zero_rows <- which(is.na(s) | s == 0);
  if (length(zero_rows) > 0) {
    warning(paste0(
      "SERa undefined for rows with zero total abundance. Pairs involving these rows will be NA. Rows: ",
      paste(zero_rows, collapse = ", ")
    ));
    X[zero_rows, ] <- NA_real_;
    s[zero_rows] <- NA_real_;
  }
  Sums <- matrix(rep(s, N), M, N);
  X <- X / Sums;
  
  
  
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

turnover_u2 <- function(X, method = NULL) {
  X <- as.matrix(X)
  
  dims <- dim(X)
  M <- dims[1]
  N <- dims[2]
  
  if (any(X < 0, na.rm = TRUE)) {
    stop("U2 overlap is undefined for negative abundances.")
  }
  
  s <- rowSums(X, na.rm = TRUE)
  zero_rows <- which(is.na(s) | s == 0)
  if (length(zero_rows) > 0) {
    warning(paste0(
      "U2 overlap is undefined for rows with zero total abundance. Pairs involving these rows will be NA. Rows: ",
      paste(zero_rows, collapse = ", ")
    ))
    X[zero_rows, ] <- NA_real_
    s[zero_rows] <- NA_real_
  }
  
  Sums <- matrix(rep(s, N), M, N)
  P <- X / Sums
  
  P_sq <- as.matrix(rowSums(P^2))
  P_sq_rep <- rep(P_sq, M)
  P1_sq <- matrix(P_sq_rep, M, M, byrow = FALSE)
  P2_sq <- matrix(P_sq_rep, M, M, byrow = TRUE)
  P12 <- P %*% t(P)
  
  denominator <- P1_sq + P2_sq + 2 * P12
  SER <- 4 * P12 / denominator
  
  if (length(zero_rows) > 0) {
    SER[zero_rows, ] <- NA_real_
    SER[, zero_rows] <- NA_real_
  }
  
  valid_rows <- setdiff(seq_len(M), zero_rows)
  if (length(valid_rows) > 0) {
    diag(SER)[valid_rows] <- 1
  }
  
  S_total <- matrix(NA_real_, M, M)
  S_common <- matrix(NA_real_, M, M)
  S_imm <- matrix(NA_real_, M, M)
  S_ext <- matrix(NA_real_, M, M)
  
  Result <- list(
    SER = SER,           # U2 overlap / similarity
    S_total = S_total,   # not defined for U2 in this implementation
    S_common = S_common, # not defined for U2 in this implementation
    S_imm = S_imm,       # not defined for U2 in this implementation
    S_ext = S_ext        # not defined for U2 in this implementation
  )
  
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


turnover_h <- function(X, method = NULL) {  # Horn turnover (overlap) index, Hill = 1
  X <- as.matrix(X)
  
  dims <- dim(X)
  M <- dims[1]
  N <- dims[2]
  
  # Horn dissimilarity is defined for non-negative abundances
  if (any(X < 0, na.rm = TRUE)) {
    stop("Horn dissimilarity is undefined for negative abundances.")
  }
  
  # Normalize rows to species relative abundances
  s <- rowSums(X, na.rm = TRUE)
  zero_rows <- which(is.na(s) | s == 0)
  if (length(zero_rows) > 0) {
    warning(paste0(
      "Horn dissimilarity is undefined for rows with zero total abundance. Pairs involving these rows will be NA. Rows: ",
      paste(zero_rows, collapse = ", ")
    ))
    X[zero_rows, ] <- NA_real_
    s[zero_rows] <- NA_real_
  }
  
  Sums <- matrix(rep(s, N), M, N)
  P <- X / Sums
  
  # Row-wise Shannon entropy
  P_logP <- P
  ok <- !is.na(P_logP) & P_logP > 0
  P_logP[ok] <- P_logP[ok] * log(P_logP[ok])
  P_logP[!ok] <- 0
  H <- -rowSums(P_logP, na.rm = TRUE)
  if (length(zero_rows) > 0) H[zero_rows] <- NA_real_
  
  # Repeat row entropies into M x M matrices
  H1 <- matrix(rep(H, M), M, M, byrow = FALSE)
  H2 <- matrix(rep(H, M), M, M, byrow = TRUE)
  
  # Pairwise entropy of the mean composition:
  # Hm[i, j] = H( (P[i, ] + P[j, ]) / 2 )
  Hm <- matrix(0, M, M)
  
  for (k in seq_len(N)) {
    mk <- outer(P[, k], P[, k], "+") / 2
    contrib <- ifelse(mk > 0 & !is.na(mk), -mk * log(mk), 0)
    Hm <- Hm + contrib
  }
  
  # Propagate NA for invalid rows
  if (length(zero_rows) > 0) {
    Hm[zero_rows, ] <- NA_real_
    Hm[, zero_rows] <- NA_real_
    H1[zero_rows, ] <- NA_real_
    H1[, zero_rows] <- NA_real_
    H2[zero_rows, ] <- NA_real_
    H2[, zero_rows] <- NA_real_
  }
  
  # Horn dissimilarity for two equally weighted assemblages
  SER <- (Hm - (H1 + H2) / 2) / log(2)
  
  # Numerical protection
  SER[SER < 0] <- 0
  SER[SER > 1] <- 1
  
  # Effective shared components in units of effective species
  S_total <- exp(Hm)
  S_common <- (1 - SER) * S_total
  S_imm <- matrix(NA_real_, M, M)
  S_ext <- matrix(NA_real_, M, M)
  
  # Clean diagonals only for valid rows
  valid_rows <- setdiff(seq_len(M), zero_rows)
  diag_vals <- exp(H)
  
  if (length(valid_rows) > 0) {
    diag(SER)[valid_rows] <- 0
    diag(S_total)[valid_rows] <- diag_vals[valid_rows]
    diag(S_common)[valid_rows] <- diag_vals[valid_rows]
  }
  
  Result <- list(
    SER = SER,                 # Horn dissimilarity
    S_total = S_total,         # effective total
    S_common = S_common,       # effective shared component
    S_imm = S_imm,             # not defined in Horn framework
    S_ext = S_ext              # not defined in Horn framework
  )
  
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
