GSVD <- function(Data, PLin = NULL, PCol = NULL) {
  # Function that performs Generalized Singular Value Decomposition
  
  # input:
  # Data - raw Matrix used for Decomposition
  # PLin - Vector with weights for lines  (in the paper, is m=rep(1/12,12))
  # PCol - Vector with weights for columns (in the paper, is a<-c(rep(0.241,6), rep(0.239,6), rep(0.275,6), rep(0.273,5),rep(0.307,6),rep(0.302,5), rep(0.417,4), rep(0.272,6), rep(0.264,5),rep(0.309,4)))
  
  # Returns:
  # d - Vector line with the singular values of the result
  # u - Line eigenvectors
  # v - Column relative eigenvectors
  
  if (is.null(PCol)) PCol <- rep(1, ncol(Data))
  
  if (is.null(PLin)) PLin <- rep(1, nrow(Data))
  
  else if (is.numeric(PLin)) PLin = PLin / sum(PLin)
  
  if (!is.numeric(PLin))
    stop("Input to 'PLin' must be of the numeric vector type. Check!")
  
  if (!is.numeric(PCol))
    stop("Input to 'PCol' must be of the numeric vector type. Check!")
  
  if (nrow(Data) != length(PLin))
    stop("The number of elements in 'Plin' must be equal to the number of lines of the 'Data' component. Check!")
  
  if (ncol(Data) != length(PCol))
    stop("The number of elements in 'PCol' must be equal to the number of columns of the 'Data' component. Check!")
  
  PLin <- as.vector(PLin)
  
  PCol <- as.vector(PCol)
  
  ncv <- min(nrow(Data)-1,ncol(Data)) # Number of valid columns
  
  AA = sweep(Data, 2, sqrt(PCol), FUN = "*")
  
  AA = sweep(AA, 1, sqrt(PLin), FUN = "*")
  
  MSVD <- svd(AA)
  d <- MSVD$d
  MU <- MSVD$u
  MV <- MSVD$v
  
  P <- diag(sqrt(1/PLin))%*%MU
  
  Q <- diag(sqrt(1/PCol))%*%MV
  
  Resp <- list(d = d[1:ncv], u = P[,1:ncv], v = Q[,1:ncv])
  
  return(Resp)
}
