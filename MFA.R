mfa <- function(data, sets, ncomps = NULL, center = TRUE, scale = TRUE) {
  # R package for Multiple Factor Analysis
  
  # Input:
  # data - data set
  # sets - list of vectors indicating the sets of variables
  # ncomps - integer, how many number of components
  # center - as R's scale() 
  # scale - as R's scale()
  
  # Return:
  # 
  # MatrixA  - Matrix with eigenvalues (Variances)
  # MatrixU  - Matrix U of the SVD of Matrix Z
  # MatrixV  - Matriz V of the SVD of Natrix Z
  # MatrixF  - Compromise Factor Score Matrix
  # MatrixEFG - Partial Factor Score Matrix
  # MatrixCCP - Matrix with Correlation of Principal Components with Groups
  # MatrixEscVar - Matrix of Partial Inertia
  
  MBQ <- function(DataQ, center=TRUE, scale=TRUE) {  
    #center, scale and get the weight of M and A
    
    # Input:
    # DataQ - raw data to be dealt with
    # center - whether centering
    # scale - whether scale
    
    # Returns:
    # MZ   - centered and scaled matrix(according to choice)
    # PLin - line weights
    # PCol - column weights
    
    MZ <- NULL    # 
    PLin <- NULL   # 
    PCol <- NULL   # 
    
    ### Begin - Center and scale data  ###
    MC <- as.matrix(DataQ) # 
    
  if(center==TRUE)
  {
    Media <- apply(MC,2,mean) # Matrix with averages by columns
    
    MC <- sweep(MC, 2, Media, FUN = "-") # Centering
  }
    
  if(scale==TRUE){
    
    Media <- apply(MC,2,mean) # Matrix with averages by columns
    
    MC_after <- sweep(MC, 2, Media, FUN = "-") # Centering
    
    SqSum <- sqrt(colSums(MC_after^2))
    
    MC <- sweep(MC, 2, SqSum, FUN = "/") # Scaling
  }
    
    ### End - Center and scale ###  
      
    PLin <- rep(1/nrow(MC),nrow(MC))
    
    Pe <- (svd(MC)$d[1])^2   # find the first eigenvalue of MC
    
    PCol <- cbind(PCol,t(rep(1/Pe,ncol(MC)))) # matrix weights of column
    
    Lista <- list(MZ=MC, PLin=PLin, PCol=PCol)
    
    return(Lista)
  }
  

  if (is.null(sets)) # Creates names for the variables if not exist
    sets <- paste("Variable", 1:10, sep = " ")
  
  ### Begin - Get the weight of the groups of variables
  NumSets = length(sets) # number of groups formed
  
  MZG   <- NULL  # general Z null array
  PLinG <- NULL  # general matrix with null line weights
  PColG <- NULL  # general matrix with null column weights
  
  j  <- 1       # initial column of the variable group
  
  k  <- sets[1] # final column of the variable group
  
  for (i in 1:NumSets) {
    
    MB   <- MBQ(data[,j:k],center=center, scale=scale)
    MZ   <- MB$MZ
    PLin <- MB$PLin
    PCol <- MB$PCol
    colnames(PCol) <- colnames(data[,j:k])

    PLinG <- PLin  # general matrix with line weights
    
    PColG <- cbind(PColG,PCol) # general matrix with column weights
    
    MZG   <- cbind(MZG,MZ)     # centered and scaled matrix
    
    j <- j + sets[i]    # begin column of next variable group
    
    k <- k + sets[i+ifelse(i!=NumSets,1,0)]  # final column of next group
    

  }
  
  PColG <- t(PColG)  
  ### End - Get the weight of the groups of variables ###
  
  ### Begin - Find the eigenvectors and eigenvalues ###
  MDS <- GSVD(MZG, PLinG, PColG) # 
  MAutoVlr  <- MDS$d  # 
  MAutoVecU <- MDS$u  # 
  MAutoVecV <- MDS$v  # 
  
  ## Matrix of variances
  MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 5))
  rownames(MEigen) <- paste("Axis", 1:length(MAutoVlr))
  colnames(MEigen) <- c("Singular value","Eigenvalue", "% cumulative eigenvalue","% Inertia","% cumulative Inertia")
  MEigen[, "Singular value"] <- MAutoVlr
  MEigen[, "Eigenvalue"] <- MAutoVlr^2
  MEigen[, "% cumulative eigenvalue"] <- cumsum(MEigen[, "Eigenvalue"])
  MEigen[, "% Inertia"] <- (MAutoVlr^2/sum(MAutoVlr^2)) * 100
  MEigen[, "% cumulative variance"] <- cumsum(MEigen[,"% Inertia"])
  
  NumAutoVlr <- length(MAutoVlr) # Number of variables
  
  NE <- length(MAutoVlr[MAutoVlr>1e-10]) # Number of significant elements
  ### End - Find the eigenvectors and eigenvalues ###
  
  ### Begin - Matrix of compromise factor score ###
  MF <-  MAutoVecU[,1:NE]%*%diag(MAutoVlr[1:NE],NE) # Matrix F - compromise factor score
  rownames(MF) <- rownames(data) # name the rows
  colnames(MF) <- paste("Axis", 1:ncol(as.matrix(MF)), sep = " ") # name the columns
  ### End - Matrix of compromise factor score ###
  
  ### Begin - Matrix of partial factor scores ###
  j  <- 1        
  
  k  <- sets[1] 
  
  LMFGrupo <- as.list(1:NumSets) 
  
  for (i in 1:NumSets) {       
    
    MFG <- NumSets * MZG[,j:k]
    
    MFG <- sweep(MFG, 2, PColG[j:k], FUN="*")
    
    LMFGrupo[[i]] <- MFG%*%MAutoVecV[j:k,] # Creates matrix of partial factor scores by group
    
    colnames(LMFGrupo[[i]]) <- paste("Axis", 1:ncol(as.matrix(LMFGrupo[[i]])), sep = " ") # name the columns
    
    j <- j + sets[i]      #Initial column of the variable Group
    
    k <- k + sets[i+ifelse(i!=NumSets,1,0)]  # Columns of the variable Group
  }
  
  names(LMFGrupo) <- paste("Group", 1:NumSets, sep = "") # name of group
  ### End - Matrix of partial factor scores ###
  
  ### Begin -  Correlation of Principal Components with Original Variables ###
  CCP <- sweep(as.matrix(MAutoVecV), 2, MAutoVlr, FUN = "*")  
  CCP <- t(CCP)
  rownames(CCP) <- paste("Axis", 1:NumAutoVlr, sep = " ")
  colnames(CCP) <- colnames(MZG)
  ### End -  Correlation of Principal Components with Original Variables ###
  
  ### Begin - Partial matrix of inertia/scores of variables ###
  CoordVar <- sweep(as.matrix(MAutoVecV), 2, sqrt(MAutoVlr), FUN = "*")  # Coordenadas das variaveis
  
  ContrVar <- sweep(as.matrix(CoordVar^2), 2, MAutoVlr, "/") # Contribuicao das variaveis
  
  ContrVar <- sweep(as.matrix(ContrVar), 1, PColG, "*")
  
  ContrGru <- matrix(data = NA, nrow = NumSets, ncol = NumAutoVlr) # Matriz com Contribuicoes dos Grupos
  
  j  <- 1      
  
  k  <- sets[1]
  
  for (i in 1:NumSets) {
    
    ContrGru[i,] <- apply(ContrVar[j:k, ], 2, sum) # Matrix with contribution of groups
    
    j <- j + sets[i]      
    
    k <- k + sets[i+ifelse(i!=NumSets,1,0)]  
    
  }
  
  EscVar <- sweep(ContrGru, 2, MAutoVlr^2, "*") # create matrix of variables scores/partial inertia
  
  colnames(EscVar) <- paste("Axis", 1:ncol(as.matrix(EscVar)), sep = " ") # name the columns
  
  #rownames(EscVar) <- c("group 1", "group 2") # name of rows ???
  ### End - Partial matrix of inertia/scores of variables ###
  
  if(is.null(ncomps))
   
    Lista <- list(MatrixPLin = PLinG,
                  MatrixPCol = PColG, MatrixZ = MZG, MatrixA = MEigen,
                  MatrixU = MAutoVecU, MatrixV = MAutoVecV, MatrixF = MF, 
                  MatrixEFG = LMFGrupo, MatrixCCP = CCP, MatrixEscVar = EscVar)

  else {   
  
  npart<-NULL
  
  j<-1
  
  k<-ncomps
  
  nparte<-as.vector(j:k)
  
  for (i in 1:NumSets) {
    
    npart <-c(npart,nparte)
    
    j <- j + NE      
    
    k <- k + NE  
    
    nparte<-as.vector(j:k)
  }
  Lista <- list(MatrixPLin = PLinG,
                MatrixPCol = PColG, MatrixZ = MZG, MatrixA = MEigen[1:ncomps,],
                MatrixU = MAutoVecU[,1:ncomps], MatrixV = MAutoVecV, MatrixF = MF[,1:ncomps], 
                MatrixEFG = as.data.frame(LMFGrupo)[,npart], MatrixCCP = CCP, MatrixEscVar = EscVar[,1:ncomps])
  }
  return(Lista)
}
