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
  # MatrixF  - Common Factor Score Matrix
  # MatrixEFG - Partial Factor Score Matrix
  # MatrixCCP - Matrix with Correlation of Principal Components with Groups
  # MatrixEscVar - Matrix of Partial Inertia
  
  MBQ <- function(DataQ,PondGeral) {  
    # balancing quantitative data
    
    # Input:
    # DataQ - data to be balanced
    # PondGeral - for frequencies data
    
    # Returns:
    # MZ   - balanced matrix
    # PLin - line weights
    # PCol - column weights
    
    MZ <- NULL    # 
    PLin <- NULL   # 
    PCol <- NULL   # 
    
    ### Begin - Center and scale data  ###
    MC <- as.matrix(DataQ) # 
    
    Media <- apply(MC,2,mean) # Matrix with averages by columns
    
    MC <- sweep(MC, 2, Media, FUN = "-") # Center
    
    SqSum <- sqrt(colSums(MC^2)/nrow(MC))
    
    MC <- sweep(MC, 2, SqSum, FUN = "/") # Scale
    ### End - Center and scale ###  
    if (TRUE) { # don't process frequence now
      MC <- as.matrix(DataQ)
      
      PLin <- rep(1,nrow(MC))
      
      NLin <- nrow(MC)
      
      SCol1 <- colSums(MC) / NLin
      
      MC <- sweep(MC, 2, SCol1, FUN = "-")
      
      SCol2 <- sqrt(colSums(MC^2)/NLin)
      
      MC <- sweep(MC, 2, SCol2, FUN = "/")
    }
    
    Pe <- GSVD(MC,PLin,rep(1,ncol(MC)))$d[1]^2   # find the first eigenvalue of MC
    
    PCol <- cbind(PCol,t(rep(1/Pe,ncol(MC)))) # matrix weights of column
    
    Lista <- list(MZ=MC, PLin=PLin, PCol=PCol)
    
    return(Lista)
  }
  

  if (is.null(sets)) # Creates names for the variables if not exist
    sets <- paste("Variable", 1:10, sep = " ")
  
  ### Begin - Balance the values of the groups of variables
  NumGrupos = length(sets) # number of groups formed
  
  MZG   <- NULL  # general Z null array
  PLinG <- NULL  # general matrix with null line weights
  PColG <- NULL  # general matrix with null column weights
  
  j  <- 1       # initial column of the variable group
  
  k  <- sets[1] # final column of the variable group
  
  for (i in 1:NumGrupos) {
    
    if (TRUE) {  # we only process quantitative data now
      MB   <- MBQ(data[,j:k],PondGeral)
      MZ   <- MB$MZ
      PLin <- MB$PLin
      PCol <- MB$PCol
      colnames(PCol) <- colnames(data[,j:k])
    }
    
    
    PLinG <- PLin  # general matrix with line weights
    
    PColG <- cbind(PColG,PCol) # general matrix with column weights
    
    MZG   <- cbind(MZG,MZ)     # balanced matrix
    
    j <- j + sets[i]    # begin column of next variable group
    
    k <- k + sets[i+ifelse(i!=NumGrupos,1,0)]  # final column of next group
    

  }
  
  PColG <- t(PColG)  
  ### End - Balance the values of the groups of variables ###
  
  ### Begin - Find the eigenvectors and eigenvalues ###
  MDS <- GSVD(MZG, PLinG, PColG) # 
  MAutoVlr  <- MDS$d  # 
  MAutoVecU <- MDS$u  # 
  MAutoVecV <- MDS$v  # 
  
  ## Matrix of variances
  MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 3))
  rownames(MEigen) <- paste("Axis", 1:length(MAutoVlr))
  colnames(MEigen) <- c("Eigenvalue", "% variance","% cumulative variance")
  MEigen[, "Eigenvalue"] <- MAutoVlr^2
  MEigen[, "% variance"] <- (MAutoVlr^2/sum(MAutoVlr^2)) * 100
  MEigen[, "% cumulative variance"] <- cumsum(MEigen[,"% variance"])
  
  NumAutoVlr <- length(MAutoVlr) # Number of car valores
  
  NE <- length(MAutoVlr[MAutoVlr>1e-10]) # Number of significant elements
  ### End - Find the eigenvectors and eigenvalues ###
  
  ### Begin - Matrix Global score ###
  MF <-  MAutoVecU[,1:NE]%*%diag(MAutoVlr[1:NE],NE) # Matrix F - global factor score
  rownames(MF) <- rownames(data) # Nomeia as linhas
  colnames(MF) <- paste("Axis", 1:ncol(as.matrix(MF)), sep = " ") # Nomeia as colunas
  ### End - Matrix Global score ###
  
  ### Begin - Matrix of factor scores by group ###
  j  <- 1        # 
  
  k  <- sets[1] # 
  
  LMFGrupo <- as.list(1:NumGrupos) # 
  
  for (i in 1:NumGrupos) {       
    
    MFG <- NumGrupos * MZG[,j:k]
    
    MFG <- sweep(MFG, 2, PColG[j:k], FUN="*")
    
    LMFGrupo[[i]] <- MFG%*%MAutoVecV[j:k,] # cria Matriz dos Escores dos Fatores por Grupo
    
    colnames(LMFGrupo[[i]]) <- paste("Axis", 1:ncol(as.matrix(LMFGrupo[[i]])), sep = " ") # Nomeia as colunas
    
    j <- j + sets[i]      # coluna inicial do Grupo de variaveis
    
    k <- k + sets[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do Grupo de variaveis  
  }
  
  names(LMFGrupo) <- paste("Group", 1:NumGrupos, sep = "") # name of group
  ### End - Matrix of factor scores by group ###
  
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
  
  ContrGru <- matrix(data = NA, nrow = NumGrupos, ncol = NumAutoVlr) # Matriz com Contribuicoes dos Grupos
  
  j  <- 1        # coluna inicial do Grupo de variaveis
  
  k  <- sets[1] # coluna final do Grupo de variaveis
  
  for (i in 1:NumGrupos) {
    
    ContrGru[i,] <- apply(ContrVar[j:k, ], 2, sum) # Matriz com Contribuicoes dos Grupos
    
    j <- j + sets[i]      # coluna inicial do Grupo de variaveis
    
    k <- k + sets[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do Grupo de variaveis  
    
  }
  
  EscVar <- sweep(ContrGru, 2, MAutoVlr^2, "*") # create matrix of variables scores/partial inertia
  
  colnames(EscVar) <- paste("Axis", 1:ncol(as.matrix(EscVar)), sep = " ") # Nomeia as colunas
  
  #rownames(EscVar) <- c("group 1", "group 2") # name of lines ???
  ### End - Partial matrix of inertia/scores of variables ###
  
  Lista <- list(MatrixPLin = PLinG,
                MatrixPCol = PColG, MatrixZ = MZG, MatrixA = MEigen,
                MatrixU = MAutoVecU, MatrixV = MAutoVecV, MatrixF = MF, 
                MatrixEFG = LMFGrupo, MatrixCCP = CCP, MatrixEscVar = EscVar)
  
  return(Lista)
}
  
  