Plot.MFA <- function(MFA,Titles = matrix(NA,1,3), PosLeg=2, BoxLeg="s", Color="s",NamArr="n") {
  # Routine to Plot Graphics of the developed MFA Method
  # Input:
  # MF - MFA function data
  # Titles - Titles for graphics. If it is not defined, it will default to default text.
  # PosLeg - 1 for caption in upper left corner
  # 2 for caption in the upper right corner - default
  # 3 for caption in lower right corner
  # 4 for caption in lower left corner
  # BoxLeg - "s" to put frame in legend - default
  # "N" does not put frame in caption
  # Color - "s" for colored graphics - default
  # "N" for black and white graphics
  # NamArr - "s" to put names points in the cloud around the
  # Centroid in the Graph corresponding to
  # Global Analysis of Individuals and Variables
  # "N" Otherwise - default
  
  # Returns:
  # Various graphics
  
  ##### begin - Information Used in Graphics #####
  # Creates Titles for graphics if they do not exist  if (!is.character(Titles[1]) || is.na(Titles[1])) Titles[1] = c("Graphic Corresponding to Global Analysis of Individuals")
  if (!is.character(Titles[2]) || is.na(Titles[2])) Titles[2] = c("Graphic Corresponding to Global Analysis of Individuals and Variables")
  if (!is.character(Titles[3]) || is.na(Titles[3])) Titles[3] = c("Graph of inertia of Variable Groups")
  
  Color  = ifelse(Color=="s","S",ifelse(Color=="n","N",Color))    # transforma em maiusculo
  BoxLeg = ifelse(BoxLeg=="s","S",ifelse(BoxLeg=="n","N",BoxLeg)) # transforma em maiusculo
  NamArr = ifelse(NamArr=="s","S",ifelse(NamArr=="n","N",NamArr)) # transforma em maiusculo
  
  if (PosLeg<1 || PosLeg>4)
    stop("Input to legend position (PosLeg) is incorrect. Check!")
  
  if (BoxLeg!="S" && BoxLeg!="N") 
    stop("Input to legend frame (BoxLeg) is incorrect. Check!")
  
  if (Color!="S" && Color!="N") 
    stop("Input to 'Color' is incorrect. Check!")
  
  if (NamArr!="S" && NamArr!="N") 
    stop("Input to 'NamArr' is incorrect. Check!")
  
  Grupos     = MFA$MatrixG  # Size of each group
  NomeGrupos = MFA$MatrixNG # Names of each group
  NomeLinhas = rownames(MFA$MatrixF) # names of line
  NumGrupos  = length(NomeGrupos) # number of groups
  cor        = 1 # 
  DescEixo1  = paste("First Principal Component (",round(MFA$MatrixA[1,2],2),"%)",sep="")
  DescEixo2  = paste("Second Principal Component (",round(MFA$MatrixA[2,2],2),"%)",sep="")
  
  if (PosLeg==1) PosLeg="topleft"     # posicao das legendas nos graficos
  if (PosLeg==2) PosLeg="topright"
  if (PosLeg==3) PosLeg="bottomright"
  if (PosLeg==4) PosLeg="bottomleft"
  
  BoxLeg = ifelse(BoxLeg=="S","o","n") # Frame in the captions, "n" without frame, "o" with frame
  
  Color_a = ifelse(Color=="S","red","black") # Colors at the points of the graphics
  Color_b = cor # For subtitles letters and their representations in the graphic
  if (Color=="S") Color_b = (cor+1):(cor+NumGrupos)
  #####   end - Information Used in Graphics  #####
  
  ##### begin - Plotting Eigenvalues #####
  mp <- barplot(MFA$MatrixA[,1],names.arg=paste(round(MFA$MatrixA[,2],2),"%",sep=""),main = "Eigenvalue")
  ##### end - Plotting Eigenvalues #####
  
  ##### begin - Plotting the Global Analysis #####
  plot(MFA$MatrixF, # cria grafico para as coordenadas principais da Analise Global
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[1],  # Titulo
       asp = 2,           # Aspecto do Grafico
       pch = 15,          # Formato dos pontos 
       cex=1,             # Tamanho dos pontos
       xlim=c(min(MFA$MatrixF[,1])-0.1,max(MFA$MatrixF[,1])+0.1), # Dimensao para as linhas do grafico
       ylim=c(min(MFA$MatrixF[,2]-0.1),max(MFA$MatrixF[,2])+0.1), # Dimensao para as colunas do grafico
       col = ifelse(Color=="S","red","black"))       # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  #LocLab(MFA$MatrixF[,1:2], NomeLinhas)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  text(MFA$MatrixF, cex=1, NomeLinhas, pos=3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  ##### end - Plotting the Global Analysis #####
  
  ##### begin - Plotting the Analysis by Group Together with the Global Analysis #####
  ## begin - Find the maximum and minimum dimensions for the columns and rows ##
  MLC <- MFA$MatrixF[,1:2]
  for (i in 1:length(MFA$MatrixEFG)) 
    MLC <- rbind(MLC,MFA$MatrixEFG[[i]][,1:2])
  maxX = max(MLC[,1]) # Dimenssoes maximas das linhas do grafico
  minX = min(MLC[,1]) # Dimenssoes minimas das linhas do grafico
  maxY = max(MLC[,2]) # Dimenssoes maximas das colunas do grafico
  minY = min(MLC[,2]) # Dimenssoes minimas das colunas do grafico
  ## end - Find the maximum and minimum dimensions for the columns and rows ##
  
  plot(MFA$MatrixF, # Creates graphs for the main coordinates of Group Analysis
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[2], # Titulo
       asp = 1,           # Aspecto do grafico
       pch = 15,          # Formato dos pontos 
       cex=1.2,           # Tamanho dos pontos
       xlim=c(minX,maxX), # Dimensao para as linhas do grafico
       ylim=c(minY,maxY), # Dimensao para as colunas do grafico
       col = Color_a)     # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  #LocLab(MFA$MatrixF[,1:2], NomeLinhas)  # Coloca os nomes dos pontos das coordenadas principais da analise global
  text(MFA$MatrixF, cex=1,NomeLinhas, pos=3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais da analise global
  ## Acrescenta no grafico da Analise Global as coordenadas principais da Analise por Grupo
  NumObserv = 4 # numero de centroides a considerar para plotagem das orbitas
  NumLinhas = nrow(MFA$MatrixEFG[[1]]) # numero de linhas
  if (NumObserv<NumLinhas) {
    Position = floor(NumLinhas/NumObserv)
    Observ = as.vector(c(rep(1,NumObserv))) # observacoes a serem plotadas orbitando os centroides
    for (i in 1:(length(Observ)-2)) {
      Observ[i+1] = Position*i
    }     
    Observ[length(Observ)]=NumLinhas # observacoes a serem plotadas orbitando os centroides
  }
  
  if (NumObserv>=NumLinhas)
    Observ = 1:NumLinhas  # observacoes a serem plotadas orbitando os centroides
  
  for (i in 1:length(MFA$MatrixEFG)) {
    if (NamArr=="N") 
      points(MFA$MatrixEFG[[i]][Observ,1:2], pch = (2 + ifelse(Color=="S",i,0)), cex = 1.2, col = 1 + ifelse(Color=="S",i,0)) # adiciona ao grafico as coordenadas principais dos Grupos
    else
      #LocLab(MFA$MatrixEFG[[i]][Observ,1:2],NomeGrupos[i], col = 1 + ifelse(Color=="S",i,0)) # Coloca os nomes dos pontos das coordenadas principais dos Grupos
      text(MFA$MatrixEFG[[i]][Observ,1:2], pos=3, cex=1, NomeGrupos[i], col = 1 + ifelse(Color=="S",i,0),xpd = TRUE) # Coloca os nomes dos pontos das coordenadas principais dos Grupos
  }
  
  ## liga os pontos de cada Analise Global com cada ponto da Analise por Grupo
  for (j in 1:length(MFA$MatrixEFG)) 
    segments(MFA$MatrixF[Observ,1], MFA$MatrixF[Observ,2], MFA$MatrixEFG[[j]][Observ,1], MFA$MatrixEFG[[j]][Observ,2], lty = cor + j, col = ifelse(Color=="S",cor + j,cor), lwd=1.5)
  
  if (NamArr=="N")
    legend(PosLeg, NomeGrupos, lty = (cor+1):(cor+NumGrupos), col = Color_b, text.col = Color_b,
           bty=BoxLeg, text.font = 6, y.intersp = 0.8,xpd = TRUE) # cria a legenda
  ##### end - Plot Analysis by Group Together with Global Analysis #####
  
  ##### begin - Plotting Correlations of Major Components with Original Variables #####
  plot(0,0, # cria grafico para as coordenadas das Correlacoes dos Componentes Principais com as Variaveis Originais
       xlab = DescEixo1, # Nomeia Eixo X
       ylab = DescEixo2, # Nomeia Eixo Y
       main = "Correlation circle", # Titulo
       asp = 1,           # Aspecto do grafico
       cex=0,             # Tamanho dos pontos
       xlim=c(-1.1,1.1),  # Dimensao para as linhas do grafico
       ylim=c(-1.1,1.1))  # Dimensao para as colunas do grafico
  
  symbols(0, 0, circles = 1, inches = FALSE, fg = 1, add = TRUE) # cria um circulo
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  j  <- 1         # coluna inicial do Grupo de variaveis
  k  <- Grupos[1] # coluna final do Grupo de variaveis
  for (i in 1:NumGrupos) {  # foi necessario criar este for para poder colocar cores diferentes para cada Grupo de variaveis
    
    arrows(0,0,MFA$MatrixCCP[1,j:k],MFA$MatrixCCP[2,j:k], lty=i, code = 2, angle = 10, col = ifelse(Color=="S",cor + i,cor)) # cria a seta apontando para cada coordenada principal
    
    if (is.null(colnames(MFA$MatrixCCP[,j:k])))
      NomeVar<- paste("Comp.", 1:Grupos[i], sep = "") # Nomeia as colunas
    else
      NomeVar<- colnames(MFA$MatrixCCP[,j:k])
    
    #LocLab(t(MFA$MatrixCCP[,j:k]), NomeVar, col = ifelse(Color=="S",cor + i,cor)) # Coloca os nomes dos pontos das coordenadas principais
    text(t(MFA$MatrixCCP[,j:k]), cex=1, pos=3, NomeVar, col = ifelse(Color=="S",cor + i,cor), xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais
    
    j <- j + Grupos[i]  # coluna inicial do Grupo de variaveis
    
    k <- k + Grupos[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do Grupo de variaveis  
    
  }
  
  legend(PosLeg, NomeGrupos, lty = cor:(cor+NumGrupos), col = Color_b, text.col = Color_b,
         bty=BoxLeg, text.font = 6, y.intersp = 0.8,xpd = TRUE) # cria a legenda
  ##### end - Plotting Correlations of Major Components with Original Variables #####
  
  ##### begin - Plotting of Partial Inertia / Variant Scores #####
  VlrMinX = ifelse(min(MFA$MatrixEscVar[,1])>0,-0.01,min(MFA$MatrixEscVar[,1])) # Valor minimo para a linha X
  VlrMinY = ifelse(min(MFA$MatrixEscVar[,2])>0,-0.01,min(MFA$MatrixEscVar[,2])) # Valor minimo para a linha Y
  VlrMaxX = 1.01 # Valor maximo para a linha X
  VlrMaxY = 1.01 # Valor maximo para a linha Y
  plot(MFA$MatrixEscVar, # cria grafico para as coordenadas Inercias Parciais/Escores das Variareis
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[3], # Titulo
       asp = 1,           # Aspecto do grafico
       pch = 15,          # Formato dos pontos 
       cex=1,             # Tamanho dos pontos
       xlim=c(VlrMinX,VlrMaxX), # Dimensao para as linhas do grafico
       ylim=c(VlrMinY,VlrMaxY), # Dimensao para as colunas do grafico
       col = Color_a)       # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  #LocLab(MFA$MatrixEscVar[,1:2],rownames(MFA$MatrixEscVar))  # Coloca os nomes dos pontos das coordenadas principais das linhas
  text(MFA$MatrixEscVar,cex=1, rownames(MFA$MatrixEscVar), pos=3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  ##### END - Plotting Partial Inertia / Variance Scores #####
}