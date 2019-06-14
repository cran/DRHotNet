DetectarEjeGrafo <-function(grafo,datos){
  ejes=c()
  for (j in c(1:length(datos$x))){
    #print(j)
    x=datos$x[j]
    y=datos$y[j]
    aux=list(x=x,y=y)
    X_aux=lpp(aux,grafo)
    eje=X_aux$data$seg
    ejes=c(ejes,eje)
  }
  return(ejes)
}

PuntoComun <- function(vertices_ejes,i,j){
  vertices_eje=vertices_ejes[i,]
  buscar=c(which((vertices_ejes[j,]==vertices_eje[1])==T),
           which((vertices_ejes[j,]==vertices_eje[2])==T))
  return(buscar)
}

VerticeComun <- function(vertices_ejes,i,j){
  vertices_eje=vertices_ejes[i,]
  buscar=c(which((vertices_ejes[j,]==vertices_eje[1])==T),
           which((vertices_ejes[j,]==vertices_eje[2])==T))
  return(vertices_ejes[j,buscar])
}

CaminoVertices <- function(Vecinos, j){
  eje=Vecinos[1,j]
  ejes_camino=c()
  vertices_camino=c()
  fin=as.numeric(Vecinos[2,j])-1
  for (k in c(1:fin)){
    ejes_camino=c(ejes_camino,Vecinos[3,which((Vecinos[1,]==eje)==T)])
    vertices_camino=c(vertices_camino,Vecinos[4,which((Vecinos[1,]==eje)==T)])
    eje=Vecinos[3,which((Vecinos[1,]==eje)==T)]
  }
  return(vertices_camino)
}

VecinosOrdenk <- function(grafo, k, i, vertices_ejes, lista_vecinos){
  resultado=c()
  for (j in c(1:grafo$lines$n)){
    lista_vecinos[[j]]=rbind(lista_vecinos[[j]],rep(1,length(lista_vecinos[[j]])))
  }
  resultado=rbind(lista_vecinos[[i]],rep(i,ncol(lista_vecinos[[i]])))
  vertices=c()
  for (v in resultado[1,]){
    vertices=c(vertices,VerticeComun(vertices_ejes,i,v))
  }
  resultado=rbind(resultado,vertices)
  if (k>1){
    for (j in c(2:k)){
      vecinos_actuales=resultado
      for (s in c(1:ncol(vecinos_actuales))){
        vecino=vecinos_actuales[1,s]
        if (length(lista_vecinos[[vecino]][1,])>1){
          buscar_i=which((lista_vecinos[[vecino]][1,]==i)==F)
        } else {
          buscar_i=1
        }
        vertices_comunes=c()
        for (l in c(1:length(buscar_i))){
          vertices_comunes=c(vertices_comunes,VerticeComun(vertices_ejes,vecino,lista_vecinos[[vecino]][1,buscar_i[l]]))
        }
        fila_aux=rep(NA,ncol(lista_vecinos[[vecino]]))
        fila_aux[buscar_i]=vertices_comunes
        fila_aux[-buscar_i]=i
        resultado=cbind(resultado,
                        rbind(lista_vecinos[[vecino]][1,],
                              rep(j,length(lista_vecinos[[vecino]][1,])),
                              rep(vecino,length(lista_vecinos[[vecino]][1,])),
                              fila_aux))
      }
    }
    if (ncol(resultado)>0){
      buscar_i=which((resultado[1,]==i)==T)
      resultado=resultado[,-buscar_i]
      if (class(resultado)=="numeric"){
        resultado=as.matrix(resultado)
      }
      buscar_duplicados=which(duplicated(resultado[1,])==T)
      if (length(buscar_duplicados)>0){
        resultado=resultado[,-buscar_duplicados]
      }
    }
  }
  rownames(resultado)=c("Neighbour","Order","Previous","UnionVertex")
  colnames(resultado)=rep("",ncol(resultado))
  return(resultado)
}


VecinosOrdenkSiguiente <- function(grafo, k, i, vertices_ejes, lista_vecinos, resultado){
  for (j in c(1:grafo$lines$n)){
    lista_vecinos[[j]]=rbind(lista_vecinos[[j]],rep(1,length(lista_vecinos[[j]])))
  }
  for (j in c(k)){
    vecinos_actuales=resultado
    for (s in c(1:ncol(vecinos_actuales))){
      vecino=vecinos_actuales[1,s]
      if (length(lista_vecinos[[vecino]][1,])>1){
        buscar_i=which((lista_vecinos[[vecino]][1,]==i)==F)
      } else {
        buscar_i=1
      }
      if (ncol(lista_vecinos[[vecino]])>0){
        vertices_comunes=c()
        for (l in c(1:length(buscar_i))){
          vertices_comunes=c(vertices_comunes,VerticeComun(vertices_ejes,vecino,lista_vecinos[[vecino]][1,buscar_i[l]]))
        }
        fila_aux=rep(NA,ncol(lista_vecinos[[vecino]]))
        fila_aux[buscar_i]=vertices_comunes
        fila_aux[-buscar_i]=i
        resultado=cbind(resultado,
                        rbind(lista_vecinos[[vecino]][1,],
                              rep(j,length(lista_vecinos[[vecino]][1,])),
                              rep(vecino,length(lista_vecinos[[vecino]][1,])),
                              fila_aux))
      }
    }
  }
  if (ncol(resultado)>0){
    buscar_i=which((resultado[1,]==i)==T)
    resultado=resultado[,-buscar_i]
    if (class(resultado)=="numeric"){
      resultado=as.matrix(resultado)
    }
    buscar_duplicados=which(duplicated(resultado[1,])==T)
    if (length(buscar_duplicados)>0){
      resultado=resultado[,-buscar_duplicados]
    }
  }
  rownames(resultado)=c("Neighbour","Order","Previous","UnionVertex")
  colnames(resultado)=rep("",ncol(resultado))
  return(resultado)
}