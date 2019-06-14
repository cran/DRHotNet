#' Computes the probability of observing a type of event in a linear network at a set of points provided
#' 
#' Given a marked point pattern lying on a linear network structure, this function uses kernel density estimation (KDE) to estimate a probability of occurrence for a type of event specified by the user through the marks of the pattern at the points of the network provided. The marks of a point pattern represent additional information of the events that are part of the pattern  
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param W - A \code{listw} object representing the neighbourhood structure that the road segments of the network form
#' @param sigma - A numeric value representing the bandwidth parameter (in meters)
#' @param x - A numeric vector representing the x-coordinates of the points of the network at which the relative probability is computed
#' @param y - A numeric vector representing the y-coordinates of the points of the network at which the relative probability is computed
#' @param mark - Mark of \code{X} that is used to characterize the type of event. The algorithm searches microzones of the network where this mark is over- or underrepresented
#' @param category_mark - A numeric/character value from the set allowed in the chosen \code{mark} to compute the relative probability in relation to it
#' @return Returns a list of length equal to \code{length(x)}. Each component of the list is itself a list, containing two numeric vectors. Given an index \code{i}, the first vector provides the relative probability of the event type being considered at the point of coordinates \code{x[i]} and \code{y[i]}, the number of events of the type of interest found at \code{sigma} meters from the point of coordinates \code{x[i]} and \code{y[i]}, and the total number of events found at \code{sigma} meters from this point (regardless of the type). The second vector contains the indexes (according to object \code{X}) of all the events used for the estimation of the corresponding relative probability, at \code{sigma} meters of distance (regardless of the type) 
#' @examples 
#' library(DRHotNet)
#' library(spatstat)
#' library(spdep)
#' library(raster)
#' library(maptools)
#' \donttest{
#' network=SampleMarkedPattern$domain
#' X=SampleMarkedPattern
#' network_lix=lixellate(network,eps=50)
#' X_lix=lpp(cbind(X$data$x,X$data$y),network_lix)
#' marks(X_lix)=marks(X)
#' sigma=100
#' collision_type="Rear-end"
#' W=NeighbourhoodMatrixNetwork(network_lix)
#' projection="+proj=utm +zone=30 ellps=WGS84 +ellps=WGS84"
#' network_lix_sp=as.SpatialLines.psp(as.psp(network_lix))
#' proj4string(network_lix_sp)=projection
#' middle_points=SpatialLinesMidPoints(network_lix_sp)@coords
#' x=middle_points[,1]
#' y=middle_points[,2]
#' relative_probabilities=RelativeProbabilityNetwork(X_lix,W,sigma,x,y, 
#'                        "Collision",collision_type)
#' }
#' @references Baddeley, A., Rubak, E., & Turner, R. (2015). Spatial point patterns: methodology and applications with R. Chapman and Hall/CRC.
#' @references Diggle, P. J. (2013). Statistical analysis of spatial and spatio-temporal point patterns. Chapman and Hall/CRC.
#' @references Kelsall, J. E., & Diggle, P. J. (1995). Kernel estimation of relative risk. Bernoulli, 1(1-2), 3-16.
#' @references McSwiggan, G., Baddeley, A., & Nair, G. (2017). Kernel density estimation on a linear network. Scandinavian Journal of Statistics, 44(2), 324-345.
#' @export
RelativeProbabilityNetwork <- function(X, W, sigma, x, y, mark, category_mark){
  
  network=X$domain
  network_sp=as.SpatialLines.psp(as.psp(network))
  vertices_segments=cbind(network$from,network$to)
  
  ### find the road segment to which each event belongs to
  segment_location=X$data$seg
  
  ### lengths segments
  length_segments=SpatialLinesLengths(network_sp)
  
  ### define kernel function (Gaussian)
  kernel=function(x) (1/(sigma*sqrt(pi)))*exp(1)^(-(x/sigma)^2)
  
  ### initialize a list to save the results
  results=list()
  ### monitorize progress
  total <- length(x)
  ### create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for (s in c(1:length(x))){
    
    ### increase progress bar
    setTxtProgressBar(pb, s)
    
    ### locate (x[s],y[s]) in the network (segment)
    i=as.numeric(DetectarEjeGrafo(network,list(x=as.numeric(x[s]),y=as.numeric(y[s]))))
    vertices_segment=vertices_segments[i,]
    
    ### find position (column index) in marks(X) for the mark of interest
    
    marksX=data.frame(marks(X))
    if (!is.null(mark)){
      find_mark=which((names(marks(X))==mark)==T)
    }
    
    ### several sums are initialized
    sum_num=0
    sum_den=0
    condition=0
    numevents=0
    all_events=c()
    numevents_num=0
    numevents_den=0
    
    ### distance from (x,y) to the extremes of the road segment they belong to
    d1=pointDistance(c(x[s],y[s]),c(network$lines$ends$x0[i],network$lines$ends$y0[i]),lonlat = F)
    d2=pointDistance(c(x[s],y[s]),c(network$lines$ends$x1[i],network$lines$ends$y1[i]),lonlat = F)
    
    k=0
    
    ### events from the own road segment are searched and added
    events=which((segment_location==i)==T)
    if (length(events)>0){
      for (e in events){
        event_distance=pointDistance(c(x[s],y[s]),
                                     c(X$data$x[e],X$data$y[e]),lonlat = F)
        ### if distance < sigma and mark conditions are fulfilled, add
        if (event_distance<sigma){
          numevents=numevents+1
          all_events=rbind(all_events,c(e,i,k,event_distance,kernel(event_distance)))
          numevents_den=numevents_den+1
          sum_den=sum_den+kernel(event_distance)
          check=marksX[e,find_mark]==category_mark
          if (length(which(check==T)==T)==length(check)){
            numevents_num=numevents_num+1
            sum_num=sum_num+kernel(event_distance)
          }
        }
      }
    }
    
    ### if d1 and d2 > sigma, no search of events is performed
    if (d1>=sigma & d2>=sigma){
      condition=1
    } else {
      k=1
      counter=0
      while (condition==0){
        counter=counter+1
        select=c()
        if (k==1){
          neighbours=VecinosOrdenk(network,k,i,vertices_segments,W$neighbours)
          if (ncol(neighbours)>0){
            select=c(1:ncol(neighbours))
            lengths=c()
            for (v in neighbours[1,]){
              find=c(which((vertices_segments[v,]==vertices_segment[1])==T),
                     which((vertices_segments[v,]==vertices_segment[2])==T))
              if (vertices_segments[v,find]==vertices_segment[1]){lengths=c(lengths,d1)}
              if (vertices_segments[v,find]==vertices_segment[2]){lengths=c(lengths,d2)}
            }
            names(lengths)=c(neighbours[3,select])
          } else {
            condition=1
          }
        } 
        if (k>=2){
          neighbours=VecinosOrdenkSiguiente(network,k,i,vertices_segments,W$neighbours,neighbours)
          select=which((neighbours[2,]==k)==T)
          ### new lengthts, from previous neighbours
          lengths=lengths[as.character(neighbours[3,select])]
        }
        
        if (length(select)>0){
          ### find events in neighbours and collect
          for (j in select){
            ### neighbours[1,j] contains neighbour index
            events=which((segment_location==neighbours[1,j])==T)
            if (length(events)>0){
              for (e in events){
                if (PuntoComun(vertices_segments,neighbours[3,j],neighbours[1,j])==1){
                  event_distance=lengths[toString(neighbours[3,j])]+
                    pointDistance(c(network$lines$ends$x0[neighbours[1,j]],network$lines$ends$y0[neighbours[1,j]]),
                                  c(X$data$x[e],X$data$y[e]),lonlat = F)
                } else {
                  event_distance=lengths[toString(neighbours[3,j])]+
                    pointDistance(c(network$lines$ends$x1[neighbours[1,j]],network$lines$ends$y1[neighbours[1,j]]),
                                  c(X$data$x[e],X$data$y[e]),lonlat = F)
                }
                ### if event_distance<sigma, add
                if (event_distance<sigma){
                  numevents=numevents+1
                  ### product vertices multiplicity
                  vertices=CaminoVertices(neighbours,j)
                  grados=vertexdegree(network)[vertices]
                  factor=1
                  for (g in grados){
                    factor=factor*(2/g)
                  }
                  all_events=rbind(all_events,c(e,i,k,event_distance,kernel(event_distance)*factor))
                  numevents_den=numevents_den+1
                  sum_den=sum_den+kernel(event_distance)*factor
                  check=marksX[e,find_mark]==category_mark
                  if (length(which(check==T)==T)==length(check)){
                    numevents_num=numevents_num+1
                    sum_num=sum_num+kernel(event_distance)*factor
                  }
                }
              }
            }
          }
          
          if (k==1){
            if (ncol(neighbours)>0){
              lengths=lengths+length_segments[neighbours[1,]]
              names(lengths)=neighbours[1,]
            }
          }
          if (k>=2){
            lengths=lengths+length_segments[neighbours[1,select]]
            names(lengths)=neighbours[1,select]
          }
          filter=which((lengths<=sigma)==F)
          lengths[filter]=Inf
          if (length(which((lengths==Inf)==T))==length(lengths)){
            condition=1
          }
          k=k+1
        } else {
          condition=1
        }
      }
    }
    ### put results of interest for point s into a list
    if (numevents>0){
      result=list()
      result[[1]]=round(c(sum_num/sum_den,numevents_num,numevents_den),4)
      result[[2]]=as.numeric(all_events[,1])
    } else{
      result=list()
      result[[1]]=c(0,0,0)
      result[[2]]=c(0)
    }
    ### storage result for point s in the general list
    names(result[[1]])=c("Relative probability","No. events category mark","No. all events")
    names(result[[2]])=NULL
    results[[s]]=result
  }
  return(results)
}
