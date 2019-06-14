#' Identifies differential risk hotspots along a linear network given a vector of relative probabilities computed over the middle points of the segments of the network
#' 
#' Given the relative probabilities of an event's occurrence along the middle points of a linear network, this function filters and groups in hotspots those segments satisfying two conditions: 1) exceeding the average of the relative probabilities for all the events in the network in \code{k} times the standard deviation of the set of probabilities, and 2) having provided \code{n} or more events of the network for the computation of their corresponding relative probability (a factor that depends on the choice of \code{sigma} when using the function \code{RelativeProbabilityNetwork}). In summary, \code{k} and \code{n} control the formation of differential risk hotspots along the network, given a set of relative probabilities covering the network. The choice of a higher value for \code{k} or \code{n} (or both) is more demanding and leads to a lower number of differential risk hotspots being detected. Users should test several values of \code{k} and \code{n} (sensitivity analysis on \code{k} and \code{n}) in order to reach reasonable choices for the research or practical purposes of their data analysis 
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param W - A \code{listw} object representing the neighbourhood structure that the road segments of the network form, possibly generated through the function \code{NeighbourhoodMatrixNetwork}
#' @param relative_probabilities - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{RelativityProbabilityNetwork}
#' @param k - A numeric value that controls the procedure of detecting differential risk hotspots (departure from average relative probability), as described above
#' @param n - A numeric value that controls the procedure of detecting differential risk hotspots (minimum size for the sample of events implicated in the computation of the relative probabilities), as described above
#' @return Returns a list that contains the differential risk hotspots found for \code{X}, given the value of \code{relative_probabilities}
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
#' k=1
#' n=20
#' hotspots=DiffHotspots_n_k(X_lix,W,relative_probabilities,k,n)
#' }
#' @export
DiffHotspots_n_k <- function(X,W,relative_probabilities,k,n){
  input=relative_probabilities
  network_lix=X$domain
  events_radius_sigma=c()
  for (e in c(1:network_lix$lines$n)){
    events_radius_sigma=c(events_radius_sigma,as.numeric(input[[e]][[1]][3]))
  }
  input_prob=c()
  for (i in c(1:length(input))){
    input_prob=c(input_prob,input[[i]][[1]][1])
  }
  ### sig contains the segments showing a high relative probability, for a large enough sample
  sig=which(((input_prob>=mean(input_prob)+k*sd(input_prob)) & (events_radius_sigma>=n))==T)
  if (length(sig)>0){
    ordered_segments=c(1:length(sig))
    ordered_segments=c(1:length(sig))
    names(ordered_segments)=sig
    hotspot=rep(0,length(sig))
    count_zeros=length(which((hotspot==0)==T))
    current_hotspot=1
    while (count_zeros>0){
      ### new hotspot contains the first segment that has not been assigned to a hotspot yet
      start=which((hotspot==0)==T)[1]
      hotspot[start]=current_hotspot
      ### neighbouring segments are searched
      aux=as.numeric(names(ordered_segments)[start])
      aux_old=c()
      while (length(aux)>0){
        neigh=c()
        for (j in c(1:length(aux))){
          neigh=c(neigh,W$neighbours[[aux[j]]])
        }
        neigh=unique(neigh)
        if (length(aux_old>0)){
          find_old=c()
          for (k in c(1:length(aux_old))){
            find_old=c(find_old,which((neigh==aux_old[k])==T))
          }
          if (length(find_old)){
            neigh=neigh[-find_old]
          }
        }
        if (length(aux)>0){
          find=c()
          for (k in c(1:length(aux))){
            find=c(find,which((neigh==aux[k])==T))
          }
          if (length(find)>0){
            neigh=neigh[-find]
          }
        }
        ### update aux_old
        aux_old=aux
        ### keep neighbouring segments that should be part of the hotspot (those also in sig)
        aux=intersect(neigh,sig)
        ### assign hotspot number
        for (m in c(1:length(aux))){
          hotspot[ordered_segments[toString(aux[m])]]=as.numeric(current_hotspot)
        }
      }
      ### increase hotspot number
      current_hotspot=current_hotspot+1
      ### find segments not assigned yet
      count_zeros=length(which((hotspot==0)==T))
    }
    ### rearrange result
    result=list()
    for (k in c(1:length(unique(hotspot)))){
      result[[k]]=as.numeric(names(ordered_segments)[hotspot==k])
    }
  } else{
    result=NULL
    message("No differential risk hotspots found for the parameters provided. Tune k or n (or both) and try again")
  }
  return(result)
}