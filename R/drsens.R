#' Performs a sensitivity analysis on the parameters \code{k} and \code{n} that are provided to \code{drhot}
#' 
#' Given a set of \code{ks} and \code{ns} parameters, this function allows the user to perform a sensitivity analysis on the parameters \code{k} and \code{n} by calling \code{drhot} for each combination of \code{k} and \code{n}
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param rel_probs - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{relpnet}
#' @param ks - A numeric vector of possible values for the \code{k} parameter that is provided to \code{drhot} 
#' @param ns - A numeric vector of possible values for the \code{n} parameter that is provided to \code{drhot} 
#' @return A matrix providing the type-specific prediction accuracy index that corresponds to the set differential risk hotspots obtained for each value of \code{k} or \code{n} provided in \code{ks} and \code{ns}, respectively. A \code{NA} value in this matrix indicates that no differential risk hotspots are found for the corresponding combination of \code{k} and \code{n}
#' @examples 
#' library(DRHotNet)
#' library(spatstat.geom)
#' library(spatstat.linnet)
#' library(spdep)
#' library(raster)
#' rel_assault <- relpnet(X = chicago, 
#' lixel_length = 50, h = 50, mark = "marks", category_mark = "assault")
#' sensitivity_analysis <- drsens(X = chicago, rel_probs = rel_assault, 
#' ks = c(1,2), ns = c(30,40))
#' @references Briz-Redon, A., Martinez-Ruiz, F., & Montes, F. (2019). Identification of differential risk hotspots for collision and vehicle type in a directed linear network. Accident Analysis & Prevention, 132, 105278.
#' @export
drsens <- function(X,rel_probs,ks,ns){
  
  network=X$domain
  
  if (rel_probs$lixel_length!=F){
    network=spatstat.linnet::lixellate(network,eps=rel_probs$lixel_length)
    #& project into the new network
    X_aux=spatstat.linnet::lpp(cbind(X$data$x,X$data$y),network)
    spatstat.geom::marks(X_aux)=spatstat.geom::marks(X)
    X=X_aux
  }
  network_lix=X$domain
  midpoints=spatstat.geom::midpoints.psp(spatstat.geom::as.psp(network_lix))
  segment_lengths=spatstat.geom::lengths_psp(spatstat.geom::as.psp(network_lix))
  
  # Midpoints as a point pattern (on the original network)
  lpp_midpoints=spatstat.linnet::lpp(midpoints,network)

  PAIs=matrix(NA,nrow=length(ks),ncol=length(ns))
  for (i in 1:length(ks)){
    # Compute distances between the middle points of those segments satisfying condition on k
    sig=which(rel_probs$probs>=mean(rel_probs$probs)+ks[i]*sd(rel_probs$probs))
    distances=spatstat.linnet::crossdist.lpp(lpp_midpoints[sig,],X)
    for (j in 1:length(ns)){
      cat(paste0("k = ",ks[i],", n = ",ns[j]),"\n")
      hotspots=drhot(X,rel_probs,ks[i],ns[j],event_distances=distances)
      if (!is.null(hotspots)){
        PAIs[i,j]=hotspots$PAI
      }
    }
  }
  rownames(PAIs)=paste0("k = ",ks)
  colnames(PAIs)=paste0("n = ",ns)
  return(PAIs)
}
