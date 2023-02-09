#' Creates the neighbourhood structure of a linear network
#' 
#' Given a linear network structure, this function creates the neighbourhood matrix ("queen" criterion) associated to it. Two segments of the network are neighbours if they share a vertex   
#' 
#' @param network - A \code{linnet} object representing a linear network structure
#' @return Returns a \code{listw} object in \code{"W"} style
#' @examples 
#' library(DRHotNet)
#' library(spatstat.geom)
#' library(spatstat.linnet)
#' library(spdep)
#' library(raster)
#' chicago_neighbourhood <- NeighbourhoodMatrixNetwork(chicago$domain)
#' class(chicago_neighbourhood)
#' chicago_neighbourhood$neighbours[[1]]
#' @export
NeighbourhoodMatrixNetwork <- function(network){
  aux=SpatialPolygons(lapply(1:network$lines$n, 
                             function(i) Polygons(list(Polygon(cbind(t(network$lines[[1]][i,c(1,3)]),
                                                                     t(network$lines[[1]][i,c(2,4)])))), 
                                                  paste0("Line",i))))
  queen=poly2nb(aux, queen=TRUE)
  W=nb2listw(queen, style="W", zero.policy=TRUE)
  return(W)
}

