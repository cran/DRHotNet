#' Plots an object obtained with \code{DiffHotspots_n_k}
#' 
#' This function plots a set of differential risk hotspots located along a linear network. An extension of the hotspots (including the kth order neighbours of the segments of the hotspots) is also plotted
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param hotspots - A set of differential risk hotspots obtained with the function \code{DiffHotspots_n_k}
#' @param order_extension - A natural number indicating a neighbourhood order to be used for constructing an extension of the differential risk hotspots. The summary is also given for the segments forming this extension 
#' @param which.plot - A numeric vector indicating which differential risk hotspots to plot (according to the way they are ordered in \code{hotspots})
#' @param eps_image - If set to \code{TRUE}, an .eps image is generated. By default it is set to \code{FALSE}
#' @examples 
#' library(DRHotNet)
#' library(spatstat.geom)
#' library(spatstat.linnet)
#' library(spdep)
#' library(raster)
#' rel_assault <- relpnet(X = chicago, 
#' lixel_length = 50, h = 50, mark = "marks", category_mark = "assault")
#' hotspots_assault <- drhot(X = chicago, rel_probs = rel_assault, 
#' k = 0.5, n = 4)
#' plothot(X = chicago, hotspots = hotspots_assault)
#' @export
plothot <- function(X, hotspots, order_extension = NULL, which.plot = NULL, eps_image=F){
  
  network=X$domain
  lixel_length=hotspots$lixel_length
  h=hotspots$h
  mark=hotspots$mark
  category_mark=hotspots$category_mark
  k=hotspots$k
  n=hotspots$n
  
  if (is.null(order_extension)){
    order_extension=round(h/lixel_length)
  }
  
  if (hotspots$lixel_length!=F){
    network=spatstat.linnet::lixellate(network,eps=hotspots$lixel_length)
    # project into the lixellized network
    X_aux=spatstat.linnet::lpp(cbind(X$data$x,X$data$y),network)
    spatstat.geom::marks(X_aux)=spatstat.geom::marks(X)
    X=X_aux
  }
  network_lix=X$domain
  
  # Create psp object
  network_lix_psp=spatstat.geom::as.psp(X)
  
  # Extract hotspots segments
  
  # Filter if specified
  
  if (!is.null(which.plot)){
    segments_hotspots <- c()
    for (j in which.plot){
      segments_hotspots <- c(segments_hotspots, hotspots[[1]][[j]])
    }
  } else {
    segments_hotspots <- c()
    for (j in 1:length(hotspots[[1]])){
      segments_hotspots <- c(segments_hotspots, hotspots[[1]][[j]])
    }
  }
  
  # Neighbourhood matrix and hotspots extension
  
  W=NeighbourhoodMatrixNetwork(network_lix)
  segments_hotspots_extension=KthOrderNeighbours(segments_hotspots,W,order=order_extension)
  
  # Plot
  
  if (eps_image){
    setEPS()
    postscript(paste0("diff_risk_hotspots_k_",gsub("\\.","_",toString(k)),"_n_",n,"_lixel_",lixel_length,
                      "_h_",h,"_type_",mark,"_",category_mark,".eps"), family="Helvetica")
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    plot(network_lix_psp, col="black", lwd=1,
         main=paste0("Differential risk hotspots '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", h = ",h, ",",
                     "\nk = ",k,", n = ",n), line=0)
    plot(network_lix_psp[segments_hotspots_extension,],
         add=T,col="#fc9272",lwd=3)
    plot(network_lix_psp[segments_hotspots,],
         add=T,col="#de2d26",lwd=3)
    dev.off()
  } else{
    par(xpd=TRUE)
    plot(network_lix_psp, col="black", lwd=1,
         main=paste0("Differential risk hotspots '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", h = ",h, ",",
                     "\nk = ",k,", n = ",n), line=0)
    plot(network_lix_psp[segments_hotspots_extension,],
         add=T,col="#fc9272",lwd=3)
    plot(network_lix_psp[segments_hotspots,],
         add=T,col="#de2d26",lwd=3)
  }
  

}
