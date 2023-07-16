#' Plots an object obtained with \code{relpnet}
#' 
#' This function plots the relative probability of occurrence of a type of event along a linear network
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param rel_probs - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{relpnet}
#' @param eps_image - If set to \code{TRUE}, an .eps image is generated. By default it is set to \code{FALSE}
#' @examples 
#' library(DRHotNet)
#' library(spatstat.geom)
#' library(spatstat.linnet)
#' library(spdep)
#' rel_assault <- relpnet(X = chicago, 
#' lixel_length = 50, h = 50, mark = "marks", category_mark = "assault")
#' plotrelp(X = chicago, rel_probs = rel_assault)
#' @export
plotrelp <- function(X, rel_probs, eps_image=F){
  
  floor_dec=function(x, level=1) round(x - 5*10^(-level-1), level)
  ceiling_dec=function(x, level=1) round(x + 5*10^(-level-1), level)
  
  network=X$domain
  lixel_length=rel_probs$lixel_length
  h=rel_probs$h
  mark=rel_probs$mark
  category_mark=rel_probs$category_mark
  
  if (rel_probs$lixel_length!=F){
    network=spatstat.linnet::lixellate(network,eps=rel_probs$lixel_length)
    # project into the lixellized network
    X_aux=spatstat.linnet::lpp(cbind(X$data$x,X$data$y),network)
    spatstat.geom::marks(X_aux)=spatstat.geom::marks(X)
    X=X_aux
  }
  network_lix=X$domain
  
  # Create psp object
  network_lix_psp=spatstat.geom::as.psp(X)
  marks(network_lix_psp)=rel_probs$probs
  
  # Plot
  
  if (eps_image){
    setEPS()
    postscript(paste0("rel_probs_lixel_",lixel_length,"_h_",h,"_type_",mark,"_",category_mark,".eps"), family="Helvetica")
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    plot(network_lix_psp,lwd=2,
         main=paste0("Relative probabilities '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", h = ",h), line=0) 

    dev.off()
  } else{
    par(xpd=TRUE)
    plot(network_lix_psp,lwd=2,
         main=paste0("Relative probabilities '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", h = ",h), line=0) 
  }
  
  
}
