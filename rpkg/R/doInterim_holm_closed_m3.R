#' Interim: conditional errors for Holm closed testing (m=3)
#' @description Interim: conditional errors for Holm closed testing (m=3)
#'
#' @param z1 A numeric vector giving first stage z-values computed from first stage p-values.
#' @param v A numeric vector giving the proportions of pre-planned measurements collected up to the interim analysis.
#' @param alpha significance level
#' @export A1 matrix of pCERs and CER  for intersections {1},{2},{3},{12},{13},{23},{123}; same as with gMCP doInterim;
#' @details eWHORM simulations
#' @author Sonja Zehetmayer
#' 

doInterim_holm_closed_m3 <- function(z1, v, alpha = 0.025) {
  # PCERs by intersection size
  pcer_k1 <- pcer_bonf_k(z1=z1,v=v,alpha=alpha,k=1) # for {i}
  pcer_k2 <- pcer_bonf_k(z1=z1,v=v,alpha=alpha,k=2) # used in pairs
  pcer_k3 <- pcer_bonf_k(z1=z1,v=v,alpha=alpha,k=3) # used in triple
  names(pcer_k1) <- names(pcer_k2) <- names(pcer_k3) <- paste0("H", 1:3)
  
  A<-matrix(rbind(c(0,0,pcer_k1[3]),#H3
                  c(0,pcer_k1[2],0),#H2
                  c(0,pcer_k2[2:3]),#(c(0,pcer_k2[2:3,2]),#H23
                  c(pcer_k1[1],0,0),#H1
                  c(pcer_k2[1],0,pcer_k2[3]),#c(pcer_k2[1,1],0,pcer_k2[1,3]), #H13
                  c(pcer_k2[1:2],0), #c(pcer_k2[2,1],pcer_k2[2,2],0), #H12
                  c(pcer_k3)),ncol=3)  #H123
  
  A1<-cbind(A,apply(A,1,sum))
  
  A1
}
