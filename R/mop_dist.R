#' MOP distance
#' @description mop_dist computes a ponderated Euclidean distance for MOP.
#'
#' @param eudist_matrix a matrix of Euclidean distance between M and G areas.
#' @param irow row from eudist_matrix where the pondarated distance will be compute
#' @param percent (numeric) percent of values sampled from te calibration region to calculate the MOP.
#' @return a vector of poderated distances for irow


mop_dist <-  function(eudist_matrix,irow,percent) {
  di <- eudist_matrix[irow, ]
  qdi <- quantile(di, probs = percent / 100, na.rm = TRUE)
  ii <-  which(di <= qdi)
  return(mean(di[ii]))
}
