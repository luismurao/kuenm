#' Extrapolation risk analysis for single comparisons
#'
#' @description kuenm_mop calculates a mobility-oriented parity layer by
#' comparing environmental values between the calibration area and the area or
#' scenario to which an ecological niche model is transferred.
#'
#' @param M.stack a RasterStack of variables representing the calibration area.
#' @param G.stack a RasterStack of variables representing the full area of interest, and areas
#' or scenarios to which models are transferred.
#' @param percent (numeric) percent of values sampled from te calibration region to calculate the MOP.
#' @param comp_each (numeric) compute distance matrix for a each fixed number of rows (default 1000).
#' @param normalized (logical) if true values of similarity are presented from 0 to 1,
#' default = TRUE.
#'
#' @return A mobility-oriented parity RasterLayer.
#'
#' @details The MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
#' of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#' @examples
#' mvars <- mvars_mop
#' gvars <- gvars_mop
#' perc <- 10
#' norm <- TRUE
#'
#' mop <- kuenm_mop(M.stack = mvars, G.stack = gvars,
#'                   percent = perc, normalized = norm)

kuenm_mop <- function(M.stack, G.stack, percent = 10,
                      comp_each=2000,normalized = TRUE) {


  mPoints <- raster::rasterToPoints(M.stack)
  m_nona <- na.omit(mPoints)
  m_naID <- attr(m_nona,"na.action")
  gPoints <- raster::rasterToPoints(G.stack)
  g_nona <- na.omit(gPoints)
  g_naID <- attr(g_nona,"na.action")

  m1 <<- m_nona[, -(1:2)]
  m2 <<- g_nona[, -(1:2)]

  if(dim(m1)[2] != dim(m2)[2]) {
    stop("Stacks must have the same dimensions")
  }

  suppressPackageStartupMessages(library("future"))
  future::plan(multiprocess)
  mop_env <- new.env()

  steps <- seq(1, dim(m2)[1], comp_each)
  kkk <- c(steps,  dim(m2)[1] + 1)
  print(kkk)
  out_index <- kuenm::plot_out(m1, m2)
  long_k <- length(kkk)

  pasos <- 1:(length(kkk) - 1)
  pasosChar <- paste0(pasos)

  for (paso in pasosChar) {
    x <- as.numeric(paso)
    mop_env[[paso]] %<-% {
      seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
      eudist <- fields::rdist(m2[seq_rdist, ], m1)
      mop_dist <- lapply(1:dim(eudist)[1], function(y){
        di <- eudist[y, ]
        qdi <- quantile(di, probs = percent / 100,
                        na.rm = TRUE)
        ii <-  which(di <= qdi)
        pond_mean <- mean(di,na.rm = TRUE)
        return(pond_mean)
      })
      mop <-unlist(mop_dist)
      return(mop)
    }
    avance <- (x / long_k) * 100
    cat("Computation progress: ", avance,"%" ,"\n")
  }

  mop_list <- as.list(mop_env)
  mop_names <- sort(as.numeric(names(mop_list)))
  mop_names <- as.character(mop_names)
  mop_vals <- unlist(mop_list[mop_names])

  if(!is.null(g_naID)){
    mop_all <- data.frame(gPoints[,1:2])
    mop_all$mop <- NA
    mop_all$mop[-g_naID] <- mop_vals
  }
  else{
    mop_all <- data.frame(gPoints[, 1:2], mop=mop_vals)
  }

  mop_max <- max(na.omit( mop_all$mop)) * 1.05
  mop_all[out_index, 3] <- mop_max


  mop2 <- unlist(mop1)
  mop_all <- data.frame(gPoints[, 1:2], mop2)
  mop_max <- max(na.omit(mop2))
  mop_max <- max(mop2)
  mop_all[out_index, 3] <- mop_max*1.05
  sp::coordinates(mop_all) <- ~x + y
  suppressWarnings({
    sp::gridded(mop_all) <- TRUE
  })

  mop_raster <- raster::raster(mop_all)
  mop_raster <- 1 - (mop_raster / mop_max)
  return(mop_raster)
}
