
#' @export
svsample_fast_cpp4SD <- function(y, draws = 1, burnin = 0, priorspec = specify_priors(), startpara, startlatent, fast_sv = get_default_fast_sv()) {
  cat("A")
  # startpara <- digest_startpara(startpara, priorspec)
  startpara <- as.list(drop(startpara))  # this way, named vectors are also accepted
  if (is.null(startpara$latent0)) {
    startpara$latent0 <- startpara$mu
  }
  
  cat("B")
  result <- .Call(`_stochvol_svsample_fast_cpp4SD`, y, draws, burnin, priorspec, startpara, startlatent, fast_sv ) #, PACKAGE = "stochvol")
  cat("C")
  if (fast_sv$store_indicators) {
    fast_sv$init_indicators <- result$indicators[NROW(result$indicators), , drop = TRUE]
  }
  result$fast_sv <- fast_sv
  result
}
