# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

arma_mahalanobis <- function(x, center, cov) {
    .Call(`_rsdm_arma_mahalanobis`, x, center, cov)
}

arma_calc_nt1 <- function(trg, ref) {
    .Call(`_rsdm_arma_calc_nt1`, trg, ref)
}

arma_calc_mic_nt1 <- function(trg, ref) {
    .Call(`_rsdm_arma_calc_mic_nt1`, trg, ref)
}

arma_calc_nt2 <- function(trg, ref) {
    .Call(`_rsdm_arma_calc_nt2`, trg, ref)
}

arma_calc_mic_nt2 <- function(trg, ref) {
    .Call(`_rsdm_arma_calc_mic_nt2`, trg, ref)
}

arma_micdet <- function(trg, ref) {
    .Call(`_rsdm_arma_micdet`, trg, ref)
}

arma_micdet_raster <- function(trg, ref) {
    .Call(`_rsdm_arma_micdet_raster`, trg, ref)
}

arma_exdet <- function(trg, ref) {
    .Call(`_rsdm_arma_exdet`, trg, ref)
}

arma_exdet_raster <- function(trg, ref, compute_mic = FALSE) {
    .Call(`_rsdm_arma_exdet_raster`, trg, ref, compute_mic)
}

