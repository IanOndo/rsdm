#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR

#include <R_ext/Print.h>
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

using std::pow;
using std::exp;
using std::log;

 //' @title covariance function
 //'
 //' Computes the covariance of the elements of a matrix
 //' @param X A numeric matrix
 //' @param norm_type Optional integer number indicating the type of normalisation to be used.
 //'  Default is \code{norm_type = 0}, i.e. normalisation is done using N-1 with N denoting the number of observations.
 //' @return A numeric matrix of the same dimensions as the input matrix X.
 //' @url https://arma.sourceforge.net/docs.html#cov
 arma::mat arma_cov(const arma::mat& X, const int norm_type = 0) {

   arma::mat out(X.n_cols, X.n_cols);

   // Degenerate cases
   if (X.n_cols == 0) {
     return out;
   } else if (X.n_rows == 0 || X.n_rows == 1) {
     out.fill(Rcpp::NumericVector::get_na());
   } else {
     out = arma::cov(X, norm_type);
   }
   return out;
 }

 //' @title Mahalanobis distance matrix
 //'
 //' @description Computes the Mahalanobis distance between each rows of a matrix
 //' @param x A numeric matrix
 //' @param center A vector of means values used to center the observations.
 //' @param cov The covariance matrix of the observations
 //' @return A numeric matrix withe same dimensions as the input matrix x.
 // [[Rcpp::export]]
 arma::vec arma_mahalanobis(arma::mat const &x,
                            arma::vec const &center,
                            arma::mat const &cov) {
   arma::mat x_cen = x.t();
   x_cen.each_col() -= center;
   arma::solve(x_cen, arma::trimatl(arma::chol(cov).t()), x_cen);
   x_cen.for_each( [](arma::mat::elem_type& val) { val = val * val; } );
   return arma::sum(x_cen, 0).t();
 }

 //' @title NT1 - Novelty index of type 1
 //'
 //' @author Ian Ondo
 //'
 //' @description Given a target and a reference dataset, computes the univariate extrapolation index
 //' as defined in Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014.
 //' @param trg A N x M numerix matrix of N observations from M variables.
 //' @param ref A N'x M numeric matrix of N' observations from the same M variables as `trg` taken as a reference dataset.
 //' @return A numeric vector of (negative) novelty values of the same length as the number of rows
 //' of the target data matrix.
 //' @name nt1
 // [[Rcpp::export]]
 arma::vec arma_calc_nt1(arma::mat &trg, arma::mat &ref){

   arma::uword i, nrows_trg = trg.n_rows, ncols_trg=trg.n_cols;

   //---------------------------------------------
   // Min/Max for each covariate in reference system
   //---------------------------------------------
   arma::mat range_arr = arma::join_cols(arma::min(ref,0),arma::max(ref,0));
   arma::rowvec diffs_ref = arma::diff(range_arr,1,0);

   //---------------------------------------------
   // Use matrix algebra to calculate univariate extrapolation (NT1)
   //---------------------------------------------
   arma::cube iud = arma::zeros<arma::cube>(nrows_trg, ncols_trg, 3);
   iud.slice(1) = trg.each_row() - range_arr.row(0);
   iud.slice(2) = range_arr.row(1) - trg.each_row();

   arma::mat rvec, UDs(nrows_trg,ncols_trg);
   for(i=0; i<nrows_trg; i++){
     rvec = arma::min(arma::mat(iud.row(i)),1);
     UDs.row(i) = arma::rowvec(rvec.t()/diffs_ref);
   }

   //---------------------------------------------
   // Compute total NT1
   //---------------------------------------------
   return arma::sum(UDs,1);
 }

 //' @title MIC1 - Most influential covariate of type 1
 //'
 //' @author Ian Ondo
 //'
 //' @description Given a target dataset, identifies the most influential variable of type 1,
 //' i.e. the variable that extrapolates the most beyond the range of reference values.
 //' @param trg A N x M numerix matrix of N observations from M variables.
 //' @param ref A N'x M numeric matrix of N' observations from the same M variables as `trg` taken as a reference dataset.
 //' @return A numeric vector of integer values denoting which variable (i.e. which column) has the highest degree of extrapolation
 //'  beyond the reference conditions
 //'  @name mic1
 // [[Rcpp::export]]
 Rcpp::IntegerVector arma_calc_mic_nt1(arma::mat &trg, arma::mat &ref){

   arma::uword i, id_min, nrows_trg = trg.n_rows, ncols_trg=trg.n_cols;

   //---------------------------------------------
   // Min/Max for each covariate in reference system
   //---------------------------------------------
   arma::mat range_arr = arma::join_cols(arma::min(ref,0),arma::max(ref,0));
   arma::rowvec diffs_ref = arma::diff(range_arr,1,0);

   //---------------------------------------------
   // Use matrix algebra to calculate univariate extrapolation (NT1)
   //---------------------------------------------
   arma::cube iud = arma::zeros<arma::cube>(nrows_trg, ncols_trg, 3);
   iud.slice(1) = trg.each_row() - range_arr.row(0);
   iud.slice(2) = range_arr.row(1) - trg.each_row();

   arma::mat rvec, UDs(nrows_trg,ncols_trg);
   for(i=0; i<nrows_trg; i++){
     rvec = arma::min(arma::mat(iud.row(i)),1);
     UDs.row(i) = arma::rowvec(rvec.t()/diffs_ref);
   }

   //---------------------------------------------
   // Compute MIC NT1
   //---------------------------------------------
   UDs.replace(0,arma::datum::nan);
   Rcpp::IntegerVector mic_nt1_rcpp(nrows_trg, NA_INTEGER);
   if(UDs.elem(arma::find_nonfinite(UDs)).eval().n_elem==UDs.n_elem) return mic_nt1_rcpp;
   arma::ivec mic_nt1(mic_nt1_rcpp.begin(), mic_nt1_rcpp.length(), false);
   arma::uvec id_finite;
   arma::rowvec r;
   for(i=0; i<nrows_trg; i++){
     r = UDs.row(i);
     id_finite = arma::find_finite(r);
     if(!id_finite.is_empty()){
       id_min = r.elem(id_finite).eval().index_min();
       mic_nt1(i) = id_finite(id_min);
     }
   }
   mic_nt1+=1;

   mic_nt1_rcpp[mic_nt1_rcpp<0] = NA_INTEGER;

   return mic_nt1_rcpp;
 }

 //' @title NT2 - Novelty of type 2
 //'
 //' @author Ian Ondo
 //'
 //' @description Given a target and reference dataset, computes the multivariate or combinatorial extrapolation index
 //' as defined in Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014.
 //' @param trg A N x M numerix matrix of N observations from M variables.
 //' @param ref A N'x M numeric matrix of N' observations from the same M variables as `trg` taken as a reference dataset.
 //' @return A numeric vector of the type 2 novelty index of the same length as the number of rows
 //' of the target data matrix.
 //' @name nt2
 // [[Rcpp::export]]
 arma::vec arma_calc_nt2(arma::mat &trg, arma::mat &ref){

   // Calculate the average and covariance matrix of the variables
   // in the reference set
   arma::vec ref_av = arma::vectorise(arma::mean(ref,0));
   arma::mat ref_cov = arma_cov(ref);

   // Calculate the mahalanobis distance of each observation to the environmental center of the reference
   // set for both the reference and the projection dataset and calculate the ratio between the two.
   arma::vec maha_ref = arma_mahalanobis(ref, ref_av, ref_cov);
   arma::vec maha_trg = arma_mahalanobis(trg, ref_av, ref_cov);
   double maha_max = arma::max(maha_ref.elem(arma::find_finite(maha_ref)));

   //---------------------------------------------
   // Compute total NT2
   //---------------------------------------------
   return maha_trg/maha_max;
 }

 //' @title MIC2 - Most influential covariate of type 2
 //'
 //' @author Ian Ondo
 //'
 //' @description Given a target dataset, identifies the most influential variable of type 2
 //' i.e. the variable with the most dissimilar correlation structure compared to a reference dataset.
 //' @param trg A N x M numerix matrix of N observations from M variables.
 //' @param ref A N'x M numeric matrix of N' observations from the same M variables as `trg` taken as a reference dataset.
 //' @return A numeric vector of integer values denoting which variable (i.e. which column) has the highest degree of extrapolation
 //' beyond the correlation structure found in the reference dataset
 //' @name mic2
 // [[Rcpp::export]]
 arma::uvec arma_calc_mic_nt2(arma::mat &trg, arma::mat &ref){

   // Calculate the average and covariance matrix of the variables
   // in the reference set
   arma::vec ref_avg = arma::vectorise(arma::mean(ref,0));
   arma::mat ref_cov = arma_cov(ref);

   // Calculate the mahalanobis distance of each observation to the environmental center of the reference
   // set for both the reference and the projection dataset and calculate the ratio between the two.
   arma::vec maha_trg, maha_trg_all = arma_mahalanobis(trg, ref_avg, ref_cov);
   arma::mat IC(trg.n_rows, trg.n_cols);
   arma::uvec id_col, id_cols = arma::regspace<arma::uvec>(0,trg.n_cols-1);
   for(arma::uword j=0; j<trg.n_cols; ++j){
     id_col = arma::find(id_cols!=j);
     maha_trg = arma_mahalanobis(trg.cols(id_col).eval(), ref_avg.elem(id_col).eval(), ref_cov.submat(id_col,id_col).eval());
     IC.col(j) = 100. * ((maha_trg_all - maha_trg)/maha_trg_all);
   }

   return arma::index_max(IC,1)+1;
 }

 //' @title MIC - Most influential covariate
 //'
 //' @description Given a target dataset, identifies the most influential covariate
 //' i.e. the covariate with the most dissimilar values compared to a reference dataset.
 //' @param trg A N x M numerix matrix of N observations from M variables.
 //' @param ref A N'x M numeric matrix of N' observations from the same M variables as `trg` taken as a reference dataset.
 //' @return A numeric vector of integer values denoting which variable (i.e. which column) has the highest novelty value
 //' @export
 // [[Rcpp::export]]
 SEXP arma_micdet(arma::mat &trg, arma::mat &ref){

   // compute most influential covariate for novelty of type 1
   Rcpp::IntegerVector mic_NT1 = arma_calc_mic_nt1(trg, ref);

   // compute most influential covariate for novelty of type 2 on remaining points (i.e. inside the range where NT1=0)
   arma::ivec arma_mic_nt1(mic_NT1.begin(), mic_NT1.length(), false);
   arma::uvec is_in_range = arma::find_finite(arma_mic_nt1);

   if(is_in_range.is_empty()){
     Rcpp::warning("arma_micdet(): cannot compute mic for novelty of type 2, because no observation falls within the range of the reference data.");
     return Rcpp::wrap(mic_NT1);
   }
   arma::mat trg_in_range = trg.rows(is_in_range);
   arma::uvec mic_NT2 = arma_calc_mic_nt2(trg_in_range, ref);
   arma_mic_nt1.elem(is_in_range) = arma::conv_to<arma::ivec>::from(mic_NT2);

   return Rcpp::wrap(mic_NT1);
 }

 //' @title MIC - Most influential covariate for rasters
 //'
 //' @author Ian Ondo
 //'
 //' @description Given a target raster, identifies the most influential variable
 //' i.e. the variable with the most dissimilar values compared to a reference raster.
 //' @param trg A RasterBrick object with M variables.
 //' @param ref A RasterBrick object with the same M variables as `trg` taken as a reference.
 //' @return A RasterLayer object with values denoting which variable (i.e. layer index) has the highest novelty value
 //' @export
 // [[Rcpp::export]]
 Rcpp::S4 arma_micdet_raster(Rcpp::S4 &trg, Rcpp::S4 &ref){

   Rcpp::S4 raster_out = Rcpp::clone(trg);
   Rcpp::S4 refdata(ref.slot("data"));
   Rcpp::S4 trgdata(raster_out.slot("data")), trgfile(raster_out.slot("file"));

   // compare number of raster layers
   arma::uword ref_nlayers = refdata.slot("nlayers"), trg_nlayers = refdata.slot("nlayers");
   if(ref_nlayers!=trg_nlayers) Rcpp::stop("arma_exdet_raster(): the reference and target rasters must have the same number of layers.");

   // compare raster layers names
   CharacterVector ref_names = refdata.slot("names"), trg_names = trgdata.slot("names");
   CharacterVector::iterator refl = ref_names.begin(), trgl=trg_names.begin();

   for(; refl != ref_names.end(); ++refl, ++trgl) if(*refl!=*trgl) Rcpp::stop("arma_micdet_raster(): layers names are not identical.");

   // get target raster info
   arma::uword ref_nrows=ref.slot("nrows"), ref_ncols=ref.slot("ncols"), trg_nrows = trg.slot("nrows"), trg_ncols = trg.slot("ncols");

   // generate matrices
   arma::mat ref_array(Rcpp::NumericVector(refdata.slot("values")).begin(), ref_nrows*ref_ncols, ref_nlayers, false, true), trg_array(Rcpp::NumericVector(trgdata.slot("values")).begin(), trg_nrows*trg_ncols, trg_nlayers, false, true);

   // find valid (finite) values
   arma::uvec ref_valid = arma::find_finite(arma::sum(ref_array,1)), trg_valid = arma::find_finite(arma::sum(trg_array,1));
   arma::mat ref_array_valid = ref_array.rows(ref_valid), trg_array_valid = trg_array.rows(trg_valid);

   // compute MIC
   int n_layers = 1;
   Rprintf("Computing most influential covariate...\n");  R_FlushConsole();
   arma::vec micdet = Rcpp::as<arma::vec>(arma_micdet(ref_array_valid, trg_array_valid));

   Rcpp::StringVector layernames(n_layers);
   Rcpp::NumericVector extent_min(n_layers), extent_max(n_layers);
   arma::mat vals = arma::mat(trg_nrows*trg_ncols, n_layers, arma::fill::value(NA_REAL));

   vals.elem(trg_valid) = micdet;
   layernames(0)="mic";
   extent_min(0)=micdet.min();
   extent_max(0)=micdet.max();

   Rcout << "Updating raster values...\n";
   R_FlushConsole();

   // update raster values
   trgdata.slot("values") = vals;

   // update other slots
   trgfile.slot("name") = "";
   trgdata.slot("min") = extent_min;
   trgdata.slot("max") = extent_max;
   trgdata.slot("inmemory") = true;
   trgdata.slot("fromdisk") = false;
   trgdata.slot("haveminmax") = true;
   trgdata.slot("names") = layernames;
   trgdata.slot("nlayers") = n_layers;
   trgdata.slot("isfactor") = true;
   Rcpp::List L(n_layers);
   L[0] = Rcpp::DataFrame::create(Named("ID") = Rcpp::Range(1,trg_nlayers), Named("Covariate") = trg_names);
   trgdata.slot("attributes") = L;

   return raster_out;
 }

 //' @title ExDet - Extrapolation Detection function
 //'
 //' @author Ian Ondo
 //'
 //' @description Given a target and reference dataset, computes the extrapolation index
 //' as defined in Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014.
 //' @param trg A N x M numerix matrix of N observations from M variables.
 //' @param ref A N'x M numeric matrix of N' observations from the same M variables as `trg` taken as a reference dataset.
 //' @return A numeric vector of novelty values of the same length as the number of rows
 //' of the target data matrix.
 //' @details Negative values denote novelty conditions of type 1, values between 0 and 1
 //' denote no novelty, i.e. conditions in the target dataset similar to the reference dataset,
 //' and values > 1 denote novelty conditions of type 2.
 //' @name exdet
 //' @references Mesgaran, M.B., Cousens, R.D. & Webber, B.L. (2014)
 //'  \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models}. Diversity & Distributions, 20: 1147-1159, DOI: 10.1111/ddi.12209
 //' @export
 // [[Rcpp::export]]
 arma::vec arma_exdet(arma::mat &trg, arma::mat &ref){

   // test for novelty of type 1
   arma::vec NT1 = arma_calc_nt1(trg, ref);

   // test for novelty of type 2 on remaining points (i.e. inside the range where NT1=0)
   arma::uvec is_in_range = arma::find(NT1==0);

   if(is_in_range.is_empty()){
     Rcpp::warning("arma_exdet(): cannot compute novelty of type 2, because no observation falls within the range of the reference data.");
     return NT1;
   }
   arma::mat trg_in_range = trg.rows(is_in_range);
   arma::vec NT2 = arma_calc_nt2(trg_in_range,ref);
   NT1.elem(is_in_range) = NT2;

   return NT1;
 }

 //' @title ExDet - Extrapolation Detection function for rasters
 //'
 //' @description Given a projection and reference raster, computes the extrapolation index
 //' as defined in Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014.
 //' @param trg A RasterBrick object of M layers of environmental conditions.
 //' @param ref A RasterBrick object with the same layers as `trg` but representing environmental conditions of reference.
 //' @param compute_mic. A logical. Should the most influential covariate be also calculated ? Default is `FALSE`.
 //' @return A RasterLayer of novelty values of the same dimensions as the `trg`.
 //' @details Negative values denote novelty conditions of type 1, values between 0 and 1
 //' denote no novelty, i.e. conditions in the target raster similar to conditions in the reference raster,
 //' and values > 1 denote novelty conditions of type 2.
 //' @references Mesgaran, M.B., Cousens, R.D. & Webber, B.L. (2014)
 //'  \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models}. Diversity & Distributions, 20: 1147-1159, DOI: 10.1111/ddi.12209
 //' @export
 // [[Rcpp::export]]
 Rcpp::S4 arma_exdet_raster(Rcpp::S4 &trg, Rcpp::S4 &ref, bool compute_mic=false){

   Rcpp::S4 raster_out = Rcpp::clone(trg);
   Rcpp::S4 refdata(ref.slot("data"));
   Rcpp::S4 trgdata(raster_out.slot("data")), trgfile(raster_out.slot("file"));

   // compare number of raster layers
   arma::uword ref_nlayers = refdata.slot("nlayers"), trg_nlayers = refdata.slot("nlayers");
   if(ref_nlayers!=trg_nlayers) Rcpp::stop("arma_exdet_raster(): the reference and target rasters must have the same number of layers.");

   // compare raster layers names
   CharacterVector ref_names = refdata.slot("names"), trg_names = trgdata.slot("names");
   CharacterVector::iterator refl = ref_names.begin(), trgl=trg_names.begin();

   for(; refl != ref_names.end(); ++refl, ++trgl) if(*refl!=*trgl) Rcpp::stop("arma_exdet_raster(): layers names are not identical.");

   // get target raster info
   arma::uword ref_nrows=ref.slot("nrows"), ref_ncols=ref.slot("ncols"), trg_nrows = trg.slot("nrows"), trg_ncols = trg.slot("ncols");

   // generate matrices
   arma::mat ref_array(Rcpp::NumericVector(refdata.slot("values")).begin(), ref_nrows*ref_ncols, ref_nlayers, false, true), trg_array(Rcpp::NumericVector(trgdata.slot("values")).begin(), trg_nrows*trg_ncols, trg_nlayers, false, true);

   // find valid (finite) values
   arma::uvec ref_valid = arma::find_finite(arma::sum(ref_array,1)), trg_valid = arma::find_finite(arma::sum(trg_array,1));
   arma::mat ref_array_valid = ref_array.rows(ref_valid), trg_array_valid = trg_array.rows(trg_valid);

   // compute extrapolation
   int n_layers = 1;
   Rcout << "Computing uni/multivariate extrapolation...\n";  R_FlushConsole();
   arma::vec micdet, exdet = arma_exdet(trg_array_valid, ref_array_valid);

   // compute mic
   if(compute_mic){
     Rprintf("Computing most influential covariate...\n");  R_FlushConsole();
     n_layers+=1;
     micdet = Rcpp::as<arma::vec>(arma_micdet(trg_array_valid, ref_array_valid));
   }

   Rcpp::StringVector layernames(n_layers);
   Rcpp::NumericVector extent_min(n_layers), extent_max(n_layers);
   Rcpp::LogicalVector is_factor(n_layers);
   arma::mat vals = arma::mat(trg_nrows*trg_ncols, n_layers, arma::fill::value(NA_REAL));

   vals.elem(trg_valid) = exdet;
   layernames(0)="exdet";
   extent_min(0)=exdet.min();
   extent_max(0)=exdet.max();
   is_factor(0)=false;

   if(n_layers>1){
     trg_valid+=trg_nrows*trg_ncols;
     vals.elem(trg_valid) = micdet;
     layernames(1)="mic";
     extent_min(1)=micdet.min();
     extent_max(1)=micdet.max();
     is_factor(1)=true;
   }

   Rcout << "Updating raster values...\n";
   R_FlushConsole();

   // update raster values
   trgdata.slot("values") = vals;

   // update other slots
   trgfile.slot("name") = "";
   trgdata.slot("min") = extent_min;
   trgdata.slot("max") = extent_max;
   trgdata.slot("inmemory") = true;
   trgdata.slot("fromdisk") = false;
   trgdata.slot("haveminmax") = true;
   trgdata.slot("names") = layernames;
   trgdata.slot("nlayers") = n_layers;

   if(compute_mic){
     trgdata.slot("isfactor") = is_factor;
     Rcpp::List L(n_layers);
     L[0] = Rcpp::DataFrame::create(Named("ID") = 1, Named("Covariate") = NA_INTEGER);
     L[1] = Rcpp::DataFrame::create(Named("ID") = Rcpp::Range(1,trg_nlayers), Named("Covariate") = trg_names);
     trgdata.slot("attributes") = L;
   }

   return raster_out;
 }

