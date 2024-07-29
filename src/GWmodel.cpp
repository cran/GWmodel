// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#define STRICT_R_HEADERS
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppEigen.h>
#include <vector>
#include <limits>
#include <functional>
#include <string>

using namespace Rcpp;
using namespace arma;
using namespace stats;

#define POWDI(x,i) pow(x,i)

#define GAUSSIAN 0
#define EXPONENTIAL 1
#define BISQUARE 2
#define TRICUBE 3
#define BOXCAR 4

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
// For debugging
void printVec(vec v) {
  int n = v.size();
  n = 10;
  
  for (int i=0; i<n; i++) {
    Rprintf("%f ", v(i));
  }
  Rprintf("\n");
}
// For debugging
void printMat(mat m) {
  uword n = m.n_rows;
  uword p = m.n_cols;
  
  n = 10;
  if (m.n_rows < n) 
  {
     n = m.n_rows;
  } 
  for (uword i=0; i<n; i++) {
    for (uword j=0; j<p; j++) {
      Rprintf("%f ", m(i, j));
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}
//-------------------------------------------Distance metric calculations
const mat mSum(2, 1, fill::ones);

//distance matrix calculation
//coords must be a matrix with 2 columns
mat coordinate_rotate(mat coords, double theta)
{
	int n = coords.n_rows;
	mat rotated_coords(n, 2);
	rotated_coords.col(0) = coords.col(0) * cos(theta) - coords.col(1) * sin(theta);
	rotated_coords.col(1) = coords.col(0) * sin(theta) + coords.col(1) * cos(theta);
	return rotated_coords;
}

//Eudclidean distance matrix
mat eu_dist_mat(mat in_locs, mat out_locs)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat eu_dist(n_in, n_out);
	int i = 0, j = 0;
	for (i = 0; i < n_in; i++)
	{
		for (j = 0; j < n_out; j++)
		{
		  eu_dist(i,j) = sum(pow(in_locs.row(i) - out_locs.row(j),2));
		}
	}
	return sqrt(eu_dist);
}
//symmetrical distance matrix
mat eu_dist_smat(mat in_locs)
{
	int n = in_locs.n_rows;
	mat eu_dist(n, n);
	for (int k = 0; k < n * n; k++)
	{
		int i = k / n, j = k % n;
		eu_dist(i, j) = sum(pow(in_locs.row(i) - in_locs.row(j), 2));
		eu_dist(j, i) = eu_dist(i, j);
	}
	return sqrt(eu_dist);
}

vec eu_dist_vec(mat in_locs, vec out_loc)
{
	int n_in = in_locs.n_rows;
	vec eu_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		eu_dist(i) = sum(pow(in_locs.row(i) - trans(out_loc), 2));
	}
	return sqrt(eu_dist);
	// mat v_span(n_in, 1, fill::ones);
	// mat m_diff = in_locs - v_span * trans(out_loc);
	// return sqrt(m_diff % m_diff * mSum);
}

//Manhattan distance matrix
mat md_dist_mat(mat in_locs, mat out_locs)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat md_dist(n_in, n_out);
	for (int i = 0; i < n_in; i++)
	{
		for (int j = 0; j < n_out; j++)
		{
			md_dist(i, j) = sum(abs(in_locs.row(i) - out_locs.row(j)));
		}
	}
	return md_dist;
}

//symmetrical distance matrix
mat md_dist_smat(mat in_locs)
{
	int n = in_locs.n_rows;
	mat md_dist(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			md_dist(i, j) = sum(abs(in_locs.row(i) - in_locs.row(j)));
			md_dist(j, i) = md_dist(i, j);
		}
	}
	return md_dist;
}
vec md_dist_vec(mat in_locs, vec out_loc)
{
	int n_in = in_locs.n_rows;
	vec md_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		md_dist(i) = sum(abs(in_locs.row(i) - trans(out_loc)));
	}
	return md_dist;
}

//Chebyshev distance matrix
mat cd_dist_mat(mat in_locs, mat out_locs)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat cd_dist(n_in, n_out);
	for (int i = 0; i < n_in; i++)
	{
		for (int j = i; j < n_out; j++)
		{
			cd_dist(i, j) = max(abs(in_locs.row(i) - out_locs.row(j)));
			cd_dist(j, i) = cd_dist(i, j);
		}
	}
	return cd_dist;
}

//symmetrical distance matrix
mat cd_dist_smat(mat in_locs)
{
	int n = in_locs.n_rows;
	mat cd_dist(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			cd_dist(i, j) = max(abs(in_locs.row(i) - in_locs.row(j)));
			cd_dist(j, i) = cd_dist(i, j);
		}
	}
	return cd_dist;
}
vec cd_dist_vec(mat in_locs, vec out_loc)
{
	int n_in = in_locs.n_rows;
	vec cd_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		cd_dist(i) = max(abs(in_locs.row(i) - trans(out_loc)));
	}
	return cd_dist;
}

//Minkowski distance matrix
mat mk_dist_mat(mat in_locs, mat out_locs, double p)
{
	int n_in = in_locs.n_rows;
	int n_out = out_locs.n_rows;
	mat mk_dist(n_in, n_out);
	for (int i = 0; i < n_in; i++)
	{
		for (int j = 0; j < n_out; j++)
		{
			mk_dist(i, j) = pow(sum(pow(abs(in_locs.row(i) - out_locs.row(j)), p)), 1.0 / p);
		}
	}
	return mk_dist;
}
//symmetrical distance matrix
mat mk_dist_smat(mat in_locs, double p)
{
	int n = in_locs.n_rows;
	mat mk_dist(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			mk_dist(i, j) = pow(sum(pow(abs(in_locs.row(i) - in_locs.row(j)), p)), 1.0 / p);
			mk_dist(j, i) = mk_dist(i, j);
		}
	}
	return mk_dist;
}
vec mk_dist_vec(mat in_locs, vec out_loc, double p)
{
	int n_in = in_locs.n_rows;
	vec mk_dist(n_in);
	for (int i = 0; i < n_in; i++)
	{
		mk_dist(i) = pow(sum(pow(abs(in_locs.row(i) - trans(out_loc)), p)), 1.0 / p);
	}
	return mk_dist;
}
//Great circle distance
// Caculate the distance with a pair of lat and long
double sp_gcdist(double lon1, double lon2, double lat1, double lat2) {
  
  double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
  double w, R, a, f, D, H1, H2;
  double lat1R, lat2R, lon1R, lon2R, DE2RA;
  
  DE2RA = M_PI/180;
  a = 6378.137;              /* WGS-84 equatorial radius in km */
    f = 1.0/298.257223563;     /* WGS-84 ellipsoid flattening factor */
    
    if (fabs(lat1 - lat2) < DBL_EPSILON) {
      if (fabs(lon1 - lon2) < DBL_EPSILON) {
        return 0.0;
        /* Wouter Buytaert bug caught 100211 */
      } else if (fabs((fabs(lon1) + fabs(lon2)) - 360.0) < DBL_EPSILON) {
        return 0.0;
      }
    }
    lat1R = lat1*DE2RA;
    lat2R = lat2*DE2RA;
    lon1R = lon1*DE2RA;
    lon2R = lon2*DE2RA;
    
    F = ( lat1R + lat2R )/2.0;
    G = ( lat1R - lat2R )/2.0;
    L = ( lon1R - lon2R )/2.0;
    
    /*
    printf("%g %g %g %g; %g %g %g\n",  *lon1, *lon2, *lat1, *lat2, F, G, L);
    */
    
    sinG2 = POWDI( sin( G ), 2 );
    cosG2 = POWDI( cos( G ), 2 );
    sinF2 = POWDI( sin( F ), 2 );
    cosF2 = POWDI( cos( F ), 2 );
    sinL2 = POWDI( sin( L ), 2 );
    cosL2 = POWDI( cos( L ), 2 );
    
    S = sinG2*cosL2 + cosF2*sinL2;
    C = cosG2*cosL2 + sinF2*sinL2;
    
    w = atan( sqrt( S/C ) );
    R = sqrt( S*C )/w;
    
    D = 2*w*a;
    H1 = ( 3*R - 1 )/( 2*C );
    H2 = ( 3*R + 1 )/( 2*S );
    
    return D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 ); 
}

// Calculate the distance vector between a point and a set of points, latitude and longitude required
vec sp_dists(mat dp, vec loc) {
  int N = dp.n_rows, j;
  vec dists(N, fill::zeros);
  double uout = loc(0), vout = loc(1);
  for (j = 0; j < N; j++) {
    dists(j) = sp_gcdist(dp(j, 0), uout, dp(j, 1), vout);
  }
  return dists;
}

// Equal to gw.dist, to be checked
// [[Rcpp::export]]
mat gw_dist(mat dp, mat rp, int focus, double p, double theta, bool longlat, bool rp_given) {
  int ndp = dp.n_rows, nrp = rp.n_rows;
  int isFocus = focus > -1;
  mat dists;
  if (p != 2 && theta != 0 && !longlat) {
    dp = coordinate_rotate(dp, theta);
    rp = coordinate_rotate(rp, theta);
  }
  if (isFocus) {
    mat prp = trans(rp.row(focus));
    if (longlat) {
      return sp_dists(dp, prp);
    } else {
      if (p == 2.0)
        return eu_dist_vec(dp, prp);
      else if(p == -1.0)
        return cd_dist_vec(dp, prp);
      else if(p == 1.0)
        return md_dist_vec(dp, prp);
      else
        return mk_dist_vec(dp, prp, p);
    }
  } else {
    if (longlat) {
      mat dists(ndp, nrp, fill::zeros);
      for (int i = 0; i < nrp; i++) {
        dists.col(i) = sp_dists(dp, trans(rp.row(i)));
      }
      return trans(dists);
    } else {
      if (p == 2.0)
        return rp_given ? eu_dist_mat(dp, rp) : eu_dist_smat(dp);
      else if (p == -1.0)
        return rp_given ? cd_dist_mat(dp, rp) : cd_dist_smat(dp);
      else if (p == 1.0)
        return rp_given ? md_dist_mat(dp, rp) : md_dist_smat(dp);
      else
        return rp_given ? mk_dist_mat(dp, rp, p) : mk_dist_smat(dp, p);
    }
  }
}
//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------Weight matrix caculation-------------------------------------------------------------------------------------
//Calucate the GW weights
double gw_weight_gaussian(double dist, double bw) {
  return exp(pow(dist, 2)/((-2)*pow(bw, 2)));
}

double gw_weight_exponential(double dist, double bw) {
  return exp(-dist/bw);
}

double gw_weight_bisquare(double dist, double bw) {
  return dist > bw ? 0 : pow(1 - pow(dist, 2)/pow(bw, 2), 2);
}

double gw_weight_tricube(double dist, double bw) {
  return dist > bw ? 0 : pow(1 - pow(dist, 3)/pow(bw, 3), 3);
}

double gw_weight_boxcar(double dist, double bw) {
  return dist > bw ? 0 : 1;
}

typedef double (*KERNEL)(double, double);
const KERNEL GWRKernel[5] = {
  gw_weight_gaussian,
  gw_weight_exponential,
  gw_weight_bisquare,
  gw_weight_tricube,
  gw_weight_boxcar
};

//gwr.weight in the R code
// [[Rcpp::export]]
mat gw_weight(mat dist, double bw, int kernel, bool adaptive) {
  const KERNEL *kerf = GWRKernel + kernel;
  int nr = dist.n_rows, nc = dist.n_cols;
  mat w(nr, nc, fill::zeros);
  if (adaptive) {
    for (int c = 0; c < nc; c++) {
      double dn = bw / nr, fixbw = 0;
      if (dn <= 1) {
        vec vdist = sort(dist.col(c));
        fixbw = vdist(int(bw) - 1);
      } else {
        fixbw = dn * max(dist.col(c));
      }
      for (int r = 0; r < nr; r++) {
        w(r, c) = (*kerf)(dist(r, c), fixbw);
      }
    }
  } else {
    for (int c = 0; c < nc; c++) {
      for (int r = 0; r < nr; r++) {
        w(r, c) = (*kerf)(dist(r, c), bw);
      }
    }
  }
  return w;
}
//------------------------------------------------------
// Geographic weights
// [[Rcpp::export]]
vec gw_weight_vec(vec vdist, double bw, int kernel, bool adaptive)
{
  const KERNEL *kerf = GWRKernel + kernel;
  int n = vdist.n_elem;
  vec wv(n, fill::zeros);
  if (adaptive) {
    for (int i = 0; i < n; i++) {
      double dn = bw / n, fixbw = 0;
      if (dn <= 1) {
        vec svdist = sort(vdist);
        fixbw = svdist(int(bw) - 1);
      } else {
        fixbw = dn * max(vdist);
      }
      wv(i) = (*kerf)(vdist(i), fixbw);
    }
  } else {
     for (int i = 0; i < n; i++) {
        wv(i) = (*kerf)(vdist(i), bw);
    }
  }
  return wv;
}
// [[Rcpp::export]]
mat gw_weight_mat(mat dist, double bw, int kernel, bool adaptive)
{
  const KERNEL *kerf = GWRKernel + kernel;
  int nr = dist.n_rows, nc = dist.n_cols;
  mat w(nr, nc, fill::zeros);
  if (adaptive) {
    for (int c = 0; c < nc; c++) {
      double dn = bw / nr, fixbw = 0;
      if (dn <= 1) {
        vec vdist = sort(dist.col(c));
        fixbw = vdist(int(bw) - 1);
      } else {
        fixbw = dn * max(dist.col(c));
      }
      for (int r = 0; r < nr; r++) {
        w(r, c) = (*kerf)(dist(r, c), fixbw);
      }
    }
  } else {
    for (int c = 0; c < nc; c++) {
      for (int r = 0; r < nr; r++) {
        w(r, c) = (*kerf)(dist(r, c), bw);
      }
    }
  }
  return w;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------GWR clalibration calucations---------------------------------------------------------------------------------
//GWR clalibration
// [[Rcpp::export]]
List gw_reg(mat x, vec y, vec w, bool hatmatrix, int focus)
{
	mat wspan(1, x.n_cols, fill::ones);
	mat xtw = trans(x % (w * wspan));
	mat xtwx = xtw * x;
	mat xtwy = trans(x) * (w % y);
	mat xtwx_inv = inv(xtwx);
	vec beta = xtwx_inv * xtwy;
	if (hatmatrix)
	{
		mat ci = xtwx_inv * xtw;
		mat s_ri = x.row(focus - 1) * ci;
		return List::create(
				Named("beta") = beta,
				Named("S_ri") = s_ri,
				Named("Ci") = ci);
	}
	else
	{
		return List::create(
				Named("beta") = beta);
	}
}

//GWR clalibration
// [[Rcpp::export]]
List gw_reg_1(mat x, vec y, vec w)
{
	mat wspan(1, x.n_cols, fill::ones);
	mat xtw = trans(x % (w * wspan));
	mat xtwx = xtw * x;
	mat xtwy = trans(x) * (w % y);
	mat xtwx_inv = inv(xtwx);
	vec beta = xtwx_inv * xtwy;
	return List::create(
			Named("beta") = beta,
			Named("xtwx_inv") = xtwx_inv);
}
// Trace of hat matrix + trace of HH' in one function. Used in beta_se

// [[Rcpp::export]]
vec trhat2(mat S)
{
	int n_obs = S.n_rows;
	double htr = 0.0;
	double htr2 = 0.0;
	vec result(2);
	for (int i = 0; i < n_obs; i++)
	{
		htr += S(i, i);
		htr2 += sum(S.row(i) % S.row(i));
	}
	result(0) = htr;
	result(1) = htr2;
	return result;
}
// Fited values
// [[Rcpp::export]]
vec fitted(mat X, mat beta)
{
	vec fitted = sum(beta % X, 1);
	return fitted;
}

// Residuals

// [[Rcpp::export]]
vec ehat(vec y, mat X, mat beta)
{
	vec fitted = sum(beta % X, 1);
	return y - fitted;
}

// Residual sum-of squares
// [[Rcpp::export]]
double rss(vec y, mat X, mat beta)
{
	vec r = ehat(y, X, beta);
	return sum(r % r);
}

//Calculate the diagnositic information
// [[Rcpp::export]]
vec gwr_diag(vec y,mat x, mat beta, mat S) {
  double ss = rss(y,x,beta);
  vec s_hat = trhat2(S);
  double n = (double) S.n_rows;
  vec result(10);
  double AIC = n*log(ss/n)+n*log(2*datum::pi)+n+s_hat(0); //AIC
  double AICc = n*log(ss/n)+n*log(2*datum::pi)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc
  double edf = n-2*s_hat(0) + s_hat(1); //edf
  double enp = 2*s_hat(0) - s_hat(1); // enp
  double yss = sum(pow(y-mean(y),2)); //yss.g
  double r2 = 1 - ss/yss;
  double r2_adj = 1-(1-r2)*(n-1)/(edf-1);
  double BIC = n * log(ss / n) + n * log(2.0 * datum::pi) + log(n) * s_hat(0);
  result(0) = AIC;
  result(1) = AICc;
  result(2) = edf;
  result(3) = enp;
  result(4) = ss;
  result(5) = r2;
  result(6) = r2_adj;
  result(7) = s_hat(0);
  result(8) = s_hat(1);
  result(9) = BIC;
  return result;
}

//Calculate the diagnositic information, no hatmatrix trace is returned
// [[Rcpp::export]]
vec gwr_diag1(vec y, mat x, mat beta, vec s_hat)
{
	double ss = rss(y, x, beta);
	// vec s_hat = trhat2(S);
	double n = (double) x.n_rows;
	// vec result(9);
	vec result(8);
	double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);																//AIC
	double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
	double edf = n - 2 * s_hat(0) + s_hat(1);																														//edf
	double enp = 2 * s_hat(0) - s_hat(1);																																// enp
	double yss = sum(pow(y - mean(y), 2));																															//yss.g
	double r2 = 1 - ss / yss;
	double r2_adj = 1 - (1 - r2) * (n - 1) / (edf - 1);
  double BIC = n * log(ss / n) + n * log(2.0 * datum::pi) + log(n) * s_hat(0);
	result(0) = AIC;
	result(1) = AICc;
	result(2) = edf;
	result(3) = enp;
	result(4) = ss;
	result(5) = r2;
	result(6) = r2_adj;
	result(7) = BIC;
	// result(8) = s_hat(1);
	return result;
	//return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc - nice tool to give to 'optimise' to find 'best' bandwidth
// [[Rcpp::export]]
double AICc(vec y, mat x, mat beta, mat S)
{
	double ss = rss(y, x, beta);
	vec s_hat = trhat2(S);
	int n = S.n_rows;
	double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
	return AICc;
	//return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}
//Could be dropped
// [[Rcpp::export]]
double AICc1(vec y, mat x, mat beta, vec s_hat)
{
  double ss = rss(y, x, beta);
  // vec s_hat = trhat2(S);
  int n = x.n_rows;
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  return AICc;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc and RSS , used for function model.selection
// [[Rcpp::export]]
vec AICc_rss(vec y, mat x, mat beta, mat S)
{
	vec result(3);
	double ss = rss(y, x, beta);
	result[0] = ss;
	vec s_hat = trhat2(S);
	int n = S.n_rows;
	double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);
	double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
	result[1] = AIC;
	result[2] = AICc;
	return result;
	//return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// [[Rcpp::export]]
vec AICc_rss1(vec y, mat x, mat beta, vec s_hat)
{
  vec result(3);
  double ss = rss(y, x, beta);
  result[0] = ss;
  int n = x.n_rows;
  double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  result[1] = AIC;
  result[2] = AICc;
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

//Caculate the i row of hatmatirx
// [[Rcpp::export]]
mat Ci_mat(mat x, vec w)
{
	return inv(trans(x) * diagmat(w) * x) * trans(x) * diagmat(w);
}

//Local R2 values for GWR
// [[Rcpp::export]]
vec gw_local_r2(mat dp, vec dybar2, vec dyhat2, bool dm_given, mat dmat, double p, double theta, bool longlat, double bw, int kernel, bool adaptive) {
  int n = dp.n_rows;
  vec localR2(n, fill::zeros);
  for (int i = 0; i < n; i++) {
    mat d = dm_given ? dmat.col(i) : gw_dist(dp, dp, i, p, theta, longlat, false);
    mat w = gw_weight(d, bw, kernel, adaptive);
    double tss = sum(dybar2 % w);
    double rss = sum(dyhat2 % w);
    localR2(i) = (tss - rss) / tss;
  }
  return localR2;
}
//BIC calculation
// [[Rcpp::export]]
double BIC(vec y, mat x, mat beta, vec s_hat)
{
  double ss = rss(y, x, beta);
  double n = (double)x.n_rows;
  double BIC = n * log(ss / n) + n * log(2 * datum::pi) + log(n) * s_hat(0);
  return BIC;
}

// GWR calibration, returns betas only
// [[Rcpp::export]]
mat gw_reg_2(mat x, vec y, vec w)
{
  mat beta;
  mat wspan(1, x.n_cols, fill::ones);
  mat xtw = trans(x % (w * wspan));
  mat xtwx = xtw * x;
  mat xtwy = trans(x) * (w % y);
  mat xtwx_inv = inv(xtwx);
  beta = xtwx_inv * xtwy;
  return beta.t();
}
// C++ version of gwr.q, used in gwr.mixed
// [[Rcpp::export]]
mat gwr_q(mat x, vec y, 
                mat dMat, double bw, int kernel, bool adaptive)
{
  int n =  dMat.n_cols;
  int m =  x.n_cols;
  mat beta(n, m);
  vec distv;
  vec w;
  for (int i = 0; i < n; i++) {
    distv = dMat.col(i);
    w = gw_weight_vec(distv, bw, kernel, adaptive);
    beta.row(i) = gw_reg_2(x, y, w);
  }
  return beta;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------Scalable GWR clalibration ---------------------------------------------------------------------------------
// Scalable GWR bandwidth optimization via CV approach
// [[Rcpp::export]]
double scgwr_loocv(vec target, mat x, vec y, int bw, int poly, mat Mx0, mat My0, mat XtX, mat XtY) {
  int n = x.n_rows, k = x.n_cols, poly1 = poly + 1;
  double b = target(0) * target(0), a = target(1) * target(1);
  vec R0 = vec(poly1) * b;
  for (int p = 1; p < poly1; p++) {
    R0(p) = pow(b, p + 1);
  }
  vec Rx(k*k*poly1, fill::zeros), Ry(k*poly1, fill::zeros);
  for (int p = 0; p < poly1; p++) {
    for (int k2 = 0; k2 < k; k2++) {
      for (int k1 = 0; k1 < k; k1++) {
        int xindex = k1*poly1*k + p*k + k2;
        Rx(xindex) = R0(p);
      }
      int yindex = p*k + k2;
      Ry(yindex) = R0(p);
    }
  }
  mat Mx = Rx * mat(1, n, fill::ones) % Mx0, My = Ry * mat(1, n, fill::ones) % My0;
  vec yhat(n, 1, fill::zeros);
  for (int i = 0; i < n; i++) {
    mat sumMx(k, k, fill::zeros);
    vec sumMy(k, fill::zeros);
    for (int k2 = 0; k2 < k; k2++) {
      for (int p = 0; p < poly1; p++) {
        for (int k1 = 0; k1 < k; k1++) {
          int xindex = k1*poly1*k + p*k + k2;
          sumMx(k1, k2) += Mx(xindex, i);
        }
        int yindex = p*k + k2;
        sumMy(k2) += My(yindex, i);
      }
    }
    sumMx = sumMx + a * XtX;
    sumMy = sumMy + a * XtY;
    if (det(sumMx) < 1e-10) {
      return 1e6;
    } else {
      mat beta = solve(sumMx, sumMy);
      yhat.row(i) = x.row(i) * beta;
    }
  }
  return sum((y - yhat) % (y - yhat));
}
//Scalable GWR C++ functions
// [[Rcpp::export]]
List scgwr_pre(mat x, vec y, int bw, int poly, double b0, mat g0, mat neighbour) {
  int n = x.n_rows;
  int k = x.n_cols;
  mat g0s(g0.n_cols + 1, g0.n_rows, fill::ones);
  mat g0t = trans(g0);
  for (int i = 0; i < bw; i++) {
    g0s.row(i + 1) = g0t.row(i);
  }
  g0s = trans(g0s);
  mat Mx0((poly + 1)*k*k, n, fill::zeros);
  mat My0((poly + 1)*k, n, fill::zeros);
  mat spanXnei(1, poly + 1, fill::ones);
  mat spanXtG(1, k, fill::ones);
  for (int i = 0; i < n; i++) {
    mat g(poly + 1, bw + 1, fill::ones);
    for (int p = 0; p < poly; p++) {
      g.row(p + 1) = pow(g0s.row(i), pow(2.0, poly/2.0)/pow(2.0, p + 1));
    }
    g = trans(g);
    g = g.rows(1, bw);
    mat xnei(bw, k, fill::zeros);
    vec ynei(bw, fill::zeros);
    for (int j = 0; j < bw; j++) {
      int inei = int(neighbour(i, j) - 1);
      xnei.row(j) = x.row(inei);
      ynei.row(j) = y(inei);
    }
    for (int k1 = 0; k1 < k; k1++) {
      mat XtG = xnei.col(k1) * spanXnei % g;
      for (int p = 0; p < (poly + 1); p++) {
        mat XtGX = XtG.col(p) * spanXtG % xnei;
        for (int k2 = 0; k2 < k; k2++) {
          int xindex = (k1 * (poly + 1) + p) * k + k2;
          Mx0(xindex, i) = sum(XtGX.col(k2));
        }
        int yindex = p * k + k1;
        vec XtGY = XtG.col(p) % ynei;
        My0(yindex, i) = sum(XtGY);
      }
    }
  }
  return List::create(
    Named("Mx0") = Mx0,
    Named("My0") = My0
  );
}
// Scalable GWR calibration
// [[Rcpp::export]]
List scgwr_reg(mat x, vec y, int bw, int poly, mat G0, mat Mx0, mat My0, mat XtX, mat XtY, mat neighbour, vec parameters) {
  int n = x.n_rows, k = x.n_cols, poly1 = poly + 1;
  double b = parameters(0), a = parameters(1);
  /**
   * Calculate Rx, Ry, and R0.
   */
  // printf("Calculate Rx, Ry, and R0 ");
  vec R0 = vec(poly1, fill::ones) * b;
  for (int p = 1; p < poly1; p++) {
    R0(p) = pow(b, p + 1);
  }
  vec Rx(k*k*poly1, fill::zeros), Ry(k*poly1, fill::zeros);
  for (int p = 0; p < poly1; p++) {
    for (int k2 = 0; k2 < k; k2++) {
      for (int k1 = 0; k1 < k; k1++) {
        Rx(k1*poly1*k + p*k + k2) = R0(p);
      }
      Ry(p*k + k2) = R0(p);
    }
  }
  /**
  * Create G0.
  */
  // printf("Create G0 ");
  mat G0s(G0.n_cols + 1, G0.n_rows, fill::ones);
  mat G0t = trans(G0);
  G0s.rows(1, bw) = G0t.rows(0, bw - 1);
  G0s = trans(G0s);
  mat G2(n, bw + 1, fill::zeros);
  for (int p = 0; p < poly; p++) {
    G2 += pow(G0s, pow(2.0, poly/2.0)/pow(2.0, p + 1)) * R0(poly - 1);
  }
  /**
   * Calculate Mx, My.
   */
  // printf("Calculate Mx, My ");
  // mat Mx00(Mx0), My00(My0);
  for (int i = 0; i < n; i++) {
    for (int k1 = 0; k1 < k; k1++) {
      for (int p = 0; p < poly1; p++) {
        for (int k2 = 0; k2 < k; k2++) {
          Mx0((k1 * (poly + 1) + p) * k + k2, i) += x(i, k1) * x(i, k2);
        }
        My0(p * k + k1, i) += x(i, k1) * y(i);
      }
    }
  }
  mat Mx = (Rx * mat(1, n, fill::ones)) % Mx0, My = (Ry * mat(1, n, fill::ones)) % My0;
  /**
   * Regression.
   */
  // printf("Regression ");
  mat Xp(bw + 1, k * poly1, fill::zeros);
  mat rowsumG(poly1, 1, fill::ones);
  mat colsumXp(1, bw + 1, fill::zeros);
  mat spanG(1, k, fill::ones);
  mat spanX(1, poly1, fill::ones);
  mat betas(n, k, fill::zeros);
  mat betasSE(n, k, fill::zeros);
  double trS = 0.0, trStS = 0.0;
  for (int i = 0; i < n; i++) {
    /**
     * Calculate G.
     */
    mat G = mat(poly1, bw + 1, fill::ones) * R0(0);
    for (int p = 0; p < poly; p++) {
      G.row(p + 1) = pow(G0s.row(i), pow(2.0, poly/2.0)/pow(2.0, p + 1)) * R0(p);
    }
    G = trans(G);
    mat g = G * rowsumG;
    /**
     * Calculate Xp.
     */
    mat xnei(bw + 1, k, fill::zeros);
    vec ynei(bw + 1, fill::zeros);
    xnei.row(0) = x.row(i);
    ynei.row(0) = y.row(i);
    for (int j = 0; j < bw; j++) {
      int inei = int(neighbour(i, j) - 1);
      xnei.row(j+1) = x.row(inei);
      ynei.row(j+1) = y(inei);
    }
    /**
     * Calculate sumMx, sumMy.
     */
    mat sumMx(k, k, fill::zeros);
    vec sumMy(k, fill::zeros);
    for (int k2 = 0; k2 < k; k2++) {
      for (int p = 0; p < poly1; p++) {
        for (int k1 = 0; k1 < k; k1++) {
          int xindex = k1*poly1*k + p*k + k2;
          sumMx(k1, k2) += Mx(xindex, i);
        }
        int yindex = p*k + k2;
        sumMy(k2) += My(yindex, i);
      }
    }
    sumMx += a * XtX;
    sumMy += a * XtY;
    mat invMx = inv(sumMx);
    betas.row(i) = trans(invMx * sumMy);
    /**
     * Calculate Diagnostic statistics, trS and trStS.
     */
    mat StS = invMx * trans(x.row(i));
    trS += det((x.row(i) * g(0, 0)) * StS);
    mat XG2X(k, k, fill::zeros);
    for (int k1 = 0; k1 < k; k1++) {
      for (int k2 = 0; k2 < k; k2++) {
        mat Gi = G2.row(i);
        XG2X(k1, k2) = sum(xnei.col(k1) % trans(Gi % Gi + 2 * a * Gi) % xnei.col(k2)) + a * a * XtX(k1, k2);
      }
    }
    mat XX = invMx * XG2X * invMx;
    mat xi = x.row(i);
    trStS += det(sum(xi * XX * trans(xi)));
    betasSE.row(i) = trans(sqrt(XX.diag()));
  }
  return List::create(
    Named("betas") = betas,
    Named("tr.S") = trS,
    Named("tr.StS") = trStS,
    Named("betas.SE") = betasSE
  );
}

//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------  Parallel GWR clalibration ---------------------------------------------------------------------------------
// GWR calibration for parallel computation
// [[Rcpp::export]]
List gw_reg_all(mat x, vec y, mat dp, bool rp_given, mat rp, bool dm_given, mat dmat, bool hatmatrix, 
                double p, double theta, bool longlat, 
                double bw, int kernel, bool adaptive,
                int ngroup, int igroup) {
  int n = rp.n_rows, k = x.n_cols;
  mat betas(n, k, fill::zeros);
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  if (hatmatrix) {
    mat betasSE(n, k, fill::zeros);
    mat s_hat(1, 2, fill::zeros);
    mat qdiag(1, n, fill::zeros);
    mat rowsumSE(n, 1, fill::ones);
    // clock_t clock0 = clock(), clock1;
    for (int i = iStart; i < iEnd; i++) {
      mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
      mat w = gw_weight(d, bw, kernel, adaptive);
      mat ws(1, k, fill::ones);
      mat xtw = trans(x %(w * ws));
      mat xtwx = xtw * x;
      mat xtwy = trans(x) * (w % y);
      mat xtwx_inv = inv(xtwx);
      betas.row(i) = trans(xtwx_inv * xtwy);
      // hatmatrix
      mat ci = xtwx_inv * xtw;
      mat si = x.row(i) * ci;
      betasSE.row(i) = trans((ci % ci) * rowsumSE);
      s_hat(0) += si(0, i);
      s_hat(1) += det(si * trans(si));
      mat onei(1, n, fill::zeros);
      onei(i) = 1;
      mat p = onei - si;
      qdiag += p % p;
    }
    return List::create(
      Named("betas") = betas,
      Named("betas.SE") = betasSE,
      Named("s_hat") = s_hat,
      Named("q.diag") = qdiag
    );
  } else {
    for (int i = iStart; i < iEnd; i++) {
      mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
      mat w = gw_weight(d, bw, kernel, adaptive);
      mat ws(1, x.n_cols, fill::ones);
      mat xtw = trans(x %(w * ws));
      mat xtwx = xtw * x;
      mat xtwy = trans(x) * (w % y);
      mat xtwx_inv = inv(xtwx);
      betas.row(i) = trans(xtwx_inv * xtwy);
    }
    return List::create(
      Named("betas") = betas
    );
  }
}

// GWR calibration for multi-threads
// [[Rcpp::export]]
#ifdef _OPENMP
List gw_reg_all_omp(mat x, vec y, mat dp, bool rp_given, mat rp, bool dm_given, mat dmat, bool hatmatrix, 
                    double p, double theta, bool longlat, 
                    double bw, int kernel, bool adaptive,
                    int threads, int ngroup, int igroup) {
  int n = rp.n_rows, k = x.n_cols;
  mat betas(n, k, fill::zeros);
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  if (hatmatrix) {
    mat betasSE(n, k, fill::zeros);
    mat s_hat(1, 2, fill::zeros);
    mat qdiag(1, n, fill::zeros);
    mat rowsumSE(n, 1, fill::ones);
    vec s_hat1(n, fill::zeros), s_hat2(n, fill::zeros);
    int thread_nums = threads > 0 ? threads : omp_get_num_procs() - 1;
    mat qdiag_all(thread_nums, n, fill::zeros);
    bool flag_error = false;
#pragma omp parallel for num_threads(thread_nums)
    for (int i = iStart; i < iEnd; i++) {
      if (!flag_error) {
        int thread_id = omp_get_thread_num();
        mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
        mat w = gw_weight(d, bw, kernel, adaptive);
        mat ws(1, k, fill::ones);
        mat xtw = trans(x %(w * ws));
        mat xtwx = xtw * x;
        mat xtwy = trans(x) * (w % y);
        try {
          mat xtwx_inv = inv(xtwx);
          betas.row(i) = trans(xtwx_inv * xtwy);
          // hatmatrix
          mat ci = xtwx_inv * xtw;
          mat si = x.row(i) * ci;
          betasSE.row(i) = trans((ci % ci) * rowsumSE);
          // s_hat(0) += si(0, i);
          // s_hat(1) += det(si * trans(si));
          s_hat1(i) = si(0, i);
          s_hat2(i) = det(si * trans(si));
          mat onei(1, n, fill::zeros);
          onei(i) = 1;
          mat p = onei - si;
          qdiag_all.row(thread_id) += p % p;
        } catch (...) {
          flag_error = true;
        }
      }
    }
    if (flag_error) {
      throw exception("Matrix seems singular");
    } else {
      s_hat(0) = sum(s_hat1);
      s_hat(1) = sum(s_hat2);
      qdiag = mat(1, thread_nums, fill::ones) * qdiag_all;
      return List::create(
        Named("betas") = betas,
        Named("betas.SE") = betasSE,
        Named("s_hat") = s_hat,
        Named("q.diag") = qdiag
      );
    }
  } else {
    bool flag_error = false;
    for (int i = iStart; i < iEnd; i++) {
      if (!flag_error) {
        mat d = dm_given ? dmat.col(i) : gw_dist(dp, rp, i, p, theta, longlat, rp_given);
        mat w = gw_weight(d, bw, kernel, adaptive);
        mat ws(1, x.n_cols, fill::ones);
        mat xtw = trans(x %(w * ws));
        mat xtwx = xtw * x;
        mat xtwy = trans(x) * (w % y);
        try {
          mat xtwx_inv = inv(xtwx);
          betas.row(i) = trans(xtwx_inv * xtwy);
        } catch (...) {
          flag_error = true;
        }
      }
    }
    if (flag_error) {
      throw exception("Matrix seems singular.");
    } else {
      return List::create(
        Named("betas") = betas
      );
    }
  }
}
#endif
// CV approach for GWR
// [[Rcpp::export]]
double gw_cv_all(mat x, vec y, mat dp, bool dm_given, mat dmat, 
                 double p, double theta, bool longlat, 
                 double bw, int kernel, bool adaptive,
                 int ngroup, int igroup) {
  int n = dp.n_rows;
  double cv = 0.0;
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  for (int i = iStart; i < iEnd; i++) {
    mat d = dm_given ? dmat.col(i) : gw_dist(dp, dp, i, p, theta, longlat, false);
    mat w = gw_weight(d, bw, kernel, adaptive);
    w(i, 0) = 0.0;
    mat ws(1, x.n_cols, fill::ones);
    mat xtw = trans(x %(w * ws));
    mat xtwx = xtw * x;
    mat xtwy = trans(x) * (w % y);
    mat xtwx_inv = inv(xtwx);
    mat betas = xtwx_inv * xtwy;
    double res = y(i) - det(x.row(i) * betas);
    cv += res * res;
  }
  return cv;
}
//OpenMP CV function
// [[Rcpp::export]]
#ifdef _OPENMP
double gw_cv_all_omp(mat x, vec y, mat dp, bool dm_given, mat dmat, 
                     double p, double theta, bool longlat, 
                     double bw, int kernel, bool adaptive,
                     int threads, int ngroup, int igroup) {
  int n = dp.n_rows;
  int thread_nums = omp_get_num_procs() - 1;
  vec cv(thread_nums, fill::zeros);
  int lgroup = floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  bool flag_error = false;
#pragma omp parallel for num_threads(thread_nums)
  for (int i = iStart; i < iEnd; i++) {
    if (!flag_error) {
      int thread_id = threads > 0 ? threads : omp_get_thread_num();
      mat d = dm_given ? dmat.col(i) : gw_dist(dp, dp, i, p, theta, longlat, false);
      mat w = gw_weight(d, bw, kernel, adaptive);
      w(i, 0) = 0.0;
      mat ws(1, x.n_cols, fill::ones);
      mat xtw = trans(x %(w * ws));
      mat xtwx = xtw * x;
      mat xtwy = trans(x) * (w % y);
      try {
        mat xtwx_inv = inv(xtwx);
        mat betas = xtwx_inv * xtwy;
        double res = y(i) - det(x.row(i) * betas);
        cv(thread_id) += res * res;
      } catch (...) {
        flag_error = true;
      }
    }
  }
  if (flag_error) {
    throw exception("Matrix seems singular.");
  }
  return sum(cv);
}
#endif
//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------  Mixed GWR clalibration ------------------------------------------------------------------------------------
//gwr.mixed fast code
// [[Rcpp::export]]
vec e_vec(int m, int n){
  vec e = linspace(0, n-1, n);
  vec ret(n, fill::zeros);
  uvec u = find(e == m);
  ret.elem(u).fill(1);
  return ret;
}

// [[Rcpp::export]]
double gwr_mixed_trace(mat x1, mat x2, vec y, 
                       mat dMat, double bw, int kernel, bool adaptive){
  int i;
  int n = x1.n_rows;
  int nc2 = x2.n_cols;
  mat mtemp, model1, model2;
  mat x3(n, nc2);
  vec y2(n);
  vec y3;
  vec hii(n, fill::zeros);
  mat m(1,n);
  double s1, s2;
  for (i = 0; i < nc2; i++) {
    mtemp = gwr_q(x1, x2.col(i), dMat, bw, kernel, adaptive);
    x3.col(i) = x2.col(i) - fitted(x1, mtemp);
  }
  
  // The following works but is slow
  for (i = 0; i < n; i++) {
    mtemp = gwr_q(x1, e_vec(i, n), dMat, bw, kernel, adaptive); // length n x nc2
    y2 = e_vec(i, n) - fitted(x1, mtemp); // length n
    model2 = gwr_q(x3, y2, dMat, 100000, 4, true);
    y3 = e_vec(i, n) - fitted(x2, model2);
    model1 = gwr_q(x1, y3, dMat.col(i), bw, kernel, adaptive); // 1 x 1 matrix
    model2 = gwr_q(x3, y2, dMat.col(i), 100000, 4, true); // n x nc2
    s1 = fitted(x1.row(i), model1)(0);  // vector with one element
    s2 = fitted(x2.row(i), model2)(0);  // vector with one element
    hii(i) = s1 + s2;
  }
  return sum(hii);
}
// [[Rcpp::export]]
List gwr_mixed_2(mat x1, mat x2, vec y, 
                       mat dMat, mat dMat_rp,
                       double bw, int kernel, bool adaptive){
  int i;
  int n = x1.n_rows;
  int nc2 = x2.n_cols;
  mat mtemp, model1, model2;
  mat x3(n, nc2);
  vec y2(n);
  vec hii(n, fill::zeros);
  mat m(1,n);  
  for (i = 0; i < nc2; i++) {
    mtemp = gwr_q(x1, x2.col(i), dMat, bw, kernel, adaptive);
    x3.col(i) = x2.col(i) - fitted(x1, mtemp);
  }
  mtemp = gwr_q(x1, y, dMat, bw, kernel, adaptive);
  y2 = y - fitted(x1, mtemp);
  model2 = gwr_q(x3, y2, dMat, 100000, 4, true);
  
  model1 = gwr_q(x1, y-fitted(x2, model2), dMat_rp, bw, kernel, adaptive);
  model2 = gwr_q(x3, y2, dMat_rp, 100000, 4, true);
  
  return List::create(
    Named("local") = model1,
    Named("global") = model2
  );
}
//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------       Generalized GWR   ------------------------------------------------------------------------------------
//Possion GWR
/*
// [[Rcpp::export]]
List gwr_poission_wt(vec y, mat x, mat W, bool verbose){
  //accuracy control
  double tol = 1e-5;
  int maxiter = 200;
  int nr = x.n_rows;
  int nc = x.n_cols;
  mat s(nr, nr, fill::zeros);
  // model calibration
  int it_count = 0;
  double llik = 1.0;
  double old_llik = 1000000.0;
  vec mu = y+0.1;
  vec nu = log(mu);
  vec y_adj;
  vec wt2(nr);
  vec wi;
  mat beta(nr, nc);
  if(verbose){
    Rprintf(" Iteration    Log-Likelihood(With bandwidth: ");
    Rprintf("%f ", bw);
    Rprintf(")=========================\n");
    wt2.ones();
    while(abs((old_llik - llik)/llik)<tol || it_count >= maxiter){
      y_adj = nu + (y-mu)/mu
      for(int i=0; i<nr; i++){
        wi= W.col(i);
        beta.row(i) = gw_reg_2(x, y, wt2*wi);
      }
      nu = fitted(x, beta);
      mu = exp(nu);
      old_llik = llik;
      llik = sum(C_dpois(y, mu, log=true));
      if(verbose){
        Rprintf("          ");
        Rprintf("%f ", it_count);
        Rprintf("       ");
        Rprintf("%f ", llik);
        Rprintf("  \n");
      }
      wt2 = mu;
      it_count++;
    }
  }
  return List::create(
    Named("wt2") = wt2,
    Named("llik") = llik,
    Named("y.adj") = y_adj
  );
}
//Binomial GWR
List gwr_binomial_wt(vec y, mat x, mat W, bool verbose){
  //accuracy control
  double tol = 1e-5;
  int maxiter = 200;
  int nr = x.n_rows;
  int nc = x.n_cols;
  mat s(nr, nr, fill::zeros);
  // model calibration
  int it_count = 0;
  double llik = 1.0;
  double old_llik = 1000000.0;
  vec mu = y+0.1;
  vec nu = log(mu);
  vec y_adj;
  vec wt2(nr);
  vec wi;
  mat beta(nr, nc);
  if(verbose){
    Rprintf(" Iteration    Log-Likelihood(With bandwidth: ");
    Rprintf("%f ", bw);
    Rprintf(")=========================\n");
    wt2.ones();
    while(abs((old_llik - llik)/llik)<tol || it_count >= maxiter){
      y_adj = nu + (y-nr*mu)/(nr*mu*(1-mu));
      for(int i=0; i<nr; i++){
        wi= W.col(i);
        beta.row(i) = gw_reg_2(x, y, wt2*wi);
      }
      nu = fitted(x, beta);
      mu = exp(nu)/(1+exp(nu));
      old_llik = llik;
      llik = sum(C_dpois(y, mu, log=true));
      if(verbose){
        Rprintf("          ");
        Rprintf("%f ", it_count);
        Rprintf("       ");
        Rprintf("%f ", llik);
        Rprintf("  \n");
      }
      wt2 = mu;
      it_count++;
    }
  }
  return List::create(
    Named("wt2") = wt2,
    Named("llik") = llik,
    Named("y.adj") = y_adj
  );
}
*/
//---------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------       Other code        ------------------------------------------------------------------------------------



//--------------------------------------       Eigen code        ------------------------------------------------------------------------------------

enum cri
{
  CVR,
  dCVR,
};

enum kern
{
  exponential,
  gaussian,
  tricube,
  boxcar,
  bisquare,
};

enum appr
{
  apprAIC,
  apprBIC,
  apprCV,
};

double bw_gwr2(const Eigen::MatrixXd &x1, const Eigen::VectorXd &y, const Eigen::MatrixXd &dp_locat, appr approach, kern kernel, bool adaptive, const Eigen::MatrixXd &dMat, bool verbose, int nlower);
double gold(std::function<double(double, const Eigen::MatrixXd &, const Eigen::VectorXd &, kern, bool, const Eigen::MatrixXd &, double, double, bool, const Eigen::MatrixXd &, bool)> fun, double lower, double upper, bool adapt_bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose);
double gwr_aic(double bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose);
double gwr_bic(double bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose);
double  gwr_cv(double bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose);
void  gw_reg_all(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &dp_locat, bool rp_given, const Eigen::MatrixXd &rp_locat, bool dm_given, const Eigen::MatrixXd &dMat, bool hatmatrix, double p, double theta, bool longlat, double bw, kern kernel, bool adaptive, Eigen::MatrixXd &betas, Eigen::MatrixXd &s_hat);
double gw_cv_all(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &dp_locat, bool dm_given, const Eigen::MatrixXd &dMat, double p, double theta, bool longlat, double bw, kern kernel, bool adaptive, int ngroup = 1, int igroup = 1);
Eigen::MatrixXd gw_dist(const Eigen::MatrixXd &dp, const Eigen::MatrixXd &rp, int focus, double p, double theta, bool longlat, bool rp_given);
Eigen::MatrixXd coordinate_rotate(const Eigen::MatrixXd &coords, double theta);
Eigen::MatrixXd gw_weight(const Eigen::MatrixXd &d, double bw, kern kernel, bool adaptive);

/*
double gw_weight_boxcar(double dist, double bw);      // Already defined?
double gw_weight_tricube(double dist, double bw);     // Already defined?
double gw_weight_bisquare(double dist, double bw);    // Already defined?
double gw_weight_exponential(double dist, double bw); // Already defined?
double gw_weight_gaussian(double dist, double bw);    // Already defined?
*/

double aic_c1(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, const Eigen::MatrixXd &betas, const Eigen::MatrixXd &s_hat);
double    bic(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, const Eigen::MatrixXd &betas, const Eigen::MatrixXd &s_hat);
void gwr_q2(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &loc, bool adaptive, bool hatmatrix, double bw, kern kernel, const Eigen::MatrixXd &dMat, const Eigen::MatrixXd &dMatSorted, int wt2, Eigen::MatrixXd &betas, Eigen::MatrixXd &S, Eigen::MatrixXd C[]);
Eigen::MatrixXd gw_weight_vec(const Eigen::VectorXd &col, const Eigen::VectorXd &colSorted, double bw, kern kernel, bool adaptive);
void gw_reg(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &w, bool hatmatrix, int focus, Eigen::MatrixXd &Ci, Eigen::MatrixXd &Si, Eigen::MatrixXd &betasi);
Eigen::VectorXd gw_fitted(const Eigen::MatrixXd &X, const Eigen::MatrixXd &beta);
void trhat2(const Eigen::MatrixXd &S, double &s_hat0, double &s_hat1);
void gwr_diag(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, const Eigen::MatrixXd &betas, const Eigen::MatrixXd &S, Rcpp::DoubleVector &result);

Eigen::VectorXd sp_dists(const Eigen::MatrixXd &dp, const Eigen::VectorXd &loc);
/* double sp_gcdist(double lon1, double lon2, double lat1, double lat2); */ // Already defined?

Eigen::VectorXd mk_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd &out_loc, double p);
Eigen::MatrixXd mk_dist_smat(const Eigen::MatrixXd &in_locs, double p);
Eigen::MatrixXd mk_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs, double p);
Eigen::VectorXd cd_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd out_loc);
Eigen::MatrixXd cd_dist_smat(const Eigen::MatrixXd &in_locs);
Eigen::MatrixXd cd_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs);
Eigen::VectorXd md_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd &out_loc);
Eigen::MatrixXd md_dist_smat(const Eigen::MatrixXd &in_locs);
Eigen::MatrixXd md_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs);
Eigen::VectorXd eu_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd &out_loc);
Eigen::MatrixXd eu_dist_smat(const Eigen::MatrixXd &in_locs);
Eigen::MatrixXd eu_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs);

// [[Rcpp::export]]
Rcpp::List new_multiscale(Eigen::MatrixXd x, Eigen::MatrixXd x1,
                          Rcpp::List dMatsParam, Eigen::MatrixXd dp_locat,
                          Eigen::VectorXd y, 
                          Eigen::Matrix<double, Eigen::Dynamic, 1> bws0, 
                          std::vector<int> var_dMat_index, bool adaptive, 
                          bool verbose, int nlower, bool hatmatrix, 
                          int max_iterations, double threshold, int max_threads,
                          std::vector<std::string> variable_names, int kerneln,
                          int approachn, int crin, int bws_reOpts)
{
  // Convert R indices (one-indexed) to C++ indices (zero-indexed)
  for (int i = 0; i < var_dMat_index.size(); i++) {
    var_dMat_index[i]--;
  }
  // Convert list of dMat matrices to vector of Eigen matrices
  // Taken from https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-March/005511.html
  std::vector<Eigen::MatrixXd> dMats(dMatsParam.size());
  for (int i = 0; i < dMatsParam.size(); i++) {
    Rcpp::NumericMatrix tmpcv = dMatsParam[i];
    double *pcv = &tmpcv(0,0);
    
    Eigen::Map<Eigen::MatrixXd> tmpmapd(pcv, tmpcv.nrow(), tmpcv.ncol());
    dMats[i] = Eigen::MatrixXd(tmpmapd.cast<double>());
  }
#ifdef _OPENMP
  int max_threads_checked = omp_get_num_procs() - 1;
  if (max_threads > 0) {
    max_threads_checked = max_threads;
  }
  omp_set_num_threads(max_threads_checked);
#endif
  
  Eigen::MatrixXd dMat, dMatSorted;
  // Create sorted dMats
  std::vector<Eigen::MatrixXd> dMatsSorted(dMats.size());
  for (int d = 0; d < dMats.size(); d++) {
    dMatSorted = dMats[d];
#pragma omp parallel for
    for (int i = 0; i < dMatSorted.cols(); i++)
    {
      std::sort(dMatSorted.col(i).begin(), dMatSorted.col(i).end());
    }
    dMatsSorted[d] = dMatSorted;
  }
  
  kern kernel = kern::bisquare;
  switch (kerneln) {
    case 0:
      kernel = kern::gaussian;
      break;
    case 1:
      kernel = kern::exponential;
      break;
    case 2:
      kernel = kern::bisquare;
      break;
    case 3:
      kernel = kern::tricube;
      break;
    case 4:
      kernel = kern::boxcar;
      break;
  }
  appr approach = appr::apprAIC;
  switch (approachn) {
    case 0:
      approach = appr::apprAIC;
      break;
    case 1:
      approach = appr::apprBIC;
      break;
    case 2:
      approach = appr::apprCV;
      break;
  }
  cri criterion = cri::dCVR;
  switch (crin) {
    case 0:
      criterion = cri::dCVR;
      break;
    case 1:
      criterion = cri::CVR;
      break;
  }
  
  Rcpp::Rcout << "Calculate the initial beta0 from the above bandwidths.\n";
  dMat = dMats[0];
  double bw_int0 = bw_gwr2(x1, y, dp_locat, approach, kernel, adaptive, dMat, verbose, nlower);
  int dp_n = dp_locat.rows(); 
  int var_n = variable_names.size();
  
  Eigen::MatrixXd *S_arrays = new Eigen::MatrixXd[var_n];
  for (int i = 0; i < var_n; i++)
  {
    S_arrays[i] = Eigen::MatrixXd::Zero(dp_n, dp_n);
  }
  Eigen::MatrixXd Beta_SE, Beta_TV, betas, S_hat;
  if (hatmatrix)
  {
    Beta_SE = Eigen::MatrixXd(dp_n, var_n);
    Beta_TV = Eigen::MatrixXd(dp_n, var_n);
  }
  Eigen::MatrixXd *C = new Eigen::MatrixXd[dp_n];
  gwr_q2(x1, y, dp_locat, adaptive, hatmatrix, bw_int0, kernel, dMat, dMatSorted, 1, betas, S_hat, C);
  
  if (hatmatrix) {
    for (int i = 0; i < var_n; i++)
    {
      for (int j = 0; j < dp_n; j++)
      {
        S_arrays[i].row(j) = x1(j, i) * C[j].row(i);
      }
    }
  }
  Rcpp::Rcout << "End of calculating the inital beta0.\n";
  Rcpp::Rcout << "Find the optimum bandwidths for each independent variable.\n";
  
  int iteration = 0;
  double **bws_vars = new double*[max_iterations+1];
  for (int i = 0; i < max_iterations+1; i++) {
    bws_vars[i] = new double[var_n];
    for (int j = 0; j < var_n; j++) {
      bws_vars[i][j] = -1.0;
    }
  }
  for (int i = 0; i < var_n; i++) {
    bws_vars[iteration][i] = bws0(i);
  }
  double criterion_val = 2.0*threshold;
  Eigen::VectorXd resid_i = y - gw_fitted(x1, betas);
  double rss0 = resid_i.transpose() * resid_i;
  double rss1 = 0.;
  auto *rss_vals = new double[max_iterations][3];
  
  bool *bw_seled = new bool[var_n];
  int *bws_change_NO = new int[var_n];
  double *bws_thresholds = new double[var_n];
  for (int i = 0; i < var_n; i++)
  {
    bw_seled[i] = false;
    bws_change_NO[i] = 0;
    bws_thresholds[i] = 0.1;
  }
  Eigen::MatrixXd beta_i, S_i;
  while ((iteration < max_iterations) && (criterion_val > threshold))
  {
    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "\tIteration " << iteration + 1 << "\n";
    for (int i = 0; i < var_n; i++)
    {
      dMat = dMats[var_dMat_index[i]];
      dMatSorted = dMatsSorted[var_dMat_index[i]];
      Eigen::MatrixXd x1_i = x1.col(i);
      Eigen::VectorXd f_i = betas.col(i).cwiseProduct(x1_i);
      Eigen::VectorXd y_i = resid_i + f_i;
      double bw_i;
      if (bw_seled[i])
      {
        bw_i = bws0(i);
      }
      else
      {
        Rcpp::Rcout << "Now select an optimum bandwidth for the variable: " << variable_names[i] << "\n";
        bw_i = bw_gwr2(x1_i, y_i, dp_locat, approach, kernel, adaptive, dMat, verbose, nlower);
        double diff = std::abs(bw_i - bws0(i));
        Rcpp::Rcout << "The newly selected bandwidth for variable: " << variable_names[i] << " is " << bw_i  << "\n";
        Rcpp::Rcout << "The bandwidth used in the last iteration is: " << bws0(i) << " and the difference in bandwidths is: " << diff << '\n';
        if (diff > bws_thresholds[i])
        {
          Rcpp::Rcout << "The bandwidth for variable " << variable_names[i] << " will be continually selected in the next iteration.\n";
          bws_change_NO[i] = 0;
        }
        else
        {
          bws_change_NO[i]++;
          if (bws_change_NO[i] < bws_reOpts)
          {
            Rcpp::Rcout << "The bandwidth for variable " << variable_names[i] << " seems to be converged for " << bws_change_NO[i] << " times. It will be continually optimized in the next " << bws_reOpts - bws_change_NO[i] << " times\n";
          }
          else // They are equal
          {
            Rcpp::Rcout << "The bandwidth for variable " << variable_names[i] << " seems to be converged and will be kept the same in the following iterations.\n";
            bw_seled[i] = true;
          }
        }
      }
      bws0(i) = bw_i;
      bws_vars[iteration+1][i] = bw_i;
      
      gwr_q2(x1_i, y_i, dp_locat, adaptive, hatmatrix, bw_i, kernel, dMat, dMatSorted, -1, beta_i, S_i, C);
      
      if (hatmatrix)
      {
        Eigen::MatrixXd S_array_i = S_arrays[i];
        S_arrays[i] = S_i * S_array_i + S_i - S_i * S_hat;
        S_hat = S_hat - S_array_i + S_arrays[i];
      }
      
      betas.col(i) = beta_i;
      resid_i = y - gw_fitted(x1, betas);
    }
    Eigen::VectorXd y_min_gw = y - gw_fitted(x1, betas);
    rss1 = y_min_gw.transpose() * y_min_gw;
    if (criterion == cri::CVR)
    {
      criterion_val = std::abs(rss1 - rss0);
      Rcpp::Rcout << "Iteration " << iteration + 1 << " change of RSS (CVR) is " << criterion_val << ".\n";
    }
    else
    {
      criterion_val = std::sqrt(std::abs(rss1 - rss0) / rss1);
      Rcpp::Rcout << "Iteration " << iteration + 1 << " change of RSS (dCVR) is " << criterion_val << ".\n";
    }
    rss0 = rss1;
    rss_vals[iteration][0] = rss0;
    rss_vals[iteration][1] = rss1;
    rss_vals[iteration][2] = criterion_val;
    iteration++;
  }
  Eigen::VectorXd yhat = gw_fitted(x1,betas);
  resid_i = y - yhat;
  Rcpp::DoubleVector S_arrays_return = Rcpp::DoubleVector(Rcpp::Dimension(var_n*dp_n*dp_n));
  long count = 0;
  for (int i = 0; i < var_n; i++) {
    for (int j = 0; j < dp_n; j++) {
      for (int k = 0; k < dp_n; k++) {
        S_arrays_return(count) = S_arrays[i](j,k);
        count++;
      }
    }
  }
  
  Rcpp::DoubleVector bws_vars_ret = Rcpp::DoubleVector(Rcpp::Dimension(iteration+1,var_n));
  for (int i = 0; i < iteration+1; i++) {
    for (int j = 0; j < var_n; j++) {
      bws_vars_ret(i,j) = bws_vars[i][j];
    }
  }
  
  Rcpp::DoubleVector gwr_diag_result = Rcpp::DoubleVector(Rcpp::Dimension(10));
  
  if (hatmatrix) {
    gwr_diag(y, x1, betas, S_hat, gwr_diag_result);
    
    // Calculate the SEs and t-values
    double RSS_gw = gwr_diag_result(4);
    double tr_Shat = gwr_diag_result(7);
    double sigma_hat11 = RSS_gw/(dp_n-tr_Shat);
    
    Eigen::MatrixXd xcoli_mat = Eigen::MatrixXd::Zero(dp_n, dp_n);
    for (int i = 0; i < var_n; i++)
    {
      xcoli_mat.diagonal() = x.col(i).cwiseInverse();
      Eigen::MatrixXd Ci = xcoli_mat * S_arrays[i];
      Beta_SE.col(i) = (Ci * Ci.transpose() * sigma_hat11).diagonal().cwiseSqrt();
      Beta_TV.col(i) = betas.col(i).array() / Beta_SE.col(i).array();
    }
  }
  
  delete [] bws_thresholds;
  delete [] bws_change_NO;
  delete [] bw_seled;
  delete [] rss_vals;
  delete [] C;
  delete [] S_arrays;
  for (int i = 0; i < max_iterations+1; i++) {
    delete [] bws_vars[i];
  }
  delete [] bws_vars;
  
  return Rcpp::List::create(Rcpp::Named("yhat", yhat),
                            Rcpp::Named("residual", resid_i),
                            Rcpp::Named("betas", betas),
                            Rcpp::Named("Shat", S_hat),
                            Rcpp::Named("S.arrays", S_arrays_return),
                            Rcpp::Named("Beta_SE", Beta_SE),
                            Rcpp::Named("Beta_TV", Beta_TV),
                            Rcpp::Named("bws0", bws0),
                            Rcpp::Named("bws.vars", bws_vars_ret),
                            Rcpp::Named("mgwr.diag", gwr_diag_result));
}

double bw_gwr2(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &dp_locat, appr approach, kern kernel, bool adaptive, const Eigen::MatrixXd &dMat, bool verbose, int nlower)
{
  int dp_n = dp_locat.rows();
  double upper, lower;
  if (adaptive)
  {
    upper = dp_n;
    lower = nlower;
  }
  else
  {
    upper = dMat.maxCoeff();
    lower = upper / 5000;
  }
  
  std::function<double(double, const Eigen::MatrixXd &, const Eigen::VectorXd &, kern, bool, const Eigen::MatrixXd &, double, double, bool, const Eigen::MatrixXd &, bool)> fun;
  if (approach == appr::apprBIC) { fun = &gwr_bic; }
  else if (approach == appr::apprAIC) { fun = &gwr_aic; }
  else { fun = &gwr_cv; } // CV
  return gold(fun, lower, upper, adaptive, x, y, kernel, adaptive, dp_locat, 2, 0, false, dMat, verbose);
}

double gold(std::function<double(double, const Eigen::MatrixXd &, const Eigen::VectorXd &, kern, bool, const Eigen::MatrixXd &, double, double, bool, const Eigen::MatrixXd &, bool)> fun, double lower, double upper, bool adapt_bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose)
{
  double eps = 1e-4;
  double R = (std::sqrt(5) - 1) / 2;
  double d = R * (upper - lower);
  double x1 = lower + d;
  double x2 = upper - d;
  if (adapt_bw)
  {
    x1 = std::floor(x1);
    x2 = std::round(x2);
  }
  double f1 = fun(x1, x, y, kernel, adaptive, dp_locat, p, theta, longlat, dMat, verbose);
  double f2 = fun(x2, x, y, kernel, adaptive, dp_locat, p, theta, longlat, dMat, verbose);
  double d1 = f2 - f1;
  double xopt = 0.;
  if (f1 < f2)
  {
    xopt = x1;
  }
  else
  {
    xopt = x2;
  }
  while ((std::abs(d) > eps) && (std::abs(d1) > eps))
  {
    d = R * d;
    if (f1 < f2)
    {
      lower = x2;
      x2 = x1;
      x1 = lower + d;
      if (adapt_bw)
        x1 = std::round(lower + d);
      f2 = f1;
      f1 = fun(x1, x, y, kernel, adaptive, dp_locat, p, theta, longlat, dMat, verbose);
    }
    else
    {
      upper = x1;
      x1 = x2;
      x2 = upper - d;
      if (adapt_bw)
        x2 = std::floor(upper - d);
      f1 = f2;
      f2 = fun(x2, x, y, kernel, adaptive, dp_locat, p, theta, longlat, dMat, verbose);
    }
    if (f1 < f2)
      xopt = x1;
    else
      xopt = x2;
    d1 = f2 - f1;
  }
  return xopt;
}

double gwr_aic(double bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose)
{
  bool DM_given = true;
  Eigen::MatrixXd betas, s_hat;
  gw_reg_all(x, y, dp_locat, false, dp_locat, DM_given, dMat, true, p, theta, longlat, bw, kernel, adaptive, betas, s_hat);
  double aicc_value = aic_c1(y, x, betas, s_hat);
  if (std::isnan(aicc_value)) { aicc_value = std::numeric_limits<double>::infinity(); }
  if (verbose)
  {
    if (adaptive)
    {
      Rcpp::Rcout << "Adaptive bandwidth (number of nearest neighbours): " << (int)bw << " AICc value: " << aicc_value << ".\n";
    }
    else
    {
      Rcpp::Rcout << "Fixed bandwidth: " << bw << " AICc value: " << aicc_value << ".\n";
    }
  }
  return aicc_value;
}

double gwr_bic(double bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose)
{
  bool DM_given = true;
  Eigen::MatrixXd betas, s_hat;
  gw_reg_all(x, y, dp_locat, false, dp_locat, DM_given, dMat, true, p, theta, longlat, bw, kernel, adaptive, betas, s_hat);
  double bic_value = bic(y, x, betas, s_hat);
  if (std::isnan(bic_value)) { bic_value = std::numeric_limits<double>::infinity(); }
  if (verbose)
  {
    if (adaptive)
    {
      Rcpp::Rcout << "Adaptive bandwidth (number of nearest neighbours): " << (int)bw << " BIC value: " << bic_value << ".\n";
    }
    else
    {
      Rcpp::Rcout << "Fixed bandwidth: " << bw << " BIC value: " << bic_value << ".\n";
    }
  }
  return bic_value;
}

double gwr_cv(double bw, const Eigen::MatrixXd &x, const Eigen::VectorXd &y, kern kernel, bool adaptive, const Eigen::MatrixXd &dp_locat, double p, double theta, bool longlat, const Eigen::MatrixXd &dMat, bool verbose)
{
  bool DM_given = true;
  Eigen::MatrixXd betas, s_hat;
  double cv_value = gw_cv_all(x, y, dp_locat, DM_given, dMat, p, theta, longlat, bw, kernel, adaptive);
  if (std::isnan(cv_value)) { cv_value = std::numeric_limits<double>::infinity(); }
  if (verbose)
  {
    if (adaptive)
    {
      Rcpp::Rcout << "Adaptive bandwidth (number of nearest neighbours): " << (int)bw << " CV value: " << cv_value << ".\n";
    }
    else
    {
      Rcpp::Rcout << "Fixed bandwidth: " << bw << " CV value: " << cv_value << ".\n";
    }
  }
  return cv_value;
}

void gw_reg_all(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &dp_locat, bool rp_given, const Eigen::MatrixXd &rp_locat, bool dm_given, const Eigen::MatrixXd &dMat, bool hatmatrix, double p, double theta, bool longlat, double bw, kern kernel, bool adaptive, Eigen::MatrixXd &betas, Eigen::MatrixXd &s_hat)
{
  int n = rp_locat.rows();
  int k = x.cols();
  betas = Eigen::MatrixXd::Zero(n, k);
  if (hatmatrix)
  {
    Eigen::MatrixXd loc_s_hat = Eigen::MatrixXd::Zero(n, 2);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      // TODO check again (use col instead?)
      //Eigen::MatrixXd d = dMat.row(i).transpose();
      Eigen::MatrixXd d = dm_given ? dMat.col(i) : gw_dist(dp_locat, rp_locat, i, p, theta, longlat, rp_given);
      Eigen::MatrixXd w = gw_weight(d, bw, kernel, adaptive);
      Eigen::MatrixXd xtw = x.cwiseProduct(w * Eigen::MatrixXd::Constant(1, k, 1)).transpose();
      Eigen::MatrixXd xtwx = xtw * x;
      Eigen::MatrixXd xtwy = x.transpose() * (w.cwiseProduct(y));
      // NOTE: More stable with ColPivHouseholderQR
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp = xtwx.colPivHouseholderQr();
      betas.row(i) = decomp.solve(xtwy).transpose();
      Eigen::MatrixXd ci = decomp.solve(xtw);
      Eigen::VectorXd si = x.row(i) * ci;
      // betasSE
      loc_s_hat(i, 0) = si(i);
      // NOTE: The resulting value does not seem to be used.
      // s_hat(0,1) += (si * si.transpose()).determinant();
      loc_s_hat(i, 1) = si.transpose() * si;
      Eigen::VectorXd onei = Eigen::VectorXd::Zero(n);
      onei(i) = 1;
      Eigen::VectorXd pp = onei - si;
      // qdiag += pp.cwiseProduct(pp);
    }
    s_hat = loc_s_hat.colwise().sum();
  }
  else
  {
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
      Eigen::MatrixXd d = dm_given ? dMat.col(i) : gw_dist(dp_locat, dp_locat, i, p, theta, longlat, rp_given);
      Eigen::MatrixXd w = gw_weight(d, bw, kernel, adaptive);
      Eigen::MatrixXd xtw = x.cwiseProduct(w * Eigen::MatrixXd::Constant(1, x.cols(), 1)).transpose();
      Eigen::MatrixXd xtwx = xtw * x;
      Eigen::MatrixXd xtwy = x.transpose() * (w.cwiseProduct(y));
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp = xtwx.colPivHouseholderQr();
      Eigen::VectorXd betasi = decomp.solve(xtwy);
      betas.row(i) = betasi.transpose();
    }
  }
}

double gw_cv_all(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &dp_locat, bool dm_given, const Eigen::MatrixXd &dMat, double p, double theta, bool longlat, double bw, kern kernel, bool adaptive, int ngroup, int igroup) 
{
  igroup--;
  int n = dp_locat.rows();
  double cv = 0.0;
  int lgroup = std::floor(((double)n) / ngroup);
  int iStart = igroup * lgroup, iEnd = (igroup + 1 < ngroup) ? (igroup + 1) * lgroup : n;
  for (int i = iStart; i < iEnd; i++) {
    Eigen::MatrixXd d = dm_given ? dMat.col(i) : gw_dist(dp_locat, dp_locat, i, p, theta, longlat, false);
    Eigen::MatrixXd w = gw_weight(d, bw, kernel, adaptive);
    w(i, 0) = 0.0;
    Eigen::MatrixXd xtw = x.cwiseProduct(w * Eigen::MatrixXd::Constant(1, x.cols(), 1)).transpose();
    Eigen::MatrixXd xtwx = xtw * x;
    Eigen::MatrixXd xtwy = x.transpose() * (w.cwiseProduct(y));
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp = xtwx.colPivHouseholderQr();
    Eigen::VectorXd betas = decomp.solve(xtwy);
    double res = y(i) - (x.row(i) * betas);
    cv += res * res;
  }
  return cv;
}

Eigen::MatrixXd gw_dist(const Eigen::MatrixXd &dp_orig, const Eigen::MatrixXd &rp_orig, int focus, double p, double theta, bool longlat, bool rp_given) {
  int ndp = dp_orig.rows(), nrp = rp_orig.rows();
  int isFocus = focus > -1;
  Eigen::MatrixXd dp, rp;
  if (p != 2 && theta != 0 && !longlat) {
    dp = coordinate_rotate(dp, theta);
    rp = coordinate_rotate(rp, theta);
  }
  else {
    dp = dp_orig;
    rp = rp_orig;
  }
  if (isFocus) {
    Eigen::VectorXd prp = rp.row(focus).transpose();
    if (longlat) {
      return sp_dists(dp, prp);
    } else {
      if (p == 2.0)
        return eu_dist_vec(dp, prp);
      else if(p == -1.0)
        return cd_dist_vec(dp, prp);
      else if(p == 1.0)
        return md_dist_vec(dp, prp);
      else
        return mk_dist_vec(dp, prp, p);
    }
  } else {
    if (longlat) {
      Eigen::MatrixXd dists = Eigen::MatrixXd::Zero(ndp, nrp);
      for (int i = 0; i < nrp; i++) {
        dists.col(i) = sp_dists(dp, rp.row(i).transpose());
      }
      return dists.transpose();
    } else {
      if (p == 2.0)
        return rp_given ? eu_dist_mat(dp, rp) : eu_dist_smat(dp);
      else if (p == -1.0)
        return rp_given ? cd_dist_mat(dp, rp) : cd_dist_smat(dp);
      else if (p == 1.0)
        return rp_given ? md_dist_mat(dp, rp) : md_dist_smat(dp);
      else
        return rp_given ? mk_dist_mat(dp, rp, p) : mk_dist_smat(dp, p);
    }
  }
}

Eigen::MatrixXd coordinate_rotate(const Eigen::MatrixXd &coords, double theta)
{
  int n = coords.rows();
  Eigen::MatrixXd rotated_coords = Eigen::MatrixXd::Zero(n, 2);
  rotated_coords.col(0) = coords.col(0) * std::cos(theta) - coords.col(1) * std::sin(theta);
  rotated_coords.col(1) = coords.col(0) * std::sin(theta) + coords.col(1) * std::cos(theta);
  return rotated_coords;
}


Eigen::MatrixXd gw_weight(const Eigen::MatrixXd &d, double bw, kern kernel, bool adaptive)
{
  std::function<double(double, double)> fun;
  switch (kernel) {
    case kern::exponential:
      fun = &gw_weight_exponential;
      break;
    case kern::gaussian:
      fun = &gw_weight_gaussian;
      break;
    case kern::tricube:
      fun = &gw_weight_tricube;
      break;
    case kern::boxcar:
      fun = &gw_weight_boxcar;
      break;
    case kern::bisquare:
      fun = &gw_weight_bisquare;
      break;
  }
  int nr = d.rows(), nc = d.cols();
  Eigen::MatrixXd w = Eigen::MatrixXd::Zero(nr, nc);
  if (adaptive)
  {
    double dn = bw / (double)nr;
    for (int c = 0; c < nc; c++)
    {
      int i = -1;
      double fixbw;
      if (dn <= 1)
      {
        // NOTE: MAYBE BAD?????
        Eigen::VectorXd vdist = d.col(c);
        std::sort(vdist.begin(), vdist.end());
        i = floor(bw) - 1; // This could maybe be -1 and we have an error
        if (i <= 0)
        {
          i = 0;
        }
        fixbw = vdist(i);
      }
      else
      {
        fixbw = dn * d.col(c).maxCoeff();
      }
      for (int r = 0; r < nr; r++)
      {
        w(r, c) = fun(d(r, c), fixbw);
      }
    }
  }
  else
  {
    for (int c = 0; c < nc; c++) {
      for (int r = 0; r < nr; r++) {
        w(r, c) = fun(d(r, c), bw);
      }
    }
  }
  return w;
}

/*
double gw_weight_gaussian(double dist, double bw)
{
  return std::exp(std::pow(dist, 2)/((-2)*std::pow(bw, 2)));
}

double gw_weight_exponential(double dist, double bw)
{
  return std::exp(-dist/bw);
}

double gw_weight_bisquare(double dist, double bw)
{
  if (bw <= 0. || dist >= bw)
  {
    return 0;
  }
  double a = 1 - (dist * dist) / (bw * bw);
  return a*a;
}

double gw_weight_tricube(double dist, double bw)
{
  return dist > bw ? 0 : std::pow(1 - std::pow(dist, 3)/std::pow(bw, 3), 3);
}

double gw_weight_boxcar(double dist, double bw)
{
  return dist > bw ? 0.0 : 1.0;
}
*/

double aic_c1(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, const Eigen::MatrixXd &betas, const Eigen::MatrixXd &s_hat)
{
  // Residual sum
  Eigen::VectorXd temp_ss = y - betas.cwiseProduct(x).rowwise().sum();
  double ss = temp_ss.transpose() * temp_ss;
  int n = x.rows();
  return n * std::log(ss / n) + n * std::log(2 * M_PI) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0)));
}

double bic(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, const Eigen::MatrixXd &betas, const Eigen::MatrixXd &s_hat)
{
  // Residual sum
  Eigen::VectorXd temp_ss = y - betas.cwiseProduct(x).rowwise().sum();
  double ss = temp_ss.transpose() * temp_ss;
  int n = x.rows();
  return n * std::log(ss / n) + n * std::log(2 * M_PI) + std::log(n) * s_hat(0);
}

void gwr_q2(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &loc, bool adaptive, bool hatmatrix, double bw, kern kernel, const Eigen::MatrixXd &dMat, const Eigen::MatrixXd &dMatSorted, int wt2, Eigen::MatrixXd &betas, Eigen::MatrixXd &S, Eigen::MatrixXd C[])
{
  int dp_n = loc.rows();
  int var_n = x.cols();
  betas = Eigen::MatrixXd::Zero(dp_n, var_n);
  S = Eigen::MatrixXd::Zero(dp_n, dp_n);
#pragma omp parallel for
  for (int i = 0; i < dMat.cols(); i++)
  {
    Eigen::MatrixXd Ci, Si, betasi;
    Eigen::MatrixXd W = gw_weight_vec(dMat.col(i), dMatSorted.col(i), bw, kernel, adaptive);
    if (wt2 != -1)
    {
      W *= wt2;
    }
    gw_reg(x, y, W, hatmatrix, i, Ci, Si, betasi);
    betas.row(i) = betasi.transpose();
    if (hatmatrix) {
      S.row(i) = Si;
      C[i] = Ci;
    }
  }
}

Eigen::MatrixXd gw_weight_vec(const Eigen::VectorXd &col, const Eigen::VectorXd &colSorted, double bw, kern kernel, bool adaptive)
{
  std::function<double(double, double)> kernfunc;
  switch (kernel) {
    case kern::exponential:
      kernfunc = &gw_weight_exponential;
      break;
    case kern::gaussian:
      kernfunc = &gw_weight_gaussian;
      break;
    case kern::tricube:
      kernfunc = &gw_weight_tricube;
      break;
    case kern::boxcar:
      kernfunc = &gw_weight_boxcar;
      break;
    case kern::bisquare:
      kernfunc = &gw_weight_bisquare;
      break;
  }
  int n = col.rows() * col.cols();
  Eigen::MatrixXd wv = Eigen::MatrixXd::Zero(n, 1);
  double dn = bw / (double)n;
  if (adaptive)
  {
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      double fixbw;
      if (dn <= 1)
      {
        int index = std::floor(bw) - 1; // NOTE: This could be -1 !!!
        if (index < 0)
        {
          index = 0;
        }
        fixbw = colSorted(index);
      }
      else
      {
        fixbw = dn * colSorted.maxCoeff();
      }
      wv(i, 0) = kernfunc(col(i), fixbw);
    }
  }
  else
  {
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
      wv(i, 0) = kernfunc(col(i), bw);
    }
  }
  return wv;
}

void gw_reg(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &w, bool hatmatrix, int focus, Eigen::MatrixXd &Ci, Eigen::MatrixXd &Si, Eigen::MatrixXd &betasi)
{
  Eigen::MatrixXd xtw = x.cwiseProduct(w * Eigen::MatrixXd::Constant(1, x.cols(), 1)).transpose();
  Eigen::MatrixXd xtwx = xtw * x;
  Eigen::MatrixXd xtwy = x.transpose() * (w.cwiseProduct(y));
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp = xtwx.colPivHouseholderQr();
  betasi = decomp.solve(xtwy);
  if (hatmatrix)
  {
    Ci = decomp.solve(xtw);
    Si = x.row(focus) * Ci;
  }
}

Eigen::VectorXd gw_fitted(const Eigen::MatrixXd &X, const Eigen::MatrixXd &beta)
{
  return beta.cwiseProduct(X).rowwise().sum();
}

void trhat2(const Eigen::MatrixXd &S, double &s_hat0, double &s_hat1)
{
  int n_obs = S.rows();
  double htr2 = 0.0;
  for (int i = 0; i < n_obs; i++)
  {
    htr2 += S.row(i) * S.row(i).transpose();
  }
  s_hat0 = S.trace();
  s_hat1 = htr2;
}

void gwr_diag(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, const Eigen::MatrixXd &betas, const Eigen::MatrixXd &S, Rcpp::DoubleVector &result) {
  Eigen::VectorXd r = y - gw_fitted(x, betas);
  double ss = r.transpose() * r;
  double s_hat0, s_hat1;
  trhat2(S, s_hat0, s_hat1);
  double n = (double) S.rows();
  double AIC = n * std::log(ss/n) + n * std::log(2*M_PI) + n + s_hat0;
  double AICc = n * std::log(ss/n) + n * std::log(2*M_PI) + n * ((n+s_hat0) / (n-2-s_hat0));
  double edf = n - 2 * s_hat0 + s_hat1;
  double enp = 2 * s_hat0 - s_hat1;
  Eigen::VectorXd y_min_ymean = y.array() - y.mean();
  double yss = y_min_ymean.transpose() * y_min_ymean;
  double r2 = 1 - ss/yss;
  double r2_adj = 1 - (1-r2) * (n-1) / (edf-1);
  double BIC = n * std::log(ss / n) + n * std::log(2.0 * M_PI) + std::log(n) * s_hat0;
  
  result(0) = AIC;
  result(1) = AICc;
  result(2) = edf;
  result(3) = enp;
  result(4) = ss;
  result(5) = r2;
  result(6) = r2_adj;
  result(7) = s_hat0;
  result(8) = s_hat1;
  result(9) = BIC;
}

Eigen::MatrixXd eu_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs)
{
  int n_in = in_locs.rows();
  int n_out = out_locs.rows();
  Eigen::MatrixXd eu_dist = Eigen::MatrixXd::Zero(n_in, n_out);
  int i = 0, j = 0;
  for (i = 0; i < n_in; i++)
  {
    for (j = 0; j < n_out; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - out_locs.row(j)).transpose();
      eu_dist(i,j) = r.transpose() * r;
    }
  }
  return eu_dist.cwiseSqrt();
}

Eigen::MatrixXd eu_dist_smat(const Eigen::MatrixXd &in_locs)
{
  int n = in_locs.rows();
  Eigen::MatrixXd eu_dist = Eigen::MatrixXd::Zero(n, n);
  for (int k = 0; k < n * n; k++)
  {
    int i = k / n, j = k % n;
    Eigen::VectorXd r = (in_locs.row(i) - in_locs.row(j)).transpose();
    eu_dist(i, j) = r.transpose() * r;
    eu_dist(j, i) = eu_dist(i, j);
  }
  return eu_dist.cwiseSqrt();
}

Eigen::VectorXd eu_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd &out_loc)
{
  int n_in = in_locs.rows();
  Eigen::VectorXd eu_dist = Eigen::VectorXd::Zero(n_in);
  for (int i = 0; i < n_in; i++)
  {
    Eigen::VectorXd r = (in_locs.row(i) - out_loc.transpose()).transpose();
    eu_dist(i) = r.transpose() * r;
  }
  return eu_dist.cwiseSqrt();
}

Eigen::MatrixXd md_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs)
{
  int n_in = in_locs.rows();
  int n_out = out_locs.rows();
  Eigen::MatrixXd md_dist = Eigen::MatrixXd::Zero(n_in, n_out);
  for (int i = 0; i < n_in; i++)
  {
    for (int j = 0; j < n_out; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - out_locs.row(j)).transpose().cwiseAbs();
      md_dist(i, j) = r.sum();
    }
  }
  return md_dist;
}

Eigen::MatrixXd md_dist_smat(const Eigen::MatrixXd &in_locs)
{
  int n = in_locs.rows();
  Eigen::MatrixXd md_dist = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - in_locs.row(j)).transpose().cwiseAbs();
      md_dist(i, j) = r.sum();
      md_dist(j, i) = md_dist(i, j);
    }
  }
  return md_dist;
}

Eigen::VectorXd md_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd &out_loc)
{
  int n_in = in_locs.rows();
  Eigen::VectorXd md_dist = Eigen::VectorXd::Zero(n_in);
  for (int i = 0; i < n_in; i++)
  {
    Eigen::VectorXd r = (in_locs.row(i) - out_loc.transpose()).transpose().cwiseAbs();
    md_dist(i) = r.sum();
  }
  return md_dist;
}

Eigen::MatrixXd cd_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs)
{
  int n_in = in_locs.rows();
  int n_out = out_locs.rows();
  Eigen::MatrixXd cd_dist = Eigen::MatrixXd::Zero(n_in, n_out);
  for (int i = 0; i < n_in; i++)
  {
    for (int j = i; j < n_out; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - out_locs.row(j)).transpose().cwiseAbs();
      cd_dist(i, j) = r.maxCoeff();
      cd_dist(j, i) = cd_dist(i, j);
    }
  }
  return cd_dist;
}

Eigen::MatrixXd cd_dist_smat(const Eigen::MatrixXd &in_locs)
{
  int n = in_locs.rows();
  Eigen::MatrixXd cd_dist = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - in_locs.row(j)).transpose().cwiseAbs();
      cd_dist(i, j) = r.maxCoeff();
      cd_dist(j, i) = cd_dist(i, j);
    }
  }
  return cd_dist;
}

Eigen::VectorXd cd_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd out_loc)
{
  int n_in = in_locs.rows();
  Eigen::VectorXd cd_dist = Eigen::VectorXd::Zero(n_in);
  for (int i = 0; i < n_in; i++)
  {
    Eigen::VectorXd r = (in_locs.row(i) - out_loc.transpose()).transpose().cwiseAbs();
    cd_dist(i) = r.maxCoeff();
  }
  return cd_dist;
}

Eigen::MatrixXd mk_dist_mat(const Eigen::MatrixXd &in_locs, const Eigen::MatrixXd &out_locs, double p)
{
  int n_in = in_locs.rows();
  int n_out = out_locs.rows();
  Eigen::MatrixXd mk_dist = Eigen::MatrixXd::Zero(n_in, n_out);
  for (int i = 0; i < n_in; i++)
  {
    for (int j = 0; j < n_out; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - out_locs.row(j)).transpose().cwiseAbs().pow(p);
      mk_dist(i, j) = std::pow(r.sum(), 1.0 / p);
    }
  }
  return mk_dist;
}

Eigen::MatrixXd mk_dist_smat(const Eigen::MatrixXd &in_locs, double p)
{
  int n = in_locs.rows();
  Eigen::MatrixXd mk_dist = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      Eigen::VectorXd r = (in_locs.row(i) - in_locs.row(j)).transpose().cwiseAbs().pow(p);
      mk_dist(i, j) = std::pow(r.sum(), 1.0 / p);
      mk_dist(j, i) = mk_dist(i, j);
    }
  }
  return mk_dist;
}

Eigen::VectorXd mk_dist_vec(const Eigen::MatrixXd &in_locs, const Eigen::VectorXd &out_loc, double p)
{
  int n_in = in_locs.rows();
  Eigen::VectorXd mk_dist = Eigen::VectorXd::Zero(n_in);
  for (int i = 0; i < n_in; i++)
  {
    Eigen::VectorXd r = (in_locs.row(i) - out_loc.transpose()).transpose().cwiseAbs().pow(p);
    mk_dist(i) = std::pow(r.sum(), 1.0 / p);
  }
  return mk_dist;
}

Eigen::VectorXd sp_dists(const Eigen::MatrixXd &dp, const Eigen::VectorXd &loc) {
  int N = dp.rows();
  Eigen::VectorXd dists = Eigen::VectorXd::Zero(N);
  double uout = loc(0), vout = loc(1);
  for (int j = 0; j < N; j++) {
    dists(j) = sp_gcdist(dp(j, 0), uout, dp(j, 1), vout);
  }
  return dists;
}

/*
double sp_gcdist(double lon1, double lon2, double lat1, double lat2) {
  double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
  double w, R, a, f, D, H1, H2;
  double lat1R, lat2R, lon1R, lon2R, DE2RA;
  
  DE2RA = M_PI/180;
  a = 6378.137;              // WGS-84 equatorial radius in km
  f = 1.0/298.257223563;     // WGS-84 ellipsoid flattening factor
  
  if (fabs(lat1 - lat2) < DBL_EPSILON) {
    if (fabs(lon1 - lon2) < DBL_EPSILON) {
      return 0.0;
      // Wouter Buytaert bug caught 100211
    } else if (std::abs((std::abs(lon1) + std::abs(lon2)) - 360.0) < DBL_EPSILON) {
      return 0.0;
    }
  }
  lat1R = lat1*DE2RA;
  lat2R = lat2*DE2RA;
  lon1R = lon1*DE2RA;
  lon2R = lon2*DE2RA;
  
  F = ( lat1R + lat2R )/2.0;
  G = ( lat1R - lat2R )/2.0;
  L = ( lon1R - lon2R )/2.0;
  
  sinG2 = std::pow( std::sin( G ), 2 );
  cosG2 = std::pow( std::cos( G ), 2 );
  sinF2 = std::pow( std::sin( F ), 2 );
  cosF2 = std::pow( std::cos( F ), 2 );
  sinL2 = std::pow( std::sin( L ), 2 );
  cosL2 = std::pow( std::cos( L ), 2 );
  
  S = sinG2*cosL2 + cosF2*sinL2;
  C = cosG2*cosL2 + sinF2*sinL2;
  
  w = std::atan( std::sqrt( S/C ) );
  R = std::sqrt( S*C )/w;
  
  D = 2*w*a;
  H1 = ( 3*R - 1 )/( 2*C );
  H2 = ( 3*R + 1 )/( 2*S );
  
  return D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 ); 
}
*/
