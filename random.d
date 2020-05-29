module gretl.random;
import std.conv;

void randInit() {
	gretl_rand_init();
}

void randFree() {
	gretl_rand_free();
}

void setSeed(long seed) { 
  gretl_rand_set_seed(seed.to!uint);
}

uint getSeed() { 
  return gretl_rand_get_seed();
}

uint uniformInt() {
	return gretl_rand_int();
}

uint uniformInt(uint k) { 
  return gretl_rand_int_max(k);
}

double runif() { 
  return gretl_rand_01();
}

double[] runif(int n, double min=0.0, double max=1.0) {
  auto result = new double[n];
  gretl_rand_uniform_minmax(result.ptr, 0, n-1, min, max);
  return result;
}

void ruinfUnsafe(double * ptr, int len, double min=0.0, double max=1.0) {
	gretl_rand_uniform_minmax(ptr, 0, len-1, min, max);
}

DoubleVector rnorm(int n, double mean=0.0, double sd=1.0) {
  auto result = DoubleVector(n);
  gretl_rand_normal_full(result.ptr, 0, n-1, mean, sd);
  return result;
}

double rnorm() { 
  return gretl_one_snormal(); 
}

DoubleVector rmvnorm(GretlMatrix mu, GretlMatrix V) {
  enforce(mu.cols == 1, "rmvnorm: mu needs to have one column");
  enforce(mu.rows == V.rows, "rmvnorm: mu and v need to have the same number of rows");
  enforce(V.rows == V.cols, "rmvnorm: v needs to be square");
  return DoubleVector(mu + chol(V)*rnorm(mu.rows));
}

DoubleVector[] rmvnorm(int n, GretlMatrix mu, GretlMatrix V) {
	DoubleVector[] result;
	result.reserve(n);
	foreach(draw; 0..n) {
		result ~= rmvnorm(mu, V);
	}
	return result;
}

int[] indexSample(int n) {
	auto result = new int[n];
	gretl_rand_int_minmax(result.ptr, n, 0, n-1);
	return result;
}

extern(C) {
	double gretl_rand_gamma_one(double shape, double scale);
	//~ GretlMatrix * inverse_wishart_matrix(GretlMatrix * S, int v, int * err);
	int gretl_rand_F(double * a, int t1, int t2, int v1, int v2); 
}

double rf(int v1, int v2) {
	double[] result = [0.0];
	gretl_rand_F(result.ptr, 0, 1, v1, v2);
	return result[0];
}

DoubleVector rf(int n, int v1, int v2) {
	auto result = DoubleVector(n);
	gretl_rand_F(result.ptr, 0, n, v1, v2);
	return result;
}

double rgamma(double shape, double scale) {
	return gretl_rand_gamma_one(shape, scale);
}

DoubleVector rgamma(int k, double shape, double scale) {
	auto result = DoubleVector(k);
	foreach(ii; 0..k) {
		result[ii] = rgamma(shape, scale);
	}
	return result;
}

/* This avoids allocation completely.
 * Main use is when you want to reuse an array.
 * Can be used with a GretlMatrix (completely avoiding the GC),
 * RMatrix (allocated by R), or DoubleMatrix (garbage collected).
 */
void rgamma(GretlMatrix gm, double shape, double scale) {
	foreach(ii; 0..gm.rows) {
		gm[ii,0] = gretl_rand_gamma_one(shape, scale);
	}
}

/*
 * S: Scale matrix
 * v: Degrees of freedom
 */
GretlMatrix * invWishartUnsafe(GretlMatrix S, int v) {
	GretlMatrix * gm;
	int * err;
	return inverse_wishart_matrix(&S, v, err);
}

DoubleMatrix invWishart(GretlMatrix S, int v) {
	return DoubleMatrix(invWishartUnsafe(S, v));
}

extern (C) {
  void gretl_rand_init();
  void gretl_rand_free();
  void gretl_rand_set_seed(uint seed);
  void gretl_rand_set_multi_seed(const GretlMatrix * seed);
  uint gretl_rand_int();
  uint gretl_rand_int_max(uint m);
	int gretl_rand_int_minmax(int *a, int n, int min, int max);
  double gretl_rand_01();
  double gretl_one_snormal();
  void gretl_rand_uniform(double * a, int t1, int t2);
  int gretl_rand_uniform_minmax(double * a, int t1, int t2, double min, double max);
  void gretl_rand_normal(double * a, int t1, int t2);
  int gretl_rand_normal_full(double * a, int t1, int t2, double mean, double sd);
  int gretl_rand_chisq(double * a, int t1, int t2, double v);
  int gretl_rand_student(double * a, int t1, int t2, double v);
  int gretl_rand_F(double * a, int t1, int t2, int v1, int v2);
  int gretl_rand_binomial(double * a, int t1, int t2, int n, double p);
  int gretl_rand_poisson(double * a, int t1, int t2, const double * m, int vec);
  int gretl_rand_weibull(double * a, int t1, int t2, double shape, double scale);
  int gretl_rand_gamma(double * a, int t1, int t2, double shape, double scale);
  double gretl_rand_gamma_one(double shape, double scale);
  int gretl_rand_GED(double * a, int t1, int t2, double nu);
  int gretl_rand_beta(double * x, int t1, int t2, double s1, double s2);
  int gretl_rand_beta_binomial(double * x, int t1, int t2, int n, double s1, double s2);
  uint gretl_rand_get_seed();
	//~ Omitting these functions since they're not in the 2017 release of Gretl or they rely on the matrix struct
  //~ GretlMatrix * halton_matrix(int m, int r, int offset, int * err);
  //~ GretlMatrix * inverse_wishart_matrix(const GretlMatrix * S, int v, int * err);
  //~ GretlMatrix * inverse_wishart_sequence(const GretlMatrix * S, int v, int replics, int * err);
  //~ void gretl_rand_set_box_muller(int s);
  //~ int gretl_rand_get_box_muller();
}
