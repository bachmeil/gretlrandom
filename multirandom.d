module gretl.multirandom;
import std.conv, std.range, std.stdio;
import gretl.matrix, gretl.vector;

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

extern(C) {
	GretlMatrix * inverse_wishart_matrix(GretlMatrix * S, int v, int * err);
}

DoubleVector rf(int n, int v1, int v2) {
	auto result = DoubleVector(n);
	gretl_rand_F(result.ptr, 0, n, v1, v2);
	return result;
}

DoubleVector rgamma(int k, double shape, double scale) {
	auto result = DoubleVector(k);
	foreach(ii; 0..k) {
		result[ii] = rgamma(shape, scale);
	}
	return result;
}

//~ /* This avoids allocation completely.
 //~ * Main use is when you want to reuse an array.
 //~ * Can be used with a GretlMatrix (completely avoiding the GC),
 //~ * RMatrix (allocated by R), or DoubleMatrix (garbage collected).
 //~ */
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

/* pdf of multivariate normal with mean mu, covariance V
 * x is the vector
 * Tested and reproduces the results in R. */
double dmvnorm(DoubleVector x, DoubleVector mu, DoubleMatrix V, bool ln=false) {
	DoubleVector dev = x - mu;
	DoubleMatrix inv_v = inv(V);
	DoubleVector z = DoubleVector(inv_v * dev);
	double quadform = 0.0;
	foreach(val; z) {
		quadform += val^^2;
	}
	double logSqrtDetSigma = 0.0;
	foreach(ii; 0..V.rows) {
		logSqrtDetSigma += log(V[ii,ii]);
	}
	if (ln) {
		return -0.5*quadform - logSqrtDetSigma - 0.5*V.rows*log(2.0*PI);
	} else {
		return exp(-0.5*quadform - logSqrtDetSigma - 0.5*V.rows*log(2.0*PI));
	}
}

extern (C) {
  void gretl_rand_set_multi_seed(const GretlMatrix * seed);
  GretlMatrix * halton_matrix(int m, int r, int offset, int * err);
  GretlMatrix * inverse_wishart_matrix(const GretlMatrix * S, int v, int * err);
  GretlMatrix * inverse_wishart_sequence(const GretlMatrix * S, int v, int replics, int * err);
  void gretl_rand_set_box_muller(int s);
  int gretl_rand_get_box_muller();
}
