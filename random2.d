module gretl.random;
import std.conv, std.range, std.stdio;

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

extern (C) {
  //~ void gretl_rand_set_multi_seed(const GretlMatrix * seed);
  //~ GretlMatrix * halton_matrix(int m, int r, int offset, int * err);
  //~ GretlMatrix * inverse_wishart_matrix(const GretlMatrix * S, int v, int * err);
  //~ GretlMatrix * inverse_wishart_sequence(const GretlMatrix * S, int v, int replics, int * err);
  //~ void gretl_rand_set_box_muller(int s);
  //~ int gretl_rand_get_box_muller();
}
