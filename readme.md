# gretlrandom

A D wrapper over Gretl's random number generators. Currently (05/29/2020) avoids any dependencies on other libraries so it's independent of those design decisions. Those functions will be provided in a different module.

# Startup and shutdown

- Must call `randInit` before calling any of these functions.
- Should call `randFree` to clean up when you're done.

```
void randInit()
void randFree()
```

- Set and get the seed with `setSeed` and `getSeed`. Should normally do this to ensure replicability. `setSeed` requires a `uint` on the C side. Since this is a rarely called function, `to!uint` is used to convert the value for convenience. Note that this operation is safe. `to` confirms that the values are within the acceptable range and throws an exception if not.

```
void setSeed(long seed)
uint getSeed()
```

# Generating integers

- `uniformInt` returns one integer over the full range of acceptable `int` values (`0` through `2^31 - 1`) if no argument is provided.
- `uniformInt` returns an integer less than the argument when called with an argument.

```
uint uniformInt(uint k) { 
  return gretl_rand_int_max(k);
}
```

# Generating uniform double values

- `runif` returns one double between 0 and 1 if called without an argument.
- `runif` returns a double array of `n` elements if called with an argument.
- `runifUnsafe` should normally not be called. It is used by other library functions generating arrays of uniform values.

```
double runif()
double[] runif(int n, double min=0.0, double max=1.0)
void ruinfUnsafe(double * ptr, int len, double min=0.0, double max=1.0)
```

# Generating normal double values

