# FortranCode
Fortran codes for strings, random number generation from various distributions, the negative binomial distribution, and approximations for sine and cosine. Create executables with

```
gfortran -o xstring_funcs.exe string_funcs.f90 xstring_funcs.f90
gfortran -o xbinomial.exe kind.f90 binomial.f90 xbinomial.f90
gfortran -o xsim_negative_binomial.exe kind.f90 binomial.f90 xsim_negative_binomial.f90
gfortran -o xtrig.exe kind.f90 trig.f90 xtrig.f90
gfortran -o xxrandom.exe kind.f90 random.f90 xxrandom.f90
```
The program `min_max.f90` from a [discussion at Fortran Discourse](https://fortran-lang.discourse.group/t/speed-of-array-intrinsics/3251) compares the speed of various ways of computing the min, max, and mean of a 1D array. Compiling on Windows with `gfortran -Ofast -march=native -flto min_max.f90` and running, sample output is

```
    time                               task
  0.0938               minval(x), maxval(x)
  0.0469                         min_max(x)
  0.1406       minval(x),maxval(x),sum(x)/n
  0.0469                    min_max_mean(x)
```
so it was faster to compute min and max in a single loop with your own function than to call the minval and maxval intrinsics.
