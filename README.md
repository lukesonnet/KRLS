# KRLS
[![Travis-CI Build Status](https://travis-ci.org/lukesonnet/KRLS.svg?branch=master)](https://travis-ci.org/lukesonnet/KRLS) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/lukesonnet/KRLS?branch=master&svg=true)](https://ci.appveyor.com/project/lukesonnet/KRLS) [![Coverage Status](https://coveralls.io/repos/github/lukesonnet/KRLS/badge.svg?branch=master)](https://coveralls.io/github/lukesonnet/KRLS?branch=master)

This package provides methods for fitting flexible functional forms for continuous and binary outcomes. This package is under development and may behave unexpectedly. It is intended to replace the KRLS package available on CRAN. Please email [Luke Sonnet](mailto:luke.sonnet@gmail.com) or leave an issue for Luke Sonnet or Chad Hazlett.

### Install latest version

You can install the latest version by running:
```R
devtools::install_github('lukesonnet/KRLS')
```

### Troubleshooting installation

This version uses `Rcpp` extensively for speed reasons. These means you need to have the right compilers on your machine.

#### Windows
If you are on Windows, you will need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/) if you haven't already. If you still are having difficulty with installing and it says that the compilation failed, try installing it without support for multiple architectures:
```R
devtools::install_github('lukesonnet/KRLS', args=c('--no-multiarch'))
```

#### Mac OSX

In order to compile the `C++` in this package, `RcppArmadillo` will require you to have compilers installed on your machine. You may already have these, but you can install them by running:

```bash
xcode-select --install
```

If you are having problems with this install on Mac OSX, specifically if you are getting errors with either `lgfortran` or `lquadmath`, then try open your Terminal and try the following:

```bash
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Also see section 2.16 [here](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf)
