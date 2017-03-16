# KRLS

This package is in alpha and may behave unexpectedly. It provides additional features and better performance relative to the KRLS package available on CRAN. Please email [Luke Sonnet](mailto:luke.sonnet@gmail.com) or leave an issue for Luke Sonnet or Chad Hazlett.

If you are having problems with this install on Mac OSX, specifically if you are getting errors with either `lgfortran` or `lquadmath`, then try open your Terminal and try the following:

    curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
    sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /

Also see section 2.16 [here](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf)
