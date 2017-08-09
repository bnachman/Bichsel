gfortran -c covlib.for
gfortran -c cov3.for
gfortran -c timdat.for
gfortran -o cov.exe cov.for covlib.o cov3.o timdat.o