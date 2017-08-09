This package contains the silicon ionizaion simulation using the Bichsel code.
The details of the Bichsel ionization calculations for silicon are described
in the paper: For the original paper (Rev. Mod. Phys. 60, p663, 1988):
  http://prola.aps.org/abstract/RMP/v60/i3/p663_1
For all applications using the content of this program, please reference to
the RMP paper above and COV/FOLD program with private communication to Hans
Bichsel. The are other related information on Hans' web page:
  http://faculty.washington.edu/hbichsel/
you may be interested in. 

Maintainer: Su Dong (last update Apr/13/2016)

Translation to C++: Maurice Garcia-Sciveres

=======================================================

Instructions:

For the Fortran code, you can do

cd bichsel/FortranVersion
sh link_cov.sh
./cov.exe

this will produce COV.OPE and COV.SPE

you can see the names of the variables in Bichsel_COV_variables.pdf

For the C++ code, you can do

cd bichsel/C++
root -b
.L bichsel-lib.C++
mymain()