export ALGENCAN=/home/gustavo/algencan-4.0.0

rm -f lovo

gfortran -c -O3 -w -fcheck=all -Wall -I$ALGENCAN/sources/algencan/inc lovo.f90 

gfortran -g -w -fcheck=all lovo.o -L$ALGENCAN/sources/algencan/lib -lalgencan -L$ALGENCAN/sources/hsl/lib -lhsl -L$ALGENCAN/sources/blas/lib -lblas sort.o subset.o -o lovo

./lovo

