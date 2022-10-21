export ALGENCAN=/home/gustavo/algencan-3.1.1

rm -f lovo

gfortran -O3 -w -fcheck=all -g lovo.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o lovo

./lovo
