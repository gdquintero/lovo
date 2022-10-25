export ALGENCAN=/home/gustavo/algencan-3.1.1

rm -f lovo3

gfortran -O3 -w -fcheck=all -g lovo3.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o lovo3

./lovo3
