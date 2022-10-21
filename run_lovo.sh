export ALGENCAN=/home/gustavo/algencan-3.1.1

rm -f ovo

gfortran -O3 -w -fcheck=all -g ovo.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o ovo

./ovo