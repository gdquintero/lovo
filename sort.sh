#! /bin/bash
#
gfortran -c -Wall sort.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
echo "Normal end of execution."
