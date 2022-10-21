#! /bin/bash
#
gfortran -c -Wall subset.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
echo "Normal end of execution."
