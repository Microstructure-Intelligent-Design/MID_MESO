#! /bin/bash
rm -f *.o
echo "> The program is compiling... Please wait..."
g++ -std=c++11 -fopenmp -c -w runfile.cpp -o runfile.o
g++ -std=c++11 -fopenmp -o run_file runfile.o ../../lib/linux/*.o -L../../lib/linux/lib -lfftw3
rm -f *.o

echo "> The runfile has been builded"

