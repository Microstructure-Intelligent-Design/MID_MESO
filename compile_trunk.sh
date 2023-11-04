#! /bin/bash
rm -f *.o
echo "> The trunk is compiling... Please wait..."
g++ -std=c++11 -fopenmp -static -c -w src/modules/license.cpp -o lib/linux/license.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/CustomSolver.cpp -o lib/linux/CustomSolver.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/ElectricField.cpp -o lib/linux/ElectricField.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/InterfaceEnergy.cpp -o lib/linux/InterfaceEnergy.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/Kinetics.cpp -o lib/linux/Kinetics.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/MagneticField.cpp -o lib/linux/MagneticField.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/Mechanics.cpp -o lib/linux/Mechanics.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/Nucleation.cpp -o lib/linux/Nucleation.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/FluidField.cpp -o lib/linux/FluidField.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/OptimizationAlgorithm.cpp -o lib/linux/OptimizationAlgorithm.o
g++ -std=c++11 -fopenmp -static -c -w src/modules/Thermodynamics.cpp -o lib/linux/Thermodynamics.o
g++ -std=c++11 -fopenmp -static -c -w src/MPF.cpp -o lib/linux/MPF.o

echo "> The trunk has been builded"

