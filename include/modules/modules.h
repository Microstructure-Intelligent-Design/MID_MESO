#pragma once
#include "Information.h"
#include "InterfaceEnergy.h"
#include "Kinetics.h"
#include "Mechanics.h"
#include "ElectricField.h"
#include "Nucleation.h"
#include "OptimizationAlgorithm.h"
#include "Thermodynamics.h"
#include "MagneticField.h"
#include "CustomSolver.h"
#include "FluidField.h"
#include "license.h"

#ifdef _WIN32
//windowsƽ̨ x86 or x68
#ifdef _WIN64
 //x64
#ifdef _DEBUG
#pragma comment(lib,"lib/x64/debug/MID.lib")
#else
#pragma comment(lib,"lib/x64/release/MID.lib")
#endif
#else
 //x86
#ifdef _DEBUG
#pragma comment(lib,"lib/x86/debug/MID.lib")
#else
#pragma comment(lib,"lib/x86/release/MID.lib")
#endif
#endif //_WIN64
#endif //_WIN32