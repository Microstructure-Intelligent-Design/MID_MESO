#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"
namespace pf {

	class ElectricField
	{
	public:
		ElectricField() {};
		~ElectricField() {};

		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information);
		void clear();
		int solve(double PotentialAccuracy, int MAXIterations);
		//void cal_electric_increment_of_oneNode(PhaseNode& node);
		void setElectricChargeDensity();

		FieldStorage_forPhaseNode* simulationField;
		Information* information;
		double* applied_electric_field;

	};

}