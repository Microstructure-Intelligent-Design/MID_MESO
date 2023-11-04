#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"

#define linear_n value2
#define Max_siteFraction_iterate_times 1000

namespace pf {


	class Thermodynamics {
	public:

		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information);
		void init_special_thermodynamicModel_in_phase(pf::PhaseEntry& phase);
		void pretreatment_thermodynamic_data_of_onePhase(pf::PhaseNode& node, int phaseIndex);
		void aftertreatment_thermodynamic_data_of_onePhase(pf::PhaseNode& node, int phaseIndex);
		void cal_thermodynamics_increment_of_oneNode(pf::PhaseNode& node);
		double cal_beta_nucleates_in_alpha(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double nucleus_radius);
		void cal_thermodynamics_data_of_onePhase_in_fixedBoundary(int phaseIndex);
		XNode get_same_ChemicalPotential_composition(pf::PhaseNode& node, pf::PhaseEntry& sample_phase, pf::PhaseEntry& aim_phase, int nucleus_radius);
		void clear();

		double cal_xLi_by_chempotential(int phaseproperty, double chempotential, double temperature, double energy_precision) {
			PhaseNode node;
			node.add_phase(0, phaseproperty, 0, 1.0);
			node[0].x.add_con(materials::Li, information->materialSystem.phases[phaseproperty].x[materials::Li].value);
			node[0].x.add_con(materials::Sn, information->materialSystem.phases[phaseproperty].x[materials::Sn].value);
			node[0].potential.add_con(materials::Li);
			node[0].potential.add_con(materials::Sn);
			node.tempValues.temperature = temperature;
			double delt_x = 0.02, precision = 0.0;
			bool is_foward_back;
			{
				node[0].x[materials::Li].value = X_EPSILON;
				node[0].x[materials::Sn].value = 1.0 - X_EPSILON;
				node[0].potential[materials::Li].chemical_part = 0.0;
				pretreatment_thermodynamic_data_of_onePhase(node, 0);
				information->materialSystem.functions.Energy(node, node[0], information->dynamicCollection);
				information->materialSystem.functions.Potential(node, node[0], information->dynamicCollection);
				aftertreatment_thermodynamic_data_of_onePhase(node, 0);
				if (node[0].potential[materials::Li].chemical_part > chempotential)
					return X_EPSILON;

				node[0].x[materials::Li].value = 1.0 - X_EPSILON;
				node[0].x[materials::Sn].value = X_EPSILON;
				node[0].potential[materials::Li].chemical_part = 0.0;
				pretreatment_thermodynamic_data_of_onePhase(node, 0);
				information->materialSystem.functions.Energy(node, node[0], information->dynamicCollection);
				information->materialSystem.functions.Potential(node, node[0], information->dynamicCollection);
				aftertreatment_thermodynamic_data_of_onePhase(node, 0);
				if (node[0].potential[materials::Li].chemical_part < chempotential)
					return 1.0 - X_EPSILON;
			}
			if (node[0].potential[materials::Li].chemical_part > chempotential)
				is_foward_back = true;
			else
				is_foward_back = false;

			do
			{
				if ((node[0].x[materials::Li].value - delt_x > X_EPSILON && is_foward_back) || ((node[0].x[materials::Li].value + delt_x) < (1.0 - X_EPSILON) && (!is_foward_back))) {
					if (is_foward_back) {
						node[0].x[materials::Li].value -= delt_x;
						node[0].x[materials::Sn].value += delt_x;
					}
					else {
						node[0].x[materials::Li].value += delt_x;
						node[0].x[materials::Sn].value -= delt_x;
					}
					node[0].potential[materials::Li].chemical_part = 0.0;
					pretreatment_thermodynamic_data_of_onePhase(node, 0);
					information->materialSystem.functions.Energy(node, node[0], information->dynamicCollection);
					information->materialSystem.functions.Potential(node, node[0], information->dynamicCollection);
					aftertreatment_thermodynamic_data_of_onePhase(node, 0);
					if (node[0].potential[materials::Li].chemical_part < chempotential && is_foward_back == true) {
						delt_x /= 1.5;
						is_foward_back = false;
					}
					else if (node[0].potential[materials::Li].chemical_part > chempotential && is_foward_back == false) {
						delt_x /= 1.5;
						is_foward_back = true;
					}
					precision = abs(node[0].potential[materials::Li].chemical_part - chempotential);
				}
				else
					delt_x /= 1.5;

			} while (precision > energy_precision);
			return node[0].x[materials::Li].value;
		}

	private:
		void do_energy_treatment(pf::PhaseEntry& phase, double molarVolume);
		FieldStorage_forPhaseNode* simulationField;
		Information* information;
		int cal_constituent_of_onePhase_by_energy_minimization(pf::PhaseNode& node, int phaseIndex);
	};
}