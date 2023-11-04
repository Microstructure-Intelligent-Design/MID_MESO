#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"

namespace pf {
	class Kinetics
	{
	public:
		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information);

		void cal_chemicalDiffusivity_of_phase(pf::PhaseNode& node, pf::PhaseEntry& phase);
		void cal_chemicalDiffusivity_of_node(pf::PhaseNode& node);
		void cal_chemicalDiffusivity_of_phase_in_fixedBoundary(int phaseIndex);
		void cal_chemicalMobility_of_phase(pf::PhaseNode& node, pf::PhaseEntry& phase);
		void cal_chemicalMobility_of_node(pf::PhaseNode& node);
		void cal_chemicalMobility_of_phase_in_fixedBoundary(int phaseIndex);

		void cal_kinetic_diffusion_term(pf::PhaseNode& node);
		void cal_kinetic_chemical_reaction_term(pf::PhaseNode& currentNode);
		void cal_kinetic_phase_transition_term(pf::PhaseNode& currentNode);

		//debug
		void write_scalar_PhaseTransitionFlux_of_phaseCon(ofstream& fout, int phaseIndex, int compIndex);
		void write_scalar_ChemicalReactionFlux_of_phaseCon(ofstream& fout, int phaseIndex, int compIndex);
		void write_scalar_DiffusionFlux_of_phaseCon(ofstream& fout, int phaseIndex, int compIndex);
		void clear();

	private:
		void cal_laplace_Grad_phaseCon_kineticCoef(pf::PhaseNode& node, int phaseIndex);
		void cal_laplace_Grad_phasePotential_kineticCoef(pf::PhaseNode& node, int phaseIndex);
		void cal_laplace_Grad_con_kineticCoef(pf::PhaseNode& node);
		void cal_laplace_Grad_potential_kineticCoef(pf::PhaseNode& node);
		FieldStorage_forPhaseNode* simulationField;
		Information* information;

	};

}
