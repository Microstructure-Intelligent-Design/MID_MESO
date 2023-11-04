#pragma once
#include "../baseTools/baseTools.h"
using namespace std;
namespace pf {

	namespace poissonEquationSolver {
		static void rhs_cal(pf::PhaseNode& node, int rhs_index) {
			node.customValues[rhs_index] = 0.0;
		}
		static void boundary(pf::PhaseNode& node, int lhs_index) {
			return;
		}
	}

	class PoissonEquationSolver
	{
	public:
		PoissonEquationSolver(FieldStorage_forPhaseNode& _phaseMesh, int LHS_index, int RHS_index, string _solver_name = "PoissonEquationSolver") {
			phaseMesh = &_phaseMesh;
			solver_name = _solver_name;
			LHS_INDEX = LHS_index;
			RHS_INDEX = RHS_index;
			setRHS_value = poissonEquationSolver::rhs_cal;
			boundary = poissonEquationSolver::boundary;
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				node->customValues.add_double(LHS_INDEX, 0.0);
				node->customValues.add_double(RHS_INDEX, 0.0);
			}
		}
		~PoissonEquationSolver() {
			clear();
		};
		void init_field(double LHS_value, double RHS_value) {
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				node->customValues.add_double(LHS_INDEX, LHS_value);
				node->customValues.add_double(RHS_INDEX, RHS_value);
			}
		}
		int solve_whole_domain(double accuracy, int MAXIterations, bool debug_solver = false, int output_step = 1000);
		int solve_whole_domain_with_average_boundary_condition(double accuracy, int MAXIterations, double average_value = 0.0, bool debug_solver = false, int output_step = 1000);
		int solve_phases_region(double accuracy, int MAXIterations, vector<int> phaseIndexes, bool debug_solver = false, int output_step = 1000);
		void set_RHS_calfunc(void(*RHS_cal)(pf::PhaseNode&, int)) {
			setRHS_value = RHS_cal;
		}
		void set_BoundaryCondition_calfunc(void(*Boundary)(pf::PhaseNode&, int)) {
			boundary = Boundary;
		}
		string solver_name;
	private:
		int LHS_INDEX, RHS_INDEX;
		FieldStorage_forPhaseNode* phaseMesh;
		void(*setRHS_value)(pf::PhaseNode&, int);
		void(*boundary)(pf::PhaseNode&, int);
		void set_RHS_value();
		void clear() {
			phaseMesh = nullptr;
			setRHS_value = nullptr;
			boundary = nullptr;
		}
	};

	namespace allenCahnEquationSolver {
		// grain_index from 0 to grain_number - 1
		static double dF_dphi_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
			return 0.0;
		}
		static double L_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
			return 0.0;
		}
		static double Source_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
			return 0.0;
		}
	}

	// partial(phi_i) / partial(t) = - Mobility_i * variation(F) / variation(phi_i) + Source
	class AllenCahnSolver {
	public:
		AllenCahnSolver(FieldStorage_forPhaseNode& _phaseMesh, int _grains_number, int _grains_start_index = SOLVER_ALLEN_CAHN, string _solver_name = "AllenCahnSolver") {
			solver_name = _solver_name;
			phaseMesh = &_phaseMesh;
			grains_number = _grains_number;
			grains_start_index = _grains_start_index;
			dF_dphi = allenCahnEquationSolver::dF_dphi_cal;
			L = allenCahnEquationSolver::L_cal;
			Source = allenCahnEquationSolver::Source_cal;
		}
		~AllenCahnSolver() {
			clear();
		};
		void set_dF_dphi_func(double(*dF_dphi_cal)(pf::PhaseNode&, int, int)) {
			dF_dphi = dF_dphi_cal;
		}
		void set_L_func(double(*L_cal)(pf::PhaseNode&, int, int)) {
			L = L_cal;
		}
		void set_Source_func(double(*Source_cal)(pf::PhaseNode&, int, int)) {
			Source = Source_cal;
		}
		// retuan MAX_PHI_VARIATION;
		double solve_one_step(int istep, double dt, bool adjust_phi_0_1 = false);
		void init_field() {
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				for (int grain = 0; grain < grains_number; grain++) {
					node->customValues.add_double((grain + grains_start_index), 0.0);
					node->customValues.add_double(-(grain + grains_start_index), 0.0);
				}
			}
		}
		string solver_name;
		int grains_number;
		int grains_start_index;
	private:
		FieldStorage_forPhaseNode* phaseMesh;
		double(*dF_dphi)(pf::PhaseNode&, int, int);
		double(*L)(pf::PhaseNode&, int, int);
		double(*Source)(pf::PhaseNode&, int, int);
		void clear() {
			phaseMesh = nullptr;
			dF_dphi = nullptr;
			L = nullptr;
			Source = nullptr;
		}
	};

	namespace cahnHilliardEquationSolver {
		// grain_index from 0 to grain_number - 1
		static double dF_dc_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return 0.0;
		}
		static double Mobility_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return 0.0;
		}
		static double Source_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
			return 0.0;
		}
	}

	// partial(c_i) / partial(t) = delt(M_i * delt(variation(F) / variation(c_i))) + Source
	class CahnHilliardSolver {
	public:
		CahnHilliardSolver(FieldStorage_forPhaseNode& _phaseMesh, string _solver_name = "CahnHilliardSolver") {
			solver_name = _solver_name;
			phaseMesh = &_phaseMesh;
			dF_dc = cahnHilliardEquationSolver::dF_dc_cal;
			Mobility = cahnHilliardEquationSolver::Mobility_cal;
			Source = cahnHilliardEquationSolver::Source_cal;
		}
		~CahnHilliardSolver() {
			clear();
		};
		void set_dF_dc_func(double(*dF_dc_cal)(pf::PhaseNode&, int, int)) {
			dF_dc = dF_dc_cal;
		}
		void set_Mobility_func(double(*Mobility_cal)(pf::PhaseNode&, int, int)) {
			Mobility = Mobility_cal;
		}
		void set_Source_func(double(*Source_cal)(pf::PhaseNode&, int, int)) {
			Source = Source_cal;
		}
		// return MAX_COMP_VARIATION
		double solve_one_step(int istep, double dt, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT, bool adjust_phi_0_1 = false);
		void init_field(int _components_number, int _comps_start_index = SOLVER_CAHN_HILLIARD) {
			components_number = _components_number;
			comps_start_index = _comps_start_index;
			for (auto node = phaseMesh->_mesh.begin(); node < phaseMesh->_mesh.end(); node++) {
				for (int comp = 0; comp < _components_number; comp++) {
					node->customValues.add_double(comp + comps_start_index, 0.0);
					node->customValues.add_double(-(comp + comps_start_index), 0.0);
					node->customValues.add_double(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD, 0.0);
					node->customValues.add_double(-(comp + BUFF_FOR_SOLVER_CAHN_HILLIARD), 0.0);
				}
			}
		}
		string solver_name;
		int components_number;
		int comps_start_index;
	private:
		FieldStorage_forPhaseNode* phaseMesh;
		double(*dF_dc)(pf::PhaseNode&, int, int);
		double(*Mobility)(pf::PhaseNode&, int, int);
		double(*Source)(pf::PhaseNode&, int, int);
		void clear() {
			phaseMesh = nullptr;
			dF_dc = nullptr;
			Mobility = nullptr;
			Source = nullptr;
		}
	};
}