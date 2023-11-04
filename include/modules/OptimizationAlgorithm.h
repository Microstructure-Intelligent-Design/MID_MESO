#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"

namespace pf{
	enum GeometricDirection {None, Surf_XY, Surf_XZ, Surf_YZ, Axis_X, Axis_Y, Axis_Z, Original_Point};
	class ConField_optimizedIndex {
	public:
		ConField_optimizedIndex(int x_increment = 0, int y_increment = 0, int z_increment = 0) {
			X_increment = x_increment;
			Y_increment = y_increment;
			Z_increment = z_increment;
		};
		ConField_optimizedIndex& operator=(const ConField_optimizedIndex& n) {
			X_increment = n.X_increment;
			Y_increment = n.Y_increment;
			Z_increment = n.Z_increment;
			return *this;
		};
		int X_increment;
		int Y_increment;
		int Z_increment;
	};
	class ConField_optimizedIndexBox {
	public:
		// angle boundary
		// alpha : angle of z-axis and the normal direction ; alpha = pow((grad_x^2 + grad_y^2) / grad_z^2, 0.5)
		// beta : angle in surf_XY, based on x-axis ; beta = abs(grad_y / grad_x)
		double alpha_down;
		double alpha_up;
		double beta_down;
		double beta_up;
		GeometricDirection ave_jump;
		//
		std::vector<ConField_optimizedIndex> indexBox;   // one direction, 6 points, the reverse direction is handled in total class
		ConField_optimizedIndexBox& operator=(const ConField_optimizedIndexBox& n) {
			alpha_down = n.alpha_down;
			alpha_up = n.alpha_up;
			beta_down = n.beta_down;
			beta_up = n.beta_up;
			ave_jump = n.ave_jump;
			indexBox = n.indexBox;
			return *this;
		};
		ConField_optimizedIndexBox() {
			alpha_down = 0;
			alpha_up = 0;
			beta_down = 0;
			beta_up = 0;
			ave_jump = GeometricDirection::None;
		};
		void do_symmetry(pf::GeometricDirection symmetry_base);
	};
	class ConField_optimizedBoxPool {
	public:
		ConField_optimizedBoxPool() {};
		~ConField_optimizedBoxPool() {};
		void init(Dimension dimension);
		ConField_optimizedIndexBox find_box_by_gradient(double grad_x, double grad_y, double grad_z);
		//ConField_optimizedIndexBox find_box_by_degree(double alpha_degree, double beta_degree);
		void clear() {
			boxPool.clear();
		};
		std::vector<ConField_optimizedIndexBox> boxPool;
		Dimension _dimension;
		ConField_optimizedBoxPool& operator=(const ConField_optimizedBoxPool& n) {
			boxPool = n.boxPool;
			return *this;
		};
	};

	struct inff {
		int phaseIndex_1;
		int phaseIndex_2;
	};
	class OptimizationAlgorithm
	{
	public:
		OptimizationAlgorithm() {};
		~OptimizationAlgorithm() {};
		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information);
		void clear();
		void optimize_phaseCon_assignment_with_cFlag(int istep, bool adjust_phi_0_1 = true);
		void optimize_bulk_increment_on_interface(int istep, int times = 2);
		// optimize : phase transition term of concentration
		void optimize_phaseCon_phase_transition_term_on_interface(int istep);

		FieldStorage_forPhaseNode* simulationField;
		Information* information;
		ConField_optimizedBoxPool ave_box_pool;
	private:
		vector<PhaseEntry*> search_for_phases_in_block_on_interface(int i, int j, int k, int phaseIndex);
		vector<double_box> auxiliary_mesh;
	};
}