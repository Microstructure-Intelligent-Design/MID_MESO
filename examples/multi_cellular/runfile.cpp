#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\multi_cellular\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
static int grains_number = 20, soft_grain = 5, radius = 12;
pf::Vector3 velosity(0.2, 0.2, 0.0);
static vector<double> sum_phi2;
static vector<pf::Vector3> sum_interface;
static double L(pf::PhaseNode& node, int grain, int grains_start_index) {
	double L = -1.0;
	return L;
}
static double dF_dphi(pf::PhaseNode& node, int grain, int grains_start_index) {
	double dF_dphi = 0.0, gamma_n = 0.0, lambda = 7.0, kappa = 60.0, miu = 40.0, ksi = 1.5e3;
	pf::Vector3 va;
	if (grain < soft_grain) {
		gamma_n = 5.0;
		va = velosity * -1;
	}
	else {
		gamma_n = 2.5;
		va = velosity;
	}
	// first term
	dF_dphi += gamma_n * node.cal_customValues_laplace(grain + pf::SOLVER_ALLEN_CAHN, 1.0, pf::FIVE_POINT);
	// second term
	double phi = node.customValues[grain + pf::SOLVER_ALLEN_CAHN], sum = 0.0;
	for (int grain2 = 0; grain2 < grains_number; grain2++)
		if (grain2 != grain) sum += node.customValues[grain2 + pf::SOLVER_ALLEN_CAHN] * node.customValues[grain2 + pf::SOLVER_ALLEN_CAHN];
	dF_dphi += -30.0 / lambda / lambda * (gamma_n * phi * (1.0 - phi) * (1.0 - 2.0 * phi) + 2.0 * kappa * phi * sum);
	// third term
	dF_dphi += -2.0 * miu / PI / radius / radius * phi * (sum_phi2[grain] - PI * radius * radius);
	// forth term
	pf::Vector3 delt_phi((node.get_neighbor_node(Direction::x_down).customValues[grain + pf::SOLVER_ALLEN_CAHN]
		- node.get_neighbor_node(Direction::x_up).customValues[grain + pf::SOLVER_ALLEN_CAHN]) / 2.0,
		(node.get_neighbor_node(Direction::y_down).customValues[grain + pf::SOLVER_ALLEN_CAHN]
			- node.get_neighbor_node(Direction::y_up).customValues[grain + pf::SOLVER_ALLEN_CAHN]) / 2.0,
		0.0/*(node.get_neighbor_node(Direction::z_down).customValues[grain + pf::SOLVER_ALLEN_CAHN]
			- node.get_neighbor_node(Direction::z_up).customValues[grain + pf::SOLVER_ALLEN_CAHN]) / 2.0*/);
	pf::Vector3 vn = va + sum_interface[grain] * 60 * kappa / lambda / lambda / ksi;
	dF_dphi += vn * delt_phi * -1;
	return dF_dphi;
}
static double Source(pf::PhaseNode& node, int grain, int grains_start_index) {
	double Source = 0.0;
	return Source;
}
static void static_phi2(pf::FieldStorage_forPhaseNode& phaseMesh) {
	for (int grain = 0; grain < grains_number; grain++) {
		sum_phi2[grain] = 0.0;
		for (auto node = phaseMesh._mesh.begin(); node < phaseMesh._mesh.end(); node++) {
			sum_phi2[grain] += node->customValues[grain + pf::SOLVER_ALLEN_CAHN] * node->customValues[grain + pf::SOLVER_ALLEN_CAHN];
		}
	}
}
static void static_interface(pf::FieldStorage_forPhaseNode& phaseMesh) {
	for (int grain = 0; grain < grains_number; grain++) {
		sum_interface[grain].set_to_zero();
		for (auto node = phaseMesh._mesh.begin(); node < phaseMesh._mesh.end(); node++) {
			pf::Vector3 delt_phi((node->get_neighbor_node(Direction::x_down).customValues[grain + pf::SOLVER_ALLEN_CAHN] 
				- node->get_neighbor_node(Direction::x_up).customValues[grain + pf::SOLVER_ALLEN_CAHN]) / 2.0,
				(node->get_neighbor_node(Direction::y_down).customValues[grain + pf::SOLVER_ALLEN_CAHN]
					- node->get_neighbor_node(Direction::y_up).customValues[grain + pf::SOLVER_ALLEN_CAHN]) / 2.0,
				0.0/*(node->get_neighbor_node(Direction::z_down).customValues[grain + pf::SOLVER_ALLEN_CAHN]
					- node->get_neighbor_node(Direction::z_up).customValues[grain + pf::SOLVER_ALLEN_CAHN]) / 2.0*/);
			double sum = 0.0;
			for (int grain2 = 0; grain2 < grains_number; grain2++)
				if (grain2 != grain) {
					sum += node->customValues[grain2 + pf::SOLVER_ALLEN_CAHN] * node->customValues[grain2 + pf::SOLVER_ALLEN_CAHN];
				}
			sum_interface[grain] += delt_phi * node->customValues[grain + pf::SOLVER_ALLEN_CAHN] * sum;
		}
	}
}
int main(int argc, char* argv[]) {
	///< main program
	{
		///< init simulation by settings
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;
		// define a poisson solver
		pf::AllenCahnSolver allenCahn(simulation.phaseMesh, grains_number, pf::SOLVER_ALLEN_CAHN);
		allenCahn.init_field();
		allenCahn.set_dF_dphi_func(dF_dphi);
		allenCahn.set_L_func(L);
		allenCahn.set_Source_func(Source);
		sum_phi2.resize(grains_number);
		sum_interface.resize(grains_number);
		simulation.information.settings.file_settings.isCustomValueOutput = true;
		simulation.information.settings.file_settings.customValue_output.add_string(0, "cells");
		int output_step = 50;
		simulation.information.settings.file_settings.file_output_step = output_step;
		// define a solver
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			static_phi2(simulation.phaseMesh);

			static_interface(simulation.phaseMesh);
			// solver cal
			double max_variation = allenCahn.solve_one_step(istep, simulation.information.settings.disperse_settings.dt);
			// data output
			if (istep % output_step == 0) {
				simulation.add_customValue_to_allnodes(0, 0.0);
				cout << "# allenCahn solver work: MAX_VARIATION of step " << to_string(istep) << " is " << to_string(max_variation) << endl;
				for(auto node = simulation.phaseMesh._mesh.begin(); node < simulation.phaseMesh._mesh.end(); node++)
					for(int grain = 0; grain < grains_number; grain++)
						if (grain < soft_grain) {
							node->customValues[0] += 2.0 * node->customValues[grain + pf::SOLVER_ALLEN_CAHN];
						}
						else {
							node->customValues[0] += node->customValues[grain + pf::SOLVER_ALLEN_CAHN];
						}
			}
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 3000;
	inf.settings.disperse_settings.Nx = 100;
	inf.settings.disperse_settings.Ny = 100;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 5e-3;
	inf.settings.disperse_settings.dx = 1.0;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "back_ground");
	inf.materialSystem.matrix_phase.set(0, 0);

	// 形核定义
	RAND_init;
	double R2_limit = (1.4 * radius) * (1.4 * radius);
	for (int grain = 0; grain < grains_number; grain++) {
		bool is_init_grain = false;
		do
		{
			is_init_grain = true;
			double rand_x = RAND_0_1 * (inf.settings.disperse_settings.Nx - 2 * radius + 6) + radius - 3, 
				rand_y = RAND_0_1 * (inf.settings.disperse_settings.Ny - 2 * radius + 6) + radius - 3;
			for (auto geo = inf.nucleationBox.geometricRegion_box.begin(); geo < inf.nucleationBox.geometricRegion_box.end(); geo++) {
				double R2 = (rand_x - geo->ellipSolid.core.x) * (rand_x - geo->ellipSolid.core.x) 
					+ (rand_y - geo->ellipSolid.core.y) * (rand_y - geo->ellipSolid.core.y);
				if (R2 < R2_limit)
					is_init_grain = false;
			}
			if (is_init_grain) {
				pf::GeometricRegion geo;
				geo.init(pf::Geometry::Geo_Ellipsoid, 0, grain);
				geo.ellipSolid.set_core(rand_x, rand_y, 0.0);
				geo.ellipSolid.set_radius(radius, radius, 0.0);
				geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + grain, 1.0);
				inf.nucleationBox.geometricRegion_box.push_back(geo);
			}
		} while (!is_init_grain);
		cout << "grain " << to_string(grain) << " has been init !" << endl;
	}
	return inf;
}