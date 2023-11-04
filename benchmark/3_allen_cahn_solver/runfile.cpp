#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_allen_cahn_solver\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
static int grains_number = 2;
static double L(pf::PhaseNode& node, int grain, int grains_start_index) {
	double L = 5.0;
	return L;
}
static double dF_dphi(pf::PhaseNode& node, int grain, int grains_start_index) {
	double A = 1.0, B = 1.0, kappa = 1.0
		, dF_dphi = 0.0, phi = node.customValues[grain + grains_start_index], drivingforce = 0.0, sum = 0.0;
	for (int index = grains_start_index; index < grains_start_index + grains_number; index++)
		if (index != (grain + grains_start_index))
			sum += node.customValues[index] * node.customValues[index];
	if (grain == 1)
		drivingforce = -1.0;
	dF_dphi = -A * phi + B * phi * phi * phi + 2.0 * phi * sum - kappa * node.cal_customValues_laplace(grain + grains_start_index, 1.0, DifferenceMethod::FIVE_POINT) + phi * (1 - phi) * drivingforce;
	return dF_dphi;
}
static double Source(pf::PhaseNode& node, int grain, int grains_start_index) {
	double Source = 0.0;
	return Source;
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
		simulation.add_customValue_to_allnodes(allenCahn.grains_start_index, 1.0);
		allenCahn.set_dF_dphi_func(dF_dphi);
		allenCahn.set_L_func(L);
		allenCahn.set_Source_func(Source);
		simulation.information.settings.file_settings.isCustomValueOutput = true;
		for (int index = allenCahn.grains_start_index; index < allenCahn.grains_start_index + grains_number; index++) {
			simulation.information.settings.file_settings.customValue_output.add_string(index, "grain" + to_string(index - allenCahn.grains_start_index));
		}
		// define a solver
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.output_in_loop(istep);

			// solver cal
			double max_variation = allenCahn.solve_one_step(istep, simulation.information.settings.disperse_settings.dt);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 1000;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 0.01;
	inf.settings.disperse_settings.dx = 1.0;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 100;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "back_ground");
	inf.materialSystem.matrix_phase.set(0, 0);

	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Ellipsoid);
	geo.ellipSolid.set_core(25.0, 25.0, 0.0);
	geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
	geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 0, 0.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	pf::GeometricRegion geo1;
	geo1.init(pf::Geometry::Geo_Ellipsoid);
	geo1.ellipSolid.set_core(25.0, 25.0, 0.0);
	geo1.ellipSolid.set_radius(10.0, 10.0, 0.0);
	geo1.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 1, 1.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo1);
	return inf;
}