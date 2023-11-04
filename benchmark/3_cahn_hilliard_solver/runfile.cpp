#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_cahn_hilliard_solver\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
static int comps_number = 1;
static double Mobility(pf::PhaseNode& node, int grain, int grains_start_index) {
	double m = 1.0;
	return m;
}
static double dF_dc(pf::PhaseNode& node, int grain, int grains_start_index) {
	double dF_dc = 0.0;
	dF_dc = node.customValues[grains_start_index + grain];
	return dF_dc;
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
		pf::CahnHilliardSolver cahnHilliard(simulation.phaseMesh);
		cahnHilliard.init_field(comps_number, pf::SOLVER_CAHN_HILLIARD);
		simulation.add_customValue_to_allnodes(cahnHilliard.comps_start_index, 0.001);
		cahnHilliard.set_dF_dc_func(dF_dc);
		cahnHilliard.set_Mobility_func(Mobility);
		cahnHilliard.set_Source_func(Source);
		simulation.information.settings.file_settings.isCustomValueOutput = true;
		for (int index = cahnHilliard.comps_start_index; index < cahnHilliard.comps_start_index + comps_number; index++) {
			simulation.information.settings.file_settings.customValue_output.add_string(index, "comp" + to_string(index - cahnHilliard.comps_start_index));
		}
		// define a solver
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			// solver cal
			double max_variation = cahnHilliard.solve_one_step(istep, simulation.information.settings.disperse_settings.dt, DifferenceMethod::FIVE_POINT);
			// data output
			simulation.output_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 5000;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 0.1;
	inf.settings.disperse_settings.dx = 1.0;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.screen_output_step = 100;
	inf.settings.file_settings.file_output_step = 100;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "back_ground");
	inf.materialSystem.matrix_phase.set(0, 0);

	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Ellipsoid);
	geo.ellipSolid.set_core(25.0, 25.0, 0.0);
	geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
	geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD + 0, 0.999);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}