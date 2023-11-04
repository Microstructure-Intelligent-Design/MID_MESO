#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\2_const_driving_force\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
int main(int argc, char* argv[]) {
	///< main program
	{
		///< init simulation by settings
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;

		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_evolution_equation(istep);

			simulation.evolve_phase_concentration_evolution_equation(istep);

			simulation.phaseFraction_assignment_and_prepare_cFlag(istep);

			simulation.phaseConcentration_assignment_with_cFlag(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

enum PHASE{Alpha, Beta};
static void energy(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
	if (p.phaseProperty == PHASE::Alpha)
		p.chemEnergyDensity = -0.1;
	else if (p.phaseProperty == PHASE::Beta)
		p.chemEnergyDensity = 0;
	return;
}
pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 10000;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 1;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-2;
	inf.settings.disperse_settings.int_width = 8;
	inf.settings.disperse_settings.x_bc = BoundaryCondition::ADIABATIC;
	// 其他功能
	inf.settings.details_settings.int_grad = Int_Gradient::Steinbach_G2009;
	inf.settings.details_settings.int_pot = Int_Potential::Steinbach_P2009;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 100;
	inf.settings.file_settings.screen_output_step = 100;
	inf.settings.file_settings.isPhiOutput = true;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(PHASE::Alpha, "alpha");
	inf.materialSystem.phases.add_Phase(PHASE::Beta, "beta");
	inf.materialSystem.matrix_phase.set(0, PHASE::Alpha);
	inf.materialSystem.functions.Energy = energy;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Polyhedron, 0, 1, PHASE::Beta, 0.0);
	geo.polyhedron.set_a_point_inside_polyhedron(Point(30.0, 0.0, 0.0));
	geo.polyhedron.add_surf(Point(25.0, 0.0, 0.0), Point(25.0, 1.0, 0.0), Point(25.0, 0.0, 1.0));
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}