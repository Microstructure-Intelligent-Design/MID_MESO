#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\2_multi_grains_juction\\")
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

			simulation.phaseFraction_assignment(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	return 1.0;
}

static double Xi_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	if ((alpha.index == 1 && beta.index == 0) || (alpha.index == 0 && beta.index == 1))
		return 1.0;
	else if((alpha.index == 2 && beta.index == 0) || (alpha.index == 0 && beta.index == 2))
		return 1.0;
	else if ((alpha.index == 2 && beta.index == 1) || (alpha.index == 1 && beta.index == 2))
		return 1.0;
	return 1.0;
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 10000;
	inf.settings.disperse_settings.Nx = 100;
	inf.settings.disperse_settings.Ny = 100;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-2;
	inf.settings.disperse_settings.y_bc = BoundaryCondition::ADIABATIC;
	// 其他功能
	inf.settings.details_settings.phi_incre_limit = 1e-2;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 1000;
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.isPhaseIndexsOutput = true;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "Grains");
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.functions.Mobility = mobility;
	inf.materialSystem.functions.Xi_ab = Xi_ab;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Polyhedron, 0, 1, 0, 0.0);
	geo.polyhedron.set_a_point_inside_polyhedron(Point(50.0, 25.0, 0.0));
	geo.polyhedron.add_surf(Point(0.0, 0.0, 0.0), Point(100.0, 0.0, 0.0), Point(100.0, 0.0, 1.0));
	geo.polyhedron.add_surf(Point(100.0, 0.0, 0.0), Point(100.0, 50.0, 0.0), Point(100.0, 50.0, 1.0));
	geo.polyhedron.add_surf(Point(100.0, 50.0, 0.0), Point(0.0, 50.0, 0.0), Point(0.0, 50.0, 1.0));
	geo.polyhedron.add_surf(Point(0.0, 50.0, 0.0), Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0));
	inf.nucleationBox.geometricRegion_box.push_back(geo);

	pf::GeometricRegion geo2(pf::Geometry::Geo_Ellipsoid, 0, 2, 0, 0.0);
	geo2.ellipSolid.set_core(50.0, 50.0, 0.0);
	geo2.ellipSolid.set_radius(20.0, 14.0, 0.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo2);

	return inf;
}