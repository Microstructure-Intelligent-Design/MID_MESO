#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\2_interface_study2D\\")
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

			simulation.evolve_phase_evolution_equation(istep, true);

			simulation.phaseFraction_assignment(istep, true);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static double Xi_abc(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma, pf::Info_DynamicCollection& inf) {
	return 0.0;
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 1000;
	inf.settings.disperse_settings.int_width = 10.0;
	inf.settings.disperse_settings.Nx = 100;
	inf.settings.disperse_settings.Ny = 100;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.y_bc = BoundaryCondition::ADIABATIC;
	inf.settings.disperse_settings.dt = 1e-2;
	// 其他功能
	inf.settings.details_settings.int_grad = Int_Gradient::Steinbach_G2009;
	inf.settings.details_settings.int_pot = Int_Potential::Steinbach_P2009;
	// 文件输入输出
	inf.settings.file_settings.isPhiOutput = true;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 100;
	inf.settings.file_settings.screen_output_step = 50;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "Grains");
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.functions.Xi_abc = Xi_abc;
	// 形核定义
	{
		GeometricRegion geo;
		geo.generate_step = 0;
		geo.geometryProperty = Geometry::Geo_Polyhedron;
		geo.phaseIndex = 1;
		geo.phaseProperty = 0;
		geo.temperature = 0.0;
		geo.polyhedron.set_a_point_inside_polyhedron(Point(25, 75, 0));
		geo.polyhedron.add_surf(Point(0, 50, 0), Point(50, 50, 0), Point(50, 50, 1));
		geo.polyhedron.add_surf(Point(50, 50, 0), Point(50, 100, 0), Point(50, 50, 1));
		geo.polyhedron.add_surf(Point(50, 100, 0), Point(0, 100, 0), Point(50, 100, 1));
		geo.polyhedron.add_surf(Point(0, 100, 0), Point(0, 50, 0), Point(0, 100, 1));
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}

	{
		GeometricRegion geo;
		geo.generate_step = 0;
		geo.geometryProperty = Geometry::Geo_Polyhedron;
		geo.phaseIndex = 2;
		geo.phaseProperty = 0;
		geo.temperature = 0.0;
		geo.polyhedron.set_a_point_inside_polyhedron(Point(75, 75, 0));
		geo.polyhedron.add_surf(Point(50, 50, 0), Point(100, 50, 0), Point(100, 50, 1));
		geo.polyhedron.add_surf(Point(100, 50, 0), Point(100, 100, 0), Point(100, 50, 1));
		geo.polyhedron.add_surf(Point(100, 100, 0), Point(50, 100, 0), Point(100, 100, 1));
		geo.polyhedron.add_surf(Point(50, 100, 0), Point(50, 50, 0), Point(50, 100, 1));
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}

	return inf;
}