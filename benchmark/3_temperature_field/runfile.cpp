#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_temperature_field\\")
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

			simulation.evolve_temperature_evolution_equation(istep);

			simulation.temperature_assignment(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static double heatSource(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
	return 0.0;
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
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 1000;
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.isTemperatureOutput = true;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "grain");
	inf.materialSystem.phases[0].heat_diffusivity = 1.0;
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.init_temperature = 0.0;
	inf.materialSystem.functions.HeatSource = heatSource;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Polyhedron, 0, 0, 0, 0.0);
	geo.polyhedron.set_a_point_inside_polyhedron(Point(50.0, 50.0, 0.0));
	geo.polyhedron.add_surf(Point(25.0, 25.0, 0.0), Point(75.0, 25.0, 0.0), Point(75.0, 25.0, 1.0));
	geo.polyhedron.add_surf(Point(75.0, 25.0, 0.0), Point(75.0, 75.0, 0.0), Point(75.0, 75.0, 1.0));
	geo.polyhedron.add_surf(Point(75.0, 75.0, 0.0), Point(25.0, 75.0, 0.0), Point(25.0, 75.0, 1.0));
	geo.polyhedron.add_surf(Point(25.0, 75.0, 0.0), Point(25.0, 25.0, 0.0), Point(25.0, 75.0, 1.0));
	geo.temperature = 1.0;
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}