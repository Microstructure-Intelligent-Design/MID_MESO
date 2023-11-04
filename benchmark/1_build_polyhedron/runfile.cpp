#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\1_build_polyhedron\\")
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

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 0;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 50;
	inf.settings.disperse_settings.x_bc = BoundaryCondition::ADIABATIC;
	inf.settings.disperse_settings.y_bc = BoundaryCondition::ADIABATIC;
	inf.settings.disperse_settings.z_bc = BoundaryCondition::ADIABATIC;
	inf.settings.disperse_settings.dt = 1e-2;
	// 其他功能
	
	// 文件输入输出
	inf.settings.file_settings.isPhiOutput = true;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "Grains");
	inf.materialSystem.matrix_phase.set(0, 0);
	// 形核定义
	{
		pf::GeometricRegion geo;
		geo.init(Geometry::Geo_Polyhedron, 0, 1, 0, 0.0);
		geo.polyhedron.set_a_point_inside_polyhedron(25, 25, 25);
		geo.polyhedron.add_surf(Point(25, 45, 25), Point(10, 25, 25), Point(25, 25, 15));//0
		geo.polyhedron.add_surf(Point(25, 45, 25), Point(10, 25, 25), Point(25, 25, 35));//1
		geo.polyhedron.add_surf(Point(25, 5, 25), Point(10, 25, 25), Point(25, 25, 15));//2
		geo.polyhedron.add_surf(Point(25, 5, 25), Point(10, 25, 25), Point(25, 25, 35));//3
		geo.polyhedron.add_surf(Point(25, 45, 25), Point(25, 25, 15), Point(40, 25, 25));//4
		geo.polyhedron.add_surf(Point(25, 45, 25), Point(25, 25, 35), Point(40, 25, 25));//5
		geo.polyhedron.add_surf(Point(25, 5, 25), Point(25, 25, 15), Point(40, 25, 25));//6
		geo.polyhedron.add_surf(Point(25, 5, 25), Point(25, 25, 35), Point(40, 25, 25));//7
		double radian[] = {0.0, AngleToRadians(-60.0), AngleToRadians(30.0)};
		geo.ellipSolid.set_rotation_radian_and_rotation_gauge(radian, RotationGauge::RG_XYZ);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	

	return inf;
}