#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_poisson_solver\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
enum ElectricField{ElectricPotential, ElecChargeDensity_Conduction};
static void rhs_cal(pf::PhaseNode& node, int rhs_index) {
	node.customValues[rhs_index] = 0.0;
}
static void boundary(pf::PhaseNode& node, int lhs_index) {
	if (node[1].phaseFraction > Simulation_Num_Cut_Off)
		node.customValues[lhs_index] = 1.0;
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
		pf::PoissonEquationSolver electricfield_solver(simulation.phaseMesh, ElectricPotential, ElecChargeDensity_Conduction);
		electricfield_solver.init_field(0.0, 0.0);
		electricfield_solver.set_BoundaryCondition_calfunc(boundary);
		electricfield_solver.set_RHS_calfunc(rhs_cal);
		// define a solver
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);
			// solver cal
			electricfield_solver.solve_whole_domain(1e-4, 10000, true, 500);

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
	inf.settings.disperse_settings.Nx = 100;
	inf.settings.disperse_settings.Ny = 100;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-2;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.isCustomValueOutput = true;
	inf.settings.file_settings.customValue_output.add_string(ElectricPotential, "electric_potential");
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "grain");
	inf.materialSystem.matrix_phase.set(0, 0);

	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Polyhedron, 0, 1, 0, 0.0);
	geo.polyhedron.set_a_point_inside_polyhedron(Point(50.0, 50.0, 0.0));
	geo.polyhedron.add_surf(Point(25.0, 25.0, 0.0), Point(75.0, 25.0, 0.0), Point(75.0, 25.0, 1.0));
	geo.polyhedron.add_surf(Point(75.0, 25.0, 0.0), Point(75.0, 75.0, 0.0), Point(75.0, 75.0, 1.0));
	geo.polyhedron.add_surf(Point(75.0, 75.0, 0.0), Point(25.0, 75.0, 0.0), Point(25.0, 75.0, 1.0));
	geo.polyhedron.add_surf(Point(25.0, 75.0, 0.0), Point(25.0, 25.0, 0.0), Point(25.0, 75.0, 1.0));
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}