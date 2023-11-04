#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\square_cavity_flow\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
static void boundary_U(pf::VectorNode& u, pf::PhaseNode& down_node, pf::PhaseNode& up_node, int Nx, int Ny, int Nz) {
	if (u._x == 0) {
		u.vals[0] = 0.0;
	}
	if (u._y == 0) {
		u.vals[0] = 0.0;
	}
	if (u._x == Nx - 1) {
		u.vals[0] = 0.0;
	}
	if (u._y == Ny - 1) {
		u.vals[0] = 1.0;
	}

}
static void boundary_V(pf::VectorNode& v, pf::PhaseNode& down_node, pf::PhaseNode& up_node, int Nx, int Ny, int Nz) {
	if (v._y == 0) {
		v.vals[0] = 0.0;
	}
	else if (v._y == Ny - 1) {
		v.vals[0] = 0.0;
	}
	else if (v._x == 0) {
		v.vals[0] = 0.0;
	}
	else if (v._x == Nx - 1) {
		v.vals[0] = 0.0;
	}
}
static void boundary_W(pf::VectorNode& w, pf::PhaseNode& down_node, pf::PhaseNode& up_node, int Nx, int Ny, int Nz) {
	w.vals[0] = 0.0;
}
static void boundary_main(pf::PhaseNode& node, int Nx, int Ny, int Nz) {
	node.velocityValues.volume_force.set_to_zero();
}

int main(int argc, char* argv[]) {
	///< main program
	{
		///< init simulation by settings
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Info_Settings& set = simulation.information.settings;
		// define a solver
		simulation.output_before_loop();
		simulation.fluidField.set_boundary_condition_for_domain_U(boundary_U);
		simulation.fluidField.set_boundary_condition_for_domain_V(boundary_V);
		simulation.fluidField.set_boundary_condition_for_domain_W(boundary_W);
		simulation.fluidField.set_boundary_condition_for_main_domain(boundary_main);
		simulation.fluidField.init_pressure_in_fluid(1.0);
		simulation.fluidField.do_boundary_condition();
		double fluid_dt = 0.1;
		for (int istep = set.disperse_settings.begin_step; istep <= set.disperse_settings.end_step; istep++) {
			simulation.fluidField.evolve_momentum_equation(fluid_dt);

			double MAX_ABS_PRESSURE_CHANGE = simulation.fluidField.do_pressure_correction(3e-4, 1.0, false, 100);

			simulation.fluidField.correcting_velocity_field(fluid_dt);

			simulation.fluidField.do_boundary_condition();

			simulation.fluidField.assign_velocity_to_main_domain();

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
	inf.settings.disperse_settings.end_step = 2000;
	inf.settings.disperse_settings.x_bc = BoundaryCondition::ADIABATIC;
	inf.settings.disperse_settings.y_bc = BoundaryCondition::ADIABATIC;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1.0;

	// 文件输入输出
	inf.settings.file_settings.file_output_step = 100;
	inf.settings.file_settings.screen_output_step = 100;
	inf.settings.file_settings.isFluidFieldOutput = true;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";

	// 材料体系定义
	inf.materialSystem.is_fluid_on = true;
	inf.materialSystem.fluidFieldMask_settings.density = 1.0;
	inf.materialSystem.fluidFieldMask_settings.viscosity = 1.0;
	inf.materialSystem.fluidFieldMask_settings.fluid_phase_box.push_back(0);
	inf.materialSystem.phases.add_Phase(0, "grain");
	inf.materialSystem.matrix_phase.set(0, 0);

	return inf;
}