#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\2_interface_dissipation\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
static double dF_dc_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
	return node.customValues[comp_index + comp_start_index];
}
static double Mobility_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
	return 1.0;
}
static double Source_cal(pf::PhaseNode& node, int comp_index, int comp_start_index) {
	return 0.0;
}
int main(int argc, char* argv[]) {
	///< main program
	{
		///< init simulation by settings
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;
		//< Cahn-Hilliard solver
		pf::CahnHilliardSolver ch_solver(simulation.phaseMesh);
		ch_solver.init_field(1);
		ch_solver.set_dF_dc_func(dF_dc_cal);
		ch_solver.set_Mobility_func(Mobility_cal);
		ch_solver.set_Source_func(Source_cal);
		simulation.add_customValue_to_allnodes(pf::SOLVER_CAHN_HILLIARD, 0.1);
		
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_evolution_equation(istep, true);

			simulation.evolve_phase_concentration_evolution_equation(istep);
			ch_solver.solve_one_step(istep, set.dt);

			simulation.phaseFraction_assignment_and_prepare_cFlag(istep);

			simulation.phaseConcentration_assignment_with_cFlag(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static Vector3 normals(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	Vector3 normal;
	normal[0] = alpha.phi_grad[0] * beta.phaseFraction - alpha.phaseFraction * beta.phi_grad[0];
	normal[1] = alpha.phi_grad[1] * beta.phaseFraction - alpha.phaseFraction * beta.phi_grad[1];
	normal[2] = alpha.phi_grad[2] * beta.phaseFraction - alpha.phaseFraction * beta.phi_grad[2];
	//normal[0] = std::sqrt(alpha.phi_grad[0] * beta.phi_grad[0]);
	//normal[1] = std::sqrt(alpha.phi_grad[1] * beta.phi_grad[1]);
	//normal[2] = std::sqrt(alpha.phi_grad[2] * beta.phi_grad[2]);
	double length = sqrt(normal * normal);
	if (length < SYS_EPSILON) {
		normal[0] = 0.0;
		normal[1] = 0.0;
		normal[2] = 0.0;
	}
	else {
		normal[0] = normal[0] / length;
		normal[1] = normal[1] / length;
		normal[2] = normal[2] / length;
	};
	return normal;
};

const int phase_property = 0, component = 0;
static void chemicalDiffusivity(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
	p.kinetics_coeff.set(component, component, 1.0);
	return;
}

static void intphaseReaction(pf::PhaseNode& n, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	double reaction_rate = 0.1;
	Vector3 norm = normals(alpha, beta, inf);
	for (auto comp1 = alpha.x.begin(); comp1 < alpha.x.end(); comp1++)
		for (auto comp2 = beta.x.begin(); comp2 < beta.x.end(); comp2++)
			if (comp1->index == comp2->index) {
				double pre = reaction_rate * norm.abs() * (comp2->value - comp1->value);
				if (pre > 0 && comp2->value > X_EPSILON && comp1->value < (1.0 - X_EPSILON)) {
					comp1->ChemicalReactionFlux += pre / alpha.phaseFraction;
					comp2->ChemicalReactionFlux -= pre / beta.phaseFraction;
				}
				else if (pre < 0 && comp1->value > X_EPSILON && comp2->value < (1.0 - X_EPSILON)) {
					comp1->ChemicalReactionFlux += pre / alpha.phaseFraction;
					comp2->ChemicalReactionFlux -= pre / beta.phaseFraction;
				}
			}
	return;
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 50000;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 1;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-3;
	inf.settings.disperse_settings.int_width = 14;
	inf.settings.disperse_settings.x_bc = BoundaryCondition::ADIABATIC;
	// 其他功能
	inf.settings.details_settings.flux_model = PhaseFluxModel::IntDiff_ConGrad;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 1000;
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.isPhiOutput = true;
	inf.settings.file_settings.isPhaseConOutput = true;
	inf.settings.file_settings.customValue_output.add_string(pf::SOLVER_CAHN_HILLIARD, "CH_COMP");
	inf.settings.file_settings.isCustomValueOutput = true;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(phase_property, "Grains");
	inf.materialSystem.phases[phase_property].x.add_nodeEntry(component, "X");
	inf.materialSystem.sys_x.add_nodeEntry(component);
	inf.materialSystem.matrix_phase.set(0, phase_property);
	inf.materialSystem.matrix_phase.x.add_con(component, 0.1);
	inf.materialSystem.functions.ChemicalDiffusivity = chemicalDiffusivity;
	inf.materialSystem.functions.InterPhasesReaction = intphaseReaction;
	inf.materialSystem.functions.Normals = normals;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Polyhedron, 0, 1, phase_property, 0.0);
	geo.polyhedron.set_a_point_inside_polyhedron(Point(30.0, 0.0, 0.0));
	geo.polyhedron.add_surf(Point(25.0, 0.0, 0.0), Point(25.0, 1.0, 0.0), Point(25.0, 0.0, 1.0));
	geo.x.add_x(component, 0.9);
	geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD, 0.9);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}