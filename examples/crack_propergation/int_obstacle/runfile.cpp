#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\crack_propergation\\")
#else
#define CPP_FILE_PATH string("")
#endif;
void settings();
pf::Information information;
int main(int argc, char* argv[]) {
	settings();
	///< main program
	{
		///< init simulation by settings
		MPF simulation;
		simulation.init_Modules(information);
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;

		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);

			//simulation.mechanics.SetEffectiveElasticConstants();
			//simulation.mechanics.initVirtualEigenstrain();
			//simulation.mechanics.Solve2();

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_evolution_equation(istep, true);

			simulation.phaseFraction_assignment(istep, false);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

//！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
//< external functions
static pf::Matrix6x6 phase_elastic_constance(int phase_property) {
	pf::Matrix6x6 C;
	switch (phase_property)
	{
	default:
		double v = 0.3, E = 210e9, Gm = E / 2.0 / (1 + v), Az = 1.0, C12 = 2 * v * Gm / (1 - 2 * v), C11 = 2 * Gm / Az + C12;
		C(0, 0) = C11;
		C(1, 1) = C11;
		C(2, 2) = C11;
		C(3, 3) = Gm;
		C(4, 4) = Gm;
		C(5, 5) = Gm;
		C(0, 1) = C12;
		C(1, 0) = C12;
		C(0, 2) = C12;
		C(2, 0) = C12;
		C(1, 2) = C12;
		C(2, 1) = C12;
		break;
	}
	return C;
}
// < external functions
//！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！

static double dfint_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
	double GC = 1.0;
	if (phase.phaseProperty == GROW_PHASE) {
		// defined by users
		// < well potential
		//double dw_dphi = 1.0 / 2.0 * phase.phaseFraction * phase.phaseFraction;
		// < obstacle potential
		double dw_dphi = 9.0 / 64.0;
		double int_incre = GC * (-2.0 * intface_width * phase.laplacian + dw_dphi / intface_width);
		double bulk_incre = (phase.phaseFraction - 1.0) * (node.mechanicalValues.Strains * node.mechanicalValues.Stresses);
		return int_incre + bulk_incre;
	}
	else {
		return 0.0;
	}
}

static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	if (alpha.index == beta.index)
		return 0.0;
	return 1.0e-12;
}

void settings() {
	// 方峙宣柊扮腎
	information.settings.disperse_settings.begin_step = 0;
	information.settings.disperse_settings.end_step = 2000;
	information.settings.disperse_settings.int_width = 1.0 * 1e-6;
	information.settings.disperse_settings.Nx = 50;
	information.settings.disperse_settings.Ny = 1;
	information.settings.disperse_settings.Nz = 1;
	information.settings.disperse_settings.y_bc = BoundaryCondition::ADIABATIC;
	information.settings.disperse_settings.dt = 1e-2;
	information.settings.disperse_settings.dx = 1e-6;
	// 凪麿孔嬬
	information.settings.details_settings.int_grad = Int_Gradient::Int_GCustom;
	information.settings.details_settings.int_pot = Int_Potential::Int_PCustom;
	// 猟周補秘 補竃
	information.settings.file_settings.isPhiOutput = true;
	information.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	information.settings.file_settings.file_output_step = 20;
	information.settings.file_settings.screen_output_step = 20;
	// 可創悶狼協吶
	information.materialSystem.phases.add_Phase(0, "Grains");
	information.materialSystem.phases[0].elastic_constants.set_to_zero();
	information.materialSystem.phases.add_Phase(GROW_PHASE, "Crack");
	information.materialSystem.phases[GROW_PHASE].elastic_constants = phase_elastic_constance(GROW_PHASE);
	information.materialSystem.matrix_phase.set(0, 0);
	information.materialSystem.functions.dfint_dphi = dfint_dphi;
	information.materialSystem.functions.Mobility = mobility;
	information.materialSystem.mechanics.effectiveElasticConstantsModel = EffectiveElasticConstantsModel::EEC_Khachaturyan;

	// 侘宰協吶
	{
		GeometricRegion geo;
		geo.generate_step = 0;
		geo.geometryProperty = Geometry::Geo_Ellipsoid;
		geo.phaseIndex = 1;
		geo.phaseProperty = GROW_PHASE;
		geo.ellipSolid.set_core(25, 0, 0);
		geo.ellipSolid.set_radius(0, 0, 0);
		information.nucleationBox.geometricRegion_box.push_back(geo);
	}
}