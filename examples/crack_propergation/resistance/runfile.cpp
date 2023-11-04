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

			//simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);

			simulation.mechanics.SetEffectiveElasticConstants();
			//simulation.mechanics.initVirtualEigenstrain();
			//simulation.mechanics.Solve2(1e-4, 10000, 1e-12, false);
			simulation.mechanics.Solve(1e-5, 1e7, 10000, false);

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
enum PHASES{Grain1, Grain2};
static pf::Matrix6x6 phase_elastic_constance(int phase_property) {
	pf::Matrix6x6 C;
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
	return C;
}
static double interpolation_func(double phi) {
	return phi * phi * (3.0 - 2.0 * phi);
}
static double dinterpolation_func_dphi(double phi) {
	return 6.0 * phi * (1.0 - phi);
}
static double Gc(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
	double gc = 0.0;
	switch (phase.phaseProperty)
	{
	case Grain1:
		gc = 1.0;
		break;
	case Grain2:
		gc = 2.0;
		break;
	default:
		gc = 1.0;
		break;
	}
	return gc;
}
static double Gc(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
	double gc = 0.0, sum_h = 0.0;
	for (auto phase = node.begin(); phase < node.end(); phase++) {
		if (phase->phaseProperty != GROW_PHASE) {
			double h = interpolation_func(phase->phaseFraction);
			sum_h += h;
			//< Gc of each phase
			gc += h * Gc(node, *phase, inf);
		}
	}
	if (sum_h < SYS_EPSILON)
		return 0.0;
	else
		return gc / sum_h;
}
static double dGc_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
	double dfunc_dPhi = dinterpolation_func_dphi(phase.phaseFraction), sum_h = 0.0, sum_Gc = 0.0, gc_phase = Gc(node, phase, inf);
	for (auto p = node.begin(); p < node.end(); p++)
		if (p->phaseProperty != GROW_PHASE) {
			double h = interpolation_func(p->phaseFraction);
			sum_h += h;
			sum_Gc += h * Gc(node, *p, inf);
		}
	if (sum_h < SYS_EPSILON)
		return 0.0;
	return dfunc_dPhi * gc_phase / sum_h - interpolation_func(phase.phaseFraction) * gc_phase / sum_h / sum_h * dfunc_dPhi;
}
static double xi_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	return Gc(node, node[Grain1], inf) * 0.1;
};
static double xi_abc(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma, pf::Info_DynamicCollection& inf) {
	return Gc(node, node[Grain1], inf) * 2.0;
};
// < external functions
//！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！

static double dfint_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
	double K = 9.0 / 64.0, int_width = 4.0e-6;
	if (phase.phaseProperty == GROW_PHASE) {
		// defined by users
		double gc = Gc(node, inf), dx = 1e-6;
		Vector3 delt_gc;
		delt_gc[0] = (Gc(node.get_neighbor_node(pf::Direction::x_down), inf) - Gc(node.get_neighbor_node(pf::Direction::x_up), inf)) / 2.0 / dx;
		delt_gc[1] = (Gc(node.get_neighbor_node(pf::Direction::y_down), inf) - Gc(node.get_neighbor_node(pf::Direction::y_up), inf)) / 2.0 / dx;
		delt_gc[2] = (Gc(node.get_neighbor_node(pf::Direction::z_down), inf) - Gc(node.get_neighbor_node(pf::Direction::z_up), inf)) / 2.0 / dx;
		double val = gc / intface_width * K - 2.0 * intface_width * (phase.phi_grad * delt_gc + gc * phase.laplacian);
		return val;
	}
	else {
		// defined by users
		double crack = 0.0;
		double abs_crackphase_grad = 0.0, crackphase_fraction = 0.0;
		for(auto p = node.begin(); p < node.end(); p++)
			if (p->phaseProperty == GROW_PHASE) {
				abs_crackphase_grad = abs(p->phi_grad * p->phi_grad);
				crackphase_fraction = p->phaseFraction;
			}
		crack = dGc_dphi(node, phase, inf) * (intface_width * abs_crackphase_grad + crackphase_fraction * K / intface_width);
		// copy from interfaceEnergy.h,		Steinbach1996
		double grad = 0.0;
		for (auto beta = node.begin(); beta < node.end(); beta++)
			if (beta->index != phase.index && beta->_flag && beta->phaseProperty != GROW_PHASE) {
				grad += 2.0 * int_width * xi_ab(node, phase, *beta, inf) * (beta->phi_grad * beta->phi_grad * phase.phaseFraction
					- phase.phi_grad * beta->phi_grad * beta->phaseFraction + phase.phaseFraction * beta->phaseFraction * beta->laplacian - beta->phaseFraction * beta->phaseFraction * phase.laplacian);
			}
		// copy from interfaceEnergy.h		Nestler_Obstacle
		double pot = 0.0;
		for (auto beta = node.begin(); beta < node.end(); beta++) {
			if (beta->index == phase.index || beta->_flag == pf_BULK/* || beta->phaseProperty == GROW_PHASE*/)
				continue;
			pot += 16.0 / int_width / PI / PI * xi_ab(node, phase, *beta, inf) * beta->phaseFraction;
			for (auto gamma = beta + 1; gamma < node.end(); gamma++) {
				if (gamma->index == phase.index || gamma->_flag == pf_BULK/* || gamma->phaseProperty == GROW_PHASE*/)
					continue;
				pot += xi_abc(node, phase, *beta, *gamma, inf) * beta->phaseFraction * gamma->phaseFraction / int_width;
			}
		}
		return crack + grad + pot;
	}
}

static void energy(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
	if (phase.phaseProperty == GROW_PHASE) {
		phase.elasEnergyDensity = -1.0 / 2.0 * (node.mechanicalValues.Strains * node.mechanicalValues.Stresses);
	}
	else {
		double crack_fraction = 0.0;
		for (auto p = node.begin(); p < node.end(); p++)
			if (p->phaseProperty == GROW_PHASE)
				crack_fraction = p->phaseFraction;
		phase.elasEnergyDensity = (crack_fraction - 1) * (crack_fraction - 1) * dinterpolation_func_dphi(phase.phaseFraction)
			* (phase_elastic_constance(phase.phaseProperty) * node.mechanicalValues.Strains * node.mechanicalValues.Strains);
	}
	return;
}

static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	if (alpha.index == beta.index)
		return 0.0;
	if (alpha.phaseProperty == GROW_PHASE || beta.phaseProperty == GROW_PHASE)
		return 1.0e-12;
	return 0.0;
}

static double interpolation_function_for_driving_force(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
	return 1.0;
}

static pf::Matrix6x6 EffectiveElasticConstants(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
	pf::Matrix6x6 elasticConstants;
	elasticConstants.set_to_zero();
	double sum_h = 0.0, crack_fraction = 0.0;
	for (auto phase = node.begin(); phase < node.end(); phase++)
		if (phase->phaseProperty != GROW_PHASE) {
			double h = interpolation_func(phase->phaseFraction);
			sum_h += h;
			elasticConstants += phase_elastic_constance(phase->phaseProperty) * h;
		}
		else
			crack_fraction = phase->phaseFraction;
	if (sum_h < SYS_EPSILON)
		return elasticConstants;
	else
		return elasticConstants * (1.0 / sum_h) * (1.0 - crack_fraction);
}

void settings() {
	// 方峙宣柊扮腎
	information.settings.disperse_settings.begin_step = 10000;
	information.settings.disperse_settings.end_step = 100000;
	information.settings.disperse_settings.int_width = 1.0 * 1e-6;  // for crack interface
	information.settings.disperse_settings.Nx = 50;
	information.settings.disperse_settings.Ny = 50;
	information.settings.disperse_settings.Nz = 1;
	information.settings.disperse_settings.y_bc = BoundaryCondition::ADIABATIC;
	information.settings.disperse_settings.dt = 0.01;
	information.settings.disperse_settings.dx = 1e-6;
	// 凪麿孔嬬
	information.settings.details_settings.int_grad = Int_Gradient::Int_GCustom;
	information.settings.details_settings.int_pot = Int_Potential::Int_PCustom;
	information.settings.details_settings.phi_incre_limit = 1e-3;
	information.settings.details_settings.difference_method = DifferenceMethod::FIVE_POINT;
	information.settings.details_settings.phases_grow_shringk_type = PHASES_GROW_SHRINK_TYPE::HALF_INT;
	// 猟周補秘 補竃
	information.settings.file_settings.NaN = 0.0;
	information.settings.file_settings.isPhiOutput = true;
	information.settings.file_settings.isMechanicalFieldOutput = true;
	information.settings.file_settings.isPhaseIndexsOutput = true;
	//-----------------------------------------------------------------------
	// adjust parameters
	information.settings.file_settings.isInterfaceEnergyOutput = true;
	information.settings.file_settings.isDrivingForceOutput = true;
	//-----------------------------------------------------------------------
	information.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	information.settings.file_settings.is_init_byDatafile = true;
	information.settings.file_settings.is_Datafile_init_by_path = false;
	information.settings.file_settings.restart_datafile_path = CPP_FILE_PATH + "DATA_init_resistance.dat";
	information.settings.file_settings.file_output_step = 500;
	information.settings.file_settings.screen_output_step = 500;
	information.settings.file_settings.is_write_datafile = true;
	information.settings.file_settings.dataFile_output_step = 1000;
	// 可創悶狼協吶
	information.materialSystem.is_mechanics_on = true;
	information.materialSystem.mechanics.effectiveElasticConstantsModel = EffectiveElasticConstantsModel::EEC_Custom;
	information.materialSystem.mechanics.rotation_gauge = RotationGauge::RG_XYX;
	information.materialSystem.mechanics.y_bc = BoundaryCondition::ADIABATIC;
	information.materialSystem.mechanics.appStrainMask[0] = true;
	information.materialSystem.mechanics.effectiveAppliedStrain[0] = 1e-4;
	double radian[] = {AngleToRadians(0.0), AngleToRadians(0.0), AngleToRadians(0.0)};
	information.materialSystem.mechanics.grainRotation.add_matrix3(0, RotationMatrix::rotationMatrix_XYX(radian));
	information.materialSystem.mechanics.grainRotation.add_matrix3(1, RotationMatrix::rotationMatrix_XYX(radian));
	information.materialSystem.mechanics.grainRotation.add_matrix3(2, RotationMatrix::rotationMatrix_XYX(radian));

	information.materialSystem.phases.add_Phase(Grain1, "Grain1");
	information.materialSystem.phases.add_Phase(Grain2, "Grain2");
	information.materialSystem.phases[Grain1].elastic_constants = phase_elastic_constance(Grain1);
	information.materialSystem.phases[Grain2].elastic_constants = phase_elastic_constance(Grain2);
	information.materialSystem.phases.add_Phase(GROW_PHASE, "Crack");
	information.materialSystem.phases[GROW_PHASE].elastic_constants.set_to_zero();
	information.materialSystem.matrix_phase.set(0, Grain1);

	information.materialSystem.functions.dfint_dphi = dfint_dphi;
	information.materialSystem.functions.Energy = energy;
	information.materialSystem.functions.Mobility = mobility;
	information.materialSystem.functions.Interpolation_function_for_driving_force = interpolation_function_for_driving_force;
	information.materialSystem.functions.EffectiveElasticConstants = EffectiveElasticConstants;

	// 侘宰協吶
	{
		GeometricRegion geo;
		geo.generate_step = 0;
		geo.geometryProperty = Geometry::Geo_Ellipsoid;
		geo.phaseIndex = 1;
		geo.phaseProperty = GROW_PHASE;
		geo.ellipSolid.set_core(25, 0, 0);
		geo.ellipSolid.set_radius(0, 3, 0);
		information.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		GeometricRegion geo;
		geo.generate_step = 0;
		geo.geometryProperty = Geometry::Geo_Polyhedron;
		geo.phaseIndex = 2;
		geo.phaseProperty = Grain2;
		geo.polyhedron.set_a_point_inside_polyhedron(25, 25, 25);
		geo.polyhedron.add_surf(pf::Point(0.0, 25.0, 0.0), pf::Point(50.0, 25.0, 0.0), pf::Point(0.0, 25.0, 1.0));
		information.nucleationBox.geometricRegion_box.push_back(geo);
	}
}