#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
#include "DATA_BASE.h"
//#include "include/DataTest.h"
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\Al_Cu_precipitation\\")
#else
#define CPP_FILE_PATH string("")
#endif;
using namespace pf;
using namespace std;
pf::Information settings();
int main(int argc, char* argv[]) {
	///< main program
	{
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.relaxation_interface(2000, 200, false);
			simulation.fill_phase_con_with_value(0, 0, materials::Al, 0.8);
			simulation.fill_phase_con_with_value(0, 0, materials::Cu, 0.2);
			simulation.fill_phase_con_with_value(0, 1, materials::Al, 0.7);
			simulation.fill_phase_con_with_value(0, 1, materials::Cu, 0.3);

			simulation.init_mesh_data(istep);

			simulation.mechanics.SetEffectiveEigenStrains();
			simulation.mechanics.SetEffectiveElasticConstants();
			simulation.mechanics.Solve(1e-3, 1e2, 1000, true);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_evolution_equation(istep, true);

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
pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dx = 1e-6;
	inf.settings.disperse_settings.int_width = 4e-6;
	inf.settings.disperse_settings.dt = 1e-3;
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 0;
	// 其他功能
	inf.settings.details_settings.con_average_range = 2;
	inf.settings.details_settings.flux_model = PhaseFluxModel::IntDiff_ConGrad;
	inf.settings.details_settings.flux_error_balance_coefficient = 0.8;
	// 文件输入输出
	inf.settings.file_settings.isDrivingForceOutput = true;
	inf.settings.file_settings.isInterfaceEnergyOutput = true;
	inf.settings.file_settings.isMechanicalFieldOutput = true;
	inf.settings.file_settings.isChemicalPotentialDivisionOutput = true;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.file_output_step = 1000;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(DATABASE::SYS_AL_CU::Get_Phase_Structure(DATABASE::SYS_AL_CU::PHASES::FCC_A1));
	inf.materialSystem.phases.add_Phase(DATABASE::SYS_AL_CU::Get_Phase_Structure(DATABASE::SYS_AL_CU::PHASES::ALCU_THETA));
	inf.materialSystem.sys_x.add_nodeEntry(materials::Al, "Al");
	inf.materialSystem.sys_x.add_nodeEntry(materials::Cu, "Cu");
	inf.materialSystem.init_temperature = 500 + 273.0;
	inf.materialSystem.R = DATABASE::R;
	inf.materialSystem.matrix_phase.set(0, DATABASE::SYS_AL_CU::FCC_A1);
	inf.materialSystem.matrix_phase.x.add_con(materials::Al, 0.8);
	inf.materialSystem.matrix_phase.x.add_con(materials::Cu, 0.2);
	inf.materialSystem.is_mechanics_on = true;
	inf.materialSystem.mechanics.rotation_gauge = RotationGauge::RG_XYX;
	double radian[] = { AngleToRadians(0.0), AngleToRadians(0.0), AngleToRadians(0.0) };
	inf.materialSystem.mechanics.grainRotation.add_matrix3(0, RotationMatrix::rotationMatrix(radian, RotationGauge::RG_XYX));
	inf.materialSystem.mechanics.grainRotation.add_matrix3(1, RotationMatrix::rotationMatrix(radian, RotationGauge::RG_XYX));
	inf.materialSystem.functions.Energy = DATABASE::SYS_AL_CU::Energy;
	inf.materialSystem.functions.Potential = DATABASE::SYS_AL_CU::Potential;
	inf.materialSystem.functions.MolarVolume = DATABASE::SYS_AL_CU::MolarVolume;
	inf.materialSystem.functions.EnergyMinimizerIterator = DATABASE::SYS_AL_CU::EnergyMinimizerIterator;
	inf.materialSystem.functions.InterPhasesReaction = DATABASE::SYS_AL_CU::IntphaseReaction;
	inf.materialSystem.functions.ChemicalDiffusivity = DATABASE::SYS_AL_CU::ChemicalDiffusivity;
	inf.materialSystem.functions.Xi_ab = DATABASE::SYS_AL_CU::xi_AB;
	inf.materialSystem.functions.Mobility = DATABASE::SYS_AL_CU::mobility;
	inf.materialSystem.functions.EffectiveEigenStrains = DATABASE::SYS_AL_CU::EffectivePhaseEigenStrains;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(Geometry::Geo_Ellipsoid, 0, 1, DATABASE::SYS_AL_CU::PHASES::ALCU_THETA, 500 + 273.0);
	geo.x.add_x(materials::Al, 0.8);
	geo.x.add_x(materials::Cu, 0.2);
	geo.ellipSolid.set_core(25.0, 25.0, 0.0);
	geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}