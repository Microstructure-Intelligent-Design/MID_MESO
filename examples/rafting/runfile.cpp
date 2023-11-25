#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
#include "DATA_BASE.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\rafting\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
int main(int argc, char* argv[]) {
	///< main program
	{
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;

		simulation.nucleations.nucleation(0);

		// relax interface
		for (int istep = 0; istep <= 5000; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.init_mesh_data(istep);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_concentration_evolution_equation(istep);

			simulation.phaseConcentration_assignment_with_cFlag(istep);

			simulation.automatically_adjust_dt(istep, 500, 50, 100, 1e-2, true);
		}
		vector<Matrix6x6> vec_M6x6;
		vec_M6x6.push_back(Get_PhaseElasticConstants(DATABASE::SYS_Elastic::PHASES::Matrix));
		vec_M6x6.push_back(Get_PhaseElasticConstants(DATABASE::SYS_Elastic::PHASES::Precipitate));
		simulation.mechanics.SetMAXElasticConstants(vec_M6x6);
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.init_mesh_data(istep);

			simulation.mechanics.SetEffectiveEigenStrains();
			simulation.mechanics.SetEffectiveElasticConstants();
			simulation.mechanics.Solve(1e-4, 1e-2, 1000, false);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_concentration_evolution_equation(istep);

			simulation.phaseConcentration_assignment_with_cFlag(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);

			simulation.automatically_adjust_dt(istep, 500, 50, 100, 1e-2, true);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}
pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.Nx = 128;
	inf.settings.disperse_settings.Ny = 128;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1;
	inf.settings.disperse_settings.dx = 1.0;
	inf.settings.disperse_settings.int_width = 4 * inf.settings.disperse_settings.dx;
	inf.settings.disperse_settings.x_bc = BoundaryCondition::PERIODIC;
	inf.settings.disperse_settings.y_bc = BoundaryCondition::PERIODIC;
	inf.settings.disperse_settings.z_bc = BoundaryCondition::PERIODIC;
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 10000;
	// 其他功能
	inf.settings.details_settings.difference_method = DifferenceMethod::NINE_POINT;
	inf.settings.details_settings.flux_model = PhaseFluxModel::IntDiff_PotentialGrad;
	inf.settings.details_settings.con_incre_limit = 1e-2;
	inf.settings.details_settings.OMP_thread_counts = 10;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.isMechanicalFieldOutput = true;
	inf.settings.file_settings.isChemicalPotentialDivisionOutput = true;
	inf.settings.file_settings.screen_output_step = 200;
	inf.settings.file_settings.file_output_step = 500;
	//inf.settings.file_settings.phases_mark;		// 标记输出，标记特殊相
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(DATABASE::SYS_Elastic::Get_Phase_Structure(DATABASE::SYS_Elastic::PHASES::Matrix));
	inf.materialSystem.matrix_phase.set(0, DATABASE::SYS_Elastic::PHASES::Matrix);
	inf.materialSystem.matrix_phase.x.add_con(DATABASE::SYS_Elastic::CON, 0.0);
	inf.materialSystem.sys_x.add_nodeEntry(DATABASE::SYS_Elastic::CON, "CON");

	inf.materialSystem.is_mechanics_on = true;
	inf.materialSystem.mechanics.effectiveElasticConstantsModel = EffectiveElasticConstantsModel::EEC_Custom;
	inf.materialSystem.mechanics.rotation_gauge = RotationGauge::RG_ZXZ;
	inf.materialSystem.mechanics.appStrainMask[0] = true;
	inf.materialSystem.mechanics.effectiveAppliedStrain[0] = 1e-2;
	Matrix3x3 matrix_rotate;
	matrix_rotate.set_to_unity();
	Matrix3x3 prec_rotate;
	prec_rotate.set_to_unity();
	inf.materialSystem.mechanics.grainRotation.add_matrix3(0, matrix_rotate);
	inf.materialSystem.mechanics.grainRotation.add_matrix3(1, prec_rotate);

	inf.materialSystem.functions.Potential = DATABASE::SYS_Elastic::Potential;
	inf.materialSystem.functions.ChemicalMobility = DATABASE::SYS_Elastic::ChemicalMobility;
	inf.materialSystem.functions.EffectiveEigenStrains = DATABASE::SYS_Elastic::EffectivePhaseEigenStrains;
	inf.materialSystem.functions.EffectiveElasticConstants = DATABASE::SYS_Elastic::EffectiveElasticConstants;
	inf.materialSystem.functions.MolarVolume = DATABASE::SYS_Elastic::MolarVolume;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Ellipsoid, 0, 0, DATABASE::SYS_Elastic::PHASES::Matrix);
	geo.ellipSolid.set_core(64, 64, 0);
	geo.ellipSolid.set_radius(15, 15, 0);
	geo.x.add_x(DATABASE::SYS_Elastic::CON, 1.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}
