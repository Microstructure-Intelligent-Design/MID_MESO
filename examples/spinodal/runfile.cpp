#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
#include "DATA_BASE.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\spinodal\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
int main(int argc, char* argv[])  {
	///< main program
	{
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;

		simulation.output_before_loop();
		double iterate_rate = 1e-10;
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			// add noise
			if (istep == set.begin_step) {
				RAND_init;
				for (auto node = simulation.phaseMesh._mesh.begin(); node < simulation.phaseMesh._mesh.end(); node++)
					(*node)[0].x[materials::Ti].value = 0.3 + RAND_0_1 * 0.01;
			}
			
			simulation.init_mesh_data(istep);

			simulation.mechanics.SetEffectiveEigenStrains();
			simulation.mechanics.SetEffectiveElasticConstants();
			if (istep == set.begin_step)
				simulation.mechanics.initVirtualEigenstrain();
			iterate_rate = simulation.mechanics.Solve2(istep, 1.0e-5, 1000, iterate_rate, false);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_phase_concentration_evolution_equation(istep);

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
	inf.settings.disperse_settings.Nx = 256;
	inf.settings.disperse_settings.Ny = 256;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1.0e-3;
	inf.settings.disperse_settings.dx = 1.0e-6;
	inf.settings.disperse_settings.int_width = 4 * inf.settings.disperse_settings.dx;
	inf.settings.disperse_settings.x_bc = BoundaryCondition::PERIODIC;
	inf.settings.disperse_settings.y_bc = BoundaryCondition::PERIODIC;
	inf.settings.disperse_settings.z_bc = BoundaryCondition::PERIODIC;
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 10000;
	// 其他功能
	inf.settings.details_settings.difference_method = DifferenceMethod::NINE_POINT;
	inf.settings.details_settings.flux_model = PhaseFluxModel::IntDiff_PotentialGrad;
	inf.settings.details_settings.con_incre_limit = 1e-3;
	inf.settings.details_settings.OMP_thread_counts = 1;

	// 文件输入输出
	inf.settings.file_settings.isMechanicalFieldOutput = true;
	inf.settings.file_settings.isChemicalPotentialDivisionOutput = true;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.screen_output_step = 100;
	inf.settings.file_settings.file_output_step = 500;
	// 材料体系定义
	inf.materialSystem.mechanics.effectiveElasticConstantsModel = EffectiveElasticConstantsModel::EEC_Khachaturyan;
	inf.materialSystem.mechanics.x_bc = BoundaryCondition::PERIODIC;
	inf.materialSystem.mechanics.y_bc = BoundaryCondition::PERIODIC;
	inf.materialSystem.mechanics.z_bc = BoundaryCondition::PERIODIC;
	inf.materialSystem.phases.add_Phase(DATABASE::SYS_TIC_ZRC::Get_Phase_Structure(DATABASE::SYS_TIC_ZRC::PHASES::FCC_A1));
	inf.materialSystem.matrix_phase.set(0, DATABASE::SYS_TIC_ZRC::PHASES::FCC_A1);
	inf.materialSystem.matrix_phase.x.add_con(materials::Ti, 0.3);
	inf.materialSystem.sys_x.add_nodeEntry(materials::Ti, "Ti");

	inf.materialSystem.is_mechanics_on = true;
	inf.materialSystem.mechanics.rotation_gauge = RotationGauge::RG_ZXZ;
	double radian[] = {0.0, 0.0, 0.0};
	inf.materialSystem.mechanics.grainRotation.add_matrix3(0, RotationMatrix::rotationMatrix(radian, RotationGauge::RG_ZXZ));
	inf.materialSystem.mechanics.grainRotation.add_matrix3(1, RotationMatrix::rotationMatrix(radian, RotationGauge::RG_ZXZ));

	inf.materialSystem.functions.Potential = DATABASE::SYS_TIC_ZRC::Potential;
	inf.materialSystem.functions.ChemicalMobility = DATABASE::SYS_TIC_ZRC::ChemicalMobility;
	inf.materialSystem.functions.EffectiveEigenStrains = DATABASE::SYS_TIC_ZRC::EffectivePhaseEigenStrains;
	inf.materialSystem.functions.MolarVolume = DATABASE::SYS_TIC_ZRC::MolarVolume;
	
	return inf;
}
