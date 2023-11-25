#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
#include "DATA_BASE.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\dendrite\\")
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

			simulation.evolve_temperature_evolution_equation(istep);

			simulation.phaseFraction_assignment(istep);

			simulation.temperature_assignment(istep);

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
	inf.settings.disperse_settings.Nx = 300;
	inf.settings.disperse_settings.Ny = 300;
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 5000;
	inf.settings.disperse_settings.dx = 0.03;
	inf.settings.disperse_settings.int_width = 4 * inf.settings.disperse_settings.dx;
	inf.settings.disperse_settings.dt = 1e-4;
	// 其他功能
	inf.settings.details_settings.int_grad = Int_Gradient::Int_GCustom;
	inf.settings.details_settings.int_pot = Int_Potential::Int_PCustom;
	inf.settings.details_settings.difference_method = DifferenceMethod::NINE_POINT;
	inf.settings.details_settings.OMP_thread_counts = 10;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.screen_output_step = 100;
	inf.settings.file_settings.file_output_step = 100;
	inf.settings.file_settings.dataFile_output_step = 100;
	inf.settings.file_settings.isTemperatureOutput = true;
	// 材料体系定义
	inf.materialSystem.init_temperature = 0.0;
	inf.materialSystem.phases.add_Phase(DATABASE::SYS_Solidification::Get_Phase_Structure(DATABASE::SYS_Solidification::Solid));
	inf.materialSystem.phases.add_Phase(DATABASE::SYS_Solidification::Get_Phase_Structure(DATABASE::SYS_Solidification::PHASES(VIRTUAL_PHASE)));
	inf.materialSystem.matrix_phase.set(0, VIRTUAL_PHASE);
	inf.materialSystem.functions.MolarVolume = DATABASE::SYS_Solidification::MolarVolume;
	inf.materialSystem.functions.Mobility = DATABASE::SYS_Solidification::mobility;
	inf.materialSystem.functions.dfint_dphi = DATABASE::SYS_Solidification::dfint_dphi;
	inf.materialSystem.functions.HeatSource = DATABASE::SYS_Solidification::heatSource;
	// 形核定义
	{
		pf::GeometricRegion geo;
		geo.init(Geometry::Geo_Ellipsoid, 0, 1, DATABASE::SYS_Solidification::Solid, 0.0);
		geo.ellipSolid.set_core(150, 150, 0);
		geo.ellipSolid.set_radius(4, 4, 0);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}

	return inf;
}