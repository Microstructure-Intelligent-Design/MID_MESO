#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_magnetic_field\\")
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

			simulation.relaxation_interface(1000, 100, false);

			simulation.init_mesh_data(istep);

			simulation.magneticField.setMagneticParameters();
			simulation.magneticField.solve_magnetizationIntensity(1e-4, 10000, true, 100);
			simulation.magneticField.solve_dMagnetizationIntensity_dPhi(1e-4, 10000, true, 100);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}
static void magnetizationIntensity(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
	n.magneticValues.mag_intensity[1] = n[0].phaseFraction;
}
static void dMagnetizationIntensity_dPhi(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
	double vec[] = { 0.0, 0.0, 0.0 };
	if (p.index == 0) {
		vec[1] = 1.0;
		n.magneticValues.dmag_intensity_dphi.add_vec(p.index, vec);
	}
	else {
		n.magneticValues.dmag_intensity_dphi.add_vec(p.index, vec);
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
	inf.settings.file_settings.file_output_step = 1000;
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.isMagneticFieldOutput = true;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "grain");
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.is_magnetics_on = true;
	inf.materialSystem.magneticFieldMask_settings.init_magnetic_potential = 0.0;
	inf.materialSystem.magneticFieldMask_settings.is_averaged = true;
	inf.materialSystem.magneticFieldMask_settings.average_value = 0.0;
	inf.materialSystem.functions.MagnetizationIntensity = magnetizationIntensity;
	inf.materialSystem.functions.dMagnetizationIntensitydPhi = dMagnetizationIntensity_dPhi;
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Ellipsoid, 0, 1, 0, 0.0);
	geo.ellipSolid.set_core(50.0, 50.0, 0.0);
	geo.ellipSolid.set_radius(20.0, 20.0, 0.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}