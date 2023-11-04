#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\2_source_term\\")
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
		for (auto node = simulation.phaseMesh._mesh.begin(); node < simulation.phaseMesh._mesh.end(); node++) {
			node->velocityValues.velocity[0] = 0.5;
			node->velocityValues.velocity[1] = 0.5;
			node->velocityValues.velocity[2] = 0.0;
		}
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);

			simulation.prepare_physical_parameters_in_mesh(istep);

			simulation.evolve_temperature_evolution_equation(istep);

			simulation.temperature_assignment(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static double heatSource(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
	return node.velocityValues.velocity * node.tempValues.temperature_grad;
}

pf::Information settings() {
	pf::Information inf;
	// ��ֵ��ɢʱ��
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 5000;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-2;

	// �ļ��������
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 100;
	inf.settings.file_settings.screen_output_step = 100;
	inf.settings.file_settings.isFluidFieldOutput = true;
	inf.settings.file_settings.isTemperatureOutput = true;
	// ������ϵ����
	inf.materialSystem.phases.add_Phase(0, "grain");
	inf.materialSystem.phases[0].heat_diffusivity = 1.0;
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.init_temperature = 0.0;
	inf.materialSystem.functions.HeatSource = heatSource;
	// �κ˶���
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Ellipsoid, 0, 0, 0, 0.0);
	geo.ellipSolid.set_core(25, 25, 0);
	geo.ellipSolid.set_radius(10, 10, 0);
	geo.temperature = 1.0;
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}