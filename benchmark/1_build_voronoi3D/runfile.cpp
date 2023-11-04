#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\1_build_voronoi3D\\")
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

			simulation.phaseFraction_assignment(istep);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static double Xi_abc(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma, pf::Info_DynamicCollection& inf) {
	return 10.0;
}
pf::Information settings() {
	pf::Information inf;
	// ��ֵ��ɢʱ��
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 2000;
	inf.settings.disperse_settings.int_width = 4;
	inf.settings.disperse_settings.Nx = 60;
	inf.settings.disperse_settings.Ny = 60;
	inf.settings.disperse_settings.Nz = 60;
	inf.settings.disperse_settings.dt = 1e-3;
	// ��������
	inf.settings.details_settings.int_grad = Int_Gradient::Steinbach_1996;
	inf.settings.details_settings.int_pot = Int_Potential::Nestler_Obstacle;
	inf.settings.details_settings.phi_incre_limit = 1e-3;
	inf.settings.details_settings.OMP_thread_counts = 5;
	// �ļ��������
	inf.settings.file_settings.file_output_step = 200;
	inf.settings.file_settings.screen_output_step = 50;
	inf.settings.file_settings.isPhiOutput = true;
	inf.settings.file_settings.is_write_datafile = false;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	// ������ϵ����
	inf.materialSystem.phases.add_Phase(0, "Grains");
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.functions.Xi_abc = Xi_abc;
	// �κ˶���
	vector<int> phase_property;
	phase_property.push_back(0);
	vector<double> phase_weight;
	phase_weight.push_back(1);
	pf::XNode x;
	inf.generate_voronoi_structure(Vector3(0, 0, 0), Vector3(60, 60, 60), 0, 12, 0, phase_property, phase_weight, x, 0.0);
	return inf;
}