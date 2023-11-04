#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_electirc_field\\")
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

			simulation.electricField.setElectricChargeDensity();

			simulation.electricField.solve(1e-4, 10000);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}


pf::Information settings() {
	pf::Information inf;
	// ��ֵ��ɢʱ��
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 0;
	inf.settings.disperse_settings.Nx = 100;
	inf.settings.disperse_settings.Ny = 100;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-2;
	// �ļ��������
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 1000;
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.isElectricFieldOutput = true;
	// ������ϵ����
	inf.materialSystem.phases.add_Phase(0, "grain");
	inf.materialSystem.phases[0].conductivity = 1.0;

	inf.materialSystem.matrix_phase.set(0, 0);

	inf.materialSystem.is_electrics_on = true;
	inf.materialSystem.electricFieldMask_settings.init_electric_potential = 0.0;
	inf.materialSystem.electricFieldMask_settings.electric_potential_phase_index.add_double(1, 1.0);
	// �κ˶���
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Polyhedron, 0, 1, 0, 0.0);
	geo.polyhedron.set_a_point_inside_polyhedron(Point(50.0, 50.0, 0.0));
	geo.polyhedron.add_surf(Point(25.0, 25.0, 0.0), Point(75.0, 25.0, 0.0), Point(75.0, 25.0, 1.0));
	geo.polyhedron.add_surf(Point(75.0, 25.0, 0.0), Point(75.0, 75.0, 0.0), Point(75.0, 75.0, 1.0));
	geo.polyhedron.add_surf(Point(75.0, 75.0, 0.0), Point(25.0, 75.0, 0.0), Point(25.0, 75.0, 1.0));
	geo.polyhedron.add_surf(Point(25.0, 75.0, 0.0), Point(25.0, 25.0, 0.0), Point(25.0, 75.0, 1.0));
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}