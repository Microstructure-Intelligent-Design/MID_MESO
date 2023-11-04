#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\1_build_polygon\\")
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

pf::Information settings() {
	pf::Information inf;
	// ��ֵ��ɢʱ��
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 0;
	inf.settings.disperse_settings.Nx = 50;
	inf.settings.disperse_settings.Ny = 50;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-2;
	// ��������
	
	// �ļ��������
	inf.settings.file_settings.isPhiOutput = true;
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	// ������ϵ����
	inf.materialSystem.phases.add_Phase(0, "Grains");
	inf.materialSystem.matrix_phase.set(0, 0);
	// �κ˶���
	{
		pf::GeometricRegion geo;
		geo.init(Geometry::Geo_Polyhedron, 0, 1, 0, 0.0);
		geo.polyhedron.set_a_point_inside_polyhedron(25, 25, 0);
		geo.polyhedron.add_surf(Point(25, 45, 0), Point(10, 25, 0), Point(10, 25, 1));
		geo.polyhedron.add_surf(Point(10, 25, 0), Point(25, 5, 0), Point(25, 5, 1));
		geo.polyhedron.add_surf(Point(25, 5, 0), Point(40, 25, 0), Point(40, 25, 1));
		geo.polyhedron.add_surf(Point(40, 25, 0), Point(25, 45, 0), Point(25, 45, 1));
		double radian[] = {0.0, AngleToRadians(30.0), 0.0};
		geo.ellipSolid.set_rotation_radian_and_rotation_gauge(radian, RotationGauge::RG_XZX);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	

	return inf;
}