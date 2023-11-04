#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("benchmark\\3_elastic_field\\")
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
		vector<Matrix6x6> stiffness;
		stiffness.push_back(simulation.information.materialSystem.phases[0].elastic_constants);
		simulation.mechanics.SetMAXElasticConstants(stiffness);
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);

			simulation.nucleations.nucleation(istep);

			simulation.init_mesh_data(istep);

			simulation.relaxation_interface(1000, 100, false);

			simulation.mechanics.SetEffectiveEigenStrains();
			simulation.mechanics.SetEffectiveElasticConstants();
			simulation.mechanics.initVirtualEigenstrain();
			//simulation.mechanics.Solve(1e-6, 0.001, 1000, true);
			simulation.mechanics.Solve2(1e-6, 1000, 0.01, true);

			simulation.output_in_loop(istep);

			simulation.data_file_write_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

static vStrain effectivePhaseEigenStrains(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
	vStrain eigenstrain;
	for (auto phase = n.begin(); phase < n.end(); phase++) {
		vStrain p_e;
		if (phase->index == 1) {
			p_e[0] = 0.01;
			p_e[1] = 0.01;
			p_e[2] = 0.01;
			double radien[] = { AngleToRadians(0.0), AngleToRadians(0.0) , AngleToRadians(0.0) };
			p_e.do_rotate(RotationMatrix::rotationMatrix(radien, RotationGauge::RG_XYX));
		}
		else if (phase->index == 0) {
			p_e[0] = 0.0;
			p_e[1] = 0.0;
			p_e[2] = 0.0;
			double radien[] = { AngleToRadians(0.0), AngleToRadians(0.0) , AngleToRadians(0.0) };
			p_e.do_rotate(RotationMatrix::rotationMatrix(radien, RotationGauge::RG_XYX));
		}
		eigenstrain += p_e * phase->phaseFraction;
	}
	return eigenstrain;
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
	inf.settings.disperse_settings.int_width = 8;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.file_settings.file_output_step = 1000;
	inf.settings.file_settings.screen_output_step = 1000;
	inf.settings.file_settings.isMechanicalFieldOutput = true;
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "grain");
	double Gm = 1.0, v = 0.3, Az = 3.0, C12 = v * Gm / (1 - 2 * v), C11 = 2 * Gm / Az + C12;
	inf.materialSystem.phases[0].elastic_constants(0, 0) = C11;
	inf.materialSystem.phases[0].elastic_constants(1, 1) = C11;
	inf.materialSystem.phases[0].elastic_constants(2, 2) = C11;
	inf.materialSystem.phases[0].elastic_constants(1, 0) = C12;
	inf.materialSystem.phases[0].elastic_constants(2, 0) = C12;
	inf.materialSystem.phases[0].elastic_constants(2, 1) = C12;
	inf.materialSystem.phases[0].elastic_constants(0, 1) = C12;
	inf.materialSystem.phases[0].elastic_constants(0, 2) = C12;
	inf.materialSystem.phases[0].elastic_constants(1, 2) = C12;
	inf.materialSystem.phases[0].elastic_constants(3, 3) = Gm;
	inf.materialSystem.phases[0].elastic_constants(4, 4) = Gm;
	inf.materialSystem.phases[0].elastic_constants(5, 5) = Gm;
	inf.materialSystem.matrix_phase.set(0, 0);
	inf.materialSystem.functions.EffectiveEigenStrains = effectivePhaseEigenStrains;
	inf.materialSystem.is_mechanics_on = true;
	inf.materialSystem.mechanics.effectiveElasticConstantsModel = EffectiveElasticConstantsModel::EEC_Khachaturyan;
	inf.materialSystem.mechanics.rotation_gauge = RotationGauge::RG_XYX;
	double radien[] = {AngleToRadians(0.0), AngleToRadians(0.0) , AngleToRadians(0.0) };
	inf.materialSystem.mechanics.grainRotation.add_matrix3(0, RotationMatrix::rotationMatrix(radien, RotationGauge::RG_XYX));
	inf.materialSystem.mechanics.grainRotation.add_matrix3(1, RotationMatrix::rotationMatrix(radien, RotationGauge::RG_XYX));
	// 形核定义
	pf::GeometricRegion geo;
	geo.init(pf::Geometry::Geo_Ellipsoid, 0, 1, 0, 0.0);
	geo.ellipSolid.set_core(50.0, 50.0, 0.0);
	geo.ellipSolid.set_radius(15.0, 15.0, 0.0);
	inf.nucleationBox.geometricRegion_box.push_back(geo);
	return inf;
}