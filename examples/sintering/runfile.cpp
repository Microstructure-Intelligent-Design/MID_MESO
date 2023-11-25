#define _CRT_SECURE_NO_WARNINGS
#include "../../include/MPF.h"
using namespace pf;
using namespace std;
#ifdef _WIN32
#define CPP_FILE_PATH string("examples\\sintering\\")
#else
#define CPP_FILE_PATH string("")
#endif;
pf::Information settings();
static int comps_number = 1, grains_number = 9;
static double Mobility(pf::PhaseNode& node, int grain, int grains_start_index) {
	double m = 0.0, con = node.customValues[grains_start_index + grain], D_vol = 0.04, D_vap = 0.002, D_surf = 16.0, D_GB = 1.6, 
		phi = con * con * con * (10.0 - 15.0 * con + 6.0 * con * con), sum = 0.0;
	for (int pIndex1 = 0; pIndex1 < grains_number; pIndex1++)
		for (int pIndex2 = 0; pIndex2 < grains_number; pIndex2++)
			if(pIndex1 != pIndex2)
				sum += node.customValues[pIndex1 + pf::SOLVER_ALLEN_CAHN] * node.customValues[pIndex2 + pf::SOLVER_ALLEN_CAHN];
	m = D_vol * phi + D_vap * (1.0 - phi) + D_surf * con * (1.0 - con) + D_GB * sum;
	return m;
}
static double dF_dc(pf::PhaseNode& node, int grain, int grains_start_index) {
	double dF_dc = 0.0, con = node.customValues[grains_start_index + grain], A = 16.0, B = 1.0, sum2 = 0.0, sum3 = 0.0, kappa = 5.0;
	for (int pIndex = 0; pIndex < grains_number; pIndex++) {
		double phi = node.customValues[pIndex + pf::SOLVER_ALLEN_CAHN];
		sum2 += phi * phi;
		sum3 += phi * phi * phi;
	}
	dF_dc = 2.0 * A * con * (1 - con) * (1 - con) - 2.0 * A * con * con * (1 - con) + B * (2.0 * con - 6.0 * sum2 + 4.0 * sum3) - 0.5 * kappa * node.cal_customValues_laplace(pf::SOLVER_CAHN_HILLIARD + grain, 0.5);
	return dF_dc;
}
static double dF_dphi_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
	double dfint_dphi = 0.0, con = node.customValues[pf::SOLVER_CAHN_HILLIARD], phi = node.customValues[pf::SOLVER_ALLEN_CAHN + grain_index], sum = 0.0, B = 1.0, kappa = 2.0;
	for (int pIndex = 0; pIndex < grains_number; pIndex++)
		sum += node.customValues[pf::SOLVER_ALLEN_CAHN + pIndex] * node.customValues[pf::SOLVER_ALLEN_CAHN + pIndex];
	dfint_dphi = B * (12.0 * (1.0 - con) * phi - 12.0 * (2.0 - con) * phi * phi + 12.0 * sum * phi) - kappa * node.cal_customValues_laplace(pf::SOLVER_ALLEN_CAHN + grain_index, 0.5);
	return dfint_dphi;
}
static double L_cal(pf::PhaseNode& node, int grain_index, int grain_start_index) {
	return 5.0;
}
int main(int argc, char* argv[]) {
	///< main program
	{
		///< init simulation by settings
		MPF simulation;
		simulation.init_Modules(settings());
		simulation.init_SimulationMesh();
		pf::Set_Disperse& set = simulation.information.settings.disperse_settings;
		// define a solver
		pf::CahnHilliardSolver cahnHilliard(simulation.phaseMesh);
		cahnHilliard.init_field(comps_number, pf::SOLVER_CAHN_HILLIARD);
		simulation.add_customValue_to_allnodes(cahnHilliard.comps_start_index, 0.001);
		cahnHilliard.set_dF_dc_func(dF_dc);
		cahnHilliard.set_Mobility_func(Mobility);
		// define a solver
		pf::AllenCahnSolver allenCahn(simulation.phaseMesh, grains_number, pf::SOLVER_ALLEN_CAHN);
		allenCahn.init_field();
		allenCahn.set_dF_dphi_func(dF_dphi_cal);
		allenCahn.set_L_func(L_cal);
		int output_step = 100;
		simulation.information.settings.file_settings.screen_output_step = output_step;
		simulation.information.settings.file_settings.file_output_step = output_step;
		simulation.information.settings.file_settings.isCustomValueOutput = true;
		for (int index = cahnHilliard.comps_start_index; index < cahnHilliard.comps_start_index + comps_number; index++)
			simulation.information.settings.file_settings.customValue_output.add_string(index, "comp" + to_string(index - cahnHilliard.comps_start_index));
		simulation.information.settings.file_settings.customValue_output.add_string(0, "grains");
		// define a solver
		simulation.output_before_loop();
		for (int istep = set.begin_step; istep <= set.end_step; istep++) {
			simulation.information.dynamicCollection.init_each_timeStep(istep, set.dt);
			simulation.add_customValue_to_allnodes(0, 0.0);
			simulation.nucleations.nucleation(istep);

			// solver cal
			double max_variation = 0.0, max_variation2 = 0.0;
			max_variation = cahnHilliard.solve_one_step(istep, simulation.information.settings.disperse_settings.dt, DifferenceMethod::FIVE_POINT, true);
			max_variation2 = allenCahn.solve_one_step(istep, simulation.information.settings.disperse_settings.dt, true);
			// data output
			if (istep % output_step == 0) {
				cout << "> cahnHilliard = " << max_variation << ", allenCahn = " << max_variation2 << endl;
				for (auto node = simulation.phaseMesh._mesh.begin(); node < simulation.phaseMesh._mesh.end(); node++)
					for (int index = 0; index < grains_number; index++)
						node->customValues[0] += node->customValues[index + pf::SOLVER_ALLEN_CAHN] * node->customValues[index + pf::SOLVER_ALLEN_CAHN];
			}
			simulation.output_in_loop(istep);
		}
		simulation.output_after_loop();
		simulation.exit_MPF();
	}
}

pf::Information settings() {
	pf::Information inf;
	// 数值离散时空
	inf.settings.disperse_settings.begin_step = 0;
	inf.settings.disperse_settings.end_step = 5000;
	inf.settings.disperse_settings.Nx = 100;
	inf.settings.disperse_settings.Ny = 100;
	inf.settings.disperse_settings.Nz = 1;
	inf.settings.disperse_settings.dt = 1e-4;
	inf.settings.disperse_settings.dx = 0.5;
	// 文件输入输出
	inf.settings.file_settings.working_folder_path = CPP_FILE_PATH + "data";
	inf.settings.details_settings.OMP_thread_counts = 10;
	// 细节
	// 材料体系定义
	inf.materialSystem.phases.add_Phase(0, "back_ground");
	inf.materialSystem.matrix_phase.set(0, 0);

	// 形核定义
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(50.0, 50.0, 0.0);
		geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD, 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 0, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(29.0, 50.0, 0.0);
		geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 1, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(71.0, 50.0, 0.0);
		geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 2, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(50.0, 29.0, 0.0);
		geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 3, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(50.0, 71.0, 0.0);
		geo.ellipSolid.set_radius(10.0, 10.0, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 4, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	double R = 5.0;
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(39.0, 39.0, 0.0);
		geo.ellipSolid.set_radius(R, R, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 5, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(61.0, 39.0, 0.0);
		geo.ellipSolid.set_radius(R, R, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 6, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(39.0, 61.0, 0.0);
		geo.ellipSolid.set_radius(R, R, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 7, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	{
		pf::GeometricRegion geo;
		geo.init(pf::Geometry::Geo_Ellipsoid);
		geo.ellipSolid.set_core(61.0, 61.0, 0.0);
		geo.ellipSolid.set_radius(R, R, 0.0);
		geo.customValues.add_double(pf::SOLVER_CAHN_HILLIARD , 0.999);
		geo.customValues.add_double(pf::SOLVER_ALLEN_CAHN + 8, 0.999);
		inf.nucleationBox.geometricRegion_box.push_back(geo);
	}
	return inf;
}