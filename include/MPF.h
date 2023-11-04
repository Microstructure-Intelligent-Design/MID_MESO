#pragma once
#include "baseTools/baseTools.h"
#include "modules/modules.h"

using namespace std;
namespace pf {
	class MPF
	{
	public:
		MPF() {};
		~MPF() {};
		void init_Modules(Information information);
		void init_SimulationMesh();
		void output_before_loop();
		void output_after_loop();
		void relaxation_interface(int relaxation_steps, int out_step, bool is_output_file = true, bool adjust_phi_0_1 = false, bool normalized_phi = false);

		//used in solve loop
		void init_mesh_data(int istep, bool init_x = true, bool init_potential = true, bool init_phase_potential = true);
		void prepare_physical_parameters_in_mesh(int istep);
		void evolve_phase_evolution_equation(int istep, bool adjust_phi_0_1 = false);
		void evolve_temperature_evolution_equation(int istep);
		void evolve_phase_concentration_evolution_equation(int istep);
		void phaseFraction_assignment_and_prepare_cFlag(int istep, bool normalized_phi = false);
		void phaseFraction_assignment(int istep, bool normalized_phi = false);
		void phaseConcentration_assignment_with_cFlag(int istep, bool adjust_con_0_1 = true);
		void temperature_assignment(int istep);
		void output_in_loop(int istep, string custom_log = "");
		void data_file_write_in_loop(int istep, bool write_now = false);

		void automatically_adjust_dt(int istep, int start_treatment_step = 1000, int delt_step = 500, double MAX_SCALE = 1e6
			, double MIN_Mob_SCALE = SYS_EPSILON, bool isReduceOutput = false);

		void write_scalar_dfint_dphi_A_B(ofstream& fout, int alphaIndex, int betaIndex);
		void write_scalar_dfbulk_dphi_A_B(ofstream& fout, int alphaIndex, int betaIndex);
		void fill_phase_con_with_value(int fill_timestep, int phaseIndex, int conIndex, double conValue);
		void fill_phase_temperature_with_value(int fill_timestep, int phaseIndex, double temperature);
		void add_customFlag_to_allnodes(int index, int flag);
		void add_customValue_to_allnodes(int index, double value);
		void add_customVec3_to_allnodes(int index, Vector3 value);

		//used in exit
		void exit_MPF();

		bool conserved_the_aim_substances(vector<int> comp, int istep, int delt_step = 500, bool change_standard_x = false);
		//< Model
		Information information;
		//< Model
		Mechanics mechanics;
		InterfaceEnergy interfaceEnergy;
		Thermodynamics thermodynamic;
		Kinetics kinetics;
		Nucleation nucleations;
		ElectricField electricField;
		MagneticField magneticField;
		OptimizationAlgorithm optimizationAlgorithm;
		FluidField fluidField;
		//< The numerical space
		FieldStorage_forPhaseNode phaseMesh;

		//used in init
		void simulation_init_default();
		void simulation_init_byDataFile();

		double_box conserved_comp;
	private:
		void calculate_energyDensity_and_potential(PhaseNode& node);
		void write_customValue_in_node(ofstream& fout, int index, string name);
		void write_customFlag_in_node(ofstream& fout, int index, string name);
		void write_customVec3_in_node(ofstream& fout, int index, string name);
		License license;
	};
};