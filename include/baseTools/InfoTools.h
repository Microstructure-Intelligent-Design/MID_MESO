#pragma once
#include "sysTool.h"
#include "PhaseNode.h"
#include "NucleationTools.h"
#include "RotationMatrix.h"

namespace pf {
	enum Note{Note_None, Note_Reactant, Note_Product};

	struct Info_DynamicCollection {
		void init_each_timeStep(int _istep, double dt) {
			istep = _istep;
			is_Nucleation = false;
			is_possible_Nucleation = false;
			real_time += dt * dt_scale;
			new_phase_index.clear();
			MAX_diffusionFlux = 0.0;
			MAX_reactionFlux = 0.0;
			MAX_phaseTransFlux = 0.0;
			MAX_phase_increment = 0.0;
			MAX_temp_increment = 0.0;
			MAX_ABS_DISPERSION = 0.0;
			MAX_ABS_dUdt = 0.0;
			MAX_ABS_dVdt = 0.0;
			MAX_ABS_dWdt = 0.0;
			MAX_ABS_dPRESSURE = 0.0;
		}
		void init_each_outputStep() {
			average_interface_phaseCon_num = 0;
			elastic_iterate_counts_once = 0;
			magneticField__iterate_counts_once = 0;
			electricField__iterate_counts_once = 0;
			fluidField_pressure__iterate_counts_once = 0;
			is_average_line_too_small = false;
			is_interface_move_too_quickly = false;
			is_DiffusionFlux_term_out_of_limit = false;
			is_PhaseTransitionFlux_term_out_of_limit = false;
			is_ChemicalReactionFlux_term_out_of_limit = false;
			is_cal_constituent_of_onePhase_by_energy_minimization_too_quickly = false;
		}
		///< used every step
		int istep;
		double real_time;
		bool is_Nucleation;
		bool is_possible_Nucleation;
		vector<int> new_phase_index;
		double MAX_phase_increment;
		double MAX_diffusionFlux;
		double MAX_reactionFlux;
		double MAX_phaseTransFlux;
		double MAX_temp_increment;
		///< used every output
		bool is_interface_move_too_quickly;
		double MAX_dphi_versus_limitPhi;

		bool is_DiffusionFlux_term_out_of_limit;
		bool is_PhaseTransitionFlux_term_out_of_limit;
		bool is_ChemicalReactionFlux_term_out_of_limit;
		double MAX_dcon_versus_limitCon;

		bool is_average_line_too_small;
		bool is_cal_constituent_of_onePhase_by_energy_minimization_too_quickly;

		double MAX_ABS_DISPERSION;
		double MAX_ABS_dUdt;
		double MAX_ABS_dVdt;
		double MAX_ABS_dWdt;
		double MAX_ABS_dPRESSURE;

		///< other data & funcs
		double dt_scale;
		double average_interface_reaction_phaseCon_rate;
		int average_interface_phaseCon_num;
		int cal_siteFraction_iterated_times;
		int cal_associated_solution_constituents_iterated_times;
		int elastic_iterate_counts_once;
		int electricField__iterate_counts_once;
		int magneticField__iterate_counts_once;
		int fluidField_pressure__iterate_counts_once;

		bool is_write_dataFiles;

		double_box user_defined_value;
		int_box user_defined_flag;

		PhaseNode inf_node;
		ConNode init_x_ratio;
		Info_DynamicCollection() {
			real_time = 0.0;
			MAX_phase_increment = 0.0;
			MAX_diffusionFlux = 0.0;
			MAX_reactionFlux = 0.0;
			MAX_phaseTransFlux = 0.0;
			MAX_temp_increment = 0.0;
			MAX_dphi_versus_limitPhi = 1.0;
			MAX_dcon_versus_limitCon = 1.0;
			MAX_ABS_DISPERSION = 0.0;
			MAX_ABS_dUdt = 0.0;
			MAX_ABS_dVdt = 0.0;
			MAX_ABS_dWdt = 0.0;
			MAX_ABS_dPRESSURE = 0.0;
			dt_scale = 1.0;
			average_interface_reaction_phaseCon_rate = 1.0;
			is_interface_move_too_quickly = false;
			is_write_dataFiles = false;
			is_DiffusionFlux_term_out_of_limit = false;
			is_PhaseTransitionFlux_term_out_of_limit = false;
			is_ChemicalReactionFlux_term_out_of_limit = false;
			is_average_line_too_small = false;
			is_cal_constituent_of_onePhase_by_energy_minimization_too_quickly = false;

			is_Nucleation = false;
			is_possible_Nucleation = false;
			average_interface_phaseCon_num = 0;
			cal_siteFraction_iterated_times = 0;
			cal_associated_solution_constituents_iterated_times = 0;
			elastic_iterate_counts_once = 0;
			electricField__iterate_counts_once = 0;
			magneticField__iterate_counts_once = 0;
			fluidField_pressure__iterate_counts_once = 0;
			istep = 0;
		}
		Info_DynamicCollection& operator=(const Info_DynamicCollection& n) {
			real_time = n.real_time;
			dt_scale = n.dt_scale;
			MAX_phase_increment = n.MAX_phase_increment;
			MAX_diffusionFlux = n.MAX_diffusionFlux;
			MAX_reactionFlux = n.MAX_reactionFlux;
			MAX_phaseTransFlux = n.MAX_phaseTransFlux;
			MAX_temp_increment = n.MAX_temp_increment;
			MAX_dphi_versus_limitPhi = n.MAX_dphi_versus_limitPhi;
			MAX_dcon_versus_limitCon = n.MAX_dcon_versus_limitCon;
			MAX_ABS_DISPERSION = n.MAX_ABS_DISPERSION;
			MAX_ABS_dUdt = n.MAX_ABS_dUdt;
			MAX_ABS_dVdt = n.MAX_ABS_dVdt;
			MAX_ABS_dWdt = n.MAX_ABS_dWdt;
			MAX_ABS_dPRESSURE = n.MAX_ABS_dPRESSURE;
			is_interface_move_too_quickly = n.is_interface_move_too_quickly;
			is_write_dataFiles = n.is_write_dataFiles;
			is_Nucleation = n.is_Nucleation;
			is_possible_Nucleation = n.is_possible_Nucleation;
			new_phase_index = n.new_phase_index;
			is_DiffusionFlux_term_out_of_limit = n.is_DiffusionFlux_term_out_of_limit;
			is_PhaseTransitionFlux_term_out_of_limit = n.is_PhaseTransitionFlux_term_out_of_limit;
			is_ChemicalReactionFlux_term_out_of_limit = n.is_ChemicalReactionFlux_term_out_of_limit;
			is_average_line_too_small = n.is_average_line_too_small;
			is_cal_constituent_of_onePhase_by_energy_minimization_too_quickly = n.is_cal_constituent_of_onePhase_by_energy_minimization_too_quickly;
			average_interface_phaseCon_num = n.average_interface_phaseCon_num;
			cal_siteFraction_iterated_times = n.cal_siteFraction_iterated_times;
			cal_associated_solution_constituents_iterated_times = n.cal_associated_solution_constituents_iterated_times;
			elastic_iterate_counts_once = n.elastic_iterate_counts_once;
			inf_node = n.inf_node;
			init_x_ratio = n.init_x_ratio;
			istep = n.istep;
			electricField__iterate_counts_once = n.electricField__iterate_counts_once;
			magneticField__iterate_counts_once = n.magneticField__iterate_counts_once;
			fluidField_pressure__iterate_counts_once = n.fluidField_pressure__iterate_counts_once;
			user_defined_value = n.user_defined_value;
			user_defined_flag = n.user_defined_flag;
			average_interface_reaction_phaseCon_rate = n.average_interface_reaction_phaseCon_rate;
			return *this;
		}
	};

	struct Info_energy_treatment {
		double_box x_limit_down;
		double_box x_limit_up;
		bool_box is_x_down_treatment;
		bool_box is_x_up_treatment;
		double treat_dx_range;
		double MAX_deltG;
		Info_energy_treatment() {
			treat_dx_range = 0;
			MAX_deltG = 0.0;
		}
		Info_energy_treatment& operator=(const Info_energy_treatment& p) {
			is_x_down_treatment = p.is_x_down_treatment;
			is_x_up_treatment = p.is_x_up_treatment;
			treat_dx_range = p.treat_dx_range;
			x_limit_down = p.x_limit_down;
			x_limit_up = p.x_limit_up;
			MAX_deltG = p.MAX_deltG;
			return *this;
		}
	};
	struct Info_Mechanics {
		// MatrixPhase: The precipitated phase is rotated based on the matrix phase
		// | 1  0  0 |
		// | 0  1  0 |
		// | 0  0  1 |
		// tensor1_matrix6     PhaseElasticConstants;                      ///< Reference elastic constants of each phase (relative to Cref), grain orientation not considered
		// PhaseElasticConstants[PrecipitatePhaseIndex]
		// tensor1_matrix6     PhaseCompliences;                           ///< Inverse of the PhaseElasticConstants, grain orientation not considered
		// PhaseCompliences[PrecipitatePhaseIndex]
		tensor1_matrix3     grainRotation;                             ///< Symmetry variants and Coordinate rotation
		// Variants[PrecipitatePhaseIndex]
		RotationGauge       rotation_gauge;

		EffectiveElasticConstantsModel effectiveElasticConstantsModel;
		//// BoundaryCondition
		// [] -> 1 - direction X, 2 - direction Y, 3 - direction Z
		// free boundary
		vector<bool> avgStrainMask;
		// apply stress
		vector<bool> loadStressMask;
		//vStress appliedStress;
		vStress effectiveAppliedStress;
		// apply strain
		vector<bool> appStrainMask;
		//vStrain appliedStrain;
		vStrain effectiveAppliedStrain;

		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		Info_Mechanics() {
			avgStrainMask.resize(3);
			loadStressMask.resize(3);
			appStrainMask.resize(3);
			rotation_gauge = RotationGauge::RG_ZXZ;
			effectiveElasticConstantsModel = EffectiveElasticConstantsModel::EEC_Khachaturyan;
			for (int ii = 0; ii < 3; ii++) {
				loadStressMask[ii] = false;
				avgStrainMask[ii] = true;
				appStrainMask[ii] = false;
			}
			x_bc = BoundaryCondition::PERIODIC;
			y_bc = BoundaryCondition::PERIODIC;
			z_bc = BoundaryCondition::PERIODIC;
		}
		Info_Mechanics& operator=(const Info_Mechanics& n) {
			//PhaseKappa = n.PhaseKappa;
			//PhaseElasticConstants = n.PhaseElasticConstants;
			//PhaseCompliences = n.PhaseCompliences;
			grainRotation = n.grainRotation;
			effectiveElasticConstantsModel = n.effectiveElasticConstantsModel;
			rotation_gauge = n.rotation_gauge;
			avgStrainMask = n.avgStrainMask;
			loadStressMask = n.loadStressMask;
			//appliedStress = n.appliedStress;
			effectiveAppliedStress = n.effectiveAppliedStress;
			appStrainMask = n.appStrainMask;
			//appliedStrain = n.appliedStrain;
			effectiveAppliedStrain = n.effectiveAppliedStrain;
			x_bc = n.x_bc;
			y_bc = n.y_bc;
			z_bc = n.z_bc;
			return *this;
		}
	};
	struct Info_NodeEntry {
		int index;
		int note;
		double value;
		double charge;
		string name;
		Info_NodeEntry() {
			index = 0;
			value = 0.0;
			charge = 0.0;
			note = Note::Note_None;
			name = "";
		}
		Info_NodeEntry& operator=(const Info_NodeEntry& ne) {
			index = ne.index;
			value = ne.value;
			charge = ne.charge;
			name = ne.name;
			note = ne.note;
			return *this;
		}
	};
	class Info_Node {
	public:
		Info_Node() {
			siteNum = 0.0;
			index = 0;
			name = "";
		}
		Info_NodeEntry& operator[](int index) {
			for (auto i = node.begin(); i < node.end(); ++i) {
				if (i->index == index) return(*i);
			}
			cout << "Info_Node error, can't find the Info_NodeEntry, index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		Info_Node& operator=(const Info_Node& n) {
			node = n.node;
			index = n.index;
			name = n.name;
			siteNum = n.siteNum;
			return *this;
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			node.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(node.size());
		};

		void add_nodeEntry(int _index, string _name = "", double _value = 0.0, int _note = Note_None, double _charge = 0) {
			for (auto i = node.begin(); i < node.end(); ++i)
				if (i->index == _index) {
					i->name = _name;
					i->value = _value;
					i->note = _note;
					i->charge = _charge;
					return;
				}
			Info_NodeEntry Ne;
			Ne.index = _index;
			Ne.name = _name;
			Ne.value = _value;
			Ne.note = _note;
			Ne.charge = _charge;
			node.push_back(Ne);
		}
		void erase(Info_NodeEntry& n) {
			for (auto i = node.begin(); i < node.end(); ++i) {
				if (i->index == n.index) node.erase(i);
				return;
			}
			cout << "Info_Node error, don't have aim NodeEntry to erase, index : " << n.index << endl;
			SYS_PROGRAM_STOP;
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<Info_NodeEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Info_NodeEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return node.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return node.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return node.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return node.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<Info_NodeEntry> node;
		int index;
		string name;
		double siteNum;
	};
	class Info_Node2 {
	public:
		Info_Node& operator[](int index) {
			for (auto i = node2.begin(); i < node2.end(); ++i) {
				if (i->index == index) return(*i);
			}
			cout << "Info_Node2 error, can't find the Info_Node, index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		Info_Node2& operator=(const Info_Node2& n) {
			node2 = n.node2;
			return *this;
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			node2.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(node2.size());
		};

		void add_Node(int _index, string _name = "", double _siteNum = 0.0) {
			for (auto i = node2.begin(); i < node2.end(); ++i)
				if (i->index == _index) {
					i->name = _name;
					i->siteNum = _siteNum;
					return;
				}
			Info_Node n;
			n.index = _index;
			n.name = _name;
			n.siteNum = _siteNum;
			node2.push_back(n);
		}
		void erase(Info_Node& n) {
			for (auto i = node2.begin(); i < node2.end(); ++i) {
				if (i->index == n.index) node2.erase(i);
				return;
			}
			cout << "Info_Node2 error, don't have aim Info_Node to erase, index : " << n.index << endl;
			SYS_PROGRAM_STOP;
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<Info_Node>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Info_Node>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return node2.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return node2.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return node2.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return node2.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<Info_Node> node2;
	};
	class Info_Phase {
	public:
		int phase_property;
		int solvent;
		double conductivity;
		double heat_diffusivity;
		string phase_name;
		ThermodynamicModel thermodynamic_model;
		LinearCompoundModel linear_comp_model;
		bool is_energy_treatment;
		Info_energy_treatment energy_treatment;
		Matrix6x6 elastic_constants;
		Info_Node x;
		Info_Node n;
		Info_Node2 sublattice;
		Info_Phase() {
			phase_property = 0;
			is_energy_treatment = false;
			solvent = SOLVENT_NONE;
			conductivity = 0.0;
			heat_diffusivity = 0.0;
			phase_name = "";
			thermodynamic_model = ThermodynamicModel::SolutionModel;
			elastic_constants.set_to_zero();
		}
		Info_Phase& operator=(const Info_Phase& p) {
			phase_property = p.phase_property;
			solvent = p.solvent;
			conductivity = p.conductivity;
			heat_diffusivity = p.heat_diffusivity;
			phase_name = p.phase_name;
			thermodynamic_model = p.thermodynamic_model;
			linear_comp_model = p.linear_comp_model;
			x = p.x;
			n = p.n;
			sublattice = p.sublattice;
			elastic_constants = p.elastic_constants;
			energy_treatment = p.energy_treatment;
			is_energy_treatment = p.is_energy_treatment;
			return *this;
		}
	};
	class Info_Phases {
	public:
		Info_Phases() {
		}
		Info_Phase& operator[](int property) {
			for (auto i = phases.begin(); i < phases.end(); ++i) {
				if (i->phase_property == property) return(*i);
			}
			cout << "Info_Phases error, can't find the Info_Phase, property : " << property << endl;
			SYS_PROGRAM_STOP;
		}
		Info_Phases& operator=(const Info_Phases& n) {
			phases = n.phases;
			return *this;
		}
		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			phases.clear();
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(phases.size());
		};

		void add_Phase(int _property, string _name) {
			for (auto i = phases.begin(); i < phases.end(); ++i)
				if (i->phase_property == _property) {
					i->phase_name = _name;
					return;
				}
			Info_Phase phase;
			phase.phase_property = _property;
			phase.phase_name = _name;
			phases.push_back(phase);
		}
		void add_Phase(Info_Phase phase) {
			phases.push_back(phase);
		}
		void erase(Info_Phase& n) {
			for (auto i = phases.begin(); i < phases.end(); ++i) {
				if (i->phase_property == n.phase_property) phases.erase(i);
				return;
			}
			cout << "Info_Phases error, don't have aim phase to erase, property : " << n.phase_property << endl;
			SYS_PROGRAM_STOP;
		}

		//std::vector<int> interphaseIndex();      

		typedef std::vector<Info_Phase>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<Info_Phase>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return phases.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return phases.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return phases.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return phases.cend(); };                           //< Constant iterator to the end of storage vector

		std::vector<Info_Phase> phases;
	};
	class Set_ElectricFieldMask
	{
	public:
		Set_ElectricFieldMask() {
			x_bc = BoundaryCondition::PERIODIC;
			y_bc = BoundaryCondition::PERIODIC;
			z_bc = BoundaryCondition::PERIODIC;
			electric_potential_x_down = 0.0;
			electric_potential_x_up = 0.0;
			electric_potential_y_down = 0.0;
			electric_potential_y_up = 0.0;
			electric_potential_z_down = 0.0;
			electric_potential_z_up = 0.0;
			init_electric_potential = 0.0;
		};
		~Set_ElectricFieldMask() {
			electric_potential_phase_index.clear();
		};
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		double init_electric_potential;
		double electric_potential_x_down;
		double electric_potential_x_up;
		double electric_potential_y_down;
		double electric_potential_y_up;
		double electric_potential_z_down;
		double electric_potential_z_up;
		double_box electric_potential_phase_index;
		Set_ElectricFieldMask& operator=(const Set_ElectricFieldMask& p) {
			x_bc = p.x_bc;
			y_bc = p.y_bc;
			z_bc = p.z_bc;
			init_electric_potential = p.init_electric_potential;
			electric_potential_x_down = p.electric_potential_x_down;
			electric_potential_x_up = p.electric_potential_x_up;
			electric_potential_y_down = p.electric_potential_y_down;
			electric_potential_y_up = p.electric_potential_y_up;
			electric_potential_z_down = p.electric_potential_z_down;
			electric_potential_z_up = p.electric_potential_z_up;
			electric_potential_phase_index = p.electric_potential_phase_index;
			return *this;
		}
	};
	class Set_MagneticFieldMask
	{
	public:
		Set_MagneticFieldMask() {
			x_bc = BoundaryCondition::PERIODIC;
			y_bc = BoundaryCondition::PERIODIC;
			z_bc = BoundaryCondition::PERIODIC;
			magnetic_potential_x_down = 0.0;
			magnetic_potential_x_up = 0.0;
			magnetic_potential_y_down = 0.0;
			magnetic_potential_y_up = 0.0;
			magnetic_potential_z_down = 0.0;
			magnetic_potential_z_up = 0.0;
			init_magnetic_potential = 0.0;
			is_averaged = false;
			average_value = 0.0;
		};
		~Set_MagneticFieldMask() {
			magnetic_potential_phase_index.clear();
		};
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		double init_magnetic_potential;
		double magnetic_potential_x_down;
		double magnetic_potential_x_up;
		double magnetic_potential_y_down;
		double magnetic_potential_y_up;
		double magnetic_potential_z_down;
		double magnetic_potential_z_up;
		bool is_averaged;
		double average_value;
		double_box magnetic_potential_phase_index;
		Set_MagneticFieldMask& operator=(const Set_MagneticFieldMask& p) {
			is_averaged = p.is_averaged;
			average_value = p.average_value;
			x_bc = p.x_bc;
			y_bc = p.y_bc;
			z_bc = p.z_bc;
			init_magnetic_potential = p.init_magnetic_potential;
			magnetic_potential_x_down = p.magnetic_potential_x_down;
			magnetic_potential_x_up = p.magnetic_potential_x_up;
			magnetic_potential_y_down = p.magnetic_potential_y_down;
			magnetic_potential_y_up = p.magnetic_potential_y_up;
			magnetic_potential_z_down = p.magnetic_potential_z_down;
			magnetic_potential_z_up = p.magnetic_potential_z_up;
			magnetic_potential_phase_index = p.magnetic_potential_phase_index;
			return *this;
		}
	};
	class Set_FluidFieldMask {
	public:
		Set_FluidFieldMask() {
			density = 1.0;
			viscosity = 1.0;
		};
		~Set_FluidFieldMask() {
			rigid_phase_box.clear();
			fluid_phase_box.clear();
		};
		vector<int> rigid_phase_box;
		vector<int> fluid_phase_box;
		double density;
		double viscosity;
		Set_FluidFieldMask& operator=(const Set_FluidFieldMask& p) {
			rigid_phase_box = p.rigid_phase_box;
			fluid_phase_box = p.fluid_phase_box;
			density = p.density;
			viscosity = p.viscosity;
			return *this;
		}
	};
	struct Set_Disperse {
		int Nx;
		int Ny;
		int Nz;
		double int_width;
		int begin_step;
		int end_step;
		double dx;
		double dt;
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		Set_Disperse() {
			Nx = 32;
			Ny = 32;
			Nz = 1;
			begin_step = 0;
			end_step = 0;
			dt = 1.0;
			dx = 1.0;
			int_width = 4.0 * dx;
			x_bc = PERIODIC;
			y_bc = PERIODIC;
			z_bc = PERIODIC;
		}
		Set_Disperse& operator=(const Set_Disperse& n) {
			Nx = n.Nx;
			Ny = n.Ny;
			Nz = n.Nz;
			int_width = n.int_width;
			begin_step = n.begin_step;
			end_step = n.end_step;
			dx = n.dx;
			dt = n.dt;
			x_bc = n.x_bc;
			y_bc = n.y_bc;
			z_bc = n.z_bc;
			return *this;
		}
	};
	struct Set_OutputFile {
		string working_folder_path;
		string restart_datafile_path;
		string log_file_name;
		string data_file_name;
		double NaN;
		int file_output_step;
		int dataFile_output_step;
		int screen_output_step;
		bool is_write_datafile;
		bool is_init_byDatafile;
		bool is_Datafile_init_by_path;
		bool isDrivingForceOutput;
		PairFlag drivingForce_index;
		bool isInterfaceEnergyOutput;
		PairFlag intEnergy_index;
		bool isChemicalPotentialDivisionOutput;
		bool isMechanicalFieldOutput;
		bool isElectricFieldOutput;
		bool isMagneticFieldOutput;
		bool isFluidFieldOutput;
		bool isEnergyDensityOutput;
		bool isPhiOutput;
		bool isPhaseIndexsOutput;
		bool isPhasePropertyOutput;
		bool isPhaseFlagsOutput;
		bool isPhaseGradientOutput;
		bool isPhaseConOutput;
		bool isPhasePotentialOutput;
		bool isTemperatureOutput;
		bool isCustomValueOutput;
		bool isGrainsValueReverse;
		string_box customValue_output;
		bool isCustomFlagOutput;
		string_box customFlag_output;
		bool isCustomVec3Output;
		string_box customVec3_output;
		Set_OutputFile() {
			working_folder_path = "";
			restart_datafile_path = "";
			log_file_name = "log";
			data_file_name = "DATA";
			dataFile_output_step = 1;
			file_output_step = 1;
			screen_output_step = 1;
			NaN = double(-1e99);
			is_init_byDatafile = false;
			is_write_datafile = false;
			isDrivingForceOutput = false;
			drivingForce_index.set(PairValueProperty::pf_QUANTITY);
			isChemicalPotentialDivisionOutput = false;
			intEnergy_index.set(PairValueProperty::pf_QUANTITY);
			isInterfaceEnergyOutput = false;
			isMechanicalFieldOutput = false;
			isElectricFieldOutput = false;
			isMagneticFieldOutput = false;
			isFluidFieldOutput = false;
			isEnergyDensityOutput = false;
			isPhiOutput = false;
			isPhaseIndexsOutput = false;
			isPhasePropertyOutput = false;
			isPhaseFlagsOutput = false;
			isPhaseConOutput = false;
			isPhasePotentialOutput = false;
			isPhaseGradientOutput = false;
			isTemperatureOutput = false;
			isGrainsValueReverse = false;
			isCustomValueOutput = false;
			customValue_output.clear();
			isCustomFlagOutput = false;
			customFlag_output.clear();
			isCustomVec3Output = false;
			customVec3_output.clear();
#ifdef _WIN32
			is_Datafile_init_by_path = false;
#else
			is_Datafile_init_by_path = true;
#endif
		}
		Set_OutputFile& operator=(const Set_OutputFile& n) {
			working_folder_path = n.working_folder_path;
			restart_datafile_path = n.restart_datafile_path;
			log_file_name = n.log_file_name;
			data_file_name = n.data_file_name;
			NaN = n.NaN;
			file_output_step = n.file_output_step;
			dataFile_output_step = n.dataFile_output_step;
			screen_output_step = n.screen_output_step;
			is_init_byDatafile = n.is_init_byDatafile;
			is_write_datafile = n.is_write_datafile;
			isDrivingForceOutput = n.isDrivingForceOutput;
			drivingForce_index = n.drivingForce_index;
			isInterfaceEnergyOutput = n.isInterfaceEnergyOutput;
			intEnergy_index = n.intEnergy_index;
			isChemicalPotentialDivisionOutput = n.isChemicalPotentialDivisionOutput;
			isMechanicalFieldOutput = n.isMechanicalFieldOutput;
			isElectricFieldOutput = n.isElectricFieldOutput;
			isMagneticFieldOutput = n.isMagneticFieldOutput;
			isFluidFieldOutput = n.isFluidFieldOutput;
			isEnergyDensityOutput = n.isEnergyDensityOutput;
			is_Datafile_init_by_path = n.is_Datafile_init_by_path;
			isPhiOutput = n.isPhiOutput;
			isPhaseGradientOutput = n.isPhaseGradientOutput;
			isPhaseIndexsOutput = n.isPhaseIndexsOutput;
			isPhasePropertyOutput = n.isPhasePropertyOutput;
			isPhaseFlagsOutput = n.isPhaseFlagsOutput;
			isPhaseConOutput = n.isPhaseConOutput;
			isPhasePotentialOutput = n.isPhasePotentialOutput;
			isTemperatureOutput = n.isTemperatureOutput;
			isGrainsValueReverse = n.isGrainsValueReverse;
			isCustomValueOutput = n.isCustomValueOutput;
			isCustomFlagOutput = n.isCustomFlagOutput;
			isCustomVec3Output = n.isCustomVec3Output;
			customValue_output = n.customValue_output;
			customFlag_output = n.customFlag_output;
			customVec3_output = n.customVec3_output;
			return *this;
		}
	};
}