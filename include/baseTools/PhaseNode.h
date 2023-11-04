#pragma once
#include"NormalNode.h"
#include"PairValueNode.h"
#include"SublatticNode.h"

using namespace std;
namespace pf {
	enum Flag { pf_BULK, pf_NEAR_INTERFACE, pf_INTERFACE }; //check in calculation
	enum CFlag { pf_NoCal, pf_COMP_NoAverage, pf_COMP_Average, pf_COMP_NewPoint, pf_COMP_DelPoint };  //will be used by concentration increment average
	struct elecNode {
		elecNode() {
			set_zero();
		}
		void set_zero() {
			elec_potential = 0.0;
			laplace_elec_potential = 0.0;
			elec_potential_grad.set_to_zero();
			effective_conductivity = 0.0;
			elec_charge_density = 0.0;
		}
		elecNode& operator=(const elecNode& n) {
			elec_potential = n.elec_potential;
			laplace_elec_potential = n.laplace_elec_potential;
			elec_potential_grad = n.elec_potential_grad;
			effective_conductivity = n.effective_conductivity;
			elec_charge_density = n.elec_charge_density;
			return *this;
		}
		double elec_potential;
		double laplace_elec_potential;
		Vector3 elec_potential_grad;
		double effective_conductivity;
		double elec_charge_density;
	};
	struct magNode {
		magNode() {
			set_zero();
		}
		void set_zero() {
			mag_potential = 0.0;
			mag_charge_density = 0.0;
			mag_intensity[0] = 0.0;
			mag_intensity[1] = 0.0;
			mag_intensity[2] = 0.0;
		}
		magNode& operator=(const magNode& n) {
			mag_potential = n.mag_potential;
			mag_charge_density = n.mag_charge_density;
			mag_intensity[0] = n.mag_intensity[0];
			mag_intensity[1] = n.mag_intensity[1];
			mag_intensity[2] = n.mag_intensity[2];
			dmag_potential_dphi = n.dmag_potential_dphi;
			dmag_intensity_dphi = n.dmag_intensity_dphi;
			dmag_charge_density_dphi = n.dmag_charge_density_dphi;
			return *this;
		}
		double mag_potential;
		Vector3 mag_intensity;
		double mag_charge_density;
		double_box dmag_potential_dphi;
		vec3_box dmag_intensity_dphi;
		double_box dmag_charge_density_dphi;
	};
	struct mechNode {
		mechNode() {
			set_zero();
		}
		void set_zero() {
			Stresses.set_to_zero();
			Strains.set_to_zero();
			EffectiveEigenStrain.set_to_zero();
		}
		mechNode& operator=(const mechNode& n) {
			Stresses = n.Stresses;
			Strains = n.Strains;
			EffectiveEigenStrain = n.EffectiveEigenStrain;
			return *this;
		}
		vStress       Stresses;                                       ///< Storage for stresses
		vStrain       Strains;                                        ///< Storage for strains
		vStrain		  EffectiveEigenStrain;
	};
	struct tempNode {
		tempNode() {
			set_zero();
		}
		void set_zero() {
			temperature = 0;
			heat_diffusivity = 0;
			temp_increment = 0;
			temperature_grad.set_to_zero();
			diffusivity_grad.set_to_zero();
			laplace = 0;
		}
		tempNode& operator=(const tempNode& n) {
			temperature = n.temperature;
			heat_diffusivity = n.heat_diffusivity;
			temp_increment = n.temp_increment;
			temperature_grad = n.temperature_grad;
			diffusivity_grad = n.diffusivity_grad;
			laplace = n.laplace;
			return *this;
		}
		double temperature;
		Vector3 temperature_grad;
		Vector3 diffusivity_grad;
		double heat_diffusivity;
		double temp_increment;
		double laplace;
	};
	struct velocityNode {
		velocityNode() {
			set_zero();
		}
		void set_zero() {
			//density = 0.0;
			pressure = 0.0;
			velocity.set_to_zero();
			volume_force.set_to_zero();
		}
		velocityNode& operator=(const velocityNode& n) {
			//density = n.density;
			pressure = n.pressure;
			velocity = n.velocity;
			volume_force = n.volume_force;
			return *this;
		}
		//double density;
		double pressure;
		Vector3 velocity;
		Vector3 volume_force;
	};

	class PhaseEntry
	{
	public:
		PhaseEntry() {
			_flag = 0;
			_cflag = 0;
			index = 0;
			phaseProperty = 0;
			phaseFraction = 0.0;
			int_increment = 0.0;
			bulk_increment = 0.0;
			chemEnergyDensity = 0.0;
			elasEnergyDensity = 0.0;
			elecEnergyDensity = 0.0;
			magEnergyDensity = 0.0;
			otherEnergyDensity = 0.0;
			laplacian = 0.0;
			phi_grad.set_to_zero();
		};
		~PhaseEntry() {
			x.clear();
			n.clear();
			potential.clear();
			sublattice.clear();
			kinetics_coeff.clear();
		};
		
		PhaseEntry& operator=(const PhaseEntry& n);
		void set(int _index, int _property, int flag = 0, double _phaseFraction = 0.0) {
			index = _index;
			phaseProperty = _property;
			_flag = flag;
			phaseFraction = _phaseFraction;
		}
		//< flags
		int _flag; //Node_Phase_Flag
		int _cflag; //phase concentration average step flag
		int index; //search for this entry
		int phaseProperty;
		//< flags
		//< data
		double phaseFraction;  // main variable
		double laplacian;
		Vector3 phi_grad;
		double int_increment;
		double bulk_increment;
		double chemEnergyDensity;
		double elasEnergyDensity;
		double elecEnergyDensity;
		double magEnergyDensity;
		double otherEnergyDensity;

		ConNode x;  // main variable
		XNode n;  //used in  associate liquid model
		ChemNode potential;
		SublatticeNode sublattice;
		PairKinetic kinetics_coeff;

		//< data
	};
	inline PhaseEntry& PhaseEntry::operator=(const PhaseEntry& node) {
		_flag = node._flag;
		_cflag = node._cflag;
		index = node.index;
		phaseProperty = node.phaseProperty;

		phaseFraction = node.phaseFraction;
		laplacian = node.laplacian;
		phi_grad = node.phi_grad;
		int_increment = node.int_increment;
		bulk_increment = node.bulk_increment;
		chemEnergyDensity = node.chemEnergyDensity;
		elasEnergyDensity = node.elasEnergyDensity;
		elecEnergyDensity = node.elecEnergyDensity;
		magEnergyDensity = node.magEnergyDensity;
		otherEnergyDensity = node.otherEnergyDensity;

		x = node.x;
		n = node.n;
		potential = node.potential;
		sublattice = node.sublattice;
		kinetics_coeff = node.kinetics_coeff;

		return *this;
	}

	class PhaseNode
	{
	public:
		PhaseNode() {
			_Phase.reserve(2);
		};
		~PhaseNode() {
			clear();
		};
		PhaseEntry& operator[](const int index);                                       //< Index operator for accessing the n's field value
		PhaseNode&  operator=(const PhaseNode& n);                                            //< Assignement operator
		typedef std::vector<PhaseEntry>::iterator iterator;                         //< Iterator over storage vector
		typedef std::vector<PhaseEntry>::const_iterator citerator;                  //< Constant iterator over storage vector
		iterator  begin() { return _Phase.begin(); };                                 //< Iterator to the begin of storage vector
		iterator  end() { return _Phase.end(); };                                   //< Iterator to the end of storage vector
		citerator cbegin() const { return _Phase.cbegin(); };                         //< Constant iterator to the begin of storage vector
		citerator cend()   const { return _Phase.cend(); };                           //< Constant iterator to the end of storage vector
		void add_phase(int _index, int _property = 0, int flag = 0, double _phasefraction = 0.0) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i)
				if (i->index == _index) {
					i->_flag = flag;
					i->phaseProperty = _property;
					i->phaseFraction = _phasefraction;
					return;
				}
			PhaseEntry entry;
			entry.index = _index;
			entry._flag = flag;
			entry.phaseProperty = _property;
			entry.phaseFraction = _phasefraction;
			_Phase.push_back(entry);
		}
		void add_phase(PhaseEntry _sample) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i)
				if (i->index == _sample.index) {
					i->phaseProperty = _sample.phaseProperty;
					i->phaseFraction = _sample.phaseFraction;
					i->x = _sample.x;
					i->sublattice = _sample.sublattice;
					i->n = _sample.n;
					i->_flag = _sample._flag;
					i->potential = _sample.potential;
					i->chemEnergyDensity = _sample.chemEnergyDensity;
					i->elasEnergyDensity = _sample.elasEnergyDensity;
					i->kinetics_coeff = _sample.kinetics_coeff;
					i->_cflag = _sample._cflag;
					return;
				}
			PhaseEntry entry;
			entry.index = _sample.index;
			entry.phaseProperty = _sample.phaseProperty;
			entry.phaseFraction = _sample.phaseFraction;
			entry.x = _sample.x;
			entry.sublattice = _sample.sublattice;
			entry.n = _sample.n;
			entry._flag = _sample._flag;
			entry.potential = _sample.potential;
			entry.chemEnergyDensity = _sample.chemEnergyDensity;
			entry.elasEnergyDensity = _sample.elasEnergyDensity;
			entry.kinetics_coeff = _sample.kinetics_coeff;
			entry._cflag = _sample._cflag;
			_Phase.push_back(entry);
		}
		bool is_phase_exist(int phaseIndex) {
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				if (phase->index == phaseIndex && phase->phaseFraction > Simulation_Num_Cut_Off)
					return true;
			return false;
		}
		void cal_x_from_phase_x() {
			for (auto con = x.begin(); con < x.end(); con++)
				con->value = 0.0;
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				for (auto con = phase->x.begin(); con < phase->x.end(); con++)
					x[con->index].value += phase->phaseFraction * con->value;
		}
		void cal_potential_from_phase_potential() {
			for (auto p = potential.begin(); p < potential.end(); p++) {
				p->value = 0.0;
				p->chemical_part = 0.0;
				p->elastic_part = 0.0;
				p->electric_part = 0.0;
				p->magnetic_part = 0.0;
				p->other_part = 0.0;
			}
			for (auto phase = (*this).begin(); phase < (*this).end(); phase++)
				for (auto p = phase->potential.begin(); p < phase->potential.end(); p++) {
					potential[p->index].value += phase->phaseFraction * p->value;
					potential[p->index].chemical_part += phase->phaseFraction * p->chemical_part;
					potential[p->index].elastic_part += phase->phaseFraction * p->elastic_part;
					potential[p->index].electric_part += phase->phaseFraction * p->electric_part;
					potential[p->index].magnetic_part += phase->phaseFraction * p->magnetic_part;
					potential[p->index].other_part += phase->phaseFraction * p->other_part;
				}
		}
		double get_ave_phaseFraction(int phaseIndex, int opt_range = 1) {
			int num = 0;
			double p_f = 0.0;
			for (int i = -opt_range; i <= +opt_range; i++)
				for (int j = -opt_range; j <= +opt_range; j++)
					for (int k = -opt_range; k <= +opt_range; k++) {
						p_f += get_long_range_node(i, j, k)[phaseIndex].phaseFraction;
						num += 1;
					}
			if (num > 0)
				p_f = p_f / num;
			else
				p_f = this->operator[](phaseIndex).phaseFraction;
			return p_f;
		}
		void normalized_phaseFraction() {
			double sum = 0.0;
			for (auto phase = _Phase.begin(); phase < _Phase.end(); phase++)
				sum += phase->phaseFraction;
			for (auto phase = _Phase.begin(); phase < _Phase.end(); phase++)
				phase->phaseFraction = phase->phaseFraction / sum;
		}
		double cal_customValues_laplace(int customValues_index, double dx = 1.0, DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT) {
			double laplace = 0.0;
			if (diff_method == DifferenceMethod::FIVE_POINT) {
				laplace = (get_neighbor_node(Direction::x_down).customValues[customValues_index]
					+ get_neighbor_node(Direction::x_up).customValues[customValues_index]
					+ get_neighbor_node(Direction::y_down).customValues[customValues_index]
					+ get_neighbor_node(Direction::y_up).customValues[customValues_index]
					+ get_neighbor_node(Direction::z_down).customValues[customValues_index]
					+ get_neighbor_node(Direction::z_up).customValues[customValues_index]
					- 6.0 * customValues[customValues_index]) / dx / dx;
			}
			else if (diff_method == DifferenceMethod::NINE_POINT) {
				laplace = (4.0 * get_neighbor_node(Direction::x_down).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::x_up).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::y_down).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::y_up).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::z_down).customValues[customValues_index]
					+ 4.0 * get_neighbor_node(Direction::z_up).customValues[customValues_index]
					+ get_long_range_node(-1, -1, 0).customValues[customValues_index]
					+ get_long_range_node(-1, 1, 0).customValues[customValues_index]
					+ get_long_range_node(1, -1, 0).customValues[customValues_index]
					+ get_long_range_node(1, 1, 0).customValues[customValues_index]
					+ get_long_range_node(-1, 0, -1).customValues[customValues_index]
					+ get_long_range_node(-1, 0, 1).customValues[customValues_index]
					+ get_long_range_node(1, 0, -1).customValues[customValues_index]
					+ get_long_range_node(1, 0, 1).customValues[customValues_index]
					+ get_long_range_node(0, -1, -1).customValues[customValues_index]
					+ get_long_range_node(0, -1, 1).customValues[customValues_index]
					+ get_long_range_node(0, 1, -1).customValues[customValues_index]
					+ get_long_range_node(0, 1, 1).customValues[customValues_index]
					- 36.0 * customValues[customValues_index]) / 6.0 / dx / dx;
			}
			else {
				cout << "Difference method define error!";
				SYS_PROGRAM_STOP;
			}
			return laplace;
		}
		Vector3 cal_customValues_gradient(int customValues_index, double dx = 1.0) {
			Vector3 grad;
			grad[0] = (get_neighbor_node(Direction::x_down).customValues[customValues_index] - get_neighbor_node(Direction::x_up).customValues[customValues_index]) / 2.0 / dx;
			grad[1] = (get_neighbor_node(Direction::y_down).customValues[customValues_index] - get_neighbor_node(Direction::y_up).customValues[customValues_index]) / 2.0 / dx;
			grad[2] = (get_neighbor_node(Direction::z_down).customValues[customValues_index] - get_neighbor_node(Direction::z_up).customValues[customValues_index]) / 2.0 / dx;
			return grad;
		}
		double cal_dfchem_dphi(PhaseEntry& alpha, PhaseEntry& beta) {
			double df_dphi = beta.chemEnergyDensity - alpha.chemEnergyDensity;
			for (auto cb = beta.x.begin(); cb < beta.x.end(); cb++)
				for (auto ca = alpha.x.begin(); ca < alpha.x.end(); ca++)
					if(cb->index == ca->index)
						df_dphi -= (cb->value - ca->value) * potential[ca->index].chemical_part;
#ifdef _DEBUG
			if (_isnan(df_dphi)) {
				cout << "DEBUG: df_dphi error !" << endl;
				SYS_PROGRAM_STOP;
			}
#endif
			return df_dphi;
		}
		double cal_dfelas_dphi(PhaseEntry& alpha, PhaseEntry& beta) {
			double df_dphi = beta.elasEnergyDensity - alpha.elasEnergyDensity;
			for (auto cb = beta.x.begin(); cb < beta.x.end(); cb++)
				for (auto ca = alpha.x.begin(); ca < alpha.x.end(); ca++)
					if (cb->index == ca->index)
						df_dphi -= (cb->value - ca->value) * potential[ca->index].elastic_part;
#ifdef _DEBUG
			if (_isnan(df_dphi)) {
				cout << "DEBUG: df_dphi error !" << endl;
				SYS_PROGRAM_STOP;
			}
#endif
			return df_dphi;
		}
		double cal_dfelec_dphi(PhaseEntry& alpha, PhaseEntry& beta) {
			double df_dphi = beta.elecEnergyDensity - alpha.elecEnergyDensity;
			for (auto cb = beta.x.begin(); cb < beta.x.end(); cb++)
				for (auto ca = alpha.x.begin(); ca < alpha.x.end(); ca++)
					if (cb->index == ca->index)
						df_dphi -= (cb->value - ca->value) * potential[ca->index].electric_part;
#ifdef _DEBUG
			if (_isnan(df_dphi)) {
				cout << "DEBUG: df_dphi error !" << endl;
				SYS_PROGRAM_STOP;
			}
#endif
			return df_dphi;
		}
		double cal_dfmag_dphi(PhaseEntry& alpha, PhaseEntry& beta) {
			double df_dphi = beta.magEnergyDensity - alpha.magEnergyDensity;
			for (auto cb = beta.x.begin(); cb < beta.x.end(); cb++)
				for (auto ca = alpha.x.begin(); ca < alpha.x.end(); ca++)
					if (cb->index == ca->index)
						df_dphi -= (cb->value - ca->value) * potential[ca->index].magnetic_part;
#ifdef _DEBUG
			if (_isnan(df_dphi)) {
				cout << "DEBUG: df_dphi error !" << endl;
				SYS_PROGRAM_STOP;
			}
#endif
			return df_dphi;
		}
		double cal_dfother_dphi(PhaseEntry& alpha, PhaseEntry& beta) {
			double df_dphi = beta.otherEnergyDensity - alpha.otherEnergyDensity;
			for (auto cb = beta.x.begin(); cb < beta.x.end(); cb++)
				for (auto ca = alpha.x.begin(); ca < alpha.x.end(); ca++)
					if (cb->index == ca->index)
						df_dphi -= (cb->value - ca->value) * potential[ca->index].other_part;
#ifdef _DEBUG
			if (_isnan(df_dphi)) {
				cout << "DEBUG: df_dphi error !" << endl;
				SYS_PROGRAM_STOP;
			}
#endif
			return df_dphi;
		}
		double cal_dftot_dphi(PhaseEntry& alpha, PhaseEntry& beta) {
			double df_dphi = beta.chemEnergyDensity + beta.elasEnergyDensity + beta.elecEnergyDensity + beta.magEnergyDensity + beta.otherEnergyDensity 
				- alpha.chemEnergyDensity - alpha.elasEnergyDensity - alpha.elecEnergyDensity - alpha.magEnergyDensity - alpha.otherEnergyDensity;
			for (auto cb = beta.x.begin(); cb < beta.x.end(); cb++)
				for (auto ca = alpha.x.begin(); ca < alpha.x.end(); ca++)
					if (cb->index == ca->index)
						df_dphi -= (cb->value - ca->value) * potential[ca->index].value;
#ifdef _DEBUG
			if (_isnan(df_dphi)) {
				cout << "DEBUG: df_dphi error !" << endl;
				SYS_PROGRAM_STOP;
			}
#endif
			return df_dphi;
		}
		void automatic_set_flag();
		//< Data
		std::vector<PhaseEntry> _Phase;
		ConNode x;
		ChemNode potential;

		elecNode electricValues;
		magNode magneticValues;
		mechNode mechanicalValues;
		tempNode tempValues;
		velocityNode velocityValues;

		double_box customValues;
		int_box customFlags;
		vec3_box customVec3s;
		//< Data

		void clear()                                                              //< Emptys the field storage. Sets flag to 0.
		{
			_Phase.clear(); 
			customValues.clear();
			customVec3s.clear();
			customFlags.clear();
			_up_x = nullptr;
			_down_x = nullptr;
			_up_y = nullptr;
			_down_y = nullptr;
			_up_z = nullptr;
			_down_z = nullptr;
		};
		int size() const                                                   //< Returns the size of storage.
		{
			return int(_Phase.size());
		};
		int phaseNum();
		void push_back(PhaseEntry& n) {
			_Phase.push_back(n);
		}
		void erase(PhaseEntry& n) {
			for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
				if (i->index == n.index) {
					_Phase.erase(i);
					return;
				}
			}
			cout << "PhaseNode error, don't have aim PhaseEntry to erase, index : " << n.index << endl;
			SYS_PROGRAM_STOP;
		}
		void connect_mesh(int x, int y, int z, PhaseNode& up_x, PhaseNode& down_x, PhaseNode& up_y, PhaseNode& down_y, PhaseNode& up_z, PhaseNode& down_z) {
			_x = x;
			_y = y;
			_z = z;
			_up_x = &up_x;
			_down_x = &down_x;
			_up_y = &up_y;
			_down_y = &down_y;
			_up_z = &up_z;
			_down_z = &down_z;
		}
		PhaseNode& get_neighbor_node(Direction _d);
		PhaseNode& get_long_range_node(int relative_x, int relative_y, int relative_z);
		int _x; // 0 - (Nx-1)
		int _y;
		int _z;
	private:
		PhaseNode* _up_x;
		PhaseNode* _down_x;
		PhaseNode* _up_y;
		PhaseNode* _down_y;
		PhaseNode* _up_z;
		PhaseNode* _down_z;
	};
	inline PhaseEntry& PhaseNode::operator[](const int index) {
		for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
			if (i->index == index) return(*i);
		}
		cout << "PhaseNode error, can't find the PhaseEntry, index : " << index << endl;
		SYS_PROGRAM_STOP;
	}
	inline int PhaseNode::phaseNum() {
		int result = 0;
		for (auto i = _Phase.begin(); i < _Phase.end(); ++i) {
			if (i->_flag != pf_BULK) result++;
		}
		if (result == 0)  ///< both pf_bulk, in the phase's bulk
			result = 1;
		return result;
	}
	inline PhaseNode& PhaseNode::operator=(const PhaseNode& n) {
		_Phase = n._Phase;
		x = n.x;
		potential = n.potential;
		electricValues = n.electricValues;
		magneticValues = n.magneticValues;
		mechanicalValues = n.mechanicalValues;
		tempValues = n.tempValues;
		velocityValues = n.velocityValues;
		customValues = n.customValues;
		customFlags = n.customFlags;
		customVec3s = n.customVec3s;
		return *this;
	}
	inline PhaseNode& PhaseNode::get_neighbor_node(Direction _d) {
		switch (_d)
		{
		case pf::x_up:
			return *_up_x;
			break;
		case pf::x_down:
			return *_down_x;
			break;
		case pf::y_up:
			return *_up_y;
			break;
		case pf::y_down:
			return *_down_y;
			break;
		case pf::z_up:
			return *_up_z;
			break;
		case pf::z_down:
			return *_down_z;
			break;
		default:
			cout << "Class get_neighbor_node parameter error !" << endl;
			SYS_PROGRAM_STOP;
			break;
		}
	}
	inline PhaseNode& PhaseNode::get_long_range_node(int relative_x, int relative_y, int relative_z) {
		if (relative_x > 0)
			return (*this).get_neighbor_node(Direction::x_up).get_long_range_node(relative_x - 1, relative_y, relative_z);
		else if (relative_x < 0)
			return (*this).get_neighbor_node(Direction::x_down).get_long_range_node(relative_x + 1, relative_y, relative_z);
		if (relative_y > 0)
			return (*this).get_neighbor_node(Direction::y_up).get_long_range_node(relative_x, relative_y - 1, relative_z);
		else if (relative_y < 0)
			return (*this).get_neighbor_node(Direction::y_down).get_long_range_node(relative_x, relative_y + 1, relative_z);
		if (relative_z > 0)
			return (*this).get_neighbor_node(Direction::z_up).get_long_range_node(relative_x, relative_y, relative_z - 1);
		else if (relative_z < 0)
			return (*this).get_neighbor_node(Direction::z_down).get_long_range_node(relative_x, relative_y, relative_z + 1);
		return (*this);
	}
	inline void PhaseNode::automatic_set_flag() {
		for (auto phase = _Phase.begin(); phase < _Phase.end(); phase++) {
			phase->_flag = pf_BULK;
			if (phase->phaseFraction >= Simulation_Num_Cut_Off && phase->phaseFraction <= (1.0 - Simulation_Num_Cut_Off)) {
				phase->_flag = pf_INTERFACE;
				if (get_neighbor_node(Direction::x_up)[phase->index].phaseFraction > (1.0 - Simulation_Num_Cut_Off) || get_neighbor_node(Direction::x_up)[phase->index].phaseFraction < Simulation_Num_Cut_Off)
					get_neighbor_node(Direction::x_up)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::x_down)[phase->index].phaseFraction > (1.0 - Simulation_Num_Cut_Off) || get_neighbor_node(Direction::x_down)[phase->index].phaseFraction < Simulation_Num_Cut_Off)
					get_neighbor_node(Direction::x_down)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::y_up)[phase->index].phaseFraction > (1.0 - Simulation_Num_Cut_Off) || get_neighbor_node(Direction::y_up)[phase->index].phaseFraction < Simulation_Num_Cut_Off)
					get_neighbor_node(Direction::y_up)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::y_down)[phase->index].phaseFraction > (1.0 - Simulation_Num_Cut_Off) || get_neighbor_node(Direction::y_down)[phase->index].phaseFraction < Simulation_Num_Cut_Off)
					get_neighbor_node(Direction::y_down)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::z_up)[phase->index].phaseFraction > (1.0 - Simulation_Num_Cut_Off) || get_neighbor_node(Direction::z_up)[phase->index].phaseFraction < Simulation_Num_Cut_Off)
					get_neighbor_node(Direction::z_up)[phase->index]._flag = pf_NEAR_INTERFACE;
				if (get_neighbor_node(Direction::z_down)[phase->index].phaseFraction > (1.0 - Simulation_Num_Cut_Off) || get_neighbor_node(Direction::z_down)[phase->index].phaseFraction < Simulation_Num_Cut_Off)
					get_neighbor_node(Direction::z_down)[phase->index]._flag = pf_NEAR_INTERFACE;
			}
			if (phase->phaseFraction < Simulation_Num_Cut_Off) {
				if (get_neighbor_node(Direction::x_up)[phase->index].phaseFraction >= Simulation_Num_Cut_Off
					|| get_neighbor_node(Direction::x_down)[phase->index].phaseFraction >= Simulation_Num_Cut_Off
					|| get_neighbor_node(Direction::y_up)[phase->index].phaseFraction >= Simulation_Num_Cut_Off
					|| get_neighbor_node(Direction::y_down)[phase->index].phaseFraction >= Simulation_Num_Cut_Off
					|| get_neighbor_node(Direction::z_down)[phase->index].phaseFraction >= Simulation_Num_Cut_Off
					|| get_neighbor_node(Direction::z_up)[phase->index].phaseFraction >= Simulation_Num_Cut_Off) {
					phase->_flag = pf_NEAR_INTERFACE;
				}
			}
			else if (phase->phaseFraction > (1.0 - Simulation_Num_Cut_Off)) {
				if (get_neighbor_node(Direction::x_up)[phase->index].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
					|| get_neighbor_node(Direction::x_down)[phase->index].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
					|| get_neighbor_node(Direction::y_up)[phase->index].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
					|| get_neighbor_node(Direction::y_down)[phase->index].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
					|| get_neighbor_node(Direction::z_up)[phase->index].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
					|| get_neighbor_node(Direction::z_down)[phase->index].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)) {
					phase->_flag = pf_NEAR_INTERFACE;
				}
			}
		}
	}

}