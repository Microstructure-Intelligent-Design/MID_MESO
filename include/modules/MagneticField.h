#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"
namespace pf {

	class MagneticField
	{
	public:
		MagneticField() {};
		~MagneticField() {};

		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information);
		void clear();
		void setMagneticParameters();
		int solve_magnetizationIntensity(double PotentialAccuracy, int MAXIterations, bool debug_solver = false, int output_step = 1000);
		int solve_dMagnetizationIntensity_dPhi(double PotentialAccuracy, int MAXIterations, bool debug_solver = false, int output_step = 1000);
		void cal_magnetic_increment_of_oneNode(PhaseNode& node);
		/*void cal_MagneticEnergyData(PhaseNode& node) {
			double vec_mag_potential[3];
			vec_mag_potential[0] = (node.get_neighbor_node(Direction::x_down).magneticValue.mag_potential - node.get_neighbor_node(Direction::x_up).magneticValue.mag_potential) / 2.0 / simulationField->dx;
			vec_mag_potential[1] = (node.get_neighbor_node(Direction::y_down).magneticValue.mag_potential - node.get_neighbor_node(Direction::y_up).magneticValue.mag_potential) / 2.0 / simulationField->dx;
			vec_mag_potential[2] = (node.get_neighbor_node(Direction::z_down).magneticValue.mag_potential - node.get_neighbor_node(Direction::z_up).magneticValue.mag_potential) / 2.0 / simulationField->dx;
			node.magEnergyDensity = 1.0 / 2.0 * VacuumPermeability_Magnetismus * (vec_mag_potential[0] * node.magneticValue.mag_intensity[0] + vec_mag_potential[1] * node.magneticValue.mag_intensity[1] + vec_mag_potential[2] * node.magneticValue.mag_intensity[2]);
			for (auto alpha = node.begin(); alpha < node.end() - 1; alpha++)
				for (auto beta = alpha + 1; beta < node.end(); beta++)
					if (alpha->_flag && beta->_flag) {
						double buff_vec[3];
						buff_vec[0] = (node.get_neighbor_node(Direction::x_down).magneticValue.dmag_potential_dphi[alpha->index] - node.get_neighbor_node(Direction::x_up).magneticValue.dmag_potential_dphi[alpha->index]) / 2.0 / simulationField->dx;
						buff_vec[1] = (node.get_neighbor_node(Direction::y_down).magneticValue.dmag_potential_dphi[alpha->index] - node.get_neighbor_node(Direction::y_up).magneticValue.dmag_potential_dphi[alpha->index]) / 2.0 / simulationField->dx;
						buff_vec[2] = (node.get_neighbor_node(Direction::z_down).magneticValue.dmag_potential_dphi[alpha->index] - node.get_neighbor_node(Direction::z_up).magneticValue.dmag_potential_dphi[alpha->index]) / 2.0 / simulationField->dx;
						double df_dphiA = 1.0 / 2.0 * VacuumPermeability_Magnetismus *
							(vec_mag_potential[0] * node.magneticValue.dmag_intensity_dphi(alpha->index, 0)
								+ vec_mag_potential[1] * node.magneticValue.dmag_intensity_dphi(alpha->index, 1)
								+ vec_mag_potential[2] * node.magneticValue.dmag_intensity_dphi(alpha->index, 2)
								+ node.magneticValue.mag_intensity[0] * buff_vec[0]
								+ node.magneticValue.mag_intensity[1] * buff_vec[1]
								+ node.magneticValue.mag_intensity[2] * buff_vec[2]);
						buff_vec[0] = (node.get_neighbor_node(Direction::x_down).magneticValue.dmag_potential_dphi[beta->index] - node.get_neighbor_node(Direction::x_up).magneticValue.dmag_potential_dphi[beta->index]) / 2.0 / simulationField->dx;
						buff_vec[1] = (node.get_neighbor_node(Direction::y_down).magneticValue.dmag_potential_dphi[beta->index] - node.get_neighbor_node(Direction::y_up).magneticValue.dmag_potential_dphi[beta->index]) / 2.0 / simulationField->dx;
						buff_vec[2] = (node.get_neighbor_node(Direction::z_down).magneticValue.dmag_potential_dphi[beta->index] - node.get_neighbor_node(Direction::z_up).magneticValue.dmag_potential_dphi[beta->index]) / 2.0 / simulationField->dx;
						double df_dphiB = 1.0 / 2.0 * VacuumPermeability_Magnetismus *
							(vec_mag_potential[0] * node.magneticValue.dmag_intensity_dphi(beta->index, 0)
								+ vec_mag_potential[1] * node.magneticValue.dmag_intensity_dphi(beta->index, 1)
								+ vec_mag_potential[2] * node.magneticValue.dmag_intensity_dphi(beta->index, 2)
								+ node.magneticValue.mag_intensity[0] * buff_vec[0]
								+ node.magneticValue.mag_intensity[1] * buff_vec[1]
								+ node.magneticValue.mag_intensity[2] * buff_vec[2]);
						node.magDF.set(alpha->index, beta->index, df_dphiB - df_dphiA);

#ifdef _DEBUG
						if (_isnan(node.magDF(alpha->index, beta->index))) {
							cout << "DEBUG: node.magDF error !" << endl;
							SYS_PROGRAM_STOP;
						}
#endif
					}
		}*/

		FieldStorage_forPhaseNode* simulationField;
		Information* information;


	};

}