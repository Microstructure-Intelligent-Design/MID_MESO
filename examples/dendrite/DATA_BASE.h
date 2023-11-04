#pragma once
#include "../../include/baseTools/baseTools.h"
// System for Li-Cu-Sn
namespace DATABASE {
	using namespace pf::materials;
	//special elememt defined lager than ElementTop
	enum Element { Bottom = pf::materials::Element::ElementTop, LI4SN };

	static void printerror(){
		cout << "Temperature, out of limit, thermodynamic calculation error!" << endl;
		SYS_PROGRAM_STOP;
	}

	static void printerror(double T) {
		cout << "Temperature, out of limit, thermodynamic calculation error!" << endl;
		SYS_PROGRAM_STOP;
	}

	const double R = 8.3144521;

	namespace SYS_Solidification {
		enum PHASES { Solid };
		enum COMPONENT{ NONE = pf::materials::ElementTop, CON};

		static pf::Info_Phase Get_Phase_Structure(PHASES phase_property) {
			pf::Info_Phase phase;
			pf::Matrix6x6 C;
			switch (phase_property)
			{
			case Solid:
				phase.phase_property = Solid;
				phase.phase_name = "Solid";
				phase.heat_diffusivity = 1.0;
				break;
			case VIRTUAL_PHASE:
				phase.phase_property = VIRTUAL_PHASE;
				phase.phase_name = "Liquid";
				phase.heat_diffusivity = 1.0;
				break;
			default:
				break;
			}
			return phase;
		}
		static double MolarVolume(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			return 1.0;
		}
		static vector<double> epsilon(pf::PhaseNode& n, int pIndex) {
			double epsilon = 0.0, depsilon = 0.0, epsilonb = 0.01, delta = 0.02, aniso = 6.0, theta0 = 0.2;
			double theta = 0.0;
			if(!Is_Equality(n[pIndex].phi_grad[0], 0.0))
				theta = atan(n[pIndex].phi_grad[1] / n[pIndex].phi_grad[0]);
			epsilon = epsilonb * (1.0 + delta * cos(aniso * (theta - theta0)));
			depsilon = -epsilonb * aniso * delta * sin(aniso * (theta - theta0));
			vector<double> vec;
			vec.push_back(epsilon);
			vec.push_back(depsilon);
			return vec;
		}
		static double dfint_dphi(pf::PhaseNode& n, pf::PhaseEntry& alpha, double iwidth, pf::Info_DynamicCollection& inf) {
			double dfdPhiA = 0.0, lap_phi = 0.0, Alpha = 0.9, gamma = 10.0, teq = 1.0, int_width = 0.03 * 4;
			double m = Alpha / PI * atan(gamma * (teq - n.tempValues.temperature));
			double vec_buff[] = {0.0, 0.0, 0.0};
			if (alpha.phaseProperty != VIRTUAL_PHASE) {
				lap_phi = n[alpha.index].laplacian;
				vector<double> vec_e_x_down = epsilon(n.get_neighbor_node(pf::Direction::x_down), alpha.index), vec_e_x_up = epsilon(n.get_neighbor_node(pf::Direction::x_up), alpha.index)
					, vec_e_y_down = epsilon(n.get_neighbor_node(pf::Direction::y_down), alpha.index), vec_e_y_up = epsilon(n.get_neighbor_node(pf::Direction::y_up), alpha.index)
					, vec_e = epsilon(n, alpha.index);
				vec_buff[0] = (vec_e_x_down[0] * vec_e_x_down[1] * n.get_neighbor_node(pf::Direction::x_down)[alpha.index].phi_grad[1]
					- vec_e_x_up[0] * vec_e_x_up[1] * n.get_neighbor_node(pf::Direction::x_up)[alpha.index].phi_grad[1]) / 2.0;
				vec_buff[1] = -(vec_e_y_down[0] * vec_e_y_down[1] * n.get_neighbor_node(pf::Direction::y_down)[alpha.index].phi_grad[0]
					- vec_e_y_up[0] * vec_e_y_up[1] * n.get_neighbor_node(pf::Direction::y_up)[alpha.index].phi_grad[0]) / 2.0;
				dfdPhiA = vec_buff[0] + vec_buff[1] - vec_e[0] * vec_e[0] * lap_phi - alpha.phaseFraction * (1.0 - alpha.phaseFraction) * (alpha.phaseFraction - 0.5 + m);
			}
			return dfdPhiA * int_width;
		}
		static double heatSource(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
			double incre_phi = 0.0, kappa = 1.8;
			for (auto alpha = n.begin(); alpha < n.end(); alpha++)
				if(alpha->phaseProperty == Solid)
					incre_phi += alpha->int_increment;
			return incre_phi * kappa;
		}
		static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 1 / 0.0003;
		}
	}

}
