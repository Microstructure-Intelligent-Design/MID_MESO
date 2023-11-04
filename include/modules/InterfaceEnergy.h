#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"

namespace pf {
	class InterfaceEnergy {
	public:
		InterfaceEnergy() {};
		~InterfaceEnergy() {};
		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information);
		void optimize_phiGradient(pf::PhaseNode& node, int phaseIndex, int optimize_range);
		void cal_interface_inrement_inNode(pf::PhaseNode& node, bool adjust_phi_0_1 = false);
		void clear();

		FieldStorage_forPhaseNode* simulationField;
		Information* information;
		Int_Gradient interface_gradient;
		Int_Potential interface_potential;
		double dfint_dphi_S2009(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double intface_width, pf::Info_DynamicCollection& inf) {
			Funcs& func = information->materialSystem.functions;
			return func.Xi_ab(node, alpha, beta, inf) * (intface_width * (beta.phaseFraction * alpha.laplacian - alpha.phaseFraction * beta.laplacian)
				+ PI * PI / 2.0 / intface_width * (alpha.phaseFraction - beta.phaseFraction));
		};
		double dfint_dphi_grad_S1996(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
			Funcs& func = information->materialSystem.functions;
			double grad = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++)
				if (beta->index != phase.index && beta->_flag) {
					grad += 2.0 * intface_width * func.Xi_ab(node, phase, *beta, inf) * (beta->phi_grad * beta->phi_grad * phase.phaseFraction
						- phase.phi_grad * beta->phi_grad * beta->phaseFraction + phase.phaseFraction * beta->phaseFraction * beta->laplacian - beta->phaseFraction * beta->phaseFraction * phase.laplacian);
				}
			return grad;
		};
		double dfint_dphi_grad_S1999(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
			Funcs& func = information->materialSystem.functions;
			double grad = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++)
				if (beta->index != phase.index && beta->_flag) {
					grad += intface_width * func.Xi_ab(node, phase, *beta, inf) * beta->laplacian;
				}
			return grad;
		};
		double dfint_dphi_pot_Nobstacle(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
			Funcs& func = information->materialSystem.functions;
			double grad = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++) {
				if (beta->index == phase.index || beta->_flag == pf_BULK)
					continue;
				grad += 16.0 / intface_width / PI / PI * func.Xi_ab(node, phase, *beta, inf) * beta->phaseFraction;
				for (auto gamma = beta + 1; gamma < node.end(); gamma++) {
					if (gamma->index == phase.index || gamma->_flag == pf_BULK)
						continue;
					grad += func.Xi_abc(node, phase, *beta, *gamma, inf) * beta->phaseFraction * gamma->phaseFraction / intface_width;
				}
			}
			return grad;
		};
		double dfint_dphi_pot_Nwell(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
			Funcs& func = information->materialSystem.functions;
			double grad = 0.0;
			for (auto beta = node.begin(); beta < node.end(); beta++) {
				if (beta->index == phase.index || beta->_flag == pf_BULK)
					continue;
				grad += 18.0 / intface_width * func.Xi_ab(node, phase, *beta, inf) * phase.phaseFraction * beta->phaseFraction * beta->phaseFraction;
				for (auto gamma = beta + 1; gamma < node.end(); gamma++) {
					if (gamma->index == phase.index || gamma->_flag == pf_BULK)
						continue;
					grad += 2.0 / intface_width * func.Xi_abc(node, phase, *beta, *gamma, inf) * phase.phaseFraction
						* beta->phaseFraction * gamma->phaseFraction * beta->phaseFraction * gamma->phaseFraction;
				}
			}
			return grad;
		};
		double dfint_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, double intface_width, pf::Info_DynamicCollection& inf) {
			Funcs& func = information->materialSystem.functions;
			if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1996(node, phase, intface_width, inf) + dfint_dphi_pot_Nobstacle(node, phase, intface_width, inf);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1996(node, phase, intface_width, inf) + dfint_dphi_pot_Nwell(node, phase, intface_width, inf);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1999(node, phase, intface_width, inf) + dfint_dphi_pot_Nobstacle(node, phase, intface_width, inf);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1999(node, phase, intface_width, inf) + dfint_dphi_pot_Nwell(node, phase, intface_width, inf);
			}
			else {
				if (func.dfint_dphi == nullptr)
					return 0.0;
				else
					return func.dfint_dphi(node, phase, intface_width, inf);
			}
		};
	};
}