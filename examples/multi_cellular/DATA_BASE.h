#pragma once
#include "../../include/baseTools/baseTools.h"
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
	namespace SYS_Cellular {
		enum PHASES { Cell, BackGround };
		static double u[] = {0.0, 0.0, 0.0};
		const int cells_numbers = 58;
		const double cell_radius = 3.0, cell_distance = 6.0;

		static double dfint_dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, double int_width, pf::Info_DynamicCollection& inf) {
			if (phase.phaseProperty == BackGround)
				return 0.0;
			double dfint_dphi = 0.0, gamma_n = 0.0, lambda = 7.0, kappa = 60.0, miu = 40.0, ksi = 1.5e3;
			if (phase.index < 6)
				gamma_n = 5.0;
			else
				gamma_n = 2.5;
			// first term
			dfint_dphi += gamma_n * (2.0 * phase.laplacian
				+ 60.0 / lambda / lambda * phase.phaseFraction * (1 - phase.phaseFraction) * (1 - 2.0 * phase.phaseFraction));
			// second term
			for (auto phi = node.begin(); phi < node.end(); phi++)
				if (phi->index != phase.index)
					dfint_dphi += 60.0 * kappa / lambda / lambda * phase.phaseFraction * phi->phaseFraction * phi->phaseFraction;
			// third term
			dfint_dphi -= 2.0 * miu / /*PI*/4.0 / cell_radius / cell_radius * phase.phaseFraction 
				* (inf.user_defined_value[phase.index] - /*PI*/4.0 * cell_radius * cell_radius);
			// forth term
			pf::Vector3 u_tot(u[0], u[1], u[2]);
			u_tot[0] += 60.0 * kappa / lambda / lambda / ksi * inf.user_defined_vec3(phase.index, 0);
			u_tot[1] += 60.0 * kappa / lambda / lambda / ksi * inf.user_defined_vec3(phase.index, 1);
			u_tot[2] += 60.0 * kappa / lambda / lambda / ksi * inf.user_defined_vec3(phase.index, 2);
			dfint_dphi -= 2.0 * (u_tot * phase.phi_grad);
			return dfint_dphi;
		}
		static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 0.5;
		}
	}

}
