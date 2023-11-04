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
	/****************** TiC-ZrC system *********************************/
	namespace SYS_TIC_ZRC {
		enum PHASES { FCC_A1};
		static double MolarVolume(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			double mV = 1.31e-5;
			return mV;
		}
		static double Energy(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			int phaseProperty = p.phaseProperty;
			double xTI = p.x[Ti].value;
			double T = N.tempValues.temperature;
			double P = 0.0, vm = MolarVolume(N, p, inf);
			switch (phaseProperty)
			{
			case FCC_A1:
				return (6.565e5 * pow(xTI, 6) - 2.01e6 * pow(xTI, 5) + 2.379e6 * pow(xTI, 4) - 1.377e6 * pow(xTI, 3)
					+ 3.636e5 * pow(xTI, 2) - 698.2 * xTI - 2.21e5) / vm;
				break;
			default:
				std::cout << "Functions error!" << std::endl;
				SYS_PROGRAM_STOP;
				break;
			}
		}

		static pf::Info_Phase Get_Phase_Structure(PHASES phase_property) {
			pf::Info_Phase phase;
			pf::Matrix6x6 matrix;
			double Gm = 200e6, v = 0.3, Az = 1.0, C12 = 2 * v * Gm / (1 - 2 * v), C11 = 2 * Gm / Az + C12;

			matrix(0, 0) = C11;
			matrix(1, 1) = C11;
			matrix(2, 2) = C11;
			matrix(3, 3) = Gm;
			matrix(4, 4) = Gm;
			matrix(5, 5) = Gm;
			matrix(0, 1) = C12;
			matrix(1, 0) = C12;
			matrix(0, 2) = C12;
			matrix(2, 0) = C12;
			matrix(1, 2) = C12;
			matrix(2, 1) = C12;
			switch (phase_property)
			{
			case FCC_A1:
				//C = Get_PhaseElasticConstants(FCC_A1);
				phase.phase_property = FCC_A1;
				phase.phase_name = "FCC_A1";
				phase.thermodynamic_model = pf::ThermodynamicModel::SolutionModel;
				phase.x.add_nodeEntry(Ti, "xTi");
				phase.elastic_constants = matrix;
				break;
			default:
				break;
			}
			return phase;
		}

		static pf::vStrain EffectivePhaseEigenStrains(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
			pf::vStrain strain;
			double eta0 = 0.0;
			double c = n[0].x[Ti].value;
			if (c < 0)
				c = 0;
			else if (c > 1)
				c = 1;
			eta0 = 0.01 * c * c * c * (10.0 - 15.0 * c + 6.0 * c * c);
			//eta0 = 0.01 - 0.01 * c * c * c * (10.0 - 15.0 * c + 6.0 * c * c);
			strain[0] = eta0;
			strain[1] = eta0;
			//strain[2] = eta0;
			return strain;
		}

		static void Potential(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			double c = p.x[Ti].value, sigma = 0.8e10; // sigma = interface energy / dx
			double lap_c = (N.get_neighbor_node(pf::Direction::x_up)[0].x[Ti].value + N.get_neighbor_node(pf::Direction::x_down)[0].x[Ti].value
				+ N.get_neighbor_node(pf::Direction::y_up)[0].x[Ti].value + N.get_neighbor_node(pf::Direction::y_down)[0].x[Ti].value
				- 4.0 * c);
			p.potential[Ti].chemical_part = (3.939e6 * pow(c, 5) - 1.005e7 * pow(c, 4) + 9.516e6 * pow(c, 3)
				- 4.131e6 * c * c + 7.272e5 * c - 698.2) / MolarVolume(N, p, inf) - 2.0 * sigma * lap_c;
			pf::vStrain strain;
			strain[0] = 0.01 * (30 * c * c - 60 * c * c * c + 30 * c * c * c * c);
			strain[1] = strain[0];
			p.potential[Ti].elastic_part = -(N.mechanicalValues.Stresses * strain);

		}

		static void ChemicalMobility(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			p.kinetics_coeff.set(Ti, Ti, 1.0e-11);
		}
	}

}
