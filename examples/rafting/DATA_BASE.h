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

	namespace SYS_Elastic {
		enum PHASES { Matrix, Precipitate };
		enum COMPONENT{ NONE = pf::materials::ElementTop, CON};

		static pf::Matrix6x6 Get_PhaseElasticConstants(PHASES phase_property) {
			pf::Matrix6x6 matrix;
			//< cubic
			if (phase_property == Matrix) {
				double Gm = 400, v = 0.3, Az = 3.0, C12 = /*2 * */v * Gm / (1 - 2 * v), C11 = 2 * Gm / Az + C12;

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
			}
			else if (phase_property == Precipitate) {
				double Gm = 200, v = 0.3, Az = 3.0, C12 = /*2 * */v * Gm / (1 - 2 * v), C11 = 2 * Gm / Az + C12;

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
			}
			return matrix;
		}
		static pf::Matrix6x6 EffectiveElasticConstants(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
			pf::Matrix6x6 elasticConstants, pC, mC, deltC, effectC;
			pC = Get_PhaseElasticConstants(Matrix);
			mC = Get_PhaseElasticConstants(Precipitate);
			deltC = pC - mC;
			effectC = (pC + mC) * 0.5;
			double c = n[0].x[CON].value;
			elasticConstants = effectC + deltC * (c * c * c * (10.0 - 15.0 * c + 6.0 * c * c) - 0.5);
			return elasticConstants;
		}
		static pf::Matrix6x6 dEffectiveElasticConstants_dc(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
			pf::Matrix6x6 elasticConstants, pC, mC, deltC;
			pC = Get_PhaseElasticConstants(Matrix);
			mC = Get_PhaseElasticConstants(Precipitate);
			deltC = pC - mC;
			double c = n[0].x[CON].value;
			elasticConstants = deltC * c * c * (10.0 - 15.0 * c + 6.0 * c * c) * 3.0 + deltC * c * c * c * (- 15.0 + 12.0 * c);
			return elasticConstants;
		}
		static pf::vStrain Get_PhaseElasticStrain(PHASES phase_property) {
			pf::vStrain strain;

			if (phase_property == Matrix) {
				strain[0] = 0.0;
				strain[1] = 0.0;
				strain[2] = 0.0;
				strain[3] = 0.0;
				strain[4] = 0.0;
				strain[5] = 0.0;
			}
			else if (phase_property == Precipitate) {
				strain[0] = 0.01;
				strain[1] = 0.01;
				strain[2] = 0.01;
				strain[3] = 0.00;
				strain[4] = 0.00;
				strain[5] = 0.00;
			}
			return strain;
		}
		static pf::Info_Phase Get_Phase_Structure(PHASES phase_property) {
			pf::Info_Phase phase;
			switch (phase_property)
			{
			case Matrix:
				phase.phase_property = Matrix;
				phase.phase_name = "Matrix";
				phase.thermodynamic_model = pf::ThermodynamicModel::SolutionModel;
				phase.x.add_nodeEntry(CON, "CON");
				break;
			case Precipitate:
				phase.phase_property = Precipitate;
				phase.phase_name = "Precipitate";
				phase.thermodynamic_model = pf::ThermodynamicModel::SolutionModel;
				phase.x.add_nodeEntry(CON, "CON");
				break;
			default:
				break;
			}
			return phase;
		}
		static pf::vStrain EffectivePhaseEigenStrains(pf::PhaseNode& n, pf::Info_DynamicCollection& inf) {
			pf::vStrain strain;
			double c = n[0].x[CON].value;
			strain = Get_PhaseElasticStrain(Precipitate) * c * c * c * (10.0 - 15.0 * c + 6.0 * c * c);
			return strain;
		}
		static pf::vStrain dEffectivePhaseEigenStrains_dc(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			pf::vStrain strain;
			double c = n[0].x[CON].value;
			strain = Get_PhaseElasticStrain(Precipitate) * c * c * (10.0 - 15.0 * c + 6.0 * c * c) * 3.0 
				+ Get_PhaseElasticStrain(Precipitate) * c * c * c * (- 15.0 + 12.0 * c);
			return strain;
		}
		static void Potential(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			//pf::Matrix6x6 delt_C = Get_PhaseElasticConstants(Precipitate) - Get_PhaseElasticConstants(Matrix);
			double Ab = 1.0, c = p.x[CON].value, k = 1.0;
			double lap_c = (N.get_neighbor_node(pf::Direction::x_up)[0].x[CON].value + N.get_neighbor_node(pf::Direction::x_down)[0].x[CON].value
				+ N.get_neighbor_node(pf::Direction::y_up)[0].x[CON].value + N.get_neighbor_node(pf::Direction::y_down)[0].x[CON].value
				- 4.0 * c);
			pf::Matrix6x6 dC_dc = dEffectiveElasticConstants_dc(N, inf);
			pf::vStrain dE_dc = dEffectivePhaseEigenStrains_dc(N, p, inf);
			p.potential[CON].chemical_part = 2.0 * Ab * c * (1 - c) * (1 - 2 * c) - 2.0 * k * lap_c;
			/*p.potential[CON].elastic_part = (delt_C * N.mechanicalValues.Strains * N.mechanicalValues.Strains * 0.5 * 30.0 * c * c * (1 - c) * (1 - c)
				- N.mechanicalValues.Stresses * Get_PhaseElasticStrain(Precipitate) * 30.0 * c * c * (1 - c) * (1 - c)) / 1.0;*/
			p.potential[CON].elastic_part = -(dE_dc * N.mechanicalValues.Stresses)
				+ 0.5 * (dC_dc * (N.mechanicalValues.Strains - N.mechanicalValues.EffectiveEigenStrain)
					* (N.mechanicalValues.Strains - N.mechanicalValues.EffectiveEigenStrain));
		}
		static double MolarVolume(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			return 1.0;
		}
		static void ChemicalMobility(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			p.kinetics_coeff.set(CON, CON, 0.01);
		}
	}

}
