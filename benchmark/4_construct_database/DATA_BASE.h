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
/*###########################################################
			SER (Standard elements references)
###########################################################*/
	namespace ghser {
		static double GHSERAL(double T) {
			double val = 0.0;
			if (T < 298.14 || T > 2900.00) {
				printerror();
			}
			else if (T >= 298.14 && T < 700.00) {
				val = -7976.15 + 137.093038 * T - 24.3671976 * T * log(T) - 0.001884662 * pow(T, 2) - 8.77664E-07 * pow(T, 3) + 74092 * pow(T, -1);
			}
			else if (T >= 700.00 && T < 933.47) {
				val = -11276.24 + 223.048446 * T - 38.5844296 * T * log(T) + 0.018531982 * pow(T, 2) - 5.764227E-06 * pow(T, 3) + 74092 * pow(T, -1);
			}
			else if (T >= 933.47 && T <= 2900.00) {
				val = -11278.361 + 188.684136 * T - 31.748192 * T * log(T) - 1.230622E+28 * pow(T, -9);
			}
			return val;
		}
		static double GLIQAL(double T) {
			double val = 0.0;
			if (T < 298.15 || T > 2900.00) {
				printerror();
			}
			else if (T >= 298.15 && T < 933.47) {
				val = +GHSERAL(T) + 11005.045 - 11.84185 * T + 7.9337E-20 * pow(T, 7);
			}
			else if (T >= 933.47 && T < 2900.00) {
				val = -795.991 + 177.430209 * T - 31.748192 * T * log(T);
			}
			return val;
		}
		static double GBCCAL(double T) {
			double val = 0.0;
			if (T < 298.14 || T > 2900.00) {
				printerror();
			}
			else if (T >= 298.14 && T < 2900.00) {
				val = 10083 - 4.813 * T + GHSERAL(T);
			}
			return val;
		}
		static double GHCPAL(double T) {
			double val = 0.0;
			if (T < 298.14 || T > 2900.00) {
				printerror();
			}
			else if (T >= 298.14 && T < 2900.00) {
				val = 5481 - 1.8 * T + GHSERAL(T);
			}
			return val;
		}

		static double GHSERCU(double T) {
			double val = 0.0;
			if (T < 298.14 || T > 3200.00) {
				printerror();
			}
			else if (T >= 298.14 && T < 1357.77) {
				val = -7770.458 + 130.485235 * T - 24.112392 * T * log(T) - 0.00265684 * pow(T, 2) + 1.29223E-07 * pow(T, 3) + 52478 * pow(T, -1);
			}
			else if (T >= 1357.77 && T <= 3200.00) {
				val = -13542.026 + 183.803828 * T - 31.38 * T * log(T) + 3.64167E+29 * pow(T, -9);
			}
			return val;
		}
		static double GLIQCU(double T) {
			double val = 0.0;
			if (T < 298.15 || T > 2900.00) {
				printerror();
			}
			else if (T >= 298.15 && T < 1357.77) {
				val = +GHSERCU(T) + 12964.735 - 9.511904 * T - 5.8489E-21 * pow(T, 7);
			}
			else if (T >= 1357.77 && T < 3200.00) {
				val = -46.545 + 173.881484 * T - 31.38 * T * log(T);
			}
			return val;
		}
		static double GBCCCU(double T) {
			double val = 0.0;
			if (T < 298.14 || T > 3200.00) {
				printerror();
			}
			else if (T >= 298.14 && T < 3200.00) {
				val = +4017 - 1.255 * T + GHSERCU(T);
			}
			return val;
		}
		static double GHCPCU(double T) {
			double val = 0.0;
			if (T < 298.14 || T > 3200.00) {
				printerror();
			}
			else if (T >= 298.14 && T < 3200.00) {
				val = +600 + 0.2 * T + GHSERCU(T);
			}
			return val;
		}
	};

	namespace SYS_AL_CU {
		enum PHASES { FCC_A1, ALCU_THETA};
		enum SUBLATTICES { _AL_CU_, _VA_, _AL_, _CU_ };
		namespace TDB {
		/***************************************************************************
		PHASE LIQUID  %  1  1.0  !
		CONSTITUENT LIQUID  :AL,CU :  !
		****************************************************************************/
			static double G_LIQUID_AL_0(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return +ghser::GLIQAL(T);
			}
			static double G_LIQUID_CU_0(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return +ghser::GLIQCU(T);
			}
			static double G_LIQUID_AL_CU_0(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return -68236.397 + 9.43 * T;
			}
			static double G_LIQUID_AL_CU_1(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return 33155.04 - 9.87 * T;
			}
			static double G_LIQUID_AL_CU_2(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return -2240.57;
			}
			//---------------------------------------------------------------------------------
			//					 Gibbs Energy of LIQUID		
			//---------------------------------------------------------------------------------
			static double G_LIQUID(pf::ConNode& x, double T, double P) {
				double xAl = x[Al].value;
				double xCu = x[Cu].value;
				double val = G_LIQUID_AL_0(T) * xAl
					+ G_LIQUID_CU_0(T) * xCu
					+ (G_LIQUID_AL_CU_0(T)
						+ G_LIQUID_AL_CU_1(T) * (xAl - xCu)
						+ G_LIQUID_AL_CU_2(T) * (xAl - xCu) * (xAl - xCu)) * xAl * xCu
					+ R * T * (xAl * log(xAl) + xCu * log(xCu));
				return val;
			}
			static double dG_LIQUID_dAL(pf::ConNode& x, double T, double P) {
				double xAl = x[Al].value;
				double xCu = x[Cu].value;
				double val = G_LIQUID_AL_0(T)
					+ (G_LIQUID_AL_CU_0(T)
						+ G_LIQUID_AL_CU_1(T) * (2.0 * xAl - xCu)
						+ G_LIQUID_AL_CU_2(T) * (3.0 * xAl - xCu) * (xAl - xCu)) * xCu
					+ R * T * (log(xAl) + 1.0);
				return val;
			}
			static double dG_LIQUID_dCU(pf::ConNode& x, double T, double P) {
				double xAl = x[Al].value;
				double xCu = x[Cu].value;
				double val = G_LIQUID_CU_0(T)
					+ (G_LIQUID_AL_CU_0(T)
						+ G_LIQUID_AL_CU_1(T) * (xAl - 2.0 * xCu)
						+ G_LIQUID_AL_CU_2(T) * (xAl - 3.0 * xCu) * (xAl - xCu)) * xAl
					+ R * T * (log(xCu) + 1.0);
				return val;
			}
	/***************************************************************************
	PHASE FCC_A1  % 2 1   1 !
	CONSTITUENT FCC_A1  :AL,CU : VA :  !
	****************************************************************************/
			static double G_FCC_A1_AL_0(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return ghser::GHSERAL(T);
			}
			static double G_FCC_A1_CU_0(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return ghser::GHSERCU(T);
			}
			static double G_FCC_A1_AL_CU_0(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return -45391.85 - 5.45 * T;
			}
			static double G_FCC_A1_AL_CU_1(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return 62616.56 - 27.06 * T;
			}
			static double G_FCC_A1_AL_CU_2(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return 20554.69 - 22.32 * T;
			}
			static double G_FCC_A1_AL_CU_3(double T) {
				if (T < 298.14 || T > 6000.00) { printerror(); }
				return -8579.66;
			}
			//---------------------------------------------------------------------------------
			//					 Gibbs Energy of FCC_A1		
			//---------------------------------------------------------------------------------

			static double G_FCC_A1(pf::ConNode& x, double T, double P) {
				double xAl = x[Al].value;
				double xCu = x[Cu].value;
				double val = G_FCC_A1_AL_0(T) * xAl
					+ G_FCC_A1_CU_0(T) * xCu
					+ (G_FCC_A1_AL_CU_0(T)
						+ G_FCC_A1_AL_CU_1(T) * (xAl - xCu)
						+ G_FCC_A1_AL_CU_2(T) * (xAl - xCu) * (xAl - xCu)
						+ G_FCC_A1_AL_CU_3(T) * (xAl - xCu) * (xAl - xCu) * (xAl - xCu)) * xAl * xCu
					+ R * T * (xAl * log(xAl) + xCu * log(xCu));
				return val;
			}
			static double dG_FCC_A1_dAL(pf::ConNode& x, double T, double P) {
				double xAl = x[Al].value;
				double xCu = x[Cu].value;
				double val = G_FCC_A1_AL_0(T)
					+ (G_FCC_A1_AL_CU_0(T)
						+ G_FCC_A1_AL_CU_1(T) * (2.0 * xAl - xCu)
						+ G_FCC_A1_AL_CU_2(T) * (3.0 * xAl - xCu) * (xAl - xCu)
						+ G_FCC_A1_AL_CU_3(T) * (4.0 * xAl - xCu) * (xAl - xCu) * (xAl - xCu)) * xCu
					+ R * T * (log(xAl) + 1.0);
				return val;
			}
			static double dG_FCC_A1_dCU(pf::ConNode& x, double T, double P) {
				double xAl = x[Al].value;
				double xCu = x[Cu].value;
				double val = G_FCC_A1_CU_0(T)
					+ (G_FCC_A1_AL_CU_0(T)
						+ G_FCC_A1_AL_CU_1(T) * (xAl - 2.0 * xCu)
						+ G_FCC_A1_AL_CU_2(T) * (xAl - 3.0 * xCu) * (xAl - xCu)
						+ G_FCC_A1_AL_CU_3(T) * (xAl - 4.0 * xCu) * (xAl - xCu) * (xAl - xCu)) * xAl
					+ R * T * (log(xCu) + 1.0);
				return val;
			}
			/***************************************************************************
			PHASE AL2CU  %  2  0.667  0.333 !
			CONSTITUENT THETA  :AL:AL,CU :  !
			****************************************************************************/
			static double G_ALCU_THETA_AL_AL_0(double T) {//G(ALCU_THETA,AL:AL;0) REF0
				if (T < 2.98150E+02 || T > 6.00000E+03) { printerror(); }
				return +ghser::GHSERAL(T) + 10083 - 4.813 * T;
			}
			static double G_ALCU_THETA_AL_CU_0(double T) {//G(ALCU_THETA,AL:CU;0) REF0
				if (T < 2.98150E+02 || T > 6.00000E+03) { printerror(); }
				return 0.667 * ghser::GHSERAL(T) + 0.333 * ghser::GHSERCU(T) - 15700 + 1.9 * T;
			}
			static double G_ALCU_THETA_AL_AL_CU_0(double T) {//G(ALCU_THETA,AL:AL,CU;0) REF0
				if (T < 2.98150E+02 || T > 6.00000E+03) { printerror(); }
				return 1000;
			}
			//---------------------------------------------------------------------------------
			//					 Gibbs Energy of THETA		
			//---------------------------------------------------------------------------------
			static double G_ALCU_THETA(pf::SublatticeNode& y, double T, double P) {//G(ALCU_THETA,AL:AL,CU;0) REF0
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double val = G_ALCU_THETA_AL_AL_0(T) * y1Al * y2Al
					+ G_ALCU_THETA_AL_CU_0(T) * y2Cu * y1Al
					+ G_ALCU_THETA_AL_AL_CU_0(T) * y1Al * y2Cu * y2Al
					+ 0.667 * R * T * y1Al * log(y1Al) + 0.333 * R * T * (y2Cu * log(y2Cu) + y2Al * log(y2Al));
				return val;
			}
			static double dG_ALCU_THETA_dAL1(pf::SublatticeNode& y, double T, double P) {//G(ALCU_THETA,AL:AL,CU;0) REF0
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double val = G_ALCU_THETA_AL_AL_0(T) * y2Al
					+ G_ALCU_THETA_AL_CU_0(T) * y2Cu
					+ G_ALCU_THETA_AL_AL_CU_0(T) * y2Cu * y2Al
					+ 0.667 * R * T * (log(y1Al) + 1.0);
				return val;
			}
			static double dG_ALCU_THETA_dAL2(pf::SublatticeNode& y, double T, double P) {//G(ALCU_THETA,AL:AL,CU;0) REF0
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double val = G_ALCU_THETA_AL_AL_0(T) * y1Al
					+ G_ALCU_THETA_AL_AL_CU_0(T) * y1Al * y2Cu
					+ 0.333 * R * T * (log(y2Al) + 1.0);
				return val;
			}
			static double dG_ALCU_THETA_dCU2(pf::SublatticeNode& y, double T, double P) {//G(ALCU_THETA,AL:AL,CU;0) REF0
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double val = G_ALCU_THETA_AL_CU_0(T) * y1Al
					+ G_ALCU_THETA_AL_AL_CU_0(T) * y1Al * y2Al
					+ 0.333 * R * T * (log(y2Cu) + 1.0);
				return val;
			}
			/***************************************************************************
			PHASE GAMA1  %  3 0.3077  0.0769   0.6154 !
			CONSTITUENT GAMA1  :AL:AL,CU:CU : !
			****************************************************************************/
			static double G_GAMA1_AL_AL_CU_0(double T) {
				if (T < 2.98150E+02 || T > 2.90000E+03) { printerror(); }
				return 0.3846 * ghser::GHSERAL(T) + 0.6154 * ghser::GHSERCU(T) - 24185 + 31 * T - 4 * T * log(T);
			}
			static double G_GAMA1_AL_CU_CU_0(double T) {
				if (T < 2.98150E+02 || T > 2.90000E+03) { printerror(); }
				return 0.3077 * ghser::GHSERAL(T) + 0.6923 * ghser::GHSERCU(T) - 20622 + 28.61 * T - 4 * T * log(T);
			}
			//---------------------------------------------------------------------------------
			//					 Gibbs Energy of THETA		
			//---------------------------------------------------------------------------------
			static double G_GAMA1(pf::SublatticeNode& y, double T, double P) {
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double y3Cu = y[_CU_][Cu].value;
				double val = G_GAMA1_AL_AL_CU_0(T) * y1Al * y2Al * y3Cu
					+ G_GAMA1_AL_CU_CU_0(T) * y1Al * y2Cu * y3Cu
					+ 0.3077 * R * T * y1Al * log(y1Al) 
					+ 0.0769 * R * T * (y2Cu * log(y2Cu) + y2Al * log(y2Al)) 
					+ 0.6154 * R * T * y3Cu * log(y3Cu);
				return val;
			}
			static double dG_GAMA1_dAL1(pf::SublatticeNode& y, double T, double P) {
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double y3Cu = y[_CU_][Cu].value;
				double val = G_GAMA1_AL_AL_CU_0(T) * y2Al * y3Cu
					+ G_GAMA1_AL_CU_CU_0(T) * y2Cu * y3Cu
					+ 0.3077 * R * T * (log(y1Al) + 1.0);
				return val;
			}
			static double dG_GAMA1_dAL2(pf::SublatticeNode& y, double T, double P) {
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double y3Cu = y[_CU_][Cu].value;
				double val = G_GAMA1_AL_AL_CU_0(T) * y1Al * y3Cu
					+ 0.0769 * R * T * (log(y2Al) + 1.0);
				return val;
			}
			static double dG_GAMA1_dCU2(pf::SublatticeNode& y, double T, double P) {
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double y3Cu = y[_CU_][Cu].value;
				double val = G_GAMA1_AL_CU_CU_0(T) * y1Al * y3Cu
					+ 0.0769 * R * T * (log(y2Cu) +1.0);
				return val;
			}
			static double dG_GAMA1_dCU3(pf::SublatticeNode& y, double T, double P) {
				double y1Al = y[_AL_][Al].value;
				double y2Cu = y[_AL_CU_][Cu].value;
				double y2Al = y[_AL_CU_][Al].value;
				double y3Cu = y[_CU_][Cu].value;
				double val = G_GAMA1_AL_AL_CU_0(T) * y1Al * y2Al
					+ G_GAMA1_AL_CU_CU_0(T) * y1Al * y2Cu
					+ 0.6154 * R * T * (log(y3Cu) + 1.0);
				return val;
			}
		}
		static pf::Info_Phase Get_Phase_Structure(PHASES phase_property) {
			pf::Info_Phase phase;
			switch (phase_property)
			{
			case FCC_A1:
				phase.phase_property = FCC_A1;
				phase.phase_name = "FCC_A1";
				phase.thermodynamic_model = pf::ThermodynamicModel::SolutionModel;
				phase.solvent = Al;
				phase.x.add_nodeEntry(Al, "Al");
				phase.x.add_nodeEntry(Cu, "Cu");
				break;
			case ALCU_THETA:
				phase.phase_property = ALCU_THETA;
				phase.phase_name = "ALCU_THETA";
				phase.thermodynamic_model = pf::ThermodynamicModel::SublatticeModel;
				phase.solvent = Al;
				phase.x.add_nodeEntry(Al, "Al");
				phase.x.add_nodeEntry(Cu, "Cu");
				phase.sublattice.add_Node(_AL_, "_AL_", 2.0 / 3.0);
				phase.sublattice.add_Node(_AL_CU_, "_AL_CU_", 1.0 / 3.0);
				phase.sublattice[_AL_].add_nodeEntry(Al, "Al");
				phase.sublattice[_AL_CU_].add_nodeEntry(Al, "Al");
				phase.sublattice[_AL_CU_].add_nodeEntry(Cu, "Cu");
				break;
			default:
				break;
			}
			return phase;
		}
		static double MolarVolume(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			return 1.0e-5;
		}
		static void Energy(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			int phaseProperty = p.phaseProperty;
			pf::ConNode& x = p.x;
			pf::SublatticeNode& y = p.sublattice;
			double T = N.tempValues.temperature, P = 0.0, vm = MolarVolume(N, p, inf);
			switch (phaseProperty)
			{
			case FCC_A1:
				p.chemEnergyDensity = TDB::G_FCC_A1(x, T, P) / vm;
				break;
			case ALCU_THETA:
				p.chemEnergyDensity = TDB::G_ALCU_THETA(y, T, P) / vm;
				break;
			default:
				std::cout << "Functions error!" << std::endl;
				SYS_PROGRAM_STOP;
				break;
			}
		}
		static void Potential(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			int phaseProperty = p.phaseProperty;
			pf::ConNode& x = p.x;
			pf::SublatticeNode& y = p.sublattice;
			double T = N.tempValues.temperature, P = 0.0, vm = MolarVolume(N, p, inf);
			switch (phaseProperty)
			{
			case FCC_A1:
				p.potential[Cu].chemical_part = TDB::dG_FCC_A1_dCU(x, T, P) / vm - TDB::dG_FCC_A1_dAL(x, T, P) / vm;
				break;
			case ALCU_THETA:
				p.potential[Cu].chemical_part = 3.0 * TDB::dG_ALCU_THETA_dCU2(y, T, P) / vm - 3.0 * TDB::dG_ALCU_THETA_dAL2(y, T, P) / vm;
				break;
			default:
				std::cout << "Functions error!" << std::endl;
				SYS_PROGRAM_STOP;
				break;
			}
		}
		static double EnergyMinimizerIterator(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			pf::ConNode& x = p.x;
			if (p.phaseProperty == ALCU_THETA) {
				pf::SublatticeNode& y = p.sublattice;
				y[_AL_][Al].value = 1.0;
				y[_AL_CU_][Al].value = (x[Al].value - y[_AL_].siteNum) / y[_AL_CU_].siteNum;
				y[_AL_CU_][Cu].value = x[Cu].value / y[_AL_CU_].siteNum;
			}
			return 0.0;
		}
		static void ChemicalDiffusivity(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			int phaseProperty = p.phaseProperty;
			pf::ConNode& x = p.x;
			double t = n.tempValues.temperature;
			switch (phaseProperty)
			{
			case FCC_A1:
				p.kinetics_coeff.set(Al, Al, 1e-11);
				p.kinetics_coeff.set(Cu, Cu, 1e-11);
				break;
			case ALCU_THETA:
				p.kinetics_coeff.set(Al, Al, 1e-11);
				p.kinetics_coeff.set(Cu, Cu, 1e-11);
				break;
			default:
				std::cout << "Functions error!" << std::endl;
				SYS_PROGRAM_STOP;
				break;
			}
		}
		//////////////////////////////////////////////////////////////////////////////////
		static void IntphaseReaction(pf::PhaseNode& n, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			double pre = 0.0, reaction_rate = 1e-4;
			for (auto comp1 = alpha.x.begin(); comp1 < alpha.x.end(); comp1++)
				for (auto comp2 = beta.x.begin(); comp2 < beta.x.end(); comp2++)
					if (comp1->index == comp2->index && comp1->index != Al) {
						pre = reaction_rate * alpha.phaseFraction * beta.phaseFraction * (beta.potential[comp2->index].value - alpha.potential[comp1->index].value);
						if (pre > 0 && comp2->value > X_EPSILON && comp1->value < (1.0 - X_EPSILON)) {
							comp1->ChemicalReactionFlux += pre / alpha.phaseFraction;
							comp2->ChemicalReactionFlux -= pre / beta.phaseFraction;
						}
						else if (pre < 0 && comp1->value > X_EPSILON && comp2->value < (1.0 - X_EPSILON)) {
							comp1->ChemicalReactionFlux += pre / alpha.phaseFraction;
							comp2->ChemicalReactionFlux -= pre / beta.phaseFraction;
						}
					}
		}
		static double xi_AB(pf::PhaseNode& n, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 1e3;
		}
		static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 1e-15;
		}
	}

}
