#pragma once
#include "../../include/baseTools/baseTools.h"
namespace TDB {
	using namespace pf::materials;

	static void printerror(){
		cout << "Temperature, out of limit, thermodynamic calculation error!" << endl;
		SYS_PROGRAM_STOP;
	}
	static void printerror(double T) {
		cout << "Temperature, out of limit, thermodynamic calculation error!" << endl;
		SYS_PROGRAM_STOP;
	}
	const double R = 8.3144521;

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
		enum SUBLATTICES { _AL_CU_, _VA_, _AL_, _CU_ };
		
		/***************************************************************************
										PHASE FCC_A1
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
		//									Gibbs Energy
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
		//---------------------------------------------------------------------------------
		//						A first order partial derivatives
		//---------------------------------------------------------------------------------
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
		//---------------------------------------------------------------------------------
		//						A Second order partial derivatives
		//---------------------------------------------------------------------------------
		static double dG_FCC_A1_dALdAL(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;
			double val = 0.0;
			return val;
		}
		static double dG_FCC_A1_dALdCU(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;
			double val = 0.0;
			return val;
		}
		static double dG_FCC_A1_dCUdAL(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;
			double val = 0.0;
			return val;
		}
		static double dG_FCC_A1_dCUdCU(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;
			double val = 0.0;
			return val;
		}
		//---------------------------------------------------------------------------------
		//						D I F F U S I O N     D A T A
		//---------------------------------------------------------------------------------
		static double Mk_FCC_Al(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;

			double MQ = 0.0;

			return exp(MQ / R / T) / R / T;
		}
		static double Mk_FCC_Cu(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;

			double MQ = 0.0;

			return exp(MQ / R / T) / R / T;
		}
		/***************************************************************************
										PHASE AL2CU
		****************************************************************************/
		static double G_ALCU_THETA_AL_AL_0(double T) {
			if (T < 2.98150E+02 || T > 6.00000E+03) { printerror(); }
			return +ghser::GHSERAL(T) + 10083 - 4.813 * T;
		}
		static double G_ALCU_THETA_AL_CU_0(double T) {
			if (T < 2.98150E+02 || T > 6.00000E+03) { printerror(); }
			return 0.667 * ghser::GHSERAL(T) + 0.333 * ghser::GHSERCU(T) - 15700 + 1.9 * T;
		}
		static double G_ALCU_THETA_AL_AL_CU_0(double T) {
			if (T < 2.98150E+02 || T > 6.00000E+03) { printerror(); }
			return 1000;
		}
		//---------------------------------------------------------------------------------
		//								Gibbs Energy	
		//---------------------------------------------------------------------------------
		static double G_ALCU_THETA(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = G_ALCU_THETA_AL_AL_0(T) * y1Al * y2Al
				+ G_ALCU_THETA_AL_CU_0(T) * y2Cu * y1Al
				+ G_ALCU_THETA_AL_AL_CU_0(T) * y1Al * y2Cu * y2Al
				+ 0.667 * R * T * y1Al * log(y1Al) + 0.333 * R * T * (y2Cu * log(y2Cu) + y2Al * log(y2Al));
			return val;
		}
		//---------------------------------------------------------------------------------
		//					A first order partial derivatives	
		//---------------------------------------------------------------------------------
		static double dG_ALCU_THETA_dAL1(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = G_ALCU_THETA_AL_AL_0(T) * y2Al
				+ G_ALCU_THETA_AL_CU_0(T) * y2Cu
				+ G_ALCU_THETA_AL_AL_CU_0(T) * y2Cu * y2Al
				+ 0.667 * R * T * (log(y1Al) + 1.0);
			return val;
		}
		static double dG_ALCU_THETA_dAL2(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = G_ALCU_THETA_AL_AL_0(T) * y1Al
				+ G_ALCU_THETA_AL_AL_CU_0(T) * y1Al * y2Cu
				+ 0.333 * R * T * (log(y2Al) + 1.0);
			return val;
		}
		static double dG_ALCU_THETA_dCU2(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = G_ALCU_THETA_AL_CU_0(T) * y1Al
				+ G_ALCU_THETA_AL_AL_CU_0(T) * y1Al * y2Al
				+ 0.333 * R * T * (log(y2Cu) + 1.0);
			return val;
		}
		//---------------------------------------------------------------------------------
		//						A Second order partial derivatives
		//---------------------------------------------------------------------------------
		static double dG_ALCU_THETA_dALdAL(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = 0.0;
			return val;
		}
		static double dG_ALCU_THETA_dALdCU(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = 0.0;
			return val;
		}
		static double dG_ALCU_THETA_dCUdAL(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = 0.0;
			return val;
		}
		static double dG_ALCU_THETA_dCUdCU(pf::SublatticeNode& y, double T, double P) {
			double y1Al = y[_AL_][Al].value;
			double y2Cu = y[_AL_CU_][Cu].value;
			double y2Al = y[_AL_CU_][Al].value;
			double val = 0.0;
			return val;
		}
		//---------------------------------------------------------------------------------
		//						D I F F U S I O N     D A T A
		//---------------------------------------------------------------------------------
		static double Mk_ALCU_THETA_Al(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;

			double MQ = 0.0;

			return exp(MQ / R / T) / R / T;
		}
		static double Mk_ALCU_THETA_Cu(pf::ConNode& x, double T, double P) {
			double xAl = x[Al].value;
			double xCu = x[Cu].value;

			double MQ = 0.0;

			return exp(MQ / R / T) / R / T;
		}

	}

}
