#pragma once
#include "FieldStorage.h"
namespace pf {
	namespace DEFULT_FUNCTIONS {

		static double xi_ab(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 1.0;
		};
		static double xi_abc(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::PhaseEntry& gamma, pf::Info_DynamicCollection& inf) {
			return 0.0;
		};
		static Vector3 normals(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			Vector3 normal;
			normal[0] = alpha.phi_grad[0] * beta.phaseFraction - alpha.phaseFraction * beta.phi_grad[0];
			normal[1] = alpha.phi_grad[1] * beta.phaseFraction - alpha.phaseFraction * beta.phi_grad[1];
			normal[2] = alpha.phi_grad[2] * beta.phaseFraction - alpha.phaseFraction * beta.phi_grad[2];
			double length = sqrt(normal * normal);
			if (length < SYS_EPSILON) {
				normal[0] = 0.0;
				normal[1] = 0.0;
				normal[2] = 0.0;
			}
			else {
				normal[0] = normal[0] / length;
				normal[1] = normal[1] / length;
				normal[2] = normal[2] / length;
			};
			return normal;
		};
		static double GHSER(pf::PhaseNode& node, int comp_index) {
			return 0.0;
		}
		static void energy(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return;
		}
		static double molarVolume(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return 1.0;
		}
		static void potential(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return;
		}
		static pf::vStrain EffectiveEigenStrains(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
			pf::vStrain strain;
			return strain;
		}
		static pf::Matrix6x6 EffectiveElasticConstants(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
			pf::Matrix6x6 elasticConstants;
			return elasticConstants;
		}
		static double energyMinimizerIterator(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return 0.0;
		}
		static void chemicalMobility(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return;
		}
		static void chemicalDiffusivity(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return;
		}
		static void innerphaseReaction(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return;
		}
		static void intphaseReaction(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return;
		}
		static double phaseSource(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 0.0;
		}
		static double heatSource(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
			return 0.0;
		}
		static double mobility(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			if (alpha.index == beta.index)
				return 0.0;
			return 1.0;
		}
		static void nucleation(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
			return;
		}
		static void magnetizationIntensity(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
			return;
		}
		static void dMagnetizationIntensity_dPhi(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
			return;
		}
		static double interpolation_function_for_driving_force(pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return sqrt(alpha.phaseFraction * beta.phaseFraction);
		}
		static double Dfint_Dphi(pf::PhaseNode& node, pf::PhaseEntry& phase, double int_width, pf::Info_DynamicCollection& inf) {
			return 0.0;
		}
		static void Return_df_dphiB_dphiA(pf::PhaseNode& node, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double df_dphiB_dphiA, pf::Info_DynamicCollection& inf) {

			return;
		}
		static void write_scalar_customization(pf::ofstream& fout, pf::FieldStorage_forPhaseNode& fs) {
			//< example
			//fout << "<DataArray type = \"Float64\" Name = \"" << "custom" <<
			//	"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			//for (int k = 0; k < fs.limit_z; k++)
			//	for (int j = 0; j < fs.limit_y; j++)
			//		for (int i = 0; i < fs.limit_x; i++) {
			//			double fix = 0.0;
			//			fout << fix << endl;
			//		}
			//fout << "</DataArray>" << endl;
		}
		static void write_vec3_customization(pf::ofstream& fout, pf::FieldStorage_forPhaseNode& fs) {
			//< example
			//string name;
			//fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			//for (int k = 0; k < fs.limit_z; k++)
			//	for (int j = 0; j < fs.limit_y; j++)
			//		for (int i = 0; i < fs.limit_x; i++) {
			//			fout << vec3[0] << " "
			//				<< vec3[1] << " "
			//				<< vec3[2] << endl;
			//		}
			//fout << "</DataArray>" << endl;
		}
	}

	class Funcs {
	public:
		double(*Ghser)(pf::PhaseNode&, int index);
		void(*Energy)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		void(*Potential)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*MolarVolume)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*EnergyMinimizerIterator)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		void(*ChemicalMobility)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		void(*ChemicalDiffusivity)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		vStrain(*EffectiveEigenStrains)(pf::PhaseNode&, pf::Info_DynamicCollection&);
		Matrix6x6(*EffectiveElasticConstants)(pf::PhaseNode&, pf::Info_DynamicCollection&);
		void(*MagnetizationIntensity)(pf::PhaseNode&, pf::Info_DynamicCollection&);
		void(*dMagnetizationIntensitydPhi)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		void(*InterPhasesReaction)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		void(*InnerPhaseReaction)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*PhaseSource)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*HeatSource)(pf::PhaseNode&, pf::Info_DynamicCollection&);
		double(*Xi_ab)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*Xi_abc)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*Mobility)(pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		double(*Interpolation_function_for_driving_force)(pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		Vector3(*Normals)(pf::PhaseEntry&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
		void(*Nucleation)(pf::PhaseNode&, pf::Info_DynamicCollection&);
		// double is int_width
		double(*dfint_dphi)(pf::PhaseNode&, pf::PhaseEntry&, double, pf::Info_DynamicCollection&);
		double (*dfint_dphi_S2009)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, double, pf::Info_DynamicCollection&);
		void(*Return_df_dphiB_dphiA)(pf::PhaseNode&, pf::PhaseEntry&, pf::PhaseEntry&, double, pf::Info_DynamicCollection&);
		// write file
		void(*write_scalar_customization)(pf::ofstream&, pf::FieldStorage_forPhaseNode&);
		void(*write_vec3_customization)(pf::ofstream&, pf::FieldStorage_forPhaseNode&);
		Funcs() {
			// parameters
			Ghser = DEFULT_FUNCTIONS::GHSER;
			Energy = DEFULT_FUNCTIONS::energy;
			Potential = DEFULT_FUNCTIONS::potential;
			MolarVolume = DEFULT_FUNCTIONS::molarVolume;
			EnergyMinimizerIterator = DEFULT_FUNCTIONS::energyMinimizerIterator;
			ChemicalMobility = DEFULT_FUNCTIONS::chemicalMobility;
			ChemicalDiffusivity = DEFULT_FUNCTIONS::chemicalDiffusivity;
			Mobility = DEFULT_FUNCTIONS::mobility;
			Interpolation_function_for_driving_force = DEFULT_FUNCTIONS::interpolation_function_for_driving_force;
			EffectiveEigenStrains = DEFULT_FUNCTIONS::EffectiveEigenStrains;
			EffectiveElasticConstants = DEFULT_FUNCTIONS::EffectiveElasticConstants;
			Nucleation = DEFULT_FUNCTIONS::nucleation;
			MagnetizationIntensity = DEFULT_FUNCTIONS::magnetizationIntensity;
			dMagnetizationIntensitydPhi = DEFULT_FUNCTIONS::dMagnetizationIntensity_dPhi;
			// evolution equation terms
			InterPhasesReaction = DEFULT_FUNCTIONS::intphaseReaction;
			InnerPhaseReaction = DEFULT_FUNCTIONS::innerphaseReaction;
			PhaseSource = DEFULT_FUNCTIONS::phaseSource;
			HeatSource = DEFULT_FUNCTIONS::heatSource;
			Xi_ab = DEFULT_FUNCTIONS::xi_ab;
			Xi_abc = DEFULT_FUNCTIONS::xi_abc;
			Normals = DEFULT_FUNCTIONS::normals;
			dfint_dphi = DEFULT_FUNCTIONS::Dfint_Dphi;
			Return_df_dphiB_dphiA = DEFULT_FUNCTIONS::Return_df_dphiB_dphiA;
			write_scalar_customization = DEFULT_FUNCTIONS::write_scalar_customization;
			write_vec3_customization = DEFULT_FUNCTIONS::write_vec3_customization;

		}
		void clear() {
			Ghser = nullptr;
			Energy = nullptr;
			Potential = nullptr;
			MolarVolume = nullptr;
			EnergyMinimizerIterator = nullptr;
			ChemicalMobility = nullptr;
			ChemicalDiffusivity = nullptr;
			InterPhasesReaction = nullptr;
			PhaseSource = nullptr;
			InnerPhaseReaction = nullptr;
			Xi_ab = nullptr;
			Xi_abc = nullptr;
			Mobility = nullptr;
			Interpolation_function_for_driving_force = nullptr;
			Normals = nullptr;
			EffectiveEigenStrains = nullptr;
			EffectiveElasticConstants = nullptr;
			Nucleation = nullptr;
			MagnetizationIntensity = nullptr;
			dMagnetizationIntensitydPhi = nullptr;
			dfint_dphi = nullptr;
			Return_df_dphiB_dphiA = nullptr;
			dfint_dphi_S2009 = nullptr;
			HeatSource = nullptr;
			write_scalar_customization = nullptr;
			write_vec3_customization = nullptr;
		}
		Funcs& operator=(const Funcs& n) {
			Ghser = n.Ghser;
			Energy = n.Energy;
			Potential = n.Potential;
			MolarVolume = n.MolarVolume;
			EnergyMinimizerIterator = n.EnergyMinimizerIterator;
			ChemicalMobility = n.ChemicalMobility;
			ChemicalDiffusivity = n.ChemicalDiffusivity;
			InterPhasesReaction = n.InterPhasesReaction;
			PhaseSource = n.PhaseSource;
			InnerPhaseReaction = n.InnerPhaseReaction;
			Xi_ab = n.Xi_ab;
			Xi_abc = n.Xi_abc;
			Mobility = n.Mobility;
			Interpolation_function_for_driving_force = n.Interpolation_function_for_driving_force;
			Normals = n.Normals;
			EffectiveEigenStrains = n.EffectiveEigenStrains;
			EffectiveElasticConstants = n.EffectiveElasticConstants;
			Nucleation = n.Nucleation;
			MagnetizationIntensity = n.MagnetizationIntensity;
			dMagnetizationIntensitydPhi = n.dMagnetizationIntensitydPhi;
			dfint_dphi = n.dfint_dphi;
			Return_df_dphiB_dphiA = n.Return_df_dphiB_dphiA;
			dfint_dphi_S2009 = n.dfint_dphi_S2009;
			HeatSource = n.HeatSource;
			write_scalar_customization = n.write_scalar_customization;
			write_vec3_customization = n.write_vec3_customization;
			return *this;
		}
	};

}


