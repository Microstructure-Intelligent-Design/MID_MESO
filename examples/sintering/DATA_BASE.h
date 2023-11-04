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

	namespace SYS_Sintering {
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
				phase.thermodynamic_model = pf::ThermodynamicModel::SolutionModel;
				phase.x.add_nodeEntry(CON, "CON");
				break;
			case VIRTUAL_PHASE:
				phase.phase_property = VIRTUAL_PHASE;
				phase.phase_name = "VIRTUAL_PHASE";
				phase.thermodynamic_model = pf::ThermodynamicModel::SolutionModel;
				phase.x.add_nodeEntry(CON, "CON");
				break;
			default:
				break;
			}
			return phase;
		}
		static void thermodynamicData(pf::PhaseNode& N, pf::Info_DynamicCollection& inf) {
			double con = N.x[CON].value, pf2 = 0.0, pf3 = 0.0, A = 16.0, B = 1.0, Kp = 5.0;
			double lap_c = (N.get_neighbor_node(pf::Direction::x_down).x[CON].value + N.get_neighbor_node(pf::Direction::x_up).x[CON].value + N.get_neighbor_node(pf::Direction::y_down).x[CON].value
				+ N.get_neighbor_node(pf::Direction::y_up).x[CON].value - 4.0 * con) / 0.5 / 0.5;
			for (auto phase = N.begin(); phase < N.end(); phase++) {
				if (phase->phaseProperty == VIRTUAL_PHASE)
					continue;
				pf2 += phase->phaseFraction * phase->phaseFraction;
				pf3 += phase->phaseFraction * phase->phaseFraction * phase->phaseFraction;
			}
			for (auto phase = N.begin(); phase < N.end(); phase++)
				phase->potential[CON].chemical_part = -0.5 * Kp * lap_c + 2.0 * A * con * (1.0 - 3.0 * con + 2.0 * con * con) + B * (2.0 * con - 6.0 * pf2 + 4.0 * pf3);
			return;
		}
		static double MolarVolume(pf::PhaseNode& N, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			return 1.0;
		}
		static void chemicalDiffusivity(pf::PhaseNode& n, pf::PhaseEntry& p, pf::Info_DynamicCollection& inf) {
			double phi = p.phaseFraction, con = n.x[CON].value, scale = phi * phi * phi * (10.0 - 15.0 * con + 6.0 * con * con), sum = 0.0;
			for (auto phase1 = n.begin(); phase1 < n.end(); phase1++)
				for (auto phase2 = n.begin(); phase2 < n.end(); phase2++)
					if (phase1->index != phase2->index && phase1->phaseProperty != VIRTUAL_PHASE && phase2->phaseProperty != VIRTUAL_PHASE)
						sum += phase1->phaseFraction * phase2->phaseFraction;
			p.kinetics_coeff.set(CON, CON, scale * 0.04 + 0.002 * (1.0 - scale) + 16.0 * con * (1.0 - con) + 1.6 * sum);
			return;
		}
		static double dfdPhi_bata_dfdPhi_alpha(pf::PhaseNode& n, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, double dx, double iwidth, pf::Info_DynamicCollection& inf) {
			double dfdPhiA = 0.0, dfdPhiB = 0.0, kphi = 2.0, B = 1.0, con = n.x[CON].value, pf2 = 0.0, lap_phi = 0.0;
			for (auto phase = n.begin(); phase < n.end(); phase++) {
				if (phase->phaseProperty == VIRTUAL_PHASE)
					continue;
				pf2 += phase->phaseFraction * phase->phaseFraction;
			}
			if (alpha.phaseProperty != VIRTUAL_PHASE) {
				lap_phi = (n.get_neighbor_node(pf::Direction::x_down)[alpha.index].phaseFraction + n.get_neighbor_node(pf::Direction::x_up)[alpha.index].phaseFraction + n.get_neighbor_node(pf::Direction::y_down)[alpha.index].phaseFraction
					+ n.get_neighbor_node(pf::Direction::y_up)[alpha.index].phaseFraction - 4.0 * alpha.phaseFraction) / dx / dx;
				dfdPhiA = -0.5 * kphi * lap_phi + B * (12.0 * (1.0 - con) * alpha.phaseFraction - 12.0 * (2.0 - con) * alpha.phaseFraction * alpha.phaseFraction + 12.0 * alpha.phaseFraction * pf2);
			}
			if (beta.phaseProperty != VIRTUAL_PHASE) {
				lap_phi = (n.get_neighbor_node(pf::Direction::x_down)[beta.index].phaseFraction + n.get_neighbor_node(pf::Direction::x_up)[beta.index].phaseFraction + n.get_neighbor_node(pf::Direction::y_down)[beta.index].phaseFraction
					+ n.get_neighbor_node(pf::Direction::y_up)[beta.index].phaseFraction - 4.0 * beta.phaseFraction) / dx / dx;
				dfdPhiB = -0.5 * kphi * lap_phi + B * (12.0 * (1.0 - con) * beta.phaseFraction - 12.0 * (2.0 - con) * beta.phaseFraction * beta.phaseFraction + 12.0 * beta.phaseFraction * pf2);
			}
			return dfdPhiB - dfdPhiA;
		}

		static double mobility(pf::PhaseNode& n, pf::PhaseEntry& alpha, pf::PhaseEntry& beta, pf::Info_DynamicCollection& inf) {
			return 5.0;
		}
	}

}
