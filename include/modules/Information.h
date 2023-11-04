#pragma once
#include "../baseTools/baseTools.h"

namespace pf {

	struct Set_Detials {
		double phi_incre_limit;
		double con_incre_limit;
		double energy_limit_iterate_precision;
		double flux_error_balance_coefficient;
		// 1 : normal phi, 0 : average phi
		int con_average_range;
		int bulk_incre_average_range;
		DifferenceMethod difference_method;
		PhaseFluxModel flux_model;
		Int_Gradient int_grad;
		Int_Potential int_pot;
		PHASES_GROW_SHRINK_TYPE phases_grow_shringk_type;
		int OMP_thread_counts;

		//<  bc[Boundary][Index_bc_axisX(y, z, limit_y)][PhaseIndex];
		vector<vector<pf::PhaseNode>> MPF_FixedBoundaryCondition;
		Set_Detials() {
			phi_incre_limit = 0.001;
			energy_limit_iterate_precision = 1.0;
			con_incre_limit = 0.001;
			flux_error_balance_coefficient = 1.0;
			con_average_range = 1;
			bulk_incre_average_range = 1;
			difference_method = pf::DifferenceMethod::FIVE_POINT;
			flux_model = pf::PhaseFluxModel::IntDiff_ConGrad;
			int_grad = Int_Gradient::Steinbach_G2009;
			int_pot = Int_Potential::Steinbach_P2009;
			phases_grow_shringk_type = PHASES_GROW_SHRINK_TYPE::NONE;
			OMP_thread_counts = 1;
			MPF_FixedBoundaryCondition.resize(6);
		}
		~Set_Detials() {
			MPF_FixedBoundaryCondition.clear();
		}
		Set_Detials& operator=(const Set_Detials& n) {
			phi_incre_limit = n.phi_incre_limit;
			energy_limit_iterate_precision = n.energy_limit_iterate_precision;
			con_incre_limit = n.con_incre_limit;
			flux_error_balance_coefficient = n.flux_error_balance_coefficient;
			con_average_range = n.con_average_range;
			bulk_incre_average_range = n.bulk_incre_average_range;
			difference_method = n.difference_method;
			int_grad = n.int_grad;
			int_pot = n.int_pot;
			flux_model = n.flux_model;
			OMP_thread_counts = n.OMP_thread_counts;
			MPF_FixedBoundaryCondition = n.MPF_FixedBoundaryCondition;
			phases_grow_shringk_type = n.phases_grow_shringk_type;
			return *this;
		}
		void automaticlly_set_FixedBoundaryCondition(Set_Disperse _disperse_settings, Info_Phases _phases, Info_Node sys_x, Funcs _funcs, Info_DynamicCollection _dynamic, PhaseEntry matrix_phase, NucleationBox nucleationBox) {
			double vec[] = { 0.0, 0.0, 0.0 };
			MPF_FixedBoundaryCondition.resize(6);
			if (_disperse_settings.x_bc == pf::FIXED) {
				int size = _disperse_settings.Ny * _disperse_settings.Nz;
				MPF_FixedBoundaryCondition[pf::Boundary::UP_X].resize(size);
				MPF_FixedBoundaryCondition[pf::Boundary::DOWN_X].resize(size);
				for (int index = 0; index < size; index++) {
					PhaseNode& up_x = MPF_FixedBoundaryCondition[pf::Boundary::UP_X][index];
					PhaseNode& down_x = MPF_FixedBoundaryCondition[pf::Boundary::DOWN_X][index];
					for (auto x = sys_x.begin(); x < sys_x.end(); x++) {
						up_x.x.add_con(x->index, x->value);
						down_x.x.add_con(x->index, x->value);
					}
					up_x.magneticValues.dmag_potential_dphi.add_double(0, 0.0);
					up_x.magneticValues.dmag_intensity_dphi.add_vec(0, vec);
					up_x.magneticValues.dmag_charge_density_dphi.add_double(0, 0.0);
					down_x.magneticValues.dmag_potential_dphi.add_double(0, 0.0);
					down_x.magneticValues.dmag_intensity_dphi.add_vec(0, vec);
					down_x.magneticValues.dmag_charge_density_dphi.add_double(0, 0.0);
					up_x.add_phase(0, matrix_phase.phaseProperty, pf_NEAR_INTERFACE, 0.0);
					down_x.add_phase(0, matrix_phase.phaseProperty, pf_NEAR_INTERFACE, 0.0);
					for (auto x = matrix_phase.x.begin(); x < matrix_phase.x.end(); x++) {
						up_x[0].x.add_con(x->index, x->value);
						up_x[0].potential.add_con(x->index, 0.0);
						down_x[0].x.add_con(x->index, x->value);
						down_x[0].potential.add_con(x->index, 0.0);
					}
					for (auto sub = _phases[matrix_phase.phaseProperty].sublattice.begin(); sub < _phases[matrix_phase.phaseProperty].sublattice.end(); sub++) {
						up_x[matrix_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						down_x[matrix_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						for (auto con = sub->begin(); con < sub->end(); con++) {
							up_x[matrix_phase.index].sublattice[sub->index].add_y(con->index);
							down_x[matrix_phase.index].sublattice[sub->index].add_y(con->index);
						}
					}
					for (auto n = _phases[matrix_phase.phaseProperty].n.begin(); n < _phases[matrix_phase.phaseProperty].n.end(); n++) {
						up_x[matrix_phase.index].n.add_x(n->index, 0.0);
						down_x[matrix_phase.index].n.add_x(n->index, 0.0);
					}

					for (auto p = nucleationBox.nucleus_box.begin(); p < nucleationBox.nucleus_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_x[p->phaseIndex].x.add_con(x->index, x->value);
							up_x[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_x[p->phaseIndex].x.add_con(x->index, x->value);
							down_x[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_x[p->phaseIndex].n.add_x(n->index, 0.0);
							down_x[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.condition_check_phase_box.begin(); p < nucleationBox.condition_check_phase_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_x[p->phaseIndex].x.add_con(x->index, x->value);
							up_x[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_x[p->phaseIndex].x.add_con(x->index, x->value);
							down_x[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_x[p->phaseIndex].n.add_x(n->index, 0.0);
							down_x[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.geometricRegion_box.begin(); p < nucleationBox.geometricRegion_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_x[p->phaseIndex].x.add_con(x->index, x->value);
							up_x[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_x[p->phaseIndex].x.add_con(x->index, x->value);
							down_x[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_x[p->phaseIndex].n.add_x(n->index, 0.0);
							down_x[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.pointSet_box.begin(); p < nucleationBox.pointSet_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_x.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_x.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_x.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_x[p->phaseIndex].x.add_con(x->index, x->value);
							up_x[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_x[p->phaseIndex].x.add_con(x->index, x->value);
							down_x[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_x[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_x[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_x[p->phaseIndex].n.add_x(n->index, 0.0);
							down_x[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
				}
			}
			if (_disperse_settings.y_bc == pf::FIXED) {
				int size = _disperse_settings.Nx * _disperse_settings.Nz;
				MPF_FixedBoundaryCondition[pf::Boundary::UP_Y].resize(size);
				MPF_FixedBoundaryCondition[pf::Boundary::DOWN_Y].resize(size);
				for (int index = 0; index < size; index++) {
					PhaseNode& up_y = MPF_FixedBoundaryCondition[pf::Boundary::UP_Y][index];
					PhaseNode& down_y = MPF_FixedBoundaryCondition[pf::Boundary::DOWN_Y][index];
					for (auto x = sys_x.begin(); x < sys_x.end(); x++) {
						up_y.x.add_con(x->index, x->value);
						down_y.x.add_con(x->index, x->value);
					}
					up_y.magneticValues.dmag_potential_dphi.add_double(0, 0.0);
					up_y.magneticValues.dmag_intensity_dphi.add_vec(0, vec);
					up_y.magneticValues.dmag_charge_density_dphi.add_double(0, 0.0);
					down_y.magneticValues.dmag_potential_dphi.add_double(0, 0.0);
					down_y.magneticValues.dmag_intensity_dphi.add_vec(0, vec);
					down_y.magneticValues.dmag_charge_density_dphi.add_double(0, 0.0);
					up_y.add_phase(0, matrix_phase.phaseProperty, pf_NEAR_INTERFACE, 0.0);
					down_y.add_phase(0, matrix_phase.phaseProperty, pf_NEAR_INTERFACE, 0.0);
					for (auto x = matrix_phase.x.begin(); x < matrix_phase.x.end(); x++) {
						up_y[0].x.add_con(x->index, x->value);
						up_y[0].potential.add_con(x->index, 0.0);
						down_y[0].x.add_con(x->index, x->value);
						down_y[0].potential.add_con(x->index, 0.0);
					}
					for (auto sub = _phases[matrix_phase.phaseProperty].sublattice.begin(); sub < _phases[matrix_phase.phaseProperty].sublattice.end(); sub++) {
						up_y[matrix_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						down_y[matrix_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						for (auto con = sub->begin(); con < sub->end(); con++) {
							up_y[matrix_phase.index].sublattice[sub->index].add_y(con->index);
							down_y[matrix_phase.index].sublattice[sub->index].add_y(con->index);
						}
					}
					for (auto n = _phases[matrix_phase.phaseProperty].n.begin(); n < _phases[matrix_phase.phaseProperty].n.end(); n++) {
						up_y[matrix_phase.index].n.add_x(n->index, 0.0);
						down_y[matrix_phase.index].n.add_x(n->index, 0.0);
					}

					for (auto p = nucleationBox.nucleus_box.begin(); p < nucleationBox.nucleus_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_y[p->phaseIndex].x.add_con(x->index, x->value);
							up_y[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_y[p->phaseIndex].x.add_con(x->index, x->value);
							down_y[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_y[p->phaseIndex].n.add_x(n->index, 0.0);
							down_y[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.condition_check_phase_box.begin(); p < nucleationBox.condition_check_phase_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_y[p->phaseIndex].x.add_con(x->index, x->value);
							up_y[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_y[p->phaseIndex].x.add_con(x->index, x->value);
							down_y[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_y[p->phaseIndex].n.add_x(n->index, 0.0);
							down_y[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.geometricRegion_box.begin(); p < nucleationBox.geometricRegion_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_y[p->phaseIndex].x.add_con(x->index, x->value);
							up_y[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_y[p->phaseIndex].x.add_con(x->index, x->value);
							down_y[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_y[p->phaseIndex].n.add_x(n->index, 0.0);
							down_y[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.pointSet_box.begin(); p < nucleationBox.pointSet_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_y.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_y.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_y.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_y[p->phaseIndex].x.add_con(x->index, x->value);
							up_y[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_y[p->phaseIndex].x.add_con(x->index, x->value);
							down_y[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_y[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_y[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_y[p->phaseIndex].n.add_x(n->index, 0.0);
							down_y[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
				}
			}
			if (_disperse_settings.z_bc == pf::FIXED) {
				int size = _disperse_settings.Nx * _disperse_settings.Ny;
				MPF_FixedBoundaryCondition[pf::Boundary::UP_Z].resize(size);
				MPF_FixedBoundaryCondition[pf::Boundary::DOWN_Z].resize(size);
				for (int index = 0; index < size; index++) {
					PhaseNode& up_z = MPF_FixedBoundaryCondition[pf::Boundary::UP_Z][index];
					PhaseNode& down_z = MPF_FixedBoundaryCondition[pf::Boundary::DOWN_Z][index];
					for (auto x = sys_x.begin(); x < sys_x.end(); x++) {
						up_z.x.add_con(x->index, x->value);
						down_z.x.add_con(x->index, x->value);
					}
					up_z.magneticValues.dmag_potential_dphi.add_double(0, 0.0);
					up_z.magneticValues.dmag_intensity_dphi.add_vec(0, vec);
					up_z.magneticValues.dmag_charge_density_dphi.add_double(0, 0.0);
					down_z.magneticValues.dmag_potential_dphi.add_double(0, 0.0);
					down_z.magneticValues.dmag_intensity_dphi.add_vec(0, vec);
					down_z.magneticValues.dmag_charge_density_dphi.add_double(0, 0.0);
					up_z.add_phase(0, matrix_phase.phaseProperty, pf_NEAR_INTERFACE, 0.0);
					down_z.add_phase(0, matrix_phase.phaseProperty, pf_NEAR_INTERFACE, 0.0);
					for (auto x = matrix_phase.x.begin(); x < matrix_phase.x.end(); x++) {
						up_z[0].x.add_con(x->index, x->value);
						up_z[0].potential.add_con(x->index, 0.0);
						down_z[0].x.add_con(x->index, x->value);
						down_z[0].potential.add_con(x->index, 0.0);
					}
					for (auto sub = _phases[matrix_phase.phaseProperty].sublattice.begin(); sub < _phases[matrix_phase.phaseProperty].sublattice.end(); sub++) {
						up_z[matrix_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						down_z[matrix_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						for (auto con = sub->begin(); con < sub->end(); con++) {
							up_z[matrix_phase.index].sublattice[sub->index].add_y(con->index);
							down_z[matrix_phase.index].sublattice[sub->index].add_y(con->index);
						}
					}
					for (auto n = _phases[matrix_phase.phaseProperty].n.begin(); n < _phases[matrix_phase.phaseProperty].n.end(); n++) {
						up_z[matrix_phase.index].n.add_x(n->index, 0.0);
						down_z[matrix_phase.index].n.add_x(n->index, 0.0);
					}

					for (auto p = nucleationBox.nucleus_box.begin(); p < nucleationBox.nucleus_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_z[p->phaseIndex].x.add_con(x->index, x->value);
							up_z[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_z[p->phaseIndex].x.add_con(x->index, x->value);
							down_z[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_z[p->phaseIndex].n.add_x(n->index, 0.0);
							down_z[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.condition_check_phase_box.begin(); p < nucleationBox.condition_check_phase_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_z[p->phaseIndex].x.add_con(x->index, x->value);
							up_z[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_z[p->phaseIndex].x.add_con(x->index, x->value);
							down_z[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_z[p->phaseIndex].n.add_x(n->index, 0.0);
							down_z[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.geometricRegion_box.begin(); p < nucleationBox.geometricRegion_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_z[p->phaseIndex].x.add_con(x->index, x->value);
							up_z[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_z[p->phaseIndex].x.add_con(x->index, x->value);
							down_z[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_z[p->phaseIndex].n.add_x(n->index, 0.0);
							down_z[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
					for (auto p = nucleationBox.pointSet_box.begin(); p < nucleationBox.pointSet_box.end(); p++) {
						if (p->phaseIndex == 0)
							continue;
						up_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						up_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						up_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_potential_dphi.add_double(p->phaseIndex, 0.0);
						down_z.magneticValues.dmag_intensity_dphi.add_vec(p->phaseIndex, vec);
						down_z.magneticValues.dmag_charge_density_dphi.add_double(p->phaseIndex, 0.0);
						up_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						down_z.add_phase(p->phaseIndex, p->phaseProperty, pf_NEAR_INTERFACE, 0.0);
						for (auto x = _phases[p->phaseProperty].x.begin(); x < _phases[p->phaseProperty].x.end(); x++) {
							up_z[p->phaseIndex].x.add_con(x->index, x->value);
							up_z[p->phaseIndex].potential.add_con(x->index, 0.0);
							down_z[p->phaseIndex].x.add_con(x->index, x->value);
							down_z[p->phaseIndex].potential.add_con(x->index, 0.0);
						}
						for (auto sub = _phases[p->phaseProperty].sublattice.begin(); sub < _phases[p->phaseProperty].sublattice.end(); sub++) {
							up_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							down_z[p->phaseIndex].sublattice.new_sub(sub->index, sub->siteNum);
							for (auto con = sub->begin(); con < sub->end(); con++) {
								up_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
								down_z[p->phaseIndex].sublattice[sub->index].add_y(con->index);
							}
						}
						for (auto n = _phases[p->phaseProperty].n.begin(); n < _phases[p->phaseProperty].n.end(); n++) {
							up_z[p->phaseIndex].n.add_x(n->index, 0.0);
							down_z[p->phaseIndex].n.add_x(n->index, 0.0);
						}
					}
				}
			}
		}
		bool set_FixedBoundaryCondition_by_direction(Set_Disperse _disperse_settings, Info_Phases _phases, Funcs _funcs, Info_DynamicCollection _dynamic, pf::Boundary _bc_direction, PhaseEntry _phase, double temperature = 273, double electric_potential = 0.0) {
			int size = 0;
			if (_bc_direction == pf::Boundary::UP_X || _bc_direction == pf::Boundary::DOWN_X) {
				if (_disperse_settings.x_bc != pf::FIXED)
					return false;
				size = _disperse_settings.Ny * _disperse_settings.Nz;
			}
			else if (_bc_direction == pf::Boundary::UP_Y || _bc_direction == pf::Boundary::DOWN_Y) {
				if (_disperse_settings.y_bc != pf::FIXED)
					return false;
				size = _disperse_settings.Nx * _disperse_settings.Nz;
			}
			else if (_bc_direction == pf::Boundary::UP_Z || _bc_direction == pf::Boundary::DOWN_Z) {
				if (_disperse_settings.z_bc != pf::FIXED)
					return false;
				size = _disperse_settings.Nx * _disperse_settings.Ny;
			}
			for (int index = 0; index < size; index++) {
				bool is_new_phase = true;
				PhaseNode& node = MPF_FixedBoundaryCondition[_bc_direction][index];
				node.tempValues.temperature = temperature;
				node.electricValues.elec_potential = electric_potential;
				for (auto phase = node.begin(); phase < node.end(); phase++)
					if (phase->index == _phase.index) {
						phase->set(_phase.index, _phase.phaseProperty, _phase._flag, _phase.phaseFraction);
						for (auto x = _phase.x.begin(); x < _phase.x.end(); x++) {
							phase->x[x->index].value = x->value;
						}
						is_new_phase = false;
						break;
					}
				if (is_new_phase) {
					node.add_phase(_phase.index, _phase.phaseProperty, _phase._flag, _phase.phaseFraction);
					for (auto x = _phase.x.begin(); x < _phase.x.end(); x++) {
						node[_phase.index].x.add_con(x->index, x->value);
						node[_phase.index].potential.add_con(x->index, 0.0);
						for (auto xx = _phase.x.begin(); xx < _phase.x.end(); xx++)
							node[_phase.index].kinetics_coeff.set(x->index, xx->index, 0.0);
					}
					for (auto sub = _phases[_phase.phaseProperty].sublattice.begin(); sub < _phases[_phase.phaseProperty].sublattice.end(); sub++) {
						node[_phase.index].sublattice.new_sub(sub->index, sub->siteNum);
						for (auto con = sub->begin(); con < sub->end(); con++)
							node[_phase.index].sublattice[sub->index].add_y(con->index);
					}
					for (auto n = _phases[_phase.phaseProperty].n.begin(); n < _phases[_phase.phaseProperty].n.end(); n++)
						node[_phase.index].n.add_x(n->index, 0.0);
				}
			}
			return true;
		}
	};
	struct Info_Settings {
		Set_Disperse disperse_settings;
		Set_OutputFile file_settings;
		Set_Detials details_settings;
		Info_Settings& operator=(const Info_Settings& n) {
			disperse_settings = n.disperse_settings;
			file_settings = n.file_settings;
			details_settings = n.details_settings;
			return *this;
		}
	};
	struct Info_MaterialSystem {
		Info_Phases phases;
		double init_temperature;
		double R;
		Info_Node sys_x;
		PhaseEntry matrix_phase;

		Funcs functions;

		bool is_mechanics_on;
		bool is_magnetics_on;
		bool is_electrics_on;
		bool is_fluid_on;
		Set_ElectricFieldMask electricFieldMask_settings;
		Set_MagneticFieldMask magneticFieldMask_settings;
		Set_FluidFieldMask fluidFieldMask_settings;
		Info_Mechanics mechanics;

		Info_MaterialSystem() {
			matrix_phase.phaseProperty = VIRTUAL_PHASE;
			is_mechanics_on = false;
			is_magnetics_on = false;
			is_electrics_on = false;
			is_fluid_on = false;
			R = 8.314;
			init_temperature = 0.0;
		}
		void clear() {
			phases.clear();
			functions.clear();
		}
		Info_MaterialSystem& operator=(const Info_MaterialSystem& n) {
			matrix_phase = n.matrix_phase;
			functions = n.functions;
			phases = n.phases;
			mechanics = n.mechanics;
			R = n.R;
			init_temperature = n.init_temperature;
			sys_x = n.sys_x;
			is_mechanics_on = n.is_mechanics_on;
			is_magnetics_on = n.is_magnetics_on;
			is_electrics_on = n.is_electrics_on;
			is_fluid_on = n.is_fluid_on;
			electricFieldMask_settings = n.electricFieldMask_settings;
			magneticFieldMask_settings = n.magneticFieldMask_settings;
			fluidFieldMask_settings = n.fluidFieldMask_settings;
			return *this;
		}
	};
	struct point_in_region_index {
		point_in_region_index(int re, int gi) {
			region = re;
			grain_index = gi;
		}
		int region;
		int grain_index;
	};
	class Information
	{
	public:
		Information() {};
		Information(const Information& _information) {
			init(_information);
		};
		void init(const Information& _information) {
			dynamicCollection = _information.dynamicCollection;
			settings = _information.settings;
			materialSystem = _information.materialSystem;
			nucleationBox = _information.nucleationBox;
			if (settings.details_settings.OMP_thread_counts >= omp_get_num_procs()) {
				settings.details_settings.OMP_thread_counts = omp_get_num_procs() - 1;
			}
			settings.file_settings.log_file_name = settings.file_settings.log_file_name;
			writer.init(settings.file_settings.working_folder_path, settings.file_settings.NaN);
			data_writer.set_path(settings.file_settings.working_folder_path);
			data_writer.set_mainName(settings.file_settings.data_file_name);
			writer.init_txt_file(settings.file_settings.log_file_name);
		};
		void statistics_information_in_phaseMesh(FieldStorage_forPhaseNode& phaseMesh);
		void normalize_phi_in_mesh(FieldStorage_forPhaseNode& phaseMesh);
		void set_init_con_ratio(FieldStorage_forPhaseNode& phaseMesh) {
			statistics_information_in_phaseMesh(phaseMesh);
			dynamicCollection.init_x_ratio = dynamicCollection.inf_node.x;
		}
		void statistics_information_phaseFraction_in_phaseMesh(FieldStorage_forPhaseNode& phaseMesh);
		void generate_voronoi_structure(Vector3 box_position, Vector3 box_size, int grain_0_index, int grain_number, int generate_step,vector<int> phases_properties, vector<double> phases_weight, XNode x, double temperature = 0.0);
		void generate_structure_from_BMP_pic(int phaseIndex, int phase_property, int generate_step, XNode x, double threshold[2], double temperature = 0.0, string fileName = "");
		string get_file_time_tail() {
			string tail = timer::return_cunrrent_time_by_string();
			for (auto i = tail.begin(); i < tail.end();) {
				if (*i == '#' || *i == '\n') {
					i = tail.erase(i);
				}
				else if (*i == ':') {
					*i = '-';
				}
				else {
					i++;
				}
			}
			return tail;
		}
		void clear() {
			materialSystem.clear();
		}

		Information& operator=(Information& n);
		Info_DynamicCollection dynamicCollection;
		Info_Settings settings;
		Info_MaterialSystem materialSystem;
		NucleationBox nucleationBox;
		//Tools
		WriteToFile writer;
		DataFile data_writer;
	};
	inline Information& Information::operator=(Information& n) {
		dynamicCollection = n.dynamicCollection;
		settings = n.settings;
		materialSystem = n.materialSystem;
		nucleationBox = n.nucleationBox;
		writer = n.writer;
		data_writer = n.data_writer;
	}

	inline void Information::statistics_information_in_phaseMesh(FieldStorage_forPhaseNode& phaseMesh) {
		double node_number = double(phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z);
		PhaseNode& statistics_node = dynamicCollection.inf_node;
		statistics_node = phaseMesh(0, 0, 0);
		statistics_node.tempValues.temperature = 0.0;
		for (auto p = statistics_node.begin(); p < statistics_node.end(); p++) {
			p->phaseFraction = 0.0;
			p->chemEnergyDensity = 0.0;
			p->elasEnergyDensity = 0.0;
			for (auto c = p->x.begin(); c < p->x.end(); c++) {
				c->value = 0.0;
				c->PhaseTransitionFlux = 0.0;
				c->ChemicalReactionFlux = 0.0;
				c->DiffusionFlux = 0.0;
			}
			for (auto c = p->potential.begin(); c < p->potential.end(); c++) {
				c->value = 0.0;
				c->chemical_part = 0.0;
				c->elastic_part = 0.0;
			}
		}
		for (auto c = statistics_node.x.begin(); c < statistics_node.x.end(); c++) {
			c->value = 0.0;
			c->PhaseTransitionFlux = 0.0;
			c->ChemicalReactionFlux = 0.0;
			c->DiffusionFlux = 0.0;
		}
		for (auto p = statistics_node.potential.begin(); p < statistics_node.potential.end(); p++) {
			p->value = 0.0;
		}

		for (int x = 0; x < phaseMesh.limit_x; x++)
			for (int y = 0; y < phaseMesh.limit_y; y++)
				for (int z = 0; z < phaseMesh.limit_z; z++) {
					PhaseNode& node = phaseMesh(x, y, z);
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						for (auto p_con = phase->x.begin(); p_con < phase->x.end(); p_con++) {
							statistics_node[phase->index].x[p_con->index].value += phase->phaseFraction * p_con->value;
						}
						for (auto p_chem = phase->potential.begin(); p_chem < phase->potential.end(); p_chem++) {
							statistics_node[phase->index].potential[p_chem->index].value += phase->phaseFraction * p_chem->value;
							statistics_node[phase->index].potential[p_chem->index].chemical_part += phase->phaseFraction * p_chem->chemical_part;
							statistics_node[phase->index].potential[p_chem->index].elastic_part += phase->phaseFraction * p_chem->elastic_part;
						}
						statistics_node[phase->index].phaseFraction += phase->phaseFraction;
						statistics_node[phase->index].chemEnergyDensity += phase->phaseFraction * phase->chemEnergyDensity;
						statistics_node[phase->index].elasEnergyDensity += phase->phaseFraction * phase->elasEnergyDensity;
					}
					for (auto comp = node.x.begin(); comp < node.x.end(); comp++) {
						statistics_node.x[comp->index].value += comp->value;
						statistics_node.x[comp->index].ChemicalReactionFlux += comp->ChemicalReactionFlux;
						statistics_node.x[comp->index].DiffusionFlux += comp->DiffusionFlux;
					}
					for (auto pot = node.potential.begin(); pot < node.potential.end(); pot++) {
						statistics_node.potential[pot->index].value += pot->value;
					}
				}

		for (auto p = statistics_node.begin(); p < statistics_node.end(); p++) {
			if (p->phaseFraction > Simulation_Num_Cut_Off) {
				p->chemEnergyDensity /= p->phaseFraction;
				p->elasEnergyDensity /= p->phaseFraction;
				for (auto p_con = p->x.begin(); p_con < p->x.end(); p_con++) {
					p_con->value = p_con->value / p->phaseFraction;
					if (std::isnan(p_con->value))
						dynamicCollection.is_write_dataFiles = false;
				}
				for (auto p_chem = p->potential.begin(); p_chem < p->potential.end(); p_chem++) {
					p_chem->value /= p->phaseFraction;
					p_chem->chemical_part /= p->phaseFraction;
					p_chem->elastic_part /= p->phaseFraction;
				}
			}
			else {
				p->chemEnergyDensity = 0.0;
				p->elasEnergyDensity = 0.0;
				for (auto x = p->x.begin(); x < p->x.end(); x++)
					x->value = 0.0;
				for (auto p_chem = p->potential.begin(); p_chem < p->potential.end(); p_chem++) {
					p_chem->value = 0.0;
					p_chem->chemical_part = 0.0;
					p_chem->elastic_part = 0.0;
				}
			}

			p->phaseFraction = p->phaseFraction / node_number;
			if (std::isnan(p->phaseFraction))
				dynamicCollection.is_write_dataFiles = false;
		}
		for (auto c = statistics_node.x.begin(); c < statistics_node.x.end(); c++) {
			c->value = c->value / node_number;
			c->ChemicalReactionFlux = c->ChemicalReactionFlux / node_number;
			c->DiffusionFlux = c->DiffusionFlux / node_number;
		}
		for (auto pot = statistics_node.potential.begin(); pot < statistics_node.potential.end(); pot++) {
			pot->value = pot->value / node_number;
		}
	}

	inline void Information::normalize_phi_in_mesh(FieldStorage_forPhaseNode& phaseMesh) {
		for (auto node = phaseMesh._mesh.begin(); node < phaseMesh._mesh.end(); node++) {
			double sum_phi = 0.0;
			for (auto phase = node->begin(); phase < node->end(); phase++)
				sum_phi += phase->phaseFraction;
			for (auto phase = node->begin(); phase < node->end(); phase++)
				phase->phaseFraction /= sum_phi;
		}
	}

	inline void Information::statistics_information_phaseFraction_in_phaseMesh(FieldStorage_forPhaseNode& phaseMesh) {
		double node_number = double(phaseMesh.limit_x * phaseMesh.limit_y * phaseMesh.limit_z);
		PhaseNode& statistics_node = dynamicCollection.inf_node;
		statistics_node = phaseMesh(0, 0, 0);
		for (auto p = statistics_node.begin(); p < statistics_node.end(); p++)
			p->phaseFraction = 0.0;
		for (auto node = phaseMesh._mesh.begin(); node < phaseMesh._mesh.end(); node++)
			for (auto phase = node->begin(); phase < node->end(); phase++)
				statistics_node[phase->index].phaseFraction += phase->phaseFraction;
		for (auto p = statistics_node.begin(); p < statistics_node.end(); p++)
			p->phaseFraction = p->phaseFraction / node_number;
	}

	inline void Information::generate_voronoi_structure(Vector3 box_position, Vector3 box_size, int grain_0_index, int grain_number, int generate_step,
		vector<int> phases_properties, vector<double> phases_weight, XNode x, double temperature) {
		Dimension dimention = Dimension::Three_Dimension;
		if (box_size[0] == 0 && box_size[1] == 0 && box_size[2] == 0)
			return;
		else if(box_size[0] == 0 || box_size[1] == 0 || box_size[2] == 0)
			dimention = Dimension::Two_Dimension;

		if (phases_properties.size() == 0) {
			phases_properties.push_back(materialSystem.matrix_phase.phaseProperty);
			phases_weight.push_back(1.0);
		}
		else if (phases_properties.size() != phases_weight.size()) {
			phases_weight.clear();
			for (unsigned int i = 0; i < phases_properties.size(); i++)
				phases_weight.push_back(1.0 / phases_properties.size());
		}
		// > normalize weight
		double sum_weight = 0.0;
		for (auto dd = phases_weight.begin(); dd < phases_weight.end(); dd++)
			sum_weight += *dd;
		for (auto dd = phases_weight.begin(); dd < phases_weight.end(); dd++)
			*dd = *dd / sum_weight;
		// > generate points
		RAND_init;
		vector<Point> points;
		for (int grain = 0; grain < grain_number; grain++) {
			double rand_x = RAND_0_1, rand_y = RAND_0_1, rand_z = RAND_0_1;
			Point p(rand_x * box_size[0] + box_position[0], rand_y * box_size[1] + box_position[1], rand_z * box_size[2] + box_position[2]);
			points.push_back(p);
			cout << "Voronoi: generate point at : (x, y, z) (" << p.x << ", " << p.y << ", " << p.x << ")" << endl;
		}
		// > points index and property
		vector<int> grains_property; // = points.size()
		for (auto grain = points.begin(); grain < points.end(); grain++) {
			double rand = RAND_0_1, sum_weight2 = 0.0;
			for (unsigned int ii = 0; ii < phases_weight.size(); ii++) {
				sum_weight2 += phases_weight[ii];
				if (sum_weight2 > rand) {
					grains_property.push_back(phases_properties[ii]);
					break;
				}
			}
		}
		// > periodic boundary condition
		vector<vector<Point>> mirror_points; // = 27 * points.size()
		int region_number = 0;
		if (dimention == Dimension::Three_Dimension)
			region_number = 27;
		else
			region_number = 9;
		mirror_points.resize(region_number);
		for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
			region->resize(grain_number);
		for (int grain = 0; grain < grain_number; grain++) {
			mirror_points[0][grain] = points[grain] + Point(0, 0, 0);
			mirror_points[1][grain] = points[grain] + Point(int(box_size[0]), 0, 0);
			mirror_points[2][grain] = points[grain] + Point(int(-box_size[0]), 0, 0);
			mirror_points[3][grain] = points[grain] + Point(0, int(box_size[1]), 0);
			mirror_points[4][grain] = points[grain] + Point(0, int(-box_size[1]), 0);
			mirror_points[5][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), 0);
			mirror_points[6][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), 0);
			mirror_points[7][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), 0);
			mirror_points[8][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), 0);
			if (dimention == Dimension::Three_Dimension) {
				mirror_points[9][grain] = points[grain] + Point(0, 0, int(box_size[2]));
				mirror_points[10][grain] = points[grain] + Point(0, 0, int(-box_size[2]));
				mirror_points[11][grain] = points[grain] + Point(int(box_size[0]), 0, int(box_size[2]));
				mirror_points[12][grain] = points[grain] + Point(int(-box_size[0]), 0, int(box_size[2]));
				mirror_points[13][grain] = points[grain] + Point(int(box_size[0]), 0, int(-box_size[2]));
				mirror_points[14][grain] = points[grain] + Point(int(-box_size[0]), 0, int(-box_size[2]));
				mirror_points[15][grain] = points[grain] + Point(0, int(box_size[1]), int(box_size[2]));
				mirror_points[16][grain] = points[grain] + Point(0, int(-box_size[1]), int(box_size[2]));
				mirror_points[17][grain] = points[grain] + Point(0, int(box_size[1]), int(-box_size[2]));
				mirror_points[18][grain] = points[grain] + Point(0, int(-box_size[1]), int(-box_size[2]));
				mirror_points[19][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), int(box_size[2]));
				mirror_points[20][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), int(box_size[2]));
				mirror_points[23][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), int(box_size[2]));
				mirror_points[24][grain] = points[grain] + Point(int(-box_size[0]), int(box_size[1]), int(-box_size[2]));
				mirror_points[26][grain] = points[grain] + Point(int(-box_size[0]), int(-box_size[1]), int(-box_size[2]));
				mirror_points[21][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), int(box_size[2]));
				mirror_points[22][grain] = points[grain] + Point(int(box_size[0]), int(box_size[1]), int(-box_size[2]));
				mirror_points[25][grain] = points[grain] + Point(int(box_size[0]), int(-box_size[1]), int(-box_size[2]));
			}
		}
		for (int region = 0; region < region_number; region++)
			for (int grain = 0; grain < grain_number; grain++) {
				Polyhedron poly(mirror_points[region][grain]);
				vector<point_in_region_index> record_points;
				point_in_region_index rp(region, grain);
				record_points.push_back(rp);
				for (unsigned int region_index = 0; region_index < mirror_points.size(); region_index++)
					for (unsigned int grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
						// Avoid inclusion points
						bool is_point_contained = false;
						for (auto re = record_points.begin(); re < record_points.end(); re++)
							if (re->region == region_index && re->grain_index == grain_index)
								is_point_contained = true;
						if (is_point_contained)
							continue;
						// prepare
						Point norm, mid_point;
						norm = poly.point_inside_polyhedron - mirror_points[region_index][grain_index];
						mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
						// Judged ipsilateral to poly center
						if (poly.check_point_inside_polyhedron(mid_point) == false)
							continue;
						// < add point
						poly.add_surf(norm, mid_point);
						// < Eliminate meaningless points in poly
						for (auto re = record_points.begin(); re < record_points.end();) {
							Point check = (points[grain] + mirror_points[re->region][re->grain_index]) / 2;
							if (poly.check_point_inside_polyhedron(check) == false) {
								re = record_points.erase(re);
							}
							else {
								++re;
							}
						}
					}
				GeometricRegion geo;
				geo.geometryProperty = Geometry::Geo_Polyhedron;
				geo.generate_step = generate_step;
				geo.polyhedron = poly;
				geo.phaseIndex = grain_0_index + grain;
				geo.phaseProperty = grains_property[grain];
				geo.temperature = temperature;
				geo.x = x;
				nucleationBox.geometricRegion_box.push_back(geo);
			}

	}

#ifdef _WIN32
	inline void Information::generate_structure_from_BMP_pic(int phaseIndex, int phase_property, int generate_step, XNode x, double threshold[2], double temperature, string fileName) {
		if (settings.disperse_settings.Nz > 1) {
			cout << "BMP init structure: Dimention is 3D, cant init structure by this methode !" << endl;
			return;
		}

		BMP24reader bmpReader;
		if (fileName == "") {
#ifdef _WIN32
			SelectFolderPath(fileName);
#else
			std::cout << "> Error, init datafile failed write filename !" << endl;
			SYS_PROGRAM_STOP;
#endif
		}
		bmpReader.safe(fileName);
		bmpReader.read(fileName);
		PointSet set;
		for (int y = 0; y < settings.disperse_settings.Ny; y++)
			for (int x = 0; x < settings.disperse_settings.Nx; x++) {
				int xx = double_to_int(x * bmpReader.bmp_width / settings.disperse_settings.Nx), yy = double_to_int(y * bmpReader.bmp_height / settings.disperse_settings.Ny);
				double graypercent = bmpReader.getGrayPercentage(xx, yy);
				if (graypercent > threshold[0] && graypercent < threshold[1])
					set.add_point(x, y, 0);
			}
		set.generate_step = generate_step;
		set.phaseIndex = phaseIndex;
		set.phaseProperty = phase_property;
		set.temperature = temperature;
		set.x = x;
		nucleationBox.pointSet_box.push_back(set);
	}
#endif
}