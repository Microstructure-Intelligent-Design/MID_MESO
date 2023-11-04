#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"

namespace pf {
	enum StraggeredGridsFlag{Minus_0_5, Plus_0_5};
	namespace fluid_boundary_condition_funcs {
		static void boundary_condition_for_domain_U(VectorNode& u, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			return;
		}
		static void boundary_condition_for_domain_V(VectorNode& v, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			return;
		}
		static void boundary_condition_for_domain_W(VectorNode& w, PhaseNode& down_node, PhaseNode& up_node, int Nx, int Ny, int Nz) {
			return;
		}
		static void boundary_condition_for_main_domain(PhaseNode& node, int Nx, int Ny, int Nz) {
			return;
		}
	};
	class FluidField
	{
	public:

		FluidField() {};
		~FluidField() { free(); };
		void init(FieldStorage_forPhaseNode& _simulationField, Information& _inf) {
			information = &_inf;
			phaseMesh = &_simulationField;
			if (!information->materialSystem.is_fluid_on)
				return;
			boundary_condition_for_domain_U = fluid_boundary_condition_funcs::boundary_condition_for_domain_U;
			boundary_condition_for_domain_V = fluid_boundary_condition_funcs::boundary_condition_for_domain_V;
			boundary_condition_for_domain_W = fluid_boundary_condition_funcs::boundary_condition_for_domain_W;
			boundary_condition_for_main_domain = fluid_boundary_condition_funcs::boundary_condition_for_main_domain;
			if (information->settings.disperse_settings.x_bc == PERIODIC) {
				U.init(information->settings.disperse_settings.Nx, information->settings.disperse_settings.Ny, information->settings.disperse_settings.Nz,
					information->settings.disperse_settings.dx, information->settings.disperse_settings.x_bc, information->settings.disperse_settings.y_bc,
					information->settings.disperse_settings.z_bc);
			}
			else {
				U.init(information->settings.disperse_settings.Nx + 1, information->settings.disperse_settings.Ny, information->settings.disperse_settings.Nz,
					information->settings.disperse_settings.dx, information->settings.disperse_settings.x_bc, information->settings.disperse_settings.y_bc,
					information->settings.disperse_settings.z_bc);
			}
			if (information->settings.disperse_settings.y_bc == PERIODIC) {
				V.init(information->settings.disperse_settings.Nx, information->settings.disperse_settings.Ny, information->settings.disperse_settings.Nz,
					information->settings.disperse_settings.dx, information->settings.disperse_settings.x_bc, information->settings.disperse_settings.y_bc,
					information->settings.disperse_settings.z_bc);
			}
			else {
				V.init(information->settings.disperse_settings.Nx, information->settings.disperse_settings.Ny + 1, information->settings.disperse_settings.Nz,
					information->settings.disperse_settings.dx, information->settings.disperse_settings.x_bc, information->settings.disperse_settings.y_bc,
					information->settings.disperse_settings.z_bc);
			}
			if (information->settings.disperse_settings.z_bc == PERIODIC) {
				W.init(information->settings.disperse_settings.Nx, information->settings.disperse_settings.Ny, information->settings.disperse_settings.Nz,
					information->settings.disperse_settings.dx, information->settings.disperse_settings.x_bc, information->settings.disperse_settings.y_bc,
					information->settings.disperse_settings.z_bc);
			}
			else {
				W.init(information->settings.disperse_settings.Nx, information->settings.disperse_settings.Ny, information->settings.disperse_settings.Nz + 1,
					information->settings.disperse_settings.dx, information->settings.disperse_settings.x_bc, information->settings.disperse_settings.y_bc,
					information->settings.disperse_settings.z_bc);
			}
			for (auto vecNode = U._mesh.begin(); vecNode < U._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value
				vecNode->vals.push_back(0.0); // increment
			}
			for (auto vecNode = V._mesh.begin(); vecNode < V._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value
				vecNode->vals.push_back(0.0); // increment
			}
			for (auto vecNode = W._mesh.begin(); vecNode < W._mesh.end(); vecNode++) {
				vecNode->vals.push_back(0.0); // value
				vecNode->vals.push_back(0.0); // increment
			}
		}

		void init_pressure_in_fluid(double pressure = 0.0);
		void do_boundary_condition();
		// return [ MAX_abs_DISPERSION, dU/dt, dV/dt, dW/dt ] in information
		void evolve_momentum_equation(double fluid_dt);
		//return MAX_abs_dPressure
		double do_pressure_correction(double accuracy = 1e-6
			, double relaxation_factor = 1.0, bool debug_solver = false, int output_step = 1000); 
		double do_pressure_correction_with_average_boundary_condition(double ave_pressure, double accuracy = 1e-6
				, double relaxation_factor = 1.0, bool debug_solver = false, int output_step = 1000);
		void correcting_velocity_field(double fluid_dt);
		void assign_velocity_to_main_domain();
		void assign_velocity_from_main_domain();

		void init_u_in_velocity_field(int x, int y, int z, double u);
		void init_v_in_velocity_field(int x, int y, int z, double v);
		void init_w_in_velocity_field(int x, int y, int z, double w);
		void set_boundary_condition_for_domain_U(void(*boundary)(pf::VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int)) {
			boundary_condition_for_domain_U = boundary;
		}
		void set_boundary_condition_for_domain_V(void(*boundary)(pf::VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int)) {
			boundary_condition_for_domain_V = boundary;
		}
		void set_boundary_condition_for_domain_W(void(*boundary)(pf::VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int)) {
			boundary_condition_for_domain_W = boundary;
		}
		void set_boundary_condition_for_main_domain(void(*boundary)(pf::PhaseNode&, int, int, int)) {
			boundary_condition_for_main_domain = boundary;
		}
		void free() {
			information = nullptr;
			phaseMesh = nullptr;
			boundary_condition_for_domain_U = nullptr;
			boundary_condition_for_domain_V = nullptr;
			boundary_condition_for_domain_W = nullptr;
			boundary_condition_for_main_domain = nullptr;
			U.free();
			V.free();
			W.free();
		}
	private:
		void(*boundary_condition_for_domain_U)(VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int);
		void(*boundary_condition_for_domain_V)(VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int);
		void(*boundary_condition_for_domain_W)(VectorNode&, pf::PhaseNode&, pf::PhaseNode&, int, int, int);
		void(*boundary_condition_for_main_domain)(PhaseNode&, int, int, int);
		Information* information;
		FieldStorage_forPhaseNode* phaseMesh;
		FieldStorage_forVector U;
		FieldStorage_forVector V;
		FieldStorage_forVector W;
	};
}