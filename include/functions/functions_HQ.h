#pragma once
#include"../MPF.h"
namespace pf {
	namespace functions_HQ {

		static void change_phaseproperty_of_grain(FieldStorage_forPhaseNode& mesh, int phaseIndex, int phaseProperty) {
			for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++)
				for (auto phase = node->begin(); phase < node->end(); phase++)
					if (phase->index == phaseIndex)
						phase->phaseProperty = phaseProperty;
		}

		static void change_x_of_grain(FieldStorage_forPhaseNode& mesh, int phaseIndex, XNode x) {
			for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++)
				for (auto phase = node->begin(); phase < node->end(); phase++)
					if (phase->index == phaseIndex) {
						phase->x = x;
						phase->potential = x;
					}
			
		}

		static void change_sublattice_structure_of_grain(FieldStorage_forPhaseNode& mesh, int phaseIndex, SublatticeNode sublattice) {
			for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++)
				for (auto phase = node->begin(); phase < node->end(); phase++)
					if (phase->index == phaseIndex)
						phase->sublattice = sublattice;
		}

		// return scale
		static void automatically_adjust_dt(double& dt_scale, double MAX_VARIATION, double LIMIT_VARIATION, int istep, int start_treatment_step = 0, int delt_step = 100, double MAX_SCALE = 1e3, double MIN_Mob_SCALE = 1e-3, bool isReduceOutput = true) {
			stringstream log;
			double delt_scale = 0.05;
			if (istep > start_treatment_step) {
				if (MAX_VARIATION > LIMIT_VARIATION) {
					dt_scale /= MAX_VARIATION / LIMIT_VARIATION + delt_scale;
					if (!isReduceOutput)
						log << "> Adjust time interval for evolving stability, scale = " << dt_scale << endl;
				}
				else if (istep % delt_step == 0 && (dt_scale * 1.1) <= MAX_SCALE) {
					dt_scale *= 1.1;
					log << "> Adjust time interval for fields evolving quickly, scale = " << dt_scale << endl;
				}
			}
			cout << log.str();
		}

		static void free_nonexistence_phase(FieldStorage_forPhaseNode& mesh) {
			Information inf;
			inf.statistics_information_phaseFraction_in_phaseMesh(mesh);
			for (auto p = inf.dynamicCollection.inf_node.begin(); p < inf.dynamicCollection.inf_node.end(); p++)
				if (p->phaseFraction < SYS_EPSILON) {
					for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++)
						node->erase((*node)[p->index]);
				}
		}

		static void merge_two_phases(FieldStorage_forPhaseNode& mesh, int exist_phase_index, int eliminate_phase_index) {
			for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++) {
				(*node)[exist_phase_index].phaseFraction += (*node)[eliminate_phase_index].phaseFraction;
				(*node)[eliminate_phase_index].phaseFraction = 0.0;
				(*node).automatic_set_flag();
			}
			for (auto node = mesh._mesh.begin(); node < mesh._mesh.end(); node++)
				node->erase((*node)[eliminate_phase_index]);
		}

		static vector<string> cal_gibbs_free_energy(pf::Information inf, int phase_property, vector<int> comps, vector<vector<double>> comps_range, double precision = 0.01) {
			if (comps.size() == 1) {
				cout << "only one component, dont need to use this function !" << endl;
				exit(0);
			}
			Thermodynamics thermodynamic;
			FieldStorage_forPhaseNode mesh;
			thermodynamic.init(mesh, inf);
			pf::PhaseNode node;
			node.add_phase(0, phase_property, 1, 1.0);
			node.tempValues.temperature = inf.materialSystem.init_temperature;
			vector<double> individual_comps_value;
			vector<string> strs;
			for (int index = 0; index < comps.size() - 1; index++) {
				individual_comps_value.push_back(comps_range[index][0]);
				node[0].x.add_con(comps[index], 0.0);
				node[0].potential.add_con(comps[index], 0.0);
				string str = "C_" + inf.materialSystem.phases[phase_property].phase_name + "_" + inf.materialSystem.sys_x[comps[index]].name + " = [";
				strs.push_back(str);
				str = "dG_" + inf.materialSystem.phases[phase_property].phase_name + "_dC_" + inf.materialSystem.sys_x[comps[index]].name + " = [";
				strs.push_back(str);
			}
			node[0].x.add_con(*(comps.end() - 1), 0.0);
			node[0].potential.add_con(*(comps.end() - 1), 0.0);
			string str = "G_" + inf.materialSystem.phases[phase_property].phase_name + " = [";
			strs.push_back(str);

			double vm = inf.materialSystem.functions.MolarVolume(node, node[0], inf.dynamicCollection);
			bool iterate = true;
			bool seperator = false;
			long times = 0;
			while (iterate)
			{
				times++;
				cout << "iterate times = " << times << endl;
				seperator = false;
				for (int index = 0; index < individual_comps_value.size(); index++) {
					if (index < (individual_comps_value.size() - 1)) {
						if (individual_comps_value[index] >= comps_range[index][1] + SYS_EPSILON) {
							individual_comps_value[index] = comps_range[index][0];
							individual_comps_value[index + 1] += precision;
						}
						else if(individual_comps_value[index] + precision >= comps_range[index][1] + SYS_EPSILON)
							seperator = true;
					}
					else if (index == (individual_comps_value.size() - 1) && individual_comps_value[index] >= comps_range[index][1] + SYS_EPSILON) {
						iterate = false;
					}
				}
				if (iterate) {
					//----------------------------------------------------------------------------------------
					double solvent = 1.0;
					for (auto val = individual_comps_value.begin(); val < individual_comps_value.end(); val++)
						solvent -= *val;
					//-----------------------------------------------------------------------------------------
					int cindex = 0;
					bool cal_thermo = true;
					for (auto c = comps.begin(); c < comps.end() - 1; c++) {
						node[0].x[*c].value = individual_comps_value[cindex];
						if (individual_comps_value[cindex] <= comps_range[cindex][0] - SYS_EPSILON || individual_comps_value[cindex] >= comps_range[cindex][1] + SYS_EPSILON)
							cal_thermo = false;
						cindex++;
					}
					node[0].x[*(comps.end() - 1)].value = solvent;
					if (solvent <= comps_range[cindex][0] || solvent >= comps_range[cindex][1])
						cal_thermo = false;
					//-----------------------------------------------------------------------------------------
					if (cal_thermo) {

						Funcs& funcs = inf.materialSystem.functions;
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							if (phase->phaseFraction > Simulation_Num_Cut_Off) {
								thermodynamic.pretreatment_thermodynamic_data_of_onePhase(node, phase->index);
								funcs.Energy(node, *phase, inf.dynamicCollection);
								funcs.Potential(node, *phase, inf.dynamicCollection);
								thermodynamic.aftertreatment_thermodynamic_data_of_onePhase(node, phase->index);
							}
						}

						int cindex = 0;
						for (auto c = comps.begin(); c < comps.end() - 1; c++) {
							if (seperator) {
								strs[2 * cindex] += to_string(node[0].x[*c].value) + " ; ";
								strs[2 * cindex + 1] += to_string(node[0].potential[*c].chemical_part * vm) + " ; ";
							}
							else {
								strs[2 * cindex] += to_string(node[0].x[*c].value) + " , ";
								strs[2 * cindex + 1] += to_string(node[0].potential[*c].chemical_part * vm) + " , ";
							}
							cindex++;
						}
						if (seperator) {
							strs[2 * cindex] += to_string(node[0].chemEnergyDensity * vm) + " ; ";
						}
						else {
							strs[2 * cindex] += to_string(node[0].chemEnergyDensity * vm) + " , ";
						}
					}
					else {
						int cindex = 0;
						for (auto c = comps.begin(); c < comps.end() - 1; c++) {
							if (seperator) {
								strs[2 * cindex] += to_string(node[0].x[*c].value) + " ; ";
								strs[2 * cindex + 1] += "NaN ; ";
							}
							else {
								strs[2 * cindex] += to_string(node[0].x[*c].value) + " , ";
								strs[2 * cindex + 1] += "NaN , ";
							}
							cindex++;
						}
						if (seperator) {
							strs[2 * cindex] += "NaN ; ";
						}
						else {
							strs[2 * cindex] += "NaN , ";
						}
					}
				}
				//-----------------------------------------------------------------------------------------
				individual_comps_value[0] += precision;
			}
			int cindex = 0;
			for (auto c = comps.begin(); c < comps.end() - 1; c++) {
				strs[2 * cindex] += " ];\n";
				strs[2 * cindex + 1] += " ];\n";
				cindex++;
			}
			strs[2 * cindex] += " ];\n";
			return strs;
		}

		static void OpenMP_test() {
#pragma omp parallel
			printf("OpenMP_test! Thread N.%d\n", omp_get_thread_num());
			return;
		}

		static void output_in_loop(pf::Information information, FieldStorage_forPhaseNode& phaseMesh, int istep, string custom_log = "") {
			int z_direction = 40;
			if (istep % information.settings.file_settings.file_output_step == 0) {
				information.statistics_information_in_phaseMesh(phaseMesh);
				PhaseNode& inf_node = information.dynamicCollection.inf_node;
				Info_Node& set_comps = information.materialSystem.sys_x;
				Info_Phases& set_phases = information.materialSystem.phases;
				double NaN = information.writer.NaN;
				// write scalar file
				ofstream fout;
				// header
				{
					string fname;
					fname = information.settings.file_settings.working_folder_path + dirSeparator + "scalar_variables_step" + to_string(istep) + ".vts";
					fout.open(fname);
					if (!fout) {
						cout << "Failed to write the vtk file..." << endl;
						fout.close();
						return;
					}
					fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
					fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
					fout << "<StructuredGrid WholeExtent=\""
						<< 0 << " " << phaseMesh.limit_x - 1 << " "
						<< 0 << " " << phaseMesh.limit_y - 1 << " "
						<< 0 << " " << z_direction - 1 << "\"> " << endl;
					fout << "<PointData Scalars= \"ScalarData\">" << endl;
				}
				// data
				{
					fout << "<DataArray type = \"Float64\" Name = \"" << "grains" <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < z_direction; k++)
						for (int j = 0; j < phaseMesh.limit_y; j++)
							for (int i = 0; i < phaseMesh.limit_x; i++) {
								double fix = 0.0;
								for (auto phase = phaseMesh(i, j, 0).begin(); phase < phaseMesh(i, j, 0).end(); phase++)
									fix += phase->phaseFraction * phase->phaseFraction;
								if (IS_NAN(fix))
									fout << NaN << endl;
								else
									fout << 1.0 - fix << endl;
							}
					fout << "</DataArray>" << endl;
					for (auto phase = inf_node.begin(); phase < inf_node.end(); phase++) {
						string name;
						name = "\"Phi" + to_string(phase->index) + "_" + set_phases[phase->phaseProperty].phase_name + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << name <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < z_direction; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i) {
									double phi = phaseMesh(i, j, 0)[phase->index].phaseFraction;
									if (IS_NAN(phi))
										fout << NaN << endl;
									else
										fout << phi << endl;
								}
						fout << "</DataArray>" << endl;
					}
					for (auto comp = set_comps.begin(); comp < set_comps.end(); comp++) {
						string name = "\"con_" + to_string(comp->index) + "_" + set_comps[comp->index].name + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << name <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < z_direction; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i) {
									double con = phaseMesh(i, j, 0).x[comp->index].value;
									if (IS_NAN(con))
										fout << NaN << endl;
									else
										fout << con << endl;
								}
						fout << "</DataArray>" << endl;
					}
					for (auto comp = set_comps.begin(); comp < set_comps.end(); comp++) {
						string name = "\"potential_" + to_string(comp->index) + "_" + set_comps[comp->index].name + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << name <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < z_direction; k++)
							for (int j = 0; j < phaseMesh.limit_y; j++)
								for (int i = 0; i < phaseMesh.limit_x; i++) {
									fout << phaseMesh(i, j, 0).potential[comp->index].value << endl;
								}
						fout << "</DataArray>" << endl;
					}
					vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
					/*for (int ele_index = 0; ele_index < 6; ele_index++)
					{
						string compname = "\"stress_" + compNameV[ele_index] + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << compname <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < z_direction; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i)
								{
									fout << phaseMesh(i, j, 0).mechanicalValues.Stresses[ele_index] << endl;
								}
						fout << "</DataArray>" << endl;
					}*/
					/*for (int ele_index = 0; ele_index < 6; ele_index++)
					{
						string compname = "\"strain_" + compNameV[ele_index] + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << compname <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < z_direction; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i)
								{
									fout << phaseMesh(i, j, 0).mechanicalValues.Strains[ele_index] << endl;
								}
						fout << "</DataArray>" << endl;
					}*/
					/*for (int ele_index = 0; ele_index < 6; ele_index++)
					{
						string compname = "\"eigenStrain_" + compNameV[ele_index] + "\" ";
						fout << "<DataArray type = \"Float64\" Name = " << compname <<
							"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
						for (int k = 0; k < z_direction; ++k)
							for (int j = 0; j < phaseMesh.limit_y; ++j)
								for (int i = 0; i < phaseMesh.limit_x; ++i)
								{
									fout << phaseMesh(i, j, 0).mechanicalValues.EffectiveEigenStrain[ele_index] << endl;
								}
						fout << "</DataArray>" << endl;
					}*/
					fout << "<DataArray type = \"Float64\" Name = \"" << "J1" <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < z_direction; ++k)
						for (int j = 0; j < phaseMesh.limit_y; ++j)
							for (int i = 0; i < phaseMesh.limit_x; ++i)
							{
								fout << phaseMesh(i, j, 0).mechanicalValues.Stresses.J1() << endl;
							}
					fout << "</DataArray>" << endl;

					fout << "<DataArray type = \"Float64\" Name = \"" << "vMises" <<
						"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
					for (int k = 0; k < z_direction; ++k)
						for (int j = 0; j < phaseMesh.limit_y; ++j)
							for (int i = 0; i < phaseMesh.limit_x; ++i)
							{
								fout << phaseMesh(i, j, 0).mechanicalValues.Stresses.Mises() << endl;
							}
					fout << "</DataArray>" << endl;
				}
				// tail
				{
				fout << "</PointData>" << endl;
				fout << "<Points>" << endl;
				fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
				for (int k = 0; k < z_direction; ++k)
					for (int j = 0; j < phaseMesh.limit_y; ++j)
						for (int i = 0; i < phaseMesh.limit_x; ++i)
						{
							fout << i << " " << j << " " << k << "\n";
						}
				fout << "</DataArray>" << endl;
				fout << "</Points>" << endl;
				fout << "</StructuredGrid>" << endl;
				fout << "</VTKFile>" << endl;
				fout.close();
				}
				
				//write vector file
				/*ofstream fout2;
				information.writer.open_vts_vec3_file(fout2, phaseMesh, istep);
				if (information.settings.file_settings.isPhaseGradientOutput)
					information.writer.write_vec3_phiGradient(fout2, phaseMesh);
				if (information.settings.file_settings.isMagneticFieldOutput)
					information.writer.write_vec3_demagnetizingIntensity(fout2, phaseMesh);
				if (information.settings.file_settings.isFluidFieldOutput)
					information.writer.write_vec3_velocity(fout2, phaseMesh);
				information.materialSystem.functions.write_vec3_customization(fout, phaseMesh);
				information.writer.close_vts_file(fout2, phaseMesh);*/
			}

			if (istep % information.settings.file_settings.screen_output_step == 0) {
				information.statistics_information_in_phaseMesh(phaseMesh);
				//////////////////////output to screen/////////////////////////
				stringstream log;
				PhaseNode& inf_node = information.dynamicCollection.inf_node;
				Info_Node& comp_inf = information.materialSystem.sys_x;
				log << "#########################################################################" << endl;
				timer::unlock_class();
				log << timer::return_cunrrent_time_by_string();
				timer::lock_class();

				log << "# > Simulation step_" << istep << " has been finished!" << endl;
				log << "# (Program Run) This " << information.settings.file_settings.screen_output_step << " steps used " << setprecision(4) << timer::interval_end() << "(s), Total " << istep << " steps used " << timer::total_duration_sec() << "(s)." << endl;
				log << "# Simulation time is " << information.dynamicCollection.real_time << " (sec.)" << endl;
				log << "# (One Step) ConcentrationField has been Averaged for " << information.dynamicCollection.average_interface_phaseCon_num << " times!" << endl;
				log << "# (One Step) ElectricField has been iterated for " << information.dynamicCollection.electricField__iterate_counts_once << " times!" << endl;
				log << "# (One Step) MagneticField has been iterated for " << information.dynamicCollection.magneticField__iterate_counts_once << " times!" << endl;
				log << "# (One Step) ElasticField has been iterated for " << information.dynamicCollection.elastic_iterate_counts_once << " times!" << endl;
				log << "# (One Step) FluidField has been iterated for " << information.dynamicCollection.fluidField_pressure__iterate_counts_once << " times!" << endl;
				log << "#================= The percentage of phase and component =================" << endl;
				for (auto p = inf_node.begin(); p < inf_node.end(); p++) {
					log << "#  Phase_" << p->index << "(" << information.materialSystem.phases[p->phaseProperty].phase_name << ")" << " : " << setprecision(5) << p->phaseFraction * 100 << " %" << endl;
					log << "   (average) Chem. energy density : " << setprecision(5) << p->chemEnergyDensity << " , Elas. energy density : " << setprecision(5) << p->elasEnergyDensity << " (J/m^3)" << endl;

					//double matter = 0.0;
					for (auto c = p->x.begin(); c < p->x.end(); c++) {
						log << "#   Con_" << comp_inf[c->index].name << " : " << c->value << endl;
						log << "    Potential : " << p->potential[c->index].value << ", Chem. pot. : " << p->potential[c->index].chemical_part << ", Elas. pot. : " << p->potential[c->index].elastic_part << " (J/m^3)" << endl;
						//matter += c->value;
					}
					//log << "#   Matter  : " << matter << "(mol/m^3)  " <<  endl;
					log << "--------------------------------------------------------------------------" << endl;
				}
				if (comp_inf.size() != 0) {
					log << "#---------------------- Total Component moleFraction ----------------------" << endl;
					for (auto c = comp_inf.begin(); c < comp_inf.end(); c++)
					{
						log << "#  Total Con_" << c->name << " : " << setprecision(8) << inf_node.x[c->index].value << endl;
					}
				}
				if (information.materialSystem.is_fluid_on) {
					log << "#============================== Fluid Field =============================" << endl;
					log << "# MAX abs Dispersion = " << setprecision(5) << information.dynamicCollection.MAX_ABS_DISPERSION << endl;
					log << "# MAX abs dUdt = " << setprecision(5) << information.dynamicCollection.MAX_ABS_dUdt << endl;
					log << "# MAX abs dVdt = " << setprecision(5) << information.dynamicCollection.MAX_ABS_dVdt << endl;
					log << "# MAX abs dWdt = " << setprecision(5) << information.dynamicCollection.MAX_ABS_dWdt << endl;
					log << "# MAX abs dPressure = " << setprecision(5) << information.dynamicCollection.MAX_ABS_dPRESSURE << endl;
				}
				log << "#====================== Phase & Con & Temp Field ========================" << endl;
				log << "# MAX Phase Increment = " << setprecision(5) << information.dynamicCollection.MAX_phase_increment << endl;
				log << "# MAX Flux : diffusion = " << setprecision(5) << information.dynamicCollection.MAX_diffusionFlux
					<< ", reaction = " << setprecision(5) << information.dynamicCollection.MAX_reactionFlux
					<< ", phaseTrans = " << setprecision(5) << information.dynamicCollection.MAX_phaseTransFlux << endl;
				log << "# MAX Temperature Increment = " << setprecision(5) << information.dynamicCollection.MAX_temp_increment << endl;
				log << "#============================== Custom Log ==============================" << endl;
				log << custom_log << endl;
				log << "#============================= D E B U G ================================" << endl;
				if (information.dynamicCollection.is_DiffusionFlux_term_out_of_limit)
					log << "# Problem : Kinetics model: Component change (due to diffusion) too quikily." << endl;
				if (information.dynamicCollection.is_PhaseTransitionFlux_term_out_of_limit)
					log << "# Problem : Kinetics model: Component change (due to phase transition) too quikily." << endl;
				if (information.dynamicCollection.is_ChemicalReactionFlux_term_out_of_limit)
					log << "# Problem : Kinetics model: Component change (due to source term) too quikily." << endl;
				if (information.dynamicCollection.is_average_line_too_small)
					log << "# Problem : Kinetics model: Average lines too small, please lengthen." << endl;
				if (information.dynamicCollection.is_interface_move_too_quickly)
					log << "# Problem : Mobility model: Interface moved too quickly, please slow down." << endl;
				if (information.dynamicCollection.is_cal_constituent_of_onePhase_by_energy_minimization_too_quickly)
					log << "# Problem : Thermodynamics model: energy minimization too quickly, please check and limit." << endl;
				log << "#########################################################################" << endl;
				log << endl << endl;
				timer::interval_begin();
				std::cout << log.str();

				information.writer.add_string_to_txt(log.str(), information.settings.file_settings.log_file_name);

				information.dynamicCollection.init_each_outputStep();
			}
		}

		enum TotalConcentrationFlag{OUT_OF_INTERFACE, NEW_IN_INTERFACE, ITERATED};
		namespace DEFULT_FUNCTIONS_FOR_TOTAL_CONCENTRATION {
			static double Mij(pf::PhaseNode& node, int comp_i, int comp_j, pf::Info_DynamicCollection& inf) {
				return 0.0;
			}
			static double Si(pf::PhaseNode& node, int comp_i, pf::Info_DynamicCollection& inf) {
				return 0.0;
			}
			static void Boundary(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
				return;
			}
			static void potential_func(pf::PhaseNode& node, int flag_index, pf::Info_DynamicCollection& inf) {
				if (flag_index != TotalConcentrationFlag::OUT_OF_INTERFACE) {
					// init phase concentration
					double sum_phase_fraction = 0.0;
					if (flag_index == TotalConcentrationFlag::NEW_IN_INTERFACE) {
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							if (phase->_flag == pf_INTERFACE) {
								sum_phase_fraction += phase->phaseFraction;
								for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
									comp->value = node.x[comp->index].value;
							}
						}
					}
					else if (flag_index == TotalConcentrationFlag::ITERATED) {
						ConNode con = node.x;
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							if (phase->_flag == pf_INTERFACE) {
								sum_phase_fraction += phase->phaseFraction;
								for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
									con[comp->index].value -= comp->value * phase->phaseFraction;
							}
						}
						for (auto phase = node.begin(); phase < node.end(); phase++)
							if (phase->_flag == pf_INTERFACE)
								for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
									comp->value += con[comp->index].value;
					}
					// iterate calculation phase concentration
					// ......


					// calculate energy density of phase and potential
					// energy(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf)
					// potential(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf)
				}
				else {
					for (auto phase = node.begin(); phase < node.end(); phase++) {
						if(phase->phaseFraction > (1.0 - Simulation_Num_Cut_Off))
							for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
								comp->value = node.x[comp->index].value;
					}
					// calculate energy density of phase and potential 
					// energy(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf)
					// potential(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf)
				}
				return;
			}
		}
		class TotalConcentration
		{
		public:
			TotalConcentration(pf::FieldStorage_forPhaseNode& _mesh, pf::Information _inf, int _flag_index = 1000) {
				phaseMesh = &_mesh;
				information = &_inf;
				flag_index = _flag_index;
				Mij = DEFULT_FUNCTIONS_FOR_TOTAL_CONCENTRATION::Mij;
				Si = DEFULT_FUNCTIONS_FOR_TOTAL_CONCENTRATION::Si;
				potential_func = DEFULT_FUNCTIONS_FOR_TOTAL_CONCENTRATION::potential_func;
				Boundary = DEFULT_FUNCTIONS_FOR_TOTAL_CONCENTRATION::Boundary;
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++)
							(*phaseMesh)(x, y, z).customFlags.add_int(_flag_index, OUT_OF_INTERFACE);
			}
			~TotalConcentration() {
				phaseMesh = nullptr;
				information = nullptr;
				Mij = nullptr;
				Si = nullptr;
				potential_func = nullptr;
				Boundary = nullptr;
			}
			void define_Mij(double (*M_ij)(pf::PhaseNode&, int, int, pf::Info_DynamicCollection&)) {
				Mij = M_ij;
			}
			void define_Si(double (*S_i)(pf::PhaseNode&, int, pf::Info_DynamicCollection&)) {
				Si = S_i;
			}
			void define_Boundary(void (*_Boundary)(pf::PhaseNode&, pf::Info_DynamicCollection&)) {
				Boundary = _Boundary;
			}
			void define_potential_func(void (*Potential_func)(pf::PhaseNode&, int, pf::Info_DynamicCollection&)) {
				potential_func = Potential_func;
			}
			int Flag_index() {
				return flag_index;
			}
			void init_total_component_with_phase_component() {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++)
							(*phaseMesh)(x, y, z).cal_x_from_phase_x();
			}
			void solve_potential() {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) 
							potential_func((*phaseMesh)(x, y, z), flag_index, information->dynamicCollection);
			}

			double solve_Cahn_Hilliard(pf::PhaseFluxModel flux_DrivingForce = pf::PhaseFluxModel::IntDiff_PotentialGrad, pf::DifferenceMethod differenceMethod = pf::DifferenceMethod::FIVE_POINT,  bool adjust_c_0_1 = false) {
				double dx = phaseMesh->dx, dt = information->settings.disperse_settings.dt * information->dynamicCollection.dt_scale;
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							if (flux_DrivingForce == pf::PhaseFluxModel::IntDiff_ConGrad) {
								for (auto comp_i = node.x.begin(); comp_i < node.x.end(); comp_i++) {
									comp_i->gradient[0] = (node.get_neighbor_node(pf::Direction::x_down).x[comp_i->index].value - node.get_neighbor_node(pf::Direction::x_up).x[comp_i->index].value) / 2.0 / dx;
									comp_i->gradient[1] = (node.get_neighbor_node(pf::Direction::y_down).x[comp_i->index].value - node.get_neighbor_node(pf::Direction::y_up).x[comp_i->index].value) / 2.0 / dx;
									comp_i->gradient[2] = (node.get_neighbor_node(pf::Direction::z_down).x[comp_i->index].value - node.get_neighbor_node(pf::Direction::z_up).x[comp_i->index].value) / 2.0 / dx;
									if (differenceMethod == pf::DifferenceMethod::FIVE_POINT) {
										comp_i->laplacian = (node.get_neighbor_node(pf::Direction::x_down).x[comp_i->index].value + node.get_neighbor_node(pf::Direction::x_up).x[comp_i->index].value
											+ node.get_neighbor_node(pf::Direction::y_down).x[comp_i->index].value + node.get_neighbor_node(pf::Direction::y_up).x[comp_i->index].value
											+ node.get_neighbor_node(pf::Direction::z_down).x[comp_i->index].value + node.get_neighbor_node(pf::Direction::z_up).x[comp_i->index].value
											- 6.0 * comp_i->value) / dx / dx;
									}
									else if (differenceMethod == pf::DifferenceMethod::NINE_POINT) {
										comp_i->laplacian = (4.0 * node.get_neighbor_node(pf::Direction::x_down).x[comp_i->index].value + 4.0 * node.get_neighbor_node(pf::Direction::x_up).x[comp_i->index].value
											+ 4.0 * node.get_neighbor_node(pf::Direction::y_down).x[comp_i->index].value + 4.0 * node.get_neighbor_node(pf::Direction::y_up).x[comp_i->index].value
											+ 4.0 * node.get_neighbor_node(pf::Direction::z_down).x[comp_i->index].value + 4.0 * node.get_neighbor_node(pf::Direction::z_up).x[comp_i->index].value
											+ node.get_long_range_node(-1, -1, 0).x[comp_i->index].value + node.get_long_range_node(-1, 1, 0).x[comp_i->index].value
											+ node.get_long_range_node(1, -1, 0).x[comp_i->index].value + node.get_long_range_node(1, 1, 0).x[comp_i->index].value
											+ node.get_long_range_node(-1, 0, -1).x[comp_i->index].value + node.get_long_range_node(-1, 0, 1).x[comp_i->index].value
											+ node.get_long_range_node(1, 0, -1).x[comp_i->index].value + node.get_long_range_node(1, 0, 1).x[comp_i->index].value
											+ node.get_long_range_node(0, -1, -1).x[comp_i->index].value + node.get_long_range_node(0, -1, 1).x[comp_i->index].value
											+ node.get_long_range_node(0, 1, -1).x[comp_i->index].value + node.get_long_range_node(0, 1, 1).x[comp_i->index].value
											- 36.0 * comp_i->value) / 6.0 / dx / dx;
									}
								}
								for (auto comp_i = node.x.begin(); comp_i < node.x.end(); comp_i++) {
									comp_i->increment = 0.0;
									for (auto comp_j = node.x.begin(); comp_j < node.x.end(); comp_j++) {
										Vector3 vec_Mij;
										double M_ij = Mij(node, comp_i->index, comp_j->index, information->dynamicCollection);
										vec_Mij[0] = (Mij(node.get_neighbor_node(pf::Direction::x_down), comp_i->index, comp_j->index, information->dynamicCollection) - Mij(node.get_neighbor_node(pf::Direction::x_up), comp_i->index, comp_j->index, information->dynamicCollection)) / 2.0 / dx;
										vec_Mij[1] = (Mij(node.get_neighbor_node(pf::Direction::y_down), comp_i->index, comp_j->index, information->dynamicCollection) - Mij(node.get_neighbor_node(pf::Direction::y_up), comp_i->index, comp_j->index, information->dynamicCollection)) / 2.0 / dx;
										vec_Mij[2] = (Mij(node.get_neighbor_node(pf::Direction::z_down), comp_i->index, comp_j->index, information->dynamicCollection) - Mij(node.get_neighbor_node(pf::Direction::z_up), comp_i->index, comp_j->index, information->dynamicCollection)) / 2.0 / dx;
										comp_i->increment += vec_Mij * comp_j->gradient + M_ij * comp_j->laplacian;
									}
									comp_i->increment += Si(node, comp_i->index, information->dynamicCollection);
								}
							}
							else if (flux_DrivingForce == pf::PhaseFluxModel::IntDiff_PotentialGrad) {
								double vm = information->materialSystem.functions.MolarVolume(node, node[0], information->dynamicCollection);
								for (auto comp_i = node.x.begin(); comp_i < node.x.end(); comp_i++) {
									node.potential[comp_i->index].gradient[0] = (node.get_neighbor_node(pf::Direction::x_down).potential[comp_i->index].value - node.get_neighbor_node(pf::Direction::x_up).potential[comp_i->index].value) / 2.0 / dx * vm;
									node.potential[comp_i->index].gradient[1] = (node.get_neighbor_node(pf::Direction::y_down).potential[comp_i->index].value - node.get_neighbor_node(pf::Direction::y_up).potential[comp_i->index].value) / 2.0 / dx * vm;
									node.potential[comp_i->index].gradient[2] = (node.get_neighbor_node(pf::Direction::z_down).potential[comp_i->index].value - node.get_neighbor_node(pf::Direction::z_up).potential[comp_i->index].value) / 2.0 / dx * vm;
									if (differenceMethod == pf::DifferenceMethod::FIVE_POINT) {
										node.potential[comp_i->index].laplacian = (node.get_neighbor_node(pf::Direction::x_down).potential[comp_i->index].value + node.get_neighbor_node(pf::Direction::x_up).potential[comp_i->index].value
											+ node.get_neighbor_node(pf::Direction::y_down).potential[comp_i->index].value + node.get_neighbor_node(pf::Direction::y_up).potential[comp_i->index].value
											+ node.get_neighbor_node(pf::Direction::z_down).potential[comp_i->index].value + node.get_neighbor_node(pf::Direction::z_up).potential[comp_i->index].value
											- 6.0 * node.potential[comp_i->index].value) / dx / dx * vm;
									}
									else if (differenceMethod == pf::DifferenceMethod::NINE_POINT) {
										node.potential[comp_i->index].laplacian = (4.0 * node.get_neighbor_node(pf::Direction::x_down).potential[comp_i->index].value + 4.0 * node.get_neighbor_node(pf::Direction::x_up).potential[comp_i->index].value
											+ 4.0 * node.get_neighbor_node(pf::Direction::y_down).potential[comp_i->index].value + 4.0 * node.get_neighbor_node(pf::Direction::y_up).potential[comp_i->index].value
											+ 4.0 * node.get_neighbor_node(pf::Direction::z_down).potential[comp_i->index].value + 4.0 * node.get_neighbor_node(pf::Direction::z_up).potential[comp_i->index].value
											+ node.get_long_range_node(-1, -1, 0).potential[comp_i->index].value + node.get_long_range_node(-1, 1, 0).potential[comp_i->index].value
											+ node.get_long_range_node(1, -1, 0).potential[comp_i->index].value + node.get_long_range_node(1, 1, 0).potential[comp_i->index].value
											+ node.get_long_range_node(-1, 0, -1).potential[comp_i->index].value + node.get_long_range_node(-1, 0, 1).potential[comp_i->index].value
											+ node.get_long_range_node(1, 0, -1).potential[comp_i->index].value + node.get_long_range_node(1, 0, 1).potential[comp_i->index].value
											+ node.get_long_range_node(0, -1, -1).potential[comp_i->index].value + node.get_long_range_node(0, -1, 1).potential[comp_i->index].value
											+ node.get_long_range_node(0, 1, -1).potential[comp_i->index].value + node.get_long_range_node(0, 1, 1).potential[comp_i->index].value
											- 36.0 * node.potential[comp_i->index].value) / 6.0 / dx / dx * vm;
									}
								}
								for (auto comp_i = node.x.begin(); comp_i < node.x.end(); comp_i++) {
									comp_i->increment = 0.0;
									for (auto comp_j = node.x.begin(); comp_j < node.x.end(); comp_j++) {
										Vector3 vec_Mij;
										double M_ij = Mij(node, comp_i->index, comp_j->index, information->dynamicCollection);
										vec_Mij[0] = (Mij(node.get_neighbor_node(pf::Direction::x_down), comp_i->index, comp_j->index, information->dynamicCollection) - Mij(node.get_neighbor_node(pf::Direction::x_up), comp_i->index, comp_j->index, information->dynamicCollection)) / 2.0 / dx;
										vec_Mij[1] = (Mij(node.get_neighbor_node(pf::Direction::y_down), comp_i->index, comp_j->index, information->dynamicCollection) - Mij(node.get_neighbor_node(pf::Direction::y_up), comp_i->index, comp_j->index, information->dynamicCollection)) / 2.0 / dx;
										vec_Mij[2] = (Mij(node.get_neighbor_node(pf::Direction::z_down), comp_i->index, comp_j->index, information->dynamicCollection) - Mij(node.get_neighbor_node(pf::Direction::z_up), comp_i->index, comp_j->index, information->dynamicCollection)) / 2.0 / dx;
										comp_i->increment += (vec_Mij * node.potential[comp_i->index].gradient + M_ij * node.potential[comp_i->index].laplacian) * vm;
									}
									comp_i->increment += Si(node, comp_i->index, information->dynamicCollection) * vm;
								}
							};
						}
				double MAX_COMP_INCREMENT = 0.0;
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							for (auto comp_i = node.x.begin(); comp_i < node.x.end(); comp_i++) {
								if (abs(dt * comp_i->increment) > MAX_COMP_INCREMENT)
									MAX_COMP_INCREMENT = abs(dt * comp_i->increment);
								comp_i->value += dt * comp_i->increment;
								if (adjust_c_0_1) {
									if (comp_i->value >= 1.0)
										comp_i->value = 1.0;
									else if (comp_i->value <= 0.0)
										comp_i->value = 0.0;
								}
							}
						}
				return MAX_COMP_INCREMENT;
			}

			void update_flag_by_phase_flag() {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							bool on_interface = false;
							for (auto phase = node.begin(); phase < node.end(); phase++)
								if (phase->_flag == pf_INTERFACE)
									on_interface = true;
							if (on_interface) {
								if (node.customFlags[flag_index] == TotalConcentrationFlag::OUT_OF_INTERFACE)
									node.customFlags[flag_index] = TotalConcentrationFlag::NEW_IN_INTERFACE;
								else
									node.customFlags[flag_index] = TotalConcentrationFlag::ITERATED;
							}
							else {
								node.customFlags[flag_index] = TotalConcentrationFlag::OUT_OF_INTERFACE;
							}
						}
			}

		private:
			pf::FieldStorage_forPhaseNode* phaseMesh;
			pf::Information* information;
			double (*Mij)(pf::PhaseNode&, int, int, pf::Info_DynamicCollection&);
			double (*Si)(pf::PhaseNode&, int, pf::Info_DynamicCollection&);
			void (*Boundary)(pf::PhaseNode&, pf::Info_DynamicCollection&);
			void (*potential_func)(pf::PhaseNode&, int, pf::Info_DynamicCollection&);
			int flag_index;
		};

		namespace DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL {
			static double Mij(pf::PhaseNode& node, int comp_i, int comp_j, pf::Info_DynamicCollection& inf) {
				return 0.0;
			}
			static double Source(pf::PhaseNode& node, int comp_i, pf::Info_DynamicCollection& inf) {
				return 0.0;
			}
			static void Boundary(pf::PhaseNode& node, pf::Info_DynamicCollection& inf) {
				return;
			}
			static double dphase_concentration_dphase_potential(pf::PhaseNode& node, pf::PhaseEntry& phase, int comp_index, int potential_index, pf::Info_DynamicCollection& inf) {
				return 0.0;
			}
			static void phase_potential(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
				return;
			}
			static void phase_concentration(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
				return;
			}
			static void phase_energy_density(pf::PhaseNode& node, pf::PhaseEntry& phase, pf::Info_DynamicCollection& inf) {
				return;
			}
		}
#define INCREMENT other_part
		class GrandPotentialFormulation
		{
		public:
			GrandPotentialFormulation(pf::FieldStorage_forPhaseNode& _mesh, pf::Information& _inf) {
				phaseMesh = &_mesh;
				information = &_inf;
				Mij = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::Mij;
				Source = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::Source;
				dphase_concentration_dphase_potential = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::dphase_concentration_dphase_potential;
				Boundary = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::Boundary;
				phase_potential = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::phase_potential;
				phase_concentration = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::phase_concentration;
				phase_energy_density = DEFULT_FUNCTIONS_FOR_GRAND_POTENTIAL::phase_energy_density;
			}

			~GrandPotentialFormulation() {
				phaseMesh = nullptr;
				information = nullptr;
				dphase_concentration_dphase_potential = nullptr;
				Boundary = nullptr;
				phase_potential = nullptr;
				phase_concentration = nullptr;
				phase_energy_density = nullptr;
			}
			void define_Mij(double (*_M_ij)(pf::PhaseNode&, int, int, pf::Info_DynamicCollection&)) {
				Mij = _M_ij;
			}
			void define_Source(double (*_Source)(pf::PhaseNode&, int, pf::Info_DynamicCollection&)) {
				Source = _Source;
			}
			void define_Boundary(void (*_Boundary)(pf::PhaseNode&, pf::Info_DynamicCollection&)) {
				Boundary = _Boundary;
			}
			void define_dphase_concentration_dphase_potential(double (*_dphase_concentration_dphase_potential)(pf::PhaseNode&, pf::PhaseEntry&, int, int, pf::Info_DynamicCollection&)) {
				dphase_concentration_dphase_potential = _dphase_concentration_dphase_potential;
			}
			void define_phase_potential(void (*_phase_potential)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&)) {
				phase_potential = _phase_potential;
			}
			void define_phase_concentration(void (*_phase_concentration)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&)) {
				phase_concentration = _phase_concentration;
			}
			void define_phase_energy_density(void (*_phase_energy_density)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&)) {
				phase_energy_density = _phase_energy_density;
			}

			void init_potential_field_by_phase_concentration() {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							for(auto phase = node.begin(); phase < node.end(); phase++)
								phase_potential(node, *phase, information->dynamicCollection);
							for (auto potential = node.potential.begin(); potential < node.potential.end(); potential++) {
								potential->value = 0.0;
								double sum_fraction = 0.0;
								for (auto phase = node.begin(); phase < node.end(); phase++)
									if (phase->phaseFraction > Simulation_Num_Cut_Off) {
										for (auto p_pot = phase->potential.begin(); p_pot < phase->potential.end(); p_pot++)
											if (p_pot->index == potential->index) {
												sum_fraction += phase->phaseFraction;
												potential->value += p_pot->value * phase->phaseFraction;
											}
									}
								if (sum_fraction < Simulation_Num_Cut_Off)
									potential->value = 0.0;
								else
									potential->value /= sum_fraction;
							}
						}
			}

			// size of phase_region = 0, means the whole domain
			double evolve_potential_field(vector<int> none_cal_region, double cut_off_value = 0.1, double error_balance_coeff = 1.0, int error_balance_range = 1, pf::DifferenceMethod differenceMethod = pf::DifferenceMethod::FIVE_POINT) {
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							double phi = calculate_phi_in_node(none_cal_region, node), G_phi = phi;
							if (phi < cut_off_value) {
								for (auto potential2 = node.potential.begin(); potential2 < node.potential.end(); potential2++)
									potential2->value = 0.0;
								continue;
							}
							else {
								if (1.0 - error_balance_coeff > 1e-3 && error_balance_range > 0) {
									double ave_phi = get_ave_phaseFraction(node, none_cal_region, error_balance_range);
									G_phi = error_balance_coeff * phi + (1 - error_balance_coeff) * ave_phi;
								}
							}
							Vector3 grad_phi;
							grad_phi = cal_laplace_grad_potential(node, none_cal_region, cut_off_value, differenceMethod);
							for (auto potential = node.potential.begin(); potential < node.potential.end(); potential++) {
								double factor = 0.0;
								for (auto phase = node.begin(); phase < node.end(); phase++)
									factor += phase->phaseFraction * dphase_concentration_dphase_potential(node, *phase, potential->index, potential->index, information->dynamicCollection);
								//factor /= phi;
								double diffusion_term = 0.0;
								for (auto potential2 = node.potential.begin(); potential2 < node.potential.end(); potential2++) {
									Vector3 vec_Mij;
									double M_ij = Mij(node, potential->index, potential2->index, information->dynamicCollection);
									vec_Mij[0] = (Mij(node.get_neighbor_node(pf::Direction::x_down), potential->index, potential2->index, information->dynamicCollection)
													- Mij(node.get_neighbor_node(pf::Direction::x_up), potential->index, potential2->index, information->dynamicCollection)) / 2.0 / phaseMesh->dx;
									vec_Mij[1] = (Mij(node.get_neighbor_node(pf::Direction::y_down), potential->index, potential2->index, information->dynamicCollection)
													- Mij(node.get_neighbor_node(pf::Direction::y_up), potential->index, potential2->index, information->dynamicCollection)) / 2.0 / phaseMesh->dx;
									vec_Mij[2] = (Mij(node.get_neighbor_node(pf::Direction::z_down), potential->index, potential2->index, information->dynamicCollection)
													- Mij(node.get_neighbor_node(pf::Direction::z_up), potential->index, potential2->index, information->dynamicCollection)) / 2.0 / phaseMesh->dx;

									diffusion_term += phi * M_ij * potential2->laplacian + (grad_phi * M_ij + vec_Mij * phi) * potential2->gradient;
								}
								potential->INCREMENT = 1.0 / factor * (diffusion_term / G_phi + Source(node, potential->index, information->dynamicCollection));
							}
						}
				double MAX_VARIATION_OF_POTENTIAL = 0.0, dt = information->settings.disperse_settings.dt * information->dynamicCollection.dt_scale;
#pragma omp parallel for
				for (int x = 0; x < phaseMesh->limit_x; x++)
					for (int y = 0; y < phaseMesh->limit_y; y++)
						for (int z = 0; z < phaseMesh->limit_z; z++) {
							PhaseNode& node = (*phaseMesh)(x, y, z);
							for (auto potential = node.potential.begin(); potential < node.potential.end(); potential++) {
								double abs_increment = abs(dt * potential->INCREMENT);
								if (abs_increment > MAX_VARIATION_OF_POTENTIAL)
									MAX_VARIATION_OF_POTENTIAL = abs_increment;
								potential->value += dt * potential->INCREMENT;
							}
						}
				return MAX_VARIATION_OF_POTENTIAL;
			}

			void cal_concentration_and_energy_density_by_potential(int istep, int step_cal_full_region) {
				if (istep % step_cal_full_region == 0) {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								for (auto phase = node.begin(); phase < node.end(); phase++) {
									phase_concentration(node, *phase, information->dynamicCollection);
									phase_energy_density(node, *phase, information->dynamicCollection);
								}
								node.cal_x_from_phase_x();
							}
				}
				else {
#pragma omp parallel for
					for (int x = 0; x < phaseMesh->limit_x; x++)
						for (int y = 0; y < phaseMesh->limit_y; y++)
							for (int z = 0; z < phaseMesh->limit_z; z++) {
								PhaseNode& node = (*phaseMesh)(x, y, z);
								for (auto phase = node.begin(); phase < node.end(); phase++)
									if (phase->_flag == pf::pf_INTERFACE) {
										phase_concentration(node, *phase, information->dynamicCollection);
										phase_energy_density(node, *phase, information->dynamicCollection);
									}
							}
				}
			}

		private:
			pf::FieldStorage_forPhaseNode* phaseMesh;
			pf::Information* information;
			double (*Mij)(pf::PhaseNode&, int, int, pf::Info_DynamicCollection&);
			double (*Source)(pf::PhaseNode&, int, pf::Info_DynamicCollection&);
			double (*dphase_concentration_dphase_potential)(pf::PhaseNode&, pf::PhaseEntry&, int, int, pf::Info_DynamicCollection&);
			void (*Boundary)(pf::PhaseNode&, pf::Info_DynamicCollection&);

			void (*phase_potential)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
			void (*phase_concentration)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
			void (*phase_energy_density)(pf::PhaseNode&, pf::PhaseEntry&, pf::Info_DynamicCollection&);
			Vector3 cal_laplace_grad_potential(pf::PhaseNode& node, vector<int> phase_region, double cut_off_value, pf::DifferenceMethod differenceMethod) {
				Vector3 phi_grad;
				phi_grad.set_to_zero();
				if (differenceMethod == DifferenceMethod::FIVE_POINT) {
					pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
					pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
					pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
					pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
					pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
					pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
					int xlength = 2;
					int ylength = 2;
					int zlength = 2;
					if (calculate_phi_in_node(phase_region, *node_upx) < cut_off_value) {
						node_upx = &node;
						xlength--;
					}
					if (calculate_phi_in_node(phase_region, *node_downx) < cut_off_value) {
						node_downx = &node;
						xlength--;
					}
					if (calculate_phi_in_node(phase_region, *node_upy) < cut_off_value) {
						node_upy = &node;
						ylength--;
					}
					if (calculate_phi_in_node(phase_region, *node_downy) < cut_off_value) {
						node_downy = &node;
						ylength--;
					}
					if (calculate_phi_in_node(phase_region, *node_upz) < cut_off_value) {
						node_upz = &node;
						zlength--;
					}
					if (calculate_phi_in_node(phase_region, *node_downz) < cut_off_value) {
						node_downz = &node;
						zlength--;
					}
					phi_grad[0] = (calculate_phi_in_node(phase_region, *node_downx) - calculate_phi_in_node(phase_region, *node_upx)) / xlength / phaseMesh->dx;
					phi_grad[1] = (calculate_phi_in_node(phase_region, *node_downy) - calculate_phi_in_node(phase_region, *node_upy)) / ylength / phaseMesh->dx;
					phi_grad[2] = (calculate_phi_in_node(phase_region, *node_downz) - calculate_phi_in_node(phase_region, *node_upz)) / zlength / phaseMesh->dx;
					for (auto chem = node.potential.begin(); chem < node.potential.end(); chem++) {

						if (xlength != 0)
							chem->gradient[0] = (node_downx->potential[chem->index].value -
								node_upx->potential[chem->index].value) / xlength / phaseMesh->dx;
						else
							chem->gradient[0] = 0.0;

						if (ylength != 0)
							chem->gradient[1] = (node_downy->potential[chem->index].value -
								node_upy->potential[chem->index].value) / ylength / phaseMesh->dx;
						else
							chem->gradient[1] = 0.0;

						if (zlength != 0)
							chem->gradient[2] = (node_downz->potential[chem->index].value -
								node_upz->potential[chem->index].value) / zlength / phaseMesh->dx;
						else
							chem->gradient[2] = 0.0;

						chem->laplacian = (node_downx->potential[chem->index].value
							+ node_upx->potential[chem->index].value
							+ node_downy->potential[chem->index].value
							+ node_upy->potential[chem->index].value
							+ node_downz->potential[chem->index].value
							+ node_upz->potential[chem->index].value - 6 * chem->value) / phaseMesh->dx / phaseMesh->dx;

					}
				}
				else if (differenceMethod == DifferenceMethod::NINE_POINT) {
					pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
					pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
					pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
					pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
					pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
					pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
					int xlength = 2;
					int ylength = 2;
					int zlength = 2;
					if (calculate_phi_in_node(phase_region, *node_upx) < cut_off_value) {
						node_upx = &node;
						xlength--;
					}
					if (calculate_phi_in_node(phase_region, *node_downx) < cut_off_value) {
						node_downx = &node;
						xlength--;
					}
					if (calculate_phi_in_node(phase_region, *node_upy) < cut_off_value) {
						node_upy = &node;
						ylength--;
					}
					if (calculate_phi_in_node(phase_region, *node_downy) < cut_off_value) {
						node_downy = &node;
						ylength--;
					}
					if (calculate_phi_in_node(phase_region, *node_upz) < cut_off_value) {
						node_upz = &node;
						zlength--;
					}
					if (calculate_phi_in_node(phase_region, *node_downz) < cut_off_value) {
						node_downz = &node;
						zlength--;
					}
					phi_grad[0] = (calculate_phi_in_node(phase_region, *node_downx) - calculate_phi_in_node(phase_region, *node_upx)) / xlength / phaseMesh->dx;
					phi_grad[1] = (calculate_phi_in_node(phase_region, *node_downy) - calculate_phi_in_node(phase_region, *node_upy)) / ylength / phaseMesh->dx;
					phi_grad[2] = (calculate_phi_in_node(phase_region, *node_downz) - calculate_phi_in_node(phase_region, *node_upz)) / zlength / phaseMesh->dx;
					pf::PhaseNode* node_upxupy = &node.get_long_range_node(1, 1, 0);
					pf::PhaseNode* node_downxdowny = &node.get_long_range_node(-1, -1, 0);
					pf::PhaseNode* node_upydownx = &node.get_long_range_node(-1, 1, 0);
					pf::PhaseNode* node_downyupx = &node.get_long_range_node(1, -1, 0);
					pf::PhaseNode* node_upxupz = &node.get_long_range_node(1, 0, 1);
					pf::PhaseNode* node_downxdownz = &node.get_long_range_node(-1, 0, -1);
					pf::PhaseNode* node_upzdownx = &node.get_long_range_node(-1, 0, 1);
					pf::PhaseNode* node_downzupx = &node.get_long_range_node(1, 0, -1);
					pf::PhaseNode* node_upzupy = &node.get_long_range_node(0, 1, 1);
					pf::PhaseNode* node_downzdowny = &node.get_long_range_node(0, -1, -1);
					pf::PhaseNode* node_upydownz = &node.get_long_range_node(0, 1, -1);
					pf::PhaseNode* node_downyupz = &node.get_long_range_node(0, -1, 1);
					if (calculate_phi_in_node(phase_region, *node_upxupy) < cut_off_value)
						node_upxupy = &node;
					if (calculate_phi_in_node(phase_region, *node_downxdowny) < cut_off_value)
						node_downxdowny = &node;
					if (calculate_phi_in_node(phase_region, *node_upydownx) < cut_off_value)
						node_upydownx = &node;
					if (calculate_phi_in_node(phase_region, *node_downyupx) < cut_off_value)
						node_downyupx = &node;
					if (calculate_phi_in_node(phase_region, *node_upxupz) < cut_off_value)
						node_upxupz = &node;
					if (calculate_phi_in_node(phase_region, *node_downxdownz) < cut_off_value)
						node_downxdownz = &node;
					if (calculate_phi_in_node(phase_region, *node_upzdownx) < cut_off_value)
						node_upzdownx = &node;
					if (calculate_phi_in_node(phase_region, *node_downzupx) < cut_off_value)
						node_downzupx = &node;
					if (calculate_phi_in_node(phase_region, *node_upzupy) < cut_off_value)
						node_upzupy = &node;
					if (calculate_phi_in_node(phase_region, *node_downzdowny) < cut_off_value)
						node_downzdowny = &node;
					if (calculate_phi_in_node(phase_region, *node_upydownz) < cut_off_value)
						node_upydownz = &node;
					if (calculate_phi_in_node(phase_region, *node_downyupz) < cut_off_value)
						node_downyupz = &node;
					for (auto chem = node.potential.begin(); chem < node.potential.end(); chem++) {
						if (xlength != 0)
							chem->gradient[0] = (node_downx->potential[chem->index].value -
								node_upx->potential[chem->index].value) / xlength / phaseMesh->dx;
						else
							chem->gradient[0] = 0.0;
						if (ylength != 0)
							chem->gradient[1] = (node_downy->potential[chem->index].value -
								node_upy->potential[chem->index].value) / ylength / phaseMesh->dx;
						else
							chem->gradient[1] = 0.0;
						if (zlength != 0)
							chem->gradient[2] = (node_downz->potential[chem->index].value -
								node_upz->potential[chem->index].value) / zlength / phaseMesh->dx;
						else
							chem->gradient[2] = 0.0;
						chem->laplacian = (4 * node_downx->potential[chem->index].value + 4 * node_upx->potential[chem->index].value
							+ 4 * node_downy->potential[chem->index].value + 4 * node_upy->potential[chem->index].value
							+ 4 * node_downz->potential[chem->index].value + 4 * node_upz->potential[chem->index].value
							+ node_upxupy->potential[chem->index].value   + node_downxdowny->potential[chem->index].value
							+ node_upydownx->potential[chem->index].value + node_downyupx->potential[chem->index].value
							+ node_upxupz->potential[chem->index].value   + node_downxdownz->potential[chem->index].value
							+ node_upzdownx->potential[chem->index].value + node_downzupx->potential[chem->index].value
							+ node_upzupy->potential[chem->index].value   + node_downzdowny->potential[chem->index].value
							+ node_upydownz->potential[chem->index].value + node_downyupz->potential[chem->index].value - 36 * chem->value) / 6.0 / phaseMesh->dx / phaseMesh->dx;
					}
				}
				return phi_grad;
			}
			double calculate_phi_in_node(vector<int> none_cal_region, pf::PhaseNode& node) {
				if (none_cal_region.size() == 0)
					return 1.0;
				else {
					double phi = 0.0;
					for (auto phaseIndex = none_cal_region.begin(); phaseIndex < none_cal_region.end(); phaseIndex++)
						phi += node[*phaseIndex].phaseFraction;
					if ((1.0 - phi) > 1.0)
						return 1.0;
					else if ((1.0 - phi) < 0.0)
						return 0.0;
				}
			}
			double get_ave_phaseFraction(pf::PhaseNode& node, vector<int> phase_region, int opt_range = 1) {
				int num = 0;
				double p_f = 0.0;
				for (int i = -opt_range; i <= +opt_range; i++)
					for (int j = -opt_range; j <= +opt_range; j++)
						for (int k = -opt_range; k <= +opt_range; k++) {
							p_f += calculate_phi_in_node(phase_region, node.get_long_range_node(i, j, k));
							num += 1;
						}
				if (num > 0)
					p_f = p_f / num;
				else
					p_f = calculate_phi_in_node(phase_region, node);
				return p_f;
			}
		};
	}
}