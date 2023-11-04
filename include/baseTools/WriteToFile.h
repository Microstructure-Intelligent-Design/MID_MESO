#pragma once
#include "sysTool.h"
#include "FieldStorage.h"
#include "NucleationTools.h"
#include "InfoTools.h"
using namespace std;

namespace pf{
	class WriteToFile
	{
	public:
		WriteToFile(string file_path, double _NaN) {
			init(file_path, _NaN);
		};
		WriteToFile() {};
		~WriteToFile() {};
		void init(string file_path, double _NaN) {
			_path = file_path; 
			NaN = _NaN;
#if defined(_WIN32)
			_mkdir(_path.c_str());
#elif defined(__linux__)
			mkdir(_path.c_str(), 0777);
			// #else more?
#endif
		}

		// open vts file
		void open_vts_scalar_file(ofstream& fout, FieldStorage_forPhaseNode& fs, int istep);
		void open_vts_vec3_file(ofstream& fout, FieldStorage_forPhaseNode& fs, int istep);
		// close vts file
		void close_vts_file(ofstream& fout, FieldStorage_forPhaseNode& fs);
		// write scalar data
		void write_scalar_phi(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Phases settings_phases);
		void write_scalar_con(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con);
		void write_scalar_phasecon(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con);
		void write_scalar_potential(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con);
		void write_scalar_phasepotential(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con);
		void write_scalar_partPotential(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con);
		void write_scalar_energyDensity(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_grains(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_grainsReserve(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_phiIndexs(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_phiProperties(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_Flags(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_temperatureField(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_electricField(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_magneticField(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_elasticField(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_fluidField(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_scalar_nucleationField(ofstream& fout, NP_FieldStorage& np, vector<int> np_phase_index, Set_OutputFile outer_settings);
		// write vector data
		void write_vec3_phiGradient(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_vec3_demagnetizingIntensity(ofstream& fout, FieldStorage_forPhaseNode& fs);
		void write_vec3_velocity(ofstream& fout, FieldStorage_forPhaseNode& fs);

		void write_nucleation_phase_field(NP_FieldStorage& np, int istep, vector<int> np_phase_index, Set_OutputFile outer_settings);

		void write_string_to_txt(string content, string _fname);
		void add_string_to_txt(string content, string _fname);
		void init_txt_file(string _fname);

		string _path;
		double NaN;
	};
	inline void WriteToFile::open_vts_scalar_file(ofstream& fout, FieldStorage_forPhaseNode& fs, int istep) {
		string fname;
		fname = _path + dirSeparator + "scalar_variables_step" + to_string(istep) + ".vts";
		
		fout.open(fname);
		if (!fout) {
			cout << "Failed to write the vtk file..." << endl;
			fout.close();
			return;
		}
		fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
		fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
		fout << "<StructuredGrid WholeExtent=\""
			<< 0 << " " << fs.limit_x - 1 << " "
			<< 0 << " " << fs.limit_y - 1 << " "
			<< 0 << " " << fs.limit_z - 1 << "\"> " << endl;
		fout << "<PointData Scalars= \"ScalarData\">" << endl;
	}
	inline void WriteToFile::open_vts_vec3_file(ofstream& fout, FieldStorage_forPhaseNode& fs, int istep) {
		string fname;
		fname = _path + dirSeparator + "vec3_variables_step" + to_string(istep) + ".vts";
		fout.open(fname);
		if (!fout) {
			cout << "Failed to write the vtk file..." << endl;
			fout.close();
			return;
		}
		fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
		fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
		fout << "<StructuredGrid WholeExtent=\""
			<< 0 << " " << fs.limit_x - 1 << " "
			<< 0 << " " << fs.limit_y - 1 << " "
			<< 0 << " " << fs.limit_z - 1 << "\"> " << endl;
		fout << "<PointData  Vectors= \"VectorData\">" << endl;
	}
	inline void WriteToFile::close_vts_file(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "</PointData>" << endl;
		fout << "<Points>" << endl;
		fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					fout << i << " " << j << " " << k << "\n";
				}
		fout << "</DataArray>" << endl;
		fout << "</Points>" << endl;
		fout << "</StructuredGrid>" << endl;
		fout << "</VTKFile>" << endl;
		fout.close();
	}

	inline void WriteToFile::write_scalar_phi(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Phases settings_phases) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto phase = sample_node.begin(); phase < sample_node.end(); phase++) {
			string name;
			name = "\"Phi" + to_string(phase->index) + "_" + settings_phases[phase->phaseProperty].phase_name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; ++k)
				for (int j = 0; j < fs.limit_y; ++j)
					for (int i = 0; i < fs.limit_x; ++i) {
						double phi = fs(i, j, k)[phase->index].phaseFraction;
						if(IS_NAN(phi))
							fout << NaN << endl;
						else
							fout << phi << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void WriteToFile::write_scalar_Flags(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto phase = sample_node.begin(); phase < sample_node.end(); phase++) {
			string name = "\"Phi" + to_string(phase->index) + "_flag\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; ++k)
				for (int j = 0; j < fs.limit_y; ++j)
					for (int i = 0; i < fs.limit_x; ++i)
					{
						fout << fs(i, j, k)[phase->index]._flag << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void WriteToFile::write_scalar_con(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto comp = sample_node.x.begin(); comp < sample_node.x.end(); comp++) {
			string name = "\"con_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; ++k)
				for (int j = 0; j < fs.limit_y; ++j)
					for (int i = 0; i < fs.limit_x; ++i) {
						double con = fs(i, j, k).x[comp->index].value;
						if(IS_NAN(con))
							fout << NaN << endl;
						else
							fout << con << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void WriteToFile::write_scalar_phasecon(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto phase = sample_node.begin(); phase < sample_node.end(); phase++) {
			for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
				string name = "\"phi" + to_string(phase->index) + "_con" + to_string(comp->index) + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << name <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < fs.limit_z; ++k)
					for (int j = 0; j < fs.limit_y; ++j)
						for (int i = 0; i < fs.limit_x; ++i) {
							if (fs(i, j, k)[phase->index].phaseFraction > Simulation_Num_Cut_Off) {
								double con = fs(i, j, k)[phase->index].x[comp->index].value;
								if(IS_NAN(con))
									fout << NaN << endl;
								else
									fout << con << endl;
							}
							else {
								fout << NaN << endl;
							}
						}
				fout << "</DataArray>" << endl;
			}
		}
	}
	inline void WriteToFile::write_scalar_potential(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto comp = sample_node.x.begin(); comp < sample_node.x.end(); comp++) {
			string name = "\"potential_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k).potential[comp->index].value << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void WriteToFile::write_scalar_phasepotential(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto phase = sample_node.begin(); phase < sample_node.end(); phase++) {
			for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
				string name = "\"phi" + to_string(phase->index) + "_pot" + to_string(comp->index) + "\" ";
				fout << "<DataArray type = \"Float64\" Name = " << name <<
					"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
				for (int k = 0; k < fs.limit_z; ++k)
					for (int j = 0; j < fs.limit_y; ++j)
						for (int i = 0; i < fs.limit_x; ++i) {
							if (fs(i, j, k)[phase->index].phaseFraction > Simulation_Num_Cut_Off) {
								double con = fs(i, j, k)[phase->index].potential[comp->index].value;
								if (IS_NAN(con))
									fout << NaN << endl;
								else
									fout << con << endl;
							}
							else {
								fout << NaN << endl;
							}
						}
				fout << "</DataArray>" << endl;
			}
		}
	}
	inline void WriteToFile::write_scalar_partPotential(ofstream& fout, FieldStorage_forPhaseNode& fs, Info_Node system_con) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto comp = sample_node.x.begin(); comp < sample_node.x.end(); comp++) {
			string name = "\"pot_chem_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k).potential[comp->index].chemical_part << endl;
					}
			fout << "</DataArray>" << endl;

			name = "\"pot_elas_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k).potential[comp->index].elastic_part << endl;
					}
			fout << "</DataArray>" << endl;

			name = "\"pot_elec_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k).potential[comp->index].electric_part << endl;
					}
			fout << "</DataArray>" << endl;

			name = "\"pot_mag_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k).potential[comp->index].magnetic_part << endl;
					}
			fout << "</DataArray>" << endl;

			name = "\"pot_other_" + to_string(comp->index) + "_" + system_con[comp->index].name + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k).potential[comp->index].other_part << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void WriteToFile::write_scalar_energyDensity(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "energyDensity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					double ed = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						ed += phase->phaseFraction * (phase->chemEnergyDensity + phase->elasEnergyDensity + phase->elecEnergyDensity + phase->magEnergyDensity + phase->otherEnergyDensity);
					if (IS_NAN(ed))
						fout << NaN << endl;
					else
						fout << ed << endl;
				}
		fout << "</DataArray>" << endl;

		fout << "<DataArray type = \"Float64\" Name = \"" << "chemEnergyDensity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					double ed = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						ed += phase->phaseFraction * (phase->chemEnergyDensity);
					if (IS_NAN(ed))
						fout << NaN << endl;
					else
						fout << ed << endl;
				}
		fout << "</DataArray>" << endl;

		fout << "<DataArray type = \"Float64\" Name = \"" << "elasEnergyDensity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					double ed = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						ed += phase->phaseFraction * (phase->elasEnergyDensity);
					if (IS_NAN(ed))
						fout << NaN << endl;
					else
						fout << ed << endl;
				}
		fout << "</DataArray>" << endl;

		fout << "<DataArray type = \"Float64\" Name = \"" << "elecEnergyDensity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					double ed = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						ed += phase->phaseFraction * (phase->elecEnergyDensity);
					if (IS_NAN(ed))
						fout << NaN << endl;
					else
						fout << ed << endl;
				}
		fout << "</DataArray>" << endl;

		fout << "<DataArray type = \"Float64\" Name = \"" << "magEnergyDensity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					double ed = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						ed += phase->phaseFraction * (phase->magEnergyDensity);
					if (IS_NAN(ed))
						fout << NaN << endl;
					else
						fout << ed << endl;
				}
		fout << "</DataArray>" << endl;

		fout << "<DataArray type = \"Float64\" Name = \"" << "otherEnergyDensity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					PhaseNode& node = fs(i, j, k);
					double ed = 0.0;
					for (auto phase = node.begin(); phase < node.end(); phase++)
						ed += phase->phaseFraction * (phase->otherEnergyDensity);
					if (IS_NAN(ed))
						fout << NaN << endl;
					else
						fout << ed << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_grains(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "grains" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					double fix = 0.0;
					for (auto phase = fs(i, j, k).begin(); phase < fs(i, j, k).end(); phase++)
						fix += phase->phaseFraction * phase->phaseFraction;
					if(IS_NAN(fix))
						fout << NaN << endl;
					else
						fout << fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_grainsReserve(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "grains_reserve" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					double fix = 0.0;
					for (auto phase = fs(i, j, k).begin(); phase < fs(i, j, k).end(); phase++)
						fix += phase->phaseFraction * phase->phaseFraction;
					if (IS_NAN(fix))
						fout << NaN << endl;
					else
						fout << 1.0 - fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_phiIndexs(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "phiIndexs" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					double fix = 0.0;
					for (auto phase = fs(i, j, k).begin(); phase < fs(i, j, k).end(); phase++)
						fix += phase->index * phase->phaseFraction;
					fout << fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_phiProperties(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "phiProperties" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					double fix = 0.0;
					for (auto phase = fs(i, j, k).begin(); phase < fs(i, j, k).end(); phase++)
						fix += phase->phaseProperty * phase->phaseFraction;
					fout << fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_temperatureField(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "temperature" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).tempValues.temperature << endl;
				}
		fout << "</DataArray>" << endl;
		fout << "<DataArray type = \"Float64\" Name = \"" << "heat_diffusivity" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).tempValues.heat_diffusivity << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_electricField(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "elec_potential" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).electricValues.elec_potential << endl;
				}
		fout << "</DataArray>" << endl;
		fout << "<DataArray type = \"Float64\" Name = \"" << "elec_charge_density" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).electricValues.elec_charge_density << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_magneticField(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "mag_potential" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).magneticValues.mag_potential << endl;
				}
		fout << "</DataArray>" << endl;
		fout << "<DataArray type = \"Float64\" Name = \"" << "mag_charge_density" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).magneticValues.mag_charge_density << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_elasticField(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };
		for (int ele_index = 0; ele_index < 6; ele_index++)
		{
			string compname = "\"stress_" + compNameV[ele_index] + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << compname <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; ++k)
				for (int j = 0; j < fs.limit_y; ++j)
					for (int i = 0; i < fs.limit_x; ++i)
					{
						fout << fs(i, j, k).mechanicalValues.Stresses[ele_index] << endl;
					}
			fout << "</DataArray>" << endl;
		}
		for (int ele_index = 0; ele_index < 6; ele_index++)
		{
			string compname = "\"strain_" + compNameV[ele_index] + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << compname <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; ++k)
				for (int j = 0; j < fs.limit_y; ++j)
					for (int i = 0; i < fs.limit_x; ++i)
					{
						fout << fs(i, j, k).mechanicalValues.Strains[ele_index] << endl;
					}
			fout << "</DataArray>" << endl;
		}
		for (int ele_index = 0; ele_index < 6; ele_index++)
		{
			string compname = "\"eigenStrain_" + compNameV[ele_index] + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << compname <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; ++k)
				for (int j = 0; j < fs.limit_y; ++j)
					for (int i = 0; i < fs.limit_x; ++i)
					{
						fout << fs(i, j, k).mechanicalValues.EffectiveEigenStrain[ele_index] << endl;
					}
			fout << "</DataArray>" << endl;
		}

		fout << "<DataArray type = \"Float64\" Name = \"" << "J1" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					fout << fs(i, j, k).mechanicalValues.Stresses.J1() << endl;
				}
		fout << "</DataArray>" << endl;

		fout << "<DataArray type = \"Float64\" Name = \"" << "vMises" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; ++k)
			for (int j = 0; j < fs.limit_y; ++j)
				for (int i = 0; i < fs.limit_x; ++i)
				{
					fout << fs(i, j, k).mechanicalValues.Stresses.Mises() << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_fluidField(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "pressure" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).velocityValues.pressure << endl;
				}
		fout << "</DataArray>" << endl;
		fout << "<DataArray type = \"Float64\" Name = \"" << "velocity_u" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).velocityValues.velocity[0] << endl;
				}
		fout << "</DataArray>" << endl;
		fout << "<DataArray type = \"Float64\" Name = \"" << "velocity_v" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).velocityValues.velocity[1] << endl;
				}
		fout << "</DataArray>" << endl;
		fout << "<DataArray type = \"Float64\" Name = \"" << "velocity_w" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).velocityValues.velocity[2] << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_scalar_nucleationField(ofstream& fout, NP_FieldStorage& np, vector<int> np_phase_index, Set_OutputFile outer_settings) {
		for (auto phaseIndex = np_phase_index.begin(); phaseIndex < np_phase_index.end(); phaseIndex++){
			string name;
			name = "\"Nucleation_Phi" + to_string(*phaseIndex) + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < np.limit_z; ++k)
				for (int j = 0; j < np.limit_y; ++j)
					for (int i = 0; i < np.limit_x; ++i) {
						bool is_phase_nucleate = false;
						NP_Node& np_node = np(i, j, k);
						for (auto nnp = np_node.conditional_phases.begin(); nnp < np_node.conditional_phases.end(); nnp++) {
							if (nnp->phaseIndex == *phaseIndex)
								is_phase_nucleate = true;
						}
						if (is_phase_nucleate)
							fout << *phaseIndex << endl;
						else
							fout << 0.0 << endl;
					}
			fout << "</DataArray>" << endl;
		}

	}

	inline void WriteToFile::write_vec3_phiGradient(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		PhaseNode& sample_node = fs(0, 0, 0);
		for (auto phase = sample_node.begin(); phase < sample_node.end(); phase++) {
			string name;
			name = "\"Phi" + to_string(phase->index) + "_Grad\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for (int k = 0; k < fs.limit_z; k++)
				for (int j = 0; j < fs.limit_y; j++)
					for (int i = 0; i < fs.limit_x; i++) {
						fout << fs(i, j, k)[phase->index].phi_grad[0] << " "
							<< fs(i, j, k)[phase->index].phi_grad[1] << " "
							<< fs(i, j, k)[phase->index].phi_grad[2] << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void WriteToFile::write_vec3_demagnetizingIntensity(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		string name = "\"demagnetizingIntensity\"";
		fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					PhaseNode& node = fs(i, j, k);
					fout << -(node.get_neighbor_node(Direction::x_down).magneticValues.mag_potential - node.get_neighbor_node(Direction::x_up).magneticValues.mag_potential) / 2.0 / fs.dx << " "
						<< -(node.get_neighbor_node(Direction::y_down).magneticValues.mag_potential - node.get_neighbor_node(Direction::y_up).magneticValues.mag_potential) / 2.0 / fs.dx << " "
						<< -(node.get_neighbor_node(Direction::z_down).magneticValues.mag_potential - node.get_neighbor_node(Direction::z_up).magneticValues.mag_potential) / 2.0 / fs.dx << "\n";
				}
		fout << "</DataArray>" << endl;
	}
	inline void WriteToFile::write_vec3_velocity(ofstream& fout, FieldStorage_forPhaseNode& fs) {
		string name;
		name = "\"velocity\" ";
		fout << "<DataArray type = \"Float64\" Name = " << name << " NumberOfComponents=\"3\" format=\"ascii\">" << endl;
		for (int k = 0; k < fs.limit_z; k++)
			for (int j = 0; j < fs.limit_y; j++)
				for (int i = 0; i < fs.limit_x; i++) {
					fout << fs(i, j, k).velocityValues.velocity[0] << " "
						<< fs(i, j, k).velocityValues.velocity[1] << " "
						<< fs(i, j, k).velocityValues.velocity[2] << endl;
				}
		fout << "</DataArray>" << endl;
	}

	inline void WriteToFile::write_string_to_txt(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to write the txt file!" << endl;
			fout.close();
			return;
		}
		fout << content << endl;
		fout.close();
	}

	inline void WriteToFile::add_string_to_txt(string content, string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname, ios::app);
		if (!fout) {
			cout << "Failed to add the string to txt file!" << endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
	}

	inline void WriteToFile::init_txt_file(string _fname) {
		string fname;
		fname = _path + dirSeparator + _fname + ".txt";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to init the txt file!" << endl;
			fout.close();
			return;
		}
		fout << "";
		fout.close();
	}

	inline void WriteToFile::write_nucleation_phase_field(NP_FieldStorage& np, int istep, vector<int> np_phase_index, Set_OutputFile outer_settings) {
		string fname;
		fname = _path + dirSeparator + "Condition_Nucleation_Region_step" + to_string(istep) + ".vts";
		ofstream fout(fname);
		if (!fout) {
			cout << "Failed to write the vtk file..." << endl;
			fout.close();
			return;
		}

		fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
		fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
		fout << "<StructuredGrid WholeExtent=\""
			<< 0 << " " << np.limit_x - 1 << " "
			<< 0 << " " << np.limit_y - 1 << " "
			<< 0 << " " << np.limit_z - 1 << "\"> " << endl;
		fout << "<PointData Scalars= \"ScalarData\">" << endl;
		for (auto phaseIndex = np_phase_index.begin(); phaseIndex < np_phase_index.end(); phaseIndex++)
		{
			string name;
			name = "\"P" + to_string(*phaseIndex) + "\" ";
			fout << "<DataArray type = \"Float64\" Name = " << name <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < np.limit_z; ++k)
				for (int j = 0; j < np.limit_y; ++j)
					for (int i = 0; i < np.limit_x; ++i) {
						bool is_phase_nucleate = false;
						NP_Node& np_node = np(i, j, k);
						for (auto nnp = np_node.conditional_phases.begin(); nnp < np_node.conditional_phases.end(); nnp++) {
							if (nnp->phaseIndex == *phaseIndex)
								is_phase_nucleate = true;
						}
						if (is_phase_nucleate)
							fout << *phaseIndex << endl;
						else
							fout << 0.0 << endl;
					}
			fout << "</DataArray>" << endl;
		}

		fout << "</PointData>" << endl;
		fout << "<Points>" << endl;
		fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
		for (int k = 0; k < np.limit_z; ++k)
			for (int j = 0; j < np.limit_y; ++j)
				for (int i = 0; i < np.limit_x; ++i)
				{
					fout << i << " " << j << " " << k << "\n";
				}
		fout << "</DataArray>" << endl;
		fout << "</Points>" << endl;
		fout << "</StructuredGrid>" << endl;
		fout << "</VTKFile>" << endl;

		fout.close();
	}

}