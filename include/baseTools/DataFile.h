#pragma once
#include "sysTool.h"
#include "FieldStorage.h"
using namespace std;

#define DataFile_NoneStep -1

namespace pf {
	struct Data_fieldStorage {
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		int Nx;
		int Ny;
		int Nz;
		double dx;
	};
	struct Data_phaseNode {
		double temperature;
		int phaseNum;
		int conNum;
		int custom_value_num;
		int custom_flag_num;
		int custom_vec_num;
	};
	struct Data_custom_value {
		int index;
		int flag;
		double value;
		double value2;
		double value3;
	};
	struct Data_phaseEntry {
		int _flag;
		int index;
		int phaseProperty;
		double phaseFraction;
		int conNum;
	};
	struct Data_nodeEntry {
		int index;
		double value;
	};
	class DataFile
	{
	public:
		DataFile() {
			mainName = "dataFile";
			format = ".dat";
			_path = "";
			DataFile_header = new char[50];
			string head = "## MPF PROGRAM FILE VERSION 1.0 ##";
			memcpy(DataFile_header, head.c_str(), 50);
			int i = 2;
		};
		~DataFile() {};
		void set_mainName(string name) {
			mainName = name;
		}
		void set_path(string path) {
			_path = path;
#if defined(_WIN32)
			_mkdir(_path.c_str());
#elif defined(__linux__)
			mkdir(_path.c_str(), 0777);
			// #else more?
#endif
		}
		bool write_dataFile(FieldStorage_forPhaseNode& phaseMesh, string mark);
		bool read_dataFile(FieldStorage_forPhaseNode& phaseMesh, string name);
		DataFile& operator=(const DataFile& n);

		void write_PhaseNode_to_dataFile(ofstream &fout, PhaseNode &n);
		void read_PhaseNode_from_dataFile(fstream& fout, PhaseNode& n);
		string _path;
		string mainName;
		string format;
		char* DataFile_header;
	};
	inline DataFile& DataFile::operator=(const DataFile& n) {
		_path = n._path;
		mainName = n.mainName;
		format = n.format;
		memcpy(DataFile_header, n.DataFile_header, 50);
	}
	inline void DataFile::write_PhaseNode_to_dataFile(ofstream& fout, PhaseNode& n) {
		Data_phaseNode phaseNode;
		Data_phaseEntry phaseEntry;
		Data_nodeEntry nodeEntry;
		phaseNode.temperature = n.tempValues.temperature;
		phaseNode.phaseNum = n.size();
		phaseNode.conNum = n.x.size();
		phaseNode.custom_flag_num = n.customFlags.size();
		phaseNode.custom_value_num = n.customValues.size();
		phaseNode.custom_vec_num = n.customVec3s.size();
		fout.write((const char*)&phaseNode, sizeof(Data_phaseNode));
		for (auto c = n.x.begin(); c < n.x.end(); c++) {
			nodeEntry.index = c->index;
			nodeEntry.value = c->value;
			fout.write((const char*)&nodeEntry, sizeof(Data_nodeEntry));
		}
		for (int i = 0; i < phaseNode.custom_value_num; i++) {
			Data_custom_value value;
			value.index = n.customValues._double_box[i].index;
			value.value = n.customValues._double_box[i].value;
			value.value2 = 0.0;
			value.value3 = 0.0;
			value.flag = 0;
			fout.write((const char*)&value, sizeof(Data_custom_value));
		}
		for (int i = 0; i < phaseNode.custom_flag_num; i++) {
			Data_custom_value value;
			value.index = n.customFlags._int_box[i].index;
			value.flag = n.customFlags._int_box[i].value;
			value.value = 0.0;
			value.value2 = 0.0;
			value.value3 = 0.0;
			fout.write((const char*)&value, sizeof(Data_custom_value));
		}
		for (int i = 0; i < phaseNode.custom_vec_num; i++) {
			Data_custom_value value;
			value.index = n.customVec3s._vec_box[i].index;
			value.flag = 0;
			value.value = n.customVec3s._vec_box[i].vec[0];
			value.value2 = n.customVec3s._vec_box[i].vec[1];
			value.value3 = n.customVec3s._vec_box[i].vec[2];
			fout.write((const char*)&value, sizeof(Data_custom_value));
		}
		///< storage for _Phase
		for (auto p = n._Phase.begin(); p < n._Phase.end(); p++) {
			phaseEntry._flag = p->_flag;
			phaseEntry.index = p->index;
			phaseEntry.phaseProperty = p->phaseProperty;
			phaseEntry.phaseFraction = p->phaseFraction;
			phaseEntry.conNum = p->x.size();
			fout.write((const char*)&phaseEntry, sizeof(Data_phaseEntry));
			///< con
			for (auto c = p->x.begin(); c < p->x.end(); c++) {
				nodeEntry.index = c->index;
				nodeEntry.value = c->value;
				fout.write((const char*)&nodeEntry, sizeof(Data_nodeEntry));
			}
		}
	}
	inline void DataFile::read_PhaseNode_from_dataFile(fstream& fin, PhaseNode& n) {
		Data_phaseNode phaseNode;
		Data_phaseEntry phaseEntry;
		Data_nodeEntry nodeEntry;
		fin.read((char*)&phaseNode, sizeof(Data_phaseNode));
		n.tempValues.temperature = phaseNode.temperature;
		n._Phase.resize(phaseNode.phaseNum);
		for (int c_num = 0; c_num < phaseNode.conNum; c_num++) {
			fin.read((char*)&nodeEntry, sizeof(Data_nodeEntry));
			n.x.add_con(nodeEntry.index, nodeEntry.value);
		}
		for (int i = 0; i < phaseNode.custom_value_num; i++) {
			Data_custom_value cus_val;
			fin.read((char*)&cus_val, sizeof(Data_custom_value));
			n.customValues.add_double(cus_val.index, cus_val.value);
		}
		for (int i = 0; i < phaseNode.custom_flag_num; i++) {
			Data_custom_value cus_val;
			fin.read((char*)&cus_val, sizeof(Data_custom_value));
			n.customFlags.add_int(cus_val.index, cus_val.flag);
		}
		for (int i = 0; i < phaseNode.custom_vec_num; i++) {
			Data_custom_value cus_val;
			fin.read((char*)&cus_val, sizeof(Data_custom_value));
			n.customVec3s.add_vec(cus_val.index, Vector3(cus_val.value, cus_val.value2, cus_val.value3));
		}
		///< storage for _Phase
		for (auto p = n._Phase.begin(); p < n._Phase.end(); p++) {
			fin.read((char*)&phaseEntry, sizeof(Data_phaseEntry));
			p->_flag = phaseEntry._flag;
			p->index = phaseEntry.index;
			p->phaseProperty = phaseEntry.phaseProperty;
			p->phaseFraction = phaseEntry.phaseFraction;
			///< con
			for (int c_num = 0; c_num < phaseEntry.conNum; c_num++) {
				fin.read((char*)&nodeEntry, sizeof(Data_nodeEntry));
				p->x.add_con(nodeEntry.index, nodeEntry.value);
			}
		}
	}
	inline bool DataFile::write_dataFile(FieldStorage_forPhaseNode& phaseMesh, string mark) {
		string fname;
		if (_path == "")
			fname = mainName + mark + format;
		else
			fname = _path + dirSeparator + mainName + mark + format;
		ofstream fout(fname, std::ios::binary);
		if (!fout) {
			cout << "Failed to write the data file!" << endl;
			fout.close();
			return false;
		}
		fout.write((const char*)DataFile_header, 50);
		{ ///< defined in sequence(same with read)
			Data_fieldStorage fieldStorage;
			fieldStorage.dx = phaseMesh.dx;
			fieldStorage.Nx = phaseMesh.limit_x;
			fieldStorage.Ny = phaseMesh.limit_y;
			fieldStorage.Nz = phaseMesh.limit_z;
			fieldStorage.x_bc = phaseMesh._BC.x_bc;
			fieldStorage.y_bc = phaseMesh._BC.y_bc;
			fieldStorage.z_bc = phaseMesh._BC.z_bc;
			fout.write((const char*)&fieldStorage, sizeof(Data_fieldStorage));
			///< storage for _mesh
			for (auto n = phaseMesh._mesh.begin(); n < phaseMesh._mesh.end(); n++)
				write_PhaseNode_to_dataFile(fout, (*n));
		}
		fout.close();
		return true;
	}
	inline bool DataFile::read_dataFile(FieldStorage_forPhaseNode& phaseMesh, string name) {
		Data_fieldStorage fieldStorage;
		char header[50];
		std::fstream fin(name, std::ios::binary | ios::in);
		if (!fin) {
			cout << "Failed to read the aim file!" << endl;
			fin.close();
			return false;
		}
		fin.read((char*)header, 50);
		if (strcmp(header, DataFile_header)) {
			cout << "File format error, check content or version!" << endl;
			fin.close();
			return false;
		}
		{  ///< defined in sequence(same with write)
			fin.read((char*)&fieldStorage, sizeof(Data_fieldStorage));
			phaseMesh.init(fieldStorage.Nx, fieldStorage.Ny, fieldStorage.Nz, fieldStorage.dx, 
				fieldStorage.x_bc, fieldStorage.y_bc, fieldStorage.z_bc);

			for (auto n = phaseMesh._mesh.begin(); n < phaseMesh._mesh.end(); n++)
				read_PhaseNode_from_dataFile(fin, (*n));
		}
		fin.close();
		return true;
	}

}