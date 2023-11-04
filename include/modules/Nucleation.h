#pragma once
#include"../baseTools/baseTools.h"
#include"Information.h"
#include"Thermodynamics.h"
#include"Kinetics.h"
#include"Mechanics.h"
using namespace std;
namespace pf {
	class Nucleation
	{
	public:
		Nucleation() {};
		~Nucleation() {};

		void init(FieldStorage_forPhaseNode& _simulationField, Information& _information, Thermodynamics& _thermodynamics, 
			Kinetics& _kinetics, Mechanics& _mechanics);
		bool add_NewPhase_byPhaseIndex(int phaseIndex, int phaseProperty);
		bool nucleation(int istep, bool normalize_phi = true);
		void clear();

	private:
		void search_region_possible_to_nucleate(int x, int y, int z, int pIndex, vector<Point>& points);
		
		NP_FieldStorage nucleation_mesh;
		NucleationBox nucleationBox;
		///< connect
		FieldStorage_forPhaseNode* simulationField;
		Information* information;
		Thermodynamics* thermodynamics;
		Kinetics* kinetics;
		Mechanics* mechanics;

		bool definiteNucleation(int istep);
		vector<int> get_nucleus_by_condition(int istep);
		void get_nucleus_by_outerFunction(int istep);
	};
	inline bool Nucleation::nucleation(int istep, bool normalize_phi) {
		vector<int> phase_index;
		information->dynamicCollection.is_Nucleation = false;
		information->dynamicCollection.new_phase_index.clear();

		if (information->nucleationBox.nucleation_property == ConditionalNucleation && istep % information->nucleationBox.nucleation_step == 0) {
			information->statistics_information_in_phaseMesh(*simulationField);
			phase_index = get_nucleus_by_condition(istep);
		}
		else if (information->nucleationBox.nucleation_property == UserDefined_Nucleated
			&& istep % information->nucleationBox.nucleation_step == 0) {
			information->statistics_information_in_phaseMesh(*simulationField);
			get_nucleus_by_outerFunction(istep);
		}
		if (normalize_phi)
			information->normalize_phi_in_mesh(*simulationField);
		information->dynamicCollection.is_Nucleation = definiteNucleation(istep);
		if (phase_index.size() != 0) {
			information->dynamicCollection.is_possible_Nucleation = true;
			information->writer.write_nucleation_phase_field(nucleation_mesh, istep, phase_index, information->settings.file_settings);
			cout << "> Nucleation module: Phases that may nucleate are marked." << endl << endl;
		}
		return information->dynamicCollection.is_Nucleation;
	}
}