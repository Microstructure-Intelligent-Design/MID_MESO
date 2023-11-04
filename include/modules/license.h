#pragma once
#include "../baseTools/baseTools.h"
namespace pf {
#ifdef _WIN32
	const string project_path = "include/license.dat";
#else
	const string project_path = "../../include/license.dat";
#endif
	class License
	{
	public:
		License();
		~License();
		bool license(string licensefile = project_path);
		void write_current_time_into_license(string licensefile = project_path);
	private:
		bool generate_a_license_dat(string licensefile);
		bool read_license(string name, ints_time& time1, ints_time& time2, ints_time& time3);
		bool write_license(string name, ints_time& time1, ints_time& time2, ints_time& time3);
		bool second_protection_read(string name, ints_time& time);
		bool second_protection_write(string name, ints_time& time);
	};
}