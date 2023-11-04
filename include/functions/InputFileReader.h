#pragma once
#include "../baseTools/sysTool.h"
namespace pf {
	using namespace std;
	class InputFileReader
	{
	public:
		InputFileReader() {
			_split = ' ';
			input_file.clear();
		}
		InputFileReader(string input_file_name, bool debug = false, int file_max_lines = 100, char split = ' ') {
			init(input_file_name);
		}
		void init(string input_file_name, bool debug = false, int file_max_lines = 100, char split = ' ') {
			input_file.clear();
			_split = split;
			std::fstream fin(input_file_name, std::ios::in);
			if (debug)
				cout << "-----------------------------------------------------------------------------------------" << endl;
			if (!fin) {
				cout << "> Failed to read the input_file:" << input_file_name << ", use default value 0 (int, double, bool), \"\" (string)" << endl;
			}
			else {
				int line = 1;
				string strline;
				while (getline(fin, strline) && line <= file_max_lines) {
					input_file.add_string(line, strline);
					line++;
				}
			}
			fin.close();
			get_valid_words();
			if (debug)
				debug_input_file_lines();
			// Remove duplicate definitions
			if(debug)
				cout << "> Reading data from input file! Following is [keyword] + [Value], the valid keyword works." << endl;
			for (auto vec1 = _valid_words.begin(); vec1 < _valid_words.end(); vec1++) {
				if (debug)
					cout << "Valid   :		[" << (*vec1)[0] << "] + [" << (*vec1)[1] << "]" << endl;
				for (auto vec2 = vec1 + 1; vec2 < _valid_words.end();) {
					if ((*vec1)[0].compare((*vec2)[0]) == 0) {
						if (debug)
							cout << "in-Valid:		[" << (*vec2)[0] << "] + [" << (*vec2)[1] << "]" << endl;
						vec2 = _valid_words.erase(vec2);
					}
					else
						++vec2;
				}
			}
			if (debug)
				cout << "-----------------------------------------------------------------------------------------" << endl;
		}
		~InputFileReader() {
			_split = ' ';
			input_file.clear();
		}
		bool read_int_value(string value_name, int& int_value) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					try
					{
						int_value = stoi((*vec)[1]);
						cout << "> Keyword: " << value_name << " is defined as " << int_value << "." << endl;
					}
					catch (double)
					{
						cout << "> Input file error! the value of Keyword: " << (*vec)[0] << " cannot be translate to int value." << endl;
						SYS_PROGRAM_STOP;
					}
					return true;
				}
			cout << "> [DEFAULT] Keyword: " << value_name << " is defined as " << int_value << "." << endl;
			return false;
		}
		bool read_double_value(string value_name, double& double_value) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					try
					{
						double_value = stod((*vec)[1]);
						cout << "> Keyword: " << value_name << " is defined as " << double_value << "." << endl;
					}
					catch (double)
					{
						cout << "> Input file error! the value of Keyword: " << (*vec)[0] << " cannot be translate to double value." << endl;
						SYS_PROGRAM_STOP;
					}
					return true;
				}
			cout << "> [DEFAULT] Keyword: " << value_name << " is defined as " << double_value << "." << endl;
			return false;
		}
		bool read_bool_value(string value_name, bool& bool_value) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					bool_value = false;
					if ((*vec)[1].compare("true") == 0 || (*vec)[1].compare("TRUE") == 0 || (*vec)[1].compare("1") == 0)
						bool_value = true;
					if(bool_value)
						cout << "> Keyword: " << value_name << " is defined as TRUE." << endl;
					else
						cout << "> Keyword: " << value_name << " is defined as FALSE." << endl;
					return true;
				}
			if (bool_value)
				cout << "> [DEFAULT] Keyword: " << value_name << " is defined as TRUE." << endl;
			else
				cout << "> [DEFAULT] Keyword: " << value_name << " is defined as FALSE." << endl;
			return false;
		}
		bool read_string_value(string value_name, string& string_value) {
			for (auto vec = _valid_words.begin(); vec < _valid_words.end(); vec++)
				if (value_name.compare((*vec)[0]) == 0) {
					string_value = (*vec)[1];
					cout << "> Keyword: " << value_name << " is defined as " << string_value << "." << endl;
					return true;
				}
			cout << "> [DEFAULT] Keyword: " << value_name << " is defined as " << string_value << "." << endl;
			return false;
		}
		void debug_input_file_lines() {
			vector<int> valie_lines;
			string valid_words = "<valid>		", invalid_words = "<in-valid>	", note_words = "<note>		";
			cout << "======================================= D E B U G =======================================" << endl;
			cout << "LINE	PROPERTY	|CONTENT" << endl;
			cout << "-----------------------------------------------------------------------------------------" << endl;
			for (int line = 1; line <= input_file.size(); line++) {
				string out = "[" + to_string(line) + "]	", str = input_file[line], equal = "=";
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else
						++s;
				}
				str = "|" + str;
				if (head == '#')
					out += note_words + str;
				else if (head == ' ' || head == '	' || head == '\n' || vec3.size() != 3)
					out += invalid_words + str;
				else if (equal.compare(vec3[1]) != 0)
					out += invalid_words + str;
				else {
					out += valid_words + str;
					valie_lines.push_back(line);
				}
				cout << out << endl;
			}
			cout << "-----------------------------------------------------------------------------------------" << endl;
			if (valie_lines.size() != _valid_words.size()) {
				cout << "Check valid words error!" << endl;
				SYS_PROGRAM_STOP;
			}
			for (int index = 0; index < valie_lines.size(); index++) {
				string out = "[" + to_string(valie_lines[index]) + "]	<valid>		|{\"" + _valid_words[index][0] + "\"}, {\"=\"}, {\"" + _valid_words[index][1] + "\"}";
				cout << out << endl;
			}
			cout << "=========================================================================================" << endl;
		}
		void get_valid_words() {
			for (int line = 1; line <= input_file.size(); line++) {
				string str = input_file[line], equal = "=";
				char head;
				head = get_first_character_of_line(str);
				vector<string> vec3 = split_string(str, _split);
				for (auto s = vec3.begin(); s < vec3.end();) {
					if ((*s).compare("") == 0)
						s = vec3.erase(s);
					else
						++s;
				}
				if (head == '#' || head == ' ' || head == '	' || head == '\n' || vec3.size() != 3)
					continue;
				else if (equal.compare(vec3[1]) != 0)
					continue;
				else {
					vector<string> valie_line;
					valie_line.push_back(vec3[0]);
					valie_line.push_back(vec3[2]);
					_valid_words.push_back(valie_line);
				}
			}
		}
		string_box get_whole_file_strings() {
			return input_file;
		}
	private:
		string_box input_file;
		// [lines]{name, value}
		vector<vector<string>> _valid_words;
		char _split;
		char get_first_character_of_line(string line) {
			if (line.compare("") != 0)
				return line.at(0);
			else
				return '0';
		}
		vector<string> split_string(string str, const char split) {
			vector<string> res;
			istringstream iss(str);
			string token;
			while (getline(iss, token, split)) {
				res.push_back(token);
			}
			return res;
		}
	};




}