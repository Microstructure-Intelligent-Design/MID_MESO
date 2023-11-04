#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<cstdio>
#include<omp.h>
#include<vector>
#include<cmath>
#include<iostream>
#include <cstring>
#include <sstream>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include<algorithm>
#include <initializer_list>
#include<complex>
#include  <stdio.h>
#include<float.h>
#ifdef _WIN32
#define SYS_PROGRAM_STOP while(true){getchar();};
#define IS_NAN(a) _isnan(a)
#include  <tchar.h>
#include  <shlobj.h>
#include  <direct.h>
#include  <Windows.h>
#include  <io.h>
#include<conio.h>
const std::string dirSeparator = "\\";                                      //< Windows style directory separator
#else
#include <sys/io.h>
#include <sys/stat.h>
#define SYS_PROGRAM_STOP std::exit(0);
#define IS_NAN(a) __isnan(a)
const std::string dirSeparator = "/";                                       //< Unix/Linux style directory separator
#endif;

#define SYS_EPSILON 1e-6
#define Tangent_Max_Number double(1e10)
#define CON_EPSILON 1e-4
#define X_EPSILON 1e-6
#define Simulation_Num_Cut_Off 1e-3
#define Is_Equality(a, b) pf::isTwoNumEquality(a, b)
#define PI double(3.1415926535897932)
#define ElementaryChargeQuantity double(1.6021892e-19)   // C
#define AvogadroConstant double(6.02214076e23)  // mol -1
#define NA double(6.02214076e23)  // mol -1
#define FaradayConstant 96485.333  // C/mol
#define VacuumPermeability_Magnetismus (PI*4e-7) // 
#define AngleToRadians(angle) double(angle/180.0*PI)

#define VA int(-1)
#define VIRTUAL_PHASE int(-1)
#define GROW_PHASE int(-2)
#define SHRINK_PHASE int(-3)
#define SOLVENT_NONE int(-1)
#define VA_str "VA"

#define RAND_init (srand((unsigned int)time(NULL)))
#define RAND_0_1 (((double)std::rand())/RAND_MAX)

namespace pf {
	using namespace std;

	enum Int_Gradient { Steinbach_1996, Steinbach_1999, Steinbach_G2009, Int_GCustom };
	enum Int_Potential { Nestler_Well, Nestler_Obstacle, Steinbach_P2009, Int_PCustom };
	enum ThermodynamicModel { SublatticeModel, SolutionModel, AssociateLiquidModel, Solution_LinearCompoundModel, Pure_Substance };
	enum PhaseFluxModel { IntDiff_PotentialGrad, IntDiff_ConGrad };
	enum DifferenceMethod { FIVE_POINT, NINE_POINT };
	enum EffectiveElasticConstantsModel { EEC_IngoSteinbach, EEC_Khachaturyan, EEC_Custom };
	enum ConcentrationEvolutionEquation { PhaseConcentrationModel, TotalConcentrationModel };
	enum NumericalLimit { DownLimit, UpLimit };
	enum Boundary { UP_X, UP_Y, UP_Z, DOWN_X, DOWN_Y, DOWN_Z };
	enum BoundaryCondition { FIXED, PERIODIC, ADIABATIC };
	enum Dimension { One_Dimension, Two_Dimension, Three_Dimension };
	enum Direction { x_up, x_down, y_up, y_down, z_up, z_down, num_of_directons };
	enum PHASES_GROW_SHRINK_TYPE {NONE, COMPLETE, BULK_LOCK, HALF_INT};
	// laplace(LHS) = RHS
	enum DEFULT_SOLVER_VALUE_INDEX { SOLVER_ALLEN_CAHN = 10000, SOLVER_CAHN_HILLIARD = 20000, BUFF_FOR_SOLVER_CAHN_HILLIARD = 30000 };

	namespace materials {
		enum Element {
			ElementBottom, H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni,
			Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce,
			Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra,
			Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Uun, Uuu, Uub, ElementTop
		};
		const vector<std::string> Element_Name = { "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
		"Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
		"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
		"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Uun", "Uuu", "Uub", "" };
	};

	static bool isTwoNumEquality(double a, double b) {
		if (std::fabs(a - b) < SYS_EPSILON)
			return true;
		else
			return false;
	}

	static int double_to_int(double a) {
		if ((a - int(a)) > 0.5)
			return int(a) + 1;
		else
			return int(a);
	}

	struct Point {
		double x;
		double y;
		double z;
		Point() {
			x = 0;
			y = 0;
			z = 0;
		}
		Point(double _x, double _y, double _z) {
			x = _x;
			y = _y;
			z = _z;
		}
		void set(double _x, double _y, double _z) {
			x = _x;
			y = _y;
			z = _z;
		}
		void do_boundary(BoundaryCondition x_bc, BoundaryCondition y_bc, BoundaryCondition z_bc, int x_limit, int y_limit, int z_limit) {
			if (x_bc == BoundaryCondition::ADIABATIC && x < 0)
				x = 0;
			else if (x_bc == BoundaryCondition::ADIABATIC && x >= x_limit)
				x = x_limit - 1;
			else if (x_bc == BoundaryCondition::FIXED && x < 0)
				x = -1;
			else if (x_bc == BoundaryCondition::FIXED && x >= x_limit)
				x = x_limit;
			else if (x_bc == BoundaryCondition::PERIODIC && x < 0)
				x += x_limit;
			else if (x_bc == BoundaryCondition::PERIODIC && x >= x_limit)
				x -= x_limit;
			else if (y_bc == BoundaryCondition::ADIABATIC && y < 0)
				y = 0;
			else if (y_bc == BoundaryCondition::ADIABATIC && y >= y_limit)
				y = y_limit - 1;
			else if (y_bc == BoundaryCondition::FIXED && y < 0)
				y = -1;
			else if (y_bc == BoundaryCondition::FIXED && y >= y_limit)
				y = y_limit;
			else if (y_bc == BoundaryCondition::PERIODIC && y < 0)
				y += y_limit;
			else if (y_bc == BoundaryCondition::PERIODIC && y >= y_limit)
				y -= y_limit;
			else if (z_bc == BoundaryCondition::ADIABATIC && z < 0)
				z = 0;
			else if (z_bc == BoundaryCondition::ADIABATIC && z >= z_limit)
				z = z_limit - 1;
			else if (z_bc == BoundaryCondition::FIXED && z < 0)
				z = -1;
			else if (z_bc == BoundaryCondition::FIXED && z >= z_limit)
				z = z_limit;
			else if (z_bc == BoundaryCondition::PERIODIC && z < 0)
				z += z_limit;
			else if (z_bc == BoundaryCondition::PERIODIC && z >= z_limit)
				z -= z_limit;
			else
				return;
			return do_boundary(x_bc, y_bc, z_bc, x_limit, y_limit, z_limit);
		}
		Point& operator=(const Point& n) {
			this->x = n.x;
			this->y = n.y;
			this->z = n.z;
			return *this;
		}
		Point operator+(const Point& n) {
			Point re;
			re.x = this->x + n.x;
			re.y = this->y + n.y;
			re.z = this->z + n.z;
			return re;
		}
		Point operator-(const Point& n) {
			Point re;
			re.x = this->x - n.x;
			re.y = this->y - n.y;
			re.z = this->z - n.z;
			return re;
		}
		Point operator*(const double n) {
			Point re;
			re.x = this->x * n;
			re.y = this->y * n;
			re.z = this->z * n;
			return re;
		}
		Point operator/(const double n) {
			Point re;
			re.x = this->x / n;
			re.y = this->y / n;
			re.z = this->z / n;
			return re;
		}
	};

	struct Point2D {
		double x;
		double y;
		Point2D() {
			x = 0;
			y = 0;
		}
		Point2D(double _x, double _y) {
			x = _x;
			y = _y;
		}
		void set(double _x, double _y) {
			x = _x;
			y = _y;
		}
		Point2D& operator=(const Point2D& n) {
			this->x = n.x;
			this->y = n.y;
			return *this;
		}
	};

	struct LinearCompoundModel {
		double m;
		double x;
		double scale_max;
		void set_parameters(double _m, double _x) {
			m = _m;
			x = _x;
		}
		LinearCompoundModel() {
			m = 1.0;
			x = 1.0;
			scale_max = 0.0;
		}
	};
	struct bool_elem {
		int index;
		bool value;
		bool_elem() {
			index = 0;
			value = 0;
		}
		bool_elem& operator=(const bool_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct bool_box {
		vector<bool_elem> _bool_box;
		typedef std::vector<bool_elem>::iterator iterator;
		typedef std::vector<bool_elem>::const_iterator citerator;
		iterator  begin() { return _bool_box.begin(); };
		iterator  end() { return _bool_box.end(); };
		bool& operator[](const int index) {
			for (auto i = _bool_box.begin(); i < _bool_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "_bool_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		bool_box& operator=(const bool_box& n) {
			_bool_box = n._bool_box;
			return *this;
		}
		void add_bool(int _index, bool _value) {
			for (auto i = _bool_box.begin(); i < _bool_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			bool_elem elem;
			elem.index = _index;
			elem.value = _value;
			_bool_box.push_back(elem);
		}
		void clear() {
			_bool_box.clear();
		}
		int size() {
			return int(_bool_box.size());
		}
	};
	struct int_elem {
		int index;
		int value;
		int_elem() {
			index = 0;
			value = 0;
		}
		int_elem& operator=(const int_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct int_box {
		vector<int_elem> _int_box;
		typedef std::vector<int_elem>::iterator iterator;
		typedef std::vector<int_elem>::const_iterator citerator;
		iterator  begin() { return _int_box.begin(); };
		iterator  end() { return _int_box.end(); };
		int& operator[](const int index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "int_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		int_box& operator=(const int_box& n) {
			_int_box = n._int_box;
			return *this;
		}
		void add_int(int _index, int _value) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			int_elem elem;
			elem.index = _index;
			elem.value = _value;
			_int_box.push_back(elem);
		}
		void clear() {
			_int_box.clear();
		}
		int size() {
			return int(_int_box.size());
		}
	};
	struct double_elem {
		int index;
		double value;
		double_elem() {
			index = 0;
			value = 0.0;
		}
		double_elem& operator=(const double_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct double_box {
		vector<double_elem> _double_box;
		typedef std::vector<double_elem>::iterator iterator;
		typedef std::vector<double_elem>::const_iterator citerator;
		iterator  begin() { return _double_box.begin(); };
		iterator  end() { return _double_box.end(); };
		double& operator[](const int index) {
			for (auto i = _double_box.begin(); i < _double_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "double_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		double_box& operator=(const double_box& n) {
			_double_box = n._double_box;
			return *this;
		}
		void add_double(int _index, double _value) {
			for (auto i = _double_box.begin(); i < _double_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			double_elem elem;
			elem.index = _index;
			elem.value = _value;
			_double_box.push_back(elem);
		}
		void clear() {
			_double_box.clear();
		}
		int size() {
			return int(_double_box.size());
		}
	};
	struct pair_flag {
		int index1;
		int index2;
		int value;
		pair_flag() {
			index1 = 0;
			index2 = 0;
			value = 0;
		}
		pair_flag& operator=(const pair_flag& n) {
			index1 = n.index1;
			index2 = n.index2;
			value = n.value;
			return *this;
		}
	};
	struct pair_flag_box {
		vector<pair_flag> _flag_box;
		typedef std::vector<pair_flag>::iterator iterator;
		typedef std::vector<pair_flag>::const_iterator citerator;
		iterator  begin() { return _flag_box.begin(); };
		iterator  end() { return _flag_box.end(); };
		int& operator()(const int index1, const int index2) {
			for (auto i = _flag_box.begin(); i < _flag_box.end(); ++i) {
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) return i->value;
			}
			cout << "int_box error, can't find the value indexs : " << index1 << ", " << index2 << endl;
			SYS_PROGRAM_STOP;
		}
		pair_flag_box& operator=(const pair_flag_box& n) {
			_flag_box = n._flag_box;
		}
		void set_flag(int index1, int index2, int flag) {
			for (auto i = _flag_box.begin(); i < _flag_box.end(); ++i)
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) {
					i->value = flag;
					return;
				}
			pair_flag elem;
			elem.index1 = index1;
			elem.index2 = index2;
			elem.value = flag;
			_flag_box.push_back(elem);
		}
		void clear() {
			_flag_box.clear();
		}
		int size() {
			return int(_flag_box.size());
		}
	};
	struct string_elem {
		int index;
		string value;
		string_elem() {
			index = 0;
			value = "";
		}
		string_elem& operator=(const string_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};
	struct string_box {
		vector<string_elem> _string_box;
		typedef std::vector<string_elem>::iterator iterator;
		typedef std::vector<string_elem>::const_iterator citerator;
		iterator  begin() { return _string_box.begin(); };
		iterator  end() { return _string_box.end(); };
		string& operator[](const int index) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			cout << "string_box error, can't find the value index : " << index << endl;
			SYS_PROGRAM_STOP;
		}
		string_box& operator=(const string_box& n) {
			_string_box = n._string_box;
			return *this;
		}
		void add_string(int _index, string _value) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			string_elem elem;
			elem.index = _index;
			elem.value = _value;
			_string_box.push_back(elem);
		}
		void clear() {
			_string_box.clear();
		}
		int size() {
			return int(_string_box.size());
		}
	};

#ifdef _WIN32
	static std::string TCHAR2STRING(TCHAR* STR)
	{
		int iLen = WideCharToMultiByte(CP_ACP, 0, STR, -1, NULL, 0, NULL, NULL);
		char* chRtn = new char[iLen * sizeof(char)];
		WideCharToMultiByte(CP_ACP, 0, STR, -1, chRtn, iLen, NULL, NULL);
		std::string str(chRtn);
		delete chRtn;
		return str;
	}
	static bool SelectFolderPath(std::string& folderName)
	{
		TCHAR szBuffer[MAX_PATH] = { 0 };
		BROWSEINFO bi;
		ZeroMemory(&bi, sizeof(BROWSEINFO));
		bi.hwndOwner = GetForegroundWindow();
		bi.pszDisplayName = szBuffer;
		bi.pidlRoot = NULL;
		bi.lpszTitle = _T("从下面选文件:");
		bi.ulFlags = BIF_BROWSEINCLUDEFILES;;
		LPITEMIDLIST idl = SHBrowseForFolder(&bi);
		if (NULL == idl){
			return false;
		}
		SHGetPathFromIDList(idl, szBuffer);
		folderName = TCHAR2STRING(szBuffer);
		return true;
		
	}
#endif;
	static char get_char_not_show() {
		char c;
#ifdef _WIN32
		c = _getch();
#else
		system("stty -echo");
		c = getchar();
		system("stty echo");
#endif
		return c;
	}
	static string get_string_from_consol(bool is_show = true, char replace_char = '*') {
		string str = "";
		int i = 0;
		if (is_show) {
			while (true) {
				char ch = getchar();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
			}
		}
		else {
			while (true) {
				char ch = get_char_not_show();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
				putchar(replace_char);
			}
		}
		putchar('\n');
		return str;
	}
}