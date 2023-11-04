#pragma once
#include "sysTool.h"
#include "NormalNode.h"
#include "BoundaryCondition.h"
#include "PairValueNode.h"
#include "SublatticNode.h"
#define X_MIN 0
#define Y_MIN 0
#define Z_MIN 0

using namespace std;
namespace pf {
	static int Index(int x, int y, int z, int limit_x, int limit_y, int limit_z) {
		return x + y * limit_x + z * limit_x * limit_y;
	}
	class FieldStorage_forPhaseNode
	{
	public:
		FieldStorage_forPhaseNode() {};
		~FieldStorage_forPhaseNode() {
			try {
				_BC.free();
				_mesh.clear();
			}
			catch(...){
				std::cout << "Destroy FieldStorage_forPhaseNode failed, problems occur." << endl;
			}
		};

		PhaseNode& operator()(const int x, const int y, const int z) ;                                       //< Index operator for accessing the n's field value
		FieldStorage_forPhaseNode&  operator=(const FieldStorage_forPhaseNode& n);

		Flag currentFlag(PhaseNode& node, int phaseIndex);      //<  seek where the node is (set in storage)
		bool is_interface(PhaseNode& node, int phaseIndex);
		void upgrade(PhaseNode& node, int phaseIndex);
		void downgrade(PhaseNode& node, int phaseIndex);
		void init(int Nx, int Ny, int Nz, double _dx, BoundaryCondition bc_x, BoundaryCondition bc_y, BoundaryCondition bc_z) {
			//<maxNum 
			_mesh.resize(Nx*Ny*Nz);
			limit_x = Nx;
			limit_y = Ny;
			limit_z = Nz;
			dx = _dx;
			_BC.init(bc_x, bc_y, bc_z, Nx, Ny, Nz);
			_dimention = Three_Dimension;
			if (Nx == 1 || Ny == 1 || Nz == 1)
				_dimention = Two_Dimension;
			if ((Nx == 1 && Ny == 1) || (Nx == 1 && Nz == 1) || (Ny == 1 && Nz == 1))
				_dimention = One_Dimension;
			for(int z = 0; z < Nz; z++)
				for (int y = 0; y < Ny; y++)
					for (int x = 0; x < Nx; x++) {
						PhaseNode* _up_x;
						PhaseNode* _down_x;
						PhaseNode* _up_y;
						PhaseNode* _down_y;
						PhaseNode* _up_z;
						PhaseNode* _down_z;
						if (bc_x == BoundaryCondition::ADIABATIC) {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
							_down_x = &_mesh[down];
						}
						else if (bc_x == BoundaryCondition::FIXED) {
							if (x != 0 && x != (Nx - 1)) {
								_up_x = &_mesh[(x + 1) + y * limit_x + z * limit_x * limit_y];
								_down_x = &_mesh[(x - 1) + y * limit_x + z * limit_x * limit_y];
							}
							else if (x == 0 && x != (Nx - 1)) {
								_up_x = &_mesh[(x + 1) + y * limit_x + z * limit_x * limit_y];
								_down_x = &_BC(Boundary::DOWN_X, y, z);
							}
							else if (x == (Nx - 1) && x != 0) {
								_up_x = &_BC(Boundary::UP_X, y, z);
								_down_x = &_mesh[(x - 1) + y * limit_x + z * limit_x * limit_y];
							}
							else {
								_up_x = &_BC(Boundary::UP_X, y, z);
								_down_x = &_BC(Boundary::DOWN_X, y, z);
							}
						}
						else {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = (Nx - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = 0 + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
							_down_x = &_mesh[down];
						}
						if (bc_y == BoundaryCondition::ADIABATIC) {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
							_down_y = &_mesh[down];
						}
						else if (bc_y == BoundaryCondition::FIXED) {
							if (y != 0 && y != (Ny - 1)) {
								_up_y = &_mesh[x + (y + 1) * limit_x + z * limit_x * limit_y];
								_down_y = &_mesh[x + (y - 1) * limit_x + z * limit_x * limit_y];
							}
							else if (y == 0 && y != (Ny - 1)) {
								_up_y = &_mesh[x + (y + 1) * limit_x + z * limit_x * limit_y];
								_down_y = &_BC(Boundary::DOWN_Y, x, z);
							}
							else if (y == (Ny - 1) && y != 0) {
								_up_y = &_BC(Boundary::UP_Y, x, z);
								_down_y = &_mesh[x + (y - 1) * limit_x + z * limit_x * limit_y];
							}
							else {
								_up_y = &_BC(Boundary::UP_Y, x, z);
								_down_y = &_BC(Boundary::DOWN_Y, x, z);
							}
						}
						else {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + (Ny - 1) * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + 0 * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
							_down_y = &_mesh[down];
						}
						if (bc_z == BoundaryCondition::ADIABATIC) {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_z = &_mesh[up];
							_down_z = &_mesh[down];
						}
						else if (bc_z == BoundaryCondition::FIXED) {
							if (z != 0 && z != (Nz - 1)) {
								_up_z = &_mesh[x + y * limit_x + (z + 1) * limit_x * limit_y];
								_down_z = &_mesh[x + y * limit_x + (z - 1) * limit_x * limit_y];
							}
							else if (z == 0 && z != (Nz - 1)) {
								_up_z = &_mesh[x + y * limit_x + (z + 1) * limit_x * limit_y];
								_down_z = &_BC(Boundary::DOWN_Z, x, y);
							}
							else if (z == (Nz - 1) && z != 0) {
								_up_z = &_BC(Boundary::UP_Z, x, y);
								_down_z = &_mesh[x + y * limit_x + (z - 1) * limit_x * limit_y];
							}
							else {
								_up_z = &_BC(Boundary::UP_Z, x, y);
								_down_z = &_BC(Boundary::DOWN_Z, x, y);
							}
						}
						else {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + (Nz - 1) * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + 0 * limit_x * limit_y;
							_up_z = &_mesh[up];
							_down_z = &_mesh[down];
						}
						_mesh[x + y * limit_x + z * limit_x * limit_y].connect_mesh(x, y, z, *_up_x, *_down_x, *_up_y, *_down_y, *_up_z, *_down_z);
					}
		}
		void add_FixedBoundaryInf(vector<vector<PhaseNode>> bc, int phaseIndex) {
			/*  vector<vector<PhaseEntry>> bc;    bc[Direction][i + j * limit_x/y][PhaseIndex]   */

			if (_BC.x_bc == FIXED) {
				try {
					for (int i = 0; i < limit_y; i++)
						for (int j = 0; j < limit_z; j++) {
							PhaseEntry upX, downX;
							int index = Index_fix_bc_axisX(i, j, limit_y);
							upX   = bc[Boundary::UP_X  ][index][phaseIndex];
							downX = bc[Boundary::DOWN_X][index][phaseIndex];
							upX.potential = upX.x;
							downX.potential = downX.x;
							for (auto x1 = upX.x.begin(); x1 < upX.x.end(); x1++)
								for (auto x2 = upX.x.begin(); x2 < upX.x.end(); x2++)
									upX.kinetics_coeff.set(x1->index, x2->index, 0.0);
							for (auto x1 = downX.x.begin(); x1 < downX.x.end(); x1++)
								for (auto x2 = downX.x.begin(); x2 < downX.x.end(); x2++)
									downX.kinetics_coeff.set(x1->index, x2->index, 0.0);
							_BC.up_x  [index].push_back(upX);
							_BC.down_x[index].push_back(downX);
							_BC.up_x  [index].tempValues     = bc[Boundary::UP_X  ][index].tempValues;
							_BC.down_x[index].tempValues     = bc[Boundary::DOWN_X][index].tempValues;
							_BC.up_x  [index].electricValues = bc[Boundary::UP_X  ][index].electricValues;
							_BC.down_x[index].electricValues = bc[Boundary::DOWN_X][index].electricValues;
							_BC.up_x  [index].magneticValues = bc[Boundary::UP_X  ][index].magneticValues;
							_BC.down_x[index].magneticValues = bc[Boundary::DOWN_X][index].magneticValues;
							_BC.up_x  [index].customValues   = bc[Boundary::UP_X  ][index].customValues;
							_BC.down_x[index].customValues   = bc[Boundary::DOWN_X][index].customValues;
							_BC.up_x  [index].x              = bc[Boundary::UP_X  ][index].x;
							_BC.down_x[index].x              = bc[Boundary::DOWN_X][index].x;
						}
				}
				catch (...) {
					std::cout << "Definition of X_BoundaryConditionInf can't be used by Program, Please Check." << endl;
				}
			}
			if (_BC.y_bc == FIXED) {
				try {
					for (int i = 0; i < limit_x; i++)
						for (int j = 0; j < limit_z; j++) {
							PhaseEntry upY, downY;
							int index = Index_fix_bc_axisY(i, j, limit_x);
							upY   = bc[Boundary::UP_Y  ][index][phaseIndex];
							downY = bc[Boundary::DOWN_Y][index][phaseIndex];
							upY.potential = upY.x;
							downY.potential = downY.x;
							for (auto x1 = upY.x.begin(); x1 < upY.x.end(); x1++)
								for (auto x2 = upY.x.begin(); x2 < upY.x.end(); x2++)
									upY.kinetics_coeff.set(x1->index, x2->index, 0.0);
							for (auto x1 = downY.x.begin(); x1 < downY.x.end(); x1++)
								for (auto x2 = downY.x.begin(); x2 < downY.x.end(); x2++)
									downY.kinetics_coeff.set(x1->index, x2->index, 0.0);
							_BC.up_y  [index].push_back(upY);
							_BC.down_y[index].push_back(downY);
							_BC.up_y  [index].tempValues     = bc[Boundary::UP_Y  ][index].tempValues;
							_BC.down_y[index].tempValues     = bc[Boundary::DOWN_Y][index].tempValues;
							_BC.up_y  [index].electricValues = bc[Boundary::UP_Y  ][index].electricValues;
							_BC.down_y[index].electricValues = bc[Boundary::DOWN_Y][index].electricValues;
							_BC.up_y  [index].magneticValues = bc[Boundary::UP_Y  ][index].magneticValues;
							_BC.down_y[index].magneticValues = bc[Boundary::DOWN_Y][index].magneticValues;
							_BC.up_y  [index].customValues   = bc[Boundary::UP_Y  ][index].customValues;
							_BC.down_y[index].customValues   = bc[Boundary::DOWN_Y][index].customValues;
							_BC.up_y  [index].x              = bc[Boundary::UP_Y  ][index].x;
							_BC.down_y[index].x              = bc[Boundary::DOWN_Y][index].x;
						}
				}
				catch (...) {
					std::cout << "Definition of Y_BoundaryConditionInf can't be used by Program, Please Check." << endl;
				}
			}
			if (_BC.z_bc == FIXED) {
				try {
					for (int i = 0; i < limit_x; i++)
						for (int j = 0; j < limit_y; j++) {
							PhaseEntry upZ, downZ;
							int index = Index_fix_bc_axisZ(i, j, limit_x);
							upZ   = bc[Boundary::UP_Z  ][index][phaseIndex];
							downZ = bc[Boundary::DOWN_Z][index][phaseIndex];
							upZ.potential = upZ.x;
							downZ.potential = downZ.x;
							for (auto x1 = upZ.x.begin(); x1 < upZ.x.end(); x1++)
								for (auto x2 = upZ.x.begin(); x2 < upZ.x.end(); x2++)
									upZ.kinetics_coeff.set(x1->index, x2->index, 0.0);
							for (auto x1 = downZ.x.begin(); x1 < downZ.x.end(); x1++)
								for (auto x2 = downZ.x.begin(); x2 < downZ.x.end(); x2++)
									downZ.kinetics_coeff.set(x1->index, x2->index, 0.0);
							_BC.up_z  [index].push_back(upZ);
							_BC.down_z[index].push_back(downZ);
							_BC.up_z  [index].tempValues     = bc[Boundary::UP_Z  ][index].tempValues;
							_BC.down_z[index].tempValues     = bc[Boundary::DOWN_Z][index].tempValues;
							_BC.up_z  [index].electricValues = bc[Boundary::UP_Z  ][index].electricValues;
							_BC.down_z[index].electricValues = bc[Boundary::DOWN_Z][index].electricValues;
							_BC.up_z  [index].magneticValues = bc[Boundary::UP_Z  ][index].magneticValues;
							_BC.down_z[index].magneticValues = bc[Boundary::DOWN_Z][index].magneticValues;
							_BC.up_z  [index].customValues   = bc[Boundary::UP_Z  ][index].customValues;
							_BC.down_z[index].customValues   = bc[Boundary::DOWN_Z][index].customValues;
							_BC.up_z  [index].x              = bc[Boundary::UP_Z  ][index].x;
							_BC.down_z[index].x              = bc[Boundary::DOWN_Z][index].x;
						}
				}
				catch (...) {
					std::cout << "Definition of Z_BoundaryConditionInf can't be used by Program, Please Check." << endl;
				}
			}
		}
		void free() {
			_BC.free();
			for (auto node = _mesh.begin(); node < _mesh.end(); node++)
				node->clear();
			_mesh.clear();
			limit_x = 0;
			limit_y = 0;
			limit_z = 0;
			dx = 0.0;
		}
		int size() {
			return int(_mesh.size());
		}

		// ******   Settings   ****** //
		BoundaryCondition_forPhaseNode _BC;
		int limit_x;
		int limit_y;
		int limit_z;
		double dx;
		Dimension _dimention;
		std::vector<PhaseNode> _mesh;
	};
	inline Flag FieldStorage_forPhaseNode::currentFlag(PhaseNode& node, int phaseIndex) {
		if (node[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off && node[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off))
			return pf_INTERFACE;
		if (node[phaseIndex].phaseFraction < Simulation_Num_Cut_Off) {
			if (node.get_neighbor_node(Direction::x_up)[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off
				|| node.get_neighbor_node(Direction::x_down)[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off
				|| node.get_neighbor_node(Direction::y_up)[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off
				|| node.get_neighbor_node(Direction::y_down)[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off
				|| node.get_neighbor_node(Direction::z_down)[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off
				|| node.get_neighbor_node(Direction::z_up)[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off)
				return pf_NEAR_INTERFACE;
		}
		else if (node[phaseIndex].phaseFraction > (1.0 - Simulation_Num_Cut_Off)) {
			if (node.get_neighbor_node(Direction::x_up)[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
				|| node.get_neighbor_node(Direction::x_down)[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
				|| node.get_neighbor_node(Direction::y_up)[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
				|| node.get_neighbor_node(Direction::y_down)[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
				|| node.get_neighbor_node(Direction::z_up)[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off)
				|| node.get_neighbor_node(Direction::z_down)[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off))
				return pf_NEAR_INTERFACE;
		}
		return pf_BULK;
	}
	inline void FieldStorage_forPhaseNode::upgrade(PhaseNode& node, int phaseIndex) {
		node[phaseIndex]._flag = pf_INTERFACE;
		if ((node._x - 1) >= 0 || _BC.x_bc == PERIODIC)
			if (node.get_neighbor_node(Direction::x_down)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::x_down)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._x + 1) < limit_x || _BC.x_bc == PERIODIC)
			if (node.get_neighbor_node(Direction::x_up)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::x_up)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._y - 1) >= 0 || _BC.y_bc == PERIODIC)
			if (node.get_neighbor_node(Direction::y_down)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::y_down)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._y + 1) < limit_y || _BC.y_bc == PERIODIC)
			if (node.get_neighbor_node(Direction::y_up)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::y_up)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._z - 1) >= 0 || _BC.z_bc == PERIODIC)
			if (node.get_neighbor_node(Direction::z_down)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::z_down)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		if ((node._z + 1) < limit_z || _BC.z_bc == PERIODIC)
			if (node.get_neighbor_node(Direction::z_up)[phaseIndex]._flag == pf_BULK)
				node.get_neighbor_node(Direction::z_up)[phaseIndex]._flag = pf_NEAR_INTERFACE;
		
	}
	inline void FieldStorage_forPhaseNode::downgrade(PhaseNode& node, int phaseIndex) {
		if (node[phaseIndex]._flag == pf_INTERFACE) {
			node[phaseIndex]._flag = pf_NEAR_INTERFACE;
			return;
		}
		bool _xu = false, _xd = false, _yu = false, _yd = false, _zu = false, _zd = false;
		if ((node._x + 1) < limit_x || _BC.x_bc == PERIODIC)
			_xu = (*this).is_interface(node.get_neighbor_node(Direction::x_up), phaseIndex);
		if ((node._x - 1) >= 0 || _BC.x_bc == PERIODIC)
			_xd = (*this).is_interface(node.get_neighbor_node(Direction::x_down), phaseIndex);
		if ((node._y + 1) < limit_y || _BC.y_bc == PERIODIC)
			_yu = (*this).is_interface(node.get_neighbor_node(Direction::y_up), phaseIndex);
		if ((node._y - 1) >= 0 || _BC.y_bc == PERIODIC)
			_yd = (*this).is_interface(node.get_neighbor_node(Direction::y_down), phaseIndex);
		if ((node._z + 1) < limit_z || _BC.z_bc == PERIODIC)
			_zu = (*this).is_interface(node.get_neighbor_node(Direction::z_up), phaseIndex);
		if ((node._z - 1) >= 0 || _BC.z_bc == PERIODIC)
			_zd = (*this).is_interface(node.get_neighbor_node(Direction::z_down), phaseIndex);
		if (_xu || _xd || _yu || _yd || _zu || _zd) {
			node[phaseIndex]._flag = pf_NEAR_INTERFACE;
		}
		else {
			node[phaseIndex]._flag = pf_BULK;
			if (node[phaseIndex].phaseFraction < Simulation_Num_Cut_Off)
				node[phaseIndex].phaseFraction = 0.0;
			if (node[phaseIndex].phaseFraction > (1 - Simulation_Num_Cut_Off))
				node[phaseIndex].phaseFraction = 1.0;
		}
	}
	inline bool FieldStorage_forPhaseNode::is_interface(PhaseNode& node, int phaseIndex) {
		if (node[phaseIndex].phaseFraction >= Simulation_Num_Cut_Off && node[phaseIndex].phaseFraction <= (1.0 - Simulation_Num_Cut_Off))
			return true;
		else
			return false;
	}
	inline PhaseNode& FieldStorage_forPhaseNode::operator()(const int x, const int y, const int z) {
#ifdef _DEBUG
		if (limit_x <= 0 || limit_y <= 0 || limit_z <= 0) {
			cout << "Settings error ! limit value error !" << endl;
			SYS_PROGRAM_STOP;
		}
#endif
		int _x = x, _y = y, _z = z;
		if (x >= 0 && x < limit_x && y >= 0 && y < limit_y && z >= 0 && z < limit_z) {
			return _mesh[x + y * limit_x + z * limit_x * limit_y];
		}
		if (x < 0) {
			if (_BC.x_bc == PERIODIC)
				_x = _x + limit_x;
			else if (_BC.x_bc == FIXED)
				return _BC(Boundary::DOWN_X, y, z);
			else if (_BC.x_bc == ADIABATIC)
				_x = 0;
		}
		if (x >= limit_x) {
			if (_BC.x_bc == PERIODIC)
				_x = _x - limit_x;
			else if (_BC.x_bc == FIXED)
				return _BC(Boundary::UP_X, y, z);
			else if (_BC.x_bc == ADIABATIC)
				_x = limit_x - 1;
		}
		if (y < 0) {
			if (_BC.y_bc == PERIODIC)
				_y = _y + limit_y;
			else if (_BC.y_bc == FIXED)
				return _BC(Boundary::DOWN_Y, x, z);
			else if (_BC.y_bc == ADIABATIC)
				_y = 0;
		}
		if (y >= limit_y) {
			if (_BC.y_bc == PERIODIC)
				_y = _y - limit_y;
			else if (_BC.y_bc == FIXED)
				return _BC(Boundary::UP_Y, x, z);
			else if (_BC.y_bc == ADIABATIC)
				_y = limit_y - 1;
		}
		if (z < 0) {
			if (_BC.z_bc == PERIODIC)
				_z = _z + limit_z;
			else if (_BC.z_bc == FIXED)
				return _BC(Boundary::DOWN_Z, x, y);
			else if (_BC.z_bc == ADIABATIC)
				_z = 0;
		}
		if (z >= limit_z) {
			if (_BC.z_bc == PERIODIC)
				_z = _z - limit_z;
			else if (_BC.z_bc == FIXED)
				return _BC(Boundary::UP_Z, x, y);
			else if (_BC.z_bc == ADIABATIC)
				_z = limit_z - 1;
		}
		return (*this)(_x, _y, _z);
	}
	inline FieldStorage_forPhaseNode& FieldStorage_forPhaseNode::operator=(const FieldStorage_forPhaseNode& n) {
		_mesh = n._mesh;
		limit_x = n.limit_x;
		limit_y = n.limit_y;
		limit_z = n.limit_z;
		_BC = n._BC;
		dx = n.dx;
		_dimention = n._dimention;
		return *this;
	}

	class FieldStorage_forMechanicNode {
	public:
		FieldStorage_forMechanicNode() {
		};
		~FieldStorage_forMechanicNode() {
			free();
		};

		void init(int _Nx, int _Ny, int _Nz, double _dx, BoundaryCondition bc_x, BoundaryCondition bc_y, BoundaryCondition bc_z) {
			//<maxNum 
			Nx = _Nx;
			Ny = _Ny;
			Nz = _Nz;
			if (bc_x != PERIODIC) {
				Nx = _Nx * 2;
			}
			if (bc_y != PERIODIC) {
				Ny = _Ny * 2;
			}
			if (bc_z != PERIODIC) {
				Nz = _Nz * 2;
			}
			_mesh.resize(Nx * Ny * Nz);
			_virtual_strain_buff.resize(Nx * Ny * Nz);
			dx = _dx;
		}
		void free() {
			_mesh.clear();
			_virtual_strain_buff.clear();
			Nx = 0;
			Ny = 0;
			Nz = 0;
			dx = 0.0;
		}
		int size() {
			return int(_mesh.size());
		}

		MechanicNode& operator()(const int x, const int y, const int z) {
			if (x >= 0 && x < Nx && y >= 0 && y < Ny && z >= 0 && z < Nz)
				return _mesh[x + y * Nx + z * Nx * Ny];
			else {
				cout << "FieldStorage_forMechanicNode error ! limit value error !" << endl;
				SYS_PROGRAM_STOP;
			}
		}

		void save_virtual_strain_to_buff() {
			for (int index = 0; index < Nx * Ny * Nz; index++)
				_virtual_strain_buff[index] = _mesh[index].VirtualEigenStrains;
		}

		void assign_virtual_strain_with_buff() {
			for (int index = 0; index < Nx * Ny * Nz; index++)
				_mesh[index].VirtualEigenStrains = _virtual_strain_buff[index];
		}

		std::vector<MechanicNode> _mesh;
		std::vector<vStrain> _virtual_strain_buff;
		int Nx;                                                                ///< System size along X direction
		int Ny;                                                                ///< System size along y direction
		int Nz;                                                                ///< System size along z direction
		double dx;
	};

	struct VectorNode {
		VectorNode() {
		}
		~VectorNode() {
			clear();
		}
		void clear() {
			_up_x = nullptr;
			_down_x = nullptr;
			_up_y = nullptr;
			_down_y = nullptr;
			_up_z = nullptr;
			_down_z = nullptr;
			vals.clear();
		}
		void connect_mesh(int x, int y, int z, VectorNode& up_x, VectorNode& down_x
			, VectorNode& up_y, VectorNode& down_y, VectorNode& up_z, VectorNode& down_z) {
			_x = x;
			_y = y;
			_z = z;
			_up_x = &up_x;
			_down_x = &down_x;
			_up_y = &up_y;
			_down_y = &down_y;
			_up_z = &up_z;
			_down_z = &down_z;
		}
		VectorNode& get_neighbor_node(Direction _d);
		VectorNode& get_long_range_node(int relative_x, int relative_y, int relative_z);
		vector<double> vals;
		int _x;
		int _y;
		int _z;
	private:
		VectorNode* _up_x;
		VectorNode* _down_x;
		VectorNode* _up_y;
		VectorNode* _down_y;
		VectorNode* _up_z;
		VectorNode* _down_z;
	};
	inline VectorNode& VectorNode::get_neighbor_node(Direction _d) {
		switch (_d)
		{
		case pf::x_up:
			return *_up_x;
			break;
		case pf::x_down:
			return *_down_x;
			break;
		case pf::y_up:
			return *_up_y;
			break;
		case pf::y_down:
			return *_down_y;
			break;
		case pf::z_up:
			return *_up_z;
			break;
		case pf::z_down:
			return *_down_z;
			break;
		default:
			cout << "Class get_neighbor_node parameter error !" << endl;
			SYS_PROGRAM_STOP;
			break;
		}
	}
	inline VectorNode& VectorNode::get_long_range_node(int relative_x, int relative_y, int relative_z) {
		if (relative_x > 0)
			return (*this).get_neighbor_node(Direction::x_up).get_long_range_node(relative_x - 1, relative_y, relative_z);
		else if (relative_x < 0)
			return (*this).get_neighbor_node(Direction::x_down).get_long_range_node(relative_x + 1, relative_y, relative_z);
		if (relative_y > 0)
			return (*this).get_neighbor_node(Direction::y_up).get_long_range_node(relative_x, relative_y - 1, relative_z);
		else if (relative_y < 0)
			return (*this).get_neighbor_node(Direction::y_down).get_long_range_node(relative_x, relative_y + 1, relative_z);
		if (relative_z > 0)
			return (*this).get_neighbor_node(Direction::z_up).get_long_range_node(relative_x, relative_y, relative_z - 1);
		else if (relative_z < 0)
			return (*this).get_neighbor_node(Direction::z_down).get_long_range_node(relative_x, relative_y, relative_z + 1);
		return (*this);
	}
	class FieldStorage_forVector
	{
	public:
		FieldStorage_forVector() {};
		~FieldStorage_forVector() {
			try {
				_mesh.clear();
			}
			catch (...) {
				std::cout << "Destroy FieldStorage_forVector failed, problems occur." << endl;
			}
		};

		VectorNode& operator()(const int x, const int y, const int z);                                       //< Index operator for accessing the n's field value
		FieldStorage_forVector& operator=(const FieldStorage_forVector& n);
		void init(int Nx, int Ny, int Nz, double _dx, BoundaryCondition bc_x, BoundaryCondition bc_y, BoundaryCondition bc_z) {
			//<maxNum 
			_mesh.resize(Nx * Ny * Nz);
			limit_x = Nx;
			limit_y = Ny;
			limit_z = Nz;
			x_bc = bc_x;
			y_bc = bc_y;
			z_bc = bc_z;
			dx = _dx;
			_dimention = Three_Dimension;
			if (Nx == 1 || Ny == 1 || Nz == 1)
				_dimention = Two_Dimension;
			if ((Nx == 1 && Ny == 1) || (Nx == 1 && Nz == 1) || (Ny == 1 && Nz == 1))
				_dimention = One_Dimension;
			for (int z = 0; z < Nz; z++)
				for (int y = 0; y < Ny; y++)
					for (int x = 0; x < Nx; x++) {
						VectorNode* _up_x;
						VectorNode* _down_x;
						VectorNode* _up_y;
						VectorNode* _down_y;
						VectorNode* _up_z;
						VectorNode* _down_z;
						if (bc_x == BoundaryCondition::ADIABATIC || bc_x == BoundaryCondition::FIXED) {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
							_down_x = &_mesh[down];
						}
						else {
							int up = (x + 1) + y * limit_x + z * limit_x * limit_y;
							int down = (x - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == 0)
								down = (Nx - 1) + y * limit_x + z * limit_x * limit_y;
							if (x == (Nx - 1))
								up = 0 + y * limit_x + z * limit_x * limit_y;
							_up_x = &_mesh[up];
							_down_x = &_mesh[down];
						}
						if (bc_y == BoundaryCondition::ADIABATIC || bc_y == BoundaryCondition::FIXED) {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
							_down_y = &_mesh[down];
						}
						else {
							int up = x + (y + 1) * limit_x + z * limit_x * limit_y;
							int down = x + (y - 1) * limit_x + z * limit_x * limit_y;
							if (y == 0)
								down = x + (Ny - 1) * limit_x + z * limit_x * limit_y;
							if (y == (Ny - 1))
								up = x + 0 * limit_x + z * limit_x * limit_y;
							_up_y = &_mesh[up];
							_down_y = &_mesh[down];
						}
						if (bc_z == BoundaryCondition::ADIABATIC || bc_z == BoundaryCondition::FIXED) {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + z * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + z * limit_x * limit_y;
							_up_z = &_mesh[up];
							_down_z = &_mesh[down];
						}
						else {
							int up = x + y * limit_x + (z + 1) * limit_x * limit_y;
							int down = x + y * limit_x + (z - 1) * limit_x * limit_y;
							if (z == 0)
								down = x + y * limit_x + (Nz - 1) * limit_x * limit_y;
							if (z == (Nz - 1))
								up = x + y * limit_x + 0 * limit_x * limit_y;
							_up_z = &_mesh[up];
							_down_z = &_mesh[down];
						}
						_mesh[x + y * limit_x + z * limit_x * limit_y].connect_mesh(x, y, z, *_up_x, *_down_x, *_up_y, *_down_y, *_up_z, *_down_z);
					}
		}
		void free() {
			for (auto node = _mesh.begin(); node < _mesh.end(); node++)
				node->clear();
			_mesh.clear();
			limit_x = 0;
			limit_y = 0;
			limit_z = 0;
			dx = 0.0;
		}
		int size() {
			return int(_mesh.size());
		}

		// ******   Settings   ****** //
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		int limit_x;
		int limit_y;
		int limit_z;
		double dx;
		Dimension _dimention;
		std::vector<VectorNode> _mesh;
	};
	inline VectorNode& FieldStorage_forVector::operator()(const int x, const int y, const int z) {
#ifdef _DEBUG
		if (limit_x <= 0 || limit_y <= 0 || limit_z <= 0) {
			cout << "Settings error ! mesh size error !" << endl;
			SYS_PROGRAM_STOP;
		}
#endif
		int _x = x, _y = y, _z = z;
		if (x >= 0 && x < limit_x && y >= 0 && y < limit_y && z >= 0 && z < limit_z) {
			return _mesh[x + y * limit_x + z * limit_x * limit_y];
		}
		if (x < 0) {
			if (x_bc == PERIODIC)
				_x = _x + limit_x;
			else
				_x = 0;
		}
		if (x >= limit_x) {
			if (x_bc == PERIODIC)
				_x = _x - limit_x;
			else
				_x = limit_x - 1;
		}
		if (y < 0) {
			if (y_bc == PERIODIC)
				_y = _y + limit_y;
			else
				_y = 0;
		}
		if (y >= limit_y) {
			if (y_bc == PERIODIC)
				_y = _y - limit_y;
			else
				_y = limit_y - 1;
		}
		if (z < 0) {
			if (z_bc == PERIODIC)
				_z = _z + limit_z;
			else
				_z = 0;
		}
		if (z >= limit_z) {
			if (z_bc == PERIODIC)
				_z = _z - limit_z;
			else
				_z = limit_z - 1;
		}
		return (*this)(_x, _y, _z);
	}
	inline FieldStorage_forVector& FieldStorage_forVector::operator=(const FieldStorage_forVector& n) {
		_mesh = n._mesh;
		dx = n.dx;
		_dimention = n._dimention;
		limit_x = n.limit_x;
		limit_y = n.limit_y;
		limit_z = n.limit_z;
		x_bc = n.x_bc;
		y_bc = n.y_bc;
		z_bc = n.z_bc;
		return *this;
	}
}
