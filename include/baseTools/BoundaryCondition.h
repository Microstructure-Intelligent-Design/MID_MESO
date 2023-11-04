#pragma once
#include "sysTool.h"
#include "NormalNode.h"
#include "SublatticNode.h"
#include "PairValueNode.h"
#include "PhaseNode.h"
using namespace std;
/*  
	 0		x           limitX
	       DOWNY
0  D *******************  U
   O *******************  P
y  W *******************  X
   N *******************
   X *******************
	 *******************
	 *******************
limitY     UPY
*/

namespace pf {
	static int Index_fix_bc_axisX(int y, int z, int limit_y) {
		return y + z * limit_y;
	}
	static int Index_fix_bc_axisY(int x, int z, int limit_x) {
		return x + z * limit_x;
	}
	static int Index_fix_bc_axisZ(int x, int y, int limit_x) {
		return x + y * limit_x;
	}
	class BoundaryCondition_forPhaseNode
	{
	public:
		BoundaryCondition_forPhaseNode() {};
		~BoundaryCondition_forPhaseNode() {};
		void init(BoundaryCondition bc_x, BoundaryCondition bc_y, BoundaryCondition bc_z, int x, int y, int z) {
			x_bc = bc_x;
			y_bc = bc_y;
			z_bc = bc_z;
			limit_x = x;
			limit_y = y;
			limit_z = z;
			if (bc_x == FIXED) {
				up_x.resize(z * y);
				down_x.resize(z * y);
				for(int k = 0; k < z; k++)
					for (int j = 0; j < y; j++) {
						int up_j = j + 1, down_j = j - 1, up_k = k + 1, down_k = k - 1;
						if (down_j < 0)
							down_j = 0;
						if (up_j > y - 1)
							up_j = y - 1;
						if (down_k < 0)
							down_k = 0;
						if (up_k > z - 1)
							up_k = z - 1;
						{
							PhaseNode& _up_x =   up_x[Index_fix_bc_axisX(j, k, limit_y)];
							PhaseNode& _down_x = up_x[Index_fix_bc_axisX(j, k, limit_y)];
							PhaseNode& _up_y =   up_x[Index_fix_bc_axisX(up_j, k, limit_y)];
							PhaseNode& _down_y = up_x[Index_fix_bc_axisX(down_j, k, limit_y)];
							PhaseNode& _up_z =   up_x[Index_fix_bc_axisX(j, up_k, limit_y)];
							PhaseNode& _down_z = up_x[Index_fix_bc_axisX(j, down_k, limit_y)];
							up_x[Index_fix_bc_axisX(j, k, limit_y)].connect_mesh(x, j, k, _up_x, _down_x, _up_y, _down_y, _up_z, _down_z);
						}
						{
							PhaseNode& _up_x =   down_x[Index_fix_bc_axisX(j, k, limit_y)];
							PhaseNode& _down_x = down_x[Index_fix_bc_axisX(j, k, limit_y)];
							PhaseNode& _up_y =   down_x[Index_fix_bc_axisX(up_j, k, limit_y)];
							PhaseNode& _down_y = down_x[Index_fix_bc_axisX(down_j, k, limit_y)];
							PhaseNode& _up_z =   down_x[Index_fix_bc_axisX(j, up_k, limit_y)];
							PhaseNode& _down_z = down_x[Index_fix_bc_axisX(j, down_k, limit_y)];
							down_x[Index_fix_bc_axisX(j, k, limit_y)].connect_mesh(-1, j, k, _up_x, _down_x, _up_y, _down_y, _up_z, _down_z);
						}
					}
			}
			if (bc_y == FIXED) {
				up_y.resize(x * z);
				down_y.resize(x * z);
				for (int k = 0; k < z; k++)
					for (int i = 0; i < x; i++) {
						int up_i = i + 1, down_i = i - 1, up_k = k + 1, down_k = k - 1;
						if (down_i < 0)
							down_i = 0;
						if (up_i > x - 1)
							up_i = x - 1;
						if (down_k < 0)
							down_k = 0;
						if (up_k > z - 1)
							up_k = z - 1;
						{
							PhaseNode& _up_x =   up_y[Index_fix_bc_axisY(up_i, k, limit_x)];
							PhaseNode& _down_x = up_y[Index_fix_bc_axisY(down_i, k, limit_x)];
							PhaseNode& _up_y =   up_y[Index_fix_bc_axisY(i, k, limit_x)];
							PhaseNode& _down_y = up_y[Index_fix_bc_axisY(i, k, limit_x)];
							PhaseNode& _up_z =   up_y[Index_fix_bc_axisY(i, up_k, limit_x)];
							PhaseNode& _down_z = up_y[Index_fix_bc_axisY(i, down_k, limit_x)];
							up_y[Index_fix_bc_axisY(i, k, limit_x)].connect_mesh(i, y, k, _up_x, _down_x, _up_y, _down_y, _up_z, _down_z);
						}
						{
							PhaseNode& _up_x =   down_y[Index_fix_bc_axisY(up_i, k, limit_x)];
							PhaseNode& _down_x = down_y[Index_fix_bc_axisY(down_i, k, limit_x)];
							PhaseNode& _up_y =   down_y[Index_fix_bc_axisY(i, k, limit_x)];
							PhaseNode& _down_y = down_y[Index_fix_bc_axisY(i, k, limit_x)];
							PhaseNode& _up_z =   down_y[Index_fix_bc_axisY(i, up_k, limit_x)];
							PhaseNode& _down_z = down_y[Index_fix_bc_axisY(i, down_k, limit_x)];
							down_y[Index_fix_bc_axisY(i, k, limit_x)].connect_mesh(i, -1, k, _up_x, _down_x, _up_y, _down_y, _up_z, _down_z);
						}
					}
			}
			if (bc_z == FIXED) {
				down_z.resize(x * y);
				up_z.resize(x * y);
				for (int j = 0; j < y; j++)
					for (int i = 0; i < x; i++) {
						int up_i = i + 1, down_i = i - 1, up_j = j + 1, down_j = j - 1;
						if (down_i < 0)
							down_i = 0;
						if (up_i > x - 1)
							up_i = x - 1;
						if (down_j < 0)
							down_j = 0;
						if (up_j > y - 1)
							up_j = y - 1;
						{
							PhaseNode& _up_x =   up_z[Index_fix_bc_axisZ(up_i, j, limit_x)];
							PhaseNode& _down_x = up_z[Index_fix_bc_axisZ(down_i, j, limit_x)];
							PhaseNode& _up_y =   up_z[Index_fix_bc_axisZ(i, up_j, limit_x)];
							PhaseNode& _down_y = up_z[Index_fix_bc_axisZ(i, down_j, limit_x)];
							PhaseNode& _up_z =   up_z[Index_fix_bc_axisZ(i, j, limit_x)];
							PhaseNode& _down_z = up_z[Index_fix_bc_axisZ(i, j, limit_x)];
							up_z[Index_fix_bc_axisZ(i, j, limit_x)].connect_mesh(i, j, z, _up_x, _down_x, _up_y, _down_y, _up_z, _down_z);
						}
						{
							PhaseNode& _up_x =   down_z[Index_fix_bc_axisZ(up_i, j, limit_x)];
							PhaseNode& _down_x = down_z[Index_fix_bc_axisZ(down_i, j, limit_x)];
							PhaseNode& _up_y =   down_z[Index_fix_bc_axisZ(i, up_j, limit_x)];
							PhaseNode& _down_y = down_z[Index_fix_bc_axisZ(i, down_j, limit_x)];
							PhaseNode& _up_z =   down_z[Index_fix_bc_axisZ(i, j, limit_x)];
							PhaseNode& _down_z = down_z[Index_fix_bc_axisZ(i, j, limit_x)];
							down_z[Index_fix_bc_axisZ(i, j, limit_x)].connect_mesh(i, j, -1, _up_x, _down_x, _up_y, _down_y, _up_z, _down_z);
						}
					}
			}
		};
		PhaseNode& operator()(Boundary b, int i, int j = 0);                                       //< Index operator for accessing the n's field value
		void free() {
			for (auto node = up_x.begin(); node < up_x.end(); node++)
				node->clear();
			for (auto node = down_x.begin(); node < down_x.end(); node++)
				node->clear();
			for (auto node = up_y.begin(); node < up_y.end(); node++)
				node->clear();
			for (auto node = down_y.begin(); node < down_y.end(); node++)
				node->clear();
			for (auto node = up_z.begin(); node < up_z.end(); node++)
				node->clear();
			for (auto node = down_z.begin(); node < down_z.end(); node++)
				node->clear();
			up_x.clear();
			down_x.clear();
			up_y.clear();
			down_y.clear();
			up_z.clear();
			down_z.clear();
		};

		size_t size(BoundaryCondition bc) {
			if (bc == UP_X) return up_x.size();
			else if (bc == UP_Y) return up_y.size();
			else if (bc == UP_Z) return up_z.size();
			else if (bc == DOWN_X) return down_x.size();
			else if (bc == DOWN_Y) return down_y.size();
			else if (bc == DOWN_Z) return down_z.size();
		}
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		int limit_x;
		int limit_y;
		int limit_z;

		std::vector<PhaseNode> up_x;
		std::vector<PhaseNode> up_y;
		std::vector<PhaseNode> up_z;
		std::vector<PhaseNode> down_x;
		std::vector<PhaseNode> down_y;
		std::vector<PhaseNode> down_z;
	};
	inline PhaseNode& BoundaryCondition_forPhaseNode::operator()(Boundary b, int i, int j) {
		if (b == UP_X) {
			if (i < 0)
				i = 0;
			else if (i >= limit_y)
				i = limit_y - 1;
			if (j < 0)
				j = 0;
			else if (j >= limit_z)
				j = limit_z - 1;

			return up_x[i + j * limit_y];
		}
		else if (b == UP_Y) {
			if (i < 0)
				i = 0;
			else if (i >= limit_x)
				i = limit_x - 1;
			if (j < 0)
				j = 0;
			else if (j >= limit_z)
				j = limit_z - 1;
			
			return up_y[i + j * limit_x]; 
		}
		else if (b == UP_Z) {
			if (i < 0)
				i = 0;
			else if (i >= limit_x)
				i = limit_x - 1;
			if (j < 0)
				j = 0;
			else if (j >= limit_y)
				j = limit_y - 1;

			return up_z[i + j * limit_x]; 
		}
		else if (b == DOWN_X) {
			if (i < 0)
				i = 0;
			else if (i >= limit_y)
				i = limit_y - 1;
			if (j < 0)
				j = 0;
			else if (j >= limit_z)
				j = limit_z - 1;

			return down_x[i + j * limit_y];
		}
		else if (b == DOWN_Y) {
			if (i < 0)
				i = 0;
			else if (i >= limit_x)
				i = limit_x - 1;
			if (j < 0)
				j = 0;
			else if (j >= limit_z)
				j = limit_z - 1;

			return down_y[i + j * limit_x]; 
		}
		else if (b == DOWN_Z) {
			if (i < 0)
				i = 0;
			else if (i >= limit_x)
				i = limit_x - 1;
			if (j < 0)
				j = 0;
			else if (j >= limit_y)
				j = limit_y - 1;

			return down_z[i + j * limit_x]; 
		}
		cout << "Boundarycondition error, don't have this boundary!" << endl;
		SYS_PROGRAM_STOP;
	}

}
