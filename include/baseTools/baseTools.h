#pragma once
#include "sysTool.h"
#include "NucleationTools.h"
#include "PhaseNode.h"
#include "FieldStorage.h"
#include "timer.h"
#include "WriteToFile.h"
#include "BoundaryCondition.h"
#include "OuterFunction.h"
#include "DataFile.h"
#include "InfoTools.h"
#include "RotationMatrix.h"
#ifdef _WIN32
#include "BMP24Reader.h"
#endif
namespace pf {
	using namespace std;

	static vector<double> getCommonTangent(vector<double> G1, vector<double> G2, vector<double> x1, vector<double> x2)
	{ //return two points of commonTangent
		vector<double> val;
		for (unsigned int i = 1; i < x1.size() - 1; i++) {
			double slopei1 = (G1[i - 1] - G1[i]) / (x1[i - 1] - x1[i]);
			double slopei2 = (G1[i] - G1[i + 1]) / (x1[i] - x1[i + 1]);
			for (unsigned int j = 1; j < x2.size() - 1; j++) {
				double slopej1 = (G2[j - 1] - G2[j]) / (x2[j - 1] - x2[j]);
				double slopej2 = (G2[j] - G2[j + 1]) / (x2[j] - x2[j + 1]);
				double slope = (G1[i] - G2[j]) / (x1[i] - x2[j]);
				if ((slope > slopei1 && slope < slopei2) || (slope < slopei1 && slope > slopei2)) {
					if ((slope > slopej1 && slope < slopej2) || (slope < slopej1 && slope > slopej2)) {
						val = { x1[i] , x2[j] };
						return val;
					}
				}
			}
		}
		cout << "Don't have Tangent!!!";
		getchar();
		return val;
	}

	static double getTangent(vector<double> point, vector<double> G2, vector<double> x2)
	{ //return one point of Tangent
		for (unsigned int i = 1; i < x2.size() - 1; i++) {
			double slope0 = (point[1] - G2[i]) / (point[0] - x2[i]);
			double slope1 = (G2[i - 1] - G2[i]) / (x2[i - 1] - x2[i]);
			double slope2 = (G2[i] - G2[i + 1]) / (x2[i] - x2[i + 1]);
			if ((slope1 > slope0 && slope2 < slope0) || (slope2 > slope0 && slope1 < slope0)) {
				return x2[i];
			}
		}
		cout << "Don't have Tangent!!!";
		getchar();
		return 0.0;
	}


	class Criterions {
	public:
		void push_back(int& n) {
			indexs.push_back(n);
		}
		vector<int> indexs;
		int step;
		Criterions&  operator=(const Criterions& n);
	};
	inline Criterions& Criterions::operator=(const Criterions& n) {
		indexs = n.indexs;
		step = n.step;
		return *this;
	}

	//二维double矢量
	struct  Vec2d
	{
		double x, y;

		Vec2d()
		{
			x = 0.0;
			y = 0.0;
		}
		Vec2d(double dx, double dy)
		{
			x = dx;
			y = dy;
		}
		void Set(double dx, double dy)
		{
			x = dx;
			y = dy;
		}
	};

	//判断点在线段上
	static bool IsPointOnLine(double px0, double py0, double px1, double py1, double px2, double py2)
	{
		bool flag = false;
		double d1 = (px1 - px0) * (py2 - py0) - (px2 - px0) * (py1 - py0);
		if ((fabs(d1) < SYS_EPSILON) && ((px0 - px1) * (px0 - px2) <= 0) && ((py0 - py1) * (py0 - py2) <= 0))
		{
			flag = true;
		}
		return flag;
	}

	//判断两线段相交
	static bool IsIntersect(double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4)
	{
		bool flag = false;
		double d = (px2 - px1) * (py4 - py3) - (py2 - py1) * (px4 - px3);
		if (d != 0)
		{
			double r = ((py1 - py3) * (px4 - px3) - (px1 - px3) * (py4 - py3)) / d;
			double s = ((py1 - py3) * (px2 - px1) - (px1 - px3) * (py2 - py1)) / d;
			if ((r >= 0) && (r <= 1) && (s >= 0) && (s <= 1))
			{
				flag = true;
			}
		}
		return flag;
	}

	//判断点在多边形内
	static bool is_Point_In_Polygon_2D(double x, double y, const vector<Vec2d>& POL)
	{
		bool isInside = false;
		int count = 0;

		//
		double minX = DBL_MAX;
		for (int i = 0; i < int(POL.size()); i++)
		{
			minX = min(minX, POL[i].x);
		}

		//
		double px = x;
		double py = y;
		double linePoint1x = x;
		double linePoint1y = y;
		double linePoint2x = minX - 10;			//取最小的X值还小的值作为射线的终点
		double linePoint2y = y;

		//遍历每一条边
		for (int i = 0; i < int(POL.size()) - 1; i++)
		{
			double cx1 = POL[i].x;
			double cy1 = POL[i].y;
			double cx2 = POL[i + 1].x;
			double cy2 = POL[i + 1].y;

			if (IsPointOnLine(px, py, cx1, cy1, cx2, cy2))
			{
				return true;
			}

			if (fabs(cy2 - cy1) < SYS_EPSILON)   //平行则不相交
			{
				continue;
			}

			if (IsPointOnLine(cx1, cy1, linePoint1x, linePoint1y, linePoint2x, linePoint2y))
			{
				if (cy1 > cy2)			//只保证上端点+1
				{
					count++;
				}
			}
			else if (IsPointOnLine(cx2, cy2, linePoint1x, linePoint1y, linePoint2x, linePoint2y))
			{
				if (cy2 > cy1)			//只保证上端点+1
				{
					count++;
				}
			}
			else if (IsIntersect(cx1, cy1, cx2, cy2, linePoint1x, linePoint1y, linePoint2x, linePoint2y))   //已经排除平行的情况
			{
				count++;
			}
		}

		if (count % 2 == 1)
		{
			isInside = true;
		}

		return isInside;
	}


	static double gauss_func_1D(double peak_distance_origin, double peak_high, double variance, double point_distance_origin) {
		double gauss_value = peak_high * exp(-(point_distance_origin - peak_distance_origin) * (point_distance_origin - peak_distance_origin)
			/ 2 / variance / variance);
		return gauss_value;
	}

	static void debug_cout(string valueName, double value) {
		cout << endl << "Debug the value: " << valueName << " = " << value << endl;
	}

}