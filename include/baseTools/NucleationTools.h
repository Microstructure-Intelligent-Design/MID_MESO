#pragma once
#include "sysTool.h"
#include "NormalNode.h"
#include "BoundaryCondition.h"
#include "RotationMatrix.h"

#define NP_AnyPhase int(-1)

#define NP_x_DownLimit value
#define NP_x_UpLimit value2

namespace pf {
	using namespace std;

	enum NucleationProperty { DefiniteNucleation, ConditionalNucleation, UserDefined_Nucleated};
	enum Geometry{ Geo_None, Geo_Ellipsoid, Geo_Polyhedron};
	enum NucleationPosition { NP_Anywhere, NP_Bulk, NP_Interface };
	enum NucleationConcentration { Default_X, Defined_X};

	struct line_func_2D {
		// temp: a * x + b * y + c = 0
		line_func_2D() {
			a = 0.0;
			b = 0.0;
			c = 0.0;
		}
		line_func_2D(Point2D p1, Point2D p2) {
			init(p1, p2);
		}
		line_func_2D(double x1, double y1, double x2, double y2) {
			init(x1, y1, x2, y2);
		}
		void init(Point2D p1, Point2D p2) {
			a = p2.y - p1.y;
			b = p1.x - p2.x;
			c = p2.x * p1.y - p2.y * p1.x;
		}
		void init(double x1, double y1, double x2, double y2) {
			a = y2 - y1;
			b = x1 - x2;
			c = x2 * y1 - y2 * x1;
		}
		bool is_point_on_line(double x, double y) {
			if (a == 0 && b == 0)
				return false;
			double val = a * x + b * y + c;
			if (Is_Equality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_above_line(double x, double y) {
			if (a == 0 && b == 0)
				return false;
			if ((a * x + b * y + c) > 0.0)
				return true;
			else
				return false;
		}
		double a, b, c;
	};

	struct surf_func_3D {
		// temp: a * x + b * y + c * z + d = 0
		surf_func_3D() {
			a = 0.0;
			b = 0.0;
			c = 0.0;
			d = 0.0;
		}
		surf_func_3D(Point p1, Point p2, Point p3) {
			init(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
		}
		surf_func_3D(Point normal, Point point) {
			a = normal.x;
			b = normal.y;
			c = normal.z;
			d = -(a * point.x + b * point.y + c * point.z);
		}
		void init(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
			Vector3 t1(x2 - x1, y2 - y1, z2 - z1);
			Vector3 t2(x3 - x1, y3 - y1, z3 - z1);
			t1 = t1.cross(t2);
			a = t1[0];
			b = t1[1];
			c = t1[2];
			d = -(a * x1 + b * y1 + c * z1);
		}
		bool is_point_on_surf(double x, double y, double z) {
			double val = a * x + b * y + c * z + d;
			if (Is_Equality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_on_surf(Point p) {
			double val = a * p.x + b * p.y + c * p.z + d;
			if (Is_Equality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_above_surf(double x, double y, double z) {
			if ((a * x + b * y + c * z + d) > 0.0)
				return true;
			else
				return false;
		}
		bool is_point_above_surf(Point p) {
			if ((a * p.x + b * p.y + c * p.z + d) > 0.0)
				return true;
			else
				return false;
		}
		surf_func_3D& operator=(const surf_func_3D& n) {
			a = n.a;
			b = n.b;
			c = n.c;
			d = n.d;
			return *this;
		}
		double a, b, c, d;
	};
	struct Ellipsoid {
		//> (x - core.x) * (x - core.x) / radius_x / radius_x + (y - core.y) * (y - core.y) / radius_y / radius_y + (z - core.z) * (z - core.z) / radius_z / radius_z = 1.0
		Ellipsoid(){
			radius_x = 0.0;
			radius_y = 0.0;
			radius_z = 0.0;
			radian[0] = 0.0;
			radian[1] = 0.0;
			radian[2] = 0.0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void set_core(double x, double y, double z) {
			core.x = x;
			core.y = y;
			core.z = z;
		}
		void set_radius(double x_radius, double y_radius, double z_radius) {
			radius_x = abs(x_radius) + SYS_EPSILON;
			radius_y = abs(y_radius) + SYS_EPSILON;
			radius_z = abs(z_radius) + SYS_EPSILON;
		}
		void set_rotation_radian_and_rotation_gauge(double _radian[3], RotationGauge _rotationGauge) {
			radian[0] = -_radian[0];
			radian[1] = -_radian[1];
			radian[2] = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		bool check_point_inside_ellipsoid(double x, double y, double z) {
			double test = (x - core.x) * (x - core.x) / radius_x / radius_x + (y - core.y) * (y - core.y) / radius_y / radius_y + (z - core.z) * (z - core.z) / radius_z / radius_z;
			if (test <= 1.0)
				return true;
			else
				return false;
		}
		bool check_point_inside_ellipsoid(Point p) {
			double test = (p.x - core.x) * (p.x - core.x) / radius_x / radius_x + (p.y - core.y) * (p.y - core.y) / radius_y / radius_y + (p.z - core.z) * (p.z - core.z) / radius_z / radius_z;
			if (test <= 1.0)
				return true;
			else
				return false;
		}
		Ellipsoid& operator=(const Ellipsoid& n) {
			this->core = n.core;
			this->radius_x = n.radius_x;
			this->radius_y = n.radius_y;
			this->radius_z = n.radius_z;
			return *this;
		}
		Point core;
		double radius_x;
		double radius_y;
		double radius_z;
		double radian[3];
		RotationGauge rotationGauge;
	};

	struct Polyhedron {
		Polyhedron(int inside_x = 0, int inside_y = 0, int inside_z = 0) {
			point_inside_polyhedron.set(inside_x, inside_y, inside_z);
			radian[0] = 0.0;
			radian[1] = 0.0;
			radian[2] = 0.0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		Polyhedron(Point inside_point) {
			point_inside_polyhedron.set(inside_point.x, inside_point.y, inside_point.z);
			radian[0] = 0.0;
			radian[1] = 0.0;
			radian[2] = 0.0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void set_rotation_radian_and_rotation_gauge(double _radian[3], RotationGauge _rotationGauge) {
			radian[0] = -_radian[0];
			radian[1] = -_radian[1];
			radian[2] = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		void add_surf(Point p1, Point p2, Point p3) {
			surf_func_3D surf(p1, p2, p3);
			surfaces.push_back(surf);
		}
		void add_surf(Point norm, Point p) {
			surf_func_3D surf(norm, p);
			surfaces.push_back(surf);
		}
		void set_a_point_inside_polyhedron(int x, int y, int z) {
			point_inside_polyhedron.set(x, y, z);
		}
		void set_a_point_inside_polyhedron(Point p) {
			point_inside_polyhedron.set(p.x, p.y, p.z);
		}
		bool check_point_inside_polyhedron(int x, int y, int z) {
			bool is_inside = true;
			for (auto surf = surfaces.begin(); surf < surfaces.end(); surf++)
				if (surf->is_point_above_surf(x, y, z) != surf->is_point_above_surf(point_inside_polyhedron) && !surf->is_point_on_surf(x, y, z))
					is_inside = false;
			return is_inside;
		}
		bool check_point_inside_polyhedron(Point p) {
			bool is_inside = true;
			for (auto surf = surfaces.begin(); surf < surfaces.end(); surf++)
				if (surf->is_point_above_surf(p) != surf->is_point_above_surf(point_inside_polyhedron) && !surf->is_point_on_surf(p))
					is_inside = false;
			return is_inside;
		}
		Polyhedron& operator=(const Polyhedron& n) {
			this->surfaces = n.surfaces;
			this->point_inside_polyhedron = n.point_inside_polyhedron;
			return *this;
		}
		vector<surf_func_3D> surfaces;
		Point point_inside_polyhedron;
		double radian[3];
		RotationGauge rotationGauge;
	};

	struct GeometricRegion {
		Geometry geometryProperty;
		Ellipsoid ellipSolid;
		Polyhedron polyhedron;
		int generate_step; 
		int phaseProperty;  
		int phaseIndex;   
		XNode x;
		double temperature;

		double_box customValues;
		int_box customFlags;
		vec3_box customVec3s;

		GeometricRegion(pf::Geometry _geometry_property = Geo_None, int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, double _temperature = 0.0) {
			init(_geometry_property, _generate_step, _phase_index, _phase_property, _temperature);
		}
		void init(pf::Geometry _geometry_property = Geo_None, int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, double _temperature = 0.0) {
			geometryProperty = _geometry_property;
			generate_step = _generate_step;
			phaseProperty = _phase_property;
			phaseIndex = _phase_index;
			temperature = _temperature;
		}
		GeometricRegion& operator=(const GeometricRegion& n) {
			this->generate_step = n.generate_step;
			this->geometryProperty = n.geometryProperty;
			this->temperature = n.temperature;
			this->phaseProperty = n.phaseProperty;
			this->phaseIndex = n.phaseIndex;
			this->ellipSolid = n.ellipSolid;
			this->polyhedron = n.polyhedron;
			this->x = n.x;
			this->customValues = n.customValues;
			this->customFlags = n.customFlags;
			this->customVec3s = n.customVec3s;
			return *this;
		}

	};

	class PointSet {
	public:
		PointSet(int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, double _temperature = 0.0) {
			init(_generate_step, _phase_index, _phase_property, _temperature);
		}
		void init(int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, double _temperature = 0.0) {
			generate_step = _generate_step;
			phaseProperty = _phase_property;
			phaseIndex = _phase_index;
			temperature = _temperature;
		}
		~PointSet() {
			points.clear();
		}
		vector<Point> points;
		int generate_step;
		int phaseProperty; 
		int phaseIndex;
		XNode x;
		double temperature;
		bool is_point_in_set(Point point, int x_limit, BoundaryCondition x_bc, int y_limit = 1, BoundaryCondition y_bc = BoundaryCondition::ADIABATIC,
			int z_limit = 1, BoundaryCondition z_bc = BoundaryCondition::ADIABATIC) {
			int x = double_to_int(point.x), y = double_to_int(point.y), z = double_to_int(point.z);
			if (x < 0) {
				if (x_bc == BoundaryCondition::PERIODIC)
					x += x_limit;
				else
					x = 0;
			}
			else if (x >= x_limit) {
				if (x_bc == BoundaryCondition::PERIODIC)
					x -= x_limit;
				else
					x = x_limit - 1;
			}
			if (y < 0) {
				if (y_bc == BoundaryCondition::PERIODIC)
					y += y_limit;
				else
					y = 0;
			}
			else if (y >= y_limit) {
				if (y_bc == BoundaryCondition::PERIODIC)
					y -= y_limit;
				else
					y = y_limit - 1;
			}
			if (z < 0) {
				if (z_bc == BoundaryCondition::PERIODIC)
					z += z_limit;
				else
					z = 0;
			}
			else if (z >= z_limit) {
				if (z_bc == BoundaryCondition::PERIODIC)
					z -= z_limit;
				else
					z = z_limit - 1;
			}

			for (auto p = points.begin(); p < points.end(); p++)
				if (p->x == x && p->y == y && p->z == z)
					return true;
			return false;
		}
		void add_point(int _x, int _y, int _z) {
			points.push_back(Point(_x, _y, _z));
		}
		PointSet& operator=(const PointSet& n) {
			this->generate_step = n.generate_step;
			this->points = n.points;
			this->temperature = n.temperature;
			this->phaseProperty = n.phaseProperty;
			this->phaseIndex = n.phaseIndex;
			this->x = n.x;
			return *this;
		}
	};

	struct Nucleus {
		int generate_step;    //range = [1, nstep)
		Point core;
		double phasefraction;   //range = (0, RADIUS_LIMIT)
		int phaseProperty;    //range = [0, phaseNum)
		int phaseIndex;     // 0 = basePhase, others[0, ~ ) = newPhase
		vector<int> contact_phases;   // phase property
		XNode x;
		Nucleus() {
			generate_step = 0;
			phasefraction = 1.0;
			phaseProperty = 0;
			phaseIndex = 0;
		}
		Nucleus& operator=(const Nucleus& n) {
			this->generate_step = n.generate_step;
			this->core = n.core;
			this->phasefraction = n.phasefraction;
			this->phaseProperty = n.phaseProperty;
			this->phaseIndex = n.phaseIndex;
			this->x = n.x;
			this->contact_phases = n.contact_phases;
			return *this;
		}
	};

	struct ConditionalPhase {
		int phaseProperty;    ///< property
		int phaseIndex;
		NucleationPosition nucleationPosition;
		vector<int> contactedPhases;   ///< property
		vector<int> inexistencePhases;
		XNode local_x_limit;  ///< local con in this range will generate wanted phase
		double global_phaseFraction_limit;   ///< aim phase nucleation will start below the global_phaseFraction 
		double possible_volume;
		double repel_radius;
		double nucleation_rate;   ///< 0.0 ~ 1.0
		NucleationConcentration nucleation_con;
		XNode defined_x;
		ConditionalPhase() {
			phaseProperty = 0;
			phaseIndex = 0;
			nucleationPosition = NP_Anywhere;
			repel_radius = 0.0;
			global_phaseFraction_limit = 1.0;
			nucleation_rate = 0.0;
			possible_volume = 0.0;
			nucleation_con = Default_X;
		}
	};

	struct NP_Node {
		vector<ConditionalPhase> conditional_phases;
		vector <pf::XNode> generate_x;
		vector <double> phaseFraction;
		int _x;
		int _y;
		int _z;
		NP_Node() {
			_x = 0;
			_y = 0;
			_z = 0;
		}
	};

	class NP_FieldStorage {
	public:
		NP_FieldStorage() {};
		~NP_FieldStorage() {};

		NP_Node& operator()(const int x, const int y, const int z);                                       //< Index operator for accessing the n's field value
		NP_FieldStorage& operator=(const NP_FieldStorage& n);
		int size() {
			return int(_mesh.size());
		}
		void init(int x, int y, int z, BoundaryCondition bc_x, BoundaryCondition bc_y, BoundaryCondition bc_z) {
			_mesh.resize(x * y * z);
			limit_num = 0;
			limit_x = x;
			limit_y = y;
			limit_z = z;
			x_bc = bc_x;
			y_bc = bc_y;
			z_bc = bc_z;
			for (int i = 0; i < x; i++)
				for (int j = 0; j < y; j++)
					for (int k = 0; k < z; k++) {
						(*this)(i, j, k)._x = i;
						(*this)(i, j, k)._y = j;
						(*this)(i, j, k)._z = k;
					}
		}
		void clear() {
			_mesh.clear();
		}
		int limit_x;
		int limit_y;
		int limit_z;
		int limit_num;
		BoundaryCondition x_bc;
		BoundaryCondition y_bc;
		BoundaryCondition z_bc;
		std::vector<NP_Node> _mesh;
	};
	inline NP_Node& NP_FieldStorage::operator()(const int x, const int y, const int z) {
		int _x = x, _y = y, _z = z;
		if (x >= 0 && x < limit_x && y >= 0 && y < limit_y && z >= 0 && z < limit_z) {
			return _mesh[x + y * limit_x + z * limit_x * limit_y];
		}
		if (x < 0)
			if (x_bc == PERIODIC)
				_x = _x + limit_x;
			else
				_x = 0;
		if (x >= limit_x)
			if (x_bc == PERIODIC)
				_x = _x - limit_x;
			else
				_x = limit_x - 1;
		if (y < 0)
			if (y_bc == PERIODIC)
				_y = _y + limit_y;
			else
				_y = 0;
		if (y >= limit_y)
			if (y_bc == PERIODIC)
				_y = _y - limit_y;
			else
				_y = limit_y - 1;
		if (z < 0)
			if (z_bc == PERIODIC)
				_z = _z + limit_z;
			else
				_z = 0;
		if (z >= limit_z)
			if (z_bc == PERIODIC)
				_z = _z - limit_z;
			else
				_z = limit_z - 1;
		return (*this)(_x, _y, _z);
	}
	inline NP_FieldStorage& NP_FieldStorage::operator=(const NP_FieldStorage& n) {
		_mesh = n._mesh;
		limit_x = n.limit_x;
		limit_y = n.limit_y;
		limit_z = n.limit_z;
		limit_num = n.limit_num;
		return *this;
	}

	struct NucleationBox {
	public:
		std::vector<pf::Nucleus> nucleus_box;
		std::vector<pf::GeometricRegion> geometricRegion_box;
		std::vector<pf::ConditionalPhase> condition_check_phase_box;
		std::vector<pf::PointSet> pointSet_box;
		///< Property
		int nucleation_property;
		///< ConditionalNucleation
		int nucleation_step;
		NucleationBox() {
			nucleation_property = NucleationProperty::DefiniteNucleation;
			nucleation_step = 500;
		};
		NucleationBox(NucleationProperty np) {
			nucleation_property = np;
			nucleation_step = 500;
		}
		NucleationBox& operator=(const NucleationBox& n) {
			nucleus_box = n.nucleus_box;
			geometricRegion_box = n.geometricRegion_box;
			condition_check_phase_box = n.condition_check_phase_box;
			pointSet_box = n.pointSet_box;
			nucleation_property = n.nucleation_property;
			nucleation_step = n.nucleation_step;
			return *this;
		}
	};
}
