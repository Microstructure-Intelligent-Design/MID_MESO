#pragma once
#include "sysTool.h"

using namespace std;
namespace pf {
    class Vector3;
    class Matrix3x3;
    class Vector6;
    class Matrix6x6;
    class vStrain;
    class vStress;
    class Matrix3x3
    {
    public:

        Matrix3x3() {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    storage[i][j] = 0.0;
        };

        Matrix3x3(const Matrix3x3& rhs) 
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    storage[i][j] = rhs(i, j);
        };

        double& operator()(const int i, const int j)
        {
#ifdef DEBUG
            if (i > 2 or j > 2 or i < 0 or j < 0)
            {
                std::cout << "Error in dMatrix3x3::operator()\n"
                    << "Access beyond storage range. (i, j) = "
                    << i << ", " << j << " > (2, 2)"
                    << "\nTerminating!!!" << std::endl;
                SYS_PROGRAM_STOP;
            }
#endif
            return storage[i][j];
        };
        double const& operator()(const int i, const int j) const
        {
#ifdef DEBUG
            if (i > 2 or j > 2 or i < 0 or j < 0)
            {
                std::cout << "Error in dMatrix3x3::operator()\n"
                    << "Access beyond storage range. (i, j) = "
                    << i << ", " << j << " > (2, 2)"
                    << "\nTerminating!!!" << std::endl;
                SYS_PROGRAM_STOP;
            }
#endif
            return storage[i][j];
        };
        void set_to_zero(void)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    storage[i][j] = 0.0;
        };
        void set_to_unity(void)
        {
            set_to_zero();
            storage[0][0] = 1.0;
            storage[1][1] = 1.0;
            storage[2][2] = 1.0;
        };
        Matrix3x3& operator=(const Matrix3x3& rhs)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    storage[i][j] = rhs(i, j);
            return *this;
        };
        Matrix3x3& operator=(const double rhs[3][3])
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    storage[i][j] = rhs[i][j];
            return *this;
        };
        bool operator==(Matrix3x3& rhs)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    if (storage[i][j] != rhs(i, j))
                        return false;
                }
            return true;
        };
        bool operator!=(Matrix3x3& rhs)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    if (storage[i][j] != rhs(i, j))
                        return true;
                }
            return false;
        };

        Matrix3x3 operator*(double m)
        {
            Matrix3x3 tmp;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp(i, j) = storage[i][j] * m;
                }
            return tmp;
        };
        Matrix3x3 operator/(double m)
        {
            Matrix3x3 tmp;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp(i, j) = storage[i][j] / m;
                }
            return tmp;
        };
        Vector3 operator*(const Vector3& rhs) const;
        Matrix3x3 operator*(const Matrix3x3& rhs) 
        {
            Matrix3x3 tmp;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                    {
                        tmp(i, j) += storage[i][k] * rhs(k, j);
                    }
            return tmp;
        };
        Matrix3x3 operator+(Matrix3x3& rhs)
        {
            Matrix3x3 tmp;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp(i, j) = storage[i][j] + rhs(i, j);
                }
            return tmp;
        };
        Matrix3x3 operator-(Matrix3x3& rhs)
        {
            Matrix3x3 tmp;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp(i, j) = storage[i][j] - rhs(i, j);
                }
            return tmp;
        };

        Matrix6x6 outer(const Matrix3x3& rhs) const;

        Matrix3x3& operator+=(Matrix3x3& rhs)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    storage[i][j] += rhs(i, j);
                }
            return *this;
        };
        Matrix3x3& operator-=(Matrix3x3& rhs)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    storage[i][j] -= rhs(i, j);
                }
            return *this;
        };
        Matrix3x3& operator*=(double m)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    storage[i][j] *= m;
                }
            return *this;
        };
        Matrix3x3& operator/=(double m)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    storage[i][j] /= m;
                }
            return *this;
        };

        inline double cal_determinant(void)
        {
            return (storage[0][0] * storage[1][1] * storage[2][2] +
                storage[0][1] * storage[1][2] * storage[2][0] +
                storage[0][2] * storage[1][0] * storage[2][1] -
                storage[0][2] * storage[1][1] * storage[2][0] -
                storage[0][1] * storage[1][0] * storage[2][2] -
                storage[0][0] * storage[1][2] * storage[2][1]);
        };
        Matrix3x3& do_invert(void)
        {
            Matrix3x3 tmp;
            double detInv = cal_determinant();

            if (detInv != 0.0) detInv = 1.0 / detInv;
            else
            {
                std::cout << "dMatrix3x3:\n" << this->print_matrix() << "Is not Invertible!" << std::endl;
                SYS_PROGRAM_STOP;
            }

            tmp(0, 0) = (storage[1][1] * storage[2][2] -
                storage[1][2] * storage[2][1]) * detInv;
            tmp(1, 0) = -(storage[1][0] * storage[2][2] -
                storage[1][2] * storage[2][0]) * detInv;
            tmp(2, 0) = (storage[1][0] * storage[2][1] -
                storage[1][1] * storage[2][0]) * detInv;
            tmp(0, 1) = -(storage[0][1] * storage[2][2] -
                storage[0][2] * storage[2][1]) * detInv;
            tmp(1, 1) = (storage[0][0] * storage[2][2] -
                storage[0][2] * storage[2][0]) * detInv;
            tmp(2, 1) = -(storage[0][0] * storage[2][1] -
                storage[0][1] * storage[2][0]) * detInv;
            tmp(0, 2) = (storage[0][1] * storage[1][2] -
                storage[1][1] * storage[0][2]) * detInv;
            tmp(1, 2) = -(storage[0][0] * storage[1][2] -
                storage[0][2] * storage[1][0]) * detInv;
            tmp(2, 2) = (storage[0][0] * storage[1][1] -
                storage[0][1] * storage[1][0]) * detInv;

            (*this) = tmp;
            return *this;
        };
        Matrix3x3 get_inverted_Matrix(void)
        {
            Matrix3x3 tmp;

            double detInv = cal_determinant();

            if (detInv != 0.0) detInv = 1.0 / detInv;
            else
            {
                std::cout << "dMatrix3x3:\n" << this->print_matrix() << "Is not Invertible!" << std::endl;
                SYS_PROGRAM_STOP;
            }

            tmp(0, 0) = (storage[1][1] * storage[2][2] -
                storage[1][2] * storage[2][1]) * detInv;
            tmp(1, 0) = -(storage[1][0] * storage[2][2] -
                storage[1][2] * storage[2][0]) * detInv;
            tmp(2, 0) = (storage[1][0] * storage[2][1] -
                storage[1][1] * storage[2][0]) * detInv;
            tmp(0, 1) = -(storage[0][1] * storage[2][2] -
                storage[0][2] * storage[2][1]) * detInv;
            tmp(1, 1) = (storage[0][0] * storage[2][2] -
                storage[0][2] * storage[2][0]) * detInv;
            tmp(2, 1) = -(storage[0][0] * storage[2][1] -
                storage[0][1] * storage[2][0]) * detInv;
            tmp(0, 2) = (storage[0][1] * storage[1][2] -
                storage[1][1] * storage[0][2]) * detInv;
            tmp(1, 2) = -(storage[0][0] * storage[1][2] -
                storage[0][2] * storage[1][0]) * detInv;
            tmp(2, 2) = (storage[0][0] * storage[1][1] -
                storage[0][1] * storage[1][0]) * detInv;

            return tmp;
        };
        Matrix3x3& do_transpose(void)
        {
            double tmp[3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[i][j] = storage[j][i];
                }
            (*this) = tmp;
            return *this;
        };
        Matrix3x3 get_transposed_Matrix(void)
        {
            Matrix3x3 TempMat;

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    TempMat(i, j) = storage[j][i];
                }
            return TempMat;
        };
        Matrix3x3& do_rotate(Matrix3x3& RotationMatrix)
        {
            Matrix3x3 TempMat;

            TempMat = (*this);

            (*this) = RotationMatrix * (TempMat * RotationMatrix.get_transposed_Matrix());

            return *this;
        };
        Matrix3x3 get_rotated_Matrix(Matrix3x3& RotationMatrix)
        {
            Matrix3x3 Out;
            Out = (*this);

            Out = RotationMatrix * (Out * RotationMatrix.get_transposed_Matrix());

            return Out;
        };
        double double_contract(Matrix3x3& rHS) const
        {
            double tmp = 0.0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp += storage[i][j] * rHS(i, j);
                }
            return tmp;
        };
        double trace(void) const
        {
            return storage[0][0] + storage[1][1] + storage[2][2];
        };
        Matrix3x3 get_sym(void) const
        {
            Matrix3x3 TempMat;

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    TempMat(i, j) = storage[i][j] + storage[j][i];
                }
            TempMat *= 0.5;

            return TempMat;
        };
        Matrix3x3 get_skew(void) const
        {
            Matrix3x3 TempMat;

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    TempMat(i, j) = storage[i][j] - storage[j][i];
                }
            TempMat *= 0.5;

            return TempMat;
        };
        double norm(void) const  /// Frobenius norm
        {
            double tmp = 0.0;
            for (int i = 0; i < 3; i++)
            {
                tmp += storage[i][0] * storage[i][0]
                    + storage[i][1] * storage[i][1]
                    + storage[i][2] * storage[i][2];
            }

            return sqrt(tmp);
        };
        std::string print_matrix(void) const
        {
            std::stringstream out;
            for (int i = 0; i < 3; i++)
            {
                out << "||" << std::setprecision(6) << std::right
                    << std::setw(8) << storage[i][0] << " "
                    << std::setw(8) << storage[i][1] << " "
                    << std::setw(8) << storage[i][2] << "||\n";
            }
            return out.str();
        };

        inline Vector6 VoigtVector() const;
        inline vStrain  VoigtStrain() const;
        inline vStress  VoigtStress() const;

        double storage[3][3];
    };

    class Vector3
    {
    public:
        Vector3() {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
        };
        Vector3(double i, double j, double k) {
            storage[0] = i;
            storage[1] = j;
            storage[2] = k;
        }

        Vector3(const Vector3& rhs)
        {
            for (int i = 0; i < 3; i++)
                storage[i] = rhs[i];
        };

        Vector3(const double vecinit[3])
        {
#ifdef DEBUG
            if (vecinit.size() != 3)
            {
                std::cout << "Error in dVector3::constructor()\n"
                    << "Initialization list size beyond storage range."
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            storage[0] = vecinit[0];
            storage[1] = vecinit[1];
            storage[2] = vecinit[2];
        }

        double& operator[](const int i)
        {
#ifdef DEBUG
            if (i > 2)
            {
                std::cout << "Error in dVector3::operator[]\n"
                    << "Access beyond storage range. i = "
                    << i << " > 2"
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            return storage[i];
        };
        double const& operator[](const int i) const
        {
#ifdef DEBUG
            if (i > 2)
            {
                std::cout << "Error in dVector3::operator[]\n"
                    << "Access beyond storage range. i = "
                    << i << " > 2"
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            return storage[i];
        };

        double getX(void) const
        {
            return storage[0];
        };

        void setX(const double newX)
        {
            storage[0] = newX;
        };

        double getY(void) const
        {
            return storage[1];
        };

        void setY(const double newY)
        {
            storage[1] = newY;
        };

        double getZ(void) const
        {
            return storage[2];
        };

        void setZ(const double newX)
        {
            storage[2] = newX;
        };

        void set_to_zero(void)
        {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
        };
        void set_to_unitX(void)
        {
            storage[0] = 1.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
        };
        void set_to_unitY(void)
        {
            storage[0] = 0.0;
            storage[1] = 1.0;
            storage[2] = 0.0;
        };
        void set_to_unitZ(void)
        {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 1.0;
        }; 
        bool operator==(const Vector3& rhs) const
        {
			for (int i = 0; i < 3; i++) {
				if (storage[i] != rhs[i])
					return false;
			}
            return true;
        };
        bool operator!=(const Vector3& rhs) const
        {
            for (int i = 0; i < 3; i++) {
                if (storage[i] != rhs[i])
                    return true;
            }
            return false;
        };
        Vector3 operator*(const double m) const
        {
            Vector3 tmp;
            tmp[0] = storage[0] * m;
            tmp[1] = storage[1] * m;
            tmp[2] = storage[2] * m;
            return tmp;
        };
        Vector3 operator/(const double m) const
        {
            Vector3 tmp;
            tmp[0] = storage[0] / m;
            tmp[1] = storage[1] / m;
            tmp[2] = storage[2] / m;
            return tmp;
        };
        double operator*(const Vector3& rhs) const
        {
            return storage[0] * rhs[0] + storage[1] * rhs[1] + storage[2] * rhs[2];
        };
        double abs() const
        {
            return sqrt(storage[0] * storage[0] +
                storage[1] * storage[1] +
                storage[2] * storage[2]);
        };
        Vector3 cross(const Vector3& rhs) const
        {
            Vector3 tmp;
            tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
            tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
            tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
            return tmp;
        };
        Vector3 operator+(const Vector3& rhs) const
        {
            Vector3 tmp;
            tmp[0] = storage[0] + rhs[0];
            tmp[1] = storage[1] + rhs[1];
            tmp[2] = storage[2] + rhs[2];
            return tmp;
        };
        Vector3 operator-(const Vector3& rhs) const
        {
            Vector3 tmp;
            tmp[0] = storage[0] - rhs[0];
            tmp[1] = storage[1] - rhs[1];
            tmp[2] = storage[2] - rhs[2];
            return tmp;
        };
        Vector3& operator*=(const double m)
        {
            storage[0] *= m;
            storage[1] *= m;
            storage[2] *= m;
            return *this;
        };
        Vector3& operator/=(const double m)
        {
            storage[0] /= m;
            storage[1] /= m;
            storage[2] /= m;
            return *this;
        };
        Vector3& operator-=(const Vector3& rhs)
        {
            storage[0] = storage[0] - rhs[0];
            storage[1] = storage[1] - rhs[1];
            storage[2] = storage[2] - rhs[2];
            return *this;
        };
        Vector3& operator+=(const Vector3& rhs)
        {
            storage[0] = storage[0] + rhs[0];
            storage[1] = storage[1] + rhs[1];
            storage[2] = storage[2] + rhs[2];
            return *this;
        };
        Vector3& operator=(const Vector3& rhs)
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            return *this;
        };
        Vector3& operator=(const double rhs[3])
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            return *this;
        };
        double length(void) const
        {
            return sqrt(storage[0] * storage[0] +
                storage[1] * storage[1] +
                storage[2] * storage[2]);
        };
        Vector3& normalize(void)
        {
            double norm = sqrt(storage[0] * storage[0] +
                storage[1] * storage[1] +
                storage[2] * storage[2]);
            if (norm != 0.0)
            {
                storage[0] /= norm;
                storage[1] /= norm;
                storage[2] /= norm;
            }
            return *this;
        };
        Vector3 normalized(void) const
        {
            double norm = sqrt(storage[0] * storage[0] +
                storage[1] * storage[1] +
                storage[2] * storage[2]);
            Vector3 tmp;
            if (norm != 0.0)
            {
                tmp[0] = storage[0] / norm;
                tmp[1] = storage[1] / norm;
                tmp[2] = storage[2] / norm;
            }
            else
            {
                tmp.set_to_zero();
            }
            return tmp;
        };
        Vector3& do_rotate(const Matrix3x3& RotationMatrix)
        {
            Vector3 Out;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                {
                    Out[i] += RotationMatrix(i, j) * storage[j];
                }
            (*this) = Out;
            return *this;
        };
        Vector3 get_rotated_vec3(const Matrix3x3& RotationMatrix) const
        {
            Vector3 Out;
            Out.set_to_zero();
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j){
                    Out[i] += RotationMatrix(i, j) * storage[j];
                }
            return Out;
        };
        Vector3 Xreflected(void) const
        {
            Vector3 Out;
            for (int i = 0; i < 3; ++i)
            {
                Out[i] = storage[i];
            }
            Out[0] *= -1.0;
            return Out;
        };
        Vector3 Yreflected(void) const
        {
            Vector3 Out;
            for (int i = 0; i < 3; ++i)
            {
                Out[i] = storage[i];
            }
            Out[1] *= -1.0;
            return Out;
        };
        Vector3 Zreflected(void) const
        {
            Vector3 Out;
            for (int i = 0; i < 3; ++i)
            {
                Out[i] = storage[i];
            }
            Out[2] *= -1.0;
            return Out;
        };
        std::string print(void) const
        {
            std::stringstream out;

            out << "(" << storage[0] << ", "
                << storage[1] << ", "
                << storage[2] << ")";
            return out.str();
        };
        double* data(void)
        {
            return storage;
        };
        const double* const_data(void) const
        {
            return storage;
        };
    protected:
    private:
        double storage[3];
    };

    inline Vector3 Matrix3x3::operator*(const Vector3& rhs) const
    {
        Vector3 tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    }

    struct vec3_elem {
        int index;
        Vector3 vec;
        vec3_elem() {
            index = 0;
            vec[0] = 0.0;
            vec[1] = 0.0;
            vec[2] = 0.0;
        }
        vec3_elem& operator=(const vec3_elem& n) {
            index = n.index;
            vec[0] = n.vec[0];
            vec[1] = n.vec[1];
            vec[2] = n.vec[2];
            return *this;
        }
    };

    struct vec3_box {
        vector<vec3_elem> _vec_box;
        typedef std::vector<vec3_elem>::iterator iterator;
        typedef std::vector<vec3_elem>::const_iterator citerator;
        iterator  begin() { return _vec_box.begin(); };
        iterator  end() { return _vec_box.end(); };
        Vector3& operator[](const int index) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                if (i->index == index) return i->vec;
            }
            cout << "vec3_box error, can't find the vec index : " << index << endl;
            SYS_PROGRAM_STOP;
        }
        double& operator()(const int index, const int index2) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                if (i->index == index) return i->vec[index2];
            }
            cout << "vec3_box error, can't find the vec index : " << index << endl;
            SYS_PROGRAM_STOP;
        }
        vec3_box& operator=(const vec3_box& n) {
            _vec_box = n._vec_box;
            return *this;
        }
        void set_value(double value = 0.0) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i) {
                i->vec[0] = value;
                i->vec[1] = value;
                i->vec[2] = value;
            }
        }
        void set_value(int _index, double value = 0.0) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = value;
                    i->vec[1] = value;
                    i->vec[2] = value;
                }
        }
        void add_vec(int _index, double* vec) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = vec[0];
                    i->vec[1] = vec[1];
                    i->vec[2] = vec[2];
                    return;
                }
            vec3_elem elem;
            elem.index = _index;
            elem.vec[0] = vec[0];
            elem.vec[1] = vec[1];
            elem.vec[2] = vec[2];
            _vec_box.push_back(elem);
        }
        void add_vec(int _index, Vector3 vec) {
            for (auto i = _vec_box.begin(); i < _vec_box.end(); ++i)
                if (i->index == _index) {
                    i->vec[0] = vec[0];
                    i->vec[1] = vec[1];
                    i->vec[2] = vec[2];
                    return;
                }
            vec3_elem elem;
            elem.index = _index;
            elem.vec[0] = vec[0];
            elem.vec[1] = vec[1];
            elem.vec[2] = vec[2];
            _vec_box.push_back(elem);
        }
        void clear() {
            _vec_box.clear();
        }
        int size() {
            return int(_vec_box.size());
        }
    };

    class Vector6
    {
    public:

        Vector6() {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
            storage[3] = 0.0;
            storage[4] = 0.0;
            storage[5] = 0.0;
        }

        Vector6(std::initializer_list<double> vecinit)
        {
#ifdef DEBUG
            if (vecinit.size() != 6)
            {
                std::cout << "Error in Vector6::constructor()\n"
                    << "Initialization list size beyond storage range."
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            int ii = 0;
            for (auto it = vecinit.begin(); it != vecinit.end(); it++)
            {
                storage[ii] = *it;
                ii += 1;
            }
        }

        double& operator[](const int i)
        {
#ifdef DEBUG
            if (i > 5)
            {
                std::cout << "Error in Vector6::operator[]\n"
                    << "Access beyond storage range. i = "
                    << i << " > 5"
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            return storage[i];
        };
        double const& operator[](const int i) const
        {
#ifdef DEBUG
            if (i > 5)
            {
                std::cout << "Error in dVector6::operator[]\n"
                    << "Access beyond storage range. i = "
                    << i << " > 5"
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            return storage[i];
        };
        bool is_nan_val_exist() {
            for (int i = 0; i < 6; i++)
				if (IS_NAN(storage[i]))
					return true;
            return false;
        }
        void set_to_zero(void)
        {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
            storage[3] = 0.0;
            storage[4] = 0.0;
            storage[5] = 0.0;
        };
        void set_to_unity(void)
        {
            storage[0] = 1.0;
            storage[1] = 1.0;
            storage[2] = 1.0;
            storage[3] = 0.0;
            storage[4] = 0.0;
            storage[5] = 0.0;
        };
        void do_abs(void)
        {
            storage[0] = abs(storage[0]);
            storage[1] = abs(storage[1]);
            storage[2] = abs(storage[2]);
            storage[3] = abs(storage[3]);
            storage[4] = abs(storage[4]);
            storage[5] = abs(storage[5]);
        };
        double trace(void) const
        {
            return storage[0] + storage[1] + storage[2];
        };
        double max_abs(void) const
        {
            // Returns maximum absolute value
            double tempmax = 0.0;
            for (int i = 0; i < 6; i++)
            {
                if (abs(storage[i]) > tempmax)
                {
                    tempmax = abs(storage[i]);
                }
            }
            return tempmax;
        };
        //    double norm(void) const  /// Frobenius norm
        //    {
        //        double tmp = storage[0] * storage[0]
        //                   + storage[1] * storage[1]
        //                   + storage[2] * storage[2]
        //                   + storage[3] * storage[3] * 2.0
        //                   + storage[4] * storage[4] * 2.0
        //                   + storage[5] * storage[5] * 2.0;
        //
        //        return sqrt(tmp);
        //    };
        Vector6 operator*(const double m) const
        {
            Vector6 tmp;
            tmp[0] = storage[0] * m;
            tmp[1] = storage[1] * m;
            tmp[2] = storage[2] * m;
            tmp[3] = storage[3] * m;
            tmp[4] = storage[4] * m;
            tmp[5] = storage[5] * m;
            return tmp;
        };
        bool operator==(const Vector6& rhs) const
        {
			for (int i = 0; i < 6; i++) {
				if (storage[i] != rhs[i])
					return false;
			}
            return true;
        };
        bool operator!=(const Vector6& rhs) const
        {
            for (int i = 0; i < 6; i++) {
                if (storage[i] != rhs[i])
                    return true;
            }
            return false;
        };
        Vector6 operator/(const double m) const
        {
            Vector6 tmp;
            tmp[0] = storage[0] / m;
            tmp[1] = storage[1] / m;
            tmp[2] = storage[2] / m;
            tmp[3] = storage[3] / m;
            tmp[4] = storage[4] / m;
            tmp[5] = storage[5] / m;
            return tmp;
        };
        Vector6& operator*=(const double m)
        {
            storage[0] *= m;
            storage[1] *= m;
            storage[2] *= m;
            storage[3] *= m;
            storage[4] *= m;
            storage[5] *= m;
            return *this;
        };

        double operator*(const Vector6& rhs)
        {
            double tmp = 0.0;
            for (int i = 0; i < 6; i++)
            {
                tmp += storage[i] * rhs[i];
            }
            return tmp;
        };

        Vector6& operator/=(const double m)
        {
            storage[0] /= m;
            storage[1] /= m;
            storage[2] /= m;
            storage[3] /= m;
            storage[4] /= m;
            storage[5] /= m;

            return *this;
        };
        Vector6 operator+(const Vector6& rhs) const
        {
            Vector6 tmp;
            tmp[0] = storage[0] + rhs[0];
            tmp[1] = storage[1] + rhs[1];
            tmp[2] = storage[2] + rhs[2];
            tmp[3] = storage[3] + rhs[3];
            tmp[4] = storage[4] + rhs[4];
            tmp[5] = storage[5] + rhs[5];
            return tmp;
        };
        Vector6& operator+=(const Vector6& rhs)
        {
            storage[0] = storage[0] + rhs[0];
            storage[1] = storage[1] + rhs[1];
            storage[2] = storage[2] + rhs[2];
            storage[3] = storage[3] + rhs[3];
            storage[4] = storage[4] + rhs[4];
            storage[5] = storage[5] + rhs[5];
            return *this;
        };
        Vector6 operator-(const Vector6& rhs) const
        {
            Vector6 tmp;
            tmp[0] = storage[0] - rhs[0];
            tmp[1] = storage[1] - rhs[1];
            tmp[2] = storage[2] - rhs[2];
            tmp[3] = storage[3] - rhs[3];
            tmp[4] = storage[4] - rhs[4];
            tmp[5] = storage[5] - rhs[5];
            return tmp;
        };
        Vector6& operator-=(const Vector6& rhs)
        {
            storage[0] = storage[0] - rhs[0];
            storage[1] = storage[1] - rhs[1];
            storage[2] = storage[2] - rhs[2];
            storage[3] = storage[3] - rhs[3];
            storage[4] = storage[4] - rhs[4];
            storage[5] = storage[5] - rhs[5];
            return *this;
        };
        Vector6& operator=(const Vector6& rhs)
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            storage[3] = rhs[3];
            storage[4] = rhs[4];
            storage[5] = rhs[5];
            return *this;
        };
        Vector6& operator=(const double rhs[6])
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            storage[3] = rhs[3];
            storage[4] = rhs[4];
            storage[5] = rhs[5];
            return *this;
        };


        //    double norm(void) const  /// Frobenius norm
        //    {
        //        double tmp = storage[0] * storage[0]
        //                   + storage[1] * storage[1]
        //                   + storage[2] * storage[2]
        //                   + storage[3] * storage[3] * 2.0
        //                   + storage[4] * storage[4] * 2.0
        //                   + storage[5] * storage[5] * 2.0;
        //        return sqrt(tmp);
        //    };

        Vector6 get_rotated_vec6(const Matrix3x3& RotationMatrix) const
        {
            double In[3][3];
            double Out[3][3];

            In[0][0] = storage[0];
            In[0][1] = storage[5];
            In[0][2] = storage[4];
            In[1][0] = storage[5];
            In[1][1] = storage[1];
            In[1][2] = storage[3];
            In[2][0] = storage[4];
            In[2][1] = storage[3];
            In[2][2] = storage[2];

            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                {
                    Out[p][q] = 0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                        {
                            Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                        }
                }
            Vector6 out;

            out[0] = Out[0][0];
            out[5] = Out[0][1];
            out[4] = Out[0][2];
            out[1] = Out[1][1];
            out[3] = Out[1][2];
            out[2] = Out[2][2];

            return out;
        };
        Vector6& do_rotate(Matrix3x3 RotationMatrix)
        {
            double In[3][3];
            double Out[3][3];

            In[0][0] = storage[0];
            In[0][1] = storage[5];
            In[0][2] = storage[4];
            In[1][0] = storage[5];
            In[1][1] = storage[1];
            In[1][2] = storage[3];
            In[2][0] = storage[4];
            In[2][1] = storage[3];
            In[2][2] = storage[2];

            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                {
                    Out[p][q] = 0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                        {
                            Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                        }
                }

            storage[0] = Out[0][0];
            storage[5] = Out[0][1];
            storage[4] = Out[0][2];
            storage[1] = Out[1][1];
            storage[3] = Out[1][2];
            storage[2] = Out[2][2];

            return *this;
        };
        double* data(void)
        {
            return storage;
        };
        const double* const_data(void) const
        {
            return storage;
        };
        Vector6& H_product(Vector6 vector)
        {
            storage[0] *= vector[0];
            storage[1] *= vector[1];
            storage[2] *= vector[2];
            storage[3] *= vector[3];
            storage[4] *= vector[4];
            storage[5] *= vector[5];

            return *this;
        };

        Matrix3x3 sym_tensor(void) const
        {
            Matrix3x3 tmp;
            tmp(0, 0) = storage[0];
            tmp(0, 1) = storage[5];
            tmp(0, 2) = storage[4];
            tmp(1, 0) = storage[5];
            tmp(1, 1) = storage[1];
            tmp(1, 2) = storage[3];
            tmp(2, 0) = storage[4];
            tmp(2, 1) = storage[3];
            tmp(2, 2) = storage[2];
            return tmp;
        };

        Matrix3x3 skew_tensor(void) const
        {
            Matrix3x3 tmp;
            tmp(0, 0) = storage[0];
            tmp(0, 1) = storage[5];
            tmp(0, 2) = storage[4];
            tmp(1, 0) = -storage[5];
            tmp(1, 1) = storage[1];
            tmp(1, 2) = storage[3];
            tmp(2, 0) = -storage[4];
            tmp(2, 1) = -storage[3];
            tmp(2, 2) = storage[2];
            return tmp;
        };

        std::string print(void) const
        {
            std::stringstream out;
            out << "< | ";
            for (int i = 0; i < 6; i++)
            {
                out << storage[i] << " " << " | ";
            }
            out << " >";
            return out.str();
        };
    protected:
        double storage[6];
    private:
    };

    class Matrix6x6
    {
    public:
        Matrix6x6() {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    storage[i][j] = 0.0;
        }
        double& operator()(const int i, const int j)
        {
#ifdef DEBUG
            if (i > 5 or j > 5)
            {
                std::cout << "Error in Matrix6x6::operator()\n"
                    << "Access beyond storage range. (i, j) = "
                    << i << ", " << j << " > (5, 5)"
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            return storage[i][j];
        };
        bool is_nan_val_exist() {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (IS_NAN(storage[i][j]))
                        return true;
            return false;
        }
        double const& operator()(const int i, const int j) const
        {
#ifdef DEBUG
            if (i > 5 or j > 5)
            {
                std::cout << "Error in Matrix6x6::operator()\n"
                    << "Access beyond storage range. (i, j) = "
                    << i << ", " << j << " > (5, 5)"
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
#endif
            return storage[i][j];
        };
        void set_to_zero(void)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    storage[i][j] = 0.0;
        };
        void set_to_unity(void)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (i == j)
                        storage[i][j] = 1.0;
                    else
                        storage[i][j] = 0.0;
        };
        
        double norm(void) const
        {
            double tmp = 0.0;
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    tmp += storage[i][j] * storage[i][j];
                }
            return sqrt(tmp);
        };
        /*double det(void) const
        {
            double determinant = 0.0;

            for (int i = 0; i < 6; i++)
            {
                double line_product = 1.0;
                for (int j = 0; j < 6; j++)
                {
                    line_product *= storage[(i+j)%6][j];
                }
                determinant += line_product;
            }

            for (int i = 0; i < 6; i++)
            {
                double line_product = 1.0;
                for (int j = 0; j < 6; j++)
                {
                    line_product *= storage[(i-j+6)%6][j];
                }
                determinant -= line_product;
            }
            std::cout << determinant << std::endl;
            return determinant;
        };

        bool is_singular(void)
        {
            if(fabs(det()) > DBL_EPSILON)
            {
                return false;
            }
            else
            {
                return true;
            }
        };*/

        Matrix6x6 operator*(const double m) const
        {
            Matrix6x6 tmp;
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    tmp(i, j) = storage[i][j] * m;
                }
            return tmp;
        };
        bool operator==(const Matrix6x6& rhs) const
        {
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 6; j++) {
					if (storage[i][j] != rhs(i, j))
						return false;
				}
            return true;
        };
        bool operator!=(const Matrix6x6& rhs) const
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++) {
                    if (storage[i][j] != rhs(i, j))
                        return true;
                }
            return false;
        };
        Matrix6x6& operator*=(const double m)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    storage[i][j] *= m;
                }
            return *this;
        };
        /*Vector6 operator*(Vector6 rhs) const
        {
            Vector6 tmp;
            tmp.set_to_zero();
            for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
            {
                tmp[i] += storage[i][j]*rhs[j];
            }
            return tmp;
        };*/
        Matrix6x6 operator*(const Matrix6x6& rhs) const
        {
            Matrix6x6 tmp;
            tmp.set_to_zero();
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    for (int k = 0; k < 6; k++)
                    {
                        tmp(i, j) += storage[i][k] * rhs(k, j);
                    }
            return tmp;
        };

        vStrain operator*(const vStress& rhs) const;

        vStress operator*(const vStrain& rhs) const;

        Matrix6x6 operator+(const Matrix6x6& rhs) const
        {
            Matrix6x6 tmp;
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    tmp(i, j) = storage[i][j] + rhs(i, j);
                }
            return tmp;
        };
        Matrix6x6 operator-(const Matrix6x6& rhs) const
        {
            Matrix6x6 tmp;
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    tmp(i, j) = storage[i][j] - rhs(i, j);
                }
            return tmp;
        };
        Matrix6x6& operator+=(const Matrix6x6& rhs)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    storage[i][j] += rhs(i, j);
                }
            return *this;
        };
        Matrix6x6& operator/=(const double m)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    storage[i][j] /= m;
                }
            return *this;
        };
        Matrix6x6& operator-=(const Matrix6x6& rhs)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    storage[i][j] -= rhs(i, j);
                }
            return *this;
        };
        Matrix6x6& operator=(const Matrix6x6& rhs)
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    storage[i][j] = rhs(i, j);
            return *this;
        };
        Matrix6x6& operator=(const double rhs[6][6])
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    storage[i][j] = rhs[i][j];
            return *this;
        };
        Matrix6x6& do_invert(void)
        {
            double Out[6][6];

            int indxc[6];
            int indxr[6];
            int ipiv[6] = { 0, 0, 0, 0, 0, 0 };
            int icol = 0;
            int irow = 0;
            double pivinv;
            double dum;
            double Uni[6][6] = { {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} };


            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    Out[i][j] = storage[i][j];
                }

            for (int i = 0; i < 6; i++)
            {
                double big = 0.0;
                for (int j = 0; j < 6; j++)
                    if (ipiv[j] != 1)
                        for (int k = 0; k < 6; k++)
                        {
                            if (ipiv[k] == 0)
                            {
                                if (fabs(Out[j][k]) >= big)
                                {
                                    big = fabs(Out[j][k]);
                                    irow = j;
                                    icol = k;
                                };
                            }
                            else if (ipiv[k] > 1)
                            {
                                std::cout << "Matrix6x6: Can Not Compute Inverse Matrix."
                                    << "Matrix:\n" << this->print()
                                    << "is Singular 1!!!" << std::endl;
                                SYS_PROGRAM_STOP;
                            }
                        };
                ++(ipiv[icol]);
                if (irow != icol)
                {
                    for (int l = 0; l < 6; l++)
                    {
                        double temp = Out[irow][l];
                        Out[irow][l] = Out[icol][l];
                        Out[icol][l] = temp;
                    };
                    for (int l = 0; l < 6; l++)
                    {
                        double temp = Uni[irow][l];
                        Uni[irow][l] = Uni[icol][l];
                        Uni[icol][l] = temp;
                    };
                };
                indxr[i] = irow;
                indxc[i] = icol;
                if (fabs(Out[icol][icol]) <= DBL_EPSILON)
                {
                    std::cout << "dMatrix6x6: Can Not Compute Inverse Matrix. Matrix:\n" << this->print() << "is Singular 2!!!" << std::endl;
                    SYS_PROGRAM_STOP;
                }
                pivinv = 1.0 / Out[icol][icol];
                Out[icol][icol] = 1.0;
                for (int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
                for (int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
                for (int ll = 0; ll < 6; ll++)
                    if (ll != icol)
                    {
                        dum = Out[ll][icol];
                        Out[ll][icol] = 0.0;
                        for (int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l] * dum;
                        for (int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l] * dum;
                    }
            }
            for (int l = 5; l >= 0; l--)
            {
                if (indxr[l] != indxc[l])
                    for (int k = 0; k < 6; k++)
                    {
                        double temp = Out[k][indxr[l]];
                        Out[k][indxr[l]] = Out[k][indxc[l]];
                        Out[k][indxc[l]] = temp;
                    };
            }

            (*this) = Out;
            return *this;
        };

        Matrix6x6 get_inverted_matrix(void) const
        {
            double Out[6][6];

            int indxc[6];
            int indxr[6];
            int ipiv[6] = { 0, 0, 0, 0, 0, 0 };
            int icol = 0;
            int irow = 0;
            double pivinv;
            double dum;
            double Uni[6][6] = { {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0} };

            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    Out[i][j] = storage[i][j];
                }

            for (int i = 0; i < 6; i++)
            {
                double big = 0.0;
                for (int j = 0; j < 6; j++)
                    if (ipiv[j] != 1)
                        for (int k = 0; k < 6; k++)
                        {
                            if (ipiv[k] == 0)
                            {
                                if (fabs(Out[j][k]) >= big)
                                {
                                    big = fabs(Out[j][k]);
                                    irow = j;
                                    icol = k;
                                };
                            }
                            else if (ipiv[k] > 1)
                            {
                                std::cout << "Matrix6x6: Can Not Compute Inverse Matrix.\n"
                                    << this->print() << "Matrix is Singular 2!!!"
                                    << std::endl;
                                SYS_PROGRAM_STOP;
                            }
                        };
                ++(ipiv[icol]);
                if (irow != icol)
                {
                    for (int l = 0; l < 6; l++)
                    {
                        double temp = Out[irow][l];
                        Out[irow][l] = Out[icol][l];
                        Out[icol][l] = temp;
                    };
                    for (int l = 0; l < 6; l++)
                    {
                        double temp = Uni[irow][l];
                        Uni[irow][l] = Uni[icol][l];
                        Uni[icol][l] = temp;
                    };
                };
                indxr[i] = irow;
                indxc[i] = icol;
                if (fabs(Out[icol][icol]) <= DBL_EPSILON)
                {
                    std::cout << "Matrix6x6:\n" << this->print()
                        << " Is not Invertible!" << std::endl;
                    SYS_PROGRAM_STOP;
                }
                pivinv = 1.0 / Out[icol][icol];
                Out[icol][icol] = 1.0;
                for (int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
                for (int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
                for (int ll = 0; ll < 6; ll++)
                    if (ll != icol)
                    {
                        dum = Out[ll][icol];
                        Out[ll][icol] = 0.0;
                        for (int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l] * dum;
                        for (int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l] * dum;
                    }
            }
            for (int l = 5; l >= 0; l--)
            {
                if (indxr[l] != indxc[l])
                    for (int k = 0; k < 6; k++)
                    {
                        double temp = Out[k][indxr[l]];
                        Out[k][indxr[l]] = Out[k][indxc[l]];
                        Out[k][indxc[l]] = temp;
                    };
            }

            Matrix6x6 tmp;
            tmp = Out;
            return tmp;
        };

        Matrix6x6& do_transpose(void)
        {
            double tmp[6][6];
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    tmp[i][j] = storage[j][i];
                }
            (*this) = tmp;
            return *this;
        };

        Matrix6x6 get_transposed_matrix(void) const
        {
            Matrix6x6 tmp;

            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    tmp(i, j) = storage[j][i];
                }
            return tmp;
        };

        Matrix6x6& do_rotate(const Matrix3x3& RotationMatrix)
        {
            double In[3][3][3][3];
            double Out[3][3][3][3];

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            In[i][j][k][l] = tensor(i, j, k, l);
                        }

            for (int m = 0; m < 3; m++)
                for (int n = 0; n < 3; n++)
                    for (int p = 0; p < 3; p++)
                        for (int q = 0; q < 3; q++)
                        {
                            Out[m][n][p][q] = 0.0;

                            for (int i = 0; i < 3; i++)
                                for (int j = 0; j < 3; j++)
                                    for (int k = 0; k < 3; k++)
                                        for (int l = 0; l < 3; l++)
                                        {
                                            //                Former version
                                            //                Out[m][n][p][q] += RotationMatrix(i,m)*
                                            //                                   RotationMatrix(j,n)*
                                            //                                   In[i][j][k][l]*
                                            //                                   RotationMatrix(k,p)*
                                            //                                   RotationMatrix(l,q);
                                            //              New version (active rotation)
                                            Out[m][n][p][q] += RotationMatrix(m, i) *
                                                RotationMatrix(n, j) *
                                                In[i][j][k][l] *
                                                RotationMatrix(p, k) *
                                                RotationMatrix(q, l);
                                        }
                        }
            int VoigtIndex[6][2] = { {0,0},{1,1},{2,2},{1,2},{0,2},{0,1} };
            for (int m = 0; m < 6; m++)
                for (int n = 0; n < 6; n++)
                {
                    int i = VoigtIndex[m][0];
                    int j = VoigtIndex[m][1];
                    int k = VoigtIndex[n][0];
                    int l = VoigtIndex[n][1];

                    storage[m][n] = Out[i][j][k][l];
                }
            return *this;
        };
        Matrix6x6 get_rotated_matrix(const Matrix3x3& RotationMatrix)
        {
            double In[3][3][3][3];
            double Out[3][3][3][3];

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            In[i][j][k][l] = tensor(i, j, k, l);
                        }

            for (int m = 0; m < 3; m++)
                for (int n = 0; n < 3; n++)
                    for (int p = 0; p < 3; p++)
                        for (int q = 0; q < 3; q++)
                        {
                            Out[m][n][p][q] = 0.0;

                            for (int i = 0; i < 3; i++)
                                for (int j = 0; j < 3; j++)
                                    for (int k = 0; k < 3; k++)
                                        for (int l = 0; l < 3; l++)
                                        {
                                            //                Former version
                                            //                Out[m][n][p][q] += RotationMatrix(i,m)*
                                            //                                   RotationMatrix(j,n)*
                                            //                                   In[i][j][k][l]*
                                            //                                   RotationMatrix(k,p)*
                                            //                                   RotationMatrix(l,q);
                                            //              New version (active rotation)
                                            Out[m][n][p][q] += RotationMatrix(m, i) *
                                                RotationMatrix(n, j) *
                                                In[i][j][k][l] *
                                                RotationMatrix(p, k) *
                                                RotationMatrix(q, l);
                                        }
                        }
            Matrix6x6 _OUT;
            int VoigtIndex[6][2] = { {0,0},{1,1},{2,2},{1,2},{0,2},{0,1} };

            for (int m = 0; m < 6; m++)
                for (int n = 0; n < 6; n++)
                {
                    int i = VoigtIndex[m][0];
                    int j = VoigtIndex[m][1];
                    int k = VoigtIndex[n][0];
                    int l = VoigtIndex[n][1];

                    _OUT(m, n) = Out[i][j][k][l];
                }
            return _OUT;
        };

        std::string print(void) const
        {
            std::stringstream out;
            for (int i = 0; i < 6; i++)
            {
                out << "||" << std::setprecision(6)
                    << std::setw(8) << storage[i][0] << " "
                    << std::setw(8) << storage[i][1] << " "
                    << std::setw(8) << storage[i][2] << " "
                    << std::setw(8) << storage[i][3] << " "
                    << std::setw(8) << storage[i][4] << " "
                    << std::setw(8) << storage[i][5] << "||\n";
            }
            return out.str();
        };
        const double& tensor(const int i, const int j, const int k, const int l) const
        {
            return storage[(i == j) ? (i) : (6 - (i + j))][(k == l) ? (k) : (6 - (k + l))];
        };
        double* data(void)
        {
            return &storage[0][0];
        };
        const double* const_data(void) const
        {
            return &storage[0][0];
        };
    protected:
    private:

        double storage[6][6];
    };

    inline Matrix6x6 Matrix3x3::outer(const Matrix3x3& rhs) const
    {
        Matrix6x6 tmp;

        tmp(0, 0) += storage[0][0] * rhs(0, 0);
        tmp(0, 5) += storage[0][0] * rhs(0, 1);
        tmp(0, 4) += storage[0][0] * rhs(0, 2);
        tmp(0, 5) += storage[0][0] * rhs(1, 0);
        tmp(0, 1) += storage[0][0] * rhs(1, 1);
        tmp(0, 3) += storage[0][0] * rhs(1, 2);
        tmp(0, 4) += storage[0][0] * rhs(2, 0);
        tmp(0, 3) += storage[0][0] * rhs(2, 1);
        tmp(0, 2) += storage[0][0] * rhs(2, 2);
        tmp(5, 0) += storage[0][1] * rhs(0, 0);
        tmp(5, 5) += storage[0][1] * rhs(0, 1);
        tmp(5, 4) += storage[0][1] * rhs(0, 2);
        tmp(5, 5) += storage[0][1] * rhs(1, 0);
        tmp(5, 1) += storage[0][1] * rhs(1, 1);
        tmp(5, 3) += storage[0][1] * rhs(1, 2);
        tmp(5, 4) += storage[0][1] * rhs(2, 0);
        tmp(5, 3) += storage[0][1] * rhs(2, 1);
        tmp(5, 2) += storage[0][1] * rhs(2, 2);
        tmp(4, 0) += storage[0][2] * rhs(0, 0);
        tmp(4, 5) += storage[0][2] * rhs(0, 1);
        tmp(4, 4) += storage[0][2] * rhs(0, 2);
        tmp(4, 5) += storage[0][2] * rhs(1, 0);
        tmp(4, 1) += storage[0][2] * rhs(1, 1);
        tmp(4, 3) += storage[0][2] * rhs(1, 2);
        tmp(4, 4) += storage[0][2] * rhs(2, 0);
        tmp(4, 3) += storage[0][2] * rhs(2, 1);
        tmp(4, 2) += storage[0][2] * rhs(2, 2);
        tmp(5, 0) += storage[1][0] * rhs(0, 0);
        tmp(5, 5) += storage[1][0] * rhs(0, 1);
        tmp(5, 4) += storage[1][0] * rhs(0, 2);
        tmp(5, 5) += storage[1][0] * rhs(1, 0);
        tmp(5, 1) += storage[1][0] * rhs(1, 1);
        tmp(5, 3) += storage[1][0] * rhs(1, 2);
        tmp(5, 4) += storage[1][0] * rhs(2, 0);
        tmp(5, 3) += storage[1][0] * rhs(2, 1);
        tmp(5, 2) += storage[1][0] * rhs(2, 2);
        tmp(1, 0) += storage[1][1] * rhs(0, 0);
        tmp(1, 5) += storage[1][1] * rhs(0, 1);
        tmp(1, 4) += storage[1][1] * rhs(0, 2);
        tmp(1, 5) += storage[1][1] * rhs(1, 0);
        tmp(1, 1) += storage[1][1] * rhs(1, 1);
        tmp(1, 3) += storage[1][1] * rhs(1, 2);
        tmp(1, 4) += storage[1][1] * rhs(2, 0);
        tmp(1, 3) += storage[1][1] * rhs(2, 1);
        tmp(1, 2) += storage[1][1] * rhs(2, 2);
        tmp(3, 0) += storage[1][2] * rhs(0, 0);
        tmp(3, 5) += storage[1][2] * rhs(0, 1);
        tmp(3, 4) += storage[1][2] * rhs(0, 2);
        tmp(3, 5) += storage[1][2] * rhs(1, 0);
        tmp(3, 1) += storage[1][2] * rhs(1, 1);
        tmp(3, 3) += storage[1][2] * rhs(1, 2);
        tmp(3, 4) += storage[1][2] * rhs(2, 0);
        tmp(3, 3) += storage[1][2] * rhs(2, 1);
        tmp(3, 2) += storage[1][2] * rhs(2, 2);
        tmp(4, 0) += storage[2][0] * rhs(0, 0);
        tmp(4, 5) += storage[2][0] * rhs(0, 1);
        tmp(4, 4) += storage[2][0] * rhs(0, 2);
        tmp(4, 5) += storage[2][0] * rhs(1, 0);
        tmp(4, 1) += storage[2][0] * rhs(1, 1);
        tmp(4, 3) += storage[2][0] * rhs(1, 2);
        tmp(4, 4) += storage[2][0] * rhs(2, 0);
        tmp(4, 3) += storage[2][0] * rhs(2, 1);
        tmp(4, 2) += storage[2][0] * rhs(2, 2);
        tmp(3, 0) += storage[2][1] * rhs(0, 0);
        tmp(3, 5) += storage[2][1] * rhs(0, 1);
        tmp(3, 4) += storage[2][1] * rhs(0, 2);
        tmp(3, 5) += storage[2][1] * rhs(1, 0);
        tmp(3, 1) += storage[2][1] * rhs(1, 1);
        tmp(3, 3) += storage[2][1] * rhs(1, 2);
        tmp(3, 4) += storage[2][1] * rhs(2, 0);
        tmp(3, 3) += storage[2][1] * rhs(2, 1);
        tmp(3, 2) += storage[2][1] * rhs(2, 2);
        tmp(2, 0) += storage[2][2] * rhs(0, 0);
        tmp(2, 5) += storage[2][2] * rhs(0, 1);
        tmp(2, 4) += storage[2][2] * rhs(0, 2);
        tmp(2, 5) += storage[2][2] * rhs(1, 0);
        tmp(2, 1) += storage[2][2] * rhs(1, 1);
        tmp(2, 3) += storage[2][2] * rhs(1, 2);
        tmp(2, 4) += storage[2][2] * rhs(2, 0);
        tmp(2, 3) += storage[2][2] * rhs(2, 1);
        tmp(2, 2) += storage[2][2] * rhs(2, 2);

        return tmp;
    }

    inline Vector6 Matrix3x3::VoigtVector() const
    {
        Vector6 tmp;
        tmp[0] = storage[0][0];
        tmp[1] = storage[1][1];
        tmp[2] = storage[2][2];
        tmp[3] = storage[1][2];
        tmp[4] = storage[0][2];
        tmp[5] = storage[0][1];
        return tmp;
    }

    class vStress : public Vector6
    {
        /*
         * Stress Voigt vector
         */
    public:
        vStress()
        {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
            storage[3] = 0.0;
            storage[4] = 0.0;
            storage[5] = 0.0;
        };
        /*vStress(vStress& rhs)
        {
            memmove(data(), rhs.data(), 6*sizeof(double));
        };*/
        vStress(const Vector6& rhs)
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            storage[3] = rhs[3];
            storage[4] = rhs[4];
            storage[5] = rhs[5];
        };
        vStress& operator=(const Vector6& rhs)
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            storage[3] = rhs[3];
            storage[4] = rhs[4];
            storage[5] = rhs[5];
            return *this;
        };
        /*vStress& operator=(vStress& rhs)
        {
            memmove(data(), rhs.data(), 6*sizeof(double));
            return *this;
        };*/
        vStress operator-(const vStress& rhs) const
        {
            vStress tmp;
            tmp[0] = storage[0] - rhs[0];
            tmp[1] = storage[1] - rhs[1];
            tmp[2] = storage[2] - rhs[2];
            tmp[3] = storage[3] - rhs[3];
            tmp[4] = storage[4] - rhs[4];
            tmp[5] = storage[5] - rhs[5];
            return tmp;
        };
        double norm(void) const  /// Frobenius norm
        {
            double tmp = storage[0] * storage[0]
                + storage[1] * storage[1]
                + storage[2] * storage[2]
                + storage[3] * storage[3] * 2.0
                + storage[4] * storage[4] * 2.0
                + storage[5] * storage[5] * 2.0;
            return sqrt(tmp);
        };

        double doublecontract(const vStress& Bstress) const  // "double-dot product"
        {
            double tmp =
                storage[0] * Bstress[0] +
                storage[1] * Bstress[1] +
                storage[2] * Bstress[2] +
                2.0 * storage[3] * Bstress[3] +
                2.0 * storage[4] * Bstress[4] +
                2.0 * storage[5] * Bstress[5];
            return tmp;
        };

        double doublecontract(const Vector6& symTensorV) const  // "double-dot product"
        {
            double tmp =
                storage[0] * symTensorV[0] +
                storage[1] * symTensorV[1] +
                storage[2] * symTensorV[2] +
                2.0 * storage[3] * symTensorV[3] +
                2.0 * storage[4] * symTensorV[4] +
                2.0 * storage[5] * symTensorV[5];
            return tmp;
        };

        inline double Pressure() const
        {
            return -(storage[0] + storage[1] + storage[2]) / 3.0;
        };

        inline double J1() const
        {
            return storage[0] + storage[1] + storage[2];
        };

        inline double Trace() const
        {
            return storage[0] + storage[1] + storage[2];
        };

        inline double Determinant() const
        {
            return storage[0] * storage[1] * storage[2] -
                storage[0] * storage[3] * storage[3] -
                storage[1] * storage[4] * storage[4] -
                storage[2] * storage[5] * storage[5] +
                storage[3] * storage[4] * storage[5] * 2;
        };

        Vector3 Invariants(void) const
        {
            Vector3 vec;
            vec[0] = Trace();
            vec[1] = 0.5 * (Trace() * Trace() -
                (storage[0] * storage[0] +
                    storage[1] * storage[1] +
                    storage[2] * storage[2] +
                    storage[3] * storage[3] * 2 +
                    storage[4] * storage[4] * 2 +
                    storage[5] * storage[5] * 2));
            vec[2] = Determinant();
            return vec;
        }

        double Mises(void) const
        {
            double vMises = 0.0;
            vMises = (storage[0] - storage[1]) * (storage[0] - storage[1]) +
                (storage[1] - storage[2]) * (storage[1] - storage[2]) +
                (storage[2] - storage[0]) * (storage[2] - storage[0]) +
                6 * (storage[3] * storage[3] + storage[4] * storage[4] + storage[5] * storage[5]);

            return sqrt(0.5 * vMises);
        }

        vStress get_rotated_matrix(const Matrix3x3& RotationMatrix) const
        {
            double In[3][3];
            double Out[3][3];

            In[0][0] = storage[0];
            In[0][1] = storage[5];
            In[0][2] = storage[4];
            In[1][0] = storage[5];
            In[1][1] = storage[1];
            In[1][2] = storage[3];
            In[2][0] = storage[4];
            In[2][1] = storage[3];
            In[2][2] = storage[2];

            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                {
                    Out[p][q] = 0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                        {
                            Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                        }
                }
            vStress _OUT;

            _OUT[0] = Out[0][0];
            _OUT[5] = Out[0][1];
            _OUT[4] = Out[0][2];
            _OUT[1] = Out[1][1];
            _OUT[3] = Out[1][2];
            _OUT[2] = Out[2][2];

            return _OUT;
        };
        vStress& do_rotate(const Matrix3x3& RotationMatrix)
        {
            double In[3][3];
            double Out[3][3];

            In[0][0] = storage[0];
            In[0][1] = storage[5];
            In[0][2] = storage[4];
            In[1][0] = storage[5];
            In[1][1] = storage[1];
            In[1][2] = storage[3];
            In[2][0] = storage[4];
            In[2][1] = storage[3];
            In[2][2] = storage[2];

            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                {
                    Out[p][q] = 0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                        {
                            Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                        }
                }

            storage[0] = Out[0][0];
            storage[5] = Out[0][1];
            storage[4] = Out[0][2];
            storage[1] = Out[1][1];
            storage[3] = Out[1][2];
            storage[2] = Out[2][2];

            return *this;
        };

        double get_tensor(const int i, const int j) const
        {
            return storage[(i == j) ? (i) : (6 - (i + j))];
        };

        Matrix3x3 tensor(void) const
        {
            Matrix3x3 tmp;
            tmp(0, 0) = storage[0];
            tmp(0, 1) = storage[5];
            tmp(0, 2) = storage[4];
            tmp(1, 0) = storage[5];
            tmp(1, 1) = storage[1];
            tmp(1, 2) = storage[3];
            tmp(2, 0) = storage[4];
            tmp(2, 1) = storage[3];
            tmp(2, 2) = storage[2];
            return tmp;
        };
    protected:
    private:
    };

    class vStrain : public Vector6
    {
        /*
         *  This is a special version of the six component Voigt vector:
         *  the off diagonal elements of the strain matrix are multiplied by 2(!)
         *  while transforming the 3x3 strain matrix to a 6-component vector,
         *  and have to be divided by two if they are transformed to a 3x3 matrix
         *  again. Therefore the transformation functions and some other functions
         *  have to be specified in different way then for the stress-type
         *  Voigt vector. The differences are specified in the functions.
         */
    public:
        vStrain()
        {
            storage[0] = 0.0;
            storage[1] = 0.0;
            storage[2] = 0.0;
            storage[3] = 0.0;
            storage[4] = 0.0;
            storage[5] = 0.0;
        };
        vStrain(const Vector6& rhs)
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            storage[3] = rhs[3];
            storage[4] = rhs[4];
            storage[5] = rhs[5];
        };
        vStrain& operator=(const Vector6& rhs)
        {
            storage[0] = rhs[0];
            storage[1] = rhs[1];
            storage[2] = rhs[2];
            storage[3] = rhs[3];
            storage[4] = rhs[4];
            storage[5] = rhs[5];
            return *this;
        };
        bool operator==(const Vector6& rhs)
        {
            for (int i = 0; i < 6; i++)
            {
                if (storage[i] != rhs[i]) { return false; };
            }
            return true;
        };
        vStrain operator-(const vStrain& rhs) const
        {
            vStrain tmp;
            tmp[0] = storage[0] - rhs[0];
            tmp[1] = storage[1] - rhs[1];
            tmp[2] = storage[2] - rhs[2];
            tmp[3] = storage[3] - rhs[3];
            tmp[4] = storage[4] - rhs[4];
            tmp[5] = storage[5] - rhs[5];
            return tmp;
        };
        double norm(void) const  /// Frobenius norm
        {
            /*
             * multiplication by 0.5 of the off diagonal elements before taking
             * square and doubling the result due to double appearance of the off
             * diagonal elements => 0.5^2 / 2 = 0.5
             */
            double tmp = storage[0] * storage[0]
                + storage[1] * storage[1]
                + storage[2] * storage[2]
                + storage[3] * storage[3] * 0.5
                + storage[4] * storage[4] * 0.5
                + storage[5] * storage[5] * 0.5;
            return sqrt(tmp);
        };

        double doublecontract(const vStrain& Bstrain) const // "double-dot product"
        {
            double tmp =
                storage[0] * Bstrain[0] +
                storage[1] * Bstrain[1] +
                storage[2] * Bstrain[2] +
                2.0 * (0.5 * storage[3]) * (0.5 * Bstrain[3]) +
                2.0 * (0.5 * storage[4]) * (0.5 * Bstrain[4]) +
                2.0 * (0.5 * storage[5]) * (0.5 * Bstrain[5]);
            return tmp;
        };

        vStrain get_rotated_matrix(const Matrix3x3& RotationMatrix)
        {
            /*
             * multiplication by 0.5 of the off diagonal elements during translating
             * to a 3x3 matrix and multiplying by two during translating back to
             * a Voigt strain vector
             */
            double In[3][3];
            double Out[3][3];

            In[0][0] = storage[0];
            In[0][1] = storage[5] * 0.5;
            In[0][2] = storage[4] * 0.5;
            In[1][0] = storage[5] * 0.5;
            In[1][1] = storage[1];
            In[1][2] = storage[3] * 0.5;
            In[2][0] = storage[4] * 0.5;
            In[2][1] = storage[3] * 0.5;
            In[2][2] = storage[2];

            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                {
                    Out[p][q] = 0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                        {
                            Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                        }
                }
            vStrain _OUT;

            _OUT[0] = Out[0][0];
            _OUT[5] = Out[0][1] * 2.0;
            _OUT[4] = Out[0][2] * 2.0;
            _OUT[1] = Out[1][1];
            _OUT[3] = Out[1][2] * 2.0;
            _OUT[2] = Out[2][2];

            return _OUT;
        };

        vStrain& do_rotate(const Matrix3x3& RotationMatrix)
        {
            /*
             * multiplication by 0.5 of the off diagonal elements during translating
             * to a 3x3 matrix and multiplying by two during translating back to
             * a Voigt strain vector
             */
            double In[3][3];
            double Out[3][3];

            In[0][0] = storage[0];
            In[0][1] = storage[5] * 0.5;
            In[0][2] = storage[4] * 0.5;
            In[1][0] = storage[5] * 0.5;
            In[1][1] = storage[1];
            In[1][2] = storage[3] * 0.5;
            In[2][0] = storage[4] * 0.5;
            In[2][1] = storage[3] * 0.5;
            In[2][2] = storage[2];

            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                {
                    Out[p][q] = 0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                        {
                            Out[p][q] += RotationMatrix(p, i) * In[i][j] * RotationMatrix(q, j);
                        }
                }

            storage[0] = Out[0][0];
            storage[5] = Out[0][1] * 2.0;
            storage[4] = Out[0][2] * 2.0;
            storage[1] = Out[1][1];
            storage[3] = Out[1][2] * 2.0;
            storage[2] = Out[2][2];

            return *this;
        };

        vStrain Ln() const;                                                     // Converts the strain to the logarithmic strain

        double get_tensor(const int i, const int j) const
        {
            /*
             * multiplication by 0.5 of the off diagonal elements during translating
             * to a 3x3 matrix
             */
            return (storage[(i == j) ? (i) : (6 - (i + j))] * ((i == j) ? (1.0) : (0.5)));
        };

        Matrix3x3 tensor(void) const
        {
            /*
             * multiplication by 0.5 of the off diagonal elements during translating
             * to a 3x3 matrix
             */
            Matrix3x3 tmp;
            tmp(0, 0) = storage[0];
            tmp(0, 1) = storage[5] * 0.5;
            tmp(0, 2) = storage[4] * 0.5;
            tmp(1, 0) = storage[5] * 0.5;
            tmp(1, 1) = storage[1];
            tmp(1, 2) = storage[3] * 0.5;
            tmp(2, 0) = storage[4] * 0.5;
            tmp(2, 1) = storage[3] * 0.5;
            tmp(2, 2) = storage[2];
            return tmp;
        };
    };

    inline vStrain Matrix3x3::VoigtStrain() const
    {
        vStrain tmp;
        tmp[0] = storage[0][0];
        tmp[1] = storage[1][1];
        tmp[2] = storage[2][2];
        tmp[3] = storage[1][2] * 2.0;
        tmp[4] = storage[0][2] * 2.0;
        tmp[5] = storage[0][1] * 2.0;
        return tmp;
    }

    inline vStress Matrix3x3::VoigtStress() const
    {
        vStress tmp;
        tmp[0] = storage[0][0];
        tmp[1] = storage[1][1];
        tmp[2] = storage[2][2];
        tmp[3] = storage[1][2];
        tmp[4] = storage[0][2];
        tmp[5] = storage[0][1];
        return tmp;
    }

    inline vStrain Matrix6x6::operator*(const vStress& rhs) const
    {
        vStrain tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    }

    inline vStress Matrix6x6::operator*(const vStrain& rhs) const
    {
        vStress tmp;
        tmp.set_to_zero();
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                tmp[i] += storage[i][j] * rhs[j];
            }
        return tmp;
    }

}