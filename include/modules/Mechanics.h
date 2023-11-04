#pragma once
#include "../baseTools/baseTools.h"
#include "Information.h"

#ifdef _WIN32
//windowsƽ̨ x86 or x68
#ifdef _WIN64
 //x64
#include "../../lib/x64/fftw3.h"
#pragma comment(lib,"lib/x64/libfftw3-3.lib")
#pragma comment(lib,"lib/x64/libfftw3f-3.lib")
#pragma comment(lib,"lib/x64/libfftw3l-3.lib")
#else
 //x86
#include "../../lib/x86/fftw3.h"
#pragma comment(lib,"lib/x86/libfftw3-3.lib")
#pragma comment(lib,"lib/x86/libfftw3f-3.lib")
#pragma comment(lib,"lib/x86/libfftw3l-3.lib")
#endif //_WIN64
#else
//unix
#include "../../lib/linux/fftw3.h"
//#ifdef __x86_64__
////x64
//#pragma comment(lib,"lib/x64/libfftw3-3.lib")
//#pragma comment(lib,"lib/x64/libfftw3f-3.lib")
//#pragma comment(lib,"lib/x64/libfftw3l-3.lib")
//#elif __i386__
////x86
//#include "../../lib/x86/fftw3.h"
//#pragma comment(lib,"lib/x86/libfftw3-3.lib")
//#pragma comment(lib,"lib/x86/libfftw3f-3.lib")
//#pragma comment(lib,"lib/x86/libfftw3l-3.lib")
//#endif
#endif //_WIN32


namespace pf {
    const std::complex< double > I(0.0, 1.0);                                       ///< sqrt(-1) declaration
    enum EvenExtensionDirection{EED_None, x_direction, y_direction, z_direction};

    class MechanicsProperty {
    public:
        MechanicsProperty() {};                                                      ///< Constructor (does nothing else)
        MechanicsProperty(Information& _inf);                                       //< Constructs and initializes class

        void Initialize(Information& _inf);                                        ///< Initializes the module, allocate the storage, assign internal variables

        vStrain AverageStrain;                                                     ///< Average strain calculated by the spectral solver.
        vStress AppliedStress;                                                     ///< Applied stress tensor
        vStrain AppliedStrain;                                                     ///< Applied strain tensor

        //Elasticity Tensors:
        Matrix6x6   AverageElasticConstants;                                       ///< Average elastic constants of the whole system
        //Matrix6x6   MAXElasticConstants;                                           ///< Maximum elastic constants of the system
        Matrix6x6   AverageCompliences;                                            ///< Inverse of the AverageElasticConstants
        Matrix6x6   MINCompliences;                                                ///< Inverse of the MAXElasticConstants
        vector<bool>     AvgStrainMask;                                            ///< Diagonal Components are 1 if Free Boundaries is selected for the corresponding direction
        vector<bool>     LoadStressMask;                                           ///< Diagonal Components are 1 if Applied Stress is selected for the corresponding direction
        vector<bool>     AppStrainMask;

        //tensor1_strain       EigenStrains;                                        ///< Reference eigenstrains for each phase field, Orientation considered if calculated with SetGrainProperties(Orientations&,..)
        tensor1_matrix6     ElasticConstants;                                    ///< Reference elastic constants for each phase field, grain orientation considered if calculated with SetGrainProperties(Orientations&,..)
        tensor1_matrix6     Compliences;                                         ///< Inverse of the ElasticConstants, grain orientation considered if calculated with SetGrainProperties(Orientations&,..)
        
        //tensor2_matrix6       Kappa;                                      ///< Coupling constants for the elasticity and diffusion for each phase field in rotated coordinates
        // Kappa[PhaseIndex][CompIndex]
        //tensor2_strain         Lambda;                                    ///< Coupling constants for the eigenstrain and diffusion for each phase field in rotated coordinates
        // Lambda[PhaseIndex][CompIndex]
        //tensor2_double          Cref;                                     /// Reference compositions for each phase and component
        // Cref[PhaseIndex][CompIndex]
        tensor1_matrix3       GrainRotation;                                ///< Symmetry variants and Coordinate rotation
        // GrainRotation[PhaseIndex]
        RotationGauge           rotation_gauge;
        EffectiveElasticConstantsModel effectiveElasticConstantsModel;
        MechanicsProperty& operator=(const MechanicsProperty& n) {
            AverageStrain = n.AverageStrain;
            AppliedStress = n.AppliedStress;
            AppliedStrain = n.AppliedStrain;
            AverageElasticConstants = n.AverageElasticConstants;
            //MAXElasticConstants = n.MAXElasticConstants;
            AverageCompliences = n.AverageCompliences;
            MINCompliences = n.MINCompliences;
            AvgStrainMask = n.AvgStrainMask;
            LoadStressMask = n.LoadStressMask;
            AppStrainMask = n.AppStrainMask;
            ElasticConstants = n.ElasticConstants;
            Compliences = n.Compliences;
            GrainRotation = n.GrainRotation;
            rotation_gauge = n.rotation_gauge;
            effectiveElasticConstantsModel = n.effectiveElasticConstantsModel;
        }
    };

	class Mechanics
	{
	public:
        Mechanics() {};
		Mechanics(FieldStorage_forPhaseNode& _phaseMesh, Information& _inf);
        ~Mechanics() {
            mechanicalField.free();
            inf = nullptr;
            phaseMesh = nullptr;
            for (int n = 0; n < 6; n++)
            {
                fftw_destroy_plan(ForwardPlanRHS[n]);

                delete[] rlRHSide[n];
                delete[] rcRHSide[n];
            }

            for (int n = 0; n < 3; n++)
            {
                fftw_destroy_plan(BackwardPlanU[n]);

                delete[] rlU[n];
                delete[] rcU[n];
                delete[] Q[n];
            }
            for (int n = 0; n < 9; n++)
            {
                fftw_destroy_plan(BackwardPlanDefGrad[n]);
                delete[] rcDefGrad[n];
                delete[] rlDefGrad[n];
            }
        }
		void Initialize(FieldStorage_forPhaseNode& _phaseMesh, Information& _inf);                                     ///< Constructor

        void ReInitialize();                                                   ///< Needed if the system has been remeshed
        void SetQ(void);                                                            ///< Sets the wave vectors

        void initGrainsProperties(int phase_index, int phase_property);
        ///////////////////////////////////////////////////////////////////////////////////////////
        // refer to OpenPhase Elasticity Khachaturyan model
        void init_mechanicalField();
        double SetEffectiveEigenStrains();
        void SetEffectiveElasticConstants();
        void SetMAXElasticConstants(vector<Matrix6x6> Cijs);

        //void cal_elastic_increment_of_oneNode(PhaseNode& node);
        ///////////////////////////////////////////////////////////////////////////////////////////

        // load stress
        int Solve(double StrainAccuracy, double StressAccuracy, int MAXIterations, bool is_dvStraindt_output = false);                                    ///< Solves elastic problem

        //return iterate rate
        void initVirtualEigenstrain();
        // apply strain
        double Solve2(double StrainAccuracy, int MAXIterations, double iterate_rate, bool is_dvStraindt_output = false, bool is_iterate_rate_auto_adjust = false);

        double*             rlU[3];												///< Displacements in real space
        FieldStorage_forMechanicNode mechanicalField;
    private:
        PhaseNode& get_meshNode_in_elasticField(int elas_x, int elas_y, int elas_z);
        void CalculateRHS(Matrix6x6 Cij);
        void CalculateRHS2(Matrix6x6 Cij);
        void ExecuteForwardFFT();
        void CalculateFourierSolution(Matrix6x6 Cij);
        void ExecuteBackwardFFT();
        void evaluate_virtualEigenstrain(Matrix6x6 Cij, Matrix6x6 Sij, double& MAXvStrainDifference, double iterate_rate);
        void SetElasticProperties1(double& MAXStrainDifference, vStress& AverageStress);
        void SetElasticProperties2(double& MAXStrainDifference, double& MAXStressDifference, int icount);
        void SetElasticBoundaryConditions(vStress TargetStress);
        void SetElasticBoundaryConditions2(Matrix6x6 Cij);
        void assignment_to_phaseField();


        Information* inf;
        FieldStorage_forPhaseNode* phaseMesh;
        MechanicsProperty MP;


        int Nz2;                                                               ///< Half of the system size along Z direction

        int rlSIZE;                                                            ///< System size (real space): Nx*Ny*Nz
        int rcSIZE;                                                            ///< System size (Fourier space): Nx*Ny*(Nz/2+1)

        double Norm;                                                                ///< Normalization coefficient: 1.0/SIZE
        double DPi_Nx;                                                              ///< 2.0*Pi/Nx constant
        double DPi_Ny;                                                              ///< 2.0*Pi/Ny constant
        double DPi_Nz;                                                              ///< 2.0*Pi/Nz constant

        Matrix6x6 C0;
        Matrix6x6 C0inverse;

        Matrix6x6 MAX_ElasticConstants;

        double*                 rlRHSide[6];                                         /// Right hand side in real space
        std::complex<double>*   rcRHSide[6];                                         /// Right hand side in Fourier space

        double*                 Q[3];                                                /// Wave vectors
        std::complex<double>*   rcU[3];												///< Displacements in reciprocal space

        double*                 rlDefGrad[9];                                        /// Deformation gradient entries in real space
        std::complex<double>*   rcDefGrad[9];                                        /// Deformation gradient entries in Fourier space

        fftw_plan               ForwardPlanRHS[6];                                               /// Forward FFT plans for RHSide
        fftw_plan               BackwardPlanDefGrad[9];                                          /// Backward FFT plans for deformation gradients
        fftw_plan               BackwardPlanU[3];                                                ///< Backward FFT plans for displacements

	};
}

