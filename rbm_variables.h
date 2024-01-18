#pragma once
// #ifndef RBM_VARIABLES_H
// #define RBM_VARIABLES_H

#include <complex>
#include <vector>
#include <string>

using pr = double;
using Complex = std::complex<pr>;
// Constants
static constexpr int ipr = 1;
static constexpr pr pr_val = 1.0;
static const Complex iunit;

class rbm_VARIABLES
{
public:
    // Variables
    pr pi, Norm;
    int nop;
    int ntrain;
    int nrbm;
    //
    // FAM training data read from binary file
    int nuv;
    Complex *twoEqp, *omegatrain;
    Complex **Xtrain, **Ytrain, **dH20, **dH02, **F20, **F02;

    // Arrays that depend on ntrain
    Complex *SFtrain, *TFtrain, *norm_eigen;
    Complex **NormKernel, **HamiltonianKernel, **NormKernelHalf, **NormKernelRegularized,
        **unitmat, **NormKernelHalfInv, **u_norm, **u_norminv;
    pr *realtemp;

    bool SMALLESTFIRST;
    int colldim, ierr;
    int *collidx, *SortedOrder, *SortedOrder_Hcoll;

    // Arrays that depend on colldim
    Complex *RBMenergy, *sqrt_norm, *RBMstrength, *RBMstrength1, *RBMstrength2;
    Complex **H_coll, **g_coll, **g_collinv, **ugn, **ugn2, **X_QRPA_RBM, **Y_QRPA_RBM;
    pr *RBMstrengthfromXY;

    // Emulator
    Complex qrpa_omega, SFomega;
    pr *sumrule_emulator;

    // FAM Namelists
    pr normcut;
    int number_of_training_points, number_of_operators, nmax_emulatorrun, VARIATION;
    Complex qrpa_omega_emulatorrun_start, qrpa_omega_emulatorrun_step;
    bool STRENGTH_EMULATOR_RUN;
    bool mirror_points[4];

    std::string fam_training_inputfile, rbm_outputfile, emulator_outputfile;

    // struct FAM_RBM_INPUT
    // {
    // pr normcut;
    // int number_of_training_points;
    // int number_of_operators;
    // std::string fam_training_inputfile;
    // std::string rbm_outputfile;
    // std::string emulator_outputfile;
    // int VARIATION;
    // bool STRENGTH_EMULATOR_RUN;
    // Complex qrpa_omega_emulatorrun_start;
    // Complex qrpa_omega_emulatorrun_step;
    // int nmax_emulatorrun;
    // bool mirror_points[9];
    // };

    // FAM_RBM_INPUT fam_rbm_input;
};

// #endif // FAMRBM_DATA_H
