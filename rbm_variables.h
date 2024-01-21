#pragma once
// #ifndef RBM_VARIABLES_H
// #define RBM_VARIABLES_H

#include <complex>
#include <vector>
#include <string>

using pr = double;
using Complex = std::complex<pr>;
// Constants
// static constexpr int ipr = 1;
// static constexpr pr pr_val = 1.0;
static const Complex iunit = Complex(0.0, 1.0);
static double pi = std::acos(-1.0);
//
class rbm_VARIABLES
{
public:
    int outp = 20, emu = 21;
    Complex **tempmat = nullptr, **tempmat2 = nullptr, **tempmat3 = nullptr;
    // Variables
    pr pi, Norm;
    int nop;
    int ntrain;
    int nrbm;
    //
    // FAM training data read from binary file
    int nuv;
    Complex *twoEqp = nullptr, *omegatrain = nullptr;
    Complex **Xtrain = nullptr, **Ytrain = nullptr, **dH20 = nullptr, **dH02 = nullptr, **F20 = nullptr, **F02 = nullptr;

    // Arrays that depend on ntrain
    Complex *SFtrain = nullptr, *TFtrain = nullptr, *norm_eigen = nullptr;
    Complex **NormKernel = nullptr, **HamiltonianKernel = nullptr, **NormKernelHalf = nullptr, **NormKernelRegularized = nullptr,
            **unitmat = nullptr, **NormKernelHalfInv = nullptr, **u_norm = nullptr, **u_norminv = nullptr;
    pr *realtemp = nullptr;

    bool SMALLESTFIRST;
    int colldim, ierr;
    int *collidx = nullptr, *SortedOrder = nullptr, *SortedOrder_Hcoll = nullptr;

    // Arrays that depend on colldim
    Complex *RBMenergy = nullptr, *sqrt_norm = nullptr, *RBMstrength = nullptr, *RBMstrength1 = nullptr, *RBMstrength2 = nullptr;
    Complex **H_coll = nullptr, **g_coll = nullptr, **g_collinv = nullptr, **ugn = nullptr, **ugn2 = nullptr, **X_QRPA_RBM = nullptr, **Y_QRPA_RBM = nullptr;
    pr **RBMstrengthfromXY = nullptr;

    // Emulator
    Complex qrpa_omega, SFomega;
    pr *sumrule_emulator = nullptr;

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
