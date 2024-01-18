#include "rbm_methods.h"
#include "ReadComplex.h"
#include "MultiDimArrayAllocate.h"
#include "MultiDimArraySetToValue.h"
#include "MultiDimArrayPrint.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

void rbm_METHODS::read_rbmfam_NAMELIST()
{
    // Default values
    normcut = 1.0e-14;
    number_of_training_points = 0;
    number_of_operators = 0;
    fam_training_inputfile = "temp";
    rbm_outputfile = "temp";
    emulator_outputfile = "temp";
    VARIATION = 1;
    STRENGTH_EMULATOR_RUN = true;
    qrpa_omega_emulatorrun_start = 0.0;
    qrpa_omega_emulatorrun_step = 0.1;
    nmax_emulatorrun = 1;
    //
    for (int i = 0; i < 4; ++i)
    {
        mirror_points[i] = true;
    }

    // Open the file for reading
    std::ifstream inputFile("rbmfam_NAMELIST.dat");

    // Check if the file is open
    if (!inputFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
    }

    inputFile >> normcut >> number_of_training_points >> number_of_operators;
    inputFile >> fam_training_inputfile >> rbm_outputfile >> emulator_outputfile;
    inputFile >> VARIATION >> std::boolalpha >> STRENGTH_EMULATOR_RUN;

    //
    qrpa_omega_emulatorrun_start = readComplex(inputFile);
    qrpa_omega_emulatorrun_step = readComplex(inputFile);

    // Read the remaining values
    inputFile >> nmax_emulatorrun;
    for (int i = 0; i < 4; ++i)
    {
        inputFile >> std::boolalpha >> mirror_points[i];
    }

    // Close the file
    inputFile.close();
    //
    nop = number_of_operators;
    ntrain = number_of_training_points;
    // Dimension of the RBM
    nrbm = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (mirror_points[i])
        {
            nrbm += ntrain;
        }
    }
}

void rbm_METHODS::read_fam_training()
{
    std::ifstream training_input(fam_training_inputfile);
    if (!training_input.is_open())
    {
        std::cerr << "Error opening file: " << fam_training_inputfile << std::endl;
        return;
    }

    std::cout << "**** reading FAM training data ****" << std::endl;

    training_input >> nuv;
    // std::cout << "qq, nuv: " << nuv << std::endl;
    //  Allocations
    //  twoEqp = new std::complex<double>[nuv];
    //  twoEqp = allocate1dArraytest(nuv);
    //  allocate1dArray<std::complex<double>>(twoEqp, nuv);
    //  for (int i = 0; i < nuv; i++)
    //  {
    //  std::cout << twoEqp[i] << std::endl;
    //  twoEqp[i] = readComplex(training_input);
    //  }
    //
    allocate1dArray<std::complex<double>>(twoEqp, nuv);
    allocate2dArray<std::complex<double>>(F20, nuv, nop);
    allocate2dArray<std::complex<double>>(F02, nuv, nop);
    //
    allocate1dArray<std::complex<double>>(omegatrain, nrbm);
    allocate2dArray<std::complex<double>>(Xtrain, nuv, nrbm);
    allocate2dArray<std::complex<double>>(Ytrain, nuv, nrbm);
    allocate2dArray<std::complex<double>>(dH20, nuv, nrbm);
    allocate2dArray<std::complex<double>>(dH02, nuv, nrbm);

    // Set to zero.
    set1dArrayToValue<std::complex<double>>(twoEqp, nuv, 0.0);
    set2dArrayToValue<std::complex<double>>(F20, nuv, nop, 0.0);
    set2dArrayToValue<std::complex<double>>(F02, nuv, nop, 0.0);
    //
    set1dArrayToValue<std::complex<double>>(omegatrain, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(Xtrain, nuv, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(Ytrain, nuv, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(dH20, nuv, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(dH02, nuv, nrbm, 0.0);

    // std::cout << "show type of twoEqp" << twoEqp << std::endl;
    readComplex1d(training_input, twoEqp, nuv);
    // for (int i = 0; i < nuv; i++)
    //{
    //     std::cout << "In main: " << twoEqp[i] << std::endl;
    // }
    //  std::cout << readComplex(training_input) << std::endl;
    //  std::cout << readComplex(training_input) << std::endl;
    for (int i = 0; i < nop; i++)
    {
        readComplexToCol(training_input, F20, nuv, i);
        readComplexToCol(training_input, F02, nuv, i);
    }
    //
    for (int i = 0; i < ntrain; i++)
    {
        omegatrain[i] = readComplex(training_input);
        readComplexToCol(training_input, Xtrain, nuv, i);
        readComplexToCol(training_input, Ytrain, nuv, i);
        readComplexToCol(training_input, dH20, nuv, i);
        readComplexToCol(training_input, dH02, nuv, i);
    }
    // print2d<std::complex<double>>(dH02, nuv, ntrain);

    // ifstream automatically close the due to its destructor.
}

void rbm_METHODS::completeData()
{
    // ----------------------------------------------------------------------------------------------------
    // Complete data
    // ----------------------------------------------------------------------------------------------------
    // dH^{20}_{\mu\nu} + X_{\mu\nu}(omega)*(E_\mu + E_\nu)
    // dH^{02}_{\mu\nu} + Y_{\mu\nu}(omega)*(E_\mu + E_\nu)
    // need to define i_rbm which shows in main.
    //
    int i_rbm;
    for (int i = 1; i <= ntrain; ++i)
    {
        for (int k = 1; k <= nuv; ++k)
        {
            dH20[-1 + k][-1 + i] += Xtrain[-1 + k][-1 + i] * twoEqp[-1 + k];
            dH02[-1 + k][-1 + i] += Ytrain[-1 + k][-1 + i] * twoEqp[-1 + k];
        }
    }
    //
    i_rbm = 0;
    if (mirror_points[-1 + 1])
    { // training is assumed to be performed in the 1st quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
        }
    }
    //
    if (mirror_points[-1 + 2])
    { // includes points in the 2nd quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
            copy2dColTo2dCol(Ytrain, Xtrain, nuv, i_rbm, -1 + i);
            copy2dColTo2dCol(Xtrain, Ytrain, nuv, i_rbm, -1 + i);
            copy2dColTo2dCol(dH02, dH20, nuv, i_rbm, -1 + i);
            copy2dColTo2dCol(dH20, dH02, nuv, i_rbm, -1 + i);
            omegatrain[i_rbm] = -omegatrain[-1 + i] + iunit * std::imag(omegatrain[-1 + i]);
        }
    }
    //
    if (mirror_points[-1 + 4])
    { // includes points in the 2nd quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
            copy2dConjColTo2dCol(Xtrain, Xtrain, nuv, i_rbm, -1 + i);
            copy2dConjColTo2dCol(Ytrain, Ytrain, nuv, i_rbm, -1 + i);
            copy2dConjColTo2dCol(dH20, dH20, nuv, i_rbm, -1 + i);
            copy2dConjColTo2dCol(dH02, dH02, nuv, i_rbm, -1 + i);
            omegatrain[i_rbm] = std::conj(omegatrain[-1 + i]);
        }
    }
    //
    if (mirror_points[-1 + 3])
    { // includes points in the 2nd quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
            copy2dMConjColTo2dCol(Ytrain, Xtrain, nuv, i_rbm, -1 + i);
            copy2dMConjColTo2dCol(Xtrain, Ytrain, nuv, i_rbm, -1 + i);
            copy2dMConjColTo2dCol(dH20, dH20, nuv, i_rbm, -1 + i);
            copy2dMConjColTo2dCol(dH02, dH02, nuv, i_rbm, -1 + i);
            omegatrain[i_rbm] = -omegatrain[-1 + i];
        }
    }
    //
    if (i_rbm != nrbm)
    {
        std::cerr << "dimension error" << std::endl;
    }
}