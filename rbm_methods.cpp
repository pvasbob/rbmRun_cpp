#include "rbm_methods.h"
#include "ReadComplex.h"
#include "MultiDimArrayAllocate.h"
#include "MultiDimArraySetToValue.h"
#include "MultiDimArrayPrint.h"
#include "MultiDimArraySort.h"
#include "DotProduct.h"

#include "MsgToScreen.h"
#include "DiagMat.h"

#include <iostream>
#include <fstream>
#include <iomanip>

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
    // print2d<std::complex<double>>(dH02, nuv, nrbm);

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
    // print2d(dH02, nuv, nrbm);
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
    // std::cout << "before 1st quad." << std::endl;
    // print2d(dH02, nuv, nrbm);
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
    // std::cout << "before 2st quad." << std::endl;
    // print2d(dH02, nuv, nrbm);
    //
    if (mirror_points[-1 + 2])
    { // includes points in the 2nd quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
            // std::cout << "i_rbm: " << i_rbm << std::endl;
            copy2dColTo2dCol(Ytrain, Xtrain, nuv, -1 + i, -1 + i_rbm);
            copy2dColTo2dCol(Xtrain, Ytrain, nuv, -1 + i, -1 + i_rbm);
            copy2dColTo2dCol(dH02, dH20, nuv, -1 + i, -1 + i_rbm);
            copy2dColTo2dCol(dH20, dH02, nuv, -1 + i, -1 + i_rbm);
            //
            // std::cout << "check: " << std::imag(omegatrain[-1 + i]) << "iunit: " << iunit << "product: " << iunit * std::imag(omegatrain[-1 + i]) << std::endl;
            //
            omegatrain[-1 + i_rbm] = -std::real(omegatrain[-1 + i]) + iunit * std::imag(omegatrain[-1 + i]);
        }
    }
    //
    // std::cout << "before 4st quad." << std::endl;
    // print2d(dH02, nuv, nrbm);
    //
    if (mirror_points[-1 + 4])
    { // includes points in the 2nd quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
            copy2dConjColTo2dCol(Xtrain, Xtrain, nuv, -1 + i, -1 + i_rbm);
            copy2dConjColTo2dCol(Ytrain, Ytrain, nuv, -1 + i, -1 + i_rbm);
            copy2dConjColTo2dCol(dH20, dH20, nuv, -1 + i, -1 + i_rbm);
            copy2dConjColTo2dCol(dH02, dH02, nuv, -1 + i, -1 + i_rbm);
            omegatrain[-1 + i_rbm] = std::conj(omegatrain[-1 + i]);
        }
    }
    //
    if (mirror_points[-1 + 3])
    { // includes points in the 2nd quadrant
        for (int i = 1; i <= ntrain; ++i)
        {
            ++i_rbm;
            copy2dMConjColTo2dCol(Ytrain, Xtrain, nuv, -1 + i, -1 + i_rbm);
            copy2dMConjColTo2dCol(Xtrain, Ytrain, nuv, -1 + i, -1 + i_rbm);
            copy2dMConjColTo2dCol(dH02, dH20, nuv, -1 + i, -1 + i_rbm);
            copy2dMConjColTo2dCol(dH20, dH02, nuv, -1 + i, -1 + i_rbm);
            omegatrain[-1 + i_rbm] = -omegatrain[-1 + i];
        }
    }
    // std::cout << "after 3st quad." << std::endl;
    // print2d(Ytrain, nuv, nrbm);
    //
    if (i_rbm != nrbm)
    {
        std::cerr << "dimension error" << std::endl;
    }
    //
    // print2d(dH02, nuv, nrbm);
}
//
void rbm_METHODS::strengthAtTraining()
{
    allocate1dArray<std::complex<double>>(SFtrain, nrbm);
    allocate1dArray<std::complex<double>>(TFtrain, nrbm);
    // no need to set to zero for SFtrain and TFtrain cause they will immediately assigned by doct product below.
    for (int i = 1; i <= nrbm; i++)
    {
        SFtrain[-1 + i] = dot2dCol2dCol(F20, Xtrain, -1 + 1, -1 + i, nuv) + dot2dCol2dCol(F02, Ytrain, -1 + 1, -1 + i, nuv);
        TFtrain[-1 + i] = dot2dConjCol2dCol(F20, Xtrain, -1 + 1, -1 + i, nuv) + dot2dConjCol2dCol(F02, Ytrain, -1 + 1, -1 + i, nuv);
    }
    // print1d<std::complex<double>>(SFtrain, nrbm);

    std::cout << "Strength functions at the training energies (SF, TF)" << std::endl;
    std::cout << "i    SF(Re, Im),   TF(Re, Im) omega(Re, Im)" << std::endl;

    for (int i = 1; i <= nrbm; ++i)
    {
        std::cout << std::setw(5) << -1 + i
                  << std::scientific << std::setprecision(10)
                  << std::setw(20) << SFtrain[-1 + i]
                  << std::setw(20) << TFtrain[-1 + i]
                  << std::setw(20) << omegatrain[-1 + i] << std::endl;
    }
}

void rbm_METHODS::kernelCalculation()
{
    //===================================
    // Allocation of matrices
    //===================================
    allocate2dArray<std::complex<double>>(NormKernel, nrbm, nrbm);
    allocate2dArray<std::complex<double>>(HamiltonianKernel, nrbm, nrbm);
    allocate2dArray<std::complex<double>>(u_norm, nrbm, nrbm);
    allocate2dArray<std::complex<double>>(u_norminv, nrbm, nrbm);
    allocate1dArray<std::complex<double>>(norm_eigen, nrbm);
    allocate1dArray<int>(collidx, nrbm);
    allocate1dArray<std::complex<double>>(sqrt_norm, nrbm);
    //
    allocate2dArray<std::complex<double>>(NormKernelHalf, nrbm, nrbm);
    allocate2dArray<std::complex<double>>(NormKernelHalfInv, nrbm, nrbm);
    //
    allocate2dArray<std ::complex<double>>(tempmat, nrbm, nrbm);
    allocate2dArray<std ::complex<double>>(tempmat2, nrbm, nrbm);
    allocate2dArray<std ::complex<double>>(tempmat3, nrbm, nrbm);
    allocate2dArray<std ::complex<double>>(NormKernelRegularized, nrbm, nrbm);
    allocate2dArray<std ::complex<double>>(unitmat, nrbm, nrbm);
    //
    allocate1dArray<double>(realtemp, nrbm);
    //
    //
    //====================================
    // set matricies to zero.
    //====================================
    set2dArrayToValue<std::complex<double>>(NormKernel, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(HamiltonianKernel, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(u_norm, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(u_norminv, nrbm, nrbm, 0.0);
    set1dArrayToValue<std::complex<double>>(norm_eigen, nrbm, 0.0);
    set1dArrayToValue<int>(collidx, nrbm, 0.0);
    set1dArrayToValue<std::complex<double>>(sqrt_norm, nrbm, 0.0);
    //
    set2dArrayToValue<std::complex<double>>(NormKernelHalf, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(NormKernelHalfInv, nrbm, nrbm, 0.0);
    //
    set2dArrayToValue<std::complex<double>>(tempmat, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(tempmat2, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(tempmat3, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(NormKernelRegularized, nrbm, nrbm, 0.0);
    set2dArrayToValue<std::complex<double>>(unitmat, nrbm, nrbm, 0.0);
    //
    set1dArrayToValue<double>(realtemp, nrbm, 0.0);
    //
    //
    //
    for (int i = 1; i <= nrbm; i++)
    {
        unitmat[-1 + i][-1 + i] = 1.0;
    }
    //=====================================
    // Kernel Calculation
    //=====================================
    //
    //
    // std::cout << "Xtrain before VARIATION: " << std::endl;
    // print2d<std::complex<double>>(Xtrain, nuv, nrbm);
    // std::cout << "Ytrain before VARIATION: " << std::endl;
    // print2d<std::complex<double>>(Ytrain, nuv, nrbm);
    //
    if (VARIATION == 1)
    {
        for (int i = 1; i <= nrbm; i++)
        {
            for (int j = 1; j <= nrbm; j++)
            {
                // for (int k = 1; k <= nuv; k++)
                //{
                //     std::cout << "Xtrain, -1+i, -1+j, -1+k: " << -1 + i << " " << -1 + j << " " << -1 + k << " " << Xtrain[-1 + k][-1 + i] << " " << Xtrain[-1 + k][-1 + j] << std::endl;
                // }
                //
                // for (int k = 1; k <= nuv; k++)
                //{
                //    std::cout << "Ytrain, -1+i, -1+j, -1+k: " << -1 + i << " " << -1 + j << " " << -1 + k << " " << Ytrain[-1 + k][-1 + i] << " " << Ytrain[-1 + k][-1 + j] << std::endl;
                //}
                // Norm Kernel
                // std::cout << "dot product Xtran Xtrain, i, j: " << i << " " << j << " " << dot2dCol2dCol(Xtrain, Xtrain, -1 + i, -1 + j, nuv) << std::endl;
                NormKernel[-1 + i][-1 + j] = dot2dConjCol2dCol(Xtrain, Xtrain, -1 + i, -1 + j, nuv) - dot2dConjCol2dCol(Ytrain, Ytrain, -1 + i, -1 + j, nuv);
                // std::cout << "Inside VARIATION: " << -1 + i << " " << -1 + j << " " << NormKernel[-1 + i][-1 + j] << std::endl;
                // Hamiltonian Kernel
                HamiltonianKernel[-1 + i][-1 + j] = dot2dConjCol2dCol(Xtrain, dH20, -1 + i, -1 + j, nuv) + dot2dConjCol2dCol(Ytrain, dH02, -1 + i, -1 + j, nuv);
            }
        }
    }
    else if (VARIATION == 2)
    {
        for (int j = 1; j <= nrbm; j++)
        {
            for (int i = 1; i <= nrbm; i++)
            {
                // Norm Kernel
                NormKernel[-1 + i][-1 + j] = dot2dConjCol2dCol(Xtrain, Xtrain, -1 + i, -1 + j, nuv) - dot2dConjCol2dCol(Ytrain, Ytrain, -1 + i, -1 + j, nuv);
                // Hamiltonian Kernel
                HamiltonianKernel[-1 + i][-1 + j] = dot2dConjCol2dCol(Xtrain, dH20, -1 + i, -1 + j, nuv) + dot2dConjCol2dCol(Ytrain, dH02, -1 + i, -1 + j, nuv);
            }
        }
    }
    else
    {
        std::cerr << "VARIATION should be 1 or 2" << std::endl;
        exit(1);
    }
    //
    // std::cout << "NormKernel: " << std::endl;
    // print2d<std::complex<double>>(NormKernel, nrbm, nrbm);
    // std::cout << "HamiltonianKernel: " << std::endl;
    // print2d<std::complex<double>>(HamiltonianKernel, nrbm, nrbm);
    //========================================================
    // Remove tiny numerical error from NormKernel
    //========================================================
    for (int i = 1; i <= nrbm; ++i)
    {
        for (int j = 1; j <= nrbm; ++j)
        {
            if (std::abs(std::real(NormKernel[-1 + i][-1 + j])) < 1.0e-15)
            {
                // NormKernel(i,j) = iunit * Aimag(NormKernel(i,j))
                NormKernel[-1 + i][-1 + j] = iunit * std::imag(NormKernel[-1 + i][-1 + j]);
            }
            if (std::abs(std::imag(NormKernel[-1 + i][-1 + j])) < 1.0e-15)
            {
                // NormKernel(i,j) = Dble(NormKernel(i,j))
                NormKernel[-1 + i][-1 + j] = std::real(NormKernel[-1 + i][-1 + j]);
            }
        }
    }
    //
    // std::cout << "NormKernel after remove tiny errors: " << std::endl;
    // print2d<std::complex<double>>(NormKernel, nrbm, nrbm);
}

void rbm_METHODS::diagNormKernel()
{
    diagGenComplexMat(nrbm, NormKernel, u_norm, norm_eigen, u_norminv, ierr);
    //
    if (ierr != 0)
    {
        std::cout << "ZGEEV ERROR: INFO = " << ierr << std::endl;
        std::cerr << "ZGEEV ERROR: INFO = " << ierr << std::endl;
        throw std::runtime_error("norm kernel diagonalization failed");
    }
}

void rbm_METHODS::sortNormEigen()
{
    SMALLESTFIRST = false;
    if (SortedOrder != nullptr)
        deallocate1dArray<int>(SortedOrder);
    allocate1dArray<int>(SortedOrder, nrbm);
    copy1dRealTo1d(norm_eigen, realtemp, nrbm, '+');
    //
    msgToScreen("realtemp:");
    print1d(realtemp, nrbm);
    //
    sort2(nrbm, realtemp, SortedOrder, SMALLESTFIRST);
    msgToScreen("SortedOrder:");
    print1d(SortedOrder, nrbm);

    if (VARIATION == 1)
    {
        // Norm kernel should be Hermitian. Remove small imaginary part
        for (int i = 1; i <= nrbm; ++i)
        {
            if (std::abs(std::imag(norm_eigen[-1 + i])) < 1.0e-15)
            {
                norm_eigen[-1 + i] = std::real(norm_eigen[-1 + i]);
            }
        }
    }

    std::cout << "Norm eigenvalues (Re, Im)" << std::endl;
    for (int i = 1; i <= nrbm; ++i)
    {
        std::cout << -1 + i + 1 << " " << std::real(norm_eigen[-1 + SortedOrder[-1 + i]]) << " "
                  << std::imag(norm_eigen[-1 + SortedOrder[-1 + i]]) << std::endl;
    }

    std::cout << "----------------------------" << std::endl;
    std::cout << "normcut = " << normcut << std::endl;
    //
    deallocate1dArray<double>(realtemp);
}

void rbm_METHODS::normEigenCutoff()
{
    int j = 0;

    std::cout << "   i      Sqrt(norm) (Re, Im) " << std::endl;
    for (int i = 1; i <= nrbm; ++i)
    {
        if (std::abs(norm_eigen[-1 + SortedOrder[-1 + i]]) < normcut)
        {
            continue; // Equivalent to Fortran's "Cycle"
        }
        else
        {
            j++;
            collidx[-1 + j] = SortedOrder[-1 + i];
        }

        sqrt_norm[-1 + j] = std::sqrt(norm_eigen[-1 + SortedOrder[-1 + i]]);
        if (VARIATION == 1 && std::real(norm_eigen[-1 + SortedOrder[-1 + i]]) < 0.0 &&
            std::abs(std::real(sqrt_norm[-1 + j])) < 1.0e-15)
        {
            sqrt_norm[-1 + j] = std::imag(sqrt_norm[-1 + j]) * iunit; // 1.0i represents the imaginary unit
        }

        std::cout << std::setw(5) << i + 1 << std::setw(20) << std::setprecision(10) << std::real(sqrt_norm[-1 + j])
                  << std::setw(20) << std::setprecision(10) << std::imag(sqrt_norm[-1 + j]) << std::endl;
    }

    colldim = j;
    std::cout << "colldim = " << colldim << std::endl;
    if (colldim == 0)
    {
        throw std::runtime_error("colldim = 0");
    }
}

void rbm_METHODS::normKernelExcludingSmallEigen()
{
    set2dArrayToValue<std::complex<double>>(NormKernelRegularized, nrbm, nrbm, 0.0);
    //
    std::cout << "nrbm, colldim: " << nrbm << " " << colldim << std::endl;
    for (int j = 1; j <= nrbm; ++j)
    {
        for (int i = 1; i <= nrbm; ++i)
        {
            for (int k = 1; k <= colldim; ++k)
            {
                NormKernelRegularized[-1 + i][-1 + j] += u_norm[-1 + i][-1 + collidx[-1 + k]] * norm_eigen[-1 + collidx[-1 + k]] * u_norminv[-1 + collidx[-1 + k]][-1 + j];
            }
        }
    }

    // check
    set2dArrayToValue<std::complex<double>>(tempmat, nrbm, nrbm, 0.0);
    for (int j = 1; j <= nrbm; ++j)
    {
        for (int i = 1; i <= nrbm; ++i)
        {
            for (int k = 1; k <= nrbm; ++k)
            {
                tempmat[-1 + i][-1 + j] += u_norm[-1 + i][-1 + k] * norm_eigen[-1 + k] * u_norminv[-1 + k][-1 + j];
            }
        }
    }
    //
    msgToScreen("HNormKernelRegularized:");
    print2d(NormKernelRegularized, nrbm, nrbm);
    msgToScreen("Htempmat:");
    print2d(tempmat, nrbm, nrbm);
    //
    std::cout << "norm kernel diagonalization check" << std::endl;
    // NOT important, implement later.
    // Print *, "norm kernel diagonalization check"
    // Print *, "max diff N - u n u^{-1} : ", maxval(Abs( NormKernel(:,:) - tempmat(:,:)))
    // Print *, "max diff N - Nreg       : ", maxval(Abs( NormKernel(:,:) - NormKernelRegularized(:,:)))
    // Print *, "norm kernel diagonalization check completed"
}

void rbm_METHODS::normKernelHalf()
{
    allocate2dArray<std::complex<double>>(NormKernelHalf, nrbm, nrbm);
    allocate2dArray<std::complex<double>>(NormKernelHalfInv, nrbm, nrbm);
    //
    for (int i = 1; i <= nrbm; ++i)
    {
        for (int j = 1; j <= nrbm; ++j)
        {
            for (int k = 1; k <= colldim; ++k)
            {
                NormKernelHalf[-1 + i][-1 + j] += u_norm[-1 + i][-1 + collidx[-1 + k]] * sqrt_norm[-1 + k] * u_norminv[-1 + collidx[-1 + k]][-1 + j];
                NormKernelHalfInv[-1 + i][-1 + j] += u_norm[-1 + i][-1 + collidx[-1 + k]] * std::pow(sqrt_norm[-1 + k], -1.0) * u_norminv[-1 + collidx[-1 + k]][-1 + j];
            }
        }
    }
    //
    msgToScreen("NormKernelHalf before tempmat:");
    print2d(NormKernelHalf, nrbm, nrbm);
    msgToScreen("NormKernelHalfInv before tempmat:");
    print2d(NormKernelHalfInv, nrbm, nrbm);
    //
    //
    //
    mult2dAnd2d<std::complex<double>>(NormKernelHalf, NormKernelHalf, tempmat, nrbm, nrbm, nrbm);
    //
    msgToScreen("NormKernelHalf: tempmat:");
    print2d(tempmat, nrbm, nrbm);
    //
    //  Below not important, implement later.
    //! check
    // Print *, "squareroot of norm kernel check"
    // Print *, "max diff N - N^{1/2} N^{1/2}: ", maxval(Abs( NormKernel(:,:) - tempmat(:,:)))
    // Print *, "squareroot of norm kernel check completed"
    deallocate2dArray(tempmat, nrbm);
}

void rbm_METHODS::Hcoll()
{
    msgToScreen("Allocation involing colldim:");
    std::cout << colldim << std::endl;
    //
    allocate2dArray(H_coll, colldim, colldim);
    allocate2dArray(g_coll, colldim, colldim);
    allocate1dArray(RBMenergy, colldim);
    allocate1dArray(SortedOrder_Hcoll, colldim);
    allocate2dArray(g_collinv, colldim, colldim);
    allocate1dArray(RBMstrength, colldim);
    allocate2dArray(ugn, nrbm, colldim);
    allocate2dArray(ugn2, colldim, nrbm);
    allocate1dArray(RBMstrength1, colldim);
    allocate1dArray(RBMstrength2, colldim);
    allocate2dArray(X_QRPA_RBM, nuv, colldim);
    allocate2dArray(Y_QRPA_RBM, nuv, colldim);
    allocate2dArray(RBMstrengthfromXY, colldim, nop);
    //
    // Memory allocation.
    set2dArrayToValue<std::complex<double>>(H_coll, colldim, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(g_coll, colldim, colldim, 0.0);
    set1dArrayToValue<std::complex<double>>(RBMenergy, colldim, 0.0);
    set1dArrayToValue<int>(SortedOrder_Hcoll, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(g_collinv, colldim, colldim, 0.0);
    set1dArrayToValue<std::complex<double>>(RBMstrength, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(ugn, nrbm, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(ugn2, colldim, nrbm, 0.0);
    set1dArrayToValue<std::complex<double>>(RBMstrength1, colldim, 0.0);
    set1dArrayToValue<std::complex<double>>(RBMstrength2, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(X_QRPA_RBM, nuv, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(Y_QRPA_RBM, nuv, colldim, 0.0);
    set2dArrayToValue<double>(RBMstrengthfromXY, colldim, nop, 0.0);
    //
    // Hcoll
    set2dArrayToValue<std::complex<double>>(H_coll, colldim, colldim, 0.0);
    for (int j = 1; j <= colldim; ++j)
    {
        for (int i = 1; i <= colldim; ++i)
        {
            for (int l = 1; l <= nrbm; ++l)
            {
                for (int k = 1; k <= nrbm; ++k)
                {
                    H_coll[-1 + i][-1 + j] += u_norminv[-1 + collidx[-1 + i]][-1 + k] * HamiltonianKernel[-1 + k][-1 + l] * (u_norm[-1 + l][-1 + collidx[-1 + j]]) / (sqrt_norm[-1 + i] * sqrt_norm[-1 + j]);
                }
            }
        }
    }
    //
    // msgToScreen("H_coll:");
    // print2d(H_coll, colldim, colldim);
}

void rbm_METHODS::diagCollectiveHamiltonian()
{
    std::cout << "**** Diagonalizing collective Hamiltonian ****" << std::endl;

    diagGenComplexMat(colldim, H_coll, g_coll, RBMenergy, g_collinv, ierr);
    //
    msgToScreen("g_coll:");
    print2d(g_coll, colldim, colldim);
    //
    msgToScreen("RBMenergy:");
    print1d(RBMenergy, colldim);
    //
    msgToScreen("g_collinv:");
    print2d(g_collinv, colldim, colldim);
    //
    if (ierr != 0)
    {
        std::cout << "Hcoll diagonalization ERROR: INFO = " << ierr << std::endl;
        std::cout << "Hcoll diagonalization ERROR: INFO = " << ierr << std::endl;
        throw std::runtime_error("Hcoll diagonalization failed");
    }
    //
    SMALLESTFIRST = true;
    if (realtemp != nullptr)
    {
        deallocate1dArray(realtemp);
    }
    allocate1dArray(realtemp, colldim);
    //
    for (int i = 0; i < colldim; i++)
    {
        realtemp[i] = std::real(RBMenergy[i]);
    }
    //
    sort2(colldim, realtemp, SortedOrder_Hcoll, SMALLESTFIRST);
    //
    msgToScreen("SortedOrder_Hcoll:");
    print1d(SortedOrder_Hcoll, colldim);
    //
    std::cout << "QROA energies" << std::endl;
    std::cout << " i    ReE      ImE" << std::endl;
    //
    for (int k = 1; k <= colldim; ++k)
    {
        std::cout << std::setw(5) << -1 + k + 1 << std::setw(20) << std::setprecision(10) << std::real(RBMenergy[-1 + SortedOrder_Hcoll[-1 + k]])
                  << std::setw(20) << std::setprecision(10) << std::imag(RBMenergy[-1 + SortedOrder_Hcoll[-1 + k]]) << std::endl;
    }
    //
    // The deallocate below can't work, for we dont know ( or very hard to track) its dimension.
    // If(Allocated(tempmat)) Deallocate(tempmat)
    allocate2dArray(tempmat, colldim, colldim);
    for (int i = 1; i <= colldim; ++i)
    {
        for (int j = 1; j <= colldim; ++j)
        {
            for (int k = 1; k <= colldim; ++k)
            {
                tempmat[-1 + i][-1 + j] += g_coll[-1 + i][-1 + SortedOrder_Hcoll[-1 + k]] * RBMenergy[-1 + SortedOrder_Hcoll[-1 + k]] * g_collinv[-1 + SortedOrder_Hcoll[-1 + k]][-1 + j];
            }
        }
    }
    //
    msgToScreen("tempmat colldim: ");
    print2d(tempmat, colldim, colldim);
}

void rbm_METHODS::calculateStrength()
{
    set2dArrayToValue<std::complex<double>>(ugn, nrbm, colldim, 0.0);
    set2dArrayToValue<std::complex<double>>(ugn2, colldim, nrbm, 0.0);
    //
    for (int i = 1; i <= nrbm; ++i)
    {
        for (int l = 1; l <= colldim; ++l)
        {
            for (int k = 1; k <= colldim; ++k)
            {
                ugn[-1 + i][-1 + l] += u_norm[-1 + i][-1 + collidx[-1 + k]] * g_coll[-1 + k][-1 + l] / sqrt_norm[-1 + k];
                ugn2[-1 + l][-1 + i] += g_collinv[-1 + l][-1 + k] * u_norminv[-1 + collidx[-1 + k]][-1 + i] / sqrt_norm[-1 + k];
            }
        }
    }
    //
    set1dArrayToValue<std::complex<double>>(RBMstrength1, colldim, 0.0);
    //
    mult1dAnd2d(SFtrain, ugn, RBMstrength1, nrbm, colldim);
    //
    set1dArrayToValue<std::complex<double>>(RBMstrength2, colldim, 0.0);
    //
    if (VARIATION == 1)
    {
        mult2dAndConj1d(ugn2, SFtrain, RBMstrength2, colldim, nrbm);
    }
    else
    {
        mult2dAnd1d(ugn2, TFtrain, RBMstrength2, colldim, nrbm);
    }
    //
    for (int i = 1; i <= colldim; ++i)
    {
        RBMstrength[-1 + i] = RBMstrength1[-1 + i] * RBMstrength2[-1 + i];
    }
    //
    msgToScreen("RBMstrength:");
    print1d(RBMstrength, colldim);
}

void rbm_METHODS::calculateQRPAXY()
{
    msgToScreen("CalXtrain");
    print2d(Xtrain, nuv, nrbm);
    msgToScreen("CalYtrain");
    print2d(Ytrain, nuv, nrbm);
    msgToScreen("Calugn");
    print2d(ugn, nrbm, colldim);
    //
    mult2dAnd2d(Xtrain, ugn, X_QRPA_RBM, nuv, nrbm, colldim);
    mult2dAnd2d(Ytrain, ugn, Y_QRPA_RBM, nuv, nrbm, colldim);
    //
    msgToScreen("X_QRPA_RBM");
    print2d(X_QRPA_RBM, nuv, colldim);
    //
    msgToScreen("Y_QRPA_RBM");
    print2d(Y_QRPA_RBM, nuv, colldim);
    //
    for (int i = 1; i <= colldim; i++)
    {
        Norm = std::real(dot2dConjCol2dCol(X_QRPA_RBM, X_QRPA_RBM, -1 + i, -1 + i, nuv) - dot2dConjCol2dCol(Y_QRPA_RBM, Y_QRPA_RBM, -1 + i, -1 + i, nuv));
        //
        std::cout << "i, Norm: " << i << " " << Norm << std::endl;
        //
        std::cout << "iunit: " << iunit << std::endl;
        //
        if (Norm > 0.0)
        {
            scale2dCol(X_QRPA_RBM, -1 + i, nuv, 1 / std::sqrt(Norm));
            scale2dCol(Y_QRPA_RBM, -1 + i, nuv, 1 / std::sqrt(Norm));
        }
        else
        {
            scale2dCol(X_QRPA_RBM, -1 + i, nuv, iunit / std::sqrt(-Norm));
            scale2dCol(Y_QRPA_RBM, -1 + i, nuv, iunit / std::sqrt(-Norm));
        }
        //
        Norm = std::real(dot2dConjCol2dCol(X_QRPA_RBM, X_QRPA_RBM, -1 + i, -1 + i, nuv) - dot2dConjCol2dCol(Y_QRPA_RBM, Y_QRPA_RBM, -1 + i, -1 + i, nuv));
    }
    //
    //
    msgToScreen("XXX_QRPA_RBM:");
    print2d(X_QRPA_RBM, nuv, colldim);
    msgToScreen("XXY_QRPA_RBM:");
    print2d(Y_QRPA_RBM, nuv, colldim);
    msgToScreen("XXF20:");
    print2d(F20, nuv, nop);
    msgToScreen("XXF02:");
    print2d(F02, nuv, nop);
    //
    //
    for (int j = 1; j <= nop; j++)
    {
        for (int i = 1; i <= colldim; i++)
        {
            RBMstrengthfromXY[-1 + i][-1 + j] = std::pow((std::abs(dot2dConjCol2dCol(X_QRPA_RBM, F20, -1 + i, -1 + j, nuv) + dot2dConjCol2dCol(Y_QRPA_RBM, F02, -1 + i, -1 + j, nuv))), 2);
            // std::cout << "Check dot: " << dot2dConjCol2dCol(X_QRPA_RBM, F20, -1 + i, -1 + j, nuv) << " " << dot2dConjCol2dCol(Y_QRPA_RBM, F02, -1 + i, -1 + j, nuv) << std::endl;
        }
    }
    //
    msgToScreen("RBMstrengthfromXY:");
    print2d(RBMstrengthfromXY, colldim, nop);
}

void rbm_METHODS::rbmOutputFile()
{
    std::ofstream outp(rbm_outputfile); // Replace "rbm_outputfile.txt" with the actual output file name

    if (outp.is_open())
    {
        outp << "Results of the RBM" << std::endl;
        outp << "1idx          2ReEnergy          3ImEnergy         4ReStr(SuuS)    "
             << "5ImStr(SuuS)    6StrfromXY(op1)    7StrfromXY(op2) ..." << std::endl;

        for (int i = 1; i <= colldim; ++i)
        {
            outp << std::setw(5) << -1 + i + 1;

            for (int j = 1; j <= nop; ++j)
            {
                outp << std::setw(20) << std::setprecision(10) << RBMenergy[-1 + SortedOrder_Hcoll[-1 + i]];
                outp << std::setw(20) << std::setprecision(10) << RBMstrength[-1 + SortedOrder_Hcoll[-1 + i]];
                outp << std::setw(20) << std::setprecision(10) << RBMstrengthfromXY[-1 + SortedOrder_Hcoll[-1 + i]][-1 + j];
            }

            outp << std::endl;
        }

        outp.close();
    }
}
