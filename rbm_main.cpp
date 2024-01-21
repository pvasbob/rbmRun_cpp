#include <iostream>
#include "rbm_methods.h"

int main()
{

    // std::cout << "Hello World!" << std::endl;
    rbm_METHODS rbm_methods;
    //!----------------------------------------------------------------------------------------------------
    //! Read Namelist
    //!----------------------------------------------------------------------------------------------------
    rbm_methods.read_rbmfam_NAMELIST();

    std::cout << "**** Parameters on training inputfile ****" << std::endl;
    std::cout << "fam_training_inputfile       = " << rbm_methods.fam_training_inputfile << std::endl;
    std::cout << "number_of_training_points    = " << rbm_methods.number_of_training_points << std::endl;
    std::cout << "number_of_operators          = " << rbm_methods.number_of_operators << std::endl;
    std::cout << std::endl;

    std::cout << "**** Parameters on RBM calculation ****" << std::endl;
    std::cout << "mirror_points                = " << rbm_methods.mirror_points << std::endl;
    std::cout << "VARIATION                    = " << rbm_methods.VARIATION << std::endl;
    std::cout << "normcut                      = " << rbm_methods.normcut << std::endl;
    std::cout << "rbm_outputfile               = " << rbm_methods.rbm_outputfile << std::endl;
    std::cout << std::endl;

    std::cout << "**** Parameters on FAM Emulator Run ****" << std::endl;
    std::cout << "STRENGTH_EMULATOR_RUN        = " << rbm_methods.STRENGTH_EMULATOR_RUN << std::endl;
    std::cout << "qrpa_omega_emulatorrun_start = " << rbm_methods.qrpa_omega_emulatorrun_start << std::endl;
    std::cout << "qrpa_omega_emulatrorun_step  = " << rbm_methods.qrpa_omega_emulatorrun_step << std::endl;
    std::cout << "nmax_emulatorrun             = " << rbm_methods.nmax_emulatorrun << std::endl;
    std::cout << "emulator_outputfile          = " << rbm_methods.emulator_outputfile << std::endl;

    //// Display the read values
    // std::cout << "main Read values from the file:" << std::endl;
    // std::cout << rbm_methods.normcut << std::endl;
    // std::cout << rbm_methods.number_of_training_points << std::endl;
    // std::cout << rbm_methods.number_of_operators << std::endl;
    // std::cout << rbm_methods.fam_training_inputfile << std::endl;
    // std::cout << rbm_methods.rbm_outputfile << std::endl;
    // std::cout << rbm_methods.emulator_outputfile << std::endl;
    // std::cout << rbm_methods.VARIATION << std::endl;
    // std::cout << rbm_methods.STRENGTH_EMULATOR_RUN << std::endl;
    //
    //// Display complex numbers
    // std::cout << "Complex numbers:" << std::endl;
    //// for (const auto &complexNumber : complexNumbers)
    ////{
    ////     std::cout << "(" << complexNumber.real() << ", " << complexNumber.imag() << ")" << std::endl;
    //// }
    //
    // std::cout << "(" << rbm_methods.qrpa_omega_emulatorrun_start.real() << ", " << rbm_methods.qrpa_omega_emulatorrun_start.imag() << ")" << std::endl;
    // std::cout << "(" << rbm_methods.qrpa_omega_emulatorrun_step.real() << ", " << rbm_methods.qrpa_omega_emulatorrun_step.imag() << ")" << std::endl;
    // std::cout << rbm_methods.nmax_emulatorrun << std::endl;
    //
    //// Display bool values
    // std::cout << "Bool values:" << std::endl;
    // for (int i = 0; i < 4; ++i)
    //{
    //     std::cout << std::boolalpha << rbm_methods.mirror_points[i] << std::endl;
    // }

    rbm_methods.read_fam_training();
    rbm_methods.completeData();
    rbm_methods.strengthAtTraining();
    rbm_methods.kernelCalculation();
    //
    rbm_methods.diagNormKernel();
    rbm_methods.sortNormEigen();
    rbm_methods.normEigenCutoff();
    rbm_methods.normKernelExcludingSmallEigen();
    rbm_methods.normKernelHalf();
    rbm_methods.Hcoll();
    rbm_methods.diagCollectiveHamiltonian();
    rbm_methods.calculateStrength();
    rbm_methods.calculateQRPAXY();
    //
    rbm_methods.rbmOutputFile();
}
