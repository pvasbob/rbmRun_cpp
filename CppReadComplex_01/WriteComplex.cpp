#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

int main()
{
    // Declare variables
    double value1 = 0.01;
    int value2 = 11;
    int value3 = 1;
    std::string str1 = "emulatortraining.dat";
    std::string str2 = "emulatorRBMoutput.dat";
    std::string str3 = "emulator.dat";
    int value4 = 1;
    bool boolValue = true;
    std::complex<double> complexNumbers[] = {{0.0, 0.5}, {0.1, 0.0}};
    int value5 = 100;
    bool boolValues[] = {false, false, false, false, true, true, true, true};

    // Open a file for writing
    std::ofstream outputFile("output_file_complex_numbers.txt");

    // Check if the file is open
    if (!outputFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Write variables to the file with parentheses
    outputFile << std::fixed << std::setprecision(2);

    outputFile << value1 << "\n";
    outputFile << value2 << "\n";
    outputFile << value3 << "\n";
    outputFile << str1 << "\n";
    outputFile << str2 << "\n";
    outputFile << str3 << "\n";
    outputFile << value4 << "\n";
    outputFile << std::boolalpha << boolValue << "\n";

    for (const auto &complexNumber : complexNumbers)
    {
        outputFile << "(" << complexNumber.real() << ",  " << std::setw(4) << complexNumber.imag() << ")\n";
    }

    outputFile << value5 << "\n";
    for (int i = 0; i < 8; ++i)
    {
        outputFile << std::boolalpha << boolValues[i] << " ";
    }

    // Close the file
    outputFile.close();

    std::cout << "Output written to file successfully." << std::endl;

    return 0;
}
