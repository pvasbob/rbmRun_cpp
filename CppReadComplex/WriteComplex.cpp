#include <iostream>
#include <fstream>
#include <complex>
#include <iomanip>

int main()
{
    // Create a complex number
    std::complex<double> myComplexNumber(1.0, 2.0);

    // Open a file for writing
    std::ofstream outputFile("complex_output.txt");

    // Check if the file is open
    if (!outputFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Write the complex number to the file with a specific format
    outputFile << "(" << std::fixed << std::setprecision(2) << myComplexNumber.real()
               << " , " << std::fixed << std::setprecision(2) << myComplexNumber.imag() << ")";

    // Close the file
    outputFile.close();

    std::cout << "Complex number written to file successfully." << std::endl;

    return 0;
}
