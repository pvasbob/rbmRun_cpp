#include "ReadComplex.h"

#include <fstream>

// Function to read a complex number from a string
std::complex<double> readComplex(std::ifstream &file)
{
    char dummy;
    double real, imag;

    file >> dummy;                 // Read the opening parenthesis
    file >> real >> dummy >> imag; // Read real, comma, imag
    file >> dummy;                 // Read the closing parenthesis
    // std::cout << real << " " << imag << std::endl;
    return std::complex<double>(real, imag);
}
