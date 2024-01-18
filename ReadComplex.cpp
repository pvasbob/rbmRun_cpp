#include "ReadComplex.h"

#include <fstream>
#include <iostream>

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

// Function to read complex array from a string and assign to array
void readComplex1d(std::ifstream &file, std::complex<double> *&array, int dim1)
{
    std::cout << "Inside readComplex1d." << std::endl;
    char dummy;
    double real, imag;
    for (int i = 0; i < dim1; i++)
    {
        // file >> dummy;
        // file >> real >> dummy >> imag;
        // file >> dummy;
        //
        // std::cout << array << std::endl;
        // std::cout << readComplex(file) << std::endl;
        array[i] = readComplex(file);
        // std::cout << real << imag << std::endl;
    }
}
