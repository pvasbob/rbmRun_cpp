#pragma once

#include <complex>

std::complex<double> readComplex(std::ifstream &file);
void readComplex1d(std::ifstream &file, std::complex<double> *&array, int dim1);
