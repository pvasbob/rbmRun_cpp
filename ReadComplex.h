#pragma once

#include <complex>

std::complex<double> readComplex(std::ifstream &file);

void readComplex1d(std::ifstream &file, std::complex<double> *&array, int dim1);

// dim1 is exactly number of rows.
void readComplexToCol(std::ifstream &file, std::complex<double> **&array, int dim1, int col_index);
