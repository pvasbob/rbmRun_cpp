#pragma once
#include <complex>

// dim1 is exactly number of rows.
void readInputToCol(std::ifstream &file, std::complex<double> **&array, int dim1, int col_index);