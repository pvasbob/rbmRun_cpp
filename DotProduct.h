#pragma once
#include <complex>

// dot product of 2d col with 2d col
std::complex<double> dot2dCol2dCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1);
// dot product of 2d conjugated col with 2d col
std::complex<double> dot2dConjCol2dCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1);
// dot product of 2d conjugated col with 2d  conjugated col
std::complex<double> dot2dConjCol2dConjCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1);
