#pragma once
#include <complex>

void diagGenComplexMat(int &dim,
                       std::complex<double> **mat,
                       std::complex<double> **pmat,
                       std::complex<double> *eigenvalue,
                       std::complex<double> **pmatinv,
                       int &info);