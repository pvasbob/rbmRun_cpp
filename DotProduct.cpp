#include "DotProduct.h"

// dot product of 2d col with 2d col
std::complex<double> dot2dCol2dCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1)
{
    // can not make this as template cause dot product makes no senese for lots of types.
    std::complex<double> sum = 0.0;
    for (int row = 0; row < dim1; row++)
    {
        sum += array1[row][col_index1] * array2[row][col_index2];
    }

    return sum;
}
//
// dot product of 2d conjugated col with 2d col
std::complex<double> dot2dConjCol2dCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1)
{
    // can not make this as template cause dot product makes no senese for lots of types.
    std::complex<double> sum = 0.0;
    for (int row = 0; row < dim1; row++)
    {
        sum += std::conj(array1[row][col_index1]) * array2[row][col_index2];
    }

    return sum;
}
//
// dot product of 2d conjugated col with 2d conjugated col
std::complex<double> dot2dConjCol2dConjCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1)
{
    // can not make this as template cause dot product makes no senese for lots of types.
    std::complex<double> sum = 0.0;
    for (int row = 0; row < dim1; row++)
    {
        sum += std::conj(array1[row][col_index1]) * std::conj(array2[row][col_index2]);
    }

    return sum;
}
