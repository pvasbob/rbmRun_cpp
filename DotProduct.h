#pragma once
#include <complex>

#include "MultiDimArraySetToValue.h"

// dot product of 2d col with 2d col
std::complex<double> dot2dCol2dCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1);
// dot product of 2d conjugated col with 2d col
std::complex<double> dot2dConjCol2dCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1);
// dot product of 2d conjugated col with 2d  conjugated col
std::complex<double> dot2dConjCol2dConjCol(std::complex<double> **&array1, std::complex<double> **&array2, int col_index1, int col_index2, int dim1);
// ==============================
// 2d and 2d  multiplication
// ==============================
template <typename T>
void mult2dAnd2d(T **&array_srcL, T **&array_srcR, T **&array_dest, int &m, int &n, int &l)
{
    // first set array_dest to zero incase assigned previously.
    set2dArrayToValue<T>(array_dest, m, l, 0.0);
    //
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < l; j++)
        {
            for (int k = 0; k < n; k++)
            {
                array_dest[i][j] += array_srcL[i][k] * array_srcR[k][j];
            }
        }
    }
}
//
// 1d and 2d multiplication.
template <typename T>
void mult1dAnd2d(T *&array_srcL, T **&array_srcR, T *&array_dest, int &m, int &n)
{
    set1dArrayToValue<T>(array_dest, n, 0.0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            array_dest[i] += array_srcL[j] * array_srcR[j][i];
        }
    }
}
//
// 2d and 1d multiplcation.
template <typename T>
void mult2dAnd1d(T **&array_srcL, T *&array_srcR, T *&array_dest, int &m, int &n)
{
    set1dArrayToValue<T>(array_dest, m, 0.0);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            array_dest[i] += array_srcL[i][j] * array_srcR[j];
        }
    }
}
//
// 2d and Conj 1d multiplcation.
template <typename T>
void mult2dAndConj1d(T **&array_srcL, T *&array_srcR, T *&array_dest, int &m, int &n)
{
    set1dArrayToValue<T>(array_dest, m, 0.0);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            array_dest[i] += array_srcL[i][j] * std::conj(array_srcR[j]);
        }
    }
}