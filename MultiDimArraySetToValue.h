#pragma once
#include <complex>

// Set 1d array to value.
template <typename T>
void set1dArrayToValue(T *&array, int dim1, T value)
{
    for (int i = 0; i < dim1; i++)
    {
        array[i] = value;
    }
}

// Set 2d array to value.
template <typename T>
void set2dArrayToValue(T **&array, int dim1, int dim2, T value)
{
    for (int i = 0; i < dim1; i++)
    {
        for (int j = 0; j < dim2; j++)
        {
            array[i][j] = value;
        }
    }
}
//
// copy col from 2d to col in another 2d
template <typename T>
void copy2dColTo2dCol(T **&array_src, T **&array_dest, int dim1, int col_src, int col_dest)
{
    for (int row = 0; row < dim1; row++)
    {
        array_dest[row][col_dest] = array_src[row][col_src];
    }
}
//
// copy conjgated col from 2d to col in another 2d
template <typename T>
void copy2dConjColTo2dCol(T **&array_src, T **&array_dest, int dim1, int col_src, int col_dest)
{
    for (int row = 0; row < dim1; row++)
    {
        array_dest[row][col_dest] = std::conj(array_src[row][col_src]);
    }
}
//
// copy minus, conjgated col from 2d to col in another 2d
template <typename T>
void copy2dMConjColTo2dCol(T **&array_src, T **&array_dest, int dim1, int col_src, int col_dest)
{
    for (int row = 0; row < dim1; row++)
    {
        array_dest[row][col_dest] = -std::conj(array_src[row][col_src]);
    }
}
