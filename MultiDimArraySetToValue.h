#pragma once
#include <complex>
#include <iostream>

// ===============================================================
// set array to value
// ===============================================================
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
// ===============================================================
// copy col to col
// ===============================================================
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
// =================================================================
// copy array to array
// =================================================================
// copy 1d to 1d
template <typename T>
void copy1dTo1d(T *&array_src, T *&array_dest, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        array_dest[i] = array_src[i];
    }
}
// copy 2d to 2d
template <typename T>
void copy2dTo2d(T **&array_src, T **&array_dest, int dim_row, int dim_col)
{
    for (int row = 0; row < dim_row; row++)
    {
        for (int col = 0; col < dim_col; col++)
        {
            array_dest[row][col] = array_src[row][col];
        }
    }
}
// copy 1d to 2d, assume row major.
template <typename T>
void copy1dTo2dRowMajor(T *&array_src, T **&array_dest, int dim_src, int dim_dest_row, int dim_dest_col)
{
    if (dim_dest_row * dim_dest_col != dim_src)
    {
        std::cerr << "Totally elements of src and dest NOT match. Quit" << std::endl;
        exit(1);
    }
    // row major.
    int ct = 0;
    for (int row = 0; row < dim_dest_row; row++)
    {
        for (int col = 0; col < dim_dest_col; col++)
        {
            array_dest[row][col] = array_src[ct];
            ct++;
        }
    }
}
// col major
template <typename T>
void copy1dTo2dColMajor(T *&array_src, T **&array_dest, int dim_src, int dim_dest_row, int dim_dest_col)
{
    if (dim_dest_row * dim_dest_col != dim_src)
    {
        std::cerr << "Totally elements of src and dest NOT match. Quit" << std::endl;
        exit(1);
    }
    // row major.
    int ct = 0;
    for (int col = 0; col < dim_dest_col; col++)
    {
        for (int row = 0; row < dim_dest_row; row++)
        {
            array_dest[row][col] = array_src[ct];
            ct++;
        }
    }
}
// / copy 2d to 1d. assume row major.
template <typename T>
void copy2dTo1dRowMajor(T **&array_src, T *&array_dest, int dim_src_row, int dim_src_col, int dim_dest)
{
    if (dim_src_row * dim_src_col != dim_dest)
    {
        std::cerr << "Totally elements of src and dest NOT match. Quit" << std::endl;
        exit(1);
    }
    // row major
    int ct = 0;
    for (int row = 0; row < dim_src_row; row++)
    {
        for (int col = 0; col < dim_src_col; col++)
        {
            array_dest[ct] = array_src[row][col];
            ct++;
        }
    }
}
//
template <typename T>
void copy2dTo1dColMajor(T **&array_src, T *&array_dest, int dim_src_row, int dim_src_col, int dim_dest)
{
    if (dim_src_row * dim_src_col != dim_dest)
    {
        std::cerr << "Totally elements of src and dest NOT match. Quit" << std::endl;
        exit(1);
    }
    // col major
    int ct = 0;
    for (int col = 0; col < dim_src_col; col++)
    {
        for (int row = 0; row < dim_src_row; row++)
        {
            array_dest[ct] = array_src[row][col];
            ct++;
        }
    }
}

// ===============================================
// set 2d diagonal elements to value.
// ===============================================
//  only have one dim for 2d becuases the concept of 'diagonal elements' only applys to squre matrix.
template <typename T>
void setDiagToValue(T **&array, int dim, T value)
{
    for (int i = 0; i < dim; i++)
    {
        array[i][i] = value;
    }
}

// ================================================
// copy 1d real part to another 1d with settted sign
// ================================================
template <typename T>
void copy1dRealTo1d(std::complex<T> *&array_src, T *&array_dest, int dim, char pm)
{
    int sign = +1;
    if (pm == '-')
        sign = -1;
    //
    for (int i = 0; i < dim; i++)
    {
        array_dest[i] = sign * real(array_src[i]);
    }
}
// =================================================
// scale a col by value
// =================================================
template <typename T>
void scale2dCol(T **&array, int col_index, int dim, std::complex<double> value)
{
    for (int i = 0; i < dim; i++)
    {
        array[i][col_index] = array[i][col_index] / value;
    }
}