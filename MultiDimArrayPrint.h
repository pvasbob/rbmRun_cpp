#pragma once
#include <iostream>

// print out elements of a col in 2d matrix.
template <typename T>
void print2dCol(T **&array, int dim1, int col_index)
{
    for (int row = 0; row < dim1; row++)
    {
        std::cout << "#row, #col " << row << " " << col_index << " " << array[row][col_index] << std::endl;
    }
}

// print out 2d elements for every col.
template <typename T>
void print2d(T **&array, int dim1, int dim2)
{
    for (int col = 0; col < dim2; col++)
    {
        print2dCol<T>(array, dim1, col);
    }
}

// print a 1d array
template <typename T>
void print1d(T *array, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        std::cout << " index: " << i << " " << array[i] << std::endl;
    }
}