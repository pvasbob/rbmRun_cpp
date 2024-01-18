#pragma once
#include <complex>
//
// 1D array memory allocation
template <typename T>
void allocate1dArray(T *&array, int dim1)
{
    array = new T[dim1];
}

// 2D array memory allocation
template <typename T>
void allocate2dArray(T **&array, int dim1, int dim2)
{
    array = new T *[dim1];
    for (int i = 0; i < dim1; i++)
    {
        array[i] = new T[dim2];
    }
}

// 1D array memory Deallocation
template <typename T>
void deallocate1dArray(T *&array, int dim1)
{
    delete[] array;
}

// 2D array memory Deallocation
template <typename T>
void deallocate2dArray(T **&array, int dim1)
{
    for (int i = 0; i < dim1; i++)
    {
        delete[] array[i];
    }
    delete[] array;
}

// std::complex<double> *allocate1dArraytest(int dim1);
// void allocate1dArray(std::complex<double> *&array, int dim1);
