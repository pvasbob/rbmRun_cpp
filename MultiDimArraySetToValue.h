#pragma once

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
