//
// 2D array memory allocation
template <typename T>
void allocate2dArray(T **array, int dim1, int dim2)
{
    array = new T *[dim1];
    for (int i = 0; i < dim1; i++)
    {
        array[i] = new T[dim2];
    }
}

// 2D array memory Deallocation
template <typename T>
void deallocate2dArray(T **array, int dim1)
{
    for (int i = 0; i < dim1; i++)
    {
        delete[] array[i];
    }
    delete[] array;
}
