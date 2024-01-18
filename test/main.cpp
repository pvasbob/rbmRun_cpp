#include <iostream>

#include "MultiDimArrayAllocate.h"

int main()
{
    // int **d2;
    //
    // d2 = new int *[10];
    // for (int row = 0; row < 10; row++)
    //{
    //    d2[row] = new int[12];
    //}
    //
    // std::cout << sizeof(d2) / sizeof(d2[0]) << std::endl;
    // std::cout << sizeof(d2[0]) / sizeof(d2[0][0]) << std::endl;

    // Example 2D array
    int rows = 5;
    int cols = 3;
    int **array = nullptr;
    // array = new int *[rows];
    // for (int i = 0; i < rows; ++i)
    //{
    //     array[i] = new int[cols];
    // }

    allocate2dArray<int>(array, rows, cols);

    // Get the number of rows
    int num_rows = sizeof(array) / sizeof(array[0]);

    std::cout << "sizeof(array): " << sizeof(array) << std::endl;
    std::cout << "sizeof(array[0]): " << sizeof(array[0]) << std::endl;
    std::cout << "sizeof(array[0][0]): " << sizeof(array[0][0]) << std::endl;

    std::cout << "Number of rows: " << num_rows << std::endl;

    // Don't forget to free the allocated memory
    for (int i = 0; i < rows; ++i)
    {
        delete[] array[i];
    }
    delete[] array;

    return 0;
}
