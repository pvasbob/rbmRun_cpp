#include "MultiDimArrayAllocate.h"

#include <iostream>

// std::complex<double> *allocate1dArraytest(int dim1)
//{
//     std::cout << "Inside allocate test" << std::endl;
//     return new std::complex<double>[dim1];
// }

void allocate1dArray(std::complex<double> *&array, int dim1)
{
    array = new std::complex<double>[dim1];
}