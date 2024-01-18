#include "ReadInputToCol.h"
#include "ReadComplex.h"

// read input to column cod_index.
void readInputToCol(std::ifstream &file, std::complex<double> **&array, int dim1, int col_index)
{
    for (int row = 0; row < dim1; row++)
    {
        // read in one complex and assign to the corresponding 2d matrix element.
        array[row][col_index] = readComplex(file);
    }
}