#include <algorithm>
#include <stdexcept>
//
#include "MultiDimArraySort.h"
#include "MsgToScreen.h"

void sort2(const int &dim, const double *vector, int *SortedOrder, const bool &SMALLESTFIRST)
{
    int SortedOrderTemp;
    bool COMPLETED;

    // Initialization of SortedOrder Array
    for (int i = 1; i <= dim; ++i)
    {
        SortedOrder[-1 + i] = i;
    }
    //
    do
    {
        COMPLETED = true;
        for (int i = 1; i <= dim - 1; ++i)
        {
            if (SMALLESTFIRST)
            {
                if (vector[-1 + SortedOrder[-1 + i]] > vector[-1 + SortedOrder[-1 + i + 1]])
                {
                    std::swap(SortedOrder[-1 + i], SortedOrder[-1 + i + 1]);
                    COMPLETED = false;
                }
            }
            else
            {
                if (vector[-1 + SortedOrder[-1 + i]] < vector[-1 + SortedOrder[-1 + i + 1]])
                {
                    std::swap(SortedOrder[-1 + i], SortedOrder[-1 + i + 1]);
                    COMPLETED = false;
                }
            }
        }
        msgToScreen("Inside SORT2:");
    } while (!COMPLETED);

    // Check
    for (int i = 1; i <= dim; ++i)
    {
        if (SortedOrder[-1 + i] > dim || SortedOrder[-1 + i] < 0)
        {
            throw std::runtime_error("SORT2 Error");
        }
    }
}
