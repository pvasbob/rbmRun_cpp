#pragma once

#include "rbm_variables.h"

class rbm_METHODS : public rbm_VARIABLES
{
public:
    void read_rbmfam_NAMELIST();
    void read_fam_training();
};
