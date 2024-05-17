#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "libranlib.h"

#ifndef ___RANDOM_H___
#define ___RANDOM_H___


class Random_Gen {
    public:
        Random_Gen();
        Random_Gen(std::string seed_string);
        ~Random_Gen();
    protected:
        long seed1, seed2;
        std::ifstream seedin;
        std::ofstream seedout;
};


#endif
