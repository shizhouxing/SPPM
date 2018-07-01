#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <map>
#include <ctime>
#include <random>
#include <stdexcept>
#include "Config.h"

namespace Utils {
    vector<string> split(char *str);
    
    // Halton Sequence
    const int prime[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97, 101, 103, 107, 109, 113
    };
    double randomQMC(int axis, long long i);
    static double random01();
    double random(double l = 0, double r = 1, int axis = -1, long long i = 0);
}

#endif
