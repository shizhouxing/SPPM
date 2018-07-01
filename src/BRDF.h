#ifndef BRDF_H
#define BRDF_H
#include "Ray.h"

class BRDF {
public:
    double specular, diffuse, refraction, rho_d, rho_s, phong_s, refractiveIndex;
    BRDF() {}
    BRDF(double specular, double diffuse, double refraction, 
        double rho_d, double rho_s, double phong_s, double refractiveIndex) {
        this->specular = specular;
        this->diffuse = diffuse;
        this->refraction = refraction;
        this->rho_d = rho_d; this->rho_s = rho_s; this->phong_s = phong_s;
        this->refractiveIndex = refractiveIndex;
    }
};

enum { 
    DIFFUSE, MIRROR, GLASS, LIGHT, MARBLE, FLOOR, WALL, DESK, STANFORD_MODEL, WATER, TEAPOT 
} BRDF_TYPES;

const BRDF BRDFs[] = {
    BRDF(0, 1, 0,       0.7, 0, 0,      0), // DIFFUSE
    BRDF(1, 0, 0,       0, 0, 0,        0), // MIRROR
    BRDF(0, 0, 1,       0, 0, 0,        1.65), // GLASS
    BRDF(0, 1, 0,       0, 0, 0,        0), // LIGHT
    BRDF(0.1, 0.9, 0,   1, 0, 50,       0 ), // MARBLE
    BRDF(0.1, 0.9, 0,   0.9, 0.1, 50,   0), // FLOOR
    BRDF(0, 1, 0,       1, 0, 0,        0), // WALL
    BRDF(0, 1, 0,       1, 0, 0,        0), // DESK
    BRDF(0, 1, 0,       0.9, 0.1, 10,   1), // STANFORD_MODEL
    BRDF(0, 0, 1,       0, 0, 0,        1.3), // WATER
    BRDF(0, 0, 1,       0, 0, 0,        1.5), // TEAPOT
};

#endif
