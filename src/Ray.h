#ifndef RAY_H
#define RAY_H
#include "Vector.h"

class Ray {
public:
    Point s;
    Vector d;
    Ray(Point s, Vector d) {
        this->s = s;
        this->d = d / sqrt(d.norm2());
    }
};

#endif
