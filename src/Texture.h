#ifndef TEXTURE_H
#define TEXTURE_H
#include "Vector.h"

class Texture {
    Color** image;
    int width, height;
public:
    Texture(char *filename) {
        import(filename);
    }
    void import(char *filename);
    Color query(double x, double y);
};

class TextureMapper {
    Texture *texture;
    Color color;
    double xx, xy, xz, xb, yx, yy, yz, yb;
public:
    TextureMapper(Texture *texture, double xx, double xy, double xz, double xb, double yx, double yy, double yz, double yb) {
        this->texture = texture;
        this->xx = xx;
        this->xy = xy;
        this->xz = xz;
        this->xb = xb;
        this->yx = yx;
        this->yy = yy;
        this->yz = yz;
        this->yb = yb;
    }
    TextureMapper(Color color) {
        this->texture = nullptr;
        this->color = color;
    }
    Color query(Point p);
};

#endif
