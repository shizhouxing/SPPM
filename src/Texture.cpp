#include "Texture.h"

void Texture::import(char *filename) {
    FILE *file = fopen(filename, "r");
    fscanf(file, "%*s%d%d%*d", &width, &height);
    image = new Color*[height];
    for (int i = 0; i < height; ++i) {
        image[i] = new Color[width];
        for (int j = 0; j < width; ++j) {
            int r, g, b;
            fscanf(file, "%d%d%d", &r, &g, &b);
            image[i][j] = Color(r / 255., g / 255., b / 255.);
        }
    }
    fclose(file);
}

Color Texture::query(double x, double y) {
    int _x = x * height;
    int _y = y * width;
    if (_x >= 0 && _x < height && _y >= 0 && _y < width) 
        return image[_x][_y];
    else
        return Color(0, 0, 0);
}

Color TextureMapper::query(Point p) {
    if (!texture) 
        return color;
    double x = xx * p.x + xy * p.y + xz * p.z + xb;
    double y = yx * p.x + yy * p.y + yz * p.z + yb;
    return texture->query(x, y);
}
