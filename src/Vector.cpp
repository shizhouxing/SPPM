#include "Vector.h"

double sqr(double x) {
    return x * x;
}

Vector operator*(const Matrix &a, const Vector &b) {
    return Vector(
        a.v[0][0] * b.x + a.v[0][1] * b.y + a.v[0][2] * b.z,
        a.v[1][0] * b.x + a.v[1][1] * b.y + a.v[1][2] * b.z,
        a.v[2][0] * b.x + a.v[2][1] * b.y + a.v[2][2] * b.z
    );
}

Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector operator*(const Vector &a, const Vector &b) {
    return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vector operator/(const Vector &a, const Vector &b) {
    return Vector(a.x / b.x, a.y / b.y, a.z / b.z);
}

Vector operator*(const Vector &a, const double &b) {
    return Vector(a.x * b, a.y * b, a.z * b);
}

Vector operator/(const Vector &a, const double &b) {
    return Vector(a.x / b, a.y / b, a.z / b);
}

Vector min(const Vector &a, const Vector &b) {
    return Vector(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

Vector max(const Vector &a, const Vector &b) {
    return Vector(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

Vector cross(const Vector &a, const Vector &b) {
    return Vector(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

bool operator<=(const Vector &a, const Vector &b) {
    return a.x <= b.x + Config::epsilon && a.y <= b.y + Config::epsilon && a.z <= b.z + Config::epsilon;
}

bool operator>=(const Vector &a, const Vector &b) {
    return a.x + Config::epsilon >= b.x && a.y + Config::epsilon >= b.y && a.z + Config::epsilon >= b.z;
}

double dot(const Vector &a, const Vector &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double det(const Vector &a, const Vector &b, const Vector &c) {
    return a.x * (b.y * c.z - b.z * c.y) 
        - b.x * (a.y * c.z - a.z * c.y) 
        + c.x * (a.y * b.z - a.z * b.y);
}

Matrix Matrix::inv() {
    Matrix T(n);
    T.v[0][0] = at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2);
    T.v[0][1] = at(1, 3) * at(3, 2) - at(1, 2) * at(3, 3);
    T.v[0][2] = at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2);
    T.v[1][0] = at(2, 3) * at(3, 1) - at(2, 1) * at(3, 3);
    T.v[1][1] = at(1, 1) * at(3, 3) - at(1, 3) * at(3, 1);
    T.v[1][2] = at(2, 1) * at(1, 3) - at(1, 1) * at(2, 3);
    T.v[2][0] = at(2, 1) * at(3, 2) - at(2, 2) * at(3, 1);
    T.v[2][1] = at(1, 2) * at(3, 1) - at(1, 1) * at(3, 2);
    T.v[2][2] = at(1, 1) * at(2, 2) - at(2, 1) * at(1, 2);
    double d = at(1, 1) * (at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2))
        - at(2, 1) * (at(1, 2) * at(3, 3) - at(1, 3) * at(3, 2))
        + at(3, 1) * (at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2));
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) 
            T.v[i][j] /= d;
    return T;
}
