#include "Utils.h"

vector<string> Utils::split(char *str) {
    vector<string> res;
    string cur = "";
    for (int i = 0; str[i]; ++i) {
        if (str[i] == ' ' || str[i] == '\n' || str[i] == '\r') {
            if (cur != "") res.push_back(cur);
            cur = "";
        }
        else 
            cur += str[i];
    }
    if (cur != "") res.push_back(cur);
    return res;
}

bool isPrime(int x) {
    for (int i = 2; i * i <= x; ++i) 
        if (x % i == 0) 
            return false;
    return x > 1;
}

double Utils::randomQMC(int axis, long long i) {
    int base = prime[axis];
    double f = 1, res = 0;
    while (i > 0) {
        f /= base;
        res += f * (i % base);
        i /= base;
    }
    return res;
}

static double Utils::random01() {
    static mt19937 *generator = nullptr;
    if (!generator) 
        generator = new mt19937(clock());
    static uniform_real_distribution<> dis(0, 1);
    return dis(*generator);
}

double Utils::random(double l, double r, int axis, long long i) {
    if (axis == -1)
        return l + random01() * (r - l);
    return l + randomQMC(axis, i) * (r - l);
}
