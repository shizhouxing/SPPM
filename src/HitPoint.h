#ifndef HITPOINT_H
#define HITPOINT_H
#include "Vector.h"
#include "BRDF.h"

class HitPoint {
public:
    Point p;
    Color weight, flux, fluxLight;
    Vector d, norm;
    int n;
    BRDF brdf;
    double r2;
    bool valid;
    HitPoint() {
        flux = fluxLight = Color(0, 0, 0);
        r2 = Config::initial_radius;
        n = 0;
        valid = false;
    }
};

class HitPointKDTreeNode {
public:
    HitPoint *hitpoint;
    Vector min, max;
    double maxr2;
    HitPointKDTreeNode *ls, *rs;
};

class HitPointKDTree {
    int n;
    HitPoint** hitpoints;
    HitPointKDTreeNode* build(int l, int r, int d);
    void del(HitPointKDTreeNode *p);
public:
    HitPointKDTreeNode *root;
    HitPointKDTree(vector<HitPoint*>* hitpoints);
    ~HitPointKDTree();
    void update(HitPointKDTreeNode *p, Point photon, Color weight, Vector d);
};

bool cmpHitPointX(HitPoint *a, HitPoint *b);
bool cmpHitPointY(HitPoint *a, HitPoint *b);
bool cmpHitPointZ(HitPoint *a, HitPoint *b);

#endif
