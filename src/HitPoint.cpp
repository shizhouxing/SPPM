#include "HitPoint.h"

HitPointKDTreeNode* HitPointKDTree::build(int l, int r, int d) {
    HitPointKDTreeNode *p = new HitPointKDTreeNode;
    p->min = Vector(1e100, 1e100, 1e100);
    p->max = p->min * (-1);
    p->maxr2 = 0;
    for (int i = l; i <= r; ++i) {
        p->min = min(p->min, hitpoints[i]->p);
        p->max = max(p->max, hitpoints[i]->p);
        p->maxr2 = max(p->maxr2, hitpoints[i]->r2);
    }
    int m = l + r >> 1;
    if (d == 0) 
        nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointX);
    else if (d == 1) 
        nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointY);
    else 
        nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointZ);
    p->hitpoint = hitpoints[m];
    if (l <= m - 1) p->ls = build(l, m - 1, (d + 1) % 3); else p->ls = nullptr;
    if (m + 1 <= r) p->rs = build(m + 1, r, (d + 1) % 3); else p->rs = nullptr;
    return p;
}

HitPointKDTree::HitPointKDTree(vector<HitPoint*>* hitpoints) {
    n = hitpoints->size();
    this->hitpoints = new HitPoint*[n];
    for (int i = 0; i < n; ++i)
        this->hitpoints[i] = (*hitpoints)[i];
    root = build(0, n - 1, 0);
}

void HitPointKDTree::del(HitPointKDTreeNode *p) {
    if (p->ls) del(p->ls);
    if (p->rs) del(p->rs);
    delete p;
}

HitPointKDTree::~HitPointKDTree() {
    if (!root) return;
    del(root);
    delete[] hitpoints;
}

void HitPointKDTree::update(HitPointKDTreeNode * p, Point photon, Color weight, Vector d) {
    if (!p) return;
    double mind = 0, maxd = 0;
	if (photon.x > p->max.x) mind += sqr(photon.x - p->max.x);
	if (photon.x < p->min.x) mind += sqr(p->min.x - photon.x);
	if (photon.y > p->max.y) mind += sqr(photon.y - p->max.y);
	if (photon.y < p->min.y) mind += sqr(p->min.y - photon.y);
	if (photon.z > p->max.z) mind += sqr(photon.z - p->max.z);
	if (photon.z < p->min.z) mind += sqr(p->min.z - photon.z);
    if (mind > p->maxr2) return;
    if (p->hitpoint->valid && (photon - p->hitpoint->p).norm2() <= p->hitpoint->r2) {
        HitPoint *hp = p->hitpoint;
        double factor = (hp->n * Config::alpha + Config::alpha) / (hp->n * Config::alpha + 1.);
        Vector dr = d - hp->norm * (2 * dot(d, hp->norm));    
        double rho = hp->brdf.rho_d + hp->brdf.rho_s * pow(dot(dr, hp->d), hp->brdf.phong_s);
        if (rho < 0) rho = 0;
        else if (rho > 1) rho = 1;
        hp->n++;
        hp->r2 *= factor;
        hp->flux = (hp->flux + hp->weight * weight * rho) * factor;  
    }  
    if (p->ls) update(p->ls, photon, weight, d);
    if (p->rs) update(p->rs, photon, weight, d);
    p->maxr2 = p->hitpoint->r2;
    if (p->ls && p->ls->hitpoint->r2 > p->maxr2)
        p->maxr2 = p->ls->hitpoint->r2;
    if (p->rs && p->rs->hitpoint->r2 > p->maxr2)
        p->maxr2 = p->rs->hitpoint->r2;
}

bool cmpHitPointX(HitPoint *a, HitPoint *b) {
    return a->p.x < b->p.x;
}

bool cmpHitPointY(HitPoint *a, HitPoint *b) {
    return a->p.y < b->p.y;
}

bool cmpHitPointZ(HitPoint *a, HitPoint *b) {
    return a->p.z < b->p.z;
}

