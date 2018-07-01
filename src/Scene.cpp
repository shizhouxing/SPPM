#include "Utils.h"
#include "Scene.h"

bool ObjectKDTreeNode::inside(Face *face) {
    Vector faceMin = face->min();
    Vector faceMax = face->max();
    return (faceMin.x < max.x || faceMin.x == max.x && faceMin.x == faceMax.x)
        && (faceMax.x > min.x || faceMax.x == min.x && faceMin.x == faceMax.x)
        && (faceMin.y < max.y || faceMin.y == max.y && faceMin.y == faceMax.y)
        && (faceMax.y > min.y || faceMax.y == min.y && faceMin.y == faceMax.y)
        && (faceMin.z < max.z || faceMin.z == max.z && faceMin.z == faceMax.z)
        && (faceMax.z > min.z || faceMax.z == min.z && faceMin.z == faceMax.z);
}

ObjectKDTreeNode* ObjectKDTree::build(int depth, int d, vector<Face*>* faces, Vector min, Vector max) {
    ObjectKDTreeNode *p = new ObjectKDTreeNode;
    p->min = min;
    p->max = max;
    Vector maxL, minR;
    if (d == 0) {
        maxL = Vector((p->min.x + p->max.x) / 2, p->max.y, p->max.z);
        minR = Vector((p->min.x + p->max.x) / 2, p->min.y, p->min.z);
    }
    else if (d == 1) {
        maxL = Vector(p->max.x, (p->min.y + p->max.y) / 2, p->max.z);
        minR = Vector(p->min.x, (p->min.y + p->max.y) / 2, p->min.z);
    }
    else {
        maxL = Vector(p->max.x, p->max.y, (p->min.z + p->max.z) / 2);
        minR = Vector(p->min.x, p->min.y, (p->min.z + p->max.z) / 2);
    }
    p->faces = new vector<Face*>;
    for (auto face : *faces)
        if (p->inside(face)) 
            p->faces->push_back(face);
         
    const int max_faces = 8;
    const int max_depth = 24;
    
    if (p->faces->size() > max_faces && depth < max_depth) {
        p->ls = build(depth + 1, (d + 1) % 3, p->faces, min, maxL);
        p->rs = build(depth + 1, (d + 1) % 3, p->faces, minR, max);
        
        vector<Face*> *faceL = p->ls->faces, *faceR = p->rs->faces;
        map<Face*, int> cnt;
        for (auto face : *faceL) cnt[face]++;
        for (auto face : *faceR) cnt[face]++;
        p->ls->faces = new vector<Face*>;
        p->rs->faces = new vector<Face*>;
        p->faces->clear();
        for (auto face : *faceL) 
            if (cnt[face] == 1) 
                p->ls->faces->push_back(face);
            else
                p->faces->push_back(face);
        for (auto face : *faceR) 
            if (cnt[face] == 1) 
                p->rs->faces->push_back(face);
    }
    else
        p->ls = p->rs = nullptr;
    return p;
}

void ObjectKDTree::getFaces(ObjectKDTreeNode *p, vector<Face*>* faces) {
    p->l = faces->size();
    for (auto face : *(p->faces))
        faces->push_back(face);
    p->r = faces->size();
    if (p->ls) getFaces(p->ls, faces);
    if (p->rs) getFaces(p->rs, faces);
}

ObjectKDTree::ObjectKDTree(vector<Face*>* faces) {
    Vector min = Vector(1e100, 1e100, 1e100);
    Vector max = min * -1;        
    for (auto face : *faces) {
        min = ::min(min, face->min());
        max = ::max(max, face->max());
    }
    root = build(1, 0, faces, min, max);
    this->faces = new vector<Face*>;
    getFaces(root, this->faces);
}

double ObjectKDTree::getCuboidIntersection(ObjectKDTreeNode *p, Ray ray) {
    if (!(ray.s >= p->min && ray.s <= p->max)) { // outside
        double t = -1e100;
        if (fabs(ray.d.x) > 0)
            t = max(t, min((p->min.x - ray.s.x) / ray.d.x, (p->max.x - ray.s.x) / ray.d.x));
        if (fabs(ray.d.y) > 0)
            t = max(t, min((p->min.y - ray.s.y) / ray.d.y, (p->max.y - ray.s.y) / ray.d.y));
        if (fabs(ray.d.z) > 0)
            t = max(t, min((p->min.z - ray.s.z) / ray.d.z, (p->max.z - ray.s.z) / ray.d.z));
        if (t < -Config::epsilon) return 1e100;
        Point pp = ray.s + ray.d * t;
        if (!(pp >= p->min && pp <= p->max)) return 1e100;
        return t;
    }
    else return -1e100;
}

void ObjectKDTree::getIntersection(
        ObjectKDTreeNode *p, Ray ray, Face* &nextFace, double &tMin, Vector &norm) {
    for (int i = 0; i < p->faces->size(); ++i) {
        pair<double, Vector> r = (*p->faces)[i]->intersect(ray);
        double t = r.first;
        if (t > 0 && t < tMin) {
            tMin = t;
            nextFace = (*p->faces)[i];
            norm = r.second;
        }
    }
    
    double tl = p->ls ? getCuboidIntersection(p->ls, ray) : 1e100;
    double tr = p->rs ? getCuboidIntersection(p->rs, ray) : 1e100;
    if (tl < tr) {
        if (tMin <= tl) return;
        if (p->ls) getIntersection(p->ls, ray, nextFace, tMin, norm);
        if (tMin <= tr) return;
        if (p->rs) getIntersection(p->rs, ray, nextFace, tMin, norm);
    }
    else {
        if (tMin <= tr) return;
        if (p->rs) getIntersection(p->rs, ray, nextFace, tMin, norm);
        if (tMin <= tl) return;
        if (p->ls) getIntersection(p->ls, ray, nextFace, tMin, norm);            
    }
}

void Scene::addObject(Object* object) {
    for (int i = 0; i < object->numFaces; ++i)
        object->faces[i]->object = object;
    objects.push_back(object);
}

Ray Scene::generateRay(long long i) {
    double alpha = Utils::random(0, 2 * M_PI, 0, i);
    Point s = sourceP +  Point(cos(alpha), 0, sin(alpha)) * sourceR;
    Vector d = sampleReflectedRay(sourceN, 0, i);
    return Ray(s + d * Config::epsilon, d);
}

Vector Scene::sampleReflectedRay(Vector norm, int depth, long long i, double s) {
    Vector u = cross(Vector(1, 0, 0), norm);
    if (u.norm2() < Config::epsilon) u = cross(Vector(0, 1, 0), norm);
    u.normalize();
    Vector v = cross(norm, u);
    v.normalize();
    double theta = Utils::random(0, 2 * M_PI, 2 * depth + 1, i);
    double phi = asin(pow(Utils::random(0, 1, 2 * depth + 2, i), 1. / (s + 1)));
    return (norm * cos(phi) + (u * cos(theta) + v * sin(theta)) * sin(phi)).normalize();
}

void Scene::trace(const Ray &ray, const Color &weight, int depth, long long i, HitPoint *hp) {
    if (depth > Config::max_tracing_depth)
        return;
    double tMin = 1e100;
    Face* nextFace = nullptr;
    Vector norm;
    
    objectKDTree->getIntersection(objectKDTree->root, ray, nextFace, tMin, norm);
    
    if (!nextFace || !hp && tMin < 1e-3) return;
    Point p = ray.s + ray.d * tMin;

    // russian roulette
    double s = BRDFs[nextFace->brdf].specular + BRDFs[nextFace->brdf].diffuse + BRDFs[nextFace->brdf].refraction;
    double action = Utils::random(0, 1) * s;
    
    Vector dr = ray.d - norm * (2 * dot(ray.d, norm));    
    
    // specular
    if (BRDFs[nextFace->brdf].specular > 0 && action <= BRDFs[nextFace->brdf].specular) {
        trace(
            Ray(p + dr * Config::epsilon, dr), 
            weight * nextFace->texture->query(p) * s, depth + 1, i, hp
        );        
        return;
    }
    action -= BRDFs[nextFace->brdf].specular;
    
    // diffuse
    if (BRDFs[nextFace->brdf].diffuse > 0 && action <= BRDFs[nextFace->brdf].diffuse) {
        if (hp) {
            hp->p = p;
            hp->weight = weight * nextFace->texture->query(p) * s;
            hp->fluxLight = hp->fluxLight + hp->weight * (nextFace->brdf == LIGHT);
            hp->brdf = BRDFs[nextFace->brdf];
            hp->norm = norm;
            if (nextFace->brdf == LIGHT) {
                hp->fluxLight = hp->fluxLight + hp->weight;
                hp->valid = false;
            }
            else
                hp->valid = true;
        }   
        else {
            double a = Utils::random();
            // phong specular
            if (a <= BRDFs[nextFace->brdf].rho_s) { 
                Vector d = sampleReflectedRay(dr, depth, i, BRDFs[nextFace->brdf].phong_s);
                trace(
                    Ray(p + d * Config::epsilon, d), 
                    weight * nextFace->texture->query(p) * s,
                    depth + 1, i, hp
                );
            }
            else {
                a -= BRDFs[nextFace->brdf].rho_s;
                hitpointsKDTree->update(hitpointsKDTree->root, p, weight, ray.d);
                Vector d = sampleReflectedRay(norm, depth, i);
                if (dot(d, norm) < 0) d = d * -1;
                if (a <= BRDFs[nextFace->brdf].rho_d) {
                    trace(
                        Ray(p + d * Config::epsilon, d), 
                        weight * nextFace->texture->query(p) * s, 
                        depth + 1, i, hp
                    );
                }   
            }
        }    
        return;
    }
    action -= BRDFs[nextFace->brdf].diffuse;
    
    // refraction
    if (BRDFs[nextFace->brdf].refraction > 0 && action <= BRDFs[nextFace->brdf].refraction) {
        if (!nextFace->object->center) 
            nextFace->object->calcCenter();
        bool incoming = dot(*(nextFace->object)->center - p, norm) < 0;
        double refractiveIndex = BRDFs[nextFace->brdf].refractiveIndex;
        if (!incoming) refractiveIndex = 1. / refractiveIndex;
        double cosThetaIn = -dot(ray.d, norm);
        double cosThetaOut2 = 1 - (1 - sqr(cosThetaIn)) / sqr(refractiveIndex);
        
        if (cosThetaOut2 >= -Config::epsilon) {
            double cosThetaOut = sqrt(cosThetaOut2);
            
            // schlick's approximation
            double R0 = sqr((1 - refractiveIndex) / (1 + refractiveIndex));
            double cosTheta = incoming ? cosThetaIn : cosThetaOut;
            double R = R0 + (1 - R0) * pow(1 - cosTheta, 5);
            
            if (Utils::random() <= R)
                trace(
                    Ray(p + dr * Config::epsilon, dr), 
                    weight * nextFace->texture->query(p) * s, depth + 1, i, hp
                );
            else {
                Vector d = ray.d / refractiveIndex
                        + norm * (cosThetaIn / refractiveIndex - cosThetaOut);
                trace(
                    Ray(p + d * Config::epsilon, d),
                    weight * nextFace->texture->query(p) * s, depth + 1, i, hp
                );                                    
            }
        }
        else { // total internal reflection
            trace(
                Ray(p + dr * Config::epsilon, dr), 
                weight * nextFace->texture->query(p) * s, depth + 1, i, hp
            );        
        }
    }
}

void Scene::initializeHitpointKDTree(vector<HitPoint*>* hitpoints) {
    if (hitpointsKDTree) 
        delete hitpointsKDTree;
    hitpointsKDTree = new HitPointKDTree(hitpoints);
    fprintf(stderr, "Hitpoint KD tree built\n");
}

void Scene::initializeObjectKDTree() {
    vector<Face*> *faces = new vector<Face*>;
    for (auto object : objects) {
        for (int i = 0; i < object->numFaces; ++i) 
            faces->push_back(object->faces[i]);
    }
    objectKDTree = new ObjectKDTree(faces);
    
    fprintf(stderr, "Object KD tree built\n");
}
