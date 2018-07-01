#ifndef CAMERA_H
#define CAMERA_H
#include "Scene.h"

class Camera {
    int w, h;    
    double fx, fy, aperture, focus;
    Color **canvas, light;
    Scene* scene;
    Point s;
    vector<HitPoint*> *hitpoints;
    void evaluateRadiance(int numRounds, int numPhotons);    
public:
    Camera(int w, int h) {
        this->w = w;
        this->h = h;
        
        canvas = new Color*[w];
        for (int i = 0; i < w; ++i) 
            canvas[i] = new Color[h];
    }
    void setScene(Scene *scene) {
        this->scene = scene;
    }
    void setPosition(Point s) {
        this->s = s;
    }
    void setLens(double fx, double fy, double aperture, double focus) {
        this->fx = fx;
        this->fy = fy;
        this->aperture = aperture;
        this->focus = focus;
    }
    void render(int numRounds, int numPhotons = 20480);
    void save(char *filename);
};

#endif
