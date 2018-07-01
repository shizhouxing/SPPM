#include "Camera.h"

void Camera::evaluateRadiance(int numRounds, int numPhotons) {
    for (int u = 0; u < w; ++u)
        for (int v = 0; v < h; ++v) {
            HitPoint *hp = (*hitpoints)[u * h + v];
            canvas[u][v] = hp->flux / (M_PI * hp->r2 * numPhotons * numRounds) + 
                light * hp->fluxLight / numRounds;
            canvas[u][v].x = pow(canvas[u][v].x, Config::gamma);
            canvas[u][v].y = pow(canvas[u][v].y, Config::gamma);
            canvas[u][v].z = pow(canvas[u][v].z, Config::gamma);
            canvas[u][v] = min(canvas[u][v], Vector(1, 1, 1));
    }   
}

void Camera::render(int numRounds, int numPhotons) {
    scene->initializeObjectKDTree();
    
    int cx = w / 2, cy = h / 2;
        
    TriangularFace focusPlane(new Point(0, 0, s.z + focus), new Point(0, 1, s.z + focus), new Point(1, 0, s.z + focus));
    
    // initialize hitpoints
    hitpoints = new vector<HitPoint*>;
    for (int u = 0; u < w; ++u) 
        for (int v = 0; v < h; ++v)
            hitpoints->push_back(new HitPoint);
    
    light = Color(1, 1, 1);
    const Color weight_init = Color(2.5, 2.5, 2.5);
    
    // SPPM
    for (int round = 0; round < numRounds; ++round) {
        fprintf(stderr, "Round %d/%d:\n", round + 1, numRounds);
        
        // ray tracing pass
        for (int u = 0; u < w; ++u) {
            fprintf(stderr, "\rRay tracing pass %.3lf%%", u * 100. / w);
            #pragma omp parallel for schedule(dynamic, 60), num_threads(8)
            for (int v = 0; v < h; ++v) {
                Point p(double(u - cx) / w / fx, double(v - cy) / h / fy, 0);
                Ray ray(s, p - s);
                double t = focusPlane.intersectPlane(ray);
                Point focusP = ray.s + ray.d * t;
                double theta = Utils::random(0, 2 * M_PI);
                ray.s = ray.s + Point(cos(theta), sin(theta), 0) * aperture;
                ray.d = focusP - ray.s;
                ray.d.normalize();
                (*hitpoints)[u * h + v]->valid = false;
                (*hitpoints)[u * h + v]->d = ray.d * -1;
                scene->trace(ray, Color(1, 1, 1), 1, (long long)round * (numPhotons + w * h) + u * h + v, (*hitpoints)[u * h + v]);
            }
        }
        fprintf(stderr, "\rRay tracing pass 100.000%%\n");
        
        scene->initializeHitpointKDTree(hitpoints);
        
        // photon tracing pass
        #pragma omp parallel for schedule(dynamic, 128), num_threads(8)
        for (int i = 0; i < numPhotons; ++i) {
            Ray ray = scene->generateRay((long long)round * numPhotons + (round + 1) * w * h + i);
            scene->trace(ray, weight_init * light, 1, (long long)round * numPhotons + i); 
        }
        fprintf(stderr, "\rPhoton tracing pass done\n");       
        
        if ((round + 1) % Config::checkpoint_interval == 0) {
            evaluateRadiance(round + 1, numPhotons);     
            char filename[100];
            sprintf(filename, "checkpoint-%d.ppm", round + 1);
            save(filename);
        }
    }    
    evaluateRadiance(numRounds, numPhotons);     
}

void Camera::save(char *filename) {
    FILE *file = fopen(filename, "w");
    fprintf(file, "P3\n%d %d\n255\n", w, h);
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) 
            fprintf(file, 
                "%d %d %d ", 
                int(canvas[j][h - i - 1].x * 255 + 0.5), 
                int(canvas[j][h - i - 1].y * 255 + 0.5), 
                int(canvas[j][h - i - 1].z * 255 + 0.5)
            );
        fprintf(file, "\n");
    }
    fclose(file);
    fprintf(stderr, "Image saved to %s\n", filename);
}
