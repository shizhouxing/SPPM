#include <cstdio>
#include "Camera.h"

Object* genWalls() {
    Texture *textureBottom = new Texture("../data/bottom.ppm");
    Texture *textureBack = new Texture("../data/back.ppm");
    Texture *textureTop = new Texture("../data/top.ppm");
    
    TextureMapper *color_back = new TextureMapper(Color(1, 1, 1));
    TextureMapper *color_front = new TextureMapper(Color(0.8, 0.2, 0.2));
    TextureMapper *color_left = new TextureMapper(Color(0.2, 0.5, 0.8));
    TextureMapper *color_right = new TextureMapper(Color(0.2, 0.8, 0.5));
    
    Object *walls = new Object;
    walls->numVertexes = 8;
    walls->numFaces = 12;
    walls->vertexes = new Point*[walls->numVertexes]{
        new Point(-1.2, 1, -2.5),
        new Point(-1.2, 1, 1.5),
        new Point(1.2, 1, 1.5),
        new Point(1.2, 1, -2.5),
        new Point(-1.2, -1, -2.5),
        new Point(-1.2, -1, 1.5),
        new Point(1.2, -1, 1.5),
        new Point(1.2, -1, -2.5)
    };
    for (int i = 0; i < walls->numVertexes; ++i) 
        *walls->vertexes[i] = *walls->vertexes[i] * 0.5;
    Point** vertexes = walls->vertexes;
    walls->faces = new Face*[walls->numFaces]{
        // back
        new TriangularFace(
            vertexes[1], vertexes[2], vertexes[5], 
            new TextureMapper(textureBack, 0, -1, 0, 0.5, 1/1.2, 0, 0, 0.5), MARBLE
        ),
        new TriangularFace(
            vertexes[2], vertexes[5], vertexes[6], 
            new TextureMapper(textureBack, 0, -1, 0, 0.5, 1/1.2, 0, 0, 0.5), MARBLE
        ),
        // top
        new TriangularFace(
            vertexes[0], vertexes[1], vertexes[2], 
            new TextureMapper(textureTop, 0, 0, 1./2, 2.5/4, 1./1.2, 0, 0, 0.5), DIFFUSE
        ),
        new TriangularFace(
            vertexes[0], vertexes[2], vertexes[3], 
            new TextureMapper(textureTop, 0, 0, 1./2, 2.5/4, 1./1.2, 0, 0, 0.5), DIFFUSE
        ),
        // bottom
        new TriangularFace(
            vertexes[4], vertexes[5], vertexes[6], 
            new TextureMapper(textureBottom, 0, 0, 1./2, 2.5/4, 1./1.2, 0, 0, 0.5), FLOOR
        ),
        new TriangularFace(
            vertexes[4], vertexes[6], vertexes[7], 
            new TextureMapper(textureBottom, 0, 0, 1./2, 2.5/4, 1./1.2, 0, 0, 0.5), FLOOR
        ),
        // left
        new TriangularFace(vertexes[0], vertexes[1], vertexes[4], color_left, WALL),
        new TriangularFace(vertexes[1], vertexes[4], vertexes[5], color_left, WALL),
        // right
        new TriangularFace(vertexes[2], vertexes[3], vertexes[6], color_right, WALL),
        new TriangularFace(vertexes[3], vertexes[6], vertexes[7], color_right, WALL),
        // front
        new TriangularFace(vertexes[0], vertexes[4], vertexes[7], color_front, WALL),
        new TriangularFace(vertexes[0], vertexes[3], vertexes[7], color_front, WALL),
    };
    return walls;
}

Object* genDesk() {
    Texture *textureBottom = new Texture("../data/bottom.ppm");
    
    TextureMapper *color = new TextureMapper(Color(1, 1, 1));
    
    Object *desk = new Object;
    desk->numVertexes = 4;
    desk->numFaces = 2;
    const double theta = 45 / 180. * M_PI;
    desk->vertexes = new Point*[desk->numVertexes]{
        new Point(-10, 0 - 10 * sin(theta), -10 * cos(theta)),
        new Point(-10, 0 + 10 * sin(theta), 10 * cos(theta)),
        new Point(10, 0 + 10 * sin(theta), 10 *  cos(theta)),
        new Point(10, 0 - 10 * sin(theta), -10 * cos(theta))
    };
    Point** vertexes = desk->vertexes;
    desk->faces = new Face*[desk->numFaces]{
        // bottom
        new TriangularFace(
            vertexes[0], vertexes[1], vertexes[2], 
            color, DESK
        ),
        new TriangularFace(
            vertexes[0], vertexes[2], vertexes[3], 
            color, DESK
        ),
    };
    return desk;
}

Object* genLight(Point p, double r) {
    Object *light = new Object;
    light->numFaces = 1;
    light->faces = new Face*[1]{
        new DiscFace(p, r, new TextureMapper(Color(1, 1, 1)), LIGHT)
    };
    return light;
}

Scene *sceneBox() {
    Object *bunny = new Object;
    bunny->importPly("../data/bunny.ply",  new TextureMapper(Color(0.4, 0.8, 0.8)), STANFORD_MODEL);    
    
    bunny->scale(
        2.8, 0, 0, 0.24,
        0, 2.8, 0, -0.09 - 0.5,
        0, 0, -2.8, 0.12
    );
    bunny->rotXZ(15 * M_PI / 180);
    bunny->center->print();
    
    Object *water = new Object;
    water->importPly("../data/water.ply",  new TextureMapper(Color(1, 1, 1)), WATER);    
    water->scale(
        1. / 5.52799 * 1.2, 0, 0, -0.6,
        0, 0.12 / (1.85354 - 1.34492), 0, -0.31731 + 0.08,
        0, 0, 1.2 / (5.59200 + 0.00456), -1.19902 + 0.75
    );
    water->center = new Point(0, -1, 0);
    water->printBox();
    
    Object *teapot = new Object;
    teapot->importBpt("../data/teapot.bpt", new TextureMapper(Color(0.1, 0.8, 0.8)), LIGHT); 
    teapot->scale(
        -0.02, 0, 0, 0,
        0, 0, 0.02, 0,
        0, 0.02, 0, 0
    );
    teapot->printBox();    
    
    Scene *scene = new Scene(Point(0.0, 0.5 - 1e-5, 0.1), 0.2, Vector(0, -1, 0));
    scene->addObject(bunny);
    scene->addObject(water);
    scene->addObject(genWalls());
    scene->addObject(genLight(Point(0, 0.5 - Config::epsilon, 0.1), 0.2));
    scene->addObject(new Sphere(
        Point(-0.32, -0.30, 0.3), 0.18, new TextureMapper(Color(1, 1, 1)), GLASS));
    scene->addObject(new Sphere(
        Point(0.42, 0.20, 0), 0.15, new TextureMapper(Color(1, 1, 1)), MIRROR));
    return scene;
}

Scene *sceneTeapot() {
    Object *water = new Object;
    water->importPly("../data/water.ply",  new TextureMapper(Color(1, 1, 1)), WATER);    
    water->scale(
        1. / 5.52799 * 1.2, 0, 0, -0.6,
        0, 0.12 / (1.85354 - 1.34492), 0, -0.31731 + -0.05,
        0, 0, 1.2 / (5.59200 + 0.00456), -1.19902 + 0.75
    );
    water->center = new Point(0, -1, 0);
    water->printBox();    
    
    Object *teapot = new Object;
    teapot->importBpt("../data/teapot.bpt", new TextureMapper(Color(0.1, 0.8, 0.8)), TEAPOT); 
    teapot->scale(
        -0.08, 0, 0, 0.02,
        0, 0, 0.08, -0.05,
        0, 0.08, 0, 0.1
    );
    teapot->printBox();    
    
    Scene *scene = new Scene(Point(0.0, 0.5 - 1e-5, 0.1), 0.2, Vector(0, -1, 0));
    scene->addObject(genWalls());
    scene->addObject(teapot);
    scene->addObject(water);
    scene->addObject(genLight(Point(0, 0.5 - Config::epsilon, 0.1), 0.2));
    return scene;
}

int main(int argc, char *argv[]) {
    Scene *scene = sceneBox();
    //Scene *scene = sceneTeapot();
    
    Camera *camera = new Camera(1024, 768);
    camera->setScene(scene);
    camera->setPosition(Point(0, 0.15, -1));
    camera->setLens(0.684, 0.811, 1e-3, 1 + 0.09);
    
    if (argc > 2)
        camera->render(stoi(argv[1]), stoi(argv[2]));
    else
        camera->render(stoi(argv[1]));
        
    camera->save("result.ppm");
    return 0;
}
