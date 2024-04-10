#include <cmath>
#include "scene.cpp"

int main() {
    int W = 512;
    int H = 512;

    // standard scene
    Sphere top = Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
    Sphere bottom = Sphere(Vector(0, -1000, 0), 940, Vector(0, 0, 1));
    Sphere front = Sphere(Vector(0, 0, 1000), 940, Vector(1, 0.75, 0.8));
    Sphere back = Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere center_sphere = Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    std::vector<Sphere> spheres = {center_sphere, top, bottom};
    Scene scene = Scene(spheres);

    //camera
    Vector Q = Vector(0, 0, 55);
    double alpha = 60 * (M_PI / 180);

    //light
    Vector S = Vector(-10, 20, 40);
    double I = 100;
 
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = Ray(Q);
            ray.compute_direction(i, j, W, H, Q, alpha);
            double t;
            Vector P;
            Vector N;
            size_t sphere_id;            
            if (scene.intersect(ray, t, P, N, sphere_id)) {
                image[(i * W + j) * 3 + 0] = 255;
                image[(i * W + j) * 3 + 1] = 255;
                image[(i * W + j) * 3 + 2] = 255;
            }
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
