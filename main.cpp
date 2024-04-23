#include <cmath>
#include "scene.cpp"

int main() {
    int W = 512;
    int H = 512;

    // standard scene
    Sphere top = Sphere(Vector(-1000, 0, 0), 940, Vector(1, 0, 0));
    Sphere bottom = Sphere(Vector(1000, 0, 0), 990, Vector(0, 0, 1));
    Sphere left = Sphere(Vector(0, -1000, 0), 940, Vector(1, 1, 0));
    Sphere right = Sphere(Vector(0, 1000, 0), 940, Vector(0, 1, 1));
    Sphere front = Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere back = Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere sphere_1 = Sphere(Vector(0, -20, 0), 10, Vector(1, 1, 1));
    sphere_1.set_mirror();
    Sphere sphere_2 = Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    sphere_2.set_refractive_index(1.5);
    Sphere sphere_3 = Sphere(Vector(0, 20, 0), 10, Vector(1, 1, 1));
    sphere_3.set_refractive_index(2);
    Sphere sphere_4 = Sphere(Vector(0, 20, 0), 9, Vector(1, 1, 1));
    sphere_4.set_refractive_index(1);
    std::vector<Sphere> spheres = {sphere_1, sphere_2, sphere_3, sphere_4, top, bottom, left, right, front, back};
    Scene scene = Scene(spheres);

    //camera
    Vector Q = Vector(0, 0, 55);
    double alpha = 60 * (M_PI / 180);

    //light
    Vector S = Vector(-20, -10, 40);
    double I = 2 * std::pow(10, 10);

    // number of rays launched
    int K = 5;
 
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = Ray(Q);
            ray.compute_direction(i, j, W, H, Q, alpha);
            double red = 0, green = 0, blue = 0;
            for (int k = 0; k < K; ++k) {
                Vector color = scene.get_color(S, I, ray, 10, 1);
                red += color[0] / K;
                green += color[1] / K;
                blue += color[2] / K;
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(red, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(green, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(blue, 1./2.2));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
