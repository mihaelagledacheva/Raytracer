#include <cmath>
#include <chrono>
#include <iostream>
#include <omp.h>
#include "scene.cpp"

//image
int W = 512;
int H = 512;

//scene
Sphere* top = new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 0, 0));
Sphere* bottom = new Sphere(Vector(1000, 0, 0), 990, Vector(0, 0, 1));
Sphere* left = new Sphere(Vector(0, -1000, 0), 940, Vector(1, 1, 0));
Sphere* right = new Sphere(Vector(0, 1000, 0), 940, Vector(0, 1, 1));
Sphere* front = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
Sphere* back = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
std::vector<Geometry*> stage = {top, bottom, left, right, front, back};

//camera
Vector Q = Vector(0, 0, 55);
double alpha = 60 * (M_PI / 180);

//light
Vector S = Vector(-20, -10, 40);
double I = 2 * std::pow(10, 10);
int D = 100;


void lab_1_1() {
    // setup
    std::vector<Geometry*> spheres(stage.begin(), stage.end());
    Sphere solid_sphere = Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    spheres.push_back(&solid_sphere);
    Scene scene = Scene(spheres);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = Ray(Q);
            ray.compute_direction(i, j, W, H, Q, alpha);
            Vector color = scene.get_color_unrefined(S, I, ray);
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 1./2.2));
        }
    }
    stbi_write_png("_results/lab_1_1.png", W, H, 3, &image[0], 0);
}

void lab_1_2(int K=1000) {
    // setup
    std::vector<Geometry*> spheres(stage.begin(), stage.end());
    Sphere reflective_sphere = Sphere(Vector(0, -10, 0), 10, Vector(1, 1, 1));
    reflective_sphere.set_mirror();
    spheres.push_back(&reflective_sphere);
    Sphere refractive_sphere = Sphere(Vector(0, 10, 0), 10, Vector(1, 1, 1));
    refractive_sphere.set_refractive_index(1.5);
    spheres.push_back(&refractive_sphere);
    Scene scene = Scene(spheres);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Ray ray = Ray(Q);
            ray.compute_direction(i, j, W, H, Q, alpha);
            double red = 0, green = 0, blue = 0;
            for (int k = 0; k < K; ++k) {
                Vector color = scene.get_color_unrefined(S, I, ray, 5);
                red += color[0] / K;
                green += color[1] / K;
                blue += color[2] / K;
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(red, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(green, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(blue, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_1_2.png", W, H, 3, &image[0], 0);
}

void lab_2_1(int K=1000) {
    // setup
    std::vector<Geometry*> spheres(stage.begin(), stage.end());
    Sphere reflective_sphere = Sphere(Vector(0, -10, 0), 10, Vector(1, 1, 1));
    reflective_sphere.set_mirror();
    spheres.push_back(&reflective_sphere);
    Sphere refractive_sphere = Sphere(Vector(0, 10, 0), 10, Vector(1, 1, 1));
    refractive_sphere.set_refractive_index(1.5);
    spheres.push_back(&refractive_sphere);
    Scene scene = Scene(spheres);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(S, I, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_2_1.png", W, H, 3, &image[0], 0);
}

void lab_2_2(int K=1000) {
    // setup
    std::vector<Geometry*> spheres(stage.begin(), stage.end());
    Sphere light = Sphere(S, 10, Vector(1, 1, 1));
    light.set_light(I);
    spheres.push_back(&light);
    Sphere solid_sphere = Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    spheres.push_back(&solid_sphere);
    Scene scene = Scene(spheres);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_2_2.png", W, H, 3, &image[0], 0);
}

void lab_2_3(int K=5000) {
    // setup
    std::vector<Geometry*> spheres(stage.begin(), stage.end());
    Sphere light = Sphere(Vector (-25, -10, -25), 5, Vector(1, 1, 1));
    light.set_light(I);
    spheres.push_back(&light);
    Sphere reflective_sphere = Sphere(Vector(0, -10, 0), 10, Vector(1, 1, 1));
    reflective_sphere.set_mirror();
    spheres.push_back(&reflective_sphere);
    Sphere refractive_sphere = Sphere(Vector(0, 10, 0), 10, Vector(1, 1, 1));
    refractive_sphere.set_refractive_index(1.5);
    spheres.push_back(&refractive_sphere);
    Scene scene = Scene(spheres);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_2_3.png", W, H, 3, &image[0], 0);
}

void lab_2_4(int K=2000) {
    // setup
    std::vector<Geometry*> spheres(stage.begin(), stage.end());
    Sphere light = Sphere(S, 5, Vector(1, 1, 1));
    light.set_light(I);
    spheres.push_back(&light);
    Sphere refractive_sphere = Sphere(Vector(0, -10, 20), 10, Vector(1, 1, 1));
    refractive_sphere.set_refractive_index(1.5);
    spheres.push_back(&refractive_sphere);
    Sphere solid_sphere = Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    spheres.push_back(&solid_sphere);
    Sphere reflective_sphere = Sphere(Vector(0, 10, -20), 10, Vector(1, 1, 1));
    reflective_sphere.set_mirror();
    spheres.push_back(&reflective_sphere);
    Scene scene = Scene(spheres);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                Vector P = Q + ray.u * D / std::abs(ray.u[2]);
                double r = sqrt(uniform(engine));
                double theta = uniform(engine) * 2 * M_PI;
                Vector Qprime = Vector(r*cos(theta), r*sin(theta), Q[2]);
                Ray rayprime = Ray(Qprime);
                Vector uprime = (P - Qprime) / (P - Qprime).norm();
                rayprime.set_direction(uprime);
                color = color + scene.get_color(light, rayprime, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_2_4.png", W, H, 3, &image[0], 0);
}

void lab_3_1(int K=32) {
    // setup
    std::vector<Geometry*> geometries(stage.begin(), stage.end());
    Sphere light = Sphere(S, 5, Vector(1, 1, 1));
    light.set_light(I);
    geometries.push_back(&light);
    TriangleMesh cat = TriangleMesh();
    cat.readOBJ("cat_model/cat.obj");
    cat.scale(0.7);
    cat.translate(Vector(-5, -15, -5));
    cat.rotate_xy(M_PI/2);
    cat.set_parameters(false, false);
    geometries.push_back(&cat);
    Scene scene = Scene(geometries);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_3_1.png", W, H, 3, &image[0], 0);
}

void lab_3_2(int K=32) {
    // setup
    std::vector<Geometry*> geometries(stage.begin(), stage.end());
    Sphere light = Sphere(S, 5, Vector(1, 1, 1));
    light.set_light(I);
    geometries.push_back(&light);
    TriangleMesh cat = TriangleMesh();
    cat.readOBJ("cat_model/cat.obj");
    cat.scale(0.7);
    cat.translate(Vector(-5, -15, -5));
    cat.rotate_xy(M_PI/2);
    cat.set_parameters(false, false);
    cat.bound();
    geometries.push_back(&cat);
    Scene scene = Scene(geometries);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_3_2.png", W, H, 3, &image[0], 0);
}

void lab_4_1(int K=32) {
    // setup
    std::vector<Geometry*> geometries(stage.begin(), stage.end());
    Sphere light = Sphere(S, 5, Vector(1, 1, 1));
    light.set_light(I);
    geometries.push_back(&light);
    TriangleMesh cat = TriangleMesh();
    cat.readOBJ("cat_model/cat.obj");
    cat.scale(0.7);
    cat.translate(Vector(-5, -15, -5));
    cat.rotate_xy(M_PI/2);
    cat.set_parameters(false, false);
    cat.bvh(100);
    geometries.push_back(&cat);
    Scene scene = Scene(geometries);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_4_1.png", W, H, 3, &image[0], 0);
}

void lab_4_2(int K=1000) {
    // setup
    std::vector<Geometry*> geometries(stage.begin(), stage.end());
    Sphere light = Sphere(S, 5, Vector(1, 1, 1));
    light.set_light(I);
    geometries.push_back(&light);
    TriangleMesh cat = TriangleMesh();
    cat.readOBJ("cat_model/cat.obj");
    cat.scale(0.7);
    cat.translate(Vector(-5, -15, -5));
    cat.rotate_xy(M_PI/2);
    cat.set_parameters(true, false);
    cat.bvh(100);
    geometries.push_back(&cat);
    Scene scene = Scene(geometries);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_4_2.png", W, H, 3, &image[0], 0);
}

void lab_4_3(int K=1000) {
    // setup
    std::vector<Geometry*> geometries(stage.begin(), stage.end());
    Sphere light = Sphere(S, 5, Vector(1, 1, 1));
    light.set_light(I);
    geometries.push_back(&light);
    TriangleMesh cat = TriangleMesh();
    cat.readOBJ("cat_model/cat.obj");
    cat.scale(0.7);
    cat.translate(Vector(-10, -15, 0));
    cat.rotate_xy(M_PI/2);
    cat.rotate_yz(M_PI/3);
    cat.set_parameters(true, true, "cat_model/cat_diff.png");
    cat.bvh(100);
    geometries.push_back(&cat);
    Scene scene = Scene(geometries);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0, 0, 0);
            for (int k = 0; k < K; ++k) {
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double x = sqrt(-2 * log (r1)) * cos(2 * M_PI * r2);
                double y = sqrt(-2 * log (r1)) * sin(2 * M_PI * r2);
                Ray ray = Ray(Q);
                ray.compute_direction(i + x, j + y, W, H, Q, alpha);
                color = color + scene.get_color(light, ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/K, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/K, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/K, 1./2.2));
        }
    }
    stbi_write_png("_results/lab_4_3.png", W, H, 3, &image[0], 0);
}

int main() {
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::seconds::rep duration;

    std::cout << "Rendering basic scene" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_1_1();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_1_1.png" << std::endl << std::endl;

    std::cout << "Implementing reflexion and refraction" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_1_2();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_1_2.png" << std::endl << std::endl;

    std::cout << "Adding indirect lighting" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_2_1();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_2_1.png" << std::endl << std::endl;

    std::cout << "Replacing point light with area (scene 1)" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_2_2();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_2_2.png" << std::endl << std::endl;
    
    std::cout << "Replacing point light with area (scene 2)" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_2_3();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_2_3.png" << std::endl << std::endl;

    std::cout << "Creating depth of field" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_2_4();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_2_4.png" << std::endl << std::endl;

    std::cout << "Rendering cat model" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_3_1();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_3_1.png" << std::endl << std::endl;

    std::cout << "Using a bounding box" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_3_2();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_3_2.png" << std::endl << std::endl;

    std::cout << "Implementing bounding volume hierarchies" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_4_1();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_4_1.png" << std::endl << std::endl;

    std::cout << "Interpolating on normals" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_4_2();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_4_2.png" << std::endl << std::endl;
    
    std::cout << "Adding texture" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_4_3();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_4_3.png" << std::endl << std::endl;

    return 0;
}
