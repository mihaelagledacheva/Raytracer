#include "ray.cpp"

class Sphere {
public:
    Vector C; //center
    double R; //radius
    Vector albedo; //color
    bool light;
    double i;
    bool mirror;
    bool transparent;
    double n_in;
    double n_out;

    explicit Sphere(Vector center, double radius, Vector color) {
        C = center;
        R = radius;
        albedo = color;
        light = false;
        mirror = false;
        transparent = false;
    }

    void set_light(double intensity) {
        light = true;
        i = intensity;
    }

    void set_mirror() {
        mirror = true;
    }

    void set_refractive_index(double refractive_index_in, double refractive_index_out = 1) {
        transparent = true;
        n_in = refractive_index_in;
        n_out = refractive_index_out;
    }

    bool intersect(Ray &ray, double &t) {
        double delta = std::pow(dot(ray.u, ray.O - C), 2) - ((ray.O - C).norm2() - std::pow(R, 2));
        if (delta < 0) {
            return false;
        } else {
            double t1 = dot(ray.u, C - ray.O) - sqrt(delta);
            double t2 = dot(ray.u, C - ray.O) + sqrt(delta);
            if (t2 < 0) {
                return false;
            } else {
                if (t1 >= 0) {
                    t = t1;
                } else {
                    t = t2;
                }
                return true;
            }
        }
    }

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N) {
        if (intersect(ray, t)) {
            P = ray.O + t * ray.u;
            N = (P - C) / (P - C).norm();
            return true;
        }
        return false;
    }
};
