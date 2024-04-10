#include "ray.cpp"

class Sphere {
public:
    Vector C; //center
    double R; //radius
    Vector albedo; //color

    explicit Sphere(Vector center, double radius, Vector color) {
        C = center;
        R = radius;
        albedo = color;
    }

    bool intersect_t(Ray ray, double &t) {
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

    // computes the point of intersection between a Ray and the sphere, if any 
    bool intersect(Ray ray, double &t, Vector &P, Vector &N) {
        if (intersect_t(ray, t)) {
            P = ray.O + t * ray.u;
            N = (P - C) / (P - C).norm();
            return true;
        }
        return false;
    }
};
