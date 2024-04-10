#include <limits>
#include "geometries.cpp"
#include <iostream>
class Scene {
public:
    std::vector<Sphere> spheres;

    explicit Scene(std::vector<Sphere> array) {
        spheres = array;
    }
    
    bool intersect_t(Ray ray, double &t) {
        t = std::numeric_limits<double>::infinity();
        bool intersection = false;
        for (auto sphere : spheres) {
            double t_sphere;
            if (sphere.intersect_t(ray, t_sphere)) {
                t = std::min(t, t_sphere);
                intersection = true;
            }
        }
        return intersection;
    }

    // computes the point of intersection between a Ray and the scene, if any 
    bool intersect(Ray ray, double &t, Vector &P, Vector &N, size_t &sphere_id) {
        t = std::numeric_limits<double>::infinity();
        bool intersection = false;
        for (size_t i = 0; i < spheres.size(); ++i) {
            double t_sphere;
            Vector P_sphere;
            Vector N_sphere;
            if (spheres[i].intersect(ray, t_sphere, P_sphere, N_sphere)) {
                if (t_sphere < t) {
                    t = t_sphere;
                    P = P_sphere;
                    N = N_sphere;
                    sphere_id = i;
                    intersection = true;
                }
            }
        }
        return intersection;
    }
};
