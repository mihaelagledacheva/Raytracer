#include <limits>
#include <random>
#include "geometries.cpp"

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

    Vector get_intensity(Vector P, Vector albedo, Vector N, double I, Vector S) {
        double d = (S - P).norm();
        Vector w = (S - P) / d;
        //visibility
        int V = 1;
        Ray ray = Ray(P + 0.0001 * N);
        ray.set_direction(w);
        double t;
        if (intersect_t(ray, t)) {
            if (t <= d) {
                V = 0;
            }
        }
        return (I / (4 * M_PI * std::pow(d, 2))) * (albedo / M_PI) * V * std::max(0., dot(N, w));
    }

    Vector get_color(Vector S, double I, Ray ray, int ray_depth, double n) {
        if (ray_depth < 0) {
            return Vector(0, 0, 0);
        }
        double t;
        Vector P;
        Vector N;
        size_t sphere_id;
        if (intersect(ray, t, P, N, sphere_id)) {
            if (spheres[sphere_id].mirror) {
                // Reflection
                Ray reflected_ray = Ray(P + 0.0001 * N);
                Vector u_r = ray.u - 2 * dot(ray.u, N) * N;
                reflected_ray.set_direction(u_r);
                return get_color(S, I, reflected_ray, ray_depth-1, n);
            } else if (spheres[sphere_id].transparent) {
                // Refraction
                Ray refracted_ray = Ray(P + 0.001 * N);
                double cos_theta_i = dot(ray.u, N);
                double sin_theta_i = sqrt(1.0 - std::pow(cos_theta_i, 2));
                N = (cos_theta_i < 0) ? N : -N;
                double n1 = (cos_theta_i < 0) ? n : spheres[sphere_id].n;
                double n2 = (cos_theta_i < 0) ? spheres[sphere_id].n : n;
                // Fresnel law
                double k0 = std::pow(n1 - n2, 2) / std::pow(n1 + n2, 2);
                double R = k0 + (1 - k0) * std::pow(1 - std::abs(dot(N, ray.u)), 5);
                double u = static_cast<double>(std::rand()) / RAND_MAX;
                if (sin_theta_i * n1 / n2 <= 1.0 && u >= R) {
                    Vector w_T = (n1 / n2) * (ray.u - dot(ray.u, N) * N);
                    Vector w_N = -N * sqrt(1 - std::pow(n1 / n2, 2) * (1 - std::pow(dot(ray.u, N), 2)));
                    refracted_ray.set_direction(w_T + w_N);
                    return get_color(S, I, refracted_ray, ray_depth-1, n2);
                } else {
                    // Total internal reflection
                    Vector u = ray.u - 2 * dot(ray.u, N) * N;
                    refracted_ray.set_direction(u);
                    return get_color(S, I, refracted_ray, ray_depth-1, n1);
                }
            } else {
                return get_intensity(P, spheres[sphere_id].albedo, N, I, S);
            }
        }
        return Vector(0, 0, 0);
    }
};
