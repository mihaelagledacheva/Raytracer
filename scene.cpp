#include <limits>
#include <random>
#include "geometries.cpp"

class Scene {
public:
    std::vector<Sphere> spheres;

    explicit Scene(std::vector<Sphere> array) {
        spheres = array;
    }
    
    bool intersect(Ray &ray, double &t) {
        t = std::numeric_limits<double>::infinity();
        bool intersection = false;
        for (auto sphere : spheres) {
            double t_sphere;
            if (sphere.intersect(ray, t_sphere)) {
                t = std::min(t, t_sphere);
                intersection = true;
            }
        }
        return intersection;
    }

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N, size_t &sphere_id) {
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

    int visibility(Vector &P, Vector &S) {
        Ray ray(P);
        Vector d = (S - P) / (S - P).norm();
        ray.set_direction(d);
        double t;
        if (intersect(ray, t)) {
            if (t <= (S - P).norm()) {
                return 0;
            }
        }
        return 1;
    }
    
    Vector get_intensity(Vector &P, Vector &albedo, Vector &N, double I, Vector &S) {
        double d = (S - P).norm();
        Vector w = (S - P) / d;
        Vector x = P + 0.0001 * N;
        int v = visibility(x, S);
        return (I / (4 * M_PI * std::pow(d, 2))) * (albedo / M_PI) * v * std::max(0., dot(N, w));
    }

    Vector random_cos(const Vector &N) {
        double r1 = uniform(engine);
        double r2 = uniform(engine);
        double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
        double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
        double z = sqrt(r2);
        Vector T1 = N;
        double min_N = std::min(std::min(std::abs(N[0]), std::abs(N[1])), std::abs(N[2]));
        for (int i = 0; i < 3; ++i) {
            if (std::abs(T1[i]) == min_N) {
                T1[i] = 0;
                T1[(i+1)%3] = N[(i+2)%3];
                T1[(i+2)%3] = -N[(i+1)%3];
                break;
            }
        }
        T1.normalize();
        Vector T2 = cross(N, T1);
        return x * T1 + y * T2 + z * N;
    }

    Vector get_color_unrefined(Vector &S, double I, Ray &ray, int ray_depth=1, bool invert=false) {
        if (ray_depth < 0) {
            return Vector(0, 0, 0);
        }
        double t;
        Vector P, N;
        size_t sphere_id;
        if (intersect(ray, t, P, N, sphere_id)) {
            if (spheres[sphere_id].mirror) {
                // Reflection
                ray.reflect(P, N);
                return get_color_unrefined(S, I, ray, ray_depth-1);
            } else if (spheres[sphere_id].transparent) {
                // Refraction
                double n1 = spheres[sphere_id].n_out;
                double n2 = spheres[sphere_id].n_in;
                Ray refracted_ray = refract(ray, P, N, n1, n2, invert);
                return get_color_unrefined(S, I, refracted_ray, ray_depth-1, invert);
            } else {
                // Diffusion
                return get_intensity(P, spheres[sphere_id].albedo, N, I, S);
            }
        }
        return Vector(0, 0, 0);
    }
    
    Vector get_color(Vector &S, double I, Ray &ray, int ray_depth=1, bool invert=false) {
        if (ray_depth < 0) {
            return Vector(0, 0, 0);
        }
        double t;
        Vector P, N;
        size_t sphere_id;
        if (intersect(ray, t, P, N, sphere_id)) {
            if (spheres[sphere_id].mirror) {
                // Reflection
                ray.reflect(P, N);
                return get_color(S, I, ray, ray_depth-1);
            } else if (spheres[sphere_id].transparent) {
                // Refraction
                double n1 = spheres[sphere_id].n_out;
                double n2 = spheres[sphere_id].n_in;
                Ray refracted_ray = refract(ray, P, N, n1, n2, invert);
                return get_color(S, I, refracted_ray, ray_depth-1, invert);
            } else {
                // Diffusion
                Vector L0 = get_intensity(P, spheres[sphere_id].albedo, N, I, S);
                // Indirect lighting
                Vector V = random_cos(N);
                Ray random_ray = Ray(P + 0.001 * N);
                random_ray.set_direction(V);
                L0 = L0 + spheres[sphere_id].albedo * get_color(S, I, random_ray, ray_depth-1);
                return L0;
            }
        }
        return Vector(0, 0, 0);
    }

    Vector get_color(Sphere &L, Ray &ray, int ray_depth=1, bool invert=false, bool diffused=false) {
        if (ray_depth < 0) {
            return Vector(0, 0, 0);
        }
        double t;
        Vector P, N;
        size_t sphere_id;
        if (intersect(ray, t, P, N, sphere_id)) {
            if (spheres[sphere_id].light) {
                if (diffused) {
                    return Vector(0, 0, 0);
                } else {
                    return Vector(1, 1, 1) * L.i / (4 * std::pow(M_PI * L.R, 2));
                }
            } else if (spheres[sphere_id].mirror) {
                // Reflection
                ray.reflect(P, N);
                return get_color(L, ray, ray_depth-1, invert);
            } else if (spheres[sphere_id].transparent) {
                // Refraction
                double n1 = spheres[sphere_id].n_out;
                double n2 = spheres[sphere_id].n_in;
                Ray refracted_ray = refract(ray, P, N, n1, n2, invert);
                return get_color(L, refracted_ray, ray_depth-1, invert);
            } else {
                // Direct lighting
                Vector D = (P - L.C) / (P - L.C).norm();
                Vector xprime = L.R * random_cos(D) + L.C;
                Vector Nprime = (xprime - L.C) / (xprime - L.C).norm();
                Vector omega_i = (xprime - P) / (xprime - P).norm();
                double pdf = std::max(dot(Nprime, D), 0.) / (M_PI * L.R * L.R);
                Vector x = P + 0.001 * N;
                Vector S = xprime + 0.001 * Nprime;
                int v = visibility(x, S);
                Vector L0 = L.i / (4 * std::pow(M_PI * L.R, 2)) * (spheres[sphere_id].albedo / M_PI) * v;
                L0 = L0 * std::max(dot(N, omega_i), 0.) * std::max(dot(Nprime, -omega_i), 0.) / ((xprime - P).norm2() * pdf);
                // Indirect lighting
                Vector V = random_cos(N);
                Ray random_ray = Ray(P + 0.001 * N);
                random_ray.set_direction(V);
                L0 = L0 + spheres[sphere_id].albedo * get_color(L, random_ray, ray_depth-1, invert, true);
                return L0;
            }
        }
        return Vector(0, 0, 0);
    }
};
