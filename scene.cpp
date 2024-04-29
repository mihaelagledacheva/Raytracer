#include <limits>
#include <random>
#include "geometries.cpp"

static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);

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
        //visibility
        int V = 1;
        Ray ray = Ray(P + 0.0001 * N);
        ray.set_direction(w);
        double t;
        if (intersect(ray, t)) {
            if (t <= d) {
                V = 0;
            }
        }
        return (I / (4 * M_PI * std::pow(d, 2))) * (albedo / M_PI) * V * std::max(0., dot(N, w));
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
    
    Vector get_color(Vector &S, double I, Ray &ray, int ray_depth, double n) {
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
                double u = uniform(engine);
                if (sin_theta_i <= n2 / n1 && u >= R) {
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
                Vector L0 = get_intensity(P, spheres[sphere_id].albedo, N, I, S);
                // Indirect lighting
                Vector V = random_cos(N);
                Ray random_ray = Ray(P + 0.001 * N);
                random_ray.set_direction(V);
                L0 = L0 + spheres[sphere_id].albedo * get_color(S, I, random_ray, ray_depth-1, n);
                return L0;
            }
        }
        return Vector(0, 0, 0);
    }

    Vector get_color(Sphere &L, Ray &ray, int ray_depth, double n, bool diffused) {
        if (ray_depth < 0) {
            return Vector(0, 0, 0);
        }
        double t;
        Vector P;
        Vector N;
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
                Ray reflected_ray = Ray(P + 0.0001 * N);
                Vector u_r = ray.u - 2 * dot(ray.u, N) * N;
                reflected_ray.set_direction(u_r);
                return get_color(L, reflected_ray, ray_depth-1, n, false);
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
                double u = uniform(engine);
                if (sin_theta_i <= n2 / n1 && u >= R) {
                    Vector w_T = (n1 / n2) * (ray.u - dot(ray.u, N) * N);
                    Vector w_N = -N * sqrt(1 - std::pow(n1 / n2, 2) * (1 - std::pow(dot(ray.u, N), 2)));
                    refracted_ray.set_direction(w_T + w_N);
                    return get_color(L, refracted_ray, ray_depth-1, n2, false);
                } else {
                    // Total internal reflection
                    Vector u = ray.u - 2 * dot(ray.u, N) * N;
                    refracted_ray.set_direction(u);
                    return get_color(L, refracted_ray, ray_depth-1, n1, false);
                }
            } else {
                // Direct lighting
                Vector D = (P - L.C) / (P - L.C).norm();
                Vector xprime = L.R * random_cos(D) + L.C;
                Vector Nprime = (xprime - L.C) / (xprime - L.C).norm();
                Vector omega_i = (xprime - P) / (xprime - P).norm();
                double pdf = std::max(dot(Nprime, D), 0.) / (M_PI * L.R * L.R);
                Vector S = P + 0.001 * N;
                int v = visibility(S, xprime);
                Vector L0 = L.i / (4 * std::pow(M_PI * L.R, 2)) * (spheres[sphere_id].albedo / M_PI) * v;
                L0 = L0 * std::max(dot(N, omega_i), 0.) * std::max(dot(Nprime, -omega_i), 0.) / ((xprime - P).norm2() * pdf);
                // Indirect lighting
                Vector V = random_cos(N);
                Ray random_ray = Ray(P + 0.001 * N);
                random_ray.set_direction(V);
                L0 = L0 + spheres[sphere_id].albedo * get_color(L, random_ray, ray_depth-1, n, true);
                return L0;
            }
        }
        return Vector(0, 0, 0);
    }
};
