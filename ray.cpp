#include <cmath>
#include "vector.cpp"

static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);

class Ray {
public:
    Vector O; //origin 
    Vector u; //direction 

    explicit Ray(Vector origin) {
        O = origin;
    }

    void set_direction(Vector direction) {
        u = direction / direction.norm();
    }

    void compute_direction(int x, int y, int W, int H, Vector &Q, double alpha) {
        double a = Q[0] + x + 0.5 - W/2;
        double b = Q[1] + y + 0.5 - H/2;
        double c = Q[2] - W/(2*tan(alpha/2));
        u = Vector(a, b, c);
        u.normalize();
    }

    void reflect(Vector &P, Vector &N) {
        O = P + 0.0001 * N;
        u = u - 2 * dot(u, N) * N;
    }
};

Ray refract(Ray ray, Vector &P, Vector &N, double n1, double n2, bool &invert) {
    Vector Normal = N;
    if (invert) {
        std::swap(n1, n2);
        Normal = -Normal;
    }
    double cos_theta_i = dot(ray.u, Normal);
    double sin_theta_i = sqrt(1.0 - std::pow(cos_theta_i, 2));
    double sin_theta_t = (n1 / n2) * sin_theta_i;
    double cos_theta_t = std::sqrt(1.0 - sin_theta_t * sin_theta_t);
    // Fresnel law
    double k0 = std::pow(n1 - n2, 2) / std::pow(n1 + n2, 2);
    double R = k0 + (1 - k0) * std::pow(1 - std::abs(dot(Normal, ray.u)), 5);
    double u = uniform(engine);
    if (sin_theta_t <= 1.0 && u >= R) {
        invert = !invert;
        Vector w = (n1 / n2) * (ray.u + cos_theta_i * Normal) - cos_theta_t * Normal;
        Ray refracted_ray = Ray(P + 0.001 * Normal);
        refracted_ray.set_direction(w);
        return refracted_ray;
    } else {
        // Total internal reflection
        Vector w = ray.u - 2 * dot(ray.u, N) * N;
        Ray refracted_ray = Ray(P + 0.001 * N);
        refracted_ray.set_direction(w);
        return refracted_ray;
    }
}
