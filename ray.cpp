#include <cmath>
#include "vector.cpp"

class Ray {
public:
    Vector O; //origin 
    Vector u; //direction 

    explicit Ray(Vector origin) {
        O = origin;
    }

    void set_direction(Vector direction) {
        u = direction;
    }

    void compute_direction(int x, int y, int W, int H, Vector Q, double alpha) {
        double a = Q[0] + x + 0.5 - W/2;
        double b = Q[1] + y + 0.5 - H/2;
        double c = Q[2] - W/(2*tan(alpha/2));
        u = Vector(a, b, c);
        u.normalize();
    }
};
