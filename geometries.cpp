#include <string>
#include "ray.cpp"
#include <iostream>
class Geometry {
public:
    Vector albedo;
    bool light;
    double i;
    bool mirror;
    bool transparent;
    double n_in;
    double n_out;
    
    Geometry() {}

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

    virtual bool intersect(Ray &ray, double &t, Vector &P, Vector &N) {
        return false;
    }
};


class Sphere : public Geometry {
public:
    Vector C; //center
    double R; //radius

    explicit Sphere(Vector center, double radius, Vector color) {
        C = center;
        R = radius;
        albedo = color;
        light = false;
        mirror = false;
        transparent = false;
    }

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N) override {
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
                P = ray.O + t * ray.u;
                N = (P - C) / (P - C).norm();
                return true;
            }
        }
    }
};


class Triangle : public Geometry {
public:
    Vector A;
    Vector B;
    Vector C;

    explicit Triangle(Vector vertix_1, Vector vertix_2, Vector vertix_3) {
        A = vertix_1;
        B = vertix_2;
        C = vertix_3;
        albedo = Vector(1, 1, 1);
        light = false;
        mirror = false;
        transparent = false;
    }

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N) override {
        Vector e1 = B - A;
        Vector e2 = C - A;
        double det = dot(e1, cross(ray.u, e2));
        double u = dot(ray.O - A, cross(ray.u, e2)) / det;
        if (u < 0 || u > 1) {
            return false;
        }
        double v = dot(ray.u, cross(ray.O - A, e1)) / det;
        if (v < 0 || u + v > 1) {
            return false;
        }
        t = dot(e2, cross(ray.O - A, e1)) / det;
        if (t < 0) {
            return false;
        }
        P = (1 - u - v) * A + u * B + v * C;
        N = cross(e1, e2) / cross(e1, e2).norm();
        return true;
    }
};


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, 
                    int ni = -1, int nj = -1, int nk = -1, 
                    int uvi = -1, int uvj = -1, int uvk = -1, 
                    int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), 
                                                          uvi(uvi), uvj(uvj), uvk(uvk), 
                                                          ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};


class TriangleMesh : public Geometry {
public:
    TriangleMesh() {
        albedo = Vector(1, 1, 1);
        light = false;
        mirror = false;
        transparent = false;
    };
    ~TriangleMesh() {};

    void scale(double factor) {
        Vector center(0, 0, 0);
        for (const auto& vertex : vertices) {
            center = center + vertex;
        }
        center = center / vertices.size();
        for (auto& vertex : vertices) {
            vertex = center + factor * (vertex - center);
        }
    }

    void translate(const Vector& translation) {
        for (auto& vertex : vertices) {
            vertex = vertex + translation;
        }
    }

    void rotate_xy(double theta) {
        for (auto& vertex : vertices) {
            double x = vertex[0] * cos(theta) - vertex[1] * sin(theta);
            double y = vertex[0] * sin(theta) + vertex[1] * cos(theta);
            vertex = Vector(x, y, vertex[2]);
        }
    }

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N) override {
        t = std::numeric_limits<double>::infinity();
        bool intersection = false;
        for (auto index : indices) {
            Vector A = vertices[index.vtxi];
            Vector B = vertices[index.vtxj];
            Vector C = vertices[index.vtxk];
            Triangle triangle = Triangle(A, B, C);
            double t_geometry;
            Vector P_geometry;
            Vector N_geometry;
            if (triangle.intersect(ray, t_geometry, P_geometry, N_geometry)) {
                if (t_geometry < t) {
                    t = t_geometry;
                    P = P_geometry;
                    N = N_geometry;
                }
                intersection = true;
            }
        }
        return intersection;
    }
    
    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
};
