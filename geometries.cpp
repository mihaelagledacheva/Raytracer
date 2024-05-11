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

    virtual bool intersect(Ray &ray, double &t, Vector &P, Vector &N, Vector &albedo) {
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

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N, Vector &albedo) override {
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
                albedo = this->albedo;
                return true;
            }
        }
    }
};

class Plane {
public:
    Vector N;
    Vector A;

    Plane() {}
    Plane(Vector point, Vector normal) {
        A = point;
        N = normal;
    }

    bool intersect(Ray &ray, double &t) {
        t = dot(A - ray.O, N) / dot(ray.u, N);
        if (t < 0) {
            return false;
        }
        return true;
    }
};

class BoundingBox {
public:
    Vector Bmin, Bmax;
    Plane X0, X1;
    Plane Y0, Y1;
    Plane Z0, Z1;

    BoundingBox() {}
    BoundingBox(Vector b, Vector B) {
        Bmin = b;
        Bmax = B;
        X0 = Plane(Vector(b[0], (b[1]+B[1])/2, (b[2]+B[2])/2), Vector(1, 0, 0));
        X1 = Plane(Vector(B[0], (b[1]+B[1])/2, (b[2]+B[2])/2), Vector(1, 0, 0));
        Y0 = Plane(Vector((b[0]+B[0])/2, b[1], (b[2]+B[2])/2), Vector(0, 1, 0));
        Y1 = Plane(Vector((b[0]+B[0])/2, B[1], (b[2]+B[2])/2), Vector(0, 1, 0));
        Z0 = Plane(Vector((b[0]+B[0])/2, (b[1]+B[1])/2, b[2]), Vector(0, 0, 1));
        Z1 = Plane(Vector((b[0]+B[0])/2, (b[1]+B[1])/2, B[2]), Vector(0, 0, 1));
    }

    bool intersect(Ray &ray) {
        double tx0, ty0, tz0;
        double tx1, ty1, tz1;
        X0.intersect(ray, tx0);
        X1.intersect(ray, tx1);
        if (tx1 < tx0) {
            std::swap(tx0, tx1);
        }
        Y0.intersect(ray, ty0);
        Y1.intersect(ray, ty1);
        if (ty1 < ty0) {
            std::swap(ty0, ty1);
        }
        Z0.intersect(ray, tz0);
        Z1.intersect(ray, tz1);
        if (tz1 < tz0) {
            std::swap(tz0, tz1);
        }
        double max0 = std::max(tx0, std::max(ty0, tz0));
        double min1 = std::min(tx1, std::min(ty1, tz1));
        if (min1 > max0) {
            return true;
        } else {
            return false;
        }
    }
};


class Triangle {
public:
    Vector A, B, C;
    Vector NA, NB, NC;
    Vector UVA, UVB, UVC;

    explicit Triangle(Vector vertix_1, Vector vertix_2, Vector vertix_3) {
        A = vertix_1;
        B = vertix_2;
        C = vertix_3;
    }

    void set_normals(Vector normal_a, Vector normal_b, Vector normal_c) {
        NA = normal_a;
        NB = normal_b;
        NC = normal_c;
    }

    void set_uvs(Vector uv_a, Vector uv_b, Vector uv_c) {
        UVA = uv_a;
        UVB = uv_b;
        UVC = uv_c;
    }

    Vector barycenter() {
        double x = (A[0] + B[0] + C[0]) / 3;
        double y = (A[1] + B[1] + C[1]) / 3;
        double z = (A[2] + B[2] + C[2]) / 3;
        return Vector(x, y, z);
    }

    bool intersect(Ray &ray, double &t, double &alpha, double &beta, double &gamma) {
        Vector e1 = B - A;
        Vector e2 = C - A;
        double det = dot(e1, cross(ray.u, e2));
        beta = dot(ray.O - A, cross(ray.u, e2)) / det;
        if (beta < 0 || beta > 1) {
            return false;
        }
        gamma = dot(ray.u, cross(ray.O - A, e1)) / det;
        if (gamma < 0 || beta + gamma > 1) {
            return false;
        }
        t = dot(e2, cross(ray.O - A, e1)) / det;
        if (t < 0) {
            return false;
        }
        alpha = 1 - beta - gamma;
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
private:
    bool bounded;
    bool shading_normal;
    bool texture;
    int channels, W, H;
    unsigned char* image_data;
    
    std::vector<Triangle> triangles;
    BoundingBox box;

    std::vector<std::vector<Triangle>> triangle_sets;
    std::vector<BoundingBox> boxes;

    BoundingBox bound(const std::vector<Triangle>& ts) {
        double x_min = ts[0].A[0];
        double x_max = ts[0].A[0];
        double y_min = ts[0].A[1];
        double y_max = ts[0].A[1];
        double z_min = ts[0].A[2];
        double z_max = ts[0].A[2];
        for (size_t i = 1; i < ts.size(); ++i) {
            x_min = std::min(ts[i].A[0], std::min(ts[i].B[0], std::min(ts[i].C[0], x_min)));
            x_max = std::max(ts[i].A[0], std::max(ts[i].B[0], std::max(ts[i].C[0], x_max)));
            y_min = std::min(ts[i].A[1], std::min(ts[i].B[1], std::min(ts[i].C[1], y_min)));
            y_max = std::max(ts[i].A[1], std::max(ts[i].B[1], std::max(ts[i].C[1], y_max)));
            z_min = std::min(ts[i].A[2], std::min(ts[i].B[2], std::min(ts[i].C[2], z_min)));
            z_max = std::max(ts[i].A[2], std::max(ts[i].B[2], std::max(ts[i].C[2], z_max)));
        }
        return BoundingBox(Vector(x_min, y_min, z_min), Vector(x_max, y_max, z_max));
    }
    
    void bvh(std::vector<Triangle> ts, BoundingBox b, int threshold) {
        if (ts.size() > 0 && ts.size() <= threshold) {
            triangle_sets.push_back(ts);
            boxes.push_back(b);
        } else if (ts.size() > threshold) {
            double x = b.Bmax[0] - b.Bmin[0];
            double y = b.Bmax[1] - b.Bmin[1];
            double z = b.Bmax[2] - b.Bmin[2];
            std::vector<Triangle> set1, set2;
            if (x >= y && x >= z) {
                double mid = b.Bmin[0] + (x / 2);
                for(int i; i < ts.size(); ++i) {
                    if (ts[i].barycenter()[0] < mid) {
                        set1.push_back(ts[i]);
                    } else {
                        set2.push_back(ts[i]);
                    }
                }
            } else if (y >= x && y >= z) {
                double mid = b.Bmin[1] + (y / 2);
                for(int i; i < ts.size(); ++i) {
                    if (ts[i].barycenter()[1] < mid) {
                        set1.push_back(ts[i]);
                    } else {
                        set2.push_back(ts[i]);
                    }
                }
            } else {
                double mid = b.Bmin[2] + (z / 2);
                for(int i; i < ts.size(); ++i) {
                    if (ts[i].barycenter()[2] < mid) {
                        set1.push_back(ts[i]);
                    } else {
                        set2.push_back(ts[i]);
                    }
                }
            }
            if (set1.size() == 0 || set2.size() == 0) {
                triangle_sets.push_back(ts);
                boxes.push_back(b);
            } else {
                bvh(set1, bound(set1), threshold);
                bvh(set2, bound(set2), threshold);
            }
        }
    }

    void intersect(Triangle& triangle, Ray& ray, double& t, Vector& P, Vector& N, Vector& albedo, bool& intersection) {
        double t_geometry;
        double alpha, beta, gamma;
        if (triangle.intersect(ray, t_geometry, alpha, beta, gamma)) {
            if (t_geometry < t) {
                t = t_geometry;
                P = alpha * triangle.A + beta * triangle.B + gamma * triangle.C;
                if (shading_normal) {
                    N = alpha * triangle.NA + beta * triangle.NB + gamma * triangle.NC;
                } else {
                    N = cross(triangle.B - triangle.A, triangle.C - triangle.A);
                }
                N.normalize();
                if (texture) {
                    Vector UV = alpha * triangle.UVA + beta * triangle.UVB + gamma * triangle.UVC;
                    double u = std::fmod(UV[0], 1.0);
                    double v = std::fmod(UV[1], 1.0);
                    if (u < 0) u += 1.0;
                    if (v < 0) v += 1.0;
                    int pixel = floor((1 - v) * (H - 1)) * W + floor(u * (W - 1));
                    if (pixel < 0 || pixel >= W * H * channels) {
                        std::cout << alpha << std::endl;
                        std::cout << beta << std::endl;
                        std::cout << gamma << std::endl;
                        std::cout << triangle.UVA[0] << std::endl;
                        std::cout << triangle.UVA[1] << std::endl;
                        std::cout << pixel << std::endl;
                    }
                    double r = std::min(1., std::pow(image_data[pixel * channels + 0] / 255.0, 2.2));
                    double g = std::min(1., std::pow(image_data[pixel * channels + 1] / 255.0, 2.2));
                    double b = std::min(1., std::pow(image_data[pixel * channels + 2] / 255.0, 2.2));
                    albedo = Vector(r, g, b);
                }
            }
            intersection = true;
        }
    }

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
        for (int i = 0; i < vertices.size(); ++i) {
            center = center + vertices[i];
        }
        center = center / vertices.size();
        for (int i = 0; i < vertices.size(); ++i) {
            vertices[i] = center + factor * (vertices[i] - center);
        }
    }

    void translate(const Vector& translation) {
        for (int i = 0; i < vertices.size(); ++i) {
            vertices[i] = vertices[i] + translation;
        }
    }

    void rotate_xy(double theta) {
        for (int i = 0; i < vertices.size(); ++i) {
            double x = vertices[i][0] * cos(theta) - vertices[i][1] * sin(theta);
            double y = vertices[i][0] * sin(theta) + vertices[i][1] * cos(theta);
            vertices[i] = Vector(x, y, vertices[i][2]);
        }
        for (int i = 0; i < normals.size(); ++i) {
            double x = normals[i][0] * cos(theta) - normals[i][1] * sin(theta);
            double y = normals[i][0] * sin(theta) + normals[i][1] * cos(theta);
            normals[i] = Vector(x, y, normals[i][2]);
            normals[i].normalize();
        }
    }

    void rotate_yz(double theta) {
        for (int i = 0; i < vertices.size(); ++i) {
            double y = vertices[i][1] * cos(theta) - vertices[i][2] * sin(theta);
            double z = vertices[i][1] * sin(theta) + vertices[i][2] * cos(theta);
            vertices[i] = Vector(vertices[i][0], y, z);
        }
        for (int i = 0; i < normals.size(); ++i) {
            double y = normals[i][1] * cos(theta) - normals[i][2] * sin(theta);
            double z = normals[i][1] * sin(theta) + normals[i][2] * cos(theta);
            normals[i] = Vector(normals[i][0], y, z);
            normals[i].normalize();
        }
    }

    void set_parameters(bool shading_normal, bool texture, const char* filename="") {
        this->shading_normal = shading_normal;
        this->texture = texture;
        if (texture) {
            image_data = stbi_load(filename, &W, &H, &channels, 0);
        }
        bounded = false;
        for (int i = 0; i < indices.size(); ++i) {
            Vector A = vertices[indices[i].vtxi];
            Vector B = vertices[indices[i].vtxj];
            Vector C = vertices[indices[i].vtxk];
            Triangle triangle = Triangle(A, B, C);
            if (shading_normal) {
                triangle.set_normals(normals[indices[i].ni], normals[indices[i].nj], normals[indices[i].nk]);
            }
            if (texture) {
                triangle.set_uvs(uvs[indices[i].uvi], uvs[indices[i].uvj], uvs[indices[i].uvk]);
            }
            triangles.push_back(triangle);
        }
    }

    void bound() {
        bounded = true;
        box = bound(triangles);
        triangle_sets.push_back(triangles);
        boxes.push_back(box);
    }

    void bvh(int threshold) {
        bounded = true;
        box = bound(triangles);
        bvh(triangles, box, threshold);
    }

    bool intersect(Ray &ray, double &t, Vector &P, Vector &N, Vector &albedo) override {
        t = std::numeric_limits<double>::infinity();
        bool intersection = false;
        if (bounded) {
            for (int i = 0; i < boxes.size(); ++i) {
                if (boxes[i].intersect(ray)) {
                    for (int j = 0; j < triangle_sets[i].size(); ++j) {
                        intersect(triangle_sets[i][j], ray, t, P, N, albedo, intersection);
                    }
                }
            }
        } else {
            for (int i = 0; i < triangles.size(); ++i) {
                intersect(triangles[i], ray, t, P, N, albedo, intersection);
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
                sscanf(line, "vt %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
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
