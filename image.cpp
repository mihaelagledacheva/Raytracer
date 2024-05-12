#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include "vector.cpp"

void lab_5_1() {
    int W, H, c;
	unsigned char *image = stbi_load("images/image_1.jpg", &W, &H, &c, 0);

    int Wprime, Hprime;
	unsigned char *target = stbi_load("images/image_2.jpg", &Wprime, &Hprime, &c, 0);
	
    std::vector<double> image_double(W*H*3);
	for (int i = 0; i < W*H*3; ++i) {
		image_double[i] = image[i];
    }
    
    for (int reductions = 0; reductions < W - Wprime; ++reductions) {
        
        std::vector<double> I(W*H);
        for (int i = 0; i < W*H; ++i) {
            I[i] = image_double[3*i] + image_double[3*i+1] + image_double[3*i+2];
        }

        std::vector<double> E(W*H);
        for (int j = 0; j < H; ++j) {
            for (int i = 0; i < W; ++i) {
                E[j*W+i] = std::abs((i < W ? I[j*W+i+1] : 0) - (i > 0 ? I[j*W+i-1] : 0)) + 
                           std::abs((j < H ? I[(j+1)*W+i] : 0) - (j > 0 ? I[(j-1)*W+i] : 0));
            }
        }

        std::vector<double> C(W*H);
        for (int j = 0; j < H; ++j) {
            for (int i = 0; i < W; ++i) {
                if (j == 0) {
                    C[j*W+i] = E[j*W+i];
                } else {
                    if (i == 0) {
                        C[j*W+i] = std::min(C[(j-1)*W+i], C[(j-1)*W+i+1]) + E[j*W+i];
                    } else if (i == W-1) {
                        C[j*W+i] = std::min(C[(j-1)*W+i-1], C[(j-1)*W+i]) + E[j*W+i];
                    } else {
                        C[j*W+i] = std::min(C[(j-1)*W+i-1], std::min(C[(j-1)*W+i], C[(j-1)*W+i+1])) + E[j*W+i];
                    }
                }
            }
        }

        int seam = 0;
        double min_pixel = C[(H-1)*W];
        for (int i = 1; i < W; ++i) {
            if (C[(H-1)*W+i] < min_pixel) {
                min_pixel = C[(H-1)*W+i];
                seam = i;
            }
        }

        for (int i = 0; i < W-1; ++i) {
            if (i >= seam) {
                for (int k = 0; k < 3; ++k) {
                    image_double[((H-1)*W+i)*3+k] = image_double[((H-1)*W+i+1)*3+k];
                }
            }
        }

        for (int j = H-2; j >= 0; --j) {
            if (seam == 0) {
                double C1 = C[j*W];
                double C2 = C[j*W+1];
                if (C1 > C2) {
                    seam = 1;
                }
            } else if (seam == W-1) {
                double C1 = C[j*W+W-2];
                double C2 = C[j*W+W-1];
                if (C1 < C2) {
                    seam = W-2;
                }
            } else {
                double C1 = C[j*W+seam-1];
                double C2 = C[j*W+seam];
                double C3 = C[j*W+seam+1];
                if (C1 < C2 && C1 < C3) {
                    seam = seam - 1;
                } else if(C3 < C1 && C3 < C2) {
                    seam = seam + 1;
                } 
            }
            for (int i = 0; i < W-1; ++i) {
                if (i >= seam) {
                    for (int k = 0; k < 3; ++k) {
                        image_double[(j*W+i)*3+k] = image_double[(j*W+i+1)*3+k];
                    }
                }
            }
        }
    }

    for (int reductions = 0; reductions < H - Hprime; ++reductions) {
        
        std::vector<double> I(W*H);
        for (int i = 0; i < W*H; ++i) {
            I[i] = image_double[3*i] + image_double[3*i+1] + image_double[3*i+2];
        }

        std::vector<double> E(W*H);
        for (int j = 0; j < H; ++j) {
            for (int i = 0; i < W; ++i) {
                E[j*W+i] = std::abs((i < W ? I[j*W+i+1] : 0) - (i > 0 ? I[j*W+i-1] : 0)) + 
                           std::abs((j < H ? I[(j+1)*W+i] : 0) - (j > 0 ? I[(j-1)*W+i] : 0));
            }
        }

        std::vector<double> C(W*H);
        for (int j = 0; j < H; ++j) {
            for (int i = 0; i < W; ++i) {
                if (i == 0) {
                    C[j*W+i] = E[j*W+i];
                } else {
                    if (j == 0) {
                        C[j*W+i] = std::min(C[j*W+i-1], C[(j+1)*W+i-1]) + E[j*W+i];
                    } else if (i == W-1) {
                        C[j*W+i] = std::min(C[(j-1)*W+i-1], C[j*W+i-1]) + E[j*W+i];
                    } else {
                        C[j*W+i] = std::min(C[(j-1)*W+i-1], std::min(C[j*W+i-1], C[(j+1)*W+i-1])) + E[j*W+i];
                    }
                }
            }
        }

        int seam = 0;
        double min_pixel = C[W-1];
        for (int j = 1; j < H; ++j) {
            if (C[j*W+W-1] < min_pixel) {
                min_pixel = C[j*W+W-1];
                seam = j;
            }
        }

        for (int j = 0; j < H-1; ++j) {
            if (j >= seam) {
                for (int k = 0; k < 3; ++k) {
                    image_double[(j*W+W-1)*3+k] = image_double[((j+1)*W+W-1)*3+k];
                }
            }
        }

        for (int i = W-2; i >= 0; --i) {
            if (seam == 0) {
                double C1 = C[i];
                double C2 = C[W+i];
                if (C1 > C2) {
                    seam = 1;
                }
            } else if (seam == H-1) {
                double C1 = C[(H-2)*W+i];
                double C2 = C[(H-1)*W+i];
                if (C1 < C2) {
                    seam = H-2;
                }
            } else {
                double C1 = C[(seam-1)*W+i];
                double C2 = C[seam*W+i];
                double C3 = C[(seam+1)*W+i];
                if (C1 < C2 && C1 < C3) {
                    seam = seam - 1;
                } else if(C3 < C1 && C3 < C2) {
                    seam = seam + 1;
                } 
            }
            for (int j = 0; j < H-1; ++j) {
                if (j >= seam) {
                    for (int k = 0; k < 3; ++k) {
                        image_double[(j*W+i)*3+k] = image_double[((j+1)*W+i)*3+k];
                    }
                }
            }
        }
    }

	std::vector<unsigned char> image_result(Wprime*Hprime*3, 0);
	for (int j = 0; j < Hprime; ++j) {
        for (int i = 0; i < Wprime; ++i) {
            for (int k = 0; k < 3; ++k) {
                image_result[(j*Wprime+i)*3+k] = image_double[(j*W+i)*3+k];
            } 
        }
    }
	stbi_write_png("_results/lab_5_1.png", Wprime, Hprime, 3, &image_result[0], 0);
}

void lab_5_2(int iterations=100) {
    int W, H, c;
	unsigned char *image = stbi_load("images/image_2.jpg", &W, &H, &c, 0);
	unsigned char *target = stbi_load("_results/lab_5_1.png", &W, &H, &c, 0);
	
    std::vector<double> image_double(W*H*3);
    std::vector<double> target_double(W*H*3);
	for (int i = 0; i < W*H*3; ++i) {
		image_double[i] = image[i];
        target_double[i] = target[i];
    }

    for (int iteration = 0; iteration < iterations; ++iteration) {
        static std::default_random_engine engine(10);
        static std::uniform_real_distribution<double> uniform(0, 1);

        double r1 = uniform(engine);
        double r2 = uniform(engine);
        double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
        double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
        double z = 1 - 2 * r2;

        Vector direction = Vector(x, y, z);

        std::vector<std::pair<double, int>> projI(W*H);
        std::vector<std::pair<double, int>> projM(W*H);

        for (int i = 0; i < W*H; ++i) {
            int index = i * 3;
            projI[i] = {dot(Vector(image_double[3*i], image_double[3*i+1], image_double[3*i+2]), direction), i};
            projM[i] = {dot(Vector(target_double[3*i], target_double[3*i+1], target_double[3*i+2]), direction), i};
        }

        std::sort(projI.begin(), projI.end(), [](auto& a, auto& b) {return a.first < b.first;});
        std::sort(projM.begin(), projM.end(), [](auto& a, auto& b) {return a.first < b.first;});

        for (int i = 0; i < W*H; ++i) {
            Vector advect = (projM[i].first - projI[i].first) * direction;
            for (int k = 0; k < 3; ++k) {
                image_double[3*projI[i].second + k] += advect[k];
            }
        }
    }

	std::vector<unsigned char> image_result(W*H*3, 0);
	for (int j = 0; j < H; ++j) {
        for (int i = 0; i < W; ++i) {
            for (int k = 0; k < 3; ++k) {
                image_result[(j*W+i)*3+k] = image_double[(j*W+i)*3+k];
            } 
        }
    }
	stbi_write_png("_results/lab_5_2.png", W, H, 3, &image_result[0], 0);
}

int main() {
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::seconds::rep duration;
    
    /*std::cout << "Image Retargeting" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_5_1();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_5_1.png" << std::endl << std::endl;*/

    std::cout << "Color Matching" << std::endl;
    start = std::chrono::steady_clock::now();
    lab_5_2();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Time elapsed: " << duration << " seconds" << std::endl;
    std::cout << "Output stored in _results/lab_5_2.png" << std::endl << std::endl;
}
