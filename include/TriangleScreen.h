// Sam Svindland
//
// Verified to run on MacOS 10.14.2
//                    Apple LLVM version 10.0.0 (clang-1000.11.45.5)
//                    cmake 3.13.2
//                    VTK 8.1.2

#ifndef TRIANGLE_SCREEN_H
#define TRIANGLE_SCREEN_H

class Triangle {
public:
    double        X[3];
    double        Y[3];
    double        Z[3];
    double colors[3][3];
    double normals[3][3];
    double shading[3];
};

class Screen {
public:
    unsigned char *buffer;
    int width, height;
};

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

#endif