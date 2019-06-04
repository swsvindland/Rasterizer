// Sam Svindland
//
// Verified to run on MacOS 10.14.2
//                    Apple LLVM version 10.0.0 (clang-1000.11.45.5)
//                    cmake 3.13.2
//                    VTK 8.1.2
#include "Matrix.h"

#ifndef CAMERA_H
#define CAMERA_H

class Camera {
public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    Matrix ViewTransform(void);
    Matrix CameraTransform(void);
    Matrix DeviceTransform(void);
};



double SineParameterize(int curFrame, int nFrames, int ramp);
Camera GetCamera(int frame, int nframes);

#endif