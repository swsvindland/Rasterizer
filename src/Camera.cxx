// Sam Svindland
//
// Verified to run on MacOS 10.14.2
//                    Apple LLVM version 10.0.0 (clang-1000.11.45.5)
//                    cmake 3.13.2
//                    VTK 8.1.2
#include "Camera.h"

double norm(double v0, double v1, double v2) {
    return sqrt(v0 * v0 + v1 * v1 + v2 * v2);
}

Matrix Camera::ViewTransform(void) {
    Matrix vt;
    vt.A[0][0] = 1 / tan(angle/2);
    vt.A[1][1] = 1 / tan(angle/2);
    vt.A[2][2] = (far+near)/(far-near);
    vt.A[2][3] = -1;
    vt.A[3][2] = (2*far*near) / (far - near);
    return vt;
}

Matrix Camera::CameraTransform(void) {
    double O[3] = {
        position[0],
        position[1],
        position[2]
    };
    double t0 = norm(O[0] - focus[0], O[1] - focus[1], O[2] - focus[2]);
    double w[3] = {
        (O[0] - focus[0]) / t0, 
        (O[1] - focus[1]) / t0, 
        (O[2] - focus[2]) / t0
    };
    double t1 = norm((up[1] * (O[2] - focus[2])) - (up[2] * (O[1] - focus[1])), (up[2] * (O[0] - focus[0])) - (up[0] * (O[2] - focus[2])), (up[0] * (O[1] - focus[1])) - (up[1] * (O[0] - focus[0])));
    double u[3] = {
        ((up[1] * (O[2] - focus[2])) - (up[2] * (O[1] - focus[1]))) / t1, 
        ((up[2] * (O[0] - focus[0])) - (up[0] * (O[2] - focus[2]))) / t1, 
        ((up[0] * (O[1] - focus[1])) - (up[1] * (O[0] - focus[0]))) / t1
    };
    double t2 = norm(((O[1] - focus[1]) * u[2]) - ((O[2] - focus[2]) * u[1]), ((O[2] - focus[2]) * u[0]) - ((O[0] - focus[0]) * u[2]),  ((O[0] - focus[0]) * u[1]) - ((O[1] - focus[1]) * u[0]));
    double v[3] = {
        (((O[1] - focus[1]) * u[2]) - ((O[2] - focus[2]) * u[1])) / t2, 
        (((O[2] - focus[2]) * u[0]) - ((O[0] - focus[0]) * u[2])) / t2,  
        (((O[0] - focus[0]) * u[1]) - ((O[1] - focus[1]) * u[0])) / t2

    };
    double t[3] = {
        0 - position[0],
        0 - position[1],
        0 - position[2],
    };
    //cerr << "U: " << u[0] << " " << u[1] << " " << u[2] << endl;
    //cerr << "V: " << v[0] << " " << v[1] << " " << v[2] << endl;
    //cerr << "W: " << w[0] << " " << w[1] << " " << w[2] << endl;
    //cerr << "O: " << O[0] << " " << O[1] << " " << O[2] << endl;
    //cerr << endl;

    Matrix ct;
    double ut = u[0]*t[0] + u[1]*t[1] + u[2]*t[2];
    double vt = v[0]*t[0] + v[1]*t[1] + v[2]*t[2];
    double wt = w[0]*t[0] + w[1]*t[1] + w[2]*t[2];
    double temp[4][4] = {{u[0], v[0], w[0], 0},
                       {u[1], v[1], w[1], 0},
                       {u[2], v[2], w[2], 0},
                       {ut, vt, wt, 1}};
    std::copy(&temp[0][0], &temp[0][0]+4*4, &ct.A[0][0]);
    return ct;
}

Matrix Camera::DeviceTransform(void) {
    double n = 1000;
    Matrix dt;
    dt.A[0][0] = n/2;
    dt.A[1][1] = n/2;
    dt.A[2][2] = 1;
    dt.A[3][0] = n/2;
    dt.A[3][1] = n/2;
    dt.A[3][3] = 1;
    return dt;
}

double SineParameterize(int curFrame, int nFrames, int ramp) {
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes) {
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}



