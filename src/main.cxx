// Sam Svindland
//
// Verified to run on MacOS 10.14.2
//                    Apple LLVM version 10.0.0 (clang-1000.11.45.5)
//                    cmake 3.13.2
//                    VTK 8.1.2

#include <iostream>
#include <math.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include "Matrix.h"
#include "Camera.h"
#include "reader.h"

using std::cerr;
using std::endl;

double max(double l []) {
    double max = 0;
    for(int i = 0; i < 3; ++i) {
        if(l[i] > max) {
            max = l[i];
        }
    }
    return max;
}
int max_index(double l []) {
    double max = 0.0;
    int index = 0;
    for(int i = 0; i < 3; ++i) {
        if(l[i] > max) {
            max = l[i];
            index = i;
        }
    }
    return index;
}

double min(double l []) {
    double min = INFINITY;
    for(int i = 0; i < 3; ++i) {
        if(l[i] < min) {
            min = l[i];
        }
    }
    return min;
}

int min_index(double l []) {
    double min = INFINITY;
    int index = 0;
    for(int i = 0; i < 3; ++i) {
        if(l[i] < min) {
            min = l[i];
            index = i;
        }
    }
    return index;
}

vtkImageData *NewImage(int width, int height) {
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, const char *filename) {
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

void renderTriangle(vtkImageData *image, float *depthBuffer, double pointY, double pointX, double leftY, double leftX, double rightY, double rightX, double pointColor[3], double leftColor[3], double rightColor[3], double pointZ, double leftZ, double rightZ, double pointShading, double leftShading, double rightShading) {
    double rowMin = ceil(fmin(fmin(pointY, leftY), rightY));
    double rowMax = floor(fmax(fmax(pointY, leftY), rightY));

    for(int r = rowMin; r <= rowMax; ++r) {
        double leftEnd;
        if(fabs(leftX - pointX) < 0.00001) {
            leftEnd = leftX;
        } else {  
            double leftM = (pointY - leftY) / (pointX - leftX);
            double leftB = leftY - (leftM * leftX);
            leftEnd = (r - leftB) / leftM;
        }
        double rightEnd;
        if(fabs(pointX - rightX) < 0.00001) {
            rightEnd = rightX;
        } else {  
            double rightM = (pointY - rightY) / (pointX - rightX);
            double rightB = rightY - (rightM * rightX);
            rightEnd = (r - rightB) / rightM;
        }

        double t0 = (r - leftY) / (pointY - leftY);
        double xLeftEnd = leftX + t0*(pointX - leftX);
        double xRightEnd = rightX + t0*(pointX - rightX);
        double zLeftEnd = leftZ + t0*(pointZ - leftZ);
        double zRightEnd = rightZ + t0*(pointZ - rightZ);
        
        double colorLeftEnd[3] = {leftColor[0] + t0*(pointColor[0] - leftColor[0]), 
                                    leftColor[1] + t0*(pointColor[1] - leftColor[1]), 
                                    leftColor[2] + t0*(pointColor[2] - leftColor[2])};
        double colorRightEnd[3] = {rightColor[0] + t0*(pointColor[0] - rightColor[0]), 
                                    rightColor[1] + t0*(pointColor[1] - rightColor[1]), 
                                    rightColor[2] + t0*(pointColor[2] - rightColor[2])};
        
        double shadingLeftEnd = leftShading + t0*(pointShading - leftShading);
        double shadingRightEnd = rightShading + t0*(pointShading - rightShading);

        for(int c = ceil(leftEnd); c <= floor(rightEnd); ++c) {
            if(r >= 0 && c >= 0 && r < 1000 && c < 1000) {
                unsigned char *buffer = (unsigned char *) image->GetScalarPointer(c,r,0);
                double t1 = (c - xLeftEnd) / (xRightEnd - xLeftEnd);
                double z = zLeftEnd + t1*(zRightEnd - zLeftEnd);
                double color[3] = {colorLeftEnd[0] + t1*(colorRightEnd[0] - colorLeftEnd[0]), 
                                    colorLeftEnd[1] + t1*(colorRightEnd[1] - colorLeftEnd[1]), 
                                    colorLeftEnd[2] + t1*(colorRightEnd[2] - colorLeftEnd[2])};
                double shading = shadingLeftEnd + t1*(shadingRightEnd - shadingLeftEnd);

                if(z > depthBuffer[(r * 1000) + c]) {
                    buffer[0] = ceil(255 * fmin(1, color[0]*shading));
                    buffer[1] = ceil(255 * fmin(1, color[1]*shading));
                    buffer[2] = ceil(255 * fmin(1, color[2]*shading));
                    depthBuffer[(r * 1000) + c] = z;
                }
            }
        }
    }
}

void managePoints(std::vector<Triangle> triangles, vtkImageData *image, float *depthBuffer) {
    for(int i = 0; i < triangles.size(); ++i) {
        double max_n = max(triangles[i].Y);
        int max_i = max_index(triangles[i].Y);
        double min_n = min(triangles[i].Y);
        int min_i = min_index(triangles[i].Y);
        
        if(fabs(min_n - triangles[i].Y[3-(max_i+min_i)]) < 0.00001) {
            int left_i = min_index(triangles[i].X);
            if(left_i != max_i) {
                renderTriangle(image, depthBuffer, max_n, triangles[i].X[max_i], triangles[i].Y[left_i], triangles[i].X[left_i], triangles[i].Y[3-(left_i+max_i)], triangles[i].X[3-(left_i+max_i)], triangles[i].colors[max_i], triangles[i].colors[left_i], triangles[i].colors[3-(left_i+max_i)], triangles[i].Z[max_i], triangles[i].Z[left_i], triangles[i].Z[3-(left_i+max_i)], triangles[i].shading[max_i], triangles[i].shading[left_i], triangles[i].shading[3-(left_i+max_i)]);
            } else {
                int right_i = max_index(triangles[i].X);
                renderTriangle(image, depthBuffer, max_n, triangles[i].X[max_i], triangles[i].Y[3-(right_i+max_i)],  triangles[i].X[3-(right_i+max_i)], triangles[i].Y[right_i], triangles[i].X[right_i], triangles[i].colors[max_i], triangles[i].colors[3-(right_i+max_i)], triangles[i].colors[right_i], triangles[i].Z[max_i], triangles[i].Z[3-(right_i+max_i)], triangles[i].Z[right_i], triangles[i].shading[max_i], triangles[i].shading[3-(right_i+max_i)], triangles[i].shading[right_i]);
            }
        } else if(fabs(max_n - triangles[i].Y[3-(max_i+min_i)]) < 0.00001) {
            int left_i = min_index(triangles[i].X);
            if(left_i != min_i) {
                renderTriangle(image, depthBuffer, min_n, triangles[i].X[min_i], triangles[i].Y[left_i], triangles[i].X[left_i], triangles[i].Y[3-(left_i+min_i)], triangles[i].X[3-(left_i+min_i)], triangles[i].colors[min_i], triangles[i].colors[left_i], triangles[i].colors[3-(left_i+min_i)], triangles[i].Z[min_i], triangles[i].Z[left_i], triangles[i].Z[3-(left_i+min_i)], triangles[i].shading[min_i], triangles[i].shading[left_i], triangles[i].shading[3-(left_i+min_i)]);
            } else {
                int right_i = max_index(triangles[i].X);
                renderTriangle(image, depthBuffer, min_n, triangles[i].X[min_i], triangles[i].Y[3-(right_i+min_i)],  triangles[i].X[3-(right_i+min_i)], triangles[i].Y[right_i], triangles[i].X[right_i], triangles[i].colors[min_i], triangles[i].colors[3-(right_i+min_i)], triangles[i].colors[right_i], triangles[i].Z[min_i], triangles[i].Z[3-(right_i+min_i)], triangles[i].Z[right_i], triangles[i].shading[min_i], triangles[i].shading[3-(right_i+min_i)], triangles[i].shading[right_i]);
            }
        } else {
            double t0 = triangles[i].X[3-(min_i+max_i)];
            double leftX, leftY, rightX, rightY;
            double leftZ = -1;
            double rightZ = -1;
            double leftColor[3];
            double rightColor[3];
            double leftShading;
            double rightShading;

            double test;
            if(fabs(triangles[i].X[max_i] - triangles[i].X[min_i]) < 0.00001) {
                test = triangles[i].X[max_i];
            } else {  
                double m = (max_n - min_n) / (triangles[i].X[max_i] - triangles[i].X[min_i]);
                double b = max_n - (m * triangles[i].X[max_i]);
                test = (triangles[i].Y[3-(min_i+max_i)] - b) / m;
            }

            if(t0 < test) {
                leftX = t0;
                leftY = triangles[i].Y[3-(min_i+max_i)];
                rightY = triangles[i].Y[3-(min_i+max_i)];

                if(fabs(triangles[i].X[max_i] - triangles[i].X[min_i]) < 0.00001) {
                    rightX = triangles[i].X[max_i];
                } else {  
                    double rightM = (max_n - min_n) / (triangles[i].X[max_i] - triangles[i].X[min_i]);
                    double rightB = max_n - (rightM * triangles[i].X[max_i]);
                    rightX = (rightY - rightB) / rightM;
                }

                double t = (rightY - triangles[i].Y[min_i]) / (triangles[i].Y[max_i] - triangles[i].Y[min_i]);
                leftZ = triangles[i].Z[3-(min_i+max_i)];
                rightZ = triangles[i].Z[min_i] + t*(triangles[i].Z[max_i] - triangles[i].Z[min_i]);

                leftColor[0] = triangles[i].colors[3-(min_i+max_i)][0];
                leftColor[1] = triangles[i].colors[3-(min_i+max_i)][1];
                leftColor[2] = triangles[i].colors[3-(min_i+max_i)][2];

                rightColor[0] = triangles[i].colors[min_i][0] + t*(triangles[i].colors[max_i][0] - triangles[i].colors[min_i][0]);
                rightColor[1] = triangles[i].colors[min_i][1] + t*(triangles[i].colors[max_i][1] - triangles[i].colors[min_i][1]);
                rightColor[2] = triangles[i].colors[min_i][2] + t*(triangles[i].colors[max_i][2] - triangles[i].colors[min_i][2]);

                leftShading = triangles[i].shading[3-(min_i+max_i)];
                rightShading = triangles[i].shading[min_i] + t*(triangles[i].shading[max_i] - triangles[i].shading[min_i]);
            } else if (t0 > test) {
                rightX = t0;
                rightY = triangles[i].Y[3-(min_i+max_i)];
                leftY = triangles[i].Y[3-(min_i+max_i)];

                if(fabs(triangles[i].X[max_i] - triangles[i].X[min_i]) < 0.00001) {
                    leftX = triangles[i].X[max_i];
                } else {  
                    double leftM = (max_n - min_n) / (triangles[i].X[max_i] - triangles[i].X[min_i]);
                    double leftB = max_n - (leftM * triangles[i].X[max_i]);
                    leftX = (leftY - leftB) / leftM;
                }

                double t = (leftY - triangles[i].Y[min_i]) / (triangles[i].Y[max_i] - triangles[i].Y[min_i]);
                rightZ = triangles[i].Z[3-(min_i+max_i)];
                leftZ = triangles[i].Z[min_i] + t*(triangles[i].Z[max_i] - triangles[i].Z[min_i]);

                rightColor[0] = triangles[i].colors[3-(min_i+max_i)][0];
                rightColor[1] = triangles[i].colors[3-(min_i+max_i)][1];
                rightColor[2] = triangles[i].colors[3-(min_i+max_i)][2];

                leftColor[0] = triangles[i].colors[min_i][0] + t*(triangles[i].colors[max_i][0] - triangles[i].colors[min_i][0]);
                leftColor[1] = triangles[i].colors[min_i][1] + t*(triangles[i].colors[max_i][1] - triangles[i].colors[min_i][1]);
                leftColor[2] = triangles[i].colors[min_i][2] + t*(triangles[i].colors[max_i][2] - triangles[i].colors[min_i][2]);

                rightShading = triangles[i].shading[3-(min_i+max_i)];
                leftShading = triangles[i].shading[min_i] + t*(triangles[i].shading[max_i] - triangles[i].shading[min_i]);
            }

            renderTriangle(image, depthBuffer, max_n, triangles[i].X[max_i], leftY, leftX, rightY, rightX, triangles[i].colors[max_i], leftColor, rightColor, triangles[i].Z[max_i], leftZ, rightZ, triangles[i].shading[max_i], leftShading, rightShading);
            renderTriangle(image, depthBuffer, min_n, triangles[i].X[min_i], leftY, leftX, rightY, rightX, triangles[i].colors[min_i], leftColor, rightColor, triangles[i].Z[min_i], leftZ, rightZ, triangles[i].shading[min_i], leftShading, rightShading);
        }
    }
}

void printMatrix(Matrix &m) {
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            cerr << m.A[i][j] << " ";
        }
        cerr << endl;
    }
    cerr << endl;
}


double calculateShading(LightingParameters lp, double *viewDirection, double *normal) {
    double ambiant = lp.Ka;

    double normL = sqrt(lp.lightDir[0] * lp.lightDir[0] + lp.lightDir[1] * lp.lightDir[1] + lp.lightDir[2] * lp.lightDir[2]);
    double normNormal = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

    double diffuse = lp.Kd * fmax(0, fabs((lp.lightDir[0]/normL) * (normal[0]/normNormal) + 
                            (lp.lightDir[1]/normL) * (normal[1]/normNormal) + 
                            (lp.lightDir[2]/normL) * (normal[2]/normNormal)));
    
    double normV = sqrt(viewDirection[0] * viewDirection[0] + viewDirection[1] * viewDirection[1] + viewDirection[2] * viewDirection[2]);
    double twoLN = 2 * ((lp.lightDir[0]) * (normal[0]) + 
                            (lp.lightDir[1]) * (normal[1]) + 
                            (lp.lightDir[2]) * (normal[2]));
    double R[3] = {twoLN * normal[0] - lp.lightDir[0], twoLN * normal[1] - lp.lightDir[1], twoLN * normal[2] - lp.lightDir[2]};    
    double normR = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    double specular = lp.Ks * pow(fmax(0, (R[0]/normR * (viewDirection[0]/normV) + 
                                            R[1]/normR * (viewDirection[1]/normV) + 
                                            R[2]/normR * (viewDirection[2]/normV))), lp.alpha);
    
    return ambiant + diffuse + specular;
}

int main() {
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1000*1000;
    for (int i = 0 ; i < npixels*3 ; i++)
    buffer[i] = 0;
    float *depthBuffer = (float *)malloc(1000*1000*sizeof(float));
    for(int i = 0; i < 1000*1000; ++i) {
        depthBuffer[i] = -1;
    }

    std::vector<Triangle> triangles = GetTriangles();

    Camera c = GetCamera(0, 1000);

    LightingParameters lp = LightingParameters();
    
    for(int i = 0; i < triangles.size(); ++i) {
        for(int j = 0; j < 3; ++j) {
            double pos[3] = {c.position[0] - triangles[i].X[j], c.position[1] - triangles[i].Y[j], c.position[2] - triangles[i].Z[j]};
            double sh = calculateShading(lp, pos, triangles[i].normals[j]);
            triangles[i].shading[j] = sh;
        }
    }

    Matrix ct = c.CameraTransform();
    Matrix vt = c.ViewTransform();
    Matrix dt = c.DeviceTransform();

    Matrix M1 = Matrix::ComposeMatrices(ct, vt);
    Matrix M2 = Matrix::ComposeMatrices(M1, dt);

    for(int i = 0; i < triangles.size(); ++i) {
        for(int j = 0; j < 3; ++j) {
            double ptIn[4] = {triangles[i].X[j], triangles[i].Y[j], triangles[i].Z[j], 1};
            double ptOut[4];
            ptOut[0] = ptIn[0]*M2.A[0][0]
                    + ptIn[1]*M2.A[1][0]
                    + ptIn[2]*M2.A[2][0]
                    + ptIn[3]*M2.A[3][0];
            ptOut[1] = ptIn[0]*M2.A[0][1]
                    + ptIn[1]*M2.A[1][1]
                    + ptIn[2]*M2.A[2][1]
                    + ptIn[3]*M2.A[3][1];
            ptOut[2] = ptIn[0]*M2.A[0][2]
                    + ptIn[1]*M2.A[1][2]
                    + ptIn[2]*M2.A[2][2]
                    + ptIn[3]*M2.A[3][2];
            ptOut[3] = ptIn[0]*M2.A[0][3]
                    + ptIn[1]*M2.A[1][3]
                    + ptIn[2]*M2.A[2][3]
                    + ptIn[3]*M2.A[3][3];
            triangles[i].X[j] = ptOut[0] / ptOut[3];
            triangles[i].Y[j] = ptOut[1] / ptOut[3];
            triangles[i].Z[j] = ptOut[2] / ptOut[3];
        }
    }
    
    managePoints(triangles, image, depthBuffer);

    char filename[32];
    sprintf(filename, "frame");
    WriteImage(image, filename);

    cerr << "Done printing" << endl;
}
