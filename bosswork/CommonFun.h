#pragma once

#include "MyOpenMesh.h"
#include "opencv2/opencv.hpp"

void writeXYZN(MyMesh &inputMesh, vector<int> &floorContour);
void writeXYZN(vector<cv::Point3d> &surfContour, const char *fileName);

void   findFloorContour(MyMesh &inputMesh, vector<int> &floorContour);
void   rotateMesh(MyMesh &inputMesh);
void   findMainO(MyMesh &inputMesh);
void   findHull(vector<cv::Point2d> cloudMat, vector<int> &edgeL);
double findContour(MyMesh &inputMesh, cv::Vec4d plane, vector<cv::Point3d> &retContour);



//�漰������������,��͹�Ľ���OK,��͹�治����,һ�㲻��
double findContourShortBreak(MyMesh &inputMesh, cv::Vec4d plane, vector<cv::Point3d> &retContour);
