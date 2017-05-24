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



//涉及到轮廓点排序,对凸的截面OK,非凸面不可以,一般不用
double findContourShortBreak(MyMesh &inputMesh, cv::Vec4d plane, vector<cv::Point3d> &retContour);
