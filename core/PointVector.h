#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>    
#include <Eigen/Dense>  
#include <math.h>
#include <sstream>
#include <set>
#include <map>

using namespace std;
using namespace Eigen;

typedef unsigned int uint;
/*
size : 68 Bytes
*/
class PointVector
{
public:
	PointVector(Vector3f s, Vector3f p1, Vector3f p2, Vector3f v);
	~PointVector();
	bool PointIs(Vector3f a);
	PointVector * Insert(Vector3f s, Vector3f p1, Vector3f p2, Vector3f v);

	static bool DeleteChain(PointVector *p);
	static void CheckOut(PointVector *p, vector<PointVector *> *s, Vector4f &c);

	void MeanCompute();

	Vector3f GiveOut() { return mvector; }/*给出该顶点下的向量均值*/
	Vector3f ThisPoint() { return mp; }   /*给出该顶点的坐标*/
	void NewModify(Vector3f v) { mp = v; }

	void Modify() { mupdate = true; }
	bool WhetherModify() { return mupdate; }

	

private:
	Vector3f mp; //顶点坐标
	vector<Vector3f> mv; //该顶点下相邻的6个面的法向量
	vector<float> mangle;//该顶点周围的角
	Vector3f mvector;

	PointVector * next=NULL;
	unsigned int mseq;

	bool mupdate = false;

	float AngleCompute(Vector3f p1, Vector3f p2, Vector3f vec);
};

