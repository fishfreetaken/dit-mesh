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

	Vector3f GiveOut() { return mvector; }/*�����ö����µ�������ֵ*/
	Vector3f ThisPoint() { return mp; }   /*�����ö��������*/
	void NewModify(Vector3f v) { mp = v; }

	void Modify() { mupdate = true; }
	bool WhetherModify() { return mupdate; }

	

private:
	Vector3f mp; //��������
	vector<Vector3f> mv; //�ö��������ڵ�6����ķ�����
	vector<float> mangle;//�ö�����Χ�Ľ�
	Vector3f mvector;

	PointVector * next=NULL;
	unsigned int mseq;

	bool mupdate = false;

	float AngleCompute(Vector3f p1, Vector3f p2, Vector3f vec);
};

