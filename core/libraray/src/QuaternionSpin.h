#pragma once
#include "PointVector.h"
//#include <Eigen/Geometry>
typedef Quaternion<float, 0> Quaternionx;

#define XACEPPTRANGE 0.02 //mm

/*
	在原始的截取的平面下，由于给出的样板并灭有在xoz平面内进行对齐，其实有在xoz平面上下都有点（可能是由于人工手动点出来的平面），
	所以所截取的平面并不是完全对齐xoz平面的，
	所以在最后平面的摆正过程中与原来的可能会出现一定偏差
*/

class QuaternionSpin
{
public:
	QuaternionSpin();
	~QuaternionSpin();

	static Vector3f VectorCompute(Vector3f a, Vector3f b, Vector3f c);
	static void ShowQuaternion(Quaternionx a);
	Vector3f ComputeQuaternion(Vector3f po);

	void SetCenterShift( Vector3f a ,int c);
	void SetTransfer(Quaternionx a) { mtransfer = a; };
	void SetHeel(float a) { mheelhight = a; };
	void SetQuaternionFromTwo(Vector3f a, Vector3f b);

	bool FindPoint(vector<Vector3f> *xoz);

	Vector3f GiveShift() { return mvshift; };

	void TransferSpin(vector<Vector3f>*in, vector<Vector3f>*out, Vector3f shift);

private:
	vector<Vector3f>*mpin, *mpout;

	vector<Quaternionx> macctransfer;

	Quaternionx mtransfer,mtransferi; //转换四元数及其逆
	
	Vector3f morigin;//morigint zero;
	Vector3f mvshift;//向原点的偏移量

	Vector3f mplowest;//xoz摆正后的触地点

	float mheelhight;

	bool InCircule(Vector3f a, vector<Vector3f> *xoz);
	bool WetherUpXZ(Vector3f cc, int *m);
	
	
};

