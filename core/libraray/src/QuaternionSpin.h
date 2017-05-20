#pragma once
#include "PointVector.h"
//#include "ManageObj.h"
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
	//QuaternionSpin(Vector4f a); //由旋转轴以及旋转角度进行构造四元数
	//QuaternionSpin(float d, float a, float b, float c);
	QuaternionSpin(Vector3f a, Vector3f b,Vector3f o, vector<Vector3f>&c, float heel);
	~QuaternionSpin();

	int QuaternionSpin::TransferSpin();

	static Vector3f VectorCompute(Vector3f a, Vector3f b, Vector3f c);
	static void ShowQuaternion(Quaternionx a);
	Vector3f ComputeQuaternion(Vector3f po);


	Vector3f ReturnShift() { return mVShift; };
	Quaternionx ReturnXZTrans() { return mXZTransfer; }
	Quaternionx ReturnXtrans() { return mXTransfer; }
	Quaternionx ReturnQuatFuse() { return mTransfer; }
	
	void TransferSpin(vector<Vector3f>*in, vector<Vector3f>*out, Vector3f shift);

private:
	vector<Vector3f> &mOutline; //input

	Quaternionx mXZTransfer, mXZTransferi; //xoz平面转换四元数及其逆
	Vector3f mVShift, mOrigen;
	Quaternionx mXTransfer, mXTransferi;
	Quaternionx mTransfer;
	
	//Vector3f morigin;//morigint zero;
	//Vector3f mvshift;//向原点的偏移量

	//Vector3f mplowest;//xoz摆正后的触地点

	float mHeelHight;
};

