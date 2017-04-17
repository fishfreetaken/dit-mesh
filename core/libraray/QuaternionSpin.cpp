#include "stdafx.h"
#include "QuaternionSpin.h"


QuaternionSpin::QuaternionSpin()
{
	morigin = Vector3f(0, 0, 0);
}

QuaternionSpin::~QuaternionSpin()
{
}

/*
	平面三个点的法向量
	求垂直于两个向量的垂直向量(单位法向量)，即该向量的旋转轴(实际输入为三个三维空间坐标点)
	这里需要注意的是空间的法向量方向问题！
*/
Vector3f QuaternionSpin::VectorCompute(Vector3f a, Vector3f b, Vector3f c) {
	Vector3f ab, ac;
	ab = b - a;
	ac = c - a;

	ab = ac.cross(ab);
	//float norm;//不一定单位化  ab.norm()
	//norm = sqrt(ab.transpose()*ab);
	//ab = ab / norm;
	return ab;
}

/*
	从a转到b
*/
void QuaternionSpin::SetQuaternionFromTwo(Vector3f a,Vector3f b) {
	mtransfer= Quaternionx::FromTwoVectors(a, b);
	mtransferi = mtransfer.inverse();
	macctransfer.push_back(mtransfer);
}

Vector3f QuaternionSpin::ComputeQuaternion(Vector3f po){
	Quaternionx  out; //旋转结果坐标
	Quaternionx p(0, po[0], po[1], po[2]);//输入坐标

	out = mtransfer*p*mtransferi;
	return Vector3f(out.x(),out.y(),out.z()); //return xyz point
}

/*for debug*/
void QuaternionSpin::ShowQuaternion(Quaternionx a){
	printf("w:%f x:%f y:%f x:%f\n", a.w(), a.x(), a.y(), a.z());
}

void QuaternionSpin::TransferSpin(vector<Vector3f>*in,vector<Vector3f>*out,Vector3f shift) {
	for (auto i : *in) {
		out->push_back(ComputeQuaternion(i)+ shift);
	}
}

/*
	1:将旋转过后的坐标与原坐标相减
	0:将传入的向量直接赋值给mvshift
*/
void QuaternionSpin::SetCenterShift(Vector3f a,int i) {
	if (i) {
		mvshift = morigin - ComputeQuaternion(a);
	}
	else {
		mvshift = a;
	}
	
}

bool QuaternionSpin::InCircule(Vector3f a,vector<Vector3f> *xoz) {
	float l = a.norm();
	int count;
	if ( l >= morigin[2]) {
		Vector3f m(sqrt((l - mheelhight)*(l+mheelhight)), 0, -mheelhight);
		mtransfer = Quaternionx::FromTwoVectors(a,m);//(0,1,0)
		mtransferi = mtransfer.inverse();
		for (auto i : *xoz){
			if (!WetherUpXZ(ComputeQuaternion(i),&count)) {
				return false;
			}
		}
		if (count) {
			macctransfer.push_back(mtransfer);
			return true;
		}
	}
	return false;
}

bool QuaternionSpin::WetherUpXZ(Vector3f cc,int *m) {
	/*if (cc.z() <-(XACEPPTRANGE+mheelhight)) {
		return false;
	}
	if ((cc.z() < (XACEPPTRANGE-mheelhight)) && (cc.z() > 0)) {
		(*m)++;
	}*/
	//两个容差线都差不多，肉眼几乎不可辨别，只是与outline相差一定角度
	if (cc.z() <-( mheelhight)) {
		return false;
	}
	if ((cc.z() < (XACEPPTRANGE - mheelhight)) && (cc.z() > 0)) {
		(*m)++;
	}
	return true;
}

bool QuaternionSpin::FindPoint(vector<Vector3f> *xoz) {
	cout << "Now is fnding Point in range..." << endl;
	for (auto i : *xoz) {
		if (InCircule(i,xoz )) {
			cout << i << endl;
			cout << ComputeQuaternion(i) << endl;
			mplowest = ComputeQuaternion(i);
			return true;
		}
	}
	return false;
}

