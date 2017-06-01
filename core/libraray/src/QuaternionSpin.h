#pragma once
#include "PointVector.h"
//#include "ManageObj.h"
//#include <Eigen/Geometry>
typedef Quaternion<float, 0> Quaternionx;

#define XACEPPTRANGE 0.02 //mm

/*
	��ԭʼ�Ľ�ȡ��ƽ���£����ڸ��������岢������xozƽ���ڽ��ж��룬��ʵ����xozƽ�����¶��е㣨�����������˹��ֶ��������ƽ�棩��
	��������ȡ��ƽ�沢������ȫ����xozƽ��ģ�
	���������ƽ��İ�����������ԭ���Ŀ��ܻ����һ��ƫ��
*/

class QuaternionSpin
{
public:
	//QuaternionSpin(Vector4f a); //����ת���Լ���ת�ǶȽ��й�����Ԫ��
	//QuaternionSpin(float d, float a, float b, float c);
	QuaternionSpin(Vector3f a, Vector3f b,Vector3f o, vector<Vector3f>&c, float heel);
	~QuaternionSpin();

	int QuaternionSpin::TransferSpin();

	static Vector3f VectorCompute(Vector3f a, Vector3f b, Vector3f c);
	static void ShowQuaternion(Quaternionx &a);
	Vector3f ComputeQuaternion(Vector3f po);


	Vector3f ReturnShift() { return mVShift; };
	Quaternionx ReturnXZTrans() { return mXZTransfer; }
	Quaternionx ReturnXtrans() { return mXTransfer; }
	Quaternionx ReturnQuatFuse() { return mTransfer; }
	
	void TransferSpin(vector<Vector3f>*in, vector<Vector3f>*out, Vector3f shift);

private:
	vector<Vector3f> &mOutline; //input

	Quaternionx mXZTransfer, mXZTransferi; //xozƽ��ת����Ԫ��������
	Vector3f mVShift, mOrigen;
	Quaternionx mXTransfer, mXTransferi;
	Quaternionx mTransfer;
	
	//Vector3f morigin;//morigint zero;
	//Vector3f mvshift;//��ԭ���ƫ����

	//Vector3f mplowest;//xoz������Ĵ��ص�

	float mHeelHight;
};

