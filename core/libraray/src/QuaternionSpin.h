#pragma once
#include "PointVector.h"
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

	Quaternionx mtransfer,mtransferi; //ת����Ԫ��������
	
	Vector3f morigin;//morigint zero;
	Vector3f mvshift;//��ԭ���ƫ����

	Vector3f mplowest;//xoz������Ĵ��ص�

	float mheelhight;

	bool InCircule(Vector3f a, vector<Vector3f> *xoz);
	bool WetherUpXZ(Vector3f cc, int *m);
	
	
};

