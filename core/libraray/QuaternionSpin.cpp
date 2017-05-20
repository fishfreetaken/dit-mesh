#include "stdafx.h"
#include "QuaternionSpin.h"

QuaternionSpin::QuaternionSpin(Vector3f a, Vector3f b,Vector3f o, vector<Vector3f>&c, float heel) :
	mOutline(c),
	mHeelHight(heel)
{
	mXZTransfer = Quaternionx::FromTwoVectors(a, b);
	mXZTransferi = mXZTransfer.inverse();
	//macctransfer.push_back(mXZTransfer);
	mOrigen = Vector3f(0, 0 ,mHeelHight);

	mVShift = mOrigen - ComputeQuaternion(o);//��ת���ƫ����
}

QuaternionSpin::~QuaternionSpin()
{
}

/*
	ƽ��������ķ�����
	��ֱ�����������Ĵ�ֱ����(��λ������)��������������ת��(ʵ������Ϊ������ά�ռ������)
	������Ҫע����ǿռ�ķ������������⣡
*/
Vector3f QuaternionSpin::VectorCompute(Vector3f a, Vector3f b, Vector3f c) {
	Vector3f ab, ac;
	ab = b - a;
	ac = c - a;
	ab = ac.cross(ab);
	return ab;
}

Vector3f QuaternionSpin::ComputeQuaternion(Vector3f po){
	Quaternionx  out; //��ת�������
	Quaternionx p(0, po[0], po[1], po[2]);//��������

	out = mXZTransfer*p*mXZTransferi;
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
	������ת��XOZƽ��
	Ȼ�����λ�Ƶ�origin��
	Ȼ�������תʹ�õװ��䵽x����
*/
int QuaternionSpin::TransferSpin() { 
	Quaternionx  out; //��ת�������  
	int sz = mOutline.size();
	Vector3f po;
	for (int i = 0; i < sz; i++) {
		po = mOutline[i];
		Quaternionx p(0, po[0], po[1], po[2]);//��������
		out = mXZTransfer*p*mXZTransferi;
		mOutline[i] = Vector3f(out.x(), out.y(), out.z())+mVShift;//�����ת��Ľ��
	}
	float li = 0; int ini=0;
	struct mvc {
		int x;
		int ith;
		Quaternionx tf;
		bool operator < (const struct mvc &m)const {
			return m.x > x;//����ǴӴ�С
		}
	};
	set<struct mvc> mis; 
	
	for (int i = 0; i < sz; i++) {
		vector<Vector3f>ft=mOutline;
		po = mOutline[i];
		li = (po - mOrigen).norm();
		if (li < mHeelHight) {
			continue;
		}
		li = sqrt((li - mHeelHight)*(li + mHeelHight));
		Quaternionx transfer = Quaternionx::FromTwoVectors((po-mOrigen),Vector3f(li,0,-mHeelHight));//(0,1,0)
		Quaternionx transferi = transfer.inverse();
		/*Quaternionx sms(0,mOutline[i][0],mOutline[i][1],mOutline[i][2]);
		out = transfer*sms*transferi;*/
		int cout = 0;
		for (int j = 0; j < sz; j++) {
			Quaternionx mm(0, ft[j][0],ft[j][1],ft[j][2]);
			out = transfer*mm*transferi;
			if (out.z() < 0) {  //-1:1322 -0.5:1330
				cout++;
				//break;
			}
			//ini++;
		}
		struct mvc git;
		git.ith = i;
		git.x = cout;
		git.tf = transfer;
		mis.insert(git);
		//if (ini) {
		//	mXTransfer = transfer;
		//	mXTransferi = transferi;
		//	//cout <<" Bottomest Point is :"<<mOutline[i]<< endl;  //for debug;
		//	mTransfer = mXTransfer*mXZTransfer;
		//	return 1;
		//}
	}
	set<struct mvc>::iterator it(mis.begin());
	mXTransfer = it->tf;
	mXTransferi = mXTransfer.inverse();
	mTransfer = mXTransfer*mXZTransfer;
	return 1;
}
//
//bool QuaternionSpin::InCircule(Vector3f a,vector<Vector3f> *xoz) {
//	float l = a.norm();
//	int count;
//	if ( l >= mOrigen[2]) {
//		Vector3f m(sqrt((l - mHeelHight)*(l+ mHeelHight)), 0, -mHeelHight);
//		mtransfer = Quaternionx::FromTwoVectors(a,m);//(0,1,0)
//		mtransferi = mtransfer.inverse();
//		for (auto i : *xoz){
//			if (!WetherUpXZ(ComputeQuaternion(i),&count)) {
//				return false;
//			}
//		}
//		if (count) {
//			macctransfer.push_back(mtransfer);
//			return true;
//		}
//	}
//	return false;
//}
//
//bool QuaternionSpin::WetherUpXZ(Vector3f cc,int *m) {
//	/*if (cc.z() <-(XACEPPTRANGE+mheelhight)) {
//		return false;
//	}
//	if ((cc.z() < (XACEPPTRANGE-mheelhight)) && (cc.z() > 0)) {
//		(*m)++;
//	}*/
//	//�����ݲ��߶���࣬���ۼ������ɱ��ֻ����outline���һ���Ƕ�
//	if (cc.z() <-( mheelhight)) {
//		return false;
//	}
//	if ((cc.z() < (XACEPPTRANGE - mheelhight)) && (cc.z() > 0)) {
//		(*m)++;
//	}
//	return true;
//}

//bool QuaternionSpin::FindPoint(vector<Vector3f> *xoz) {
//	cout << "Now is fnding Point in range..." << endl;
//	for (auto i : *xoz) {
//		if (InCircule(i,xoz )) {
//			cout << i << endl;
//			cout << ComputeQuaternion(i) << endl;
//			mplowest = ComputeQuaternion(i);
//			return true;
//		}
//	}
//	return false;
//}

