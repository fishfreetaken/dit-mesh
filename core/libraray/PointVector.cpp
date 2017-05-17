#include "stdafx.h"
#include "PointVector.h"
#define M_PI 3.1415926535897932384 

PointVector::PointVector(Vector3f s,Vector3f p1,Vector3f p2,Vector3f v)
{
	mp = s;
	AngleCompute(p1, p2, v);
}

PointVector::~PointVector()
{
}

/*�����µ㣬�������ӽ���Ƭ��������������ʽ����������*/
PointVector * PointVector::Insert(Vector3f s, Vector3f p1, Vector3f p2, Vector3f v) {
	if (!PointIs(s)) {
		if (next != NULL) {
			next->Insert(s, p1, p2, v);
		}
		else {
			next = new PointVector(s,p1,p2,v);
			return next;
		}
	}
	else {
		AngleCompute(p1, p2, v);
		return this;
	}
}

/*������Ӧ���ڵ����з���Ҫ��ĵ�*/
void PointVector::CheckOut(PointVector *p,vector<PointVector *> *s,Vector4f &c) {
	
	Vector3f b = p->ThisPoint();
	Vector4f a(b[0],b[1],b[2], 1);
	if ((a.transpose()*c) <= 0) {
		s->push_back(p);
		p->Modify();
	}
	if (p->next != NULL) {
		CheckOut(p->next, s, c);
	}
}

/*�ͷ�����*/
bool PointVector::DeleteChain(PointVector *p) {
	if (p->next != NULL) {
		if (DeleteChain(p->next)) {
			delete p->next;
			return true;
		}
	}
	else {
		return true;
	}
}

/*�Ƿ�Ϊ�ö���*/
bool PointVector::PointIs(Vector3f p) {
	if (mp == p) {
		return true;
	}
	return false;
}

/*���㲢�洢�н�*/
float PointVector::AngleCompute(Vector3f p1,Vector3f p2,Vector3f vec) {
	Vector3f ab;
	Vector3f ac;
	float c, angle;
	ab = p1 - mp;
	ac = p2 - mp;
	c = (ab.transpose()*ac);
	angle= (sqrt(ab.transpose()*ab)*sqrt(ac.transpose()*ac));
	c = c / angle;
	angle = acos(c);
	mangle.push_back(angle);
	mv.push_back(vec);
	return 0;
}

/*����ö����������ֵ*/
void PointVector::MeanCompute() {
	Vector3f mean(0,0,0);
	float sum=0;
	for (int i = 0; i < mangle.size(); i++) {
		sum += mangle[i];
	}
	for (int i = 0; i < mv.size(); i++) {
		mean += (mangle[i] / sum)*mv[i];
	}
	mvector=mean ;
	mvector=mvector / (sqrt(mvector.transpose()*mvector));//��һ��
}
