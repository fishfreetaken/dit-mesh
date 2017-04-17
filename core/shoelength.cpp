// shoelength.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
//#include <GL/glut.h>
//#include <flann/flann.h>

#include "PointVector.h"
#include "TriMesh.h"
#include "TriCesh.h"
#include "OutLine.h"
#include "ManageObj.h"
#include "QuaternionSpin.h";
#include "AddLength.h"
#include "MyOpenMesh.h"
// ----------------------------------------------------------------------------

#define stepforward 5 //Ь鸳�����5mm 
#define mlength 147
#define blength 238
#define tlength  91
#define outputname "outline.xyz"
#define openfilename "shoestl.obj"

int CAddLength() {
	/*
		������������㣺
		Pt1 28.649047  0.477465  90.453272
		Pt2 136.864537 -1.992595 157.171762
		Pt3 237.816803 1.965069  14.085081
	*/
	Vector3f p[3];
	p[0] << 28.649047, 0.477465, 90.453272;
	p[1] << 136.864537, -1.992595, 157.171762;
	p[2] << 237.816803, 1.965069, 14.085081;

	/*
		����ΪֺΧƽ�涨��
		Pt1 149.726660 -45.038730 7.791545
		Pt2 175.241229 -4.591933 34.809070
		Pt3 162.523123 36.191162 4.322560
	*/
	Vector3f p2[3];
	p2[0] << 149.726660, -45.038730, 7.791545;
	p2[1] << 175.241229, -4.591933, 34.809070;
	p2[2] << 162.523123, 36.191162, 4.322560;

	vector<Vector3f>allpoint;
	vector<Vector3f>outlinep;
	vector<Vector3i>alltri;

	ManageObj mobj(openfilename, &allpoint, &alltri);
	if (mobj.OpenObj()) {
		cout << "open OBJ file failed!" << endl;
		return 0;
	}
	mobj.ReadObj();
	OutLine ole(&allpoint);
	ole.SetCoe(p[0],p[1],p[2]);
	ole.Extract(alltri);
	ole.GiveOutlinePoints(&outlinep, 4);
	// cout << ole.CentralPoint(outlinep) << endl;

	ManageObj::OutFilePointObj(&outlinep, "outline3.obj");

	/*ManageObj mobj("outline51.obj", &allpoint);
	if (mobj.OpenObj()) {
		cout << "open OBJ file failed!" << endl;
		return 0;
	}
	mobj.ReadObj();

	AddLength addl(&allpoint,Vector3f(237.4680,4.0618,11.7566), 5);
	
	addl.SetCoe(p2);
	addl.OutlienPointSelect();
	vector<float> statis;

	vector<Vector3f> fuse;
	if (addl.IteratorDiff(&statis)) {
		cout << "out after move file!" << endl;
		addl.LinePointMove();
		fuse = addl.GiveOutlineSelect();
		ManageObj::OutFilePointObj(&fuse, "lengthM3.obj");
		cout << statis.size() << endl;
		ManageObj::OutFilePointObj(&statis, "differ3.obj");
	}*/
	return 1;
}

void Quaerniont() {
	vector<Vector3f>allpoint;
	vector<Vector3f>outpoint;
	vector<Vector3f>midpoint;
	ManageObj file("cctvt.obj", &allpoint);
	file.OpenObj();
	file.ReadObj();
	cout << allpoint.size() << endl;
	/*
		����������;
	*/
	Vector3f p[3];
	p[0] << 63.463964,  -57.425816, 42.609455;//(28.2512570,0,91.3379820);
	p[1] << 166.733872, -124.976311, 13.108240;
	p[2] << 166.188756, -68.167633, -156.102804;
	//Vector3f origin(0,0,90.06);//��ǵ��ԭ��
	QuaternionSpin a;
	a.SetQuaternionFromTwo(QuaternionSpin::VectorCompute(p[0], p[2], p[1]), Vector3f(0, 1, 0));
	a.SetCenterShift(p[0],1);//
	a.SetHeel(90.06);//���ø���
	a.TransferSpin(&allpoint,&midpoint,a.GiveShift());
	if (a.FindPoint(&midpoint)) {
		a.TransferSpin(&midpoint,&outpoint, Vector3f(28.2512570, 0, 91.3379820));
		//ManageObj::OutFilePointObj(&outpoint, "transfershift4.obj");
	}
}

#define STEPADDLENGHT 5 //mm
#define HEELHIGHT 90.06 //mm

string increname(int i) {
	string head = "cutoutline-"; string suffix = "-ext.obj";
	char num[3];
	sprintf(num, "%d", i);
	string acf(num);
	head += acf + suffix;
	return head;
}

int main(int argc, char* argv[])
{
	//Quaerniont();
	//CAddLength();

	/*
		������������㣺
		Pt1 237.816803 1.965069  14.085081
		Pt2 49.4055, -0.9980, 184.6677
		Pt3 28.649047  0.477465  90.453272
	*/

	MyMesh::Point p[3];
	MyMesh::Point p2[3];
	//p2[0] = MyMesh::Point(237.816803, 1.965069, 14.085081);
	p2[0] = MyMesh::Point(237.5670, 1.2265, 11.8482);
	p2[1] = MyMesh::Point(49.4055, -0.9980, 184.6677);
	p2[2] = MyMesh::Point(28.649047, 0.477465, 90.453272);
	
	//float coeabc = OutLine::GiveCoe(p2[0], p2[1], p2[2], &coe);  //??

	//Vector3f turnpoint(49.2967,-2.4764,184.5850); //?

	p[0] = MyMesh::Point(149.726660, -45.038730, 7.791545);
	p[1] = MyMesh::Point(175.241229, -4.591933, 34.809070);
	p[2] = MyMesh::Point(162.523123, 36.191162, 4.322560);
	
	vector<Vector3f> outline;
	MyOpenMesh ios;

	MyMesh::VertexHandle vertex[3];
	ios.ReadStlfile("shoestl.stl");
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	SurfaceCoe sfc(vertex, ios.mesh);

	SurfaceCoe *meta;
	if (sfc.Init(1)) {
		/*sfc.OutlineEigen(&outline);// ��ԭ�е�outlien�������
		ManageObj::OutFilePointObj(&outline, "outline2.obj");
		outline.clear();*/

		ios.FindNearest(p[0], p[1], p[2], vertex);		//��ʼ�����λ��
		meta = sfc.FindMetara(vertex[0], vertex[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ
		//meta->OutlineEigen(&outline);				
		//ManageObj::OutFilePointObj(&outline, "outline2zhizhi.obj");
		//outline.clear();

		meta->AllocateXCoe(stepforward);			//���ղ�ͬ��ʽ���õ���ϵ��
		//meta->OutlineExpansion();	//outline ��Χ�� private����

		//meta->OutlineEigen(&outline);
		//ManageObj::OutFilePointObj(&outline, "outline2zhizhi2.obj");
		//outline.clear();
	}
	vector<MySurCutArry> arry;
	sfc.OutCutOutline(stepforward / meta->ReturnLength(),arry);  //5������/��������
	vector<SurfaceCoe*> cutout(arry.size());
	//Vector3f axi = sfc.AxieCut(HEELHIGHT); //�����õ������ߺ����ģ�
	Vector3f axi = sfc.TempVector();//ֱ�Ӹ���������

	ManageObj arps("3.txt");
	vector<MyMesh::Point> arp;
	arps.ReadMeshPoints(arp);

	if (!axi.norm()) {
		cout << "axi erro zero!" << endl;
		return 0;
	}
	axi.normalize();
	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(axi, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
			cutout[i]->AllocateXCoe2();//��ɢ
			//cutout[i]->OutlineEigen(&outline);
			//cout << outline.size() << endl; //���outline
			//ManageObj::OutFilePointObj(&outline, increname(i).c_str());
			//outline.clear();
		}
	}
	ios.ShoeExpansion(cutout,meta);
	
	ios.WriteStlfile("shoestlext-large.stl",1);

	//�ͷ�new������������ڴ�
	for (int i = 0; i < arry.size(); i++) { 
		delete cutout[i]; 
	}
	
	/*vector<Vector3f> outline;
	MyOpenMesh ios;
	ios.ReadStlfile("shoestl.stl");

	ios.MoveVertex(3);
	ios.WriteStlfile("movelen3.stl",1);*/
	//ios.BottomVertex(&outline);
	//ManageObj::OutFilePointObj(&outline, "bottomshoedown.obj");
	system("pause");
	return 0;
}


