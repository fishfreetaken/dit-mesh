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

#define openfilename "shoestl.obj"

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
	string head = "Metara-Outline"; string suffix = ".obj";
	char num[3];
	sprintf(num, "%d", i);
	string acf(num);
	head += acf + suffix;
	return head;
}

void BackExpansion(MyOpenMesh&ios, SurfaceCoe&sfc, SurfaceCoe*back, SurfaceCoe*wist,float exp) {
	vector<Vector3f>outline;
	vector<MySurCutArry> arry = sfc.OutCutOutline(exp, wist, back,sfc,ios);
	vector<SurfaceCoe*> cutout(arry.size());

	ManageObj arps("55.txt");
	vector<MyMesh::Point> arp;
	if (!arps.ReadMeshPoints(arp)) {
		cout << "Back Read file" << endl;
		return;
	}

	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
			cutout[i]->AllocateXCoe();// (meta, p2[0], p2[2]);//��ɢ

			cutout[i]->OutlineEigen(&outline);
			ManageObj::OutFilePointObj(&outline, increname(i).c_str());
			outline.clear();
		}
	}
	ios.ShoeExpansionWist(cutout);
	ios.WriteStlfile("back-shoestlext-3-26-gauss-large.stl", 1);
	cout << "Back stl file output over!" << endl;
}


void WistExpansion(MyOpenMesh&ios, SurfaceCoe&sfc,SurfaceCoe*meta, SurfaceCoe*back, SurfaceCoe*wist, float exp) {
	vector<Vector3f>outline;
	vector<MySurCutArry> arrywist = sfc.OutCutOutline(3, meta, wist, back);
	
	ManageObj wistarps("44.txt");
	vector<MyMesh::Point> arpwist;
	if (!wistarps.ReadMeshPoints(arpwist)) {
		cout << "Wist Read file!" << endl;
		return;
	}

	vector<SurfaceCoe*> cutoutwist(arrywist.size());
	for (int i = 0; i < cutoutwist.size(); i++) {
		cutoutwist[i] = new SurfaceCoe(arrywist[i].n, arrywist[i].a, arrywist[i].x, ios.mesh);
		if (cutoutwist[i]->Init()) {
			cutoutwist[i]->InitTwoPoints(arpwist[i * 2], arpwist[i * 2 + 1]);
			cutoutwist[i]->AllocateXCoe();

			cutoutwist[i]->OutlineEigen(&outline);
			ManageObj::OutFilePointObj(&outline, increname(i).c_str());
			outline.clear();
		}
	}
	ios.ShoeExpansionWist(cutoutwist);
	ios.WriteStlfile("wist-shoestlext-31-ext-large.stl", 1);

	for (int i = 0; i < arrywist.size(); i++) {
		delete cutoutwist[i];
	}
	cout << "Wist expansion over!" << endl;
}
bool MetaraTopRefine(vector<MyMesh::Point>&aps) {
	if (aps.size() < 3) {
		cout << "aps num is not enough!" << endl;
		return 0;
	}
	vector<MyMesh::Point>sm;
	sm.push_back(aps[0]);
	for (int i = 1; i < aps.size(); i++) {
		//cout << aps[i][1] * sm[sm.size() - 1][1] <<" "<< abs(aps[i][0] - sm[sm.size() - 1][0]) <<endl;
		if ((aps[i][1] * sm[sm.size() - 1][1]) < 0) {
			sm.push_back(aps[i]);
		}
		else if (abs(aps[i][0] - sm[sm.size() - 1][0]) > 2) {
			sm.push_back(aps[i]);
		}
		if (sm.size() >= 3) {
			break;
		}
	}
	if (sm.size() < 3) {
		cout << "sm num is not enough!" << endl;
		return 0;
	}
	aps.clear();
	for (auto i : sm) {
		aps.push_back(i);
	}
	return 1;
}

void MetaraExpansion(MyOpenMesh&ios, SurfaceCoe&sfc, SurfaceCoe*meta, float exp) {
	vector<Vector3f>outline;
	
	vector<MySurCutArry> arry = sfc.OutCutOutline(exp, meta, sfc.TempVector());
	vector<SurfaceCoe*> cutout(arry.size());

	vector<MyMesh::Point> arp;
	vector<MyMesh::Point> aps;//��ȡ���ƽ������

	/*ManageObj arps("33.txt");
	if (!arps.ReadMeshPoints(arp)) {
		cout << "Metara Read File Eroor" << endl;
		return;
	}*/
	//SurfacePure *uptop = new SurfacePure(Vector3f(125.350765, -5.539230, 162.495641), Vector3f(98.406893, 9.889621, 171.593580), Vector3f(64.520314, 0.782637, 181.256004));

	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			//cutout[i]->InitMidEndPoint(arp,123.6,aps);//184.6*2/3
			cutout[i]->InitMidEndPoint(arp);//���㷨
			//cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
			//cutout[i]->AllocateXCoe(uptop);// (meta, p2[0], p2[2]);//��ɢ

			//cutout[i]->OutlineEigenf(increname(i).c_str());

			/*cutout[i]->OutlineEigen(&outline);
			ManageObj::OutFilePointObj(&outline, increname(i).c_str());
			outline.clear();*/
		}
	}
	//ManageObj::OutFileOutlinePointXyz(&arp, "pair.xyz");//output arp for debug;
	for (int i = arry.size() - 1; i > 0; i--) { //��ʼ��������
		cutout[i]->InitMidTopPoint(aps);
		if (aps.size() >= 3) {
			break;
		}
	}
	if (aps.size() < 3) {
		cout << "aps not enough!" << endl;
		return;
	}
	/*if (!MetaraTopRefine(aps)) {
		cout << "Wan Kou Surface is not Satisfied!" << endl;
		return;
	}*/

	SurfacePure *uptopt = new SurfacePure(aps[0], aps[1], aps[2]);
	for (int i = 0; i < arry.size(); i++) {
		cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
		cutout[i]->AllocateXCoe(uptopt);// (meta, p2[0], p2[2]);//��ɢ
	}
	//ManageObj arpst("mt.txt");
	//vector<MyMesh::Point> arpp;
	//arpst.ReadMeshPoints(arpp);
	//ios.ShoeExpansion(cutout, arpp); //??debug error!

	//ios.ShoeExpansion(cutout,uptop);//just for debug
	ios.ShoeExpansion(cutout);
	ios.WriteStlfile("metara-shoestlext-auto-large3.stl", 1);

	for (int i = 0; i < arry.size(); i++) {
		delete cutout[i];
	}
	cout << "Metara Expansion Over!" << endl;
}

void MoveLength(MyOpenMesh&ios, SurfaceCoe&sfc, SurfaceCoe*meta, float exp) {
	cout << "move add lenth..." << endl;
	float ex=sfc.FindAddLenth(meta,exp);

	ios.ShoeAddLength(sfc.ReturnStartPoint(), meta, exp);
	
	ios.WriteStlfile("len-shoestl-ext.stl", 1);
	cout << "Add Length OVER!" << endl;
}

int main(int argc, char* argv[])
{
	//Quaerniont();

	/* ������������� */
	MyMesh::Point p2[3];
	//p2[0] = MyMesh::Point(237.816803, 1.965069, 14.085081);
	p2[0] = MyMesh::Point(237.5670, 1.2265, 11.8482);
	p2[1] = MyMesh::Point(49.4055, -0.9980, 184.6677);
	p2[2] = MyMesh::Point(28.649047, 0.477465, 90.453272);
	
	//float coeabc = OutLine::GiveCoe(p2[0], p2[1], p2[2], &coe);  //??
	/*��Χ*/
	MyMesh::Point p[3];
	p[0] = MyMesh::Point(149.726660, -45.038730, 7.791545);
	p[1] = MyMesh::Point(175.241229, -4.591933, 34.809070);
	p[2] = MyMesh::Point(162.523123, 36.191162, 4.322560);
	
	vector<Vector3f> outline;
	MyOpenMesh ios;

	MyMesh::VertexHandle vertex[3];
	ios.ReadStlfile("shoestl.stl");

	/*���µ�˼·����ʹ����������ֵ����Լ���ײ�У�����ϵķ������д���*/
	//vector<MyMesh::Point> bottomline;//�ײ������ߵ���ȡ
	//vector<MyOutBottom>  botall;
	//botall=ios.ShoeBottomLine(bottomline);//��Ҫ�����ͬЬ��ƥ��ʹ�����⣬һ���潫������λ�����������յıȽ�����Ľ����1.2~1.3֮�䣬һ����������
	//ios.ShoeBottomLine(botall);
	//ManageObj::OutFilePointObj(bottomline, "bottomoutlien.obj");
	//ManageObj::OutFilePointAna(botall,"botana.obj");
	//ios.FindFloorContour(bottomline);
	//ManageObj::OutFilePointObj(bottomline, "bottom-z.obj");

	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	SurfaceCoe sfc(vertex, ios.mesh); //��������

	SurfaceCoe *meta,*wist,*back;
	if (!sfc.Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}
	//sfc.OutlineEigen(&outline);// ��ԭ�е�outlien�������
	//ManageObj::OutFilePointObj(&outline, "outline.obj");
	//outline.clear();

	ios.FindNearest(p[0], p[1], p[2], vertex);		//��ʼ�����λ��
	meta = sfc.FindMetara(vertex[0], vertex[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ

	//MoveLength(ios, sfc, meta, 3);

	MetaraExpansion(ios, sfc, meta, 5);  //��Χ�ӷ�

	wist=sfc.FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	back = sfc.FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ

	WistExpansion(ios, sfc, meta, back, wist, 3); //��Χ�ӷ�

	BackExpansion(ios, sfc, back, wist, 3); //��Χ�ӷ�

	//back->OutlineEigen(&outline);// ��ԭ�е�outlien�������
	//ManageObj::OutFilePointObj(&outline, "outlineback-non.obj");
	//outline.clear();

	//wist->OutlineEigen(&outline);// ��ԭ�е�outlien�������
	//ManageObj::OutFilePointObj(&outline, "outlinewist-non.obj");
	//outline.clear();

	/*meta->OutlineEigen(&outline);				
	ManageObj::OutFilePointObj(&outline, "outlinezhi-non.obj");
	outline.clear();*/

	meta->AllocateXCoe(stepforward);			//���ղ�ͬ��ʽ���õ���ϵ��

	system("pause");
	return 0;
}


