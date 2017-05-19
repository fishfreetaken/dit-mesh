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

string increname(int i,string ghead) {

	string head = "-Outline"; string suffix = ".obj";
	char num[3];
	sprintf(num, "%d", i);
	string acf(num);
	head =ghead+head+acf + suffix;
	return head;
}

bool BackExpansion(MyOpenMesh&ios, MyMesh::VertexHandle *vertex, SurfaceCoe*back, SurfaceCoe*wist, SurfaceCoe* &sfc, float exp) {
	vector<Vector3f>outline;
	vector<MySurCutArry> arry = sfc->OutCutOutline(exp, wist, back,ios);
	vector<SurfaceCoe*> cutout(arry.size());
	//vector<MyMesh::Point> arp;
	//ManageObj arps("55.txt");
	//if (!arps.ReadMeshPoints(arp)) {
	//	cout << "Back Read file" << endl;
	//	return;
	//}
	//for (int i = 0; i < arry.size(); i++) {
	//	cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
	//	if (cutout[i]->Init()) {
	//		cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
	//		cutout[i]->AllocateXCoe();// (meta, p2[0], p2[2]);//��ɢ

	//		cutout[i]->OutlineEigen(&outline);
	//		ManageObj::OutFilePointObj(&outline, increname(i).c_str());
	//		outline.clear();
	//	}
	//}
	MyMesh::Point arpmid, arpend;
	for (int i = 0; i < arry.size(); i++ ) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitMidEndPoint(arpmid, arpend);
			cutout[i]->InitTwoPoints(arpmid, arpend);
			cutout[i]->AllocateXCoe();// (meta, p2[0], p2[2]);//��ɢ
		}
		
	}
	ios.ShoeExpansionWist(cutout);
	ios.WriteStlfile("back-PFHext-large3.stl", 1);
	cout << "Back stl file output over!" << endl;

	delete back;
	delete sfc;
	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "Wist sfc init error!" << endl;
		return 0;
	}
	back = sfc->FindWaistLine(back); //��ֱ��x��Ѱ����Χ
	return 1;
}


bool WistExpansion(MyOpenMesh&ios, MyMesh::VertexHandle *vertex,SurfaceCoe*meta, SurfaceCoe*back, SurfaceCoe*&wist, SurfaceCoe* &sfc, float exp) {
	vector<Vector3f>outline;
	vector<MySurCutArry> arrywist = sfc->OutCutOutline(exp, meta, wist, back);
	vector<SurfaceCoe*> cutoutwist(arrywist.size());
	//vector<MyMesh::Point> arpwist;
	/*ManageObj wistarps("44.txt");
	if (!wistarps.ReadMeshPoints(arpwist)) {
		cout << "Wist Read file!" << endl;
		return;
	}
	for (int i = 0; i < cutoutwist.size(); i++) {
		cutoutwist[i] = new SurfaceCoe(arrywist[i].n, arrywist[i].a, arrywist[i].x, ios.mesh);
		if (cutoutwist[i]->Init()) {
			cutoutwist[i]->InitTwoPoints(arpwist[i * 2], arpwist[i * 2 + 1]);
			cutoutwist[i]->AllocateXCoe();

			cutoutwist[i]->OutlineEigen(&outline);
			ManageObj::OutFilePointObj(&outline, increname(i,"Wist").c_str());
			outline.clear();
		}
	}*/
	MyMesh::Point arpmid, arpend;
	for (int i = 0; i < cutoutwist.size(); i++) {
		cutoutwist[i] = new SurfaceCoe(arrywist[i].n, arrywist[i].a, arrywist[i].x, ios.mesh);
		if (cutoutwist[i]->Init()) {
			cutoutwist[i]->InitMidEndPoint(arpmid, arpend);
			cutoutwist[i]->InitTwoPoints(arpmid, arpend);
			cutoutwist[i]->AllocateXCoe();
		}
	}
	ios.ShoeExpansionWist(cutoutwist);
	ios.WriteStlfile("wist-PFHext-ext-large3.stl", 1);

	for (int i = 0; i < arrywist.size(); i++) {
		delete cutoutwist[i];
	}
	cout << "Wist expansion over!" << endl;

	delete wist;
	delete sfc;
	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "Wist sfc init error!" << endl;
		return 0;
	}
	wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, MyMesh::VertexHandle *vertex,MyMesh::VertexHandle *vertex2, SurfaceCoe* &meta,SurfaceCoe* &sfc,float exp) {
	vector<Vector3f>outline; //for debug output
	vector<MySurCutArry> arry = sfc->OutCutOutline(exp, meta, sfc->TempVector());
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
			cutout[i]->InitMidEndPoint(arp);//���㷨
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
		return 0;
	}

	SurfacePure *uptopt = new SurfacePure(aps[0], aps[1], aps[2]);
	for (int i = 0; i < arry.size(); i++) {
		cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
		cutout[i]->AllocateXCoe(uptopt);// (meta, p2[0], p2[2]);//��ɢ
	}

	ios.ShoeExpansion(cutout);
	//ios.WriteStlfile("metara-PFHext-auto-large5.stl", 1);

	for (int i = 0; i < arry.size(); i++) {
		delete cutout[i];
	}
	delete meta;
	cout << "Metara Expansion Over!" << endl;
	delete sfc;
	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "Metara sfc init error!" << endl;
		return 0;
	}

	meta= sfc->FindMetara(vertex2[0], vertex2[2]);	//��Χ����Ҫ����Ѱ����Χ���Ѿ���ʼ��������Χ
	return 1;
}

bool MoveLength(MyOpenMesh&ios, MyMesh::VertexHandle *vertex, SurfaceCoe*meta, SurfaceCoe* &sfc, float exp) {
	cout << "move add lenth..." << endl;

	float ex=sfc->FindAddLenth(meta,exp);
	ios.ShoeAddLength(sfc->ReturnStartPoint(), meta, ex);

	ios.WriteStlfile("len-PFH-ext5.stl", 1);
	cout << "Add Length OVER!" << endl;

	delete sfc;
	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "add len sfc init error!" << endl;
		return  0;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	//Quaerniont();

	/* ������������� */
	MyMesh::Point p2[3];
	/*p2[0] = MyMesh::Point(237.5670, 1.2265, 11.8482);
	p2[1] = MyMesh::Point(49.4055, -0.9980, 184.6677);
	p2[2] = MyMesh::Point(28.649047, 0.477465, 90.453272);*/
	//PFH-10
	p2[0] = MyMesh::Point(251.237072, - 1.049236, 12.656890);
	p2[1] = MyMesh::Point(3.140964, 0.610319, 101.681542);
	p2[2] = MyMesh::Point(3.739703, - 2.403674, 10.501304);

	
	//float coeabc = OutLine::GiveCoe(p2[0], p2[1], p2[2], &coe);  //??
	/*��Χ*/
	MyMesh::Point p[3];
	/*p[0] = MyMesh::Point(149.726660, -45.038730, 7.791545);
	p[1] = MyMesh::Point(175.241229, -4.591933, 34.809070);
	p[2] = MyMesh::Point(162.523123, 36.191162, 4.322560);*/
	//PFH-10-BZ-6
	p[0] = MyMesh::Point(154.442764, -45.383193, 5.547377);
	p[1] = MyMesh::Point(0,0,0);
	p[2] = MyMesh::Point(165.761880,35.007070, 5.553459);

	vector<Vector3f> outline;
	MyOpenMesh ios;

	MyMesh::VertexHandle vertex[3], vertex2[3];
	//ios.ReadStlfile("shoestl.stl");
	ios.ReadStlfile("PFH-10-BZ-6.stl");

	/*���µ�˼·����ʹ����������ֵ����Լ���ײ�У�����ϵķ������д���*/
	//vector<MyMesh::Point> bottomline;//�ײ������ߵ���ȡ
	//vector<MyOutBottom>  botall;
	//botall=ios.ShoeBottomLine(bottomline);//��Ҫ�����ͬЬ��ƥ��ʹ�����⣬һ���潫������λ�����������յıȽ�����Ľ����1.2~1.3֮�䣬һ����������
	//ios.ShoeBottomLine(botall);
	//ManageObj::OutFilePointObj(bottomline, "bottomoutlien.obj");
	//ManageObj::OutFilePointAna(botall,"botana.obj");
	//ios.FindFloorContour(bottomline);
	//ManageObj::OutFilePointObj(bottomline, "bottom-z.obj");
	//ios.InitTriPoint(vertex);//�漰����ת���꣬���Ӧ�ò�̫����
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	ios.FindNearest(p[0], p[1], p[2], vertex2);		//��ʼ�����λ��

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back;

	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}
	//sfc.OutlineEigen(&outline);// ��ԭ�е�outlien�������
	//ManageObj::OutFilePointObj(&outline, "outline.obj");
	//outline.clear();
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ
	
	judge=MoveLength(ios, vertex, meta,sfc,5);
	if (!judge) {
		cout << "Main error 0" << endl;
		return -1;
	}

	judge=MetaraExpansion(ios,vertex,vertex2,meta,sfc,5);  //��Χ�ӷ�
	if (!judge) {
		cout << "Main error 1" << endl;
		return -1;
	}

	wist=sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ

	judge=WistExpansion(ios, vertex, meta, back, wist,sfc, 3); //��Χ�ӷ�
	if (!judge) {
		cout << "Main error 2" << endl;
		return -1;
	}

	judge=BackExpansion(ios, vertex, back, wist,sfc, 3);		 //��Χ�ӷ�
	if (!judge) {
		cout << "Main error 3" << endl;
		return -1;
	}

	//back->OutlineEigen(&outline);// ��ԭ�е�outlien�������
	//ManageObj::OutFilePointObj(&outline, "outlineback-non.obj");
	//outline.clear();

	//wist->OutlineEigen(&outline);// ��ԭ�е�outlien�������
	//ManageObj::OutFilePointObj(&outline, "outlinewist-non.obj");
	//outline.clear();

	/*meta->OutlineEigen(&outline);				
	ManageObj::OutFilePointObj(&outline, "outlinezhi-non.obj");
	outline.clear();*/

	//meta->AllocateXCoe(stepforward);			//���ղ�ͬ��ʽ���õ���ϵ��

	system("pause");
	return 0;
}


