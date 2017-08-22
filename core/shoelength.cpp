// shoelength.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
//#include <GL/glut.h>
//#include <flann/flann.h>

#include "PointVector.h"
#include "ManageObj.h"
#include "QuaternionSpin.h";
#include "MyOpenMesh.h"
// ----------------------------------------------------------------------------

#define OUTPUTFILET 1 //���Ƶ��Ե�ʱ���Ƿ�����м������ļ�

int Straighten(MyMesh::VertexHandle *vertex, MyOpenMesh&ios, SurfaceCoe* &sfc,float heel) {
	MyMesh::Point mp[3];
	mp[0] = ios.mesh.point(*vertex);
	mp[1] = ios.mesh.point(*(vertex + 1));
	mp[2] = ios.mesh.point(*(vertex + 2));
	
	Vector3f p[3],ab ,ac;
	p[0] = ios.EigenTransfer(mp[0]);
	p[1] = ios.EigenTransfer(mp[1]);
	p[2] = ios.EigenTransfer(mp[2]);

	ab = p[1] - p[0];
	ac = p[2] - p[0];

	ab = ac.cross(ab);
	vector<Vector3f>outline;
	sfc->OutlineEigen(&outline);
	if (!outline.size()) {
		return 0;
	}

	QuaternionSpin m(ab,Vector3f(0,1,0),p[2],outline, heel);  //Ĭ��p[2]Ϊԭ�� ��һ������Ǻ����
	if (!m.TransferSpin()) {
		cout << "Not Find Bottomest Point!" << endl;
		return 0;
	}
	ios.ShoeSpin(m.ReturnQuatFuse(), m.ReturnShift());//��ת

	ios.WriteStlfile("Shoe-Spin-Model.stl", 1);

	cout << "ģ�Ͱ������������ģ�ͣ�" << endl;
}

string increname(int i,string ghead) {

	string head = "-Outline"; string suffix = ".obj";
	char num[3];
	sprintf(num, "%d", i);
	string acf(num);
	head =ghead+head+acf + suffix;
	return head;
}

bool ToeExpansion(MyOpenMesh&ios, SurfaceCoe*sfc, SurfaceCoe*toe,float dist,float exp,string name)
{
	vector<Vector3f> outline;
	//cout << "���ڶ�ֺΧ����..." << endl;
	string cnoutline = name;
	string cnstl = name;
	cnoutline += ".obj";
	cnstl += ".stl";
	
	float range = 7;//ǰ�����ƽ������mm//ԭʼĬ�ϵ���9mm,��֪�����
	SurfaceCoe *toea, *toeb;
	toea = sfc->SfcMoveXLen(toe,dist+ range*2);//������ǰ��7��������
	toeb = sfc->SfcMoveXLen(toe, dist- range);//��
	toe->AllocateXCoe(exp);//����Χ������ָ��Χ��

	toe->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline, cnoutline.c_str());

	ios.ShoeExpansionWist3(toea, toe, toeb, range);

	delete toea;
	delete toeb;

#ifdef OUTPUTFILET
	//ios.WriteStlfile("Toe-ext-Large.stl", 1);
	ios.WriteStlfile(cnstl.c_str(), 1);
#endif // OUTPUTFILET
	//ios.WriteStlfile("Last-ext-Large.stl", 1);
	//cout << "ֺΧ���ν���!" << endl;
	return 1;
}

bool BackExpansion(MyOpenMesh&ios, SurfaceCoe*back, SurfaceCoe*wist, SurfaceCoe* &sfc, float exp) {
	cout << "���ڶԱ�Χ����..." << endl;
	MyMesh::Point pcc;
	int iith = sfc->UpOneInch(back->ReturnIth(2), pcc);
	MyMesh::VertexHandle start =sfc->FindNearest(pcc);
	SurfaceCoe *metc = new SurfaceCoe(back->ReturnCoe(), start, 0, ios.mesh);
	metc->Init();
	metc->SetMIth(iith);

	//back->InitMidEndPoint();
	back->AllocateXCoe(exp);
	ios.ShoeExpansionWist2(wist,back,metc);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Back-ext-Large.stl", 1);
#endif // OUTPUTFILET
	//ios.WriteStlfile("Last-ext-Large.stl", 1);
	cout << "��Χ���ν���!" << endl;
	delete metc;
	return 1;
}

bool WistExpansion(MyOpenMesh&ios,SurfaceCoe*&meta, SurfaceCoe*&back, SurfaceCoe*&wist, float exp) {
	cout << "���ڶ���Χ����..." << endl;

	wist->AllocateXCoe(exp);

	ios.ShoeExpansionWist2(meta, wist, back);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET

	cout << "��Χ���ν���!" << endl;
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, SurfaceCoe* &meta,SurfaceCoe* &sfc,float heelhight,float exp) {
	vector<Vector3f>outline; //for debug output
	cout << "���ڶ���Χ���б���..." << endl;
	int ithmet=0;

	vector<MySurCutArry> arry =sfc->OutCutOutline(exp, meta, heelhight, &ithmet);
	//vector<MySurCutArry> arry = sfc->OutCutOutline(exp, meta, heelhight, ithmet);
	ithmet--;
	if (ithmet <= 0) {
		cout << "meta line error!" << endl;
	}
	
	vector<SurfaceCoe*> cutout(arry.size());

	/*for (int i = MIDDLEFILTERMOVE + 1; i < arry.size(); i++) {
		arry[i].x += exp > 0 ? -0.0005 : 0.0005;
	}*/

	vector<MyMesh::Point> arp;
	vector<MyMesh::Point> aps;//��ȡ���ƽ������

	string ccc = "cutout.obj";
	
	//string pp = "abcdefghijklmnopqrstuvwxyz0123456";
	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitMidEndPoint(arp);//���㷨
			//cout << cutout[i]->ReturnLength() << endl;

			/*cutout[i]->OutlineEigen(&outline);
			ccc = to_string(i) + ccc;
			ManageObj::OutFilePointObj(&outline,ccc.c_str());
			ccc = "cutout.obj";*/
		}
	}
	ccc = "cutout2.obj";

	/*cutout[arry.size()-1]->InitMidTopPoint(aps,1.0/3.0);//����ǳ�ʼ��һ����ڵ�ƽ��
	cutout[arry.size() - 2]->InitMidTopPoint(aps,2.0/3.0);
	cutout[arry.size() - 3]->InitMidTopPoint(aps,1.0/4.0);
	if (aps.size() < 3) {
		cout << "aps not enough!" << endl;
		return 0;
	}*/

	//SurfacePure *uptopt = new SurfacePure(aps[0], aps[1], aps[2]);
	float exteli = 0;//���������ǰ�����չϵ������һ��

	float maxf=abs(meta->DistSurface(sfc->ReturnStartPoint())); //��ʼ�㵽��Χ����Զ����
	
	for (int i = 0; i < arry.size(); i++) {
		cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
		//if (!cutout[i]->AllocateXCoe(uptopt)) {
		//	cout << "Metara "<<i<< " Allocate Failed!" << endl;
		//	return 0;
		//}// (meta, p2[0], p2[2]);//��ɢ
		//float mmn=abs(meta->DistSurface(cutout[i]->ReturnStartPoint()));

		if (i < ithmet) {
			/*if (!cutout[i]->AllocateXCoe(2*mmn/maxf-1,i)) {
				cout << "Metara " << i << " Allocate Failed!" << endl;
				return 0;
			}*/
			
			//if (!cutout[i]->AllocateXCoe()) { //�ȼӺ��ټӿ�
			//	cout << "Metara " << i << " Allocate Failed!" << endl;
			//	return 0;
			//}
		}
		else {
			cutout[i]->AllocateXCoe(&exteli);
		}

		//cutout[i]->OutlineEigenM(&outline);
		//ccc = pp[i] + ccc;
		//ManageObj::OutFilePointObj(&outline, ccc.c_str());
		////ManageObj::OutFilePointObj(outline, ccc.c_str());
		//ccc = "cutout2.obj";
	}
	float metali = cutout[ithmet]->ReturnExtensionli();
	for (int i = 0; i < ithmet; i++) {
		float mmn = abs(meta->DistSurface(cutout[i]->ReturnStartPoint()));
		if (!cutout[i]->AllocateXCoe(2 * mmn / maxf - 1, metali*(1 - mmn / maxf))) {
			cout << "Metara " << i << " Allocate Failed!" << endl;
			return 0;
		}
	}

	ios.ShoeExpansion(cutout, sfc,meta);//�������meta�������ּ����㷨
	//ios.ShoeExpansion(cutout, sfc);  //for debug;

	//ManageObj::OutFilePointObj(&outline, "endbottom.obj");

	ios.MetaraFileter(cutout[ithmet]);//��Χ��ƽ������һ����Ҫ

#ifdef OUTPUTFILET
	ios.WriteStlfile("Metara-ext-Large.stl", 1);
#endif // OUTPUTFILET
	
	for (int i = 0; i < arry.size(); i++) {
		delete cutout[i];
	}
	cout << "��Χ���ν���!" << endl;
	return 1;
}

bool MoveLength(MyOpenMesh&ios, SurfaceCoe*meta, SurfaceCoe* &sfc, float exp) {
	//cout << "move add lenth..." << endl;
	cout << "Ь鸼ӳ�..." << endl;
	float ex=sfc->FindAddLenth(meta,exp);
	ios.ShoeAddLength(sfc->ReturnStartPoint(), meta, ex);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Len-ext-Large.stl", 1);
#endif // OUTPUTFILET
	cout << "�ӳ����ν���!" << endl;
	return 1;
}

int LastDeformation(string in, string out, MyMesh::Point *psfc,MyMesh::Point *pmeta,float heelhight,float * giveout,vector<vector<float>>&vv_toe) {
	///* ������������� */
	MyMesh::Point p2[3];
	/////*��Χ*/
	MyMesh::Point p[3];

	p2[0] = *psfc;
	p2[1] = *(psfc + 1);
	p2[2] = *(psfc + 2);

	p[0] = *pmeta;
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = *(pmeta + 1);

	vector<Vector3f> outline;
	MyOpenMesh ios;
	//ios.ReadStlfile(instlfilename.c_str());

	MyMesh::VertexHandle vertex[3], vertex2[3], vertex_toe[3];//
	//cout << "Now is Reading STL file!" << endl;
	//ios.ReadStlfile("PF-010-60-6-YX-YP.stl");
	ios.ReadStlfile(in.c_str());

	//float heelhight = 60;
	//float giveout[] = { 0,5,5,3 }; //���ȣ���Χ����Χ����Χ, ֺΧ1λ�ã�ֺΧ1���Σ�ֺΧ2λ�ã�ֺΧ2���Ρ�����
	//float toechange[] = { 179,2 };

	ios.FindNearest(p2[0], p2[1], p2[2], vertex);
	ios.FindNearest(p[0], p[1], p[2], vertex2);		//��ʼ�����λ�� p[0], p[1], p[2],

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back,*toe;

	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}
	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"outline.obj");*/

	if (0) {	//һ��Ҫ����
		if (!heelhight) {
			cout << "heel hight not ini!" << endl;
			return 0;
		}
		//cout << "ģ�Ͱ���" << endl;
		Straighten(vertex, ios, sfc, heelhight); //���һ��Ҫ��������	��
		delete sfc; //����֮����Ҫ���³�ʼ��sfc��������
		sfc = new SurfaceCoe(vertex, ios.mesh); //��������
		if (!sfc->Init(1)) {
			cout << "sfc init error!" << endl;
			return 0;
		}
	}

	//����ֺΧ����Χ����Χ����Χ��Ȼ���ٱ��Σ�

	MyMesh::VertexHandle vertex_wist[3], vertex_back[3];

	vector<float> vlentoe;  //�ų�����ǰ�ͱ��κ�����⣻
	vector<SurfaceCoe *> vsfcoe;

	vector<MyMesh::VertexHandle> vvh_toe;

	//string mg[3] = {"toe1line.obj","toe2line.obj","toe3line.obj"};
	string mg2[3] = { "toe1line2.obj","toe2line2.obj","toe3line2.obj" };
	int cc = 0;
	for (auto i : vv_toe) {
		toe = sfc->FindToeBottomPoint(i[0]);
		if (toe == NULL) {
			cout << "mis a toe line error!" << endl;
			continue;
		}
		vlentoe.push_back(toe->ReturnLength());
		toe->ReturnTriPoint(vertex_toe, ios);

		vvh_toe.push_back(vertex_toe[0]);
		vvh_toe.push_back(vertex_toe[1]);
		vvh_toe.push_back(vertex_toe[2]);

		/*toe->OutlineEigen(&outline);
		ManageObj::OutFilePointObj(&outline, mg[cc++].c_str());*/

		vsfcoe.push_back(toe);
	}

	//float flentoe[2];
	//SurfaceCoe *toe178; // 2
	//SurfaceCoe *toe190; // 1

	//SUFACECOETOE toe178, toe190;

	//if (vv_toe.size() == 1) { //ֻ����������ɣ�һ����ֻ����һ��ֺΧ�ߣ���һ������������ֺΧ��
	//	//vv_toe[0][1]==178?
	//	float f_toe1, f_toe2;
	//	f_toe1= vv_toe[0][1] == 190 ?190:17
	//	toe = sfc->FindToeBottomPoint(vv_toe[0][1]);
	//	if (toe == NULL) {
	//		cout << "mis a toe line error!" << endl;
	//	}

	//	vlentoe.push_back(toe->ReturnLength());
	//	toe->ReturnTriPoint(vertex_toe, ios);
	//}


	float vlenmeta = 0;
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ
	meta->ReturnTriPoint(vertex2, ios);//����hander�������³�ʼ����Χ
	vlenmeta = meta->ReturnLength();

	/*meta->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"metaralineA.obj");*/

	wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ

	float lenwist = wist->ReturnLength();
	float lenback = back->ReturnLength();

	wist->ReturnTriPoint(vertex_wist, ios);
	back->ReturnTriPoint(vertex_back, ios);

	if (*giveout) {
		judge = MoveLength(ios, meta, sfc, giveout[0]);
		if (!judge) {
			cout << "Main error 0" << endl;
			return -1;
		}
		delete sfc;
		delete meta;
		sfc = new SurfaceCoe(vertex, ios.mesh); //��������
		if (!sfc->Init(1)) {
			cout << "add len sfc init error!" << endl;
			return  0;
		}
		meta = new SurfaceCoe(vertex2, ios.mesh);
		meta->Init(0);
		meta->SetMIth(sfc);
		meta->ReturnLength();
		//vlenmeta = meta->ReturnLength()-vlenmeta;
		//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//���³�ʼ����Χ
	}

	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"outline-len.obj");*/

	if (*(giveout+1)) {
		//giveout[1] -= vlenmeta;//������������Ĺ�����
		judge = MetaraExpansion(ios, meta, sfc, heelhight, giveout[1]);  //��Χ�ӷ�
		if (!judge) {
			cout << "Main error 1" << endl;
			return -1;
		}
		delete meta;
		delete sfc;
		delete wist;
		delete back;
		sfc = new SurfaceCoe(vertex, ios.mesh); //��������
		if (!sfc->Init(1)) {
			cout << "Metara sfc init error!" << endl;
			return 0;
		}
		meta = new SurfaceCoe(vertex2, ios.mesh);
		meta->Init(0);
		meta->SetMIth(sfc);

		meta->OutlineEigen(&outline);
		ManageObj::OutFilePointObj(&outline, "metaralineB.obj");

		wist = new SurfaceCoe(vertex_wist, ios.mesh);
		wist->Init(0);
		wist->SetMIth(sfc);
		back = new SurfaceCoe(vertex_back, ios.mesh);
		back->Init(0);
		back->SetMIth(sfc);

		for (int j = 0; j < vsfcoe.size();j++) {
			delete vsfcoe[j];
			vertex_toe[0] = vvh_toe[j * 3 + 0];
			vertex_toe[1] = vvh_toe[j * 3 + 1];
			vertex_toe[2] = vvh_toe[j * 3 + 2];

			toe = new SurfaceCoe(vertex_toe, ios.mesh);
			toe->Init(0);
			toe->SetMIth(sfc);
			vsfcoe[j] = toe;
		}
	}

	string cctoename[] = {"Toe-ext-Large-190","Toe-ext-Large-178"};
	for (int j = 0; j < vlentoe.size();j++) {
		cout << "���ڶ�ֺΧ" << j << "����..." << endl;
		vlentoe[j] += (vv_toe[j][1]- vsfcoe[j]->ReturnLength());
		if (vlentoe[j] < 0.005) {
			continue;
		}
		ToeExpansion(ios, sfc, vsfcoe[j], vv_toe[j][0], vlentoe[j],cctoename[j]);
		cout << "ֺΧ" << j << "���ν���" << endl;
	}

	for (int j = 0; j < vv_toe.size(); j++) {
		delete vsfcoe[j];
	}

	if (*(giveout+2)) {
		lenwist = wist->ReturnLength() - lenwist;
		lenwist = giveout[2] - lenwist;
		if (abs(lenwist)>0.005) {//����һ�����ȷ�Χ
			judge = WistExpansion(ios, meta, back, wist, lenwist); //��Χ�ӷ�
			if (!judge) {
				cout << "Main error 2" << endl;
				return -1;
			}
			delete wist;
			delete meta;
			delete back;
			delete sfc;
			sfc = new SurfaceCoe(vertex, ios.mesh); //��������
			if (!sfc->Init(1)) {
				cout << "Wist sfc init error!" << endl;
				return 0;
			}
			//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ
			meta = new SurfaceCoe(vertex2, ios.mesh);
			meta->Init(0);
			meta->SetMIth(sfc);
			wist = new SurfaceCoe(vertex_wist, ios.mesh);
			wist->Init(0);
			wist->SetMIth(sfc);
			back = new SurfaceCoe(vertex_back, ios.mesh);
			back->Init(0);
			back->SetMIth(sfc);
		}
	}

	if (*(giveout+3)) {
		lenback = back->ReturnLength() - lenback;
		lenback = giveout[3] - lenback;
		if (abs(lenback)>0.005) {
			judge = BackExpansion(ios, back, wist, sfc, lenback);		 //��Χ�ӷ�
			if (!judge) {
				cout << "Main error 3" << endl;
				return -1;
			}
		}
	}

	delete sfc;
	delete meta;
	delete wist;
	delete back;

	meta = new SurfaceCoe(vertex2, ios.mesh);
	meta->Init(0);
	meta->SetMIth(sfc);

	//for (int j = 0; j < vsfcoe.size(); j++) {
	//	//delete vsfcoe[j];
	//	vertex_toe[0] = vvh_toe[j * 3 + 0];
	//	vertex_toe[1] = vvh_toe[j * 3 + 1];
	//	vertex_toe[2] = vvh_toe[j * 3 + 2];

	//	toe = new SurfaceCoe(vertex_toe, ios.mesh);
	//	toe->Init(0);
	//	toe->SetMIth(sfc);
	//	vsfcoe[j] = toe;
	//}

	ios.WriteStlfile(out.c_str(),1);

	return 0;
}

int main(int argc, char* argv[])
{
	//while (1) {
		//cout << "������TXT���������ı��ļ�·����" << endl;
		//string strs= "C:/Users/47108/Desktop/giv/901.txt";
		string strs = "C:/Users/47108/Desktop/giv/20170718/tsi.txt";
		//cin >> strs;
		
		fstream file;
		file.open(strs.c_str());

		string s;

		string inputfile, outfile;

		getline(file, inputfile);
		getline(file, outfile);

		MyMesh::Point Axial[3];
		MyMesh::Point Metara[2];
		for (int i = 0; i < 3; i++) {
			getline(file, s);
			sscanf(s.c_str(),"%f %f %f\n", &Axial[i][0], &Axial[i][1], &Axial[i][2]);
		}
		for (int i = 0; i < 2; i++) {
			getline(file, s);
			sscanf(s.c_str(), "%f %f %f\n", &Metara[i][0], &Metara[i][1], &Metara[i][2]);
		}
		float hight;
		float param[4];
		getline(file, s);
		sscanf(s.c_str(), "%f\n", &hight);

		getline(file, s);
		sscanf(s.c_str(), "%f %f %f %f\n", param, param+1, param+2, param+3);

		vector<vector<float>> toe;
		float toef[2];
		while(getline(file, s)) {
			sscanf(s.c_str(), "%f %f\n", toef,toef+1);
			vector<float> temp;
			temp.push_back(toef[0]);
			temp.push_back(toef[1]);
			toe.push_back(temp);
		}

		//fscanf()
		LastDeformation(inputfile,outfile, Axial, Metara, hight,param,toe);
		cout << "���ν�������������" << endl;
	//}
		system("pause");
	return 0;
}


