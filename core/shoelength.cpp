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


	/*int pos =sf.find_last_of('.');
	string outf;
	if (pos < 0) {
		outf += "Shoe-Spin-Model.stl";
	}
	else {
		outf += sf.substr(0, pos);
		outf += "-spin.stl";
	}
	ios.WriteStlfile(outf.c_str(), 1);*/
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

bool ToeExpansion(MyOpenMesh&ios, SurfaceCoe*sfc, SurfaceCoe*toe,float dist,float exp)
{
	//vector<Vector3f> outline;
	cout << "���ڶ�ֺΧ����..." << endl;
	
	float range = 9;
	SurfaceCoe *toea, *toeb;
	toea = sfc->SfcMoveXLen(toe,dist+ range);//������ǰ��7��������
	toeb = sfc->SfcMoveXLen(toe, dist- range);
	toe->AllocateXCoe(exp);//����Χ������ָ��Χ��
	ios.ShoeExpansionWist3(toea, toe, toeb, range);

	/*toe->OutlineEigenM(&outline);
	ManageObj::OutFilePointObj(&outline, "toeaf-outline.obj");
	outline.clear();

	toeb->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline, "toeb-outline.obj");
	outline.clear();

	toea->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline, "toea-outline.obj");
	outline.clear();*/

	delete toea;
	delete toeb;

#ifdef OUTPUTFILET
	ios.WriteStlfile("Toe-ext-Large.stl", 1);
#endif // OUTPUTFILET
	//ios.WriteStlfile("Last-ext-Large.stl", 1);
	cout << "ֺΧ���ν���!" << endl;
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

	/*vector<Vector3f>outline;
	wist->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"wist-exp-large.obj");
	outline.clear();*/

	wist->AllocateXCoe(exp);

	//wist->OutlineEigenM(&outline);
	//ManageObj::OutFilePointObj(&outline, "wist-exp-M-large.obj");

	ios.ShoeExpansionWist2(meta, wist, back);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET

	cout << "��Χ���ν���!" << endl;
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, SurfaceCoe* &meta,SurfaceCoe* &sfc,float heelhight,float exp) {
	//vector<Vector3f>outline; //for debug output
	vector<MySurCutArry> arry = sfc->OutCutOutline(exp, meta, heelhight);
	vector<SurfaceCoe*> cutout(arry.size());

	vector<MyMesh::Point> arp;
	vector<MyMesh::Point> aps;//��ȡ���ƽ������

	/*ManageObj arps("33.txt");
	if (!arps.ReadMeshPoints(arp)) {
		cout << "Metara Read File Eroor" << endl;
		return;
	}*/
	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitMidEndPoint(arp);//���㷨
		}
	}
	cutout[arry.size()-1]->InitMidTopPoint(aps,1.0/3.0);
	cutout[arry.size() - 2]->InitMidTopPoint(aps,2.0/3.0);
	cutout[arry.size() - 3]->InitMidTopPoint(aps,1.0/4.0);
	if (aps.size() < 3) {
		cout << "aps not enough!" << endl;
		return 0;
	}
	vector<Vector4f> moutline;
	SurfacePure *uptopt = new SurfacePure(aps[0], aps[1], aps[2]);
	for (int i = 0; i < arry.size(); i++) {
		cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
		if (!cutout[i]->AllocateXCoe(uptopt)) {
			cout << "Metara "<<i<< " Allocate Failed!" << endl;
			return 0;
		}// (meta, p2[0], p2[2]);//��ɢ

		//cutout[i]->OutlineEigenM(&outline);//������ɢ�������
		//cutout[i]->OutlineEigen(&outline);
		//ManageObj::OutFilePointObj(&outline, increname(i, "Metara-E-").c_str());
		//outline.clear();

		/*cutout[i]->OutlineEigen(&moutline);
		ManageObj::OutFilePointObj(&moutline, increname(i, "Metara-Normal-").c_str());
		moutline.clear();*/
	}
	ios.ShoeExpansion(cutout);
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

	float ex=sfc->FindAddLenth(meta,exp);
	ios.ShoeAddLength(sfc->ReturnStartPoint(), meta, ex);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Len-ext-Large.stl", 1);
#endif // OUTPUTFILET
	cout << "�ӳ����ν���!" << endl;
	return 1;
}

int main(int argc, char* argv[])
{
	//string sts = "3.36";
	//string sts2 = "206.381052 0.234966 10.934925";
	//float ap[3];
	//
	//regex e("([0-9]+\\s)(\\.?)([0-9][0-9])");
	//smatch sm;
	//regex_match(sts.cbegin(),sts.cend(),sm,e);
	//for (int i = 0; i < sm.size(); i++) {
	//	cout << sm[i] << endl;
	//}

	//if (argc != 2) {
	//	printf("��ָ��Ь鸳�ʼ�������ļ�·����\n");
	//	return 0;
	//}
	//ManageObj indextxt(argv[1]);

	////ManageObj indextxt("C:\\Users\\47108\\Desktop\\InputFile.txt");
	//if (indextxt.OpenObj()) {
	//	printf("��ʼ���ļ���ȡʧ�ܣ���������ļ�·�������Ƿ���ȷ��\n");
	//	return 0;
	//}
	//string instlfilename;
	//string outstlfilename;
	///* ������������� */
	MyMesh::Point p2[3];
	///*��Χ*/
	MyMesh::Point p[3];
	//p[2] = MyMesh::Point(0, 0, 0);
	//float heelhight; //����
	//float giveout[] = { 0,0,0,0 }; //���ȣ���Χ����Χ����Χ
	//int whstraighten=1;
	//
	//indextxt.ReadInitTxtFile(instlfilename, p2,p,heelhight,giveout,whstraighten,outstlfilename);
	//printf("Ь鸶�ȡ����·��: %s\n",instlfilename.c_str());
	//printf("���� : %fmm\n", heelhight);
	//printf("�������ӣ� %fmm\n", giveout[0]);
	//printf("��Χ���ӣ� %fmm\n", giveout[1]);
	//printf("��Χ���ӣ� %fmm\n", giveout[2]);
	//printf("��Χ���ӣ� %fmm\n", giveout[3]);
	//printf("Ь��������·��: %s\n", outstlfilename.c_str());
	//if (whstraighten) {
	//	printf("ģ�Ͱ�������������ģ�������\n");
	//}
	//else {
	//	printf("û�б�Ҫ����ģ�͡�\n");
	//}
	//cout << endl;

	/*p2[0] = MyMesh::Point(237.5670, 1.2265, 11.8482);
	p2[1] = MyMesh::Point(49.4055, -0.9980, 184.6677);
	p2[2] = MyMesh::Point(28.649047, 0.477465, 90.453272);*/
	//PFH-10
	/*p2[0] = MyMesh::Point(251.237072, - 1.049236, 12.656890);
	p2[1] = MyMesh::Point(3.140964, 0.610319, 101.681542);
	p2[2] = MyMesh::Point(3.739703, - 2.403674, 10.501304);*/

	////PFH-30
	/*p2[0] = MyMesh::Point(234.318844 ,0.371156, 12.771103);
	p2[1] = MyMesh::Point(4.988612, 0.331804, 150.494285);
	p2[2] = MyMesh::Point(2.467105, - 3.981125, 59.876513);*/

	//PF-010-60-6-YX-YP.stl
	p2[0] = MyMesh::Point(227.634315, 1.743002, 10.427064);
	p2[1] = MyMesh::Point(21.450806, 1.800634, 130.709625);
	p2[2] = MyMesh::Point(1.305391, -0.475013, 59.906918);

	////PF - 027 - 90 - 6 - YX.stl
	//p2[0] = MyMesh::Point(206.381052, 0.234966, 10.934925);
	//p2[1] = MyMesh::Point(31.406775, 0.118433, 157.515320);
	//p2[2] = MyMesh::Point(0.720481, 0.177762, 90.235443);

	//float coeabc = OutLine::GiveCoe(p2[0], p2[1], p2[2], &coe);  //??
	
	/*p[0] = MyMesh::Point(149.726660, -45.038730, 7.791545);
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = MyMesh::Point(162.523123, 36.191162, 4.322560);*/

	//PFH-10-BZ-6
	/*p[0] = MyMesh::Point(154.442764, -45.383193, 5.547377);
	p[1] = MyMesh::Point(0,0,0);
	p[2] = MyMesh::Point(165.761880,35.007070, 5.553459);*/

	////PFH-60
	/*p[0] = MyMesh::Point(147.624629, 35.107334, 5.667529);
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = MyMesh::Point(138.070815, -44.152159, 8.595856);*/

	//PF-010-60-6-YX-YP.stl
	p[0] = MyMesh::Point(149.144440, 34.511948, 3.422032);
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = MyMesh::Point(135.451709, -44.498840, 5.291002);

	//PF-027-90-6-YX.stl
	/*p[0] = MyMesh::Point(128.3734, 32.8419, 2.9396);
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = MyMesh::Point(119.4374, -45.9890, 6.9330);*/

	vector<Vector3f> outline;
	MyOpenMesh ios;
	//ios.ReadStlfile(instlfilename.c_str());

	MyMesh::VertexHandle vertex[3], vertex2[3];//
	cout << "Now is Reading STL file!" << endl;
	ios.ReadStlfile("PF-010-60-6-YX-YP.stl");
	//ios.ReadStlfile("PF-027-90-6-YX.stl");
	
	//ios.ReadStlfile("shoestl.stl");
	//ios.ReadStlfile("PFH-10-BZ-6.stl");
	//ios.ReadStlfile("PFH-60-BZ-6.stl");
	//ios.ReadStlfile(stlfilename.c_str());
	float heelhight = 60;
	float giveout[] = { 0,5,5,3 }; //���ȣ���Χ����Χ����Χ, ֺΧ1λ�ã�ֺΧ1���Σ�ֺΧ2λ�ã�ֺΧ2���Ρ�����
	float toechange[] = { 179,2 };

	/*���µ�˼·����ʹ����������ֵ����Լ���ײ�У�����ϵķ������д���*/
	//vector<MyMesh::Point> bottomline;//�ײ������ߵ���ȡ
	//vector<MyOutBottom>  botall;
	//botall=ios.ShoeBottomLine(bottomline);//��Ҫ�����ͬЬ��ƥ��ʹ�����⣬һ���潫������λ�����������յıȽ�����Ľ����1.2~1.3֮�䣬һ����������
	//ios.ShoeBottomLine(botall);
	//ManageObj::OutFilePointObj(bottomline, "bottomoutlien.obj");
	//ManageObj::OutFilePointAna(botall,"botana.obj");
	//ios.FindFloorContour(bottomline);
	//ManageObj::OutFilePointObj(bottomline, "bottom-z.obj");
	//ios.InitTriPoint(vertex);//ʹ���������㣬x������Զ��z������ߣ�x������������õ��ĵ㲢�����룬���һ�Ҫ����ģ�Ͱ��������⣬��һ�����ظ��Թ��������Բ����ã����ĸ������룡��������Ŀ���ģ���Զ�����ʶ������⣨���ӣ�
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	ios.FindNearest(p[0], p[1], p[2], vertex2);		//��ʼ�����λ�� p[0], p[1], p[2],

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back,*toe=NULL;

	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}

	if (0) {
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
	
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ
	meta->ReturnTriPoint(vertex2,ios);//����hander�������³�ʼ����Χ
	//meta->SetMIth(sfc);

	/*meta->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline, "meta-outline.obj");
	outline.clear();*/

	//wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	//back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ

	if (giveout[0]) {
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
		//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//���³�ʼ����Χ
	}
	
	MyMesh::VertexHandle vertex_wist[3], vertex_back[3], vertex_toe[3];

	wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ

	float lentoe = 0;
	if (toechange[0]) {
		toe = sfc->FindToeBottomPoint(toechange[0]);
		lentoe = toe->ReturnLength();
		toe->ReturnTriPoint(vertex_toe, ios);
	}

	float lenwist = wist->ReturnLength();
	float lenback = back->ReturnLength();
	
	wist->ReturnTriPoint(vertex_wist,ios);
	back->ReturnTriPoint(vertex_back,ios);
	
	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"sfc-outline.obj");
	outline.clear();*/

	if (giveout[1]) {
		judge = MetaraExpansion(ios, meta, sfc,heelhight ,giveout[1]);  //��Χ�ӷ�
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
		//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//��Χ����Ҫ����Ѱ����Χ���Ѿ���ʼ��������Χ

		//wist=sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
		//back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ
		wist = new SurfaceCoe(vertex_wist, ios.mesh);
		wist->Init(0);
		wist->SetMIth(sfc);
		back = new SurfaceCoe(vertex_back, ios.mesh);
		back->Init(0);
		back->SetMIth(sfc);

		if (toechange[0]) {
			delete toe;
			toe = new SurfaceCoe(vertex_toe, ios.mesh);
			toe->Init(0);
			toe->SetMIth(sfc);
		}
	}

	if (toechange[0]) {
		lentoe = toe->ReturnLength() - lentoe;
		lentoe = toechange[1] - lentoe;
		ToeExpansion(ios, sfc, toe, toechange[0], lentoe);
	}

	if (giveout[2]) {
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
			//wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
			//back = sfc->FindWaistLine(wist);
			wist = new SurfaceCoe(vertex_wist, ios.mesh);
			wist->Init(0);
			wist->SetMIth(sfc);
			back = new SurfaceCoe(vertex_back, ios.mesh);
			back->Init(0);
			back->SetMIth(sfc);
		}
	}
	
	if (giveout[3]) {
		lenback = back->ReturnLength() - lenback;
		lenback = giveout[3] - lenback;
		if (abs(lenback)>0.005) {
			judge = BackExpansion(ios, back, wist, sfc, lenback);		 //��Χ�ӷ�
			if (!judge) {
				cout << "Main error 3" << endl;
				return -1;
			}
			/*delete sfc;
			delete wist;
			delete back;
			meta = new SurfaceCoe(vertex2, ios.mesh);
			meta->Init(0);

			wist = new SurfaceCoe(vertex_wist, ios.mesh);
			wist->Init(0);
			back = new SurfaceCoe(vertex_back, ios.mesh);
			back->Init(0);*/
		}
	}

	//ios.WriteStlfile(outstlfilename.c_str(), 1);

	

	delete sfc;
	delete meta;
	delete wist;
	delete back;

	system("pause");//For Debug;
	return 0;
}


