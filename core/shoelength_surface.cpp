// shoelength.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
//#include <afxwin.h>
//#include <GL/glut.h>
//#include <flann/flann.h>

#include "PointVector.h"
#include "ManageObj.h"
#include "QuaternionSpin.h";
#include "MyOpenMesh.h"
#include "ShoeLastDeformation.h"
// ----------------------------------------------------------------------------

//#define OUTPUTFILET 1 //���Ƶ��Ե�ʱ���Ƿ�����м������ļ�

int Straighten(MyMesh::VertexHandle *vertex, MyOpenMesh&ios, SurfaceCoe* &sfc,float heel, string sf) 
{
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
	if (!m.TransferSpin()) 
	{
		cout << "Not Find Bottomest Point!" << endl;
		return 0;
	}
	ios.ShoeSpin(m.ReturnQuatFuse(), m.ReturnShift());//��ת


	int pos =sf.find_last_of('.');
	string outf;
	if (pos < 0) 
	{
		outf += "Shoe-Spin-Model.stl";
	}
	else 
	{
		outf += sf.substr(0, pos);
		outf += "-spin.stl";
	}
	//ios.WriteStlfile("Shoe-Spin-Model.stl", 1);
	ios.WriteStlfile(outf.c_str(), 1);

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

bool BackExpansion(MyOpenMesh&ios, SurfaceCoe*back, SurfaceCoe*wist, SurfaceCoe* &sfc, float exp) {
	cout << "���ڶԱ�Χ����..." << endl;
	MyMesh::Point pcc;
	int iith = sfc->UpOneInch(back->ReturnIth(2), pcc);
	MyMesh::VertexHandle start =sfc->FindNearest(pcc);
	SurfaceCoe *metc = new SurfaceCoe(back->ReturnCoe(), start, 0, ios.mesh);
	metc->Init();
	metc->SetMIth(iith);

	back->InitMidEndPoint();
	back->AllocateXCoe(exp);
	ios.ShoeExpansionWist2(wist,back,metc);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET
	//ios.WriteStlfile("Last-ext-Large.stl", 1);
	cout << "��Χ���ν���!" << endl;
	return 1;
}

bool WistExpansion(MyOpenMesh&ios,SurfaceCoe*&meta, SurfaceCoe*&back, SurfaceCoe*&wist, float exp) 
{
	cout << "���ڶ���Χ����..." << endl;
	wist->InitMidEndPoint();
	wist->AllocateXCoe(exp);
	//cout << wist->ReturnLength() << endl;
	ios.ShoeExpansionWist2(meta, wist, back);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET

	cout << "��Χ���ν���!" << endl;
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, SurfaceCoe* &meta,SurfaceCoe* &sfc,float heelhight,float exp) {
	vector<Vector3f>outline; //for debug output
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

			/*cutout[i]->OutlineEigen(&outline);
			ManageObj::OutFilePointObj(&outline, increname(i,"Metara-").c_str());
			outline.clear();*/
		}
	}
	//ManageObj::OutFileOutlinePointXyz(&arp, "pair.xyz");//output arp for debug;

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

void main()
{
	while(1)
	{
		string featureFileName;// = "e:\\Exchange\\ChenDong\\Software\\shoestla\\shoestla.txt";
		string inStlFilename;// = "e:\\Exchange\\ChenDong\\Software\\shoestla\\shoestla.stl";
		string outStlFilename;// = "e:\\Exchange\\ChenDong\\Software\\shoestla\\shoestla2.stl";


		float dLength, pawGirth,  waistGirth, backGirth;
		std::cout << "\n***��ӭʹ��DITЬ鸱��β������***\n" << endl;
		std::cout << "(1) ������STLģ��·�� >>";
		cin >> inStlFilename;
		std::cout << "(2) ������Ь��������ļ�·��>> ";
		cin >> featureFileName;
		std::cout << "(3) �����볤�ȱ�����>>";
		cin >> dLength;
		std::cout << "(4) ��������Χ������>>";
		cin >> pawGirth;
		std::cout << "(5) ��������Χ������>>";
		cin >> waistGirth;
		std::cout << "(6) �����뱳Χ������>>";
		cin >> backGirth;

		int pos = inStlFilename.find_last_of('.');
		string outf;
		if (pos > 0)
		{
			outStlFilename += inStlFilename.substr(0, pos);
			outStlFilename += "-deformed.stl";
		}
		cout << "inStlFilename" <<inStlFilename<<endl;
		cout <<"outStlFilename"<< outStlFilename << endl;
		ShoeLastDeformation(inStlFilename, outStlFilename, featureFileName, dLength, pawGirth, waistGirth, backGirth);

	}
}

int ShoeLastDeformation(std::string inStlFilename, std::string outStlFilename,
						std::string featureFileName,
						float dLength, float pawGirth, float  waistGirth, float backGirth)
{
	//�ж��ļ��Ƿ����
	ManageObj indextxt((char*)featureFileName.c_str());
	if (indextxt.OpenObj()) 
	{
		printf("��ʼ���ļ���ȡʧ�ܣ���������ļ�·�������Ƿ���ȷ��\n");
		return 0;
	}

	/* ������������� */
	MyMesh::Point p2[3];
	/*��Χ*/
	MyMesh::Point p[3];
	p[2] = MyMesh::Point(0, 0, 0);
	float heelhight=0; //����
	float giveout[4] = {0}; // = { 0,0,0,0 }; //���ȣ���Χ����Χ����Χ
	giveout[0] = dLength; // ����
	giveout[1] = pawGirth;
	giveout[2] = waistGirth;
	giveout[3] = backGirth;
	//int whstraighten=1;
	
	// ��ȡЬ�������
	indextxt.ReadInitTxtFile2(p2, p, heelhight);

	printf("Ь鸶�ȡ����·��: %s\n", inStlFilename.c_str());
	printf("���� : %fmm\n", heelhight);
	printf("�������ӣ� %fmm\n", giveout[0]);
	printf("��Χ���ӣ� %fmm\n", giveout[1]);
	printf("��Χ���ӣ� %fmm\n", giveout[2]);
	printf("��Χ���ӣ� %fmm\n", giveout[3]);
	printf("Ь��������·��: %s\n", outStlFilename.c_str());
	//if (whstraighten) 
	//{
	//	printf("ģ�Ͱ�������������ģ�������\n");
	//}
	//else 
	//{
	//	printf("û�б�Ҫ����ģ�͡�\n");
	//}
	cout << endl;
	
	vector<Vector3f> outline;
	// ��ȡSTL�ļ�
	MyOpenMesh ios;
	ios.ReadStlfile(inStlFilename.c_str());

	MyMesh::VertexHandle vertex[3], vertex2[3];//

	/*���µ�˼·����ʹ����������ֵ����Լ���ײ�У�����ϵķ������д���*/

	//ios.InitTriPoint(vertex);//ʹ���������㣬x������Զ��z������ߣ�x������������õ��ĵ㲢�����룬���һ�Ҫ����ģ�Ͱ��������⣬��һ�����ظ��Թ��������Բ����ã����ĸ������룡��������Ŀ���ģ���Զ�����ʶ������⣨���ӣ�
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	ios.FindNearest(p[0], p[2], p[1], vertex2);		//��ʼ�����λ�� p[0], p[1], p[2],

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back;

	sfc = new SurfaceCoe(vertex, ios.mesh); //��������
	if (!sfc->Init(1)) 
	{
		cout << "sfc init error!" << endl;
		return 0;
	}

	if (1) 
	{
		if (!heelhight) 
		{
			cout << "heel hight not ini!" << endl;
			return 0;
		}
		//cout << "ģ�Ͱ���" << endl;
		Straighten(vertex, ios, sfc, heelhight, outStlFilename); //���һ��Ҫ��������	��
		delete sfc; //����֮����Ҫ���³�ʼ��sfc��������
		sfc = new SurfaceCoe(vertex, ios.mesh); //��������
		if (!sfc->Init(1)) 
		{
			cout << "sfc init error!" << endl;
			return 0;
		}
	}
	
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ

	if (giveout[0])
	{
		judge = MoveLength(ios, meta, sfc, giveout[0]);
		if (!judge) 
		{
			cout << "Main error 0" << endl;
			return -1;
		}
		delete sfc;
		delete meta;
		sfc = new SurfaceCoe(vertex, ios.mesh); //��������
		if (!sfc->Init(1)) 
		{
			cout << "add len sfc init error!" << endl;
			return  0;
		}
		meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//���³�ʼ����Χ
	}

	wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
	back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ
	float lenwist = wist->ReturnLength();
	float lenback = back->ReturnLength();
	
	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"sfc-outline.obj");
	outline.clear();*/

	if (giveout[1]) //��Χ�ӷ�
	{
		judge = MetaraExpansion(ios, meta, sfc,heelhight ,giveout[1]);  
		if (!judge) {
			cout << "Main error 1" << endl;
			return -1;
		}
		delete meta;
		delete sfc;
		sfc = new SurfaceCoe(vertex, ios.mesh); //��������
		if (!sfc->Init(1)) {
			cout << "Metara sfc init error!" << endl;
			return 0;
		}
		meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//��Χ����Ҫ����Ѱ����Χ���Ѿ���ʼ��������Χ
		wist=sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
		back = sfc->FindWaistLine(wist); //��ֱ��x��Ѱ�ұ�Χ
	}
	
	if (giveout[2]) 
	{
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
			meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//Ѱ����Χ���Ѿ���ʼ��������Χ
			wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
			back = sfc->FindWaistLine(wist);
		}
	}
	
	if (giveout[3]) 
	{
		lenback = back->ReturnLength() - lenback;
		lenback = giveout[3] - lenback;
		if (abs(lenback)>0.005) {
			judge = BackExpansion(ios, back, wist, sfc, lenback);		 //��Χ�ӷ�
			if (!judge) 
			{
				cout << "Main error 3" << endl;
				return -1;
			}
			/*delete sfc;
			delete wist;
			delete back;*/
			//sfc = new SurfaceCoe(vertex, ios.mesh); //��������
			//if (!sfc->Init(1)) {
			//	cout << "Wist sfc init error!" << endl;
			//	return 0;
			//}
			//wist = sfc->FindWaistLine(meta); //��ֱ��x��Ѱ����Χ
			//back = sfc->FindWaistLine(wist);
			//cout << back->ReturnLength() << endl;
		}
	}

	ios.WriteStlfile(outStlFilename.c_str(), 1);

	delete sfc;
	delete meta;
	delete wist;
	delete back;

	//system("pause");//For Debug;
	return 0;
}


