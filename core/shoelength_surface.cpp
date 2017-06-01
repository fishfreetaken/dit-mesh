// shoelength.cpp : 定义控制台应用程序的入口点。
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

//#define OUTPUTFILET 1 //控制调试的时候是否输出中间生成文件

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

	QuaternionSpin m(ab,Vector3f(0,1,0),p[2],outline, heel);  //默认p[2]为原点 有一点必须是后跟点
	if (!m.TransferSpin()) 
	{
		cout << "Not Find Bottomest Point!" << endl;
		return 0;
	}
	ios.ShoeSpin(m.ReturnQuatFuse(), m.ReturnShift());//旋转


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

	cout << "模型摆正并输出摆正模型！" << endl;
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
	cout << "正在对背围变形..." << endl;
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
	cout << "背围变形结束!" << endl;
	return 1;
}

bool WistExpansion(MyOpenMesh&ios,SurfaceCoe*&meta, SurfaceCoe*&back, SurfaceCoe*&wist, float exp) 
{
	cout << "正在对腰围变形..." << endl;
	wist->InitMidEndPoint();
	wist->AllocateXCoe(exp);
	//cout << wist->ReturnLength() << endl;
	ios.ShoeExpansionWist2(meta, wist, back);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET

	cout << "腰围变形结束!" << endl;
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, SurfaceCoe* &meta,SurfaceCoe* &sfc,float heelhight,float exp) {
	vector<Vector3f>outline; //for debug output
	vector<MySurCutArry> arry = sfc->OutCutOutline(exp, meta, heelhight);
	vector<SurfaceCoe*> cutout(arry.size());

	vector<MyMesh::Point> arp;
	vector<MyMesh::Point> aps;//提取碗口平面三点

	/*ManageObj arps("33.txt");
	if (!arps.ReadMeshPoints(arp)) {
		cout << "Metara Read File Eroor" << endl;
		return;
	}*/

	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitMidEndPoint(arp);//新算法

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
		}// (meta, p2[0], p2[2]);//扩散

		//cutout[i]->OutlineEigenM(&outline);//到处扩散后的轮廓
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
	cout << "掌围变形结束!" << endl;
	return 1;
}

bool MoveLength(MyOpenMesh&ios, SurfaceCoe*meta, SurfaceCoe* &sfc, float exp) {
	//cout << "move add lenth..." << endl;

	float ex=sfc->FindAddLenth(meta,exp);
	ios.ShoeAddLength(sfc->ReturnStartPoint(), meta, ex);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Len-ext-Large.stl", 1);
#endif // OUTPUTFILET
	cout << "加长变形结束!" << endl;
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
		std::cout << "\n***欢迎使用DIT鞋楦变形测试软件***\n" << endl;
		std::cout << "(1) 请输入STL模型路径 >>";
		cin >> inStlFilename;
		std::cout << "(2) 请输入鞋楦特征点文件路径>> ";
		cin >> featureFileName;
		std::cout << "(3) 请输入长度变形量>>";
		cin >> dLength;
		std::cout << "(4) 请输入掌围变形量>>";
		cin >> pawGirth;
		std::cout << "(5) 请输入腰围变形量>>";
		cin >> waistGirth;
		std::cout << "(6) 请输入背围变形量>>";
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
	//判断文件是否存在
	ManageObj indextxt((char*)featureFileName.c_str());
	if (indextxt.OpenObj()) 
	{
		printf("初始化文件读取失败，请检查你的文件路径配置是否正确。\n");
		return 0;
	}

	/* 纵向横切面三点 */
	MyMesh::Point p2[3];
	/*掌围*/
	MyMesh::Point p[3];
	p[2] = MyMesh::Point(0, 0, 0);
	float heelhight=0; //跟高
	float giveout[4] = {0}; // = { 0,0,0,0 }; //长度，掌围，腰围，背围
	giveout[0] = dLength; // 长度
	giveout[1] = pawGirth;
	giveout[2] = waistGirth;
	giveout[3] = backGirth;
	//int whstraighten=1;
	
	// 读取鞋楦特征点
	indextxt.ReadInitTxtFile2(p2, p, heelhight);

	printf("鞋楦读取配置路径: %s\n", inStlFilename.c_str());
	printf("跟高 : %fmm\n", heelhight);
	printf("长度增加： %fmm\n", giveout[0]);
	printf("掌围增加： %fmm\n", giveout[1]);
	printf("腰围增加： %fmm\n", giveout[2]);
	printf("背围增加： %fmm\n", giveout[3]);
	printf("鞋楦输出配置路径: %s\n", outStlFilename.c_str());
	//if (whstraighten) 
	//{
	//	printf("模型摆正并将摆正后模型输出。\n");
	//}
	//else 
	//{
	//	printf("没有必要摆正模型。\n");
	//}
	cout << endl;
	
	vector<Vector3f> outline;
	// 读取STL文件
	MyOpenMesh ios;
	ios.ReadStlfile(inStlFilename.c_str());

	MyMesh::VertexHandle vertex[3], vertex2[3];//

	/*最新的思路可以使用向量方差值最大，以及最底层校验相结合的方法进行处理*/

	//ios.InitTriPoint(vertex);//使用特征三点，x轴向最远，z轴向最高，x轴向最近，所得到的点并不理想，而且还要面临模型摆正的问题，有一定的重复性工作，所以不适用，重心更不理想！可以深入的考虑模型自动摆正识别的问题（复杂）
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	ios.FindNearest(p[0], p[2], p[1], vertex2);		//初始化最近位置 p[0], p[1], p[2],

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back;

	sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
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
		//cout << "模型摆正" << endl;
		Straighten(vertex, ios, sfc, heelhight, outStlFilename); //这个一定要给出跟高	！
		delete sfc; //摆正之后需要重新初始化sfc中轴横截面
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) 
		{
			cout << "sfc init error!" << endl;
			return 0;
		}
	}
	
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围

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
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) 
		{
			cout << "add len sfc init error!" << endl;
			return  0;
		}
		meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//重新初始化掌围
	}

	wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
	back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围
	float lenwist = wist->ReturnLength();
	float lenback = back->ReturnLength();
	
	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"sfc-outline.obj");
	outline.clear();*/

	if (giveout[1]) //掌围加肥
	{
		judge = MetaraExpansion(ios, meta, sfc,heelhight ,giveout[1]);  
		if (!judge) {
			cout << "Main error 1" << endl;
			return -1;
		}
		delete meta;
		delete sfc;
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) {
			cout << "Metara sfc init error!" << endl;
			return 0;
		}
		meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//扩围后需要重新寻找掌围，已经初始化过的掌围
		wist=sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
		back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围
	}
	
	if (giveout[2]) 
	{
		lenwist = wist->ReturnLength() - lenwist;
		lenwist = giveout[2] - lenwist;
		if (abs(lenwist)>0.005) {//定义一个精度范围
			judge = WistExpansion(ios, meta, back, wist, lenwist); //腰围加肥
			if (!judge) {
				cout << "Main error 2" << endl;
				return -1;
			}
			delete wist;
			delete meta;
			delete back;
			delete sfc;
			sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
			if (!sfc->Init(1)) {
				cout << "Wist sfc init error!" << endl;
				return 0;
			}
			meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围
			wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
			back = sfc->FindWaistLine(wist);
		}
	}
	
	if (giveout[3]) 
	{
		lenback = back->ReturnLength() - lenback;
		lenback = giveout[3] - lenback;
		if (abs(lenback)>0.005) {
			judge = BackExpansion(ios, back, wist, sfc, lenback);		 //背围加肥
			if (!judge) 
			{
				cout << "Main error 3" << endl;
				return -1;
			}
			/*delete sfc;
			delete wist;
			delete back;*/
			//sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
			//if (!sfc->Init(1)) {
			//	cout << "Wist sfc init error!" << endl;
			//	return 0;
			//}
			//wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
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


