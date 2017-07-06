// shoelength.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
//#include <GL/glut.h>
//#include <flann/flann.h>

#include "PointVector.h"
#include "ManageObj.h"
#include "QuaternionSpin.h";
#include "MyOpenMesh.h"
// ----------------------------------------------------------------------------

#define OUTPUTFILET 1 //控制调试的时候是否输出中间生成文件

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

	QuaternionSpin m(ab,Vector3f(0,1,0),p[2],outline, heel);  //默认p[2]为原点 有一点必须是后跟点
	if (!m.TransferSpin()) {
		cout << "Not Find Bottomest Point!" << endl;
		return 0;
	}
	ios.ShoeSpin(m.ReturnQuatFuse(), m.ReturnShift());//旋转


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

bool ToeExpansion(MyOpenMesh&ios, SurfaceCoe*sfc, SurfaceCoe*toe,float dist,float exp)
{
	//vector<Vector3f> outline;
	cout << "正在对趾围变形..." << endl;
	
	float range = 9;
	SurfaceCoe *toea, *toeb;
	toea = sfc->SfcMoveXLen(toe,dist+ range);//控制在前后7毫米左右
	toeb = sfc->SfcMoveXLen(toe, dist- range);
	toe->AllocateXCoe(exp);//将该围度扩大到指定围长
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
	cout << "趾围变形结束!" << endl;
	return 1;
}

bool BackExpansion(MyOpenMesh&ios, SurfaceCoe*back, SurfaceCoe*wist, SurfaceCoe* &sfc, float exp) {
	cout << "正在对背围变形..." << endl;
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
	cout << "背围变形结束!" << endl;
	delete metc;
	return 1;
}

bool WistExpansion(MyOpenMesh&ios,SurfaceCoe*&meta, SurfaceCoe*&back, SurfaceCoe*&wist, float exp) {
	cout << "正在对腰围变形..." << endl;

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

	cout << "腰围变形结束!" << endl;
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, SurfaceCoe* &meta,SurfaceCoe* &sfc,float heelhight,float exp) {
	//vector<Vector3f>outline; //for debug output
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
	//	printf("请指定鞋楦初始化配置文件路径。\n");
	//	return 0;
	//}
	//ManageObj indextxt(argv[1]);

	////ManageObj indextxt("C:\\Users\\47108\\Desktop\\InputFile.txt");
	//if (indextxt.OpenObj()) {
	//	printf("初始化文件读取失败，请检查你的文件路径配置是否正确。\n");
	//	return 0;
	//}
	//string instlfilename;
	//string outstlfilename;
	///* 纵向横切面三点 */
	MyMesh::Point p2[3];
	///*掌围*/
	MyMesh::Point p[3];
	//p[2] = MyMesh::Point(0, 0, 0);
	//float heelhight; //跟高
	//float giveout[] = { 0,0,0,0 }; //长度，掌围，腰围，背围
	//int whstraighten=1;
	//
	//indextxt.ReadInitTxtFile(instlfilename, p2,p,heelhight,giveout,whstraighten,outstlfilename);
	//printf("鞋楦读取配置路径: %s\n",instlfilename.c_str());
	//printf("跟高 : %fmm\n", heelhight);
	//printf("长度增加： %fmm\n", giveout[0]);
	//printf("掌围增加： %fmm\n", giveout[1]);
	//printf("腰围增加： %fmm\n", giveout[2]);
	//printf("背围增加： %fmm\n", giveout[3]);
	//printf("鞋楦输出配置路径: %s\n", outstlfilename.c_str());
	//if (whstraighten) {
	//	printf("模型摆正并将摆正后模型输出。\n");
	//}
	//else {
	//	printf("没有必要摆正模型。\n");
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
	float giveout[] = { 0,5,5,3 }; //长度，掌围，腰围，背围, 趾围1位置，趾围1变形，趾围2位置，趾围2变形。。。
	float toechange[] = { 179,2 };

	/*最新的思路可以使用向量方差值最大，以及最底层校验相结合的方法进行处理*/
	//vector<MyMesh::Point> bottomline;//底部轮廓线的提取
	//vector<MyOutBottom>  botall;
	//botall=ios.ShoeBottomLine(bottomline);//需要解决不同鞋子匹配使用问题，一方面将向量单位化，对于最终的比较理想的结果在1.2~1.3之间，一定存在杂项
	//ios.ShoeBottomLine(botall);
	//ManageObj::OutFilePointObj(bottomline, "bottomoutlien.obj");
	//ManageObj::OutFilePointAna(botall,"botana.obj");
	//ios.FindFloorContour(bottomline);
	//ManageObj::OutFilePointObj(bottomline, "bottom-z.obj");
	//ios.InitTriPoint(vertex);//使用特征三点，x轴向最远，z轴向最高，x轴向最近，所得到的点并不理想，而且还要面临模型摆正的问题，有一定的重复性工作，所以不适用，重心更不理想！可以深入的考虑模型自动摆正识别的问题（复杂）
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	ios.FindNearest(p[0], p[1], p[2], vertex2);		//初始化最近位置 p[0], p[1], p[2],

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back,*toe=NULL;

	sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}

	if (0) {
		if (!heelhight) {
			cout << "heel hight not ini!" << endl;
			return 0;
		}
		//cout << "模型摆正" << endl;
		Straighten(vertex, ios, sfc, heelhight); //这个一定要给出跟高	！
		delete sfc; //摆正之后需要重新初始化sfc中轴横截面
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) {
			cout << "sfc init error!" << endl;
			return 0;
		}
	}
	
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围
	meta->ReturnTriPoint(vertex2,ios);//返回hander用于重新初始化掌围
	//meta->SetMIth(sfc);

	/*meta->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline, "meta-outline.obj");
	outline.clear();*/

	//wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
	//back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围

	if (giveout[0]) {
		judge = MoveLength(ios, meta, sfc, giveout[0]);
		if (!judge) {
			cout << "Main error 0" << endl;
			return -1;
		}
		delete sfc;
		delete meta;
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) {
			cout << "add len sfc init error!" << endl;
			return  0;
		}
		meta = new SurfaceCoe(vertex2, ios.mesh);
		meta->Init(0);
		meta->SetMIth(sfc);
		//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//重新初始化掌围
	}
	
	MyMesh::VertexHandle vertex_wist[3], vertex_back[3], vertex_toe[3];

	wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
	back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围

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
		judge = MetaraExpansion(ios, meta, sfc,heelhight ,giveout[1]);  //掌围加肥
		if (!judge) {
			cout << "Main error 1" << endl;
			return -1;
		}
		delete meta;
		delete sfc;
		delete wist;
		delete back;
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) {
			cout << "Metara sfc init error!" << endl;
			return 0;
		}
		meta = new SurfaceCoe(vertex2, ios.mesh);
		meta->Init(0);
		meta->SetMIth(sfc);
		//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//扩围后需要重新寻找掌围，已经初始化过的掌围

		//wist=sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
		//back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围
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
			//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围
			meta = new SurfaceCoe(vertex2, ios.mesh);
			meta->Init(0);
			meta->SetMIth(sfc);
			//wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
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
			judge = BackExpansion(ios, back, wist, sfc, lenback);		 //背围加肥
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


