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

bool ToeExpansion(MyOpenMesh&ios, SurfaceCoe*sfc, SurfaceCoe*toe,float dist,float exp,string name)
{
	vector<Vector3f> outline;
	//cout << "正在对趾围变形..." << endl;
	string cnoutline = name;
	string cnstl = name;
	cnoutline += ".obj";
	cnstl += ".stl";
	
	float range = 7;//前后控制平滑距离mm//原始默认的是9mm,不知道这个
	SurfaceCoe *toea, *toeb;
	toea = sfc->SfcMoveXLen(toe,dist+ range*2);//控制在前后7毫米左右
	toeb = sfc->SfcMoveXLen(toe, dist- range);//后
	toe->AllocateXCoe(exp);//将该围度扩大到指定围长

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
	//cout << "趾围变形结束!" << endl;
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

	wist->AllocateXCoe(exp);

	ios.ShoeExpansionWist2(meta, wist, back);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET

	cout << "腰围变形结束!" << endl;
	return 1;
}

bool MetaraExpansion(MyOpenMesh&ios, SurfaceCoe* &meta,SurfaceCoe* &sfc,float heelhight,float exp) {
	vector<Vector3f>outline; //for debug output
	cout << "正在对掌围进行变形..." << endl;
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
	vector<MyMesh::Point> aps;//提取碗口平面三点

	string ccc = "cutout.obj";
	
	//string pp = "abcdefghijklmnopqrstuvwxyz0123456";
	for (int i = 0; i < arry.size(); i++) {
		cutout[i] = new SurfaceCoe(arry[i].n, arry[i].a, arry[i].x, ios.mesh);
		if (cutout[i]->Init()) {
			cutout[i]->InitMidEndPoint(arp);//新算法
			//cout << cutout[i]->ReturnLength() << endl;

			/*cutout[i]->OutlineEigen(&outline);
			ccc = to_string(i) + ccc;
			ManageObj::OutFilePointObj(&outline,ccc.c_str());
			ccc = "cutout.obj";*/
		}
	}
	ccc = "cutout2.obj";

	/*cutout[arry.size()-1]->InitMidTopPoint(aps,1.0/3.0);//这个是初始化一个碗口的平面
	cutout[arry.size() - 2]->InitMidTopPoint(aps,2.0/3.0);
	cutout[arry.size() - 3]->InitMidTopPoint(aps,1.0/4.0);
	if (aps.size() < 3) {
		cout << "aps not enough!" << endl;
		return 0;
	}*/

	//SurfacePure *uptopt = new SurfacePure(aps[0], aps[1], aps[2]);
	float exteli = 0;//这个用来和前面的扩展系数保持一致

	float maxf=abs(meta->DistSurface(sfc->ReturnStartPoint())); //起始点到掌围面最远距离
	
	for (int i = 0; i < arry.size(); i++) {
		cutout[i]->InitTwoPoints(arp[i * 2], arp[i * 2 + 1]);
		//if (!cutout[i]->AllocateXCoe(uptopt)) {
		//	cout << "Metara "<<i<< " Allocate Failed!" << endl;
		//	return 0;
		//}// (meta, p2[0], p2[2]);//扩散
		//float mmn=abs(meta->DistSurface(cutout[i]->ReturnStartPoint()));

		if (i < ithmet) {
			/*if (!cutout[i]->AllocateXCoe(2*mmn/maxf-1,i)) {
				cout << "Metara " << i << " Allocate Failed!" << endl;
				return 0;
			}*/
			
			//if (!cutout[i]->AllocateXCoe()) { //先加厚，再加宽
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

	ios.ShoeExpansion(cutout, sfc,meta);//这个用来meta进行区分检验算法
	//ios.ShoeExpansion(cutout, sfc);  //for debug;

	//ManageObj::OutFilePointObj(&outline, "endbottom.obj");

	ios.MetaraFileter(cutout[ithmet]);//掌围线平滑，有一定必要

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
	cout << "鞋楦加长..." << endl;
	float ex=sfc->FindAddLenth(meta,exp);
	ios.ShoeAddLength(sfc->ReturnStartPoint(), meta, ex);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Len-ext-Large.stl", 1);
#endif // OUTPUTFILET
	cout << "加长变形结束!" << endl;
	return 1;
}

int LastDeformation(string in, string out, MyMesh::Point *psfc,MyMesh::Point *pmeta,float heelhight,float * giveout,vector<vector<float>>&vv_toe) {
	///* 纵向横切面三点 */
	MyMesh::Point p2[3];
	/////*掌围*/
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
	//float giveout[] = { 0,5,5,3 }; //长度，掌围，腰围，背围, 趾围1位置，趾围1变形，趾围2位置，趾围2变形。。。
	//float toechange[] = { 179,2 };

	ios.FindNearest(p2[0], p2[1], p2[2], vertex);
	ios.FindNearest(p[0], p[1], p[2], vertex2);		//初始化最近位置 p[0], p[1], p[2],

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back,*toe;

	sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}
	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"outline.obj");*/

	if (0) {	//一定要摆正
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

	//先找趾围，掌围，腰围，背围，然后再变形；

	MyMesh::VertexHandle vertex_wist[3], vertex_back[3];

	vector<float> vlentoe;  //脚长变形前和变形后的问题；
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

	//if (vv_toe.size() == 1) { //只有两种情况吧，一种是只增加一条趾围线，另一种是增加两条趾围线
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
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围
	meta->ReturnTriPoint(vertex2, ios);//返回hander用于重新初始化掌围
	vlenmeta = meta->ReturnLength();

	/*meta->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"metaralineA.obj");*/

	wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
	back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围

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
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) {
			cout << "add len sfc init error!" << endl;
			return  0;
		}
		meta = new SurfaceCoe(vertex2, ios.mesh);
		meta->Init(0);
		meta->SetMIth(sfc);
		meta->ReturnLength();
		//vlenmeta = meta->ReturnLength()-vlenmeta;
		//meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//重新初始化掌围
	}

	/*sfc->OutlineEigen(&outline);
	ManageObj::OutFilePointObj(&outline,"outline-len.obj");*/

	if (*(giveout+1)) {
		//giveout[1] -= vlenmeta;//加上这个好像变的过长了
		judge = MetaraExpansion(ios, meta, sfc, heelhight, giveout[1]);  //掌围加肥
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
		cout << "正在对趾围" << j << "变形..." << endl;
		vlentoe[j] += (vv_toe[j][1]- vsfcoe[j]->ReturnLength());
		if (vlentoe[j] < 0.005) {
			continue;
		}
		ToeExpansion(ios, sfc, vsfcoe[j], vv_toe[j][0], vlentoe[j],cctoename[j]);
		cout << "趾围" << j << "变形结束" << endl;
	}

	for (int j = 0; j < vv_toe.size(); j++) {
		delete vsfcoe[j];
	}

	if (*(giveout+2)) {
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
			judge = BackExpansion(ios, back, wist, sfc, lenback);		 //背围加肥
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
		//cout << "请输入TXT变形配置文本文件路径：" << endl;
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
		cout << "变形结束，请检查结果！" << endl;
	//}
		system("pause");
	return 0;
}


