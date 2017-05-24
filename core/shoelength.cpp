// shoelength.cpp : 定义控制台应用程序的入口点。
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

#ifdef OUTPUTFILET
	ios.WriteStlfile("Shoe-Spin-MS.stl", 1);
#endif // OUTPUTFILET
	cout << "Shoe Spin Over" << endl;
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
	MyMesh::Point pcc;
	int iith = sfc->UpOneInch(back->ReturnIth(2), pcc);
	MyMesh::VertexHandle start =sfc->FindNearest(pcc);
	SurfaceCoe *metc = new SurfaceCoe(back->ReturnCoe(), start, 0, ios.mesh);
	metc->Init();
	metc->SetMIth(iith);

	back->InitMidEndPoint();
	back->AllocateXCoe(exp);
	ios.ShoeExpansionWist2(wist,back,metc);
	ios.WriteStlfile("Last-ext-Large.stl", 1);
	cout << "Back stl file output over!" << endl;

	return 1;
}

bool WistExpansion(MyOpenMesh&ios,SurfaceCoe*&meta, SurfaceCoe*&back, SurfaceCoe*&wist, float exp) {
	wist->InitMidEndPoint();
	wist->AllocateXCoe(exp);
	//cout << wist->ReturnLength() << endl;
	ios.ShoeExpansionWist2(meta, wist, back);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Wist-ext-Large.stl", 1);
#endif // OUTPUTFILET

	cout << "Wist expansion over!" << endl;
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
	cout << "Metara Expansion Over!" << endl;
	return 1;
}

bool MoveLength(MyOpenMesh&ios, SurfaceCoe*meta, SurfaceCoe* &sfc, float exp) {
	cout << "move add lenth..." << endl;

	float ex=sfc->FindAddLenth(meta,exp);
	ios.ShoeAddLength(sfc->ReturnStartPoint(), meta, ex);

#ifdef OUTPUTFILET
	ios.WriteStlfile("Len-ext-Large.stl", 1);
#endif // OUTPUTFILET
	cout << "Add Length OVER!" << endl;
	return 1;
}

int main(int argc, char* argv[])
{

	/* 纵向横切面三点 */
	MyMesh::Point p2[3];
	/*p2[0] = MyMesh::Point(237.5670, 1.2265, 11.8482);
	p2[1] = MyMesh::Point(49.4055, -0.9980, 184.6677);
	p2[2] = MyMesh::Point(28.649047, 0.477465, 90.453272);*/
	//PFH-10
	/*p2[0] = MyMesh::Point(251.237072, - 1.049236, 12.656890);
	p2[1] = MyMesh::Point(3.140964, 0.610319, 101.681542);
	p2[2] = MyMesh::Point(3.739703, - 2.403674, 10.501304);*/

	////PFH-30
	p2[0] = MyMesh::Point(234.318844 ,0.371156, 12.771103);
	p2[1] = MyMesh::Point(4.988612, 0.331804, 150.494285);
	p2[2] = MyMesh::Point(2.467105, - 3.981125, 59.876513);

	//float coeabc = OutLine::GiveCoe(p2[0], p2[1], p2[2], &coe);  //??
	/*掌围*/
	MyMesh::Point p[3];
	/*p[0] = MyMesh::Point(149.726660, -45.038730, 7.791545);
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = MyMesh::Point(162.523123, 36.191162, 4.322560);*/

	//PFH-10-BZ-6
	/*p[0] = MyMesh::Point(154.442764, -45.383193, 5.547377);
	p[1] = MyMesh::Point(0,0,0);
	p[2] = MyMesh::Point(165.761880,35.007070, 5.553459);*/

	////PFH-60
	p[0] = MyMesh::Point(147.624629, 35.107334, 5.667529);
	p[1] = MyMesh::Point(0, 0, 0);
	p[2] = MyMesh::Point(138.070815, -44.152159, 8.595856);

	vector<Vector3f> outline;
	MyOpenMesh ios;

	MyMesh::VertexHandle vertex[3], vertex2[3];//
	ios.ReadStlfile("shoestl.stl");
	//ios.ReadStlfile("PFH-10-BZ-6.stl");
	//ios.ReadStlfile("PFH-60-BZ-6.stl");
	float heelhight = 90.06;
	float giveout[] = { 5,5,3,3 }; //长度，掌围，腰围，背围

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
	ios.FindNearest(p[0], p[1], p[2], vertex2);		//初始化最近位置

	bool judge = false;
	SurfaceCoe *sfc, *meta, *wist, *back;

	sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}

	if (!heelhight) {
		cout << "heel hight not ini!" << endl;
		return 0;
	}
	Straighten(vertex, ios, sfc, heelhight); //这个一定要给出跟高	！
	delete sfc; //摆正之后需要重新初始化sfc中轴横截面
	sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
	if (!sfc->Init(1)) {
		cout << "sfc init error!" << endl;
		return 0;
	}
	
	meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围

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
		meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//重新初始化掌围
	}

	wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
	back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围
	float lenwist = wist->ReturnLength();
	float lenback = back->ReturnLength();
	
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
		sfc = new SurfaceCoe(vertex, ios.mesh); //中轴横截面
		if (!sfc->Init(1)) {
			cout << "Metara sfc init error!" << endl;
			return 0;
		}
		meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//扩围后需要重新寻找掌围，已经初始化过的掌围
		wist=sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
		back = sfc->FindWaistLine(wist); //垂直于x轴寻找背围
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
			meta = sfc->FindMetara(vertex2[0], vertex2[2]);	//寻找掌围，已经初始化过的掌围
			wist = sfc->FindWaistLine(meta); //垂直于x轴寻找腰围
			back = sfc->FindWaistLine(wist);
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
	delete sfc;
	delete meta;
	delete wist;
	delete back;
	//back->OutlineEigen(&outline);// 将原有的outlien进行输出
	//ManageObj::OutFilePointObj(&outline, "outlineback-non.obj");
	//outline.clear();

	//wist->OutlineEigen(&outline);// 将原有的outlien进行输出
	//ManageObj::OutFilePointObj(&outline, "outlinewist-non.obj");
	//outline.clear();

	/*meta->OutlineEigen(&outline);				
	ManageObj::OutFilePointObj(&outline, "outlinezhi-non.obj");
	outline.clear();*/
	//meta->AllocateXCoe(stepforward);			//按照不同方式配置递增系数

	system("pause");
	return 0;
}


