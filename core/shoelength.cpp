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


#define stepforward 5 //鞋楦长增加5mm 
#define mlength 147
#define blength 238
#define tlength  91
#define outputname "outline.xyz"
#define openfilename "shoestl.obj"


/*已知三个点求平面方程系数 A,B,C,D*/
Vector4f computecoe(Vector3f *p) {
	Vector4f x;
	x[0] = ((*(p + 1))[1] - (*(p + 0))[1])*((*(p + 2))[2] - (*(p + 0))[2]) - ((*(p + 1))[2] - (*(p + 0))[2])*((*(p + 2))[1] - (*(p + 0))[1]);
	x[1] = ((*(p + 1))[2] - (*(p + 0))[2])*((*(p + 2))[0] - (*(p + 0))[0]) - ((*(p + 1))[0] - (*(p + 0))[0])*((*(p + 2))[2] - (*(p + 0))[2]);
	x[2] = ((*(p + 1))[0] - (*(p + 0))[0])*((*(p + 2))[1] - (*(p + 0))[1]) - ((*(p + 1))[1] - (*(p + 0))[1])*((*(p + 2))[0] - (*(p + 0))[0]);
	x[3] = 0 - ((x[0] * ((*p)[0])) + (x[1] * ((*p)[1])) + (x[2] * ((*p)[2])));
	return x;
}
Vector4f ccomputecoe(Vector3f *p) {
	Vector3f ab = *(p+1) - *p;
	Vector3f ac = *(p+2) - *p;
	ab = ac.cross(ab);
	float d;
	Vector3f zero(0, 0, 0);
	d = ab.transpose()*(zero - *p);

	return Vector4f(ab[0], ab[1], ab[2], d);
}

unsigned int unTriangles;
int cpyint(const char*& p)
{
	int cpy;
	char *memwriter;
	memwriter = (char*)&cpy;
	memcpy(memwriter, p, 4);
	p += 4;
	return cpy;
}
float cpyfloat(const char*& p)
{
	float cpy;
	char *memwriter;
	memwriter = (char*)&cpy;
	memcpy(memwriter, p, 4);
	p += 4;
	return cpy;
}
/*已知平面三个点求平面法向量*/
Vector3f VectorCompute(Vector3f a, Vector3f b, Vector3f c) {
	Vector3f ab, ac;
	ab = b - a;
	ac = c - a;
	
	ab = ac.cross(ab);
	float norm;
	norm=sqrt(ab.transpose()*ab);
	ab = ab / norm;
	return ab;
}

PointVector *root;
vector<TriMesh> trilist;
bool ReadBinary(const char *buffer)
{
	const char* p = buffer;
	vector<Vector3f>pointList;
	char name[80];

	memcpy(name, p, 80);
	cout << name << endl;
	p += 80;
	unTriangles = cpyint(p); //357000
	cout << "total number triangles: " << unTriangles << endl;

	Vector3f vec(cpyfloat(p), cpyfloat(p), cpyfloat(p)); //三角面片的法向量
	Vector3f a(cpyfloat(p), cpyfloat(p), cpyfloat(p));   //三个顶点
	Vector3f b(cpyfloat(p), cpyfloat(p), cpyfloat(p));
	Vector3f c(cpyfloat(p), cpyfloat(p), cpyfloat(p));

	root = new PointVector(a,b,c,vec); //init

	trilist.push_back(TriMesh( root, root->Insert(b, c, a, vec), root->Insert(c, b, a, vec),vec,p));
	//root->Insert(b, c, a, vec);
	//root->Insert(c, b, a, vec);
	p += 2;//跳过尾部标志

	cout << "reading file..." << endl;
	for (unsigned int i = 0; i < unTriangles; i++)
	{
		if (!(i % 1000)) {
			cout <<"current step : "<< i << endl;
		}

		//p += 12;//跳过头部法向量
		Vector3f vec2(cpyfloat(p), cpyfloat(p), cpyfloat(p)); //三角面片的法向量
		Vector3f a2(cpyfloat(p), cpyfloat(p), cpyfloat(p));   //三个顶点
		Vector3f b2(cpyfloat(p), cpyfloat(p), cpyfloat(p));
		Vector3f c2(cpyfloat(p), cpyfloat(p), cpyfloat(p));
		
		/*cout << vec2 << endl;
		cout << "---------------" << endl;
		cout << VectorCompute(a2,b2,c2) << endl; //debug
		cout << "-----------------------" << endl;*/

		trilist.push_back(TriMesh(root->Insert(a2,b2,c2,vec2), root->Insert(b2, c2, a2, vec2), root->Insert(c2, b2, a2, vec2), vec2, p));
		/*root->Insert(a2, b2, c2, vec2);
		root->Insert(b2, a2, c2, vec2);
		root->Insert(c2, b2, a2, vec2);*/
		
		p += 2;//跳过尾部标志
	}
	cout << "trilist NO. : "<< trilist.size() << endl;
	return true;
}

/*
	x法向直接乘系数
*/
int JustFast(Vector3f *a, Vector4f cc) {
	Vector4f sa((*a)[2], (*a)[1], (*a)[0], 1);
	Vector3f sb(cc[0], cc[1], cc[2]);
	float lms = sa.transpose()*cc;
	if (sa.transpose()*cc < 0) {
		(a)[0] += ((a)[0] / mlength)*stepforward;
		return 1;
	}
	return 0;
}
struct mpair {
	Vector3f a;
	float s;
	unsigned int l;//点索引
};
vector<struct mpair> select;
int JustFast2(Vector3f a, Vector4f cc,vector<struct mpair> *sel,unsigned int th,float * max) {
	Vector4f sa((a)[0], (a)[1], (a)[2], 1);
	Vector3f sb(cc[0], cc[1], cc[2]);
	float lms = sa.transpose()*cc;
	if (lms < 0) {
		struct mpair mj;
		mj.l = th;
		mj.a = a;
		mj.s = (0 - lms) / sqrt(sb.transpose()*sb);//距离平面的距离
		if (mj.s > *max) {
			*max = mj.s;
		}
		
		sel->push_back(mj);
		return 1;
	}
	return 0;
}
vector<TriCesh> trilis;

bool ReadBinary2(const char *buffer,Vector4f cc) {
	const char* p = buffer;
	vector<Vector3f>pointList;
	char name[80];

	memcpy(name, p, 80);
	cout << name << endl;
	p += 80;
	unTriangles = cpyint(p); //357000
	cout << "total number triangles: " << unTriangles << endl;


	cout << "reading file..." << endl;
	for (unsigned int i = 0; i < unTriangles; i++)
	{
		if (!(i % 10000)) {
			cout << "current step : " << i << endl;
		}

		Vector3f vec(cpyfloat(p), cpyfloat(p), cpyfloat(p)); //三角面片的法向量
		Vector3f a(cpyfloat(p), cpyfloat(p), cpyfloat(p));   //三个顶点
		Vector3f b(cpyfloat(p), cpyfloat(p), cpyfloat(p));
		Vector3f c(cpyfloat(p), cpyfloat(p), cpyfloat(p));
		
		if ((JustFast(&a, cc))||(JustFast(&b, cc))||(JustFast(&c, cc))) {
			vec = VectorCompute(a, b, c);
		}
		//trilis.push_back(TriCesh(a, b, c, vec, p));

		p += 2;//跳过尾部标志
	}
	//cout << "Read over NO. is :" << trilis.size() << endl;
	return true;
}

void PointMove(){
	cout << "point moving..." << endl;
	vector<PointVector *> slectpoint;
	Vector3f p[3];
	Vector4f coe;
	p[0] << 149.726660 ,-45.038730 , 7.791545;
	p[1] << 175.241229 ,- 4.591933, 34.809070;
	p[2] << 162.523123 ,36.191162, 4.322560;
	coe=ccomputecoe(p);
	PointVector::CheckOut(root, &slectpoint, coe);
	cout << "show the select number : "<<slectpoint.size() << endl;
	for (int i = 0; i < slectpoint.size(); i++) {
		Vector3f a(slectpoint[i]->GiveOut()[0] * stepforward, 0, 0);
		Vector3f b = slectpoint[i]->ThisPoint();
		slectpoint[i]->NewModify(b + a);
	}
}

void OutStlDoc() {
	cout << "outputing to file..." << endl;
	FILE *fw;
	fopen_s(&fw, outputname, "wb");
	char head[80] = "This document is create by CD. in DIT.(2017 March 16)";
	fwrite(head,1,80,fw);
	fwrite(&unTriangles, 1, 4, fw);
	//unTriangles 

	for (unsigned int i = 0; i < trilis.size(); i++) {
		//trilist[i].UpdateVec();
		trilis[i].InputFile(fw);
	}
	fclose(fw);
}

vector<string>outfile;
float dmax = 0;
bool ReadASCII(const char *buffer, Vector4f cc)
{
	double  d[3];
	string s;
	stringstream ss(buffer);

	cout << "ReadASCII" << endl;

	while (getline(ss, s)) {
		if (s[0] == 'v')
		{
			int m = 2;
			for (int i = 0; i < 3; i++) {
				char st[20];
				int k = m;
				int j = 0;
				for (; (s[k] != ' ')&&(s[k]!='\0'); k++) {
					st[j] = s[k];
					j++;
				}
				m = ++k;
				d[i] = atof(st);
			}
			JustFast2(Vector3f(d[0], d[1], d[2]), cc, &select, outfile.size(),&dmax);
		}
		
		outfile.push_back(s);
	}
	return true;
}

string transff(Vector3f a) {
	char st[3][30];
	string s="v ";
	for (int i = 0; i < 3; i++) {
		sprintf(st[i], "%lf", a[i]);
		int j = 0;
		while (st[i][j] != '\0') {
			s += st[i][j];
			j++;
		}
		if (i != 2) {
			s += ' ';
		}
	}
	return s;
}

void DataReplace() {
	cout << dmax  << endl;
	cout << "Data Replacing... :select data is :" <<select.size()<<" "<< endl;
	for (int i = 0; i < select.size(); i++) {
		if (!(i % 1000)) {
			cout << "now process is : " << i << endl;
		}
		select[i].a[0]+=(select[i].s / dmax)*stepforward;
		outfile[select[i].l]= transff(select[i].a);
	}
}

void OutObjDoc() {
	cout << "Output to file :outfile is: " << outfile.size()<< endl;
	FILE *fp;
	fopen_s(&fp, outputname, "w");
	for (unsigned int i = 0; i < outfile.size() ; i++) {
		char mv[100];
		strcpy(mv, outfile[i].c_str());
		mv[outfile[i].size()] = '\n';
		fwrite(mv,1,outfile[i].size()+1,fp);
	}
	fclose(fp);
}

void ReadObjAscii() {
	/*
		三点为趾围平面定点
		Pt1 149.726660 -45.038730 7.791545
		Pt2 175.241229 -4.591933 34.809070
		Pt3 162.523123 36.191162 4.322560
	*/
	cout << "point moving..." << endl;
	vector<PointVector *> slectpoint;
	Vector3f p[3];
	Vector4f coe;
	p[0] << 149.726660, -45.038730, 7.791545;
	p[1] << 175.241229, -4.591933, 34.809070;
	p[2] << 162.523123, 36.191162, 4.322560;
	coe = ccomputecoe(p);// 计算平面方程系数

	FILE * pFile;
	long lSize;
	char* buffer;
	size_t result;

	fopen_s(&pFile, openfilename, "r");
	if (pFile == NULL)
	{
		fputs("File error", stderr);
		exit(1);
	}

	fseek(pFile, 0, SEEK_END);
	lSize = ftell(pFile);
	rewind(pFile);

	buffer = (char*)malloc(sizeof(char)*lSize);
	if (buffer == NULL)
	{
		fputs("Memory error", stderr);
		exit(2);
	}

	result = fread(buffer, 1, lSize, pFile);
	if (result != lSize)
	{
		fputs("Reading error", stderr);
		exit(3);
	}
	ios::sync_with_stdio(false);
	ReadASCII(buffer, coe);
	ios::sync_with_stdio(true);
	free(buffer);
	fclose(pFile);
}

int CAddLength() {
	/*
		纵向横切面三点：
		Pt1 28.649047  0.477465  90.453272
		Pt2 136.864537 -1.992595 157.171762
		Pt3 237.816803 1.965069  14.085081
	*/
	Vector3f p[3];
	p[0] << 28.649047, 0.477465, 90.453272;
	p[1] << 136.864537, -1.992595, 157.171762;
	p[2] << 237.816803, 1.965069, 14.085081;

	/*
		三点为趾围平面定点
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
		横切面三点;
	*/
	Vector3f p[3];
	p[0] << 63.463964,  -57.425816, 42.609455;//(28.2512570,0,91.3379820);
	p[1] << 166.733872, -124.976311, 13.108240;
	p[2] << 166.188756, -68.167633, -156.102804;
	//Vector3f origin(0,0,90.06);//标记点的原点
	QuaternionSpin a;
	a.SetQuaternionFromTwo(QuaternionSpin::VectorCompute(p[0], p[2], p[1]), Vector3f(0, 1, 0));
	a.SetCenterShift(p[0],1);//
	a.SetHeel(90.06);//设置跟高
	a.TransferSpin(&allpoint,&midpoint,a.GiveShift());
	if (a.FindPoint(&midpoint)) {
		a.TransferSpin(&midpoint,&outpoint, Vector3f(28.2512570, 0, 91.3379820));
		//ManageObj::OutFilePointObj(&outpoint, "transfershift4.obj");
	}
}

#define STEPADDLENGHT 5 //mm
#define HEELHIGHT 90.06 //mm

string increname(int i) {
	string head = "cutoutline-"; string suffix = ".obj";
	char num[3];
	sprintf(num, "%d", i);
	string acf(num);
	head += acf + suffix;
	return head;
}

int main(int argc, char* argv[])
{
	//Quaerniont();
	//GoOutLine();
	//CAddLength();

	/*
		纵向横切面三点：
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
	
	float movestep;
	vector<Vector3f> outline;
	MyOpenMesh ios;

	MyMesh::VertexHandle vertex[3];
	ios.ReadStlfile("shoestl.stl");
	ios.FindNearest(p2[0],p2[1],p2[2],vertex);
	SurfaceCoe sfc(vertex, ios.mesh);
	SurfaceCoe *meta;
	if (sfc.Init(1)) {
		ios.FindNearest(p[0], p[1], p[2], vertex);
		 meta= sfc.FindMetara(vertex[0], vertex[2]);
		
		//sfc.OutlineEigen(&gooutline);
		//ManageObj::OutFilePointObj(&gooutline, "outline4.obj");
	}
	vector<MySurCutArry> arry;
	sfc.OutlineXCoe((1+5/ meta->ReturnLength()),arry);
	vector<SurfaceCoe*> cutout(arry.size());
	//Vector3f axi = sfc.AxieCut(HEELHIGHT);
	Vector3f axi = sfc.TempVector();

	if (axi.norm()) {
		axi.normalize();
		for (int i = 0; i < arry.size(); i++) {
			cutout[i] = new SurfaceCoe(axi, arry[i].a, arry[i].x, ios.mesh);
			if (cutout[i]->Init()) {
				cutout[i]->OutlineEigen(&outline);
				cout << outline.size() << endl;
				ManageObj::OutFilePointObj(&outline, increname(i).c_str());
				outline.clear();
			}
		}
	}

	//ios.GoOutline2(p[0], p[1], p[2]);
	//ios.OutlineEigen(&outline);
	//ManageObj::OutFilePointObj(&outline, "outlinezhizhi.obj");
	/*movestep=ios.MetaraOutlineSort(5.0);
	if (movestep) {
		cout <<"Add length :"<<movestep << endl;
		if (ios.OutlineEigen(&outline)) {
			ManageObj::OutFilePointObj(&outline, "outlinezhizhi3.obj");
		}
	}*/
	//vector<Vector3f> outline;
	//MyOpenMesh ios;
	//ios.ReadStlfile("shoestl.stl");
	//ios.VertexNormals();
	////ios.MoveVertex(3);
	////ios.WriteStlfile("movelen3.stl");
	//ios.BottomVertex(&outline);
	//ManageObj::OutFilePointObj(&outline, "bottomshoedown.obj");
	system("pause");
	return 0;
}


