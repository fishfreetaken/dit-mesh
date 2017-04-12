#pragma once

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PointVector.h"

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

#define ADDDIFFERENCE  0.0005

#define OUTLINEINTERVAL  10 //MM 弧长

#define MAXIMUMX 9999

#define CUTINTERVEL 10 //mm

#define CUTSECTION 22  //将跖趾围到高跟处分成22段 21条截面  第一条切线忽略，从arry 1 开始
#define CUTSECTION2 17  //将跖趾围到鞋尖点处分成17段 16条截面 arr 总共37条截面   第一条切线忽略，从arry
  
class MyOpenMesh
{
public:
	MyOpenMesh(float a);
	MyOpenMesh() {};

	~MyOpenMesh() {};

	struct OutNoraml {
		MyMesh::Point a;
		MyMesh::Normal n;
		float x=0; //系数
		float d=0;
	};
	
	MyMesh  mesh;

	void ReadStlfile(char * argg);
	void WriteStlfile(char *argg,int i);

	void SetOutlineVector(vector<Vector3f> *a) { mOutline = a; };

	void MoveVertex(float len); //沿着法向量移动指定长度 for test
	void ReleaseVertexNormals();
	void BottomVertex(vector<Vector3f>*a); //提取鞋楦底部点

	float  MetaraOutlineSort(float tar); //实际上就是跖趾围进行加肥

	void GoOutline2(Vector3f a, Vector3f b, Vector3f c);
	bool GoOutline(Vector3f a, Vector3f b, Vector3f c, vector<struct OutNoraml>*lm);

	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline 输出
	
	void CrossVector(vector<struct OutNoraml> *v); //跖趾围加肥

	void MyOpenMesh::FindMetaraPoint(vector<struct OutNoraml> *v);

	static Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a);
	static MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a);

	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, Vector3f *p);
	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, MyMesh::VertexHandle *p);
private:
	
	OpenMesh::IO::Options opt;

	float mHeelHight=0;//跟高

	/*纵向横截面*/
	vector<Vector3f> *mOutline;	//用来进行摆正以及鞋楦头加长
	vector<struct OutNoraml> mOutline2;

	//vector<struct SurfaceCoe> mOutlineArray;

	//Vector3f mPointA;			//鞋跟高处一点
	//Vector3f mPointB;			//上底板平面的最凸出一点
	//MyMesh::Point mPointC; //纵向截取面的鞋楦尖端一点
	Vector3f mPointC;
	Vector3f mPointE; //鞋楦中轴线上的起始点，鞋楦跟高outline 离PointA最近的点(x)出沿着弧长走，l=(x-5)*0.4+30;
	MyMesh::VertexHandle mVertexStart; //outline指定起点
	MyMesh::VertexHandle mVertexMid;
	MyMesh::VertexHandle mVertexEnd;
	MyMesh::VertexHandle mVertexEE;

	Vector4f mCoe; //纵向横切面的平面方程参数
	Vector3f mCoeABC;

	Vector3f mMiddleAix; //中轴线

	vector<Vector4f> mSurfaceArray;


	void SetSurfaceEquation(Vector3f a, Vector3f b, Vector3f c);

	float DistSurface(MyMesh::Point a);

	int PointJudge(float a, float b);

	void AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b);

	int NextHalfEdgeJudge(MyMesh::HalfedgeHandle heh);

	void IterationHalfEdge(MyMesh::HalfedgeHandle heh); //迭代，可能堆栈的内存分配不够，可以尽量避免使用迭代 换用while循环

	float TotalLengh(vector<struct OutNoraml>a);
	
	/*
		将鞋楦分为若干个部分：
		1 鞋楦底部；
		2 鞋楦楦口（最顶层平面处）；
		3 跖趾围平面截取的前半前部分--鞋楦加长和鞋楦加肥（鞋楦面增加）；
		4 跖趾围和腰围之间的鞋楦面部分；
		5 腰围和背围之间的鞋楦面部分；
		6 背围以后的部分；
	*/
	//vector<MyMesh::VertexHandle> mMetaraFront;	//跖趾围前半部分（两部分，前半部分是鞋楦面，后半部分是底板），用来加长；

	//vector<MyMesh::VertexHandle> mMetaraMid;		//跖趾围与碗口最凸点与中轴线确定的平面部分；

	//vector<MyMesh::VertexHandle> mMetaraEnd;		//碗口最凸点与中轴线垂直的平面往后所有部分；

	//vector<MyMesh::VertexHandle> mShoeBottom;		//鞋楦底板；
	//vector<MyMesh::VertexHandle> mWristTop;		//碗口顶部部分；
};

typedef struct MyOpenMesh::OutNoraml MyOutNormal;


class SurfaceCoe {
public:
	SurfaceCoe(MyMesh::VertexHandle *a, MyMesh &b);					//三点确定一个平面
	SurfaceCoe(Vector3f bf, MyMesh::VertexHandle df,float x,MyMesh &b);			//一条直线加一个点确定一个平面,中间点即为起始点
	SurfaceCoe(MyMesh::VertexHandle start, MyMesh::Point end, MyMesh &d);	//三点确定一个平面，首先给出两个点，另一个点可变，给除初始位置点和end点
	
	~SurfaceCoe() {};

	struct CutArry {
		MyMesh::VertexHandle a;
		float x=1;
	};

	bool Init(int a);
	bool Init();

	void InitTwoPoints(MyMesh::Point a, MyMesh::Point b) { mVertexMid = a; mVertexEnd = b; CoquerMidEnd();} //用于一条直线与一个点确定一个平面

	void OutlineExpansion(float tar);

	Vector3f AxieCut(float heelhight);//根据中轴线将鞋楦进行切割分切 (等会儿再切！！！)  给出中轴线向量
	Vector3f TempVector(); //临时用手点出来点作为中轴线
	void OutlineXCoe(float a, vector<struct CutArry> &arryx);   //给出沿横切线的分割点；

	SurfaceCoe *FindMetara(MyMesh::VertexHandle end, MyMesh::VertexHandle mid);
	
	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline 输出
	void SetMidPoint(MyMesh::Point a) { mVertexMid = a; }
	float ReturnLength() { return mLength; }
	
private:
	MyMesh &mesh;

	Vector4f mCoe;
	Vector3f mCoeABC;
	MyMesh::Point mVertexStart; //outline指定起点
	MyMesh::Point mVertexMid;
	MyMesh::Point mVertexEnd;
	MyMesh::VertexHandle mHandleBegin;

	float mX=1; //需要进行增放的比例系数（>1） 如果不变 保持为1，增大则>1 缩小则<1;
	int mIth[3] = {0,0,0};  //start: 0,  mid  end  metara
	float mLen[3] = { 0,0,0 };
	float mLength=0;
	float mExtension = 0;

	vector<MyOutNormal> mOutline2;

	void AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b);
	void OutlineXCoe();
	void OutlineXCoe2();

	void OutlineRefine();
	void CoquerMidEnd();

	int NextHalfEdgeJudge(MyMesh::HalfedgeHandle heh) {
		MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);		//指向点
		MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);   //出发点
		float ac = PointJudge(DistSurface(mesh.point(vertex_s)), DistSurface(mesh.point(vertex_i)));
		if (ac != -1) {
			return 1;
		}
		return 0;
	};

	float DistSurface(MyMesh::Point a) {
		//Vector4f sa(a[0], a[1], a[2], 1);
		return (a[0] * mCoe[0] + a[1] * mCoe[1] + a[2] * mCoe[2] + mCoe[3]);
		//return sa.dot(mCoe); // mCoeABC.norm());
	};
	float DistPoints(MyMesh::Point a, MyMesh::Point b) {
		MyMesh::Point c(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
		return c.norm();
	}

	int PointJudge(float a, float b) {
		int i[2];
		if ((a == 0) || (b ==0)) {
			return 0;
		}
		i[0] = a < 0 ? -1 : 1;
		i[1] = b < 0 ? -1 : 1;

		return (i[0] + i[1]) != 0 ? -1 : 1;  //两个值符号相同 -1 //符号相反 1
	};

	void IterationHalfEdge(MyMesh::HalfedgeHandle heh) {
		//偏移至下一条半边，以vertex为基点索引该边对应顶点
		MyMesh::HalfedgeHandle heh_next = heh;//=mesh.next_halfedge_handle(heh);
		MyMesh::VertexHandle vertex_i, vertex_s;
		while ((vertex_i != mHandleBegin)&&(vertex_s != mHandleBegin)) {
			vertex_i = mesh.to_vertex_handle(heh_next);	//指向点
			vertex_s = mesh.from_vertex_handle(heh_next);  //出发点
			float ac = PointJudge(DistSurface(mesh.point(vertex_s)), DistSurface(mesh.point(vertex_i)));
			if (ac == 0) {
				AddOutlinePoint(vertex_i, vertex_s);
				heh_next = mesh.opposite_halfedge_handle(heh_next);
				heh_next = mesh.next_halfedge_handle(heh_next);
				while (!NextHalfEdgeJudge(heh_next)) {
					heh_next = mesh.next_halfedge_handle(heh_next);
					heh_next = mesh.opposite_halfedge_handle(heh_next);
					heh_next = mesh.next_halfedge_handle(heh_next);
				}
			}
			else if (ac == 1) {
				AddOutlinePoint(vertex_i, vertex_s);
				heh_next = mesh.opposite_halfedge_handle(heh_next);
				heh_next = mesh.next_halfedge_handle(heh_next);
			}
			else {  //-1
				heh_next = mesh.next_halfedge_handle(heh_next);
			}
		}
	};


	float TotalLengh(vector<MyOutNormal>a) {
		float s = 0;
		if (!a.size()) {
			cout << "error total length!" << endl;
			return 0;
		}
		MyMesh::Point k = a[0].a;
		for (int i = 1; i < a.size(); i++) {
			s += (a[i].a - k).norm();
			k = a[i].a;
		}
		return s;
	};
};

typedef struct SurfaceCoe::CutArry MySurCutArry;


