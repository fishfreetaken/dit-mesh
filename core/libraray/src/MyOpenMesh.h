#pragma once

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PointVector.h"

/*
	在进行半边迭代求取轮廓线的时候，发现迭代需要消耗一定的内存资源，由于vs所分配的资源有限，所以进行大迭代的情况下很容易出现资源不足导致程序崩溃，所以尽量使用while或者for进行处理
*/
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
typedef Quaternion<float, 0> Quaternionx;

#define ADDDIFFERENCE  0.0005

#define OUTLINEINTERVAL  10 //MM 弧长

#define MAXIMUMX 9999

#define CUTINTERVEL 10 //mm

#define CUTSECTION 22  //将跖趾围到高跟处分成22段 21条截面  第一条切线忽略，从arry 1 开始
#define CUTSECTION2 15  //将跖趾围到鞋尖点处分成15段 14条截面 arr 总共35条截面   第一条切线忽略，从arry
#define CUTSECTION3  33  //32个截面
#define WISTSECTION 8
#define WISTSECTION2 5

#define GAUSSIONFILTERNUM 3

#define ONEINCHLEN 25.4 //mm

typedef class SurfaceCoe SurFC;
class MyOpenMesh
{
public:
	MyOpenMesh(float a);
	MyOpenMesh() {};

	~MyOpenMesh() {};

	struct OutNoraml {  //很重要。两个class都用到了
		MyMesh::Point a; //其实只要增量就可以了，这个可以在以后省掉，这个用来调试看结果的！,扩散后的坐标
		MyMesh::Normal n; //扩散的法向量
		float x=0;  //扩增系数
		float d=0;	//记录起始点到该点的距离
		MyMesh::Point f=MyMesh::Point(0,0,0); //记录相对于原来的递增量
		//MyMesh::Point pro = MyMesh::Point(0, 0, 0); //记录增长的比率
		//float k = 0;//扩散增量
	};

	struct OutBottom {
		MyMesh::VertexHandle n;
		MyMesh::Point a;
		MyMesh::Normal s;
		//int i=0;
		float x = 0;
		//bool operator < (const struct OutBottom &m)const {
		//	return x > m.x;
		//	//return i < a.i;
		//}
		bool operator < (const struct OutBottom &m)const {
			return x > m.x;
			//return i < a.i;
		}
	};
	
	MyMesh  mesh;

	void ReadStlfile(char * argg);
	void WriteStlfile(char *argg,int i);

	void MoveVertex(float len); //沿着法向量移动指定长度 for test
	void ReleaseVertexNormals();
	void BottomVertex(vector<Vector3f>*a); //提取鞋楦底部点

	static Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a);
	static MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a);

	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, MyMesh::VertexHandle *p);
	MyMesh::VertexHandle FindNearest(MyMesh::Point a);

	void ShoeExpansion(vector<SurfaceCoe *> &arr);
	void ShoeExpansionWist(vector<SurfaceCoe *> &arr);
	vector<MyMesh::Point> ShoeExpansion(vector<SurfaceCoe*> &arr, SurfaceCoe* met, vector<MyMesh::Point>&css,int ith); //debug

	set<struct OutBottom> ShoeBottomLine(vector<MyMesh::Point>&a);

	int TotalVertex() {
		return mesh.n_vertices();
	}
private:
	
	OpenMesh::IO::Options opt;
	
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

	//vector<MyMesh::VertexHandle> mShoeBottom;		//鞋楦底板；//这个应该是需要的
	//vector<MyMesh::VertexHandle> mWristTop;		//碗口顶部部分；
};

typedef struct MyOpenMesh::OutBottom MyOutBottom;
typedef struct MyOpenMesh::OutNoraml MyOutNormal;


class SurfaceCoe {
public:
	SurfaceCoe(MyMesh::VertexHandle *a, MyMesh &b);					//三点确定一个平面
	SurfaceCoe(MyMesh::Point mid, MyMesh::Point end, MyMesh::VertexHandle c, MyMesh &d);  //三个点确定一个平面
	SurfaceCoe(Vector3f bf, MyMesh::VertexHandle df, float x, MyMesh &b);			//一条直线加一个点确定一个平面,中间点即为起始点
	SurfaceCoe(MyMesh::VertexHandle start, MyMesh::Point end, MyMesh &d);	//三点确定一个平面，首先给出两个点，另一个点可变，给除初始位置点和end点
	
	~SurfaceCoe() {};

	struct CutArry {
		MyMesh::VertexHandle a;
		float x = 1;
		Vector3f n; //这个可以用来表示该界面的法向向量
	};

	bool Init(int cmd);
	bool Init();
	void InitTwoPoints(MyMesh::Point a, MyMesh::Point b) { mVertexMid = a; mVertexEnd = b; CoquerMidEnd(); } //用于一条直线与一个点确定一个平面

	Vector3f AxieCut(float heelhight);//根据中轴线将鞋楦进行切割分切 (等会儿再切！！！)  给出中轴线向量
	Vector3f TempVector(); //临时用手点出来点作为中轴线
	//void OutCutOutline(float a, vector<struct CutArry> &arryx);   //给出沿横切线的分割点；
	vector<struct CutArry> OutCutOutline(float exp,SurfaceCoe *a, Vector3f axi);   //给出沿横切线的分割点；
	vector<struct CutArry> OutCutOutline(float exp, SurfaceCoe *meta, SurfaceCoe *metb, SurfaceCoe *metc);//腰围增加给出横切面

	SurfaceCoe *FindMetara(MyMesh::VertexHandle end, MyMesh::VertexHandle mid); //给出起始和end点,沿着中轴线处进行寻找
	SurfaceCoe* FindWaistLine(SurfaceCoe *met);

	int UpOneInch(int ith, MyMesh::Point &a);
	MyMesh::VertexHandle FindNearest(MyMesh::Point a);

	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline 输出
	bool OutlineEigen(vector<Vector4f> *a);
	bool OutlineEigenf(vector<Vector3f> *a);
	void SetMidPoint(MyMesh::Point a) { mVertexMid = a; }
	
	float AllocateXCoe(float a);		//初始化增量系数，从底边作为起始点  修改后加入float进行调试
	float AllocateXCoe();// (SurfaceCoe*met, MyMesh::Point a, MyMesh::Point b);			//初始化增量系数，两种不同类型，从横切轮廓中轴点出作为起始点

	float DistSurface(MyMesh::Point a) {
		return (a[0] * mCoe[0] + a[1] * mCoe[1] + a[2] * mCoe[2] + mCoe[3]); //Vector4f sa(a[0], a[1], a[2], 1);//return sa.dot(mCoe); // mCoeABC.norm());
	};
	MyMesh::Point FindNearestPoint(MyMesh::Point a, float &s);

	float ReturnExtension() { return mExtension; }
	float ReturnLength() { return mLength; }
	Vector3f ReturnCoe() { return mCoeABC; }
	MyMesh::Point ReturnStartPoint() { return mVertexStart; }
	int ReturnIth(int i) { return mIth[i]; }
	void SetMIth(int i) { mIth[2] = i; }
	MyMesh::VertexHandle ReturnVertexHandle() { return mHandleBegin; }
private:
	MyMesh &mesh;

	Vector4f mCoe;
	Vector3f mCoeABC;
	MyMesh::Point mVertexStart; //outline指定起点
	MyMesh::Point mVertexMid;
	MyMesh::Point mVertexEnd;
	MyMesh::VertexHandle mHandleBegin;

	float mX = 1; //需要进行增放的比例系数（>1） 如果不变 保持为1，增大则>1 缩小则<1;
	int mIth[3] = { 0,0,0 };  //start: 0,  mid  end  (metara 保留，不太必要) 第三个用来保存在sfc主横截轮廓的起始点位置
	float mLen[3] = { 0,0,0 };
	float mLength = 0;
	float mExtension = 0;
	float mExtensionli = 0; //for debug;

	vector<MyOutNormal> mOutline2;
	void AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b);//添加outline主体的点

	void OutlineRefine();		//平滑outline
	void CoquerMidEnd();		//寻找中值点以及终值点
	float OutlineExpansion();	//变形扩围，放入相应的allocatexcoe里面

	int NextHalfEdgeJudge(MyMesh::HalfedgeHandle heh) {
		MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);		//指向点
		MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);   //出发点
		float ac = PointJudge(DistSurface(mesh.point(vertex_s)), DistSurface(mesh.point(vertex_i)));
		if (ac != -1) {
			return 1;
		}
		return 0;
	};
	Vector3f NextHalfEdgeJudge2(MyMesh::HalfedgeHandle heh) {
		MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);		//指向点
		MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);   //出发点
		MyMesh::Point as = mesh.point(vertex_s);
		MyMesh::Point ai = mesh.point(vertex_i);
		if (PointJudge(DistSurface(as), DistSurface(ai)) != -1) {
			Vector3f cac = MyOpenMesh::EigenTransfer(ai);// (ca[0], ca[1], ca[2]);
			Vector3f cbc = MyOpenMesh::EigenTransfer(as);// (cb[0], cb[1], cb[2]);
			Vector3f v = cac - cbc;
			float t;
			t = (0 - (mCoeABC.dot(cac) + mCoe[3])) / (mCoeABC.dot(v));
			return Vector3f(v[0] * t + cac[0], v[1] * t + cac[1], v[2] * t + cac[2]);
		}
		return Vector3f(9999, 9999, 9999);
	};

	float DistPoints(MyMesh::Point a, MyMesh::Point b) {
		MyMesh::Point c(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
		return c.norm();
	}

	int PointJudge(float a, float b) {
		int i[2];
		if ((a == 0) || (b == 0)) {
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
		while ((vertex_i != mHandleBegin) && (vertex_s != mHandleBegin)) {
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

	void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)
	{
		int nCenter;//数组中心点  
		double dDis;//数组中一点到中心点距离  

		//中间变量  
		double dValue;
		double dSum;
		dSum = 0;

		*pnWidowSize = 1 + 2 * ceil(3 * sigma);// [-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数  

		nCenter = (*pnWidowSize) / 2;
		//生成高斯数据  
		for (int i = 0; i<(*pnWidowSize); i++)
		{
			dDis = (double)(i - nCenter);
			dValue = exp(-(1 / 2)*dDis*dDis / (sigma*sigma)) / (sqrt(2 * 3.1415926)*sigma);
			(*pdKernel)[i] = dValue;
			dSum += dValue;
		}
		//归一化  
		for (int i = 0; i<(*pnWidowSize); i++)
		{
			(*pdKernel)[i] /= dSum;
		}
	}

	//用高斯滤波器平滑
	//pGray 输入；pResult 输出；
	void GaussianSmooth(vector<float>&pGray, vector<float>&pResult, double sigma)
	{
		int x, y,i;
		int sz = pGray.size();

		int nWindowSize;//高斯滤波器长度  

		int nLen;//窗口长度  

		double *pdKernel;//一维高斯滤波器  

		double dDotMul;//高斯系数与图像数据的点乘  

		double dWeightSum;//滤波系数总和  

		nWindowSize = 1 + 2 * ceil(3 * sigma);// [-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数  

		if ((pdKernel = (double *)malloc(nWindowSize * sizeof(double))) == NULL)
		{
			printf("malloc memory for pdKernel,failed!!");
			exit(0);
		}

		//产生一维高斯数据  
		CreatGauss(sigma, &pdKernel, &nWindowSize);

		nLen = nWindowSize / 2;

		//x方向滤波  
		for (int x = 0; x<sz; x++)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for (int i = (-nLen); i <= nLen; i++)
			{
				//判断是否在图像内部  
				if ((i + x) >= 0 && (i + x)<sz)
				{
					dDotMul += (double)(pGray[i + x] * pdKernel[nLen + i]);
					dWeightSum += pdKernel[nLen + i];
				}
			}
			pResult.push_back(dDotMul / dWeightSum);
		}
		free(pdKernel);
	}

	struct cmp {
		float a;
		int i;
		bool operator <(const struct cmp  &b) const
		{
			return b.a > a;
		}
	};
	//void AxieSurfacePoint(MyMesh::Point a, MyMesh::Point b, MyMesh::Point cs); //计算点在直线上的投影

};

typedef struct SurfaceCoe::CutArry MySurCutArry;

