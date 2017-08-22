#pragma once

// -------------------- OpenMesh
//#include "OpenMesh/Core/IO/MeshIO.hh"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
//#include "PointVector.h"
#include"QuaternionSpin.h"
#include <math.h>
/*
	�ڽ��а�ߵ�����ȡ�����ߵ�ʱ�򣬷��ֵ�����Ҫ����һ�����ڴ���Դ������vs���������Դ���ޣ����Խ��д����������º����׳�����Դ���㵼�³�����������Ծ���ʹ��while����for���д���
*/
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
typedef Quaternion<float, 0> Quaternionx;

#define ADDDIFFERENCE  0.0005

#define OUTLINEINTERVAL  10 //MM ����

#define MAXIMUMX 9999

//#define CUTINTERVEL 10 //mm

//#define CUTSECTION 22  //����ֺΧ���߸����ֳ�22�� 21������  ��һ�����ߺ��ԣ���arry 1 ��ʼ
//#define CUTSECTION2 15  //����ֺΧ��Ь��㴦�ֳ�15�� 14������ arr �ܹ�35������   ��һ�����ߺ��ԣ���arry
#define CUTSECTION3  33  //32������
//#define WISTSECTION 16
//#define WISTSECTION2 12 //6

#define GAUSSIONFILTERNUM 3

#define ONEINCHLEN 25.4 //mm

#define ITERATIONCISHU 9
#define ITERATIONCISHUPOINT 3//��Χ����Χ�ĸ�˹�˲�����Χ 4

#define TOPOFFSET 7.5

#define DIFWINDOW 7 //�����󷽲�Ĵ��ڴ�С
//#define SHOEBOTTOMLL 2 //�ڵײ���ȡ�Ĺ����У���������ʹ��set����vector,vector���ܻ��ܶ�

#define ADDLENGTHSTEP   0.055 //0.01

#define MIDDLEFILTERMOVE 17 //ƽ����cutoutƽ����ʼ��

#define FILETERDIFFERMETARA 0.09

#define COEINVERWIDEN  0.5  //����ӿ�ϵ����

class SurfacePure 
{
public:
	SurfacePure(Vector3f a, Vector3f b, Vector3f c) :
		start(a),
		mid(b),
		end(c)
	{
		init();
	}
	SurfacePure(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c)
	{
		start = Vector3f(a[0], a[1], a[2]);
		mid = Vector3f(b[0], b[1], b[2]);
		end = Vector3f(c[0], c[1], c[2]);
		init();
	}
	void init() 
	{
		Vector3f ab = start - mid;
		ab = (start - end).cross(ab);
		mCoeABC = Vector3f(ab[0], ab[1], ab[2]);

		float d = ab.dot(Vector3f(0, 0, 0) - start) / mCoeABC.norm();
		mCoeABC.normalize();
		//mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);
		mCoe.data()[0] = mCoeABC[0];
		mCoe.data()[1] = mCoeABC[1];
		mCoe.data()[2] = mCoeABC[2];
		mCoe.data()[3] = d;
	}

	float DistSurface(MyMesh::Point a) {
		return (a[0] * mCoe[0] + a[1] * mCoe[1] + a[2] * mCoe[2] + mCoe[3]); //Vector4f sa(a[0], a[1], a[2], 1);//return sa.dot(mCoe); // mCoeABC.norm());
	};
private:
	Vector4f mCoe;
	Vector3f mCoeABC;

	Vector3f start;
	Vector3f mid;
	Vector3f end;
};

typedef class SurfaceCoe SurFC;
class MyOpenMesh
{
public:
	MyOpenMesh(float a);
	MyOpenMesh() {};

	~MyOpenMesh() {};

	struct OutNoraml {  //����Ҫ������class���õ���
		MyMesh::Point a; //��ʵֻҪ�����Ϳ����ˣ�����������Ժ�ʡ��������������Կ�����ģ�,��ɢ�������
		MyMesh::Normal n; //��ɢ�ķ�����
		//MyMesh::Normal nf; //��¼�õ��ԭʼ������

		float x=0;  //����ϵ��
		float d=0;	//��¼��ʼ�㵽�õ�ľ���
		MyMesh::Point f=MyMesh::Point(0,0,0); //��¼�����ԭ���ĵ�����
		MyMesh::Point m= MyMesh::Point(0, 0, 0); //��¼�ƶ�������

		float q = 0;//���������¼�õ㵽ǰƽ�����ƽ��ľ��룻0
		float h = 0;//���������¼�õ㵽��ƽ�����ƽ����룻 1
		
		//MyMesh::Point pro = MyMesh::Point(0, 0, 0); //��¼�����ı���
		//float k = 0;//��ɢ����
	};

	struct OutBottom {
		MyMesh::VertexHandle n;
		MyMesh::Point a;
		MyMesh::Normal s;
		float x = 1;
		bool operator < (const struct OutBottom &m)const {
			//return x > m.x;
			return n.idx()< m.n.idx();
		}
	};
	struct OutBottomLine {
		MyMesh::Normal s;   //������׼��--������Ŀ�ƽ��
		MyMesh::Point a;
		float x;			//��׼�������ĳ���
		int ith;			//mOutline2 ϵ�����õ�Ϊ�ڼ�����
		bool operator < (const struct OutBottomLine &m)const {
			return m.x < x; //�Ӵ�С
		}
	};
	
	MyMesh  mesh;

	void ReadStlfile(const char * argg);
	void WriteStlfile(const char *argg,int i);

	void MoveVertex(float len); //���ŷ������ƶ�ָ������ for test
	void ReleaseVertexNormals();
	void BottomVertex(vector<Vector3f>*a); //��ȡЬ鸵ײ���
	
	static Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a);
	static MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a);

	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, MyMesh::VertexHandle *p);
	MyMesh::VertexHandle FindNearest(MyMesh::Point a);

	//void ShoeExpansion(vector<SurfaceCoe *> &arr, SurfaceCoe*sfc,vector<Vector3f>&giv);
	void ShoeExpansion(vector<SurfaceCoe *> &arr, SurfaceCoe*sfc);
	void ShoeExpansion(vector<SurfaceCoe *> &arr, SurfaceCoe*sfc, SurfaceCoe*meta);
	MyMesh::Point VertexFilter(MyMesh::VertexHandle vh);//������ռ�����Ƚϴ�
	MyMesh::Point VertexFilter2(MyMesh::VertexHandle vh);//Զ����ռ�����Ƚϴ�

	void ShoeExpansionWist(vector<SurfaceCoe *> &arr);
	void ShoeExpansionWist(SurfaceCoe*meta, SurfaceCoe*metb, SurfaceCoe*metc);
	void ShoeExpansionWist2(SurfaceCoe*meta, SurfaceCoe*metb, SurfaceCoe*metc);
	void ShoeExpansionWist3(SurfaceCoe*meta, SurfaceCoe*metb, SurfaceCoe*metc,float exp);
	void ShoeExpansion(vector<SurfaceCoe*> &arr, vector<MyMesh::Point>&css); //debug
	
	void ShoeAddLength(MyMesh::Point a, SurfaceCoe*met, float exp);

	void ShoeSpin(Quaternionx &a, Vector3f b); //

	void MetaraFileter(SurfaceCoe* meta);//��������Χ����һ����ƽ���˲�

	

private:
	OpenMesh::IO::Options opt;
	void BotIteration(set<struct OutBottom>&arr, int idx, int iver); //�����ڵײ����е�����ȡ�������ٶ�̫��
	void BotIteration(set<int>&arr,int idx,int iver);
	MyMesh::Point MyOpenMesh::GaussFilter(vector<struct OutBottom>&vm, float wind);
	/*
		��Ь鸷�Ϊ���ɸ����֣�
		1 Ь鸵ײ���
		2 Ь�鸿ڣ����ƽ�洦����
		3 ��ֺΧƽ���ȡ��ǰ��ǰ����--Ь鸼ӳ���Ь鸼ӷʣ�Ь������ӣ���
		4 ��ֺΧ����Χ֮���Ь��沿�֣�
		5 ��Χ�ͱ�Χ֮���Ь��沿�֣�
		6 ��Χ�Ժ�Ĳ��֣�
	*/
	//vector<MyMesh::VertexHandle> mMetaraFront;	//��ֺΧǰ�벿�֣������֣�ǰ�벿����Ь��棬��벿���ǵװ壩�������ӳ���

	//vector<MyMesh::VertexHandle> mMetaraMid;		//��ֺΧ�������͹����������ȷ����ƽ�沿�֣�

	//vector<MyMesh::VertexHandle> mMetaraEnd;		//�����͹���������ߴ�ֱ��ƽ���������в��֣�

	//vector<MyMesh::VertexHandle> mShoeBottom;		//Ь鸵װ壻//���Ӧ������Ҫ��
	//vector<MyMesh::VertexHandle> mWristTop;		//��ڶ������֣�
	void TailGaussionFilter(vector<MyMesh::VertexHandle>& vm,int count);
	void LaplacianFilter(vector<MyMesh::VertexHandle>& vm, int count);//ƽ��ֵ
	void LaplacianFilterDist(vector<MyMesh::VertexHandle>& vm, int count);//���ݾ������ƽ��ֵ
	MyMesh::Point  GaussionArroundVertex(vector<MyMesh::Point>& sm,float sigma);

	void getFaceCentroid(std::vector<MyMesh::Point> &centroid);
	void updateVertexPosition(
		std::vector<MyMesh::Normal> &filtered_normals, 
		int iteration_number, bool fixed_boundary);
	void updateVertexPosition( vector<MyMesh::VertexHandle>&vmv,
		std::vector<MyMesh::Normal> &filtered_normals,
		int iteration_number, bool fixed_boundary);

};

typedef struct MyOpenMesh::OutBottom MyOutBottom;
typedef struct MyOpenMesh::OutNoraml MyOutNormal;
typedef struct MyOpenMesh::OutBottomLine MyBotOutLine;

class SurfaceCoe 
{
public:
	SurfaceCoe(MyMesh::VertexHandle *a, MyMesh &b);					//����ȷ��һ��ƽ��
	SurfaceCoe(MyMesh::Point mid, MyMesh::Point end, MyMesh::VertexHandle c, MyMesh &d);  //������ȷ��һ��ƽ��
	SurfaceCoe(Vector3f bf, MyMesh::VertexHandle df, float x, MyMesh &b);			//һ��ֱ�߼�һ����ȷ��һ��ƽ��,�м�㼴Ϊ��ʼ��
	SurfaceCoe(MyMesh::VertexHandle start, MyMesh::Point end, MyMesh &d);	//����ȷ��һ��ƽ�棬���ȸ��������㣬��һ����ɱ䣬������ʼλ�õ��end��
	~SurfaceCoe() {};

	struct CutArry 
	{
		MyMesh::VertexHandle a;
		float x;// = 1;
		Vector3f n; //�������������ʾ�ý���ķ�������
		CutArry()
		{
			x = 1;
		}
	};

	bool Init(int cmd);
	bool Init();
	void InitTwoPoints(MyMesh::Point a, MyMesh::Point b) { mVertexMid = a; mVertexEnd = b; CoquerMidEnd(); } //����һ��ֱ����һ����ȷ��һ��ƽ��

	Vector3f AxieCut(float heelhight);//���������߽�Ь鸽����и���� (�Ȼ�����У�����)  ��������������
	

	//void OutCutOutline(float a, vector<struct CutArry> &arryx);   //�����غ����ߵķָ�㣻
	vector<struct CutArry> OutCutOutline(float exp,SurfaceCoe *a, float h,int*imet);   //�����غ����ߵķָ�㣻
	vector<struct CutArry> OutCutOutline(float exp, SurfaceCoe *a,float h,int& imet);   //�����غ����ߵķָ�㣻�Զ����������ߵĵ�

	SurfaceCoe *FindMetara(MyMesh::VertexHandle end, MyMesh::VertexHandle mid); //������ʼ��end��,���������ߴ�����Ѱ��
	SurfaceCoe* FindWaistLine(SurfaceCoe *met);
	SurfaceCoe* FindToeBottomPoint(float a);//���������ֺΧ��Ь鸵װ���������루ab-4��gb-4��hb-4������λֺΧ�����꣬������Բο��ĵ�˵��
	SurfaceCoe* SfcMoveXLen(SurfaceCoe *toe, float x);

	float FindAddLenth(SurfaceCoe *met,float ext);

	int UpOneInch(int ith, MyMesh::Point &a);
	
	MyMesh::VertexHandle FindNearest(MyMesh::Point a);
	
	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline ���
	bool OutlineEigenM(vector<Vector3f> *a); //output vector<Vector3f> outline ���
	bool OutlineEigenaf(vector<Vector3f> *a); //output vector<Vector3f> outline ���
	bool OutlineEigen(vector<Vector4f> *a);
	bool OutlineEigenf(const char *a);
	void SetMidPoint(MyMesh::Point a) { mVertexMid = a; }
	MyMesh::VertexHandle PointExchange(MyMesh::Point a){
		if (mVertexStart[1] > 0) {
			mVertexMid = mVertexStart;
			mVertexEnd = mVertexEnd;
		}
		else {
			mVertexMid = mVertexEnd;
			mVertexEnd = mVertexStart;
		}
		if (mVertexStart[1] == 0) {
			cout << "Point Exchange Error!" << endl;
		}
		MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
		MyMesh::VertexIter	v_it_s;
		float min = MAXIMUMX;
		float lin;
		MyMesh::Point pp;
		for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		{
			pp = mesh.point(*v_it);
			lin = (pp - a).norm();
			if (lin < min) {
				min = lin;
				v_it_s = v_it;
			}
		}
		mVertexStart = mesh.point(*v_it_s);
		mHandleBegin = mesh.vertex_handle(v_it_s->idx());
		Init(0);
		return mHandleBegin;
	}
	
	float AllocateXCoe();
	float AllocateXCoe(float *a);		//��ʼ������ϵ�����ӵױ���Ϊ��ʼ��  �޸ĺ����float���е���
	float AllocateXCoe(float &a);
	float AllocateXCoe(float a, float c); // (SurfaceCoe*met, MyMesh::Point a, MyMesh::Point b);			//��ʼ������ϵ�������ֲ�ͬ���ͣ��Ӻ���������������Ϊ��ʼ��
	float AllocateXCoe(SurfacePure*met);  //��������������
	float TopSlide(SurfacePure*met);//�������ж����ƽ������


	float DistSurface(MyMesh::Point a) {
		return (a[0] * mCoe[0] + a[1] * mCoe[1] + a[2] * mCoe[2] + mCoe[3]); //Vector4f sa(a[0], a[1], a[2], 1);//return sa.dot(mCoe); // mCoeABC.norm());
	};

	int FindNearestOutline(MyMesh::Point a);

	float ReturnExtension() { return mExtension; }
	float ReturnLength() { 	return mLength;}
	Vector3f ReturnCoe() { return mCoeABC; }
	MyMesh::Point ReturnStartPoint() { return mVertexStart; }
	MyMesh::Point ReturnMidPoint() { return mVertexMid; }
	MyMesh::Point ReturnEndPoint() { return mVertexEnd; }

	int ReturnIth(int i) { return mIth[i]; }
	void SetMIth(int i) { mIth[2] = i; }

	void SetMIth(SurfaceCoe*sfc) {
		mIth[2]=sfc->FindNearestOutline(mVertexStart);
	}

	MyMesh::VertexHandle ReturnVertexHandle() { return mHandleBegin; }
	void ReturnTriPoint(MyMesh::VertexHandle *vr, MyOpenMesh&ios) {
		ios.FindNearest(mVertexStart,mVertexMid,mVertexEnd, vr);
	}
	
	void InitMidEndPoint(vector<MyMesh::Point>&a);
	
	void InitMidEndPoint(MyMesh::Point&mid, MyMesh::Point&end) { 
		vector<MyMesh::Point>mv; 
		InitMidEndPoint(mv); 
		mid = mv[0]; 
		end = mv[1]; 
	}
	void InitMidEndPoint() {
		vector<MyMesh::Point>mv;
		InitMidEndPoint(mv);
		InitTwoPoints(mv[0], mv[1]);
	}

	int InitMidTopPoint(vector<MyMesh::Point>&a);
	void InitMidTopPoint(vector<MyMesh::Point>&fw, float x);

	void CalculateLen() { mLength = TotalLengh(mOutline2); }
	int ReturnMoutline2Len() { return mOutline2.size(); }

	void CoquerMidEnd();		//Ѱ����ֵ���Լ���ֵ��

	MyMesh::Point FindNearestPoint(MyMesh::Point a, float &s);
	MyMesh::Point FindNearestPoint(MyMesh::Point a, float &s, int i,float m);//���int i�����ж� 0ǰƽ�滹�� 1��ƽ��ģ�
	//MyMesh::Point FindNearestPoint(MyMesh::VertexHandle a, float &s);
	MyMesh::Point FindNearestPoint(MyMesh::Point a, MyMesh::Normal nf, float &s);

	MyMesh::Point ReturnSpecificPoint(int i) {
		if (i >= mOutline2.size()) {
			cout << "Specific error Point!" << endl;
			return MyMesh::Point(0, 0, 0);
		}
		return mOutline2[i].a;
	}

	SurfacePure* LastBottomSurf(SurfaceCoe* sfc) {
		if (!mOutline2.size()) {
			cout << "SurfacePure is zeror" << endl;
			return NULL;
		}
		//MyMesh::Point ssp = mOutline2[(mIth[0] + mIth[1]) / 2].a;

		//int iif=sfc->FindNearestOutline(ssp);
		//iif -= 7;//����������ľ���

		//ssp = sfc->ReturnSpecificPoint(iif);

		MyMesh::Point ssp=sfc->ReturnSpecificPoint(sfc->ReturnIth(1));
		//cout << (mIth[0] + mIth[1]) / 2 << endl;
		//cout << ssp << endl;
		SurfacePure *met = new SurfacePure(mOutline2[mIth[0]].a, mOutline2[mIth[1]].a,ssp );

		return met;
	}

	void initQH(SurfacePure *f, int ii) {  //give up!
		float ss = 0;
		for (int i = 0; i < mOutline2.size(); i++) {
			ss = f->DistSurface(mOutline2[i].a);
			if (ii) {
				mOutline2[i].h = ss;
			}
			else {
				mOutline2[i].q = ss;
			}
		}
	}

	float ReturnExtensionli() {
		return mExtensionli;
	}

private:
	MyMesh &mesh;

	Vector4f mCoe;
	Vector3f mCoeABC;
	MyMesh::Point mVertexStart; //outlineָ�����
	MyMesh::Point mVertexMid;
	MyMesh::Point mVertexEnd;
	MyMesh::VertexHandle mHandleBegin;

	float mX = 1; //��Ҫ�������ŵı���ϵ����>1�� ������� ����Ϊ1��������>1 ��С��<1;
	int mIth[3];// = { 0 };  //start: 0,  mid  end  (metara ��������̫��Ҫ) ����������������sfc�������������ʼ��λ��
	float mLen[3];// = { 0 }; //start-mid  mid-end  end-start distance
	float mLength = 0;
	float mExtension = 0;
	float mExtensionli = 0; //for debug;

	std::vector<MyOutNormal> mOutline2;
	void AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b);//���outline����ĵ�

	void OutlineRefine();		//ƽ��outline
	
	float OutlineExpansion();	//������Χ��������Ӧ��allocatexcoe����
	float OutlineExpansion(float *a);
	float OutlineExpansion(float a);
	//void InitMidEndPoint();

	int NextHalfEdgeJudge(MyMesh::HalfedgeHandle heh) {
		MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);		//ָ���
		MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);   //������
		float ac = PointJudge(DistSurface(mesh.point(vertex_s)), DistSurface(mesh.point(vertex_i)));
		if (ac != -1) {
			return 1;
		}
		return 0;
	};
	Vector3f NextHalfEdgeJudge2(MyMesh::HalfedgeHandle heh) {
		MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);		//ָ���
		MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);   //������
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

		return (i[0] + i[1]) != 0 ? -1 : 1;  //����ֵ������ͬ -1 //�����෴ 1
	};

	void IterationHalfEdge(MyMesh::HalfedgeHandle heh) {
		//ƫ������һ����ߣ���vertexΪ���������ñ߶�Ӧ����
		MyMesh::HalfedgeHandle heh_next = heh;//=mesh.next_halfedge_handle(heh);
		MyMesh::VertexHandle vertex_i, vertex_s;
		while ((vertex_i != mHandleBegin) && (vertex_s != mHandleBegin)) {
			vertex_i = mesh.to_vertex_handle(heh_next);	//ָ���
			vertex_s = mesh.from_vertex_handle(heh_next);  //������
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
		int nCenter;//�������ĵ�  
		double dDis;//������һ�㵽���ĵ����  

		//�м����  
		double dValue;
		double dSum;
		dSum = 0;

		*pnWidowSize = 1 + 2 * ceil(3 * sigma);// [-3*sigma,3*sigma] �������ݣ��Ḳ�Ǿ��󲿷��˲�ϵ��  

		nCenter = (*pnWidowSize) / 2;
		//���ɸ�˹����  
		for (int i = 0; i<(*pnWidowSize); i++)
		{
			dDis = (double)(i - nCenter);
			dValue = exp(-(1 / 2)*dDis*dDis / (sigma*sigma)) / (sqrt(2 * 3.1415926)*sigma);
			(*pdKernel)[i] = dValue;
			dSum += dValue;
		}
		//��һ��  
		for (int i = 0; i<(*pnWidowSize); i++)
		{
			(*pdKernel)[i] /= dSum;
		}
	}

	//�ø�˹�˲���ƽ��
	//pGray ���룻pResult �����
	void GaussianSmooth(vector<float>&pGray, vector<float>&pResult, double sigma)
	{
		int x, y,i;
		int sz = pGray.size();

		int nWindowSize;//��˹�˲�������  

		int nLen;//���ڳ���  

		double *pdKernel;//һά��˹�˲���  

		double dDotMul;//��˹ϵ����ͼ�����ݵĵ��  

		double dWeightSum;//�˲�ϵ���ܺ�  

		nWindowSize = 1 + 2 * ceil(3 * sigma);// [-3*sigma,3*sigma] �������ݣ��Ḳ�Ǿ��󲿷��˲�ϵ��  

		if ((pdKernel = (double *)malloc(nWindowSize * sizeof(double))) == NULL)
		{
			printf("malloc memory for pdKernel,failed!!");
			exit(0);
		}

		//����һά��˹����  
		CreatGauss(sigma, &pdKernel, &nWindowSize);

		nLen = nWindowSize / 2;

		//x�����˲�  
		for (int x = 0; x<sz; x++)
		{
			dDotMul = 0;
			dWeightSum = 0;
			for (int i = (-nLen); i <= nLen; i++)
			{
				//�ж��Ƿ���ͼ���ڲ�  
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
	//void AxieSurfacePoint(MyMesh::Point a, MyMesh::Point b, MyMesh::Point cs); //�������ֱ���ϵ�ͶӰ
	float CrossPointAxi(MyMesh::Point a,MyMesh::Point b) 
	{
		//float cc = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (a.norm()*b.norm());
		float cc = (a|b) / (a.norm()*b.norm());
		return cc;
	}

	/*
		Ѱ��outline������ĵ����ţ�
	*/
	int findOutlineNearestIth(MyMesh::Point m)
	{
		int t=0;
		float ini_dist = (m - mOutline2[0].a).norm();
		float min_dist = ini_dist;
		for (int i = 1; i < mOutline2.size(); i++) 
		{
			ini_dist=(m - mOutline2[i].a).norm();
			if (ini_dist<min_dist) {
				min_dist = ini_dist;
				t = i;
			}
		}
		//findPoint = mOutline2[t].a;
		return t;
	}

	float arcLengthFromStartPoint(int axi) {
		if ((axi < 0) || ((axi+1) > mOutline2.size())) {
			cout << "Number beyong Erro" << endl;
			return 0;
		}
		float len = 0;

		for (int i = 1; i <=axi; i++) {
			len+=(mOutline2[i].a- mOutline2[i-1].a).norm();
		}
		return len;
	}
	
};

typedef struct SurfaceCoe::CutArry MySurCutArry;

struct SUFACECOETOE { //�������������Χ�ȵĳ�ʼ������һ�������Ľṹ�壻
	SurfaceCoe *toe;
	MyMesh::VertexHandle vertex_handle[3];
	float pos = 0; //ֺΧ�������λ��

	float len = 0;
};

