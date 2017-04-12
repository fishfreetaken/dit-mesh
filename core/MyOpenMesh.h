#pragma once

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PointVector.h"

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

#define ADDDIFFERENCE  0.0005

#define OUTLINEINTERVAL  10 //MM ����

#define MAXIMUMX 9999

#define CUTINTERVEL 10 //mm

#define CUTSECTION 22  //����ֺΧ���߸����ֳ�22�� 21������  ��һ�����ߺ��ԣ���arry 1 ��ʼ
#define CUTSECTION2 17  //����ֺΧ��Ь��㴦�ֳ�17�� 16������ arr �ܹ�37������   ��һ�����ߺ��ԣ���arry
  
class MyOpenMesh
{
public:
	MyOpenMesh(float a);
	MyOpenMesh() {};

	~MyOpenMesh() {};

	struct OutNoraml {
		MyMesh::Point a;
		MyMesh::Normal n;
		float x=0; //ϵ��
		float d=0;
	};
	
	MyMesh  mesh;

	void ReadStlfile(char * argg);
	void WriteStlfile(char *argg,int i);

	void SetOutlineVector(vector<Vector3f> *a) { mOutline = a; };

	void MoveVertex(float len); //���ŷ������ƶ�ָ������ for test
	void ReleaseVertexNormals();
	void BottomVertex(vector<Vector3f>*a); //��ȡЬ鸵ײ���

	float  MetaraOutlineSort(float tar); //ʵ���Ͼ�����ֺΧ���мӷ�

	void GoOutline2(Vector3f a, Vector3f b, Vector3f c);
	bool GoOutline(Vector3f a, Vector3f b, Vector3f c, vector<struct OutNoraml>*lm);

	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline ���
	
	void CrossVector(vector<struct OutNoraml> *v); //��ֺΧ�ӷ�

	void MyOpenMesh::FindMetaraPoint(vector<struct OutNoraml> *v);

	static Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a);
	static MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a);

	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, Vector3f *p);
	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, MyMesh::VertexHandle *p);
private:
	
	OpenMesh::IO::Options opt;

	float mHeelHight=0;//����

	/*��������*/
	vector<Vector3f> *mOutline;	//�������а����Լ�Ь�ͷ�ӳ�
	vector<struct OutNoraml> mOutline2;

	//vector<struct SurfaceCoe> mOutlineArray;

	//Vector3f mPointA;			//Ь���ߴ�һ��
	//Vector3f mPointB;			//�ϵװ�ƽ�����͹��һ��
	//MyMesh::Point mPointC; //�����ȡ���Ь鸼��һ��
	Vector3f mPointC;
	Vector3f mPointE; //Ь��������ϵ���ʼ�㣬Ь鸸���outline ��PointA����ĵ�(x)�����Ż����ߣ�l=(x-5)*0.4+30;
	MyMesh::VertexHandle mVertexStart; //outlineָ�����
	MyMesh::VertexHandle mVertexMid;
	MyMesh::VertexHandle mVertexEnd;
	MyMesh::VertexHandle mVertexEE;

	Vector4f mCoe; //����������ƽ�淽�̲���
	Vector3f mCoeABC;

	Vector3f mMiddleAix; //������

	vector<Vector4f> mSurfaceArray;


	void SetSurfaceEquation(Vector3f a, Vector3f b, Vector3f c);

	float DistSurface(MyMesh::Point a);

	int PointJudge(float a, float b);

	void AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b);

	int NextHalfEdgeJudge(MyMesh::HalfedgeHandle heh);

	void IterationHalfEdge(MyMesh::HalfedgeHandle heh); //���������ܶ�ջ���ڴ���䲻�������Ծ�������ʹ�õ��� ����whileѭ��

	float TotalLengh(vector<struct OutNoraml>a);
	
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

	//vector<MyMesh::VertexHandle> mShoeBottom;		//Ь鸵װ壻
	//vector<MyMesh::VertexHandle> mWristTop;		//��ڶ������֣�
};

typedef struct MyOpenMesh::OutNoraml MyOutNormal;


class SurfaceCoe {
public:
	SurfaceCoe(MyMesh::VertexHandle *a, MyMesh &b);					//����ȷ��һ��ƽ��
	SurfaceCoe(Vector3f bf, MyMesh::VertexHandle df,float x,MyMesh &b);			//һ��ֱ�߼�һ����ȷ��һ��ƽ��,�м�㼴Ϊ��ʼ��
	SurfaceCoe(MyMesh::VertexHandle start, MyMesh::Point end, MyMesh &d);	//����ȷ��һ��ƽ�棬���ȸ��������㣬��һ����ɱ䣬������ʼλ�õ��end��
	
	~SurfaceCoe() {};

	struct CutArry {
		MyMesh::VertexHandle a;
		float x=1;
	};

	bool Init(int a);
	bool Init();

	void InitTwoPoints(MyMesh::Point a, MyMesh::Point b) { mVertexMid = a; mVertexEnd = b; CoquerMidEnd();} //����һ��ֱ����һ����ȷ��һ��ƽ��

	void OutlineExpansion(float tar);

	Vector3f AxieCut(float heelhight);//���������߽�Ь鸽����и���� (�Ȼ�����У�����)  ��������������
	Vector3f TempVector(); //��ʱ���ֵ��������Ϊ������
	void OutlineXCoe(float a, vector<struct CutArry> &arryx);   //�����غ����ߵķָ�㣻

	SurfaceCoe *FindMetara(MyMesh::VertexHandle end, MyMesh::VertexHandle mid);
	
	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline ���
	void SetMidPoint(MyMesh::Point a) { mVertexMid = a; }
	float ReturnLength() { return mLength; }
	
private:
	MyMesh &mesh;

	Vector4f mCoe;
	Vector3f mCoeABC;
	MyMesh::Point mVertexStart; //outlineָ�����
	MyMesh::Point mVertexMid;
	MyMesh::Point mVertexEnd;
	MyMesh::VertexHandle mHandleBegin;

	float mX=1; //��Ҫ�������ŵı���ϵ����>1�� ������� ����Ϊ1��������>1 ��С��<1;
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
		MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);		//ָ���
		MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);   //������
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

		return (i[0] + i[1]) != 0 ? -1 : 1;  //����ֵ������ͬ -1 //�����෴ 1
	};

	void IterationHalfEdge(MyMesh::HalfedgeHandle heh) {
		//ƫ������һ����ߣ���vertexΪ���������ñ߶�Ӧ����
		MyMesh::HalfedgeHandle heh_next = heh;//=mesh.next_halfedge_handle(heh);
		MyMesh::VertexHandle vertex_i, vertex_s;
		while ((vertex_i != mHandleBegin)&&(vertex_s != mHandleBegin)) {
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
};

typedef struct SurfaceCoe::CutArry MySurCutArry;


