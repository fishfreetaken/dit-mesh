#pragma once

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "PointVector.h"

/*
	�ڽ��а�ߵ�����ȡ�����ߵ�ʱ�򣬷��ֵ�����Ҫ����һ�����ڴ���Դ������vs���������Դ���ޣ����Խ��д����������º����׳�����Դ���㵼�³�����������Ծ���ʹ��while����for���д���
*/
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
typedef Quaternion<float, 0> Quaternionx;

#define ADDDIFFERENCE  0.0005

#define OUTLINEINTERVAL  10 //MM ����

#define MAXIMUMX 9999

#define CUTINTERVEL 10 //mm

#define CUTSECTION 22  //����ֺΧ���߸����ֳ�22�� 21������  ��һ�����ߺ��ԣ���arry 1 ��ʼ
#define CUTSECTION2 15  //����ֺΧ��Ь��㴦�ֳ�15�� 14������ arr �ܹ�35������   ��һ�����ߺ��ԣ���arry
#define CUTSECTION3  33  //32������

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
		float x=0;  //����ϵ��
		float d=0;	//��¼��ʼ�㵽�õ�ľ���
		MyMesh::Point f=MyMesh::Point(0,0,0); //��¼�����ԭ���ĵ�����
		//MyMesh::Point pro = MyMesh::Point(0, 0, 0); //��¼�����ı���
		//float k = 0;//��ɢ����
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

	void MoveVertex(float len); //���ŷ������ƶ�ָ������ for test
	void ReleaseVertexNormals();
	void BottomVertex(vector<Vector3f>*a); //��ȡЬ鸵ײ���

	static Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a);
	static MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a);

	void FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, MyMesh::VertexHandle *p);

	void ShoeExpansion(vector<SurfaceCoe *> &arr, SurfaceCoe* met);
	void ShoeExpansion(vector<SurfaceCoe*> &arrx, SurfaceCoe* met, int ith);
	vector<MyMesh::Point> ShoeExpansion(vector<SurfaceCoe*> &arr, SurfaceCoe* met, vector<MyMesh::Point>&css); //debug

	set<struct OutBottom> ShoeBottomLine(vector<MyMesh::Point>&a);

	int TotalVertex() {
		return mesh.n_vertices();
	}
private:
	
	OpenMesh::IO::Options opt;
	
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
};

typedef struct MyOpenMesh::OutBottom MyOutBottom;
typedef struct MyOpenMesh::OutNoraml MyOutNormal;


class SurfaceCoe {
public:
	SurfaceCoe(MyMesh::VertexHandle *a, MyMesh &b);					//����ȷ��һ��ƽ��
	SurfaceCoe(Vector3f bf, MyMesh::VertexHandle df, float x, MyMesh &b);			//һ��ֱ�߼�һ����ȷ��һ��ƽ��,�м�㼴Ϊ��ʼ��
	SurfaceCoe(MyMesh::VertexHandle start, MyMesh::Point end, MyMesh &d);	//����ȷ��һ��ƽ�棬���ȸ��������㣬��һ����ɱ䣬������ʼλ�õ��end��

	~SurfaceCoe() {};

	struct CutArry {
		MyMesh::VertexHandle a;
		float x = 1;
		Vector3f n; //�������������ʾ�ý���ķ�������
	};

	bool Init(int cmd);
	bool Init();
	void InitTwoPoints(MyMesh::Point a, MyMesh::Point b) { mVertexMid = a; mVertexEnd = b; CoquerMidEnd(); } //����һ��ֱ����һ����ȷ��һ��ƽ��

	Vector3f AxieCut(float heelhight);//���������߽�Ь鸽����и���� (�Ȼ�����У�����)  ��������������
	Vector3f TempVector(); //��ʱ���ֵ��������Ϊ������
	void OutCutOutline(float a, vector<struct CutArry> &arryx);   //�����غ����ߵķָ�㣻
	vector<struct CutArry> OutCutOutline(float exp,SurfaceCoe *a, Vector3f axi,int &ith);   //�����غ����ߵķָ�㣻

	SurfaceCoe *FindMetara(MyMesh::VertexHandle end, MyMesh::VertexHandle mid); //������ʼ��end��,���������ߴ�����Ѱ��

	bool OutlineEigen(vector<Vector3f> *a); //output vector<Vector3f> outline ���
	bool OutlineEigen(vector<Vector4f> *a);
	void SetMidPoint(MyMesh::Point a) { mVertexMid = a; }
	
	float AllocateXCoe(float a);		//��ʼ������ϵ�����ӵױ���Ϊ��ʼ��  �޸ĺ����float���е���
	float AllocateXCoe2();			//��ʼ������ϵ�������ֲ�ͬ���ͣ��Ӻ���������������Ϊ��ʼ��

	float DistSurface(MyMesh::Point a) {
		return (a[0] * mCoe[0] + a[1] * mCoe[1] + a[2] * mCoe[2] + mCoe[3]); //Vector4f sa(a[0], a[1], a[2], 1);//return sa.dot(mCoe); // mCoeABC.norm());
	};
	MyMesh::Point FindNearestPoint(MyMesh::Point a, float &s);

	float ReturnExtension() { return mExtension; }
	float ReturnLength() { return mLength; }
	Vector3f ReturnCoe() { return mCoeABC; }
private:
	MyMesh &mesh;

	Vector4f mCoe;
	Vector3f mCoeABC;
	MyMesh::Point mVertexStart; //outlineָ�����
	MyMesh::Point mVertexMid;
	MyMesh::Point mVertexEnd;
	MyMesh::VertexHandle mHandleBegin;

	float mX = 1; //��Ҫ�������ŵı���ϵ����>1�� ������� ����Ϊ1��������>1 ��С��<1;
	int mIth[3] = { 0,0,0 };  //start: 0,  mid  end  metara
	float mLen[3] = { 0,0,0 };
	float mLength = 0;
	float mExtension = 0;

	vector<MyOutNormal> mOutline2;
	void AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b);//���outline����ĵ�

	void OutlineRefine();		//ƽ��outline
	void CoquerMidEnd();		//Ѱ����ֵ���Լ���ֵ��
	float OutlineExpansion();	//������Χ��������Ӧ��allocatexcoe����

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
};

typedef struct SurfaceCoe::CutArry MySurCutArry;

