#pragma once
#include "PointVector.h"
#include "MyOpenMesh.h"

struct VertexTri {
	uint i;  //���
	uint p[3];  //����������������
	uint v[3];  //�����������ڵ����������α��
	Vector3f a; //����ķ�����
};
typedef struct VertexTri VertexTriMesh;

class ManageObj
{
public:
	ManageObj(char * name, vector<Vector3f>*p, vector<Vector3i>*t);
	ManageObj(char * name, vector<Vector3f>*p);
	ManageObj(char *name);

	class BaseHandle
	{
	public:

		explicit BaseHandle(int _idx = -1) : idx_(_idx) {}

		/// Get the underlying index of this handle
		int idx() const { return idx_; }

		/// The handle is valid iff the index is not equal to -1.
		bool is_valid() const { return idx_ != -1; }

		/// reset handle to be invalid
		void reset() { idx_ = -1; }
		/// reset handle to be invalid
		void invalidate() { idx_ = -1; }

		bool operator==(const BaseHandle& _rhs) const {
			return (this->idx_ == _rhs.idx_);
		}

		bool operator!=(const BaseHandle& _rhs) const {
			return (this->idx_ != _rhs.idx_);
		}

		bool operator<(const BaseHandle& _rhs) const {
			return (this->idx_ < _rhs.idx_);
		}


		// this is to be used only by the iterators
		void __increment() { ++idx_; }
		void __decrement() { --idx_; }

		void __increment(int amount) { idx_ += amount; }
		void __decrement(int amount) { idx_ -= amount; }

	private:

		int idx_;
	};

	struct VertexHandle : public BaseHandle
	{
		explicit VertexHandle(int _idx = -1) : BaseHandle(_idx) {}
	};
	class CmpVec
	{
		public:

			CmpVec(float _eps = FLT_MIN) : eps_(_eps) {}

			bool operator()(const Vector3f& _v0, const Vector3f& _v1) const
			{
				if (fabs(_v0[0] - _v1[0]) <= eps_)
				{
					if (fabs(_v0[1] - _v1[1]) <= eps_)
					{
						return (_v0[2] < _v1[2] - eps_);
					}
					else return (_v0[1] < _v1[1] - eps_);
				}
				else return (_v0[0] < _v1[0] - eps_);
			}

		private:
			float eps_;
	};


	~ManageObj();//���������ر��ļ���
	int OpenObj();
	
	void ReadObj();
	void ReadXyz();
	void ReadStl();
	void ReadMeshPoints(vector<MyMesh::Point>& ab);
	void ReadMeshPoints2(vector<MyMesh::Point>& ab);

	static void OutFileOutlinePointXyz(vector<Vector3f>* vp, char * outfilename);
	static void OutFileOutlinePointXyz(vector<MyMesh::Point>* vp, char * outfilename);

	static void OutFilePointObj(vector<Vector3f>* vp,const char * outfilenam);
	static void OutFilePointObj(vector<float> *a, const char * outfilenam);
	static void OutFileVectorFloat(vector<float>&a, char * outfilename);
	static void OutFilePointObj(vector<Vector4f>*vp, const char * outfilenam);
	static void OutFilePointObj(vector<MyMesh::Point>&mp, const char * outfilename);

	static void OutFilePointAna(set<MyOutBottom>&mp, const char * outfilename);
	

	//void ObjLine();
private:
	char *mrfilename;
	char *mbuffer;
	vector<Vector3f> *mvpoint;
	vector<Vector3i> *mptri; //һ��������Ƭ������������

	vector<VertexTriMesh> *mMeshTri; //���¶�����mesh��������

	Vector3f StringFloat2(string s);
	Vector3f StringToFloatV(string s);
	Vector3i StringToIntV(string s);

	MyMesh::Point StringFloat3(string s);
	MyMesh::Point StringFloat4(string s);

	void FreeBuf() { free(mbuffer); };

	map<Vector3f, uint, CmpVec>  mVMap;
	//map<Vector3f, uint, CmpVec>::iterator mVMapIt;

	int cpyint(const char*& p)
	{
		int cpy;
		memcpy((char*)&cpy, p, 4);
		p += 4;
		return cpy;
	}
	float cpyfloat(const char*& p)
	{
		float cpy;
		/*char *writer;
		writer = (char*)&cpy;*/
		memcpy((char*)&cpy, p, 4);
		p += 4;
		return cpy;
	}
};






