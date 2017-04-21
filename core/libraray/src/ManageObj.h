#pragma once
#include "PointVector.h"
#include "MyOpenMesh.h"

struct VertexTri {
	uint i;  //编号
	uint p[3];  //三个顶点的索引编号
	uint v[3];  //该三角形相邻的三个三角形编号
	Vector3f a; //该面的法向量
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


	~ManageObj();//可以用来关闭文件夹
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
	

	//void ObjLine();
private:
	char *mrfilename;
	char *mbuffer;
	vector<Vector3f> *mvpoint;
	vector<Vector3i> *mptri; //一个三角面片的三个顶点编号

	vector<VertexTriMesh> *mMeshTri; //重新定义了mesh的三角形

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






