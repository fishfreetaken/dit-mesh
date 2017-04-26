#include "stdafx.h"
#include "MyOpenMesh.h"

MyOpenMesh::MyOpenMesh(float a)
{
}

void MyOpenMesh::ReadStlfile(char * argg) {
	// read mesh from stdin
	if (!OpenMesh::IO::read_mesh(mesh, argg, opt))
	{
		std::cerr << "Error: Cannot read mesh from " << std::endl;
		return;
	}
	//vertex_normals init!
	mesh.request_vertex_normals();

	// assure we have vertex normals
	if (!mesh.has_vertex_normals())
	{
		std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
		return ;
	}
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// we need face normals to update the vertex normals
		mesh.request_face_normals();
		// let the mesh update the normals
		mesh.update_normals();
		// dispose the face normals, as we don't need them anymore
		mesh.release_face_normals();
	}
}

void MyOpenMesh::WriteStlfile(char *argg,int i){
	if (!OpenMesh::IO::write_mesh(mesh, argg ,i)) //0 ascii 1 binary
	{
		std::cerr << "Error: cannot write mesh to " << argg << std::endl;
	}
}

void MyOpenMesh::FindNearest(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c, MyMesh::VertexHandle *p) {  //2
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	v_it_s[3];
	float min[3] = { 9999,9999,9999 };
	float lin;
	MyMesh::Point point[3], pp;
	cout << "Finding Nearest Vertex..." << endl;
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		pp = mesh.point(*v_it);
		lin = (pp - a).norm();
		if (lin < min[0]) {
			min[0] = lin;
			v_it_s[0] = v_it;
		}
		lin = (pp - b).norm();
		if (lin < min[1]) {
			min[1] = lin;
			v_it_s[1] = v_it;
		}
		lin = (pp - c).norm();
		if (lin < min[2]) {
			min[2] = lin;
			v_it_s[2] = v_it;
		}
	}
	*(p) = mesh.vertex_handle(v_it_s[0]->idx());
	*(p+1) = mesh.vertex_handle(v_it_s[1]->idx());
	*(p+2) = mesh.vertex_handle(v_it_s[2]->idx());
}

Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a){
	return Vector3f(a[0], a[1], a[2]);
}
MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a) {
	return MyMesh::Point(a[0], a[1], a[2]);
}

void MyOpenMesh::MoveVertex(float len) {
	// move all vertices one unit length along it's normal direction
	cout << "Moving poits..." << endl;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin();
		v_it != mesh.vertices_end(); ++v_it)
	{
		mesh.set_point(*v_it, mesh.point(*v_it) + mesh.normal(*v_it)*len); //移动指定距离；也是设置点位置的一种方式
	 }
}

void MyOpenMesh::ReleaseVertexNormals() {
	mesh.release_vertex_normals();
	if (mesh.has_vertex_normals())
	{
		std::cerr << "Ouch! ERROR! Shouldn't have any vertex normals anymore!\n";
	}
}

void MyOpenMesh::BottomVertex(vector<Vector3f>*a) {
	MyMesh::Normal m;
	MyMesh::Point cc;
	cout << "extracting bottom points...." << endl;
	/*if (!mHeelHight) {
		cout << "heel hight no init!" << endl;
	}*/
	for (MyMesh::VertexIter v_it = mesh.vertices_begin();
		v_it != mesh.vertices_end(); ++v_it)
	{
		//std::cout << "Vertex Bottom #" << *v_it << endl; //加调试输出太占用时间
		m = mesh.normal(*v_it); //>0是鞋底，0<是鞋楦表面
		if (m[2] > 0) { 
			cc = mesh.point(*v_it);
			//if (cc[2] <= (mHeelHight + 6)) {//使用跟高来将鞋底进行分离
			//	a->push_back(Vector3f(cc[0], cc[1], cc[2]));
			//}
		}
	}
}

set<MyOutBottom> MyOpenMesh::ShoeBottomLine(vector<MyMesh::Point>&a) {
	struct vmmt {
		MyMesh::VertexHandle n;
		bool operator < (const struct vmmt &a)const {
			return n.idx() < a.n.idx();
		}
	};
	struct vmt {
		MyMesh::Point a;
		float x;
	}git;
	cout << "Out the Bottom Line" << endl;
	set<MyOutBottom>arr;
	vector<struct vmt> xxt;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		MyMesh::Normal es(0, 0, 0),dif(0,0,0);

		//set<struct vmmt> xt;
		//struct vmmt xxx;
		//xxx.n = mesh.vertex_handle(v_it->idx());
		////xxx.i = v_it->idx();
		//xt.insert(xxx);
		//for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		//{
		//	for (MyMesh::VertexVertexIter vv_it1 = mesh.vv_iter(mesh.vertex_handle(vv_it->idx())); vv_it1.is_valid(); ++vv_it1)
		//	{
		//		for (MyMesh::VertexVertexIter vv_it2 = mesh.vv_iter(mesh.vertex_handle(vv_it1->idx())); vv_it2.is_valid(); ++vv_it2)
		//		{
		//			xxx.n = mesh.vertex_handle(vv_it2->idx());
		//			xt.insert(xxx);
		//		}
		//		xxx.n = mesh.vertex_handle(vv_it1->idx());
		//		xt.insert(xxx);
		//	}
		//	xxx.n = mesh.vertex_handle(vv_it->idx());
		//	xt.insert(xxx);
		//}
		git.a = mesh.point(*v_it);
		git.x = 1;
		xxt.push_back(git);

		for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			for (MyMesh::VertexVertexIter vv_it1 = mesh.vv_iter(mesh.vertex_handle(vv_it->idx())); vv_it1.is_valid(); ++vv_it1)
			{
				for (MyMesh::VertexVertexIter vv_it2 = mesh.vv_iter(mesh.vertex_handle(vv_it1->idx())); vv_it2.is_valid(); ++vv_it2)
				{
					for (MyMesh::VertexVertexIter vv_it3 = mesh.vv_iter(mesh.vertex_handle(vv_it2->idx())); vv_it3.is_valid(); ++vv_it3)
					{
						git.a = mesh.point(*vv_it2);
						git.x = 0.1;
						xxt.push_back(git);
					}
					git.a = mesh.point(*vv_it2);
					git.x = 0.2;
					xxt.push_back(git);
				}
				git.a = mesh.point(*vv_it1);
				git.x = 0.6;
				xxt.push_back(git);
			}
			git.a = mesh.point(*vv_it);
			git.x = 0.8;
			xxt.push_back(git);
		}
		
		if (!xxt.size()) {
			continue;
		}
		float xii = 0;
		for (auto i:xxt) {
			es += i.a;
		}
		es = es / xxt.size();
		for (auto i:xxt)
		{
			dif[0]+= ((i.a[0] - es[0])*i.x) * ((i.a[0] - es[0])*i.x);
			dif[1]+= ((i.a[1] - es[1])*i.x) * ((i.a[1] - es[1])*i.x);
			dif[2]+= ((i.a[2] - es[2])*i.x) * ((i.a[2] - es[2])*i.x);
		}
		dif[0] = sqrt(dif[0]);
		dif[1] = sqrt(dif[1]);
		dif[2] = sqrt(dif[2]);
		xxt.clear();
		//if (!xt.size()) {
		//	continue;
		//}
		////cout << xt.size() << endl;
		//for (auto i : xt) {
		//	es += mesh.normal(i.n);
		//}
		//es = es / xt.size();
		//for (auto i:xt)
		//{
		//	dif[0]+= (mesh.normal(i.n) - es)[0] * (mesh.normal(i.n) - es)[0];
		//	dif[1]+= (mesh.normal(i.n) - es)[1] * (mesh.normal(i.n) - es)[1];
		//	dif[2]+= (mesh.normal(i.n) - es)[2] * (mesh.normal(i.n) - es)[2];
		//}
		//dif[0] = sqrt(dif[0]);
		//dif[1] = sqrt(dif[1]);
		//dif[2] = sqrt(dif[2]);

		MyOutBottom stt;
		stt.n = mesh.vertex_handle(v_it->idx());
		stt.x = dif.norm();
		stt.s = mesh.normal(*v_it);
		stt.a = mesh.point(*v_it);
		if ((stt.x >=4.79)&&(stt.x<4.8)) {
			a.push_back(mesh.point(*v_it));
		}	
		arr.insert(stt);
	}
	return arr;
}

vector<MyMesh::Point> MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arr, SurfaceCoe* met,vector<MyMesh::Point>&css) {
	MyMesh::Point  p,p1,p2;// n为递增向量
	float s1, s2; float lin=0;
	vector<MyMesh::Point>  result;
	for (int j=0;j<css.size();j++)// p :css)
	{
		p = css[j];
		if (met->DistSurface(p)<0) {
			for (int i = 0; i < CUTSECTION - 1; i++) {
				lin = arr[i]->DistSurface(p);
				if ( lin< 0) {
					if (i == (CUTSECTION - 2)) {
						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					continue;
				}else {
					if (lin > 0) {
						p2 = arr[i]->FindNearestPoint(p, s2);
						p1 = i ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1);
						p+= p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					}else {
						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					break;
				}
			}
		}else if (met->DistSurface(p) > 0) {
			for (int i = (CUTSECTION-1); i < arr.size(); i++) {
				lin = arr[i]->DistSurface(p);
				if (lin> 0) {
					if (i == (arr.size() - 1)) {
						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					continue;
				}
				else {
					if (lin < 0) {
						p2 = arr[i]->FindNearestPoint(p, s2);
						p1 = i>(CUTSECTION-1) ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1);
						p+= p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					}else {
						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					break;
				}
			}
		}else {
			p+=met->FindNearestPoint(p,s1); //s==0;
		}
		result.push_back(p);
	}
	return result;
}

void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arr, SurfaceCoe* met) {
	//MyMesh::Normal m;
	cout << "Now is shoe Expansing..." << endl;
	MyMesh::Point p, p1, p2;// n为递增向量
	float s1, s2; float lin = 0;
	int k = mesh.n_vertices();
	int oc = 0;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		oc++;
		if (!(oc % 10000)) {
			cout << "now is :" << (oc * 100 / k) << "%" << endl;
		}
		/*if (mesh.normal(*v_it)[2] > 0) {
			continue;
		}*/
		p = mesh.point(*v_it);
		if (met->DistSurface(p)<0) {
			for (int i = 0; i < CUTSECTION - 1; i++) {
				lin = arr[i]->DistSurface(p);
				if (lin< 0) {
					if (i == (CUTSECTION - 2)) {
						p += arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					continue;
				}
				else {
					if (lin > 0) {
						p2 = arr[i]->FindNearestPoint(p, s2);
						p1 = i ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1);
						p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					}
					else {
						p += arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					break;
				}
			}
		}
		else if (met->DistSurface(p) > 0) {
			for (int i = (CUTSECTION - 1); i < arr.size(); i++) {
				lin = arr[i]->DistSurface(p);
				if (lin> 0) {
					if (i == (arr.size() - 1)) {
						p += arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					continue;
				}
				else {
					if (lin < 0) {
						p2 = arr[i]->FindNearestPoint(p, s2);
						p1 = i>(CUTSECTION - 1) ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1);
						p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					}
					else {
						p += arr[i]->FindNearestPoint(p, s1); //s==0;
					}
					break;
				}
			}
		}else {
			p += met->FindNearestPoint(p, s1); //s==0;
		}
		//if (p.norm() > 999) {  //测试是否有坏点
		//	cout<<"error points"<<endl;
		//}
		mesh.set_point(*v_it, p);
	}
}


void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arrx, SurfaceCoe* met,int ith) {
	//MyMesh::Normal m;
	if (ith == -1) {
		cout << "Shoe Expanding init error!" << endl;
	}
	vector<SurfaceCoe*>arr;
	for (int i = 0; i <= arrx.size(); i++) {
		if (i == ith) {
			arr.push_back(met);
			continue;
		}
		arr.push_back(arrx[i]);
	}
	cout << "Now is shoe Expansing..." << endl;
	MyMesh::Point p, p1, p2;// n为递增向量
	float s1, s2; float lin = 0;
	int k = mesh.n_vertices();
	int oc = 0;
	int arrlen = arr.size();
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		oc++;
		if (!(oc % 10000)) {
			cout << "now is :" << (oc * 100 / k) << "%" << endl;
		}
		p = mesh.point(*v_it);
		//int len = arr.size();
		int i = arrlen /2;
		while (1) {
			if (arr[i]->DistSurface(p) > 0) {
				if (arr[i-1]->DistSurface(p) > 0) {
					i -= i / 2;
					if (i == 1) {
						p += arr[i-1]->FindNearestPoint(p, s1); //s==0;
						break;
					}
				}else if(arr[i-1]->DistSurface(p)<0){
					p2 = arr[i]->FindNearestPoint(p, s2);
					p1 = arr[i - 1]->FindNearestPoint(p, s1);
					p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					break;
				}else {
					p += arr[i]->FindNearestPoint(p, s1); //s==0;
					break;
				}
			}
			else if (arr[i]->DistSurface(p) < 0) {
				if (arr[i + 1]->DistSurface(p) < 0) {
					i += (arrlen - i) / 2;
					if (i == (arrlen - 1)) {
						p += arr[i - 1]->FindNearestPoint(p, s1); //s==0;
						break;
					}
				}else if (arr[i - 1]->DistSurface(p)>0) {
					p2 = arr[i]->FindNearestPoint(p, s2);
					p1 = arr[i - 1]->FindNearestPoint(p, s1);
					p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					break;
				}else {
					p += arr[i]->FindNearestPoint(p, s1); //s==0;
					break;
				}
			}else {
				p += arr[i]->FindNearestPoint(p, s1); //s==0;
				break;
			}
		}
		
		if (p.norm() > 999) {  //测试是否有坏点
			cout<<"error points"<<endl;
		}
		mesh.set_point(*v_it, p);
	}
}

SurfaceCoe::SurfaceCoe(MyMesh::VertexHandle *vertex, MyMesh &b):
	mesh(b)
{
	mVertexStart = mesh.point(*vertex);
	mVertexMid = mesh.point(*(vertex + 1));
	mVertexEnd = mesh.point(*(vertex + 2));
	mHandleBegin = *vertex;
}
SurfaceCoe::SurfaceCoe(MyMesh::VertexHandle a, MyMesh::Point c, MyMesh &d):
	mesh(d),
	mHandleBegin(a)
{
	mVertexStart = mesh.point(a);
	mVertexEnd = c;
}
SurfaceCoe::SurfaceCoe(Vector3f bf, MyMesh::VertexHandle df, float x, MyMesh &b):
	mesh(b),
	mCoeABC(bf),
	mHandleBegin(df),
	mX(x)
{
	mVertexStart = mesh.point(df);
	float d = mCoeABC.dot(Vector3f((0 - mVertexStart[0]), (0 - mVertexStart[1]), (0 - mVertexStart[2])));
	mCoe= Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);
}

bool SurfaceCoe::Init() {//(MyMesh::VertexHandle *vertex,MyMesh &mmesh){
	MyMesh::VertexHandle vertex_i;
	MyMesh::HalfedgeHandle heh;
	MyMesh::Point st;
	for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mHandleBegin); voh_it.is_valid(); ++voh_it) {
		heh = mesh.halfedge_handle(voh_it->idx());
		Vector3f jug = NextHalfEdgeJudge2(mesh.next_halfedge_handle(heh));
		if (jug != Vector3f(9999, 9999, 9999)) {
			if (mVertexStart[1]< jug[1]) {
				heh = mesh.next_halfedge_handle(heh);
				IterationHalfEdge(heh);
				return true;
			}
		}
	}
	return false;
}

bool SurfaceCoe::Init(int cmd){//(MyMesh::VertexHandle *vertex,MyMesh &mmesh){
	mOutline2.clear();
	mLength = 0;

	MyOutNormal abc;
	abc.a = mVertexStart;
	abc.n = mesh.normal(mHandleBegin);
	mOutline2.push_back(abc);

	Vector3f a= MyOpenMesh::EigenTransfer(mVertexStart);
	Vector3f b= MyOpenMesh::EigenTransfer(mVertexMid);
	Vector3f c= MyOpenMesh::EigenTransfer(mVertexEnd);
	Vector3f ab = a - b;
	ab = (a-c).cross(ab);
	mCoeABC = Vector3f(ab[0], ab[1], ab[2]);
	
	float d = ab.dot(Vector3f(0, 0, 0) - a)/ mCoeABC.norm();
	mCoeABC.normalize();
	mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);
	//float d = ab.dot(Vector3f(0, 0, 0) - a);
	//mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);

	MyMesh::VertexHandle vertex_i;
	MyMesh::HalfedgeHandle heh;
	MyMesh::Point st;
	bool ini = 0;
	for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mHandleBegin); voh_it.is_valid(); ++voh_it) {
		heh = mesh.halfedge_handle(voh_it->idx());
		//vertex_i = mesh.to_vertex_handle(heh);
		//st = mesh.point(vertex_i);
		//if (st[2] > mVertexStart[2]) { //沿着z轴方向往上走（正方向）
		//	if (NextHalfEdgeJudge(mesh.next_halfedge_handle(heh))) {
		//		ini = 1;
		//		heh = mesh.next_halfedge_handle(heh);
		//		IterationHalfEdge(heh);
		//	}
		//}
		Vector3f jug = NextHalfEdgeJudge2(mesh.next_halfedge_handle(heh));
		if (jug != Vector3f(9999, 9999, 9999)) {
			if (mVertexStart[2]< jug[2]) {
				heh = mesh.next_halfedge_handle(heh);
				IterationHalfEdge(heh);
				ini = 1;
			}
		}
	}
	if (ini) {
		if (cmd) {
			OutlineRefine();
		}
		CoquerMidEnd();
	}
	return ini;
}

void SurfaceCoe::CoquerMidEnd() {
	float ins, mid = MAXIMUMX, end = MAXIMUMX;
	MyMesh::Point k = mVertexStart;
	for (int i = 1; i < mOutline2.size(); i++) {
		mLength += (mOutline2[i].a - k).norm();
		mOutline2[i].d = mLength;
		ins = DistPoints(mOutline2[i].a, mVertexMid);
		k = mOutline2[i].a;
		if (ins < mid) {
			mid = ins;
			mIth[0] = i;
			mLen[0] = mLength;
		}
		ins = DistPoints(mOutline2[i].a, mVertexEnd);
		if (ins < end) {
			end = ins;
			mIth[1] = i;
			mLen[1] = mLength - mLen[0];
		}
	}
	mLen[2] = mLength - mLen[1] - mLen[0];
}

void SurfaceCoe::OutlineRefine() {
 	MyMesh::Point k;
	k = mOutline2[0].a;
	int j = 0;
	for (int i = mOutline2.size() / 8; i < mOutline2.size(); i++) {
		if (DistPoints(mOutline2[i].a, k) < 0.003) {
			j = i;
			cout << i << endl;
			break;
		}
	}
	vector<MyOutNormal> st;
	for (int i = 0; i < j; i++) {
		st.push_back(mOutline2[i]);
	}
	mOutline2 = st;
}

void SurfaceCoe::AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b) {
	Vector3f cac = MyOpenMesh::EigenTransfer(mesh.point(a));// (ca[0], ca[1], ca[2]);
	Vector3f cbc = MyOpenMesh::EigenTransfer(mesh.point(b));// (cb[0], cb[1], cb[2]);
	Vector3f v = cac - cbc;
	float t;
	t = (0 - (mCoeABC.dot(cac) + mCoe[3])) / (mCoeABC.dot(v));
	MyOpenMesh::OutNoraml mc;
	cbc = Vector3f(v[0] * t + cac[0], v[1] * t + cac[1], v[2] * t + cac[2]);

	mc.a = MyMesh::Point(cbc[0], cbc[1], cbc[2]);
	mc.n = (mesh.normal(a) + mesh.normal(b)) / 2;//平面法向量投影
	Vector3f nv(mc.n[0], mc.n[1], mc.n[2]);
	Vector3f ln = nv - mCoeABC*(nv.dot(mCoeABC));
	ln.normalize();
	mc.n = MyMesh::Normal(ln[0], ln[1], ln[2]);
	mOutline2.push_back(mc);
}

//float SurfaceCoe::AllocateXCoe(float tar) {
//	if ((!mLen[0]) || (!mLen[1])) {
//		cout << "zero error!" << endl;
//		return -1;
//	}
//	for (int i = 0; i <= mIth[0]; i++) {
//		mOutline2[i].x = mOutline2[i].d / mLen[0];
//	}
//	for (int i = mIth[1] - 1; i >= mIth[0]; i--) {
//		mOutline2[i].x = abs(mOutline2[i].d - mLen[0]) / mLen[1];  //1
//	}
//
//	mExtension = tar > 0 ? -tar : tar;
//	return OutlineExpansion();
//}

float SurfaceCoe::AllocateXCoe(float tar) {
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error!" << endl;
		return -1;
	}
	float mm = mLen[0] > mLen[1] ? mLen[0] : mLen[1];

	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = mOutline2[i].d / mm; // mLen[0];
	}
	for (int i = mIth[1]; i >= mIth[0]; i--) {
		mOutline2[i].x = 1 - abs(mOutline2[i].d - mLen[0]) / mm;  // mLen[1];  //1
	}
	vector<float>input, output;
	for (int i = 0; i <=mIth[1]; i++) {
		input.push_back(mOutline2[i].x);
	}
	GaussianSmooth(input, output, 2);
	for (int i = 0; i <= mIth[1]; i++) {
		mOutline2[i].x = output[i];
	}

	mExtension = tar > 0 ? -tar : tar;
	return OutlineExpansion();
}

//float SurfaceCoe::AllocateXCoe2() {
//	if ((!mLen[0]) || (!mLen[1])) {
//		cout << "zero error2!" << endl;
//		return -1;
//	}
//	for (int i = 0; i <= mIth[0]; i++) {
//		mOutline2[i].x = (1 - mOutline2[i].d / mLen[0]);// (mLen[0] / mLen[2]);
//	}
//	for (int i = mOutline2.size()-1; i >= mIth[1]; i--) {
//		mOutline2[i].x = (mOutline2[i].d - mLen[1] - mLen[0]) / mLen[2];  //2
//	}
//
//	mExtension = mX >0 ? -mLength*mX : mLength*mX; //由扩散系数转换为扩散距离
//
//	return OutlineExpansion();
//}

float SurfaceCoe::AllocateXCoe2() {
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error2!" << endl;
		return -1;
	}
	float mm = mLen[0] > mLen[2] ? mLen[0] : mLen[2];
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = (1 - mOutline2[i].d / mm);//mLen[0]
	}
	for (int i = mOutline2.size() - 1; i >= mIth[1]; i--) {
		mOutline2[i].x = (mOutline2[i].d - mLen[1] - mLen[0]) / mm;  //mLen[2]
	}

	vector<float>input,output;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		input.push_back(mOutline2[i].x);
	}
	for (int i = 0; i <= mIth[0]; i++) {
		input.push_back(mOutline2[i].x);
	}
	GaussianSmooth(input,output,2);

	for (int i = 0; i < mOutline2.size() - mIth[1]; i++) {
		mOutline2[i + mIth[1]].x = output[i];
	}
	for (int i = 0; i <=mIth[0] ; i++) {
		mOutline2[i].x = output[i + mOutline2.size() - mIth[1]];
	}

	mExtension = mX >0 ? -mLength*mX : mLength*mX; //由扩散系数转换为扩散距离

	return OutlineExpansion();
}

//float SurfaceCoe::OutlineExpansion() {  //(float tar)//一般如果分布之间差异过大，导致收敛摆动过大，收敛效果并不好
//	vector<MyOutNormal> bso = mOutline2;
//	float s = 0, li = mExtension*10, pp = 0;
//	while (abs(abs(s - mLength) - abs(mExtension)) > ADDDIFFERENCE) {
//		bso = mOutline2;
//		for (int j = 0; j < bso.size(); j++) {
//			bso[j].a += bso[j].n*bso[j].x*li;
//		}
//		s = TotalLengh(bso);
//		pp = abs(mExtension)/(s - mLength);
//		li *= pp > 0 ? pp : pp<-1? -pp: -1/pp;
//		//li *= pp > 0 ? pp : pp < -1 ? -pp : 1 - pp;
//	}
//	for (int j = 0; j < mOutline2.size(); j++) {
//		mOutline2[j].f = mOutline2[j].n*mOutline2[j].x*li;
//		mOutline2[j].a += mOutline2[j].f;
//	}
//	return li;
//}

float SurfaceCoe::OutlineExpansion() {  //(float tar) //这个是取半值
	vector<MyOutNormal> bso = mOutline2;
	float sout = mExtension * 20;
	float s = 0, li = sout/2, pp = 0;
	float aa = 0, bb = 1, cc = 0.5;
	while (abs(abs(s - mLength) - abs(mExtension)) > ADDDIFFERENCE) {
		bso = mOutline2;
		for (int j = 0; j < bso.size(); j++) {
			bso[j].a += bso[j].n*bso[j].x*li;
		}
		s = TotalLengh(bso);
		pp = abs(mExtension) / (s - mLength);
		if ((pp > 1)||(pp<0)) {
			aa = cc;
			cc = cc+(bb - cc) / 2;
			li = sout*cc;
		}
		else if (pp < 1) {
			bb = cc;
			cc =aa+(cc - aa) / 2;
			li = sout*cc;
		}
		else {
			break;
		}
	}
	for (int j = 0; j < mOutline2.size(); j++) {
		mOutline2[j].f = mOutline2[j].n*mOutline2[j].x*li;
		mOutline2[j].a += mOutline2[j].f;
	}
	return li;
}

bool SurfaceCoe::OutlineEigen(vector<Vector3f> *a) {
	if (!mOutline2.size()) {
		cout << "outline no points" << endl;
		return false;
	}
	for (auto i : mOutline2) {
		a->push_back(Vector3f(i.a[0], i.a[1], i.a[2]));
	}
	return true;
}

bool SurfaceCoe::OutlineEigen(vector<Vector4f> *a) {
	if (!mOutline2.size()) {
		cout << "outline no points" << endl;
		return false;
	}
	for (auto i : mOutline2) {
		a->push_back(Vector4f(i.a[0], i.a[1], i.a[2],i.x));
	}
	return true;
}

SurfaceCoe* SurfaceCoe::FindMetara(MyMesh::VertexHandle start, MyMesh::VertexHandle end) {
	SurfaceCoe sfc(start, mesh.point(end), mesh);
	cout << "Now is finding Metara..." << endl;
	float dd = MAXIMUMX ,ss;
	for (int i = mIth[0]*0.13; i < mIth[0]*0.5; i++) {
		sfc.SetMidPoint(mOutline2[i].a);
		if (sfc.Init(0)) {
			ss = sfc.ReturnLength();
			if (ss < dd) {
				dd = ss;
				mIth[2] = i;
			}
		}
	}
	cout << "j:" << mIth[2] << endl;
	SurfaceCoe *ret= new SurfaceCoe(start, mesh.point(end), mesh);
	ret->SetMidPoint(mOutline2[mIth[2]].a);
	ret->Init(0);
	return ret;
}

Vector3f SurfaceCoe::AxieCut(float heelhight) {
	if (heelhight < 5) {
		cout << "heelhight is too loaw!" << endl;
		return Vector3f(0,0,0);
	}
	float aim = (heelhight - 5)*0.4 + 30;
	if (abs(mOutline2[mIth[1]].d - mOutline2[mIth[0]].d)>aim) {
		cout << "heelhight is too long error!" << endl;
		return Vector3f(0, 0, 0);
	}
	float ss = 0;
	Vector3f axie;
	MyMesh::Point k = mOutline2[mIth[1]].a;
	for (int i = mIth[1]-1; i > mIth[0]; i--) {
		ss += (k - mOutline2[i].a).norm();
		k = mOutline2[i].a;
		if (abs(ss - aim) < 0.8) {
			axie=Vector3f(mOutline2[i].a[0]-mVertexStart[0], mOutline2[i].a[1] - mVertexStart[1], mOutline2[i].a[2] - mVertexStart[2]);
			axie.normalize();
			break;
		}
	}
	return axie;
}

void SurfaceCoe::OutCutOutline(float x ,vector<MySurCutArry> &arryx) {
	struct CutArry2 {
		MyMesh::Point a;
		float x = 1;
	};
	vector<struct CutArry2> arry;
	int interv=(mIth[0] - mIth[2])/ CUTSECTION; //把后一段分成22段，给出21个截面
	for (int i = 1; i < CUTSECTION; i++) {
		struct CutArry2 st;
		st.a = mOutline2[mIth[2]+i*interv].a;
		st.x = x;
		arry.push_back(st);
	}
	interv = mIth[2] / CUTSECTION2; //把后一段分成17段，给出16个截面
	for (int i = 1; i < CUTSECTION2; i++) {
		struct CutArry2 st;
		st.a = mOutline2[mIth[2] - i*interv].a;
		st.x = x*(CUTSECTION2-i)/ CUTSECTION2;
		arry.push_back(st);
	}

	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	*v_it_s = (MyMesh::VertexIter*)malloc(sizeof(MyMesh::VertexIter)*arry.size());
	float * minar =(float*)malloc(sizeof(float)*arry.size());
	float lin;
	for (int i = 0; i < arry.size(); i++) {
		*(minar + i) = MAXIMUMX;
	}
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		for (int i = 0; i < arry.size();i++) {
			lin = (arry[i].a - mesh.point(*v_it)).norm();
			if (*(minar + i) > lin) {
				*(v_it_s + i) = v_it;
				*(minar + i) = lin;
			}
		}
	}
	
	for (int i = 0; i < arry.size(); i++) {
		MySurCutArry stt;
		stt.a = mesh.vertex_handle((*(v_it_s + i))->idx());
		stt.x = arry[i].x;
		arryx.push_back(stt);
	}
	free(minar);
	free(v_it_s);
}

vector<MySurCutArry> SurfaceCoe::OutCutOutline(float exp,SurfaceCoe *met,Vector3f axi,int &ith) {
	float xi = exp / met->ReturnLength();
	struct CutArry2 {
		MyMesh::Point a;
		float x = 1;
		Vector3f n; 
	};
	float ju,max0,max1;
	max0 = met->DistSurface(mVertexStart);
	max1 =abs(met->DistSurface(mVertexMid));

	Vector3f coemet = met->ReturnCoe();
	float thert=acos(coemet.dot(axi));
	Vector3f mix = coemet.cross(axi);

	//Vector4f cmix(thert,mix[0],mix[1],mix[2]);
	float lin;
	bool ini = 0;
	vector<struct CutArry2> arry;
	Quaternionx out;
	int interv = mIth[0] / CUTSECTION3; //把第一段分成37份
	for (int i = 1; i < CUTSECTION3; i++) {
		struct CutArry2 st;
		st.a = mOutline2[i*interv].a;
		ju=met->DistSurface(st.a);
		if (ju >0) {
			st.x = xi*(1-ju/max0);
			lin = thert*(ju / max0 );
		}else if(ju<0){
			if (!ini) {
				ini = true;
				ith = i-1;
			}
			st.x = xi;
			lin = thert*(abs(ju) / max1);
		}else {
			continue;
		}
		Quaternionx mtransfer(cos(lin),sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0,coemet[0],coemet[1],coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(),out.y(),out.z());
		arry.push_back(st);
	}

	vector<MySurCutArry>arryx;
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	*v_it_s = (MyMesh::VertexIter*)malloc(sizeof(MyMesh::VertexIter)*arry.size());
	float * minar = (float*)malloc(sizeof(float)*arry.size());

	for (int i = 0; i < arry.size(); i++) {
		*(minar + i) = MAXIMUMX;
	}
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		for (int i = 0; i < arry.size(); i++) {
			lin = (arry[i].a - mesh.point(*v_it)).norm();
			if (*(minar + i) > lin) {
				*(v_it_s + i) = v_it;
				*(minar + i) = lin;
			}
		}
	}

	for (int i = 0; i < arry.size(); i++) {
		MySurCutArry stt;
		stt.a = mesh.vertex_handle((*(v_it_s + i))->idx());
		stt.x = arry[i].x;
		stt.n = arry[i].n;
		arryx.push_back(stt);
	}
	free(minar);
	free(v_it_s);
	return  arryx;
}

Vector3f SurfaceCoe::TempVector() {
	return Vector3f(mVertexStart[0]-38.5335, mVertexStart[1]+0.4859, mVertexStart[2]-157.5532).normalized();
	//return Vector3f(38.5335 - mVertexStart[0], -0.4859 - mVertexStart[1], 157.5532 - mVertexStart[2]).normalized();
}

MyMesh::Point SurfaceCoe::FindNearestPoint(MyMesh::Point a, float &s) {
	int k = 0, n = 0, m = 0;
	s = abs(DistSurface(a));
	float ss = 0, mig = MAXIMUMX;
	for (int i = 0; i < mOutline2.size(); i++) {
		ss = DistPoints(a, mOutline2[i].a);
		if (ss < mig) {
			mig = ss;
			k = i;
		}
	}
	if (!mOutline2[k].x) {
		return MyMesh::Point(0, 0, 0);
	}
	MyMesh::Point fb, fc, fd, j1, j2, p;
	fb = mOutline2[k].a;
	float t;
	for (int i = 1; i < 4; i++) {
		m = k - i;
		m = m < 0 ? mOutline2.size()+m : m;

		fc = mOutline2[m].a;
		fd = fb - fc;
		t = ((a[0] - fb[0])*fd[0] + (a[1] - fb[1])*fd[1] + (a[2] - fb[2])*fd[2]) / (fd[0] * fd[0] + fd[1] * fd[1] + fd[2] * fd[2]); //两直线的交点
		p = t*fd + fb;
		j1 = fb - p;
		j2 = fc - p;
		if ((j1[0] * j2[0] + j1[1] * j2[1] + j1[2] * j2[2]) < 0) {
			float s1 = j1.norm(), s2 = j2.norm();
			return (mOutline2[k].f*s2 / (s1 + s2) + mOutline2[m].f*s1 / (s1 + s2));
		}

		n = k + i;
		n = n >= mOutline2.size() ? n- mOutline2.size() : n;

		fc = mOutline2[n].a;
		fd = fb - fc;
		t = ((a[0] - fb[0])*fd[0] + (a[1] - fb[1])*fd[1] + (a[2] - fb[2])*fd[2]) / (fd[0] * fd[0] + fd[1] * fd[1] + fd[2] * fd[2]); //Vector3f s = t*d + a; 	t = (c-a).dot(d) / (d.dot(d));
		p = t*fd + fb;
		j1 = fb - p;
		j2 = fc - p;
		if ((j1[0] * j2[0] + j1[1] * j2[1] + j1[2] * j2[2]) < 0) {
			float s1 = j1.norm(), s2 = j2.norm();
			return (mOutline2[k].f*s2 / (s1 + s2) + mOutline2[n].f*s1 / (s1 + s2));
		}
	}
	return mOutline2[k].f;//MyMesh::Point(0,0,0);  // mOutline2[k].f
}
