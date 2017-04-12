#include "stdafx.h"
#include "MyOpenMesh.h"

MyOpenMesh::MyOpenMesh(float a):
	mHeelHight(a)
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

void MyOpenMesh::GoOutline2(Vector3f a, Vector3f b, Vector3f c) {
	Vector3f p[3];
	MyMesh::Point ca(a[0],a[1],a[2]);
	MyMesh::Point cb(b[0], b[1], b[2]);
	MyMesh::Point cc(c[0], c[1], c[2]);
	FindNearest(ca,cb,cc,p);
	SetSurfaceEquation(p[0], p[1],p[2]);

	MyMesh::VertexHandle vertex_i;
	MyMesh::HalfedgeHandle heh;
	MyMesh::Point st, at;
	at = mesh.point(mVertexStart);
	int init = 0;
	for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mVertexStart); voh_it.is_valid(); ++voh_it) {
		heh = mesh.halfedge_handle(voh_it->idx());
		vertex_i = mesh.to_vertex_handle(heh);
		st = mesh.point(vertex_i);
		if (st[2] > at[2]) { //z轴方向往上走
			if (NextHalfEdgeJudge(mesh.next_halfedge_handle(heh))) {
				init = 1;
				heh = mesh.next_halfedge_handle(heh);
				IterationHalfEdge(heh);
			}
		}
	}
	if (!init) {
		cout << "start point not init!" << endl;
	}
}

bool MyOpenMesh::GoOutline(Vector3f a, Vector3f b, Vector3f c,vector<struct OutNoraml>*lm) {
	//Vector3f p[3];//FindNearest(ca, cb, cc, p);在外部做
	SetSurfaceEquation(a, b, c);

	MyMesh::VertexHandle vertex_i;
	MyMesh::HalfedgeHandle heh;
	MyMesh::Point st, at;
	at = mesh.point(mVertexStart);
	bool init = 0;
	for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mVertexStart); voh_it.is_valid(); ++voh_it) {
		heh = mesh.halfedge_handle(voh_it->idx());
		vertex_i = mesh.to_vertex_handle(heh);
		st = mesh.point(vertex_i);
		if (st[2] > at[2]) { //z轴方向往上走
			if (NextHalfEdgeJudge(mesh.next_halfedge_handle(heh))) {
				init = 1;
				heh = mesh.next_halfedge_handle(heh);
				IterationHalfEdge(heh);
			}
		}
	}
	if (!init) {
		cout << "Start point not Init!" << endl;
	}
	return init;
}

void MyOpenMesh::FindNearest(MyMesh::Point a,MyMesh::Point b,MyMesh::Point c,Vector3f *p) {  //2
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	v_it_s[3];
	float min[3] = {9999,9999,9999};
	float lin;
	MyMesh::Point point[3],pp;

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		pp = mesh.point(*v_it);
		lin = (pp - a).norm();
		if (lin<min[0]) {
			min[0] = lin;
			v_it_s[0] = v_it;
		}
		lin = (pp - b).norm();
		if (lin<min[1]) {
			min[1] = lin;
			v_it_s[1] = v_it;
		}
		lin = (pp- c).norm();
		if (lin<min[2]) {
			min[2] = lin;
			v_it_s[2] = v_it;
		}
	}
	mVertexStart = mesh.vertex_handle(v_it_s[2]->idx());
	mVertexMid = mesh.vertex_handle(v_it_s[1]->idx());
	mVertexEnd = mesh.vertex_handle(v_it_s[0]->idx());
	point[0]=mesh.point(mVertexStart); point[1] = mesh.point(mVertexMid); point[2] = mesh.point(mVertexEnd);
	*p = Vector3f(point[0][0], point[0][1], point[0][2]);
	*(p+1)= Vector3f(point[1][0], point[1][1], point[1][2]);
	*(p+2)= Vector3f(point[2][0], point[2][1], point[2][2]);
	MyOutNormal abc;
	abc.a = point[0];
	abc.n = mesh.normal(mVertexStart);
	abc.x = 0;
	mOutline2.push_back(abc);
	return;//Vector3f(point[0], point[1], point[2]);
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

int MyOpenMesh::NextHalfEdgeJudge(MyMesh::HalfedgeHandle heh) {
	MyMesh::VertexHandle vertex_i = mesh.to_vertex_handle(heh);	//指向点
	MyMesh::VertexHandle vertex_s = mesh.from_vertex_handle(heh);  //出发点
	float ac = PointJudge(DistSurface(mesh.point(vertex_s)), DistSurface(mesh.point(vertex_i)));
	if (ac != -1) {
		return 1;
	}
	return 0;
}

void MyOpenMesh::IterationHalfEdge(MyMesh::HalfedgeHandle heh) {
	//偏移至下一条半边，以vertex为基点索引该边对应顶点
	MyMesh::HalfedgeHandle heh_next = heh;//=mesh.next_halfedge_handle(heh);
	MyMesh::VertexHandle vertex_i, vertex_s;
	while (vertex_i.idx() != mVertexStart.idx()) {
		vertex_i = mesh.to_vertex_handle(heh_next);	//指向点
		vertex_s = mesh.from_vertex_handle(heh_next);  //出发点
		float ac = PointJudge(DistSurface(mesh.point(vertex_s)), DistSurface(mesh.point(vertex_i)));
		if (ac == 0) {
			AddOutlinePoint(vertex_i, vertex_s);
			heh_next = mesh.next_halfedge_handle(heh_next);
			heh_next = mesh.opposite_halfedge_handle(heh_next);
			heh_next = mesh.next_halfedge_handle(heh_next);
			heh_next = mesh.next_halfedge_handle(heh_next);
		}
		else if (ac == 1) {
			AddOutlinePoint(vertex_i,vertex_s);
			heh_next = mesh.opposite_halfedge_handle(heh_next);
			heh_next = mesh.next_halfedge_handle(heh_next);
		}
		else {  //-1
			heh_next = mesh.next_halfedge_handle(heh_next);
		}
	}
}

int MyOpenMesh::PointJudge(float a, float b) {
	int i[2];
	if ((a == 0) || (b == 0)) {
		return 0;
	}
	i[0] = a < 0 ? -1 : 1;
	i[1] = b < 0 ? -1 : 1;

	return (i[0] + i[1]) != 0 ? -1 : 1;  //两个值符号相同 -1 //符号相反 1
}

void MyOpenMesh::AddOutlinePoint(MyMesh::VertexHandle a, MyMesh::VertexHandle b) {
	MyMesh::Point ca, cb;

	ca = mesh.point(a);
	cb = mesh.point(b);
	Vector3f cac(ca[0],ca[1],ca[2]);
	Vector3f cbc(cb[0],cb[1],cb[2]);
	Vector3f v = cac - cbc;
	float t;
	t = (0 - (mCoeABC.dot(cac) + mCoe[3])) / (mCoeABC.dot(v));
	MyOpenMesh::OutNoraml mc;
	cbc = Vector3f(v[0] * t + ca[0], v[1] * t + ca[1], v[2] * t + ca[2]);
	mc.a = MyMesh::Point(cbc[0],cbc[1],cbc[2]);
	mc.n = (mesh.normal(a) + mesh.normal(b)) / 2;//平面法向量投影
	Vector3f nv(mc.n[0], mc.n[1], mc.n[2]);
	Vector3f nf ( mCoeABC[0],mCoeABC[1],mCoeABC[2]);
	nf.normalize();
	Vector3f ln =nv- nf*(nv.dot(nf));
	ln.normalize();
	mc.n = MyMesh::Normal(ln[0],ln[1],ln[2]);
	mOutline2.push_back(mc);
}
Vector3f MyOpenMesh::EigenTransfer(MyMesh::Point a){
	return Vector3f(a[0], a[1], a[2]);
}
MyMesh::Point MyOpenMesh::MeshTransfer(Vector3f a) {
	return MyMesh::Point(a[0], a[1], a[2]);
}

void MyOpenMesh::SetSurfaceEquation(Vector3f a, Vector3f b, Vector3f c) {
	cout << "Openmesh Set Coe" << endl;
	Vector3f ab = b - a;
	Vector3f ac = c - a;
	ab = ac.cross(ab);
	float d;
	Vector3f zero(0, 0, 0);
	d = ab.dot(zero - a);
	mCoe = Vector4f(ab[0], ab[1], ab[2], d);
	mCoeABC = Vector3f(ab[0], ab[1], ab[2]);
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
	if (!mHeelHight) {
		cout << "heel hight no init!" << endl;
	}
	for (MyMesh::VertexIter v_it = mesh.vertices_begin();
		v_it != mesh.vertices_end(); ++v_it)
	{
		//std::cout << "Vertex Bottom #" << *v_it << endl; //加调试输出太占用时间
		m = mesh.normal(*v_it); //>0是鞋底，0<是鞋楦表面
		if (m[2] > 0) { 
			cc = mesh.point(*v_it);
			if (cc[2] <= (mHeelHight + 6)) {//使用跟高来将鞋底进行分离
				a->push_back(Vector3f(cc[0], cc[1], cc[2]));
			}
		}
	}
}

bool MyOpenMesh::OutlineEigen(vector<Vector3f> *a) {
	if (!mOutline2.size()) {
		cout << "outline no points" << endl;
		return false;
	}
	for (auto i : mOutline2) {
		a->push_back(Vector3f(i.a[0],i.a[1],i.a[2]));
	}
	return true;
}

float MyOpenMesh::DistSurface(MyMesh::Point a) {
	Vector4f sa(a[0], a[1], a[2], 1);
	return (sa.dot(mCoe) / mCoeABC.norm());
}

float MyOpenMesh::MetaraOutlineSort(float tar) {
	MyMesh::Point a, b, c;
	a = mesh.point(mVertexStart);
	b = mesh.point(mVertexMid);
	c = mesh.point(mVertexEnd);
	float sa = 0, sb = 0, af = 0, bf = 0;
	MyMesh::Point k = mOutline2[0].a;
	int jk = 0;
	int ja, jb;
	bool inverse = false;
	for (int i = 1; i < mOutline2.size(); i++) {
		sb += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		if (mOutline2[i].a == b) {
			jk++;
			ja = i;
			sa = sb;
			cout << "outline find mid!: "<<i << endl;
		}
		if (mOutline2[i].a == c) {
			jb = i;
			inverse = jk ? true : false;
			sb = sb-sa;
			cout << "outline find end!" << i << endl;
			break;
		}
	}
	k= mOutline2[0].a;
	for (int i = 1; i <= ja; i++) {
		af += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		mOutline2[i].x = af / sa;
	}
	k = c;
	for (int i = jb - 1; i >= ja; i--) {
		bf += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		mOutline2[i].x = bf / sb;
	}

	vector<MyOutNormal> bso = mOutline2;
	float s = 0, source = 0, li = 0;
	source = TotalLengh(bso);
	li = tar>0 ? -tar : tar;

	while (abs(abs(s - source) - abs(tar)) > ADDDIFFERENCE) {
		bso = mOutline2;
		for (int j = 0; j < b.size(); j++) {
			bso[j].a += bso[j].n*bso[j].x*li;
		}
		s = TotalLengh(bso);
		li *= abs(tar) / abs(s - source);
	}
	for (int j = 0; j < mOutline2.size(); j++) {
		mOutline2[j].a += mOutline2[j].n*mOutline2[j].x*li;
	}
	return tar>0 ? -li : li;
}

float MyOpenMesh::TotalLengh(vector<MyOutNormal>a) {
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
}

//38.8228 -1.6320  158.2785
void MyOpenMesh::CrossVector(vector<MyOutNormal> *v) {
	MyMesh::Point ee = MyMesh::Point(49.2967, -2.4764, 184.5850);//mesh.point(mVertexEE);  //move end! 49.2967,-2.4764,184.5850
	MyMesh::Point pc = MyMesh::Point(38.8228,-1.6320,158.2785)-mesh.point(mVertexStart); //中轴线向量
	Vector3f pcc(pc[0], pc[1], pc[2]);

	MyMesh::Point k = (*v)[0].a;
	MyMesh::Point bb = mesh.point(mVertexMid);
	MyMesh::Point sk;
	float ss = 0;
	for (int i = 0; i < v->size();i++) {
		sk = (*v)[i].a;
		ss += (sk - k).norm();
	}

	for (int i = 1; i < v->size(); i++) {
		sk = (*v)[i].a;
		ss += ( sk- k).norm();
		k = sk;
		if (sk != ee) {
			if (abs(ss - OUTLINEINTERVAL) < 0.07) {
				cout << "push one coe!" << endl;
				mSurfaceArray.push_back(Vector4f(pcc[0],pcc[1],pcc[2],(0-1)*pcc.dot(Vector3f(sk[0], sk[1],sk[2]))));
				ss = 0;
			}
		}
		else {
			cout << "make end!" << endl;
			break;
		}
	}
}

void MyOpenMesh::FindMetaraPoint(vector<MyOutNormal> *v) {
	MyMesh::Point mv;
	MyMesh::Point start = mesh.point(mVertexStart);
	MyMesh::Point end = mesh.point(mVertexMid);

	for (int i = 0; i < v->size(); i++) {
		mv = (*v)[i].a;
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
	MyMesh::Point st, at;
	bool ini = 0;
	for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mHandleBegin); voh_it.is_valid(); ++voh_it) {
		heh = mesh.halfedge_handle(voh_it->idx());
		vertex_i = mesh.to_vertex_handle(heh);
		st = mesh.point(vertex_i);
		if (st[2] > mVertexStart[2]) { // 沿着z轴方向往上走
			if (NextHalfEdgeJudge(mesh.next_halfedge_handle(heh))) {
				ini = 1;
				heh = mesh.next_halfedge_handle(heh);
				IterationHalfEdge(heh);
			}
		}
	}
	return ini;
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
	Vector3f ab = b - a;
	ab = (c-a).cross(ab);
	mCoeABC = Vector3f(ab[0], ab[1], ab[2]);
	
	float d = ab.dot(Vector3f(0, 0, 0) - a)/ mCoeABC.norm();
	mCoeABC.normalize();
	mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);
	//float d = ab.dot(Vector3f(0, 0, 0) - a);
	//mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);

	MyMesh::VertexHandle vertex_i;
	MyMesh::HalfedgeHandle heh;
	MyMesh::Point st, at;
	bool ini = 0;
	for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mHandleBegin); voh_it.is_valid(); ++voh_it) {
		heh = mesh.halfedge_handle(voh_it->idx());
		vertex_i = mesh.to_vertex_handle(heh);
		st = mesh.point(vertex_i);
		if (st[2] > mVertexStart[2]) { // 沿着z轴方向往上走
			if (NextHalfEdgeJudge(mesh.next_halfedge_handle(heh))) {
				ini = 1;
				heh = mesh.next_halfedge_handle(heh);
				IterationHalfEdge(heh);
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

void SurfaceCoe::OutlineXCoe(){
	MyMesh::Point k;
	k = mVertexStart;
	float af=0;
	for (int i = 1; i <= mIth[0]; i++) {
		af += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		mOutline2[i].x = af / mLen[0];
	}
	k = mVertexMid; af = 0;
	for (int i = mIth[1]-1; i >= mIth[0]; i--) {
		af += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		mOutline2[i].x = af / mLen[1];
	}
}

void SurfaceCoe::OutlineXCoe2() {
	MyMesh::Point k;
	k = mVertexStart;
	float af = 0;
	for (int i = 1; i <= mIth[0]; i++) {
		af += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		mOutline2[i].x = af / mLen[0];
	}
	k = mVertexMid; af = 0;
	for (int i = mOutline2.size()-1; i >= mIth[1]; i--) {
		af += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		mOutline2[i].x = af / mLen[1];
	}
}

void SurfaceCoe::OutlineExpansion(float tar) {
	vector<MyOutNormal> bso = mOutline2;
	float s = 0, li = 0;
	li = tar>0 ? -tar : tar;

	while (abs(abs(s - mLength) - abs(tar)) > ADDDIFFERENCE) {
		bso = mOutline2;
		for (int j = 0; j < bso.size(); j++) {
			bso[j].a += bso[j].n*bso[j].x*li;
		}
		s = TotalLengh(bso);
		li *= abs(tar) / abs(s - mLength);
	}
	for (int j = 0; j < mOutline2.size(); j++) {
		mOutline2[j].a += mOutline2[j].n*mOutline2[j].x*li;
	}
	mExtension= tar>0 ? -li : li;
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

SurfaceCoe* SurfaceCoe::FindMetara(MyMesh::VertexHandle start, MyMesh::VertexHandle end) {
	SurfaceCoe sfc(start, mesh.point(end), mesh);
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
	//SurfaceCoe sfc(axie,mHandleBegin);
}

void SurfaceCoe::OutlineXCoe(float x ,vector<MySurCutArry> &arryx) {
	MyMesh::Point k = mOutline2[mIth[2]].a; float ss = 0;
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
	interv = mIth[2] / CUTSECTION2; //把后一段分成22段，给出21个截面
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

Vector3f SurfaceCoe::TempVector() {
	return Vector3f(mVertexStart[0]-38.5335, mVertexStart[1]+0.4859, mVertexStart[2]-157.5532);
}
