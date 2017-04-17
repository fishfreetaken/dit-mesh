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

//void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arr, SurfaceCoe* met) {
//
//	MyMesh::Normal m;
//	MyMesh::Point p,p1,p2;// n为递增向量
//	float s1, s2; float lin=0;
//	int k = mesh.n_vertices();
//	int oc = 0;
//	for (MyMesh::VertexIter v_it = mesh.vertices_begin();v_it != mesh.vertices_end(); ++v_it)
//	{
//		oc++;
//		if (!(oc%10000)) {
//			cout << "now is :" << (oc*100/k)<<"%"<< endl;
//		}
//		if (mesh.normal(*v_it)[2] > 0) {
//			continue;
//		}
//		p = mesh.point(*v_it);
//		//if (p[1] < 1) { //for debug;
//		//	int ccd = 000;
//		//}
//		if (met->DistSurface(p)<0) {
//			for (int i = 0; i < CUTSECTION - 1; i++) {
//				lin = arr[i]->DistSurface(p);
//				if ( lin< 0) {
//					if (i == (CUTSECTION - 2)) {
//						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
//					}
//					continue;
//				}else {
//					if (lin > 0) {
//						p2 = arr[i]->FindNearestPoint(p, s2);
//						p1 = i ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1);
//						p+= p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
//					}else {
//						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
//					}
//					break;
//				}
//			}
//		}else if (met->DistSurface(p) > 0) {
//			for (int i = (CUTSECTION-1); i < arr.size(); i++) {
//				lin = arr[i]->DistSurface(p);
//				if (lin> 0) {
//					if (i == (arr.size() - 1)) {
//						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
//					}
//					continue;
//				}
//				else {
//					if (lin < 0) {
//						p2 = arr[i]->FindNearestPoint(p, s2);
//						p1 = i>(CUTSECTION-1) ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1);
//						p+= p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
//					}else {
//						p+= arr[i]->FindNearestPoint(p, s1); //s==0;
//					}
//					break;
//				}
//			}
//		}else {
//			p+=met->FindNearestPoint(p,s1); //s==0;
//		}
//		mesh.set_point(*v_it, p);
//	}
//}
void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arr, SurfaceCoe* met) {
	//MyMesh::Normal m;
	MyMesh::Point p, p1, p2,pro;// n为递增向量
	float s1, s2,sa=0,sb=0; float lin = 0;
	int k = mesh.n_vertices();
	int oc = 0;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		oc++;
		if (!(oc % 10000)) {
			cout << "now is :" << (oc * 100 / k) << "%" << endl;
		}
		if (mesh.normal(*v_it)[2] > 0) {
			continue;
		}
		p = mesh.point(*v_it);
		//m = mesh.normal(*v_it);
		//m.normalize();
		if (met->DistSurface(p)<0) {
			for (int i = 0; i < CUTSECTION - 1; i++) {
				lin = arr[i]->DistSurface(p);
				if (lin< 0) {
					if (i == (CUTSECTION - 2)) {
						pro= arr[i]->FindNearestPoint(p, s1,sa); //s==0;
						p += pro*sa;
					}
					continue;
				}
				else {
					if (lin > 0) {
						p2 = arr[i]->FindNearestPoint(p, s2,sb);
						p1 = i ? arr[i - 1]->FindNearestPoint(p, s1,sa) : met->FindNearestPoint(p, s1,sa);
						pro= p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
						pro.normalize();
						sa= sa*(s2 / (s1 + s2)) + sb*(s1 / (s1 + s2));
						p += sa*pro;
					}
					else {
						pro= arr[i]->FindNearestPoint(p, s1,sa); //s==0;
						p += sa*pro;
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
						pro= arr[i]->FindNearestPoint(p, s1,sa); //s==0;
						p += sa*pro;
					}
					continue;
				}
				else {
					if (lin < 0) {
						p2 = arr[i]->FindNearestPoint(p, s2,sb);
						p1 = i>(CUTSECTION - 1) ? arr[i - 1]->FindNearestPoint(p, s1) : met->FindNearestPoint(p, s1,sb);
						pro= p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
						pro.normalize();
						sa = sa*(s2 / (s1 + s2)) + sb*(s1 / (s1 + s2));
						p += sa*pro;
					}
					else {
						pro= arr[i]->FindNearestPoint(p, s1,sa); //s==0;
						p += pro*sa;
					}
					break;
				}
			}
		}
		else {
			pro= met->FindNearestPoint(p, s1,sa); //s==0;
			p += pro*sa;
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

void SurfaceCoe::AllocateXCoe(float tar){
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error!" << endl;
		return;
	}
	for (int i = 1; i <= mIth[0]; i++) {
		mOutline2[i].x = mOutline2[i].d / mLen[0];
	}
	for (int i = mIth[1]-1; i >= mIth[0]; i--) {
		mOutline2[i].x = abs(mOutline2[i].d-mLen[0]-mLen[1]) / mLen[1];
	}
	mExtension = tar>0 ? -tar : tar;
	OutlineExpansion();
}

void SurfaceCoe::AllocateXCoe2() {
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error2!" << endl;
		return;
	}
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = (mLen[0]-mOutline2[i].d) / mLen[0];
	}
	for (int i = mOutline2.size()-1; i >= mIth[1]; i--) {
		mOutline2[i].x = (mOutline2[i].d-mLen[1]-mLen[0]) / mLen[2];
	}
	mExtension = mX >0 ? -mLength*mX : mLength*mX; //由扩散系数转换为扩散距离
	OutlineExpansion();
}

void SurfaceCoe::OutlineExpansion() {  //(float tar)
	vector<MyOutNormal> bso = mOutline2;
	float s = 0, li = mExtension,pp=0;
	while (abs(abs(s - mLength) - abs(mExtension)) > ADDDIFFERENCE) {
		bso = mOutline2;
		for (int j = 0; j < bso.size(); j++) {
			bso[j].a += bso[j].n*bso[j].x*li;
		}
		s = TotalLengh(bso);
		pp = abs(mExtension)/(s - mLength);
		li *= pp > 0 ? pp : pp<-1? -pp: -1/pp;
	}
	MyMesh::Point aa,bb;
	for (int j = 0; j < mOutline2.size(); j++) {
		mOutline2[j].f = mOutline2[j].n*mOutline2[j].x*li;
		mOutline2[j].a += mOutline2[j].f;
		mOutline2[j].k = li*mOutline2[j].x;
		/*bb = mOutline2[j].a;
		aa= mOutline2[j].n*mOutline2[j].x*li;
		bb[0] = bb[0] ? bb[0] : 1; 这个也不太行，理论上都不通过！
		bb[1] = bb[1] ? bb[1] : 1;
		bb[2] = bb[2] ? bb[2] : 1;
		mOutline2[j].pro = MyMesh::Point(aa[0] / bb[0], aa[1] / bb[1], aa[2] / bb[2]);*/
	}
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

Vector3f SurfaceCoe::TempVector() {
	return Vector3f(mVertexStart[0]-38.5335, mVertexStart[1]+0.4859, mVertexStart[2]-157.5532);
}

MyMesh::Point SurfaceCoe::FindNearestPoint(MyMesh::Point a,float &s) {
	float ss = 0,mig= MAXIMUMX;
	int k = 0;
	for (int i = 0; i < mOutline2.size();i++) {
		ss = DistPoints(a, mOutline2[i].a);
		if (ss < mig) {
			mig = ss;
			k = i;
		}
	}
	s = mig;
	/*if (mOutline2[k].x == 0) {
		return MyMesh::Point(0, 0, 0);
	}*/
	return mOutline2[k].n;  // mOutline2[k].f
}
MyMesh::Point SurfaceCoe::FindNearestPoint(MyMesh::Point a, float &s,float &sc) {
	float ss = 0, mig = MAXIMUMX;
	int k = 0;
	for (int i = 0; i < mOutline2.size(); i++) {
		ss = DistPoints(a, mOutline2[i].a);
		if (ss < mig) {
			mig = ss;
			k = i;
		}
	}
	s = mig;
	/*if (mOutline2[k].x == 0) {
	return MyMesh::Point(0, 0, 0);
	}*/
	sc=mOutline2[k].k;
	return mOutline2[k].n;  // mOutline2[k].f
}