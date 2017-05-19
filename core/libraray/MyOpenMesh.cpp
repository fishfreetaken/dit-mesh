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
	MyMesh::Point pp;
	cout << "Finding Tri Nearest Vertex..." << endl;
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

void MyOpenMesh::InitTriPoint(MyMesh::VertexHandle *p) {
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	v_it_s[3];
	float min[3] = { 0,0,9999 };

	MyMesh::Point pp, es(0, 0, 0);
	cout << "Init Tri Nearest Vertex..." << endl;

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		pp = mesh.point(*v_it);
		if (min[0] < pp[0]) { //x轴向最远的
			min[0] = pp[0];
			v_it_s[0] = v_it;
		}
		if (min[1] < pp[2]) { //z轴向最高的
			min[1] = pp[2];
			v_it_s[1] = v_it;
		}
		if (min[2] > pp[0]) { //x轴向最近的
			min[2] = pp[0];
			v_it_s[2] = v_it;
		}
		//es += pp; //for test
	}
	*(p) = mesh.vertex_handle(v_it_s[0]->idx());
	*(p + 1) = mesh.vertex_handle(v_it_s[1]->idx());
	*(p + 2) = mesh.vertex_handle(v_it_s[2]->idx());

	/*cout << "P1: " << mesh.point(*p) << endl;
	cout << "P2: " << mesh.point(*(p + 1)) << endl;
	cout << "P3: " << mesh.point(*(p + 2)) << endl;*/
	//int init = 0;
	MyMesh::Point tri[3];
	tri[0] = mesh.point(*p);
	tri[1] = mesh.point(*(p + 1));
	tri[2] = mesh.point(*(p + 2));
	//if (abs(tri[0][1]) > 0.5) {
		tri[0][1] = 0;
		//init++;
	//}
	//if (abs(tri[1][1]) > 0.5) {
		tri[1][1] = 0;
		///init++;
	//}
	//if (abs(tri[2][1]) > 0.5) {
		tri[2][1] = 0;
	//	init++;
	//}
	//if (!init) {
	//	return;
	//}
	FindNearest(tri[0], tri[1], tri[2], p);

	cout << "P1: " << mesh.point(*p) << endl;
	cout << "P2: " << mesh.point(*(p + 1)) << endl;
	cout << "P3: " << mesh.point(*(p + 2)) << endl;
	return;

	/*es = es / mesh.n_vertices(); //重心位置-
	cout <<"ES: "<<es << endl;	
	float lin,cp=0;
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it){
		lin = (es - mesh.point(*v_it)).norm();
		if (lin > cp) {
			cp = lin;
			v_it_s[0] = v_it;
		}
	}
	cout << "PX: " << mesh.point(mesh.vertex_handle(v_it_s[0]->idx())) << endl;*/
}

MyMesh::VertexHandle MyOpenMesh::FindNearest(MyMesh::Point a) {  //2
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	v_it_s;
	float min = 9999;
	float lin;
	MyMesh::Point point[3], pp;
	cout << "Finding Signle Nearest Vertex..." << endl;
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		pp = mesh.point(*v_it);
		lin = (pp - a).norm();
		if (lin < min) {
			min = lin;
			v_it_s = v_it;
		}
	}
	return mesh.vertex_handle(v_it_s->idx());
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
void MyOpenMesh::BotIteration(set<MyOutBottom>&arr, int idx ,int iver) {
	if (iver > ITERATIONCISHU) {
		return;
	}
	iver++;
	MyMesh::VertexVertexIter vv_it = mesh.vv_iter(mesh.vertex_handle(idx));
	for (; vv_it.is_valid(); ++vv_it)
	{
		BotIteration(arr,vv_it->idx(),iver);
		MyOutBottom git;
		git.s = mesh.normal(*vv_it);
		git.s.normalize();
		git.n = mesh.vertex_handle(vv_it->idx());
		arr.insert(git);
	}
	return;
}

void MyOpenMesh::ShoeBottomLine(vector<MyOutBottom>&arrx) {
	cout << "Out the Bottom Line" << endl;
	//vector<MyOutBottom>arrx;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		MyMesh::Normal es(0, 0, 0), dif(0, 0, 0);
		
#ifdef SHOEBOTTOMLL
		set<MyOutBottom>arr; //使用vector整体提升的很高很高，之间差别特别明显
#else
		vector<MyOutBottom>arr; //使用vector整体提升的很高很高，之间差别特别明显
#endif // DEBUG

		MyOutBottom git;

		for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			for (MyMesh::VertexVertexIter vv_it1 = mesh.vv_iter(mesh.vertex_handle(vv_it->idx())); vv_it1.is_valid(); ++vv_it1)
			{
				for (MyMesh::VertexVertexIter vv_it2 = mesh.vv_iter(mesh.vertex_handle(vv_it1->idx())); vv_it2.is_valid(); ++vv_it2)
				{
					for (MyMesh::VertexVertexIter vv_it3 = mesh.vv_iter(mesh.vertex_handle(vv_it2->idx())); vv_it3.is_valid(); ++vv_it3)
					{
						for (MyMesh::VertexVertexIter vv_it4 = mesh.vv_iter(mesh.vertex_handle(vv_it3->idx())); vv_it4.is_valid(); ++vv_it4)
						{
							git.s = mesh.normal(*vv_it4);
							git.n = mesh.vertex_handle(vv_it4->idx());
#ifdef SHOEBOTTOMLL
							arr.insert(git);
#else
							arr.push_back(git);
#endif // DEBUG
							
						}
						git.s = mesh.normal(*vv_it3);
						git.n = mesh.vertex_handle(vv_it3->idx());
#ifdef SHOEBOTTOMLL
						arr.insert(git);
#else
						arr.push_back(git);
#endif // DEBUG
					}
					git.s = mesh.normal(*vv_it2);
					git.n = mesh.vertex_handle(vv_it2->idx());
#ifdef SHOEBOTTOMLL
					arr.insert(git);
#else
					arr.push_back(git);
#endif // DEBUG
				}
				git.s = mesh.normal(*vv_it1);
				git.n = mesh.vertex_handle(vv_it1->idx());
#ifdef SHOEBOTTOMLL
				arr.insert(git);
#else
				arr.push_back(git);
#endif // DEBUG
			}
			git.s = mesh.normal(*vv_it);
			git.n = mesh.vertex_handle(vv_it->idx());
#ifdef SHOEBOTTOMLL
			arr.insert(git);
#else
			arr.push_back(git);
#endif // DEBUG
		}

		for (auto i : arr) {
			es += i.s;
		}
		es = es / arr.size();
		for (auto i : arr)
		{
			dif[0] += ((i.s[0] - es[0])*i.x) * ((i.s[0] - es[0])*i.x);
			dif[1] += ((i.s[1] - es[1])*i.x) * ((i.s[1] - es[1])*i.x);
			dif[2] += ((i.s[2] - es[2])*i.x) * ((i.s[2] - es[2])*i.x);
		}
		dif[0] = sqrt(dif[0]);
		dif[1] = sqrt(dif[1]);
		dif[2] = sqrt(dif[2]);

		MyOutBottom stt;
		stt.n = mesh.vertex_handle(v_it->idx());
		stt.s = dif;
		stt.a = mesh.point(*v_it);
		arrx.push_back(stt);
	}
}

void MyOpenMesh::FindFloorContour(vector<MyMesh::Point> &out) {
	vector<int> floorContour;
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	v_it_s[3];
	MyMesh::Point curPoint;
	double minX, maxX, minY, maxY;
	cout << "Finding Floor Contour..." << endl;
	v_it = mesh.vertices_begin();
	minX = maxX = mesh.point(*v_it)[0];
	minY = maxY = mesh.point(*v_it)[1];
	for (; v_it != v_end; v_it++)
	{
		curPoint = mesh.point(*v_it);
		if (curPoint[0] > maxX)
			maxX = curPoint[0];
		if (curPoint[0] < minX)
			minX = curPoint[0];
		if (curPoint[1] > maxY)
			maxY = curPoint[1];
		if (curPoint[1] < minY)
			minY = curPoint[1];
		//mesh.set_point(*v_it, curPoint);
		//v_it->idx();
	}
	double xSampleS = 0.5;//0.05;
	int xSampleL = int((maxX - minX) / xSampleS + 1);
	vector<vector<int>> xSamplePtV(xSampleL);
	v_it = mesh.vertices_begin();
	for (; v_it != v_end; v_it++)
	{
		curPoint = mesh.point(*v_it);
		int xIndex = int((curPoint[0] - minX) / xSampleS);
		xSamplePtV[xIndex].push_back(v_it->idx());
	}
	floorContour.clear();
	for (int i = 0; i < xSampleL; i++)
	{
		if (xSamplePtV[i].size() == 0)
			continue;
		double minZ = mesh.point(mesh.vertex_handle(xSamplePtV[i][0]))[2];
		for (unsigned int j = 1; j < xSamplePtV[i].size(); j++)
		{
			curPoint = mesh.point(mesh.vertex_handle(xSamplePtV[i][j]));
			if (curPoint[2] < minZ)
				minZ = curPoint[2];
		}
		minY = maxY = mesh.point(mesh.vertex_handle(xSamplePtV[i][0]))[1];
		int left = 0;
		int right = 0;
		for (unsigned int j = 0; j < xSamplePtV[i].size(); j++)
		{
			curPoint = mesh.point(mesh.vertex_handle(xSamplePtV[i][j]));
			//floorContour.push_back(curPoint);

			if ((curPoint[2] - minZ) > 5)
				continue;
			if (curPoint[1] < minY)
			{
				minY = curPoint[1];
				left = j;
			}
			if (curPoint[1] > maxY)
			{
				maxY = curPoint[1];
				right = j;
			}
		}
		floorContour.push_back(xSamplePtV[i][left]);
		floorContour.push_back(xSamplePtV[i][right]);
	}
	//vector<cv::Point2d> cloudMat;
	//vector<int> edgeL;
	/*for (unsigned int i = 0; i < floorContour.size(); i++)
	{
		curPoint = mesh.point(mesh.vertex_handle(floorContour[i]));
		cloudMat.push_back(cv::Point2d(curPoint[0], curPoint[1]));
	}*/

	for (auto i : floorContour) {
		out.push_back(mesh.point(mesh.vertex_handle(i)));
	}
}

void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arr, vector<MyMesh::Point>&css) {
	MyMesh::Point  p,p1,p2;// n为递增向量
	float s1, s2; int sta=0,end=0;
	vector<MyMesh::Point>  result;

	float arrlen = arr.size();
	for (int j=0;j<css.size();j++)// p :css)
	{
		p = css[j];
		int i = arrlen / 2;
		end = arrlen;
		sta = 0;
		while (1) {
			if (arr[i]->DistSurface(p) > 0) {
				if (arr[i - 1]->DistSurface(p) > 0) {
					if (i == 1) {
						p += arr[i - 1]->FindNearestPoint(p, s1); //s==0;
						break;
					}
					end = i;
					i -= (i - sta) / 2;
				}
				else if (arr[i - 1]->DistSurface(p)<0) {
					p2 = arr[i]->FindNearestPoint(p, s2);
					p1 = arr[i - 1]->FindNearestPoint(p, s1);
					p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					break;
				}
				else {
					p += arr[i - 1]->FindNearestPoint(p, s1); //s==0;
					break;
				}
			}
			else if (arr[i]->DistSurface(p) < 0) {
				if (arr[i + 1]->DistSurface(p) < 0) {
					if (i == (arrlen - 2)) {  //1
						p += arr[i + 1]->FindNearestPoint(p, s1); //s==0;
						break;
					}
					sta = i;
					i += (end - i) / 2;
				}
				else if (arr[i + 1]->DistSurface(p)>0) {
					p2 = arr[i]->FindNearestPoint(p, s2);
					p1 = arr[i + 1]->FindNearestPoint(p, s1);
					p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					cout << s1 << " : " << p1 << " || " << s2 << " : " << p2 << " || "<< p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2)) << endl;
					break;
				}
				else {
					p += arr[i + 1]->FindNearestPoint(p, s1); //s==0;
					break;
				}
			}
			else {
				p += arr[i]->FindNearestPoint(p, s1); //s==0;
				break;
			}
		}
		result.push_back(p);
	}
	for (int i = 0; i < css.size(); i++) {
		cout << result[i] - css[i] << " : "<< (result[i] - css[i]).norm()<<" : ";
		if (i >= 1) {
			cout << (result[i]-result[i-1]).norm();
		}
		cout << endl;
	}
	
	return;
}

void MyOpenMesh::ShoeAddLength(MyMesh::Point start, SurfaceCoe*met ,float exp) {
	MyMesh::Point p;
	float max = abs(met->DistSurface(start));
	float lin;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p = mesh.point(*v_it);
		lin = met->DistSurface(p);
		if (lin <= 0) {
			continue;
		}
		p[0] += exp*lin / max; //x轴方向移动
		mesh.set_point(*v_it, p);
	}
}

void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arrx) {

	vector<SurfaceCoe*>&arr=arrx;
	cout << "Now is shoe Expansing..." << endl;
	MyMesh::Point p, p1, p2;// n为递增向量
	float s1, s2;

	int k = mesh.n_vertices();
	int oc = 0;
	int arrlen = arr.size();

	int sta = 0, end = 0;

	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		oc++;
		if (!(oc % 10000)) {
			cout << "now is :" << (oc * 100 / k) << "%" << endl;
		}
		p = mesh.point(*v_it);

		int i = arrlen / 2;
		end = arrlen;
		sta = 0;
		while (1) {
			if (arr[i]->DistSurface(p) > 0) {
				if (arr[i-1]->DistSurface(p) > 0) {
					if (i == 1) {
						p += arr[i - 1]->FindNearestPoint(p, s1); //s==0;
						break;
					}
					end = i;
					i -= (i-sta) / 2;
				}else if(arr[i-1]->DistSurface(p)<0){
					p2 = arr[i]->FindNearestPoint(p, s2);
					p1 = arr[i - 1]->FindNearestPoint(p, s1);
					p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					break;
				}else {
					p += arr[i-1]->FindNearestPoint(p, s1); //s==0;
					break;
				}
			}
			else if (arr[i]->DistSurface(p) < 0) {
				if (arr[i + 1]->DistSurface(p) < 0) {
					if (i == (arrlen - 2)) {  //1
						p += arr[i + 1]->FindNearestPoint(p, s1); //s==0;
						break;
					}
					sta = i;
					i += (end - i) / 2;
				}else if (arr[i + 1]->DistSurface(p)>0) {
					p2 = arr[i]->FindNearestPoint(p, s2);
					p1 = arr[i + 1]->FindNearestPoint(p, s1);
					p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
					break;
				}else {
					p += arr[i+1]->FindNearestPoint(p, s1); //s==0;
					break;
				}
			}else {
				p += arr[i]->FindNearestPoint(p, s1); //s==0;
				break;
			}
		}
		mesh.set_point(*v_it, p);
	}
}

void MyOpenMesh::ShoeExpansionWist(vector<SurfaceCoe*> &arr) {
	cout << "Now is shoe wist Expansing..." << endl;
	MyMesh::Point p, p1, p2;// n为递增向量
	float s1, s2;

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
		if (arr[0]->DistSurface(p) >= 0) {
			continue;
		}
		if (arr[arr.size()-1]->DistSurface(p) <= 0) {
			continue;
		}

		for (int i = 1; i < arrlen; i++) {
			if (arr[i]->DistSurface(p) > 0) {
				p2 = arr[i]->FindNearestPoint(p, s2);
				p1 = arr[i - 1]->FindNearestPoint(p, s1);
				p += p1*(s2 / (s1 + s2)) + p2*(s1 / (s1 + s2));
				break;
			}
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
SurfaceCoe::SurfaceCoe(MyMesh::Point mid, MyMesh::Point end, MyMesh::VertexHandle c, MyMesh &d) :
	mesh(d),
	mHandleBegin(c),
	mVertexMid(mid),
	mVertexEnd(end)
{
	mVertexStart = mesh.point(c);
}


bool SurfaceCoe::Init() {//(MyMesh::VertexHandle *vertex,MyMesh &mmesh){
	MyOutNormal abc;
	abc.a = mVertexStart;
	abc.n = mesh.normal(mHandleBegin);
	Vector3f nv(abc.n[0], abc.n[1], abc.n[2]);
	Vector3f ln = nv - mCoeABC*(nv.dot(mCoeABC));
	ln.normalize();
	abc.n = MyMesh::Normal(ln[0], ln[1], ln[2]);

	mOutline2.push_back(abc);

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

	MyOutNormal abc;
	abc.a = mVertexStart;
	abc.n = mesh.normal(mHandleBegin);
	Vector3f nv(abc.n[0], abc.n[1], abc.n[2]);
	Vector3f ln = nv - mCoeABC*(nv.dot(mCoeABC));
	ln.normalize();
	abc.n = MyMesh::Normal(ln[0], ln[1], ln[2]);

	mOutline2.push_back(abc);

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
	/*if (ini) {
		if (cmd) {
			OutlineRefine();
		}
		CoquerMidEnd();
	}*/
	if (ini) {
		switch (cmd) {
			case 1:
				OutlineRefine();
				CoquerMidEnd();
				break;
			case 2:
				//NULL
				mLength= TotalLengh(mOutline2);
				break;
			default:
				CoquerMidEnd();
				break;
		}
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
	int ini = 0;
	for (int i = mOutline2.size() / 8; i < mOutline2.size(); i++) {
		if (DistPoints(mOutline2[i].a, k) < 0.003) {
			j = i;
			//cout << i << endl;
			ini = 1;
			break;
		}
	}
	if (!ini) {
		return;
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
	float s1 = DistPoints(mc.a, mesh.point(a));
	float s2 = DistPoints(mc.a, mesh.point(b));
	mc.n = mesh.normal(a)*(s2/(s1+s2)) + mesh.normal(b)*(s1/(s1+s2));//平面法向量投影
	//mc.n = (mesh.normal(a) + mesh.normal(b)) / 2;//平面法向量投影
	Vector3f nv(mc.n[0], mc.n[1], mc.n[2]);
	Vector3f ln = nv - mCoeABC*(nv.dot(mCoeABC));
	ln.normalize();
	mc.n = MyMesh::Normal(ln[0], ln[1], ln[2]);
	mOutline2.push_back(mc);
}

float SurfaceCoe::AllocateXCoe(float tar) {
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error!" << endl;
		return -1;
	}
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = (1 - mOutline2[i].d / mLen[0]);//mLen[0]
	}
	for (int i = mOutline2.size() - 1; i >= mIth[1]; i--) {
		mOutline2[i].x = (mOutline2[i].d - mLen[1] - mLen[0]) / mLen[2];  //mLen[2]
	}

	vector<float>input, output;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		input.push_back(mOutline2[i].x);
	}
	for (int i = 0; i <= mIth[0]; i++) {
		input.push_back(mOutline2[i].x);
	}
	GaussianSmooth(input, output, GAUSSIONFILTERNUM);

	for (int i = 0; i < mOutline2.size() - mIth[1]; i++) {
		mOutline2[i + mIth[1]].x = output[i];
	}
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = output[i + mOutline2.size() - mIth[1]];
	}
	
	mExtension = tar > 0 ? -tar : tar;
	return OutlineExpansion();
}

float SurfaceCoe::TopSlide(SurfacePure*met) {
	if (abs(met->DistSurface(mVertexStart)) > TOPOFFSET) {
		return 0;
	}

	vector<float>input, output;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		if (abs(met->DistSurface(mOutline2[i].a)) <= TOPOFFSET) {
			mOutline2[i].x = abs(met->DistSurface(mOutline2[i].a))/ TOPOFFSET;
		}
		input.push_back(mOutline2[i].x);
	}
	for (int i = 0; i <= mIth[0]; i++) {
		if (abs(met->DistSurface(mOutline2[i].a)) <= TOPOFFSET) {
			mOutline2[i].x = abs(met->DistSurface(mOutline2[i].a))/ TOPOFFSET;
		}
		input.push_back(mOutline2[i].x);
	}
	GaussianSmooth(input, output,3);

	for (int i = 0; i < mOutline2.size() - mIth[1]; i++) {
		mOutline2[i + mIth[1]].x = output[i];
	}
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = output[i + mOutline2.size() - mIth[1]];
	}
	return 1;
}

float SurfaceCoe::AllocateXCoe(){//(SurfaceCoe*met,MyMesh::Point aa,MyMesh::Point bb) {
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error2!" << endl;
		return 0;
	}
	//float mm = mLen[0] > mLen[2] ? mLen[0] : mLen[2]; //???
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = (1 - mOutline2[i].d / mLen[0]);//mLen[0]
	}
	for (int i = mOutline2.size() - 1; i >= mIth[1]; i--) {
		mOutline2[i].x = (mOutline2[i].d - mLen[1] - mLen[0]) / mLen[2];  //mLen[2]
	}

	vector<float>input,output;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		input.push_back(mOutline2[i].x);
	}
	for (int i = 0; i <= mIth[0]; i++) {
		input.push_back(mOutline2[i].x);
	}
	GaussianSmooth(input,output, GAUSSIONFILTERNUM);

	for (int i = 0; i < mOutline2.size() - mIth[1]; i++) {
		mOutline2[i + mIth[1]].x = output[i];
	}
	for (int i = 0; i <=mIth[0] ; i++) {
		mOutline2[i].x = output[i + mOutline2.size() - mIth[1]];
	}

	mExtension = mLength*mX; //由扩散系数转换为扩散距离

	return OutlineExpansion();
}

float SurfaceCoe::AllocateXCoe(SurfacePure*met) {//(SurfaceCoe*met,MyMesh::Point aa,MyMesh::Point bb) {
	if ((!mLen[0]) || (!mLen[1])) {
		cout << "zero error2!" << endl;
		return 0;
	}
	//float mm = mLen[0] > mLen[2] ? mLen[0] : mLen[2]; //???
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = (1 - mOutline2[i].d / mLen[0]);//mLen[0]
	}
	for (int i = mOutline2.size() - 1; i >= mIth[1]; i--) {
		mOutline2[i].x = (mOutline2[i].d - mLen[1] - mLen[0]) / mLen[2];  //mLen[2]
	}

	vector<float>input, output;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		input.push_back(mOutline2[i].x);
	}
	for (int i = 0; i <= mIth[0]; i++) {
		input.push_back(mOutline2[i].x);
	}
	GaussianSmooth(input, output, GAUSSIONFILTERNUM);

	for (int i = 0; i < mOutline2.size() - mIth[1]; i++) {
		mOutline2[i + mIth[1]].x = output[i];
	}
	for (int i = 0; i <= mIth[0]; i++) {
		mOutline2[i].x = output[i + mOutline2.size() - mIth[1]];
	}
	mExtension = mLength*mX; //由扩散系数转换为扩散距离

	TopSlide(met);
	mExtension = mLength*mX;

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
	//li += (met->re() - li)*(abs(met->DistSurface(mVertexStart)) / ls);
	for (int j = 0; j < mOutline2.size(); j++) {
		mOutline2[j].f = mOutline2[j].n*mOutline2[j].x*li;
		//mOutline2[j].a += mOutline2[j].f;//其实都是未移动的点
		mOutline2[j].m = mOutline2[j].a + mOutline2[j].f;//用于输出调试
	}
	mExtensionli = li;
	return li;
}
class CircleIndex {
public:
	CircleIndex(MyMesh::Normal a,int b) :
	mN(a),
	Ith(b)
	{}
	MyMesh::Normal RtNormal() { return mN; }
	int RtIth() { return Ith; }

	CircleIndex *next;
	CircleIndex *up;
private:
	MyMesh::Normal mN;
	int Ith;
};

void SurfaceCoe::InitMidEndPoint(vector<MyMesh::Point>&fw) {
	int window = 7;
	float upstep = 7;//默认

	MyMesh::Point bua(0,0,0), bub(0,0,0);
	int bot = 0; float lin = mOutline2[0].a[2];
	float topv= mOutline2[0].a[2];
	for (int i = 1; i < mOutline2.size(); i++) {
		if (lin > mOutline2[i].a[2]) {
			lin = mOutline2[i].a[2];
			bot = i;
		}
		if (topv < mOutline2[i].a[2]) {
			topv= mOutline2[i].a[2];
		}
	}//默认情况下bot应该处于500左右范围内的值(该值远远大于window窗口滤波范围)
	upstep = abs((topv - lin) / 7);
	//cout << upstep << endl;
	set<MyBotOutLine> sum1,sum2;//for debug
	lin = 0; int ith=0;
	for (int i = bot; i < mOutline2.size(); i++) {
		if (abs(mOutline2[i].a[2] - mOutline2[bot].a[2]) > upstep) {
			break;
		}
		vector<MyMesh::Normal>sm;
		sm.push_back(mOutline2[i].n);
		for (int j = 1; j < window; j++) {
			sm.push_back(mOutline2[i - j].n);
			sm.push_back(mOutline2[i + j].n);
		}
		MyMesh::Normal es(0, 0, 0), dif(0, 0, 0);
		for (auto s : sm) {
			es += s;
		}
		es = es / (window * 2);
		for (auto s : sm) {
			dif[0] += ((s[0] - es[0])) * ((s[0] - es[0]));
			dif[1] += ((s[1] - es[1])) * ((s[1] - es[1]));
			dif[2] += ((s[2] - es[2])) * ((s[2] - es[2]));
		}
		dif[0] = sqrt(dif[0]);
		dif[1] = sqrt(dif[1]);
		dif[2] = sqrt(dif[2]);
		MyBotOutLine git;
		git.s = dif;
		git.x = dif.norm();
		git.ith = i;
		sum1.insert(git);
		if (lin < dif.norm()) {
			lin = dif.norm();
			ith = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return;
	}
	if (mOutline2[ith].a[1]>0) {
		bua = mOutline2[ith].a;
	}
	else {
		bub = mOutline2[ith].a;
	}//step one over!

	lin = 0; ith = 0;
	for (int i = bot; i >=0; i--) {
		if (abs(mOutline2[i].a[2] - mOutline2[bot].a[2]) > upstep) {
			break;
		}
		vector<MyMesh::Normal>sm;
		sm.push_back(mOutline2[i].n);
		for (int j = 1; j < window; j++) {
			sm.push_back(mOutline2[i - j].n);
			sm.push_back(mOutline2[i + j].n);
		}
		MyMesh::Normal es(0, 0, 0), dif(0, 0, 0);
		for (auto i : sm) {
			es += i;
		}
		es = es / (window * 2);
		for (auto i : sm) {
			dif[0] += ((i[0] - es[0])) * ((i[0] - es[0]));
			dif[1] += ((i[1] - es[1])) * ((i[1] - es[1]));
			dif[2] += ((i[2] - es[2])) * ((i[2] - es[2]));
		}
		dif[0] = sqrt(dif[0]);
		dif[1] = sqrt(dif[1]);
		dif[2] = sqrt(dif[2]);
		MyBotOutLine git;
		git.s = dif;
		git.x = dif.norm();
		git.ith = i;
		sum2.insert(git);
		if (lin < dif.norm()) {
			lin = dif.norm();
			ith = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return;
	}
	if (mOutline2[ith].a[1]>0) {
		bua = mOutline2[ith].a;
	}
	else {
		bub = mOutline2[ith].a;
	}//step two is on running;
	if (bua.norm()) {
		fw.push_back(bua);
	}
	else {
		fw.push_back(mOutline2[bot-mOutline2.size()/8].a);
	}
	if (bub.norm()) {
		fw.push_back(bub);
	}
	else {
		for (int i = bot+1; i < mOutline2.size(); i++) {
			if (abs(mOutline2[bot].a[2] - mOutline2[i].a[2]) > 1) {
				ith = i;
				break;
			}
		}
		fw.push_back(mOutline2[ith].a);
	}
}

void SurfaceCoe::InitMidTopPoint(vector<MyMesh::Point>&fw) {
	int window = 7;
	float upstep = 7;//默认

	MyMesh::Point bua(0,0,0), bub(0,0,0);
	float lin = mOutline2[0].a[2];
	float topv= mOutline2[0].a[2];
	for (int i = 1; i < mOutline2.size(); i++) {
		if (lin < mOutline2[i].a[2]) {
			lin = mOutline2[i].a[2];
			//bot = i;
		}
		if (topv > mOutline2[i].a[2]) {
			topv= mOutline2[i].a[2];
		}
	}//默认情况下bot应该处于500左右范围内的值(该值远远大于window窗口滤波范围)
	upstep = abs((topv - lin) / 8);
	//cout << upstep << endl;

	set<MyBotOutLine> sum1, sum2;//for debug
	vector<MyBotOutLine>tf;
	for (int i = mOutline2.size() * 2 / 3; i < mOutline2.size(); i++) {
		MyBotOutLine gt;
		gt.ith = i;
		gt.s = mOutline2[i].n;
		gt.a = mOutline2[i].a;
		tf.push_back(gt);
	}
	for (int i = 0; i < mOutline2.size()/3; i++) {
		MyBotOutLine gt;
		gt.ith = i;
		gt.s = mOutline2[i].n;
		gt.a = mOutline2[i].a;
		tf.push_back(gt);
	}
	int bot = 0; lin = tf[0].a[2];
	for (int i = 0; i < tf.size(); i++) {
		if (lin < tf[i].a[2]) {
			lin = tf[i].a[2];
			bot = i;
		}
	}
	
	lin = 0; int ith=0;
	for (int i = bot; i < tf.size(); i++) {
		if (abs(tf[i].a[2] - tf[bot].a[2]) > upstep) {
			break;
		}
		vector<MyMesh::Normal>sm;
		sm.push_back(mOutline2[i].n);
		for (int j = 1; j < window; j++) {
			sm.push_back(tf[i - j].s);
			sm.push_back(tf[i + j].s);
		}
		MyMesh::Normal es(0, 0, 0), dif(0, 0, 0);
		for (auto s : sm) {
			es += s;
		}
		es = es / (window * 2);
		for (auto s : sm) {
			dif[0] += ((s[0] - es[0])) * ((s[0] - es[0]));
			dif[1] += ((s[1] - es[1])) * ((s[1] - es[1]));
			dif[2] += ((s[2] - es[2])) * ((s[2] - es[2]));
		}
		dif[0] = sqrt(dif[0]);
		dif[1] = sqrt(dif[1]);
		dif[2] = sqrt(dif[2]);
		MyBotOutLine git;
		git.s = dif;
		git.x = dif.norm();
		git.ith = i;
		sum1.insert(git);
		if (lin < dif.norm()) {
			lin = dif.norm();
			ith = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return;
	}
	if (tf[ith].a[1]>0) {
		bua = tf[ith].a;
	}
	else {
		bub = tf[ith].a;
	}//step one over!

	lin = 0; ith = 0;
	for (int i = bot; i >=0; i--) {
		if (abs(tf[i].a[2] - tf[bot].a[2]) > upstep) {
			break;
		}
		vector<MyMesh::Normal>sm;
		sm.push_back(mOutline2[i].n);
		for (int j = 1; j < window; j++) {
			sm.push_back(tf[i - j].s);
			sm.push_back(tf[i + j].s);
		}
		MyMesh::Normal es(0, 0, 0), dif(0, 0, 0);
		for (auto i : sm) {
			es += i;
		}
		es = es / (window * 2);
		for (auto i : sm) {
			dif[0] += ((i[0] - es[0])) * ((i[0] - es[0]));
			dif[1] += ((i[1] - es[1])) * ((i[1] - es[1]));
			dif[2] += ((i[2] - es[2])) * ((i[2] - es[2]));
		}
		dif[0] = sqrt(dif[0]);
		dif[1] = sqrt(dif[1]);
		dif[2] = sqrt(dif[2]);
		MyBotOutLine git;
		git.s = dif;
		git.x = dif.norm();
		git.ith = i;
		sum2.insert(git);
		if (lin < dif.norm()) {
			lin = dif.norm();
			ith = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return;
	}
	if (tf[ith].a[1]>0) {
		bua = tf[ith].a;
	}
	else {
		bub = tf[ith].a;
	}//step two is on running;
	if (bua.norm()) {
		fw.push_back(bua);
	}
	if (bub.norm()) {
		fw.push_back(bub);
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
bool SurfaceCoe::OutlineEigenf(const char  *filename) {
	if (!mOutline2.size()) {
		cout << "outline no points" << endl;
		return false;
	}
	FILE *fp;
	fopen_s(&fp, filename, "w");
	for (unsigned int i = 0; i < mOutline2.size(); i++) {
		fprintf(fp, "%d %f,%f,%f  %f,%f,%f\n",i,mOutline2[i].n[0], mOutline2[i].n[1],mOutline2[i].n[2], mOutline2[i].a[0],mOutline2[i].a[1],mOutline2[i].a[2]);
	}
	fclose(fp);
	return true;
}

bool SurfaceCoe::OutlineEigen(vector<Vector4f> *a) {
	if (!mOutline2.size()) {
		cout << "outline no points" << endl;
		return false;
	}
	for (auto i : mOutline2) {
		a->push_back(Vector4f(i.f[0], i.f[1], i.f[2],i.x));
	}
	return true;
}

SurfaceCoe* SurfaceCoe::FindMetara(MyMesh::VertexHandle start, MyMesh::VertexHandle end) {
	SurfaceCoe sfc(start, mesh.point(end), mesh);
	cout << "Now is finding Metara..." << endl;
	float dd = MAXIMUMX, ss;
	int metara;
	for (int i = mIth[0] * 0.13; i < mIth[0] * 0.5; i++) {
		sfc.SetMidPoint(mOutline2[i].a);
		if (sfc.Init(0)) {
			ss = sfc.ReturnLength();
			if (ss < dd) {
				dd = ss;
				metara = i;   //mLen[2]=i;
			}
		}
	}
	cout << "Metara ith:" << metara << endl;
	vector<MyMesh::Point> fm;
	sfc.InitMidEndPoint(fm);
	if (fm.size() < 2) {
		cout << "metara point error" << endl;
		return NULL;
	}
	cout << fm[0] << endl;
	cout << fm[1] << endl;
	MyMesh::VertexHandle sfcstartp = FindNearest(mOutline2[metara].a);
	
	SurfaceCoe *ret = new SurfaceCoe(fm[0], fm[1], sfcstartp, mesh);
	//SurfaceCoe *ret = new SurfaceCoe(MyMesh::Point(160.563766,34.669678,2.144665),MyMesh::Point(146.947586,-44.111824,3.824428),sfcstartp,mesh);//手动给出底边的两个点
	ret->Init(0);
	ret->SetMIth(metara);
	return ret;
}

SurfaceCoe* SurfaceCoe::FindWaistLine(SurfaceCoe *met) {  //腰围和背围可以使用同一个函数，其原理相同
	int ccmid = 0; float ss = MAXIMUMX, lin = 0;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		lin = abs(met->DistSurface(mOutline2[i].a));
		if (lin<ss) {
			ss = lin;
			ccmid = i;
		}
	}
	MyMesh::Point mdp, enp,of;
	int iith = UpOneInch(met->ReturnIth(2),of);
	MyMesh::VertexHandle start = FindNearest(of);
	

	float min = MAXIMUMX; int mid = 0;
	for (int i = mIth[1]+10; i < ccmid ; i++) {
		mdp = enp = mOutline2[i].a;
		enp[1] += 1;
		SurfaceCoe sfc(mdp, enp, start ,mesh);
		if (sfc.Init(2)) {
			lin = sfc.ReturnLength();
			if (lin <min ) {
				min = lin;
				mid = i;   //mLen[2]=i;
			}
		}
	}
	cout << "Waist Ith:" << mid << endl;
	mdp = enp = mOutline2[mid].a;
	enp[1] += 1;
	SurfaceCoe* ret= new SurfaceCoe(enp,mdp, start, mesh);
	ret->SetMIth(iith);
	ret->Init(2);
	return ret;
}

float SurfaceCoe::FindAddLenth(SurfaceCoe *met,float ext) {
	float min = MAXIMUMX, lin; int ith = 0,len=mOutline2.size();
	float max = abs(met->DistSurface(mOutline2[0].a));
	struct ctt {
		float i;
		MyMesh::Point a;
	};
	vector<struct ctt> lm;
	float source=0;
	for (int i = mIth[1]; i < len; i++) {
		lin = met->DistSurface(mOutline2[i].a);
		if (lin >= 0) {
			lin = abs(lin);
			if (lin < min) {
				min = lin;
				ith = i;
			}
			struct ctt st;
			st.i = lin / max;
			st.a = mOutline2[i].a;
			lm.push_back(st);
			if (lm.size() > 1) {
				source+=DistPoints(st.a,lm[lm.size()-2].a);
			}
		}
	}
	float ex = 2 * ext;//for save
	float a = 0;
	float b =0.5;
	float c = 1;
	while (abs(abs(lin -source)- ext)>ADDDIFFERENCE){
		lin = 0;
		vector<struct ctt> lmc=lm;
		for (int i = 0; i < lmc.size(); i++) {
			lmc[i].a[0] += lmc[i].i*ex*b;
			if (i >=1) {
				lin += DistPoints(lmc[i].a, lmc[i-1].a);
			}
		}
		/*MyMesh::Point cc = lmc[0].a;
		for (int i = 1; i < lmc.size(); i++) {
			lin +=DistPoints(cc, lmc[i].a);
			cc = lmc[i].a;
		}*/
		float pp = abs(lin - source) - ext;
		if (pp > 0) {
			c = b;
			b = (b + a) / 2;// c - (c - a) / 2
		}
		else if(pp<0){
			a = b;
			b = (c + b) / 2;
		}
		else {
			break;
		}
	}
	return b*ex;
}

int SurfaceCoe::UpOneInch(int ith, MyMesh::Point &af) {
	MyMesh::Point k = mOutline2[ith].a;
	MyMesh::Point b;
	float dis=0;
	for (int i = ith+1; i < mIth[0]; i++) {
		dis += (mOutline2[i].a - k).norm();
		k = mOutline2[i].a;
		if (dis >= ONEINCHLEN) {
			af= mOutline2[i].a;
			return i;
		}
	}
	return -1;
}

//SurfaceCoe* SurfaceCoe::FindBackLine() {}

MyMesh::VertexHandle SurfaceCoe::FindNearest(MyMesh::Point a) {  //2
	MyMesh::VertexIter  v_it, v_end(mesh.vertices_end());
	MyMesh::VertexIter	v_it_s;
	float min = MAXIMUMX;
	float lin;
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		lin = (mesh.point(*v_it) - a).norm();
		if (lin < min) {
			min = lin;
			v_it_s = v_it;
		}
	}
	return  mesh.vertex_handle(v_it_s->idx());
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

vector<MySurCutArry> SurfaceCoe::OutCutOutline(float exp,SurfaceCoe *met,Vector3f axi) {
	float xi = exp>0 ? -exp / met->ReturnLength():exp/met->ReturnLength();
	struct CutArry2 {
		MyMesh::Point a;
		float x = 1;
		Vector3f n; 
	};
	float ju,max0,max1;
	max0 = met->DistSurface(mVertexStart);
	max1 =abs(met->DistSurface(mVertexMid));

	Vector3f coemet = met->ReturnCoe();
	//float thert=acos(coemet.dot(axi));//原版的，忘记除以2了
	float thert = acos(coemet.dot(axi))/2;
	Vector3f mix = coemet.cross(axi);
	mix.normalize();

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
			//st.x = xi;
			lin = thert*(ju / max0 );
		}else if(ju<0){
			if (!ini) {
				ini = true;
				struct CutArry2 sts;
				sts.x = xi;
				sts.n = met->ReturnCoe();
				sts.a = met->ReturnStartPoint();
				arry.push_back(sts);
			}
			//st.x = xi;
			st.x = xi*(1 - abs(ju) / max1);
			lin = thert*(abs(ju) / max1);
		}else {
			continue;
		}
		Quaternionx mtransfer(cos(lin),sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0,coemet[0],coemet[1],coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(),out.y(),out.z());
		//st.n.normalize();//意义不大，相差不是很多！！
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

vector<MySurCutArry> SurfaceCoe::OutCutOutline(float exp, SurfaceCoe *meta, SurfaceCoe *metb, MyOpenMesh&ios) { //背围
	cout << "back cut out..." << endl;
	float xi = exp>0 ? -exp / metb->ReturnLength() : exp / metb->ReturnLength();
	struct CutArry2 {
		MyMesh::Point a;
		float x = 1;
		float xa = 1;
		float xb = 1;
		Vector3f n;
	};
	vector<float> outfile;
	outfile.push_back(0);
	float test[] = { 0.05,0.15,0.25,0.6,0.85,0.95 };

	MyMesh::Point pcc;
	int iith=UpOneInch(metb->ReturnIth(2),pcc);
	MyMesh::VertexHandle start = FindNearest(pcc);
	SurfaceCoe *metc = new SurfaceCoe(metb->ReturnCoe(),start,0,ios.mesh);
	metc->Init();
	metc->SetMIth(iith);

	float ju, max0, max1;
	max0 = abs(metb->DistSurface(meta->ReturnStartPoint())); //修改为以到平面b为基准点
	max1 = abs(metb->DistSurface(metc->ReturnStartPoint()));

	Vector3f coemet = metb->ReturnCoe();
	Vector3f axi = meta->ReturnCoe();
	float thert = acos(coemet.dot(axi)) / 2;
	Vector3f mix = coemet.cross(axi);
	mix.normalize();//这个很重要，一定要单位化

	float lin;
	bool ini = 0;
	vector<struct CutArry2> arry;
	Quaternionx out;
	int interv = (metb->ReturnIth(2) - meta->ReturnIth(2)) / WISTSECTION2; //把第一段分成5份
																		  //printf("interv:%d\n", interv);
	struct CutArry2 st;

	//int mma = WISTSECTION - 1;
	int mma = 0;
	if (metb->ReturnIth(2) - meta->ReturnIth(2) - WISTSECTION2*interv){
		st.a = mOutline2[(metb->ReturnIth(2) + meta->ReturnIth(2) - (WISTSECTION2 - 1)*interv) / 2].a;
		ju = abs(metb->DistSurface(st.a)) / max0;
		lin = thert*ju;
		
		st.xa = ju*ju;
		st.x = ju;
		st.xb = 1 - (1 - ju)*(1 - ju);
		st.x = xi*(1 - (ju + st.xa*(1 - ju) + st.xb*ju) / 2);
		//cout << st.xa << " : " << ju << " : " << st.xb << " : " << (ju + st.xa*(1 - ju) + st.xb*ju) / 2 << " : " << st.x << endl;
		
#ifdef SWITCHOPEN
		st.x = (1 - ju)*xi;
#endif // 

#ifdef SWITCHOPEN2
		st.x = test[mma]*xi;
#endif // 
		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
		arry.push_back(st);
		mma++;

		outfile.push_back(1-ju);
	}

	for (int i = WISTSECTION2 - 1; i >= 1; i--) {
		st.a = mOutline2[metb->ReturnIth(2) - i*interv].a;
		ju = abs(metb->DistSurface(st.a)) / max0;

		lin = thert*ju;

		st.xa = ju*ju;
		st.x = ju;
		st.xb = 1 - (1 - ju)*(1 - ju);
		st.x = xi*(1 - (ju + st.xa*(1 - ju) + st.xb*ju) / 2);
		cout << st.xa << " : " << ju << " : " << st.xb << " : " << (ju + st.xa*(1 - ju) + st.xb*ju) / 2 << " : " << st.x << endl;

#ifdef SWITCHOPEN
		st.x = (1 - ju)*xi;
#endif // 

#ifdef SWITCHOPEN2
		st.x = test[mma] * xi;
#endif // 
		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
		arry.push_back(st);
		mma++;

		outfile.push_back(1-ju);
	}
	outfile.push_back(1);

	interv = (metc->ReturnIth(2) - metb->ReturnIth(2)) / WISTSECTION2; //把第一段分成5份
	//printf("interv:%d\n", interv);
	//int mmb = WISTSECTION2 - 1;
	int mmb = 0;
	for (int i = 1; i < WISTSECTION2; i++) {
		st.a = mOutline2[metb->ReturnIth(2) + i*interv].a;
		ju = abs(metb->DistSurface(st.a)) / max1;
		lin = thert*ju;
		
		st.x = 1 - ju;
		st.xb = st.x*st.x;
		st.xa = 1 - (1 - st.x)*(1 - st.x);
		//cout << st.xa << " : " << st.x << " : " << st.xb << " : " << (st.x + st.xa*st.x + st.xb*(1 - st.x)) / 2 << endl;
		st.x = xi*(st.x + st.xa*st.x + st.xb*(1 - st.x)) / 2;

#ifdef SWITCHOPEN
		st.x = (1 - ju)*xi;
#endif // 

#ifdef SWITCHOPEN2
		st.x = (1-test[mmb]) * xi;
#endif // 

		st.n = metc->ReturnCoe();
		arry.push_back(st);
		mmb++;

		outfile.push_back(1-ju);
	}
	
	if (metc->ReturnIth(2) - metb->ReturnIth(2) - WISTSECTION2*interv) {
		st.a = mOutline2[(metc->ReturnIth(2) + metb->ReturnIth(2) + (WISTSECTION2 - 1)*interv) / 2].a;
		ju = abs(metb->DistSurface(st.a)) / max1;
		lin = thert*ju;
		
		st.x = 1 - ju;
		st.xb = st.x*st.x;
		st.xa = 1 - (1 - st.x)*(1 - st.x);
		cout << st.xa << " : " << st.x << " : " << st.xb << " : " << (st.x + st.xa*st.x + st.xb*(1 - st.x)) / 2 << endl;
		st.x = xi*(st.x + st.xa*st.x + st.xb*(1 - st.x)) / 2;

#ifdef SWITCHOPEN
		st.x = (1 - ju)*xi;
#endif // 
#ifdef SWITCHOPEN2
		st.x = (1-test[mmb]) * xi;
#endif // 

		st.n = metc->ReturnCoe();
		arry.push_back(st);
		mmb++;

		outfile.push_back(1-ju);
	}

#ifdef DEBUGBACK
	FILE *fp;
	fopen_s(&fp, "cutout2.txt", "w");
	for (auto i : outfile) {
		fprintf(fp, "%f\n", i);
	}
	fclose(fp);
#endif // DEBUG

	outfile.push_back(0);
	vector<float>gaussout;
	GaussianSmooth(outfile,gaussout,0.2);
	/*for (int i = 0; i < gaussout.size(); i++) {
	printf("back-gauss %d: %f \n", i, gaussout[i]);
	}*/

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

	MySurCutArry stt;
	stt.a = meta->ReturnVertexHandle();
	stt.x = 0;
	stt.n = meta->ReturnCoe();
	arryx.push_back(stt);
	int jj = 0;
	for (; jj < mma; jj++) {
		stt.a = mesh.vertex_handle((*(v_it_s + jj))->idx());
		stt.x = gaussout[jj + 1] * xi;//arry[jj].x;
		stt.n = arry[jj].n;
		arryx.push_back(stt);
	}
	stt.a = metb->ReturnVertexHandle();
	stt.x = gaussout[jj+1]*xi;
	stt.n = metb->ReturnCoe();
	arryx.push_back(stt);
	for (; jj < mma + mmb; jj++) {
		stt.a = mesh.vertex_handle((*(v_it_s + jj))->idx());
		stt.x = gaussout[jj + 2] * xi;//arry[jj].x;
		stt.n = arry[jj].n;
		arryx.push_back(stt);
	}
	stt.a = metc->ReturnVertexHandle();
	stt.x = 0;
	stt.n = metc->ReturnCoe();
	arryx.push_back(stt);

	free(minar);
	free(v_it_s);
	delete metc;
	return arryx;
}

vector<MySurCutArry> SurfaceCoe::OutCutOutline(float exp, SurfaceCoe *meta,SurfaceCoe *metb,SurfaceCoe *metc) { //腰围
	cout << "wist cut out..." << endl;
	vector<float>gauss;
	gauss.push_back(0);

	float xi = exp>0 ? -exp / metb->ReturnLength() : exp / metb->ReturnLength();
	struct CutArry2 {
		MyMesh::Point a;
		float x = 1;
		float xa = 1;
		float xb = 1;
		Vector3f n;
	};
	float ju, max0, max1;
	max0 = abs(metb->DistSurface(meta->ReturnStartPoint())); //修改为以到平面b为基准点
	max1 = abs(metb->DistSurface(metc->ReturnStartPoint()));

	Vector3f coemet = metb->ReturnCoe();
	Vector3f axi = meta->ReturnCoe();
	float thert = acos(coemet.dot(axi))/2;
	Vector3f mix = coemet.cross(axi);
	mix.normalize();//这个很重要，一定要单位化

	float lin;
	bool ini = 0;
	vector<struct CutArry2> arry;
	Quaternionx out;
	int interv = (metb->ReturnIth(2)-meta->ReturnIth(2))/ WISTSECTION; //把第一段分成5份
	//printf("interv:%d\n", interv);
	struct CutArry2 st;

	int mma = 0;
	if (metb->ReturnIth(2) - meta->ReturnIth(2) - WISTSECTION*interv) {
		st.a = mOutline2[(metb->ReturnIth(2) + meta->ReturnIth(2) - (WISTSECTION - 1)*interv) / 2].a;
		ju = abs(metb->DistSurface(st.a)) / max0;
		lin = thert*ju;
		//st.x = xi*(1 - ju);

		st.xa = ju*ju;
		st.x = ju;
		st.xb = 1 - (1 - ju)*(1 - ju);
		
		st.x = xi*(1-(ju + st.xa*(1 - ju) + st.xb*ju)/2);

		gauss.push_back(1 - ju);
		//cout << st.xa << " : " << ju << " : " << st.xb << " : " << (ju + st.xa*(1 - ju) + st.xb*ju) / 2 <<" : "<<st.x<< endl;
		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
		arry.push_back(st);
		mma++;
	}

	for (int i = WISTSECTION-1; i >=1; i--) {
		st.a = mOutline2[metb->ReturnIth(2)-i*interv].a;
		ju = abs(metb->DistSurface(st.a))/max0;

		lin = thert*ju;
		//st.x = xi*(1-ju);

		st.xa = ju*ju;
		st.x = ju;
		st.xb = 1 - (1 - ju)*(1 - ju);
		st.x = xi*(1-(ju + st.xa*(1-ju) + st.xb*ju)/2);
		//cout << st.xa << " : " << ju << " : " << st.xb << " : " << (ju + st.xa*(1 - ju) + st.xb*ju) / 2 << " : " << st.x << endl;
		gauss.push_back(1 - ju);

		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
		arry.push_back(st);
		mma++;
	}

	gauss.push_back(1);

	interv = (metc->ReturnIth(2) - metb->ReturnIth(2)) / WISTSECTION2; //把第二段分成5份
	//printf("interv:%d\n", interv);
	coemet = metb->ReturnCoe();
	axi = metc->ReturnCoe();
	thert = acos(coemet.dot(axi))/2;
	//mix = coemet.cross(axi);
	//mix.normalize();//其实等于-1就行
	mix = Vector3f(0, -1, 0);

	int mmb = 0;
	for (int i = 1; i < WISTSECTION2; i++) {
		st.a = mOutline2[metb->ReturnIth(2)+i*interv].a;
		ju = abs(metb->DistSurface(st.a))/max1;
		lin = thert*ju;
		//cout << ju << " : ";
		//ju = 1 - ju*ju;
		//st.x = xi*(1-ju);// *ju;
		//cout << ju << " : ";
		//cout << 1-ju << endl;
		st.x = 1 - ju;
		st.xb = st.x*st.x;
		st.xa = 1 - (1 - st.x)*(1 - st.x);
		//cout << st.xa << " : " << st.x << " : " << st.xb << " : " << (st.x + st.xa*st.x + st.xb*(1 - st.x)) / 2 << endl;
		st.x = xi*(st.x + st.xa*st.x + st.xb*(1 - st.x)) / 2;
		gauss.push_back(1 - ju);

		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
		arry.push_back(st);
		mmb++;
	}
	
	if ( metc->ReturnIth(2) - metb->ReturnIth(2) - WISTSECTION2*interv) {
		st.a = mOutline2[(metc->ReturnIth(2) + metb->ReturnIth(2) + (WISTSECTION2 - 1)*interv) / 2].a;
		ju = abs(metb->DistSurface(st.a)) / max1;
		lin = thert*ju;
		//cout << ju << " : ";
		//ju = 1 - ju*ju;
		//st.x = xi*(1-ju);// *ju;
		//cout << ju << " : ";
		//cout << 1-ju << endl;
		st.x = 1 - ju;
		st.xb = st.x*st.x;
		st.xa = 1 - (1-st.x)*(1-st.x);
		cout << st.xa << " : " << st.x << " : " << st.xb << " : " << (st.x + st.xa*st.x + st.xb*(1-st.x))/2 << endl;
		st.x = xi*(st.x + st.xa*st.x + st.xb*(1 - st.x)) /2;
		gauss.push_back(1 - ju);

		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
		arry.push_back(st);
		mmb++;
	}

	gauss.push_back(0); //just for smooth
	
	vector<float>gaussout;
	GaussianSmooth(gauss, gaussout, 0.2);
	/*for (int i = 0; i < gaussout.size(); i++) {
		printf("gauss %d: %f \n", i, gaussout[i]);
	}*/
	//gaussout = gauss;

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

	MySurCutArry stt;
	stt.a = meta->ReturnVertexHandle();
	stt.x = 0;//gaussout[0] * xi; 
	stt.n = meta->ReturnCoe();
	arryx.push_back(stt);
	int jj = 0;
	for (; jj < mma; jj++) {
		stt.a = mesh.vertex_handle((*(v_it_s + jj))->idx());
		stt.x = xi*gaussout[jj + 1];//arry[jj].x;
#ifdef STTWISTX
		stt.x = arry[jj].x;
#endif // 
		stt.n = arry[jj].n;
		arryx.push_back(stt);
	}
	stt.a = metb->ReturnVertexHandle();
	stt.x = gaussout[jj+1]*xi;
#ifdef STTWISTX
	stt.x = xi;
#endif // 
	stt.n = metb->ReturnCoe();
	arryx.push_back(stt);
	for (; jj < mma +mmb;jj++) {
		stt.a = mesh.vertex_handle((*(v_it_s + jj))->idx());
		stt.x = xi*gaussout[jj + 2]; //arry[jj].x;
#ifdef STTWISTX
		stt.x = arry[jj].x;
#endif // 

		stt.n = arry[jj].n;
		arryx.push_back(stt);
	}
	stt.a = metc->ReturnVertexHandle();
	stt.x = 0;// gaussout[gaussout.size() - 1] * xi; ///0;
	stt.n = metc->ReturnCoe();
	arryx.push_back(stt);

	free(minar);
	free(v_it_s);

	return arryx;
}

Vector3f SurfaceCoe::TempVector() {
	//return Vector3f(mVertexStart[0]-38.5335, mVertexStart[1]+0.4859, mVertexStart[2]-157.5532).normalized();
	return Vector3f(mVertexStart[0] - mVertexEnd[0], mVertexStart[1]- mVertexEnd[1], mVertexStart[2] -mVertexEnd[2]).normalized();
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
