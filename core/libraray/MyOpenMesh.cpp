#include "stdafx.h"
#include "MyOpenMesh.h"

struct CutArry2
{
	MyMesh::Point a;
	float x = 1;
	Vector3f n;
};

MyOpenMesh::MyOpenMesh(float a)
{
}

void MyOpenMesh::ReadStlfile(const char * argg) 
{
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

void MyOpenMesh::WriteStlfile(const char *argg,int i){
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
	//cout << "Finding Tri Nearest Vertex..." << endl;
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

void MyOpenMesh::ShoeSpin(Quaternionx &af, Vector3f sf) {
	MyMesh::Point p;
	Quaternionx afi = af.inverse();
	Quaternionx out;
	MyMesh::Point mshift = MyMesh::Point(sf[0], sf[1], sf[2]);
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p = mesh.point(*v_it);
		Quaternionx ss(0, p[0], p[1], p[2]);
		out = af*ss*afi;
		p = MyMesh::Point(out.x(), out.y(), out.z())+mshift;
		mesh.set_point(*v_it,p);
	}
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

void MyOpenMesh::ShoeExpansion(vector<SurfaceCoe*> &arrx) { //这个需要考虑变的方向

	vector<SurfaceCoe*>&arr=arrx;
	cout << "正在对掌围进行变形..." << endl;
	MyMesh::Point p, p1, p2;// n为递增向量
	float s1, s2;
	MyMesh::Normal nf;

	/*int k = mesh.n_vertices();
	int oc = 0;*/
	int arrlen = arr.size();

	int sta = 0, end = 0;

	vector<MyMesh::VertexHandle>select; //选取最后不连续点进行滤波；

	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
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
						/*nf = mesh.normal(*v_it);
						nf.normalized();
						p += nf*(arr[i + 1]->FindNearestPointF(p, s1));*/
						select.push_back(mesh.vertex_handle(v_it->idx()));
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
	TailGaussionFilter(select, 6);//滤波五次
	TotalGaussionFilter(2);
}

void MyOpenMesh::TotalGaussionFilter(int count) { //全局滤波
	float sigma = 1;
	int inc = count;
	MyMesh::Point p;
	MyMesh::VertexHandle h;
	while (inc--) {
		for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			p = mesh.point(*v_it);
			h = mesh.vertex_handle(v_it->idx());
			vector<MyMesh::Point> vtemp;
			vtemp.push_back(p);
			for (MyMesh::VertexVertexIter vv_it = mesh.vv_begin(h); vv_it.is_valid(); ++vv_it)
			{
				vtemp.push_back(mesh.point(mesh.vertex_handle(vv_it->idx())));
			}
			p = GaussionArroundVertex(vtemp, sigma);
			mesh.set_point(h, p);
		}
	}
}

void MyOpenMesh::TailGaussionFilter(vector<MyMesh::VertexHandle>& vm,int count)
{
	MyMesh::Point p;
	float sigma = 1;
	int num = count;
	while (num--) {
		for (auto i : vm) {
			p = mesh.point(i);
			vector<MyMesh::Point> vtemp;
			vtemp.push_back(p);
			for (MyMesh::VertexVertexIter vv_it = mesh.vv_begin(i); vv_it.is_valid(); ++vv_it)
			{
				vtemp.push_back(mesh.point(mesh.vertex_handle(vv_it->idx())));
			}
			p = GaussionArroundVertex(vtemp, sigma);
			mesh.set_point(i, p);
		}
	}
}

MyMesh::Point MyOpenMesh::GaussionArroundVertex(vector<MyMesh::Point>& pGray,float sigma) {
	int sz = pGray.size();
	float dWeightSum = 0;//滤波系数总和  

	MyMesh::Point dDotMul(0, 0, 0);//高斯系数与图像数据的点乘  
	float dDis, dValue;

	for (int i = 0; i<sz; i++)
	{
		dDis = (pGray[0] - pGray[i]).norm();//距离
		dValue = exp(-(1 / 2)*dDis*dDis / (sigma*sigma)) / (sqrt(2 * 3.1415926)*sigma);

		dDotMul += dValue*pGray[i];//叠加的点坐标
		dWeightSum += dValue;
	}
	return dDotMul / dWeightSum;
}

void MyOpenMesh::ShoeExpansionWist(vector<SurfaceCoe*> &arr) {
	//cout << "Now is shoe wist Expansing..." << endl;
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

void MyOpenMesh::ShoeExpansionWist(SurfaceCoe*meta, SurfaceCoe*metb, SurfaceCoe*metc){
	//cout << "Now is shoe wist2 Expansing..." << endl;
	MyMesh::Point p,p1,p2; float lin,s1,s2;
	//int cis = 4;//迭代次数4次
	float sigma = 1;
	SurfaceCoe *smc;
	map<int,MyMesh::Point>select;
	float ran = 6; //mm 扩散缓冲区域
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p = mesh.point(*v_it);
		lin=meta->DistSurface(p);
		if (lin >=0) {
			continue;
		}
		if ((lin <0) && (lin >(0-ran))) {
			select[v_it->idx()] = MyMesh::Point(0,0,0);
			continue;
		}
		lin = metc->DistSurface(p);
		if (lin <= -ran) {
			continue;
		}
		if ((lin > -ran) && (lin <0)) {
			select[v_it->idx()] = MyMesh::Point(0, 0, 0);
			continue;
		}
		smc = metb->DistSurface(p) > 0 ? meta : metc;
		p1 = metb->FindNearestPoint(p, s1);
		p2 = smc->FindNearestPoint(p, s2);
		select[v_it->idx()] = p1*(s2 / (s1 + s2));// + p2*(s1 / (s1 + s2));
	}
	map<int, MyMesh::Point>::iterator it(select.begin()),it_s;
	for (; it != select.end(); it++) {
		set<int> cvm;
		BotIteration(cvm, it->first, 0);
		vector<MyOutBottom> smv;
		for (auto i : cvm) {
			it_s=select.find(i);
			MyOutBottom git;
			if (it_s == select.end()) {
				git.a = mesh.point(mesh.vertex_handle(i));
				git.s = MyMesh::Point(0,0,0);
				smv.push_back(git);
				continue;
			}
			git.a= mesh.point(mesh.vertex_handle(i));
			git.s = select[i];
			smv.push_back(git);
		}
		it->second = GaussFilter(smv, sigma);
	}
	for (it = select.begin(); it != select.end(); it++) {
		p = mesh.point(mesh.vertex_handle(it->first))+it->second;
		mesh.set_point(mesh.vertex_handle(it->first),p);
	}
}

void MyOpenMesh::ShoeExpansionWist2(SurfaceCoe*meta, SurfaceCoe*metb, SurfaceCoe*metc) {
	//cout << "Now is shoe wist2 Expansing..." << endl;
	MyMesh::Point p, p1, p2; float lin, s1, s2;
	//int cis = 4;//迭代次数4次
	float sigma = 1;
	SurfaceCoe *smc;
	map<int, MyMesh::Point>select;
	float ran = 6; //mm 扩散缓冲区域
	struct scc {
		int i;//main
		vector<int> veh;
	};
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		p = mesh.point(*v_it);
		lin = meta->DistSurface(p);
		if (lin >= 0) {
			continue;
		}
		if ((lin <0) && (lin >(0 - ran))) {
			select[v_it->idx()] = MyMesh::Point(0, 0, 0);
			continue;
		}
		lin = metc->DistSurface(p);
		if (lin <= -ran) {
			continue;
		}
		if ((lin > -ran) && (lin <0)) {
			select[v_it->idx()] = MyMesh::Point(0, 0, 0);
			continue;
		}
		smc = metb->DistSurface(p) > 0 ? meta : metc;
		p1 = metb->FindNearestPoint(p, s1);
		p2 = smc->FindNearestPoint(p, s2);
		select[v_it->idx()] = p1*(s2 / (s1 + s2));// + p2*(s1 / (s1 + s2));
	}

	map<int, MyMesh::Point>::iterator it, it_s;
	for (int i = 0; i < 5; i++) { //迭代5次滤波
		it = select.begin();
		for (; it != select.end(); it++) {
			//set<int> cvm;
			//BotIteration(cvm, it->first, 0);
			vector<MyOutBottom> smv;
			for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(mesh.vertex_handle(it->first)); vv_it.is_valid(); ++vv_it)
			{
				it_s = select.find(vv_it->idx());
				MyOutBottom git;
				if (it_s == select.end()) {
					git.a = mesh.point(mesh.vertex_handle(vv_it->idx()));
					git.s = MyMesh::Point(0, 0, 0);
					smv.push_back(git);
					continue;
				}
				git.a = mesh.point(mesh.vertex_handle(vv_it->idx()));
				git.s = select[vv_it->idx()];
				smv.push_back(git);
			}
			it->second = GaussFilter(smv, sigma);
		}
	}
	for (it = select.begin(); it != select.end(); it++) {
		p = mesh.point(mesh.vertex_handle(it->first)) + it->second;
		mesh.set_point(mesh.vertex_handle(it->first), p);
	}
}

MyMesh::Point MyOpenMesh::GaussFilter(vector<MyOutBottom>&pGray,float sigma) {
	int sz = pGray.size();
	float dWeightSum=0;//滤波系数总和  

	MyMesh::Normal dDotMul(0, 0, 0);//高斯系数与图像数据的点乘  
	float dDis, dValue;
	for (int i = 0; i<sz; i++)
	{
		dDis = (pGray[0].a - pGray[i].a).norm();//距离
		dValue = exp(-(1 / 2)*dDis*dDis / (sigma*sigma)) / (sqrt(2 * 3.1415926)*sigma);

		dDotMul += dValue*pGray[i].s;//叠加的点坐标
		dWeightSum += dValue;
	}
	return dDotMul / dWeightSum;
}

void MyOpenMesh::BotIteration(set<int>&arr, int idx, int iver) {
	if (iver > ITERATIONCISHUPOINT) {
		return;
	}
	iver++;
	MyMesh::VertexVertexIter vv_it = mesh.vv_iter(mesh.vertex_handle(idx));
	for (; vv_it.is_valid(); ++vv_it)
	{
		BotIteration(arr, vv_it->idx(), iver);

		arr.insert(vv_it->idx());
	}
	return;
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
	mCoe.data()[0] = mCoeABC[0];
	mCoe.data()[1] = mCoeABC[1];
	mCoe.data()[2] = mCoeABC[2];
	mCoe.data()[3] = d;
	//mCoe= Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);
}
SurfaceCoe::SurfaceCoe(MyMesh::Point mid, MyMesh::Point end, MyMesh::VertexHandle cf, MyMesh &d) :
	mesh(d),
	mHandleBegin(cf),
	mVertexMid(mid),
	mVertexEnd(end)
{
	mVertexStart = mesh.point(cf);

	Vector3f a = MyOpenMesh::EigenTransfer(mVertexStart);
	Vector3f b = MyOpenMesh::EigenTransfer(mVertexMid);
	Vector3f c = MyOpenMesh::EigenTransfer(mVertexEnd);
	Vector3f ab = a - b;
	ab = (a - c).cross(ab);
	//mCoeABC = Vector3f(ab[0], ab[1], ab[2]);

	mCoeABC.data()[0] = ab[0];
	mCoeABC.data()[1] = ab[1];
	mCoeABC.data()[2] = ab[2];

	float df = ab.dot(Vector3f(0, 0, 0) - a) / mCoeABC.norm();
	mCoeABC.normalize();
	//mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);

	mCoe.data()[0] = mCoeABC[0];
	mCoe.data()[1] = mCoeABC[1];
	mCoe.data()[2] = mCoeABC[2];
	mCoe.data()[3] = df;
}


bool SurfaceCoe::Init() {//(MyMesh::VertexHandle *vertex,MyMesh &mmesh){
	MyOutNormal abc;
	abc.a = mVertexStart;
	abc.n = mesh.normal(mHandleBegin);
	//abc.nf = abc.n;//for debug;
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

	mCoeABC.data()[0] = ab[0];
	mCoeABC.data()[1] = ab[1];
	mCoeABC.data()[2] = ab[2];
	
	float d = ab.dot(Vector3f(0, 0, 0) - a)/ mCoeABC.norm();
	mCoeABC.normalize();
	//mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);

	mCoe.data()[0] = mCoeABC[0];
	mCoe.data()[1] = mCoeABC[1];
	mCoe.data()[2] = mCoeABC[2];
	mCoe.data()[3] = d;
	//float d = ab.dot(Vector3f(0, 0, 0) - a);
	//mCoe = Vector4f(mCoeABC[0], mCoeABC[1], mCoeABC[2], d);

	MyOutNormal abc;
	abc.a = mVertexStart;
	abc.n = mesh.normal(mHandleBegin);
	//abc.nf = abc.n;//for debug;
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
	if (!ini) {
		for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(mHandleBegin); voh_it.is_valid(); ++voh_it) {
			heh = mesh.halfedge_handle(voh_it->idx());
			Vector3f jug = NextHalfEdgeJudge2(mesh.next_halfedge_handle(heh));
			if (jug != Vector3f(9999, 9999, 9999)) {
				if (mVertexStart[1]< jug[1]) {
					heh = mesh.next_halfedge_handle(heh);
					IterationHalfEdge(heh);
					ini = 1;
				}
			}
		}
	}
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
	mLength = 0;
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
	//mc.nf = mc.n;
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
	
	mExtension = -tar;
	return OutlineExpansion();
}
//
//float SurfaceCoe::TopSlide(SurfacePure*met) {
//	if (abs(met->DistSurface(mVertexStart)) > TOPOFFSET) {
//		return 0;
//	}
//
//	vector<float>input, output;
//	for (int i = mIth[1]; i < mOutline2.size(); i++) {
//		if (abs(met->DistSurface(mOutline2[i].a)) <= TOPOFFSET) {
//			mOutline2[i].x = abs(met->DistSurface(mOutline2[i].a))/ TOPOFFSET;
//		}
//		input.push_back(mOutline2[i].x);
//	}
//	for (int i = 0; i <= mIth[0]; i++) {
//		if (abs(met->DistSurface(mOutline2[i].a)) <= TOPOFFSET) {
//			mOutline2[i].x = abs(met->DistSurface(mOutline2[i].a))/ TOPOFFSET;
//		}
//		input.push_back(mOutline2[i].x);
//	}
//	GaussianSmooth(input, output,3);
//
//	for (int i = 0; i < mOutline2.size() - mIth[1]; i++) {
//		mOutline2[i + mIth[1]].x = output[i];
//	}
//	for (int i = 0; i <= mIth[0]; i++) {
//		mOutline2[i].x = output[i + mOutline2.size() - mIth[1]];
//	}
//	return 1;
//}
float SurfaceCoe::TopSlide(SurfacePure*met) {
	if (abs(met->DistSurface(mVertexStart)) > 1) { //判断是否为碗口顶板线，如果是进一步做处理
		return 0;
	}
	float window = GAUSSIONFILTERNUM + 1;
	int sz = mOutline2.size();
	//for (int i = 0; i < mIth[0]*2/3; i++) {// 1/2 or 1/3
	//	mOutline2[i].x = 1;//mLen[0]
	//}
	//for (int i = sz - 1; i > (mIth[1]*2+sz)/3; i--) {
	//	mOutline2[i].x = 1;  //mLen[2]
	//}
	for (int i = 0; i < mIth[0]*1/2; i++) {// 1/2 or 1/3
		mOutline2[i].x = 1;//mLen[0]
	}
	for (int i = sz - 1; i > (mIth[1]+sz)/2; i--) {
		mOutline2[i].x = 1;  //mLen[2]
	}

	//滤波两次，增加平滑效果
	vector<float>input, output;
	for (int i = mIth[1]; i < mOutline2.size(); i++) {
		input.push_back(mOutline2[i].x);
	}
	for (int i = 0; i <= mIth[0]; i++) {
		input.push_back(mOutline2[i].x);
	}

	GaussianSmooth(input, output, window);
	input.clear();
	input = output;
	output.clear();
	GaussianSmooth(input, output, window);
	input.clear();
	input = output;
	output.clear();
	GaussianSmooth(input, output, window);
	input.clear();
	input = output;
	output.clear();
	GaussianSmooth(input, output, window);
	input.clear();
	input = output;
	output.clear();
	GaussianSmooth(input, output, window);
	input.clear();
	input = output;
	output.clear();
	GaussianSmooth(input, output, window);
	input.clear();
	input = output;
	output.clear();
	GaussianSmooth(input, output, window);

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

	TopSlide(met);

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
	//mExtension = mLength*mX; //由扩散系数转换为扩散距离

	//TopSlide(met);
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
	float sout = mExtension * 10;
	float s = 0, li = sout/2, pp = 0;
	float aa = 0, bb = 1, cc = 0.5;
	while (abs(abs(s - mLength) - abs(mExtension)) > ADDDIFFERENCE) {
		bso = mOutline2;
		for (int j = 0; j < bso.size(); j++) {
			bso[j].a += bso[j].n*bso[j].x*li;
		}
		s = TotalLengh(bso);
		pp = -mExtension/(s - mLength);
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
	MyMesh::Point ac=mOutline2[0].m;
	float lii = 0;
	for (auto i : mOutline2) {
		lii+=(ac - i.m).norm();
		ac = i.m;
	}
	mLength = lii;
	mExtensionli = li;
	return li;
}
//class CircleIndex {
//public:
//	CircleIndex(MyMesh::Normal a,int b) :
//	mN(a),
//	Ith(b)
//	{}
//	MyMesh::Normal RtNormal() { return mN; }
//	int RtIth() { return Ith; }
//
//	CircleIndex *next;
//	CircleIndex *up;
//private:
//	MyMesh::Normal mN;
//	int Ith;
//};

void SurfaceCoe::InitMidEndPoint(vector<MyMesh::Point>&fw) {
	int window = DIFWINDOW;
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
	//cout << "topv: " << topv << " " << lin << endl;
	//cout << upstep << endl;
	//set<MyBotOutLine> sum1,sum2;//for debug
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
		/*MyBotOutLine git;
		git.s = dif;
		git.x = dif.norm();
		git.ith = i;
		sum1.push_back(git);*/
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
		/*MyBotOutLine git;
		git.s = dif;
		git.x = dif.norm();
		git.ith = i;
		sum2.push_back(git);*/
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

int SurfaceCoe::InitMidTopPoint(vector<MyMesh::Point>&fw) {
	int window = DIFWINDOW;
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

	//set<MyBotOutLine> sum1,sum2;//for debug;
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
		sm.push_back(tf[i].s);
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

		if (lin < dif.norm()) {
			lin = dif.norm();
			ith = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return 0;
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
		sm.push_back(tf[i].s);
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

		if (lin < dif.norm()) {
			lin = dif.norm();
			ith = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return 0;
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
	return 1;
}

void SurfaceCoe::InitMidTopPoint(vector<MyMesh::Point>&fw, float x) {
	int window = DIFWINDOW;
	float upstep = 7;//默认

	float lin = mOutline2[0].a[2];
	float topv = mOutline2[0].a[2];
	for (int i = 1; i < mOutline2.size(); i++) {
		if (lin < mOutline2[i].a[2]) {
			lin = mOutline2[i].a[2];
		}
		if (topv > mOutline2[i].a[2]) {
			topv = mOutline2[i].a[2];
		}
	}//默认情况下bot应该处于500左右范围内的值(该值远远大于window窗口滤波范围)
	upstep = abs((topv - lin) / 8);

	vector<MyBotOutLine>tf;
	for (int i = mOutline2.size() * 2 / 3; i < mOutline2.size(); i++) {
		MyBotOutLine gt;
		gt.ith = i;
		gt.s = mOutline2[i].n;
		gt.a = mOutline2[i].a;
		tf.push_back(gt);
	}
	for (int i = 0; i < mOutline2.size() / 3; i++) {
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

	lin = 0; int ith[2] = {0};
	for (int i = bot; i < tf.size(); i++) {
		if (abs(tf[i].a[2] - tf[bot].a[2]) > upstep) {
			break;
		}
		vector<MyMesh::Normal>sm;
		sm.push_back(tf[i].s);
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

		if (lin < dif.norm()) {
			lin = dif.norm();
			ith[1] = i;
		}
	}
	if (!ith) {
		cout << "ith is zero!" << endl;
		return;
	}

	lin = 0; //ith[] = 0;
	for (int i = bot; i >= 0; i--) {
		if (abs(tf[i].a[2] - tf[bot].a[2]) > upstep) {
			break;
		}
		vector<MyMesh::Normal>sm;
		sm.push_back(tf[i].s);
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

		if (lin < dif.norm()) {
			lin = dif.norm();
			ith[0] = i;
		}
	}
	/*cout << ith[0]<<" : "<< ith[1] << endl;
	cout << tf[ith[0]].a << endl;
	cout << tf[ith[1]].a << endl;
	cout << ith[0] + (ith[1] - ith[0])*x << endl;
	cout << tf[ith[0] + (ith[1] - ith[0])*x].a << endl;*/
	fw.push_back(tf[ith[0]+(ith[1] - ith[0])*x].a);
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
bool SurfaceCoe::OutlineEigenM(vector<Vector3f> *a) {
	if (!mOutline2.size()) {
		cout << "outline no points" << endl;
		return false;
	}
	for (auto i : mOutline2) {
		a->push_back(Vector3f(i.m[0], i.m[1], i.m[2]));
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
	SurfaceCoe *sfc=new SurfaceCoe(start, mesh.point(end), mesh);
	//cout << "Now is finding Metara..." << endl;
	float dd = MAXIMUMX, ss;
	int metara;
	for (int i = mIth[0] * 0.13; i < mIth[0]*0.85 ; i++) {
		sfc->SetMidPoint(mOutline2[i].a);
		if (sfc->Init(0)) {
			ss = sfc->ReturnLength();
			if (ss < dd) {
				dd = ss;
				metara = i;   //mLen[2]=i;
			}
		}
	}
	MyMesh::VertexHandle sfcstartp=sfc->PointExchange(mOutline2[metara].a);
	//cout << "Find Metara ith:" << metara << endl;
	vector<MyMesh::Point> fm;
	sfc->InitMidEndPoint(fm);
	if (fm.size() < 2) {
		cout << "metara point error" << endl;
		return NULL;
	}
	//sfc->GiveMidEndPoint(fm[0],fm[1]);

	//MyMesh::VertexHandle sfcstartp = FindNearest(mOutline2[metara].a);
	SurfaceCoe *ret = new SurfaceCoe(fm[0], fm[1], sfcstartp, mesh);
	//SurfaceCoe *ret = new SurfaceCoe(MyMesh::Point(160.563766,34.669678,2.144665),MyMesh::Point(146.947586,-44.111824,3.824428),sfcstartp,mesh);//手动给出底边的两个点
	ret->Init();
	ret->CoquerMidEnd();
	ret->SetMIth(metara);
	//sfc->SetMIth(metara);
	delete sfc;
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
	MyMesh::Point mdp, enp, of;
	int iith = UpOneInch(met->ReturnIth(2),of);
	MyMesh::VertexHandle start = FindNearest(of);

	float min = MAXIMUMX; int mid = 0;
	for (int i = mIth[1]+10; i < ccmid ; i++) {
		mdp = enp = mOutline2[i].a;
		enp[1] += 1;
		SurfaceCoe sfc(mdp, enp, start ,mesh);
		if (sfc.Init()) {
			sfc.CalculateLen();
			lin = sfc.ReturnLength();
			if (lin <min ) {
				min = lin;
				mid = i;   //mLen[2]=i;
			}
		}
	}
	//cout << "Waist Ith:" << mid << endl;
	if (!mid) {
		cout << "mid error!" << endl;
	}
	mdp = enp = mOutline2[mid].a;
	enp[1] += 1;
	SurfaceCoe* ret= new SurfaceCoe(enp,mdp, start, mesh);
	ret->SetMIth(iith);
	ret->Init();
	ret->InitMidEndPoint();
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
			af= mOutline2[i-1].a;
			return i-1;
		}
	}
	return -1;
}

SurfaceCoe* SurfaceCoe::SfcMoveXLen(SurfaceCoe *toe,float x) 
{
	if ((!x) || (x<0)) {
		cout << "Toe x is Zero or minus!" << endl;
		return NULL;
	}
	MyMesh::Point pstart = mOutline2[mIth[1]].a;
	float tdist = 0; int  iith;
	for (int i = toe->ReturnIth(2); i < mOutline2.size(); i++) {
		tdist += (pstart - mOutline2[i].a).norm();
		if (tdist > x) {
			iith = i - 1; //定出跖围在鞋楦底部的起始点
			break;
		}
	}
	MyMesh::VertexHandle start = FindNearest(mOutline2[iith].a);
	SurfaceCoe *metc = new SurfaceCoe(toe->ReturnCoe(), start, 0, mesh);
	metc->Init();
	metc->SetMIth(iith);
	return metc;
}

SurfaceCoe* SurfaceCoe::FindToeBottomPoint(float dis)
{
	if ((!dis)||(dis<0)) {
		cout << "Toe dis is Zero!" << endl;
		return NULL;
	}
	MyMesh::Point pstart = mOutline2[mIth[1]].a;
	float tdist = 0; int  iith;
	for (int i = mIth[1]+1; i < mOutline2.size();i++) {
		tdist += (pstart - mOutline2[i].a).norm();
		if (tdist > dis) {
			iith = i - 1; //定出跖围在鞋楦底部的起始点
			break;
		}
	}

	MyMesh::Point mdp, enp;
	MyMesh::VertexHandle start = FindNearest(mOutline2[iith].a);

	float min = MAXIMUMX; int mid = 0;
	float templen = 0;
	for (int i = 10; i < mIth[0]-50; i++) {
		mdp = enp = mOutline2[i].a;
		enp[1] += 1;
		SurfaceCoe sfc(mdp, enp, start, mesh);
		if (sfc.Init()) {
			sfc.CalculateLen();
			templen = sfc.ReturnLength();
			if (templen <min) {
				min = templen;
				mid = i;   //这个是定出底部，在线段上半部分找一个最小围度点；
			}
		}
	}
	if (!mid) {
		cout << "mid error!" << endl;
	}
	start= FindNearest(mOutline2[mid].a);
	mdp = enp = mOutline2[iith].a;
	enp[1] += 1;
	SurfaceCoe* ret = new SurfaceCoe(enp, mdp, start, mesh);
	ret->SetMIth(iith);
	ret->Init();
	ret->InitMidEndPoint();
	return ret;
}


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

int SurfaceCoe::FindNearestOutline(MyMesh::Point a) {
	float lin = MAXIMUMX,ch;
	int ith = 0;
	for (int i = 0; i < mOutline2.size(); i++) {
		ch = (a-mOutline2[i].a).norm();
		if (lin > ch) {
			lin = ch;
			ith = i;
		}
	}
	return ith;
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
	printf("inverv: %d\n", interv);
	for (int i = 1; i < CUTSECTION3; i++) 
	{
		struct CutArry2 st;
		st.a = mOutline2[i*interv].a;
		cout <<i<<" : "<< i*interv << endl;
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

vector<MySurCutArry> SurfaceCoe::OutCutOutline(float exp, SurfaceCoe *met,float heelh) //改用距离做分量
{ 
	float xi = exp>0 ? -exp / met->ReturnLength() : exp / met->ReturnLength();
	//中轴横截点
	float szz = (mVertexMid - mVertexEnd).norm();
	float ssst = heelh*0.4 + szz*0.25;
	//printf("szzz:%f  ssst:%f\n", szz, ssst);
	MyMesh::Point axip;
	if (ssst >= szz) 
	{
		axip = mVertexMid;
	}
	else 
	{
		axip = (mVertexMid - mVertexEnd)*(heelh*(0.4)/szz + 0.25) + mVertexEnd;
	}
	Vector3f axi(mVertexStart[0] - axip[0], mVertexStart[1] - axip[1], mVertexStart[2] - axip[2]);
	axi.normalize();

	float ju, max0, max1;
	max0 = met->DistSurface(mVertexStart);
	max1 = abs(met->DistSurface(axip)); //vertexMid or vertexend

	Vector3f coemet = met->ReturnCoe();
	//float thert=acos(coemet.dot(axi));//原版的，忘记除以2了
	float thert = acos(coemet.dot(axi)) / 2;
	Vector3f mix = coemet.cross(axi);
	mix.normalize();

	float lin;
	bool ini = 0;
	float slen = mLen[0];

	if (!slen) {
		cout << "cut out len error!" << endl;
	}
	MyMesh::Point ack = mOutline2[0].a;
	float interv = slen / CUTSECTION3;//原则上分成32份，实际上取32个点，最后一个预留出来余量；
	//printf("slen: %f interv: %f\n", slen,interv);
	
	vector<struct CutArry2> arry;
	Quaternionx out;
	float linslen = 0; int ith = 1;
	for (int i = 1; i < CUTSECTION3; i++) {
		while (linslen < interv) {
			linslen += (mOutline2[ith].a-ack).norm();
			ack = mOutline2[ith].a;
			ith++;
		}
		ith--;
		ack = mOutline2[ith - 1].a;
		if (ith >= (mIth[0]-10)) 
		{
			break;
		}
	//	printf("i:%d ith:%d lislen:%f\n",i, ith, linslen);
		linslen = 0;
		struct CutArry2 st;
		st.a = mOutline2[ith].a;
		ju = met->DistSurface(st.a);
		if (ju >0) {
			st.x = xi*(1 - ju / max0);
			//st.x = xi;
			lin = thert*(ju / max0);
		}
		else if (ju<0) {
			if (!ini) {
				ini = true;
				struct CutArry2 sts;
				sts.x = xi;
				sts.n = met->ReturnCoe();
				sts.a = met->ReturnStartPoint();
				arry.push_back(sts);
			}
			st.x = xi;
			//st.x = xi*(1 - abs(ju) / max1);
			lin = thert*(abs(ju) / max1);
		}
		else {
			continue;
		}
		Quaternionx mtransfer(cos(lin), sin(lin)*mix[0], sin(lin)*mix[1], sin(lin)*mix[2]);
		Quaternionx scoemet(0, coemet[0], coemet[1], coemet[2]);
		out = mtransfer*scoemet*(mtransfer.inverse());
		st.n = Vector3f(out.x(), out.y(), out.z());
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

float SurfaceCoe::FindNearestPointF(MyMesh::Point a, float &s) {
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
		return 0;
	}
	MyMesh::Point fb, fc, fd, j1, j2, p;
	fb = mOutline2[k].a;
	float t;
	for (int i = 1; i < 4; i++) {
		m = k - i;
		m = m < 0 ? mOutline2.size() + m : m;

		fc = mOutline2[m].a;
		fd = fb - fc;
		t = ((a[0] - fb[0])*fd[0] + (a[1] - fb[1])*fd[1] + (a[2] - fb[2])*fd[2]) / (fd[0] * fd[0] + fd[1] * fd[1] + fd[2] * fd[2]); //两直线的交点
		p = t*fd + fb;
		j1 = fb - p;
		j2 = fc - p;
		if ((j1[0] * j2[0] + j1[1] * j2[1] + j1[2] * j2[2]) < 0) {
			float s1 = j1.norm(), s2 = j2.norm();
			return (mOutline2[k].x*s2 / (s1 + s2) + mOutline2[m].x*s1 / (s1 + s2))*mExtensionli;
		}

		n = k + i;
		n = n >= mOutline2.size() ? n - mOutline2.size() : n;

		fc = mOutline2[n].a;
		fd = fb - fc;
		t = ((a[0] - fb[0])*fd[0] + (a[1] - fb[1])*fd[1] + (a[2] - fb[2])*fd[2]) / (fd[0] * fd[0] + fd[1] * fd[1] + fd[2] * fd[2]); //Vector3f s = t*d + a; 	t = (c-a).dot(d) / (d.dot(d));
		p = t*fd + fb;
		j1 = fb - p;
		j2 = fc - p;
		if ((j1[0] * j2[0] + j1[1] * j2[1] + j1[2] * j2[2]) < 0) {
			float s1 = j1.norm(), s2 = j2.norm();
			return (mOutline2[k].x*s2 / (s1 + s2) + mOutline2[n].x*s1 / (s1 + s2))*mExtensionli;
		}
	}
	return mOutline2[k].x*mExtensionli;//使用该参数加肥变形出来的模型极其不连续，效果并不理想只能在边界附近进行小部分变形适用
}
