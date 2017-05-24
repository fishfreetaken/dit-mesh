#include "stdafx.h"
#include "CommonFun.h"

void outputMat(cv::Mat1d &outMat)
{
	for (int i = 0; i < outMat.rows; i++)
	{
		for (int j = 0; j < outMat.cols; j++)
		{
			printf("%10f", outMat[i][j]);
		}
		printf("\n");
	}
}


void writeXYZN(MyMesh &inputMesh, vector<int> &floorContour)
{
	MyMesh::Point curPoint;
	FILE * fp = fopen("data/zzz.xyzn", "w");
	for (unsigned int i = 0; i < floorContour.size(); i++)
	{
		curPoint = inputMesh.point(inputMesh.vertex_handle(floorContour[i]));
		fprintf(fp, "%f %f %f 0 0 1\n", curPoint[0], curPoint[1], curPoint[2]);
	}
	fclose(fp);
}
void writeXYZN(vector<cv::Point3d> &surfContour, const char *fileName)
{
	MyMesh::Point curPoint;
	
	FILE * fp = fopen(fileName, "w");

	for (unsigned int i = 0; i < surfContour.size(); i++)
		fprintf(fp, "%10f %10f %10f 0 0 1\n", surfContour[i].x, surfContour[i].y, surfContour[i].z);

	fclose(fp);
}
void writeXYZ(vector<cv::Point3d> &surfContour)
{
	MyMesh::Point curPoint;
	FILE * fp = fopen("zzz.xyz", "w");
	fprintf(fp, "3, %zd\n", surfContour.size());
	for (unsigned int i = 0; i < surfContour.size(); i++)
	{
		fprintf(fp, "%10f %10f %10f\n", surfContour[i].x, surfContour[i].y, surfContour[i].z);
	}
	fclose(fp);
}

double findContour(MyMesh &inputMesh, cv::Vec4d A, vector<cv::Point3d> &surfContour)
{
	surfContour.clear();
	cv::Point3d   newPt;
	MyMesh::Point curPoint;
	vector<double> disToPlane(inputMesh.n_vertices());

	//遍历顶点，计算顶点到平面距离
	for (MyMesh::VertexIter v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		disToPlane[v_it->idx()] = A[0] * curPoint[0] + A[1] * curPoint[1] + A[2] * curPoint[2] + A[3]; //不失一般性,未做正则化
	}
	//遍历棱，如有棱的两端位于平面两侧，则求该棱与平面的交点，作为轮廓点起始点, break
	MyMesh::HalfedgeHandle startI;
	MyMesh::HalfedgeIter v_it;
	MyMesh::VertexHandle fPointH, tPointH;

	for (v_it = inputMesh.halfedges_begin(); v_it != inputMesh.halfedges_end(); v_it++)
	{
		fPointH = inputMesh.from_vertex_handle(*v_it);
		tPointH = inputMesh.to_vertex_handle(*v_it);
		double a = disToPlane[fPointH.idx()];
		double b = disToPlane[tPointH.idx()];
		if (a * b < 0)
		{
			MyMesh::Point fp = inputMesh.point(fPointH);
			MyMesh::Point tp = inputMesh.point(tPointH);
			double dt = a / (a - b);
			newPt.x = fp[0] + dt * (tp[0] - fp[0]);
			newPt.y = fp[1] + dt * (tp[1] - fp[1]);
			newPt.z = fp[2] + dt * (tp[2] - fp[2]);
			surfContour.push_back(newPt);
			startI = *v_it;
			break;
		}
	}
	if (v_it == inputMesh.halfedges_end())//没有相交点
		return 0;
	MyMesh::HalfedgeHandle he = startI; //he为迭代量
	he = inputMesh.next_halfedge_handle(he);
	he = inputMesh.next_halfedge_handle(he);
	he = inputMesh.next_halfedge_handle(he);
	while (1)
	{
		he = inputMesh.opposite_halfedge_handle(he);
		he = inputMesh.next_halfedge_handle(he);
		if (he == startI)
			break;
		fPointH = inputMesh.from_vertex_handle(he);
		tPointH = inputMesh.to_vertex_handle(he);
		MyMesh::Point fp = inputMesh.point(fPointH);
		MyMesh::Point tp = inputMesh.point(tPointH);
		double a = disToPlane[fPointH.idx()];
		double b = disToPlane[tPointH.idx()];
		if (a * b < 0)
		{
			MyMesh::Point fp = inputMesh.point(fPointH);
			MyMesh::Point tp = inputMesh.point(tPointH);
			double dt = a / (a - b);
			newPt.x = fp[0] + dt * (tp[0] - fp[0]);
			newPt.y = fp[1] + dt * (tp[1] - fp[1]);
			newPt.z = fp[2] + dt * (tp[2] - fp[2]);
			surfContour.push_back(newPt);
			continue;
		}
		he = inputMesh.next_halfedge_handle(he);
		if (he == startI)
			break;
		fPointH = inputMesh.from_vertex_handle(he);
		tPointH = inputMesh.to_vertex_handle(he);
		fp = inputMesh.point(fPointH);
		tp = inputMesh.point(tPointH);
		a = disToPlane[fPointH.idx()];
		b = disToPlane[tPointH.idx()];
		if (a * b < 0)
		{
			MyMesh::Point fp = inputMesh.point(fPointH);
			MyMesh::Point tp = inputMesh.point(tPointH);
			double dt = a / (a - b);
			newPt.x = fp[0] + dt * (tp[0] - fp[0]);
			newPt.y = fp[1] + dt * (tp[1] - fp[1]);
			newPt.z = fp[2] + dt * (tp[2] - fp[2]);
			surfContour.push_back(newPt);
			continue;
		}
		break;
	}
	double contourLength = 0; //计算轮廓周长
	for (unsigned int i = 0; i < surfContour.size() - 1; i++)
	{
		double dx = surfContour[i + 1].x - surfContour[i].x;
		double dy = surfContour[i + 1].y - surfContour[i].y;
		double dz = surfContour[i + 1].z - surfContour[i].z;
		contourLength += sqrt(dx * dx + dy * dy + dz * dz);
	}
	double dx = surfContour[0].x - surfContour[surfContour.size() - 1].x;
	double dy = surfContour[0].y - surfContour[surfContour.size() - 1].y;
	double dz = surfContour[0].z - surfContour[surfContour.size() - 1].z;
	contourLength += sqrt(dx * dx + dy * dy + dz * dz);
	return contourLength;
}


//如涉及到轮廓点排序,对凸的截面OK,非凸面不可以,一般不用
double findContourShortBreak(MyMesh &inputMesh, cv::Vec4d A, vector<cv::Point3d> &retContour)
{
	vector<cv::Point3d> surfContour;
	MyMesh::Point curPoint;
	vector<double> disToPlane(inputMesh.n_vertices());
	int vIndex = 0;
	//遍历顶点，计算顶点到平面距离
	for (MyMesh::VertexIter v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		disToPlane[v_it->idx()] = A[0] * curPoint[0] + A[1] * curPoint[1] + A[2] * curPoint[2] + A[3];
	}
	//遍历棱，如有棱的两端位于平面两侧，则求该棱与平面的交点，作为轮廓点，正反棱问题是否存在？
	for (MyMesh::HalfedgeIter v_it = inputMesh.halfedges_begin(); v_it != inputMesh.halfedges_end(); v_it++)
	{
		MyMesh::VertexHandle	fPointH = inputMesh.from_vertex_handle(*v_it);
		MyMesh::VertexHandle	tPointH = inputMesh.to_vertex_handle(*v_it);
		if (disToPlane[fPointH.idx()] * disToPlane[tPointH.idx()] <= 0)
		{
			MyMesh::Point fp = inputMesh.point(fPointH);
			MyMesh::Point tp = inputMesh.point(tPointH);
			double dt = -(disToPlane[fPointH.idx()]) / (disToPlane[tPointH.idx()] - disToPlane[fPointH.idx()]);
			cv::Point3d newPt;
			newPt.x = fp[0] + dt * (tp[0] - fp[0]);
			newPt.y = fp[1] + dt * (tp[1] - fp[1]);
			newPt.z = fp[2] + dt * (tp[2] - fp[2]);
			surfContour.push_back(newPt);
		}
	}
	if (surfContour.size()< 3)
		return -1;
	//一般来说到此结束，为了求轮廓周长，需要把轮廓点排序
	//求中心，根据角度依次排序
	double sumX = 0, sumY = 0, sumZ = 0;
	for (unsigned int i = 0; i < surfContour.size(); i++)
	{
		sumX += surfContour[i].x;
		sumY += surfContour[i].y;
		sumZ += surfContour[i].z;
	}
	cv::Point3d centerP = cv::Point3d(sumX / surfContour.size(), sumY / surfContour.size(), sumZ / surfContour.size());

	vector<double> surfTheta(surfContour.size());
	cv::Point3d  initV = surfContour[0] - centerP;
	initV /= sqrt(initV.x * initV.x + initV.y * initV.y + initV.z * initV.z);
	surfTheta[0] = 0;
	for (unsigned int i = 0; i < surfContour.size(); i++)
	{
		cv::Point3d  curV = surfContour[i] - centerP;
		surfTheta[i] = acos(curV.ddot(initV) / sqrt(curV.x * curV.x + curV.y * curV.y + curV.z * curV.z));
		cv::Point3d  tmpNorm = initV.cross(curV);
		if (tmpNorm.ddot(cv::Point3d(A[0], A[1], A[2])) < 0)  //与平面法向量相反，需要调整
			surfTheta[i] = CV_2PI - surfTheta[i];
	}

	vector<int> dst;
	cv::sortIdx(surfTheta, dst, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);
	retContour.resize(surfContour.size());//调整轮廓点顺序，依次计算周长
	for (unsigned int i = 0; i < surfContour.size(); i++)
		retContour[i] = surfContour[dst[i]];


	double contourLength = 0; //计算轮廓周长
	for (unsigned int i = 0; i < surfContour.size() - 1; i++)
	{
		double dx = retContour[i + 1].x - retContour[i].x;
		double dy = retContour[i + 1].y - retContour[i].y;
		double dz = retContour[i + 1].z - retContour[i].z;
		contourLength += sqrt(dx * dx + dy * dy + dz * dz);
	}
	double dx = retContour[0].x - retContour[surfContour.size() - 1].x;
	double dy = retContour[0].y - retContour[surfContour.size() - 1].y;
	double dz = retContour[0].z - retContour[surfContour.size() - 1].z;
	contourLength += sqrt(dx * dx + dy * dy + dz * dz);


	//旋转到XOY平面，不管平移量
	cv::Point3d planeN = cv::Point3d(A[0], A[1], A[2]);
	planeN /= sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
	cv::Point3d zN = cv::Point3d(0, 0, 1);
	double theta = acos(planeN.dot(zN));
	cv::Point3d newZ = planeN.cross(zN);
	newZ /= sqrt(newZ.x * newZ.x + newZ.y * newZ.y + newZ.z * newZ.z);
	newZ *= theta;
	cv::Mat1d rotatM = cv::Mat1d::eye(3, 3);
	cv::Mat1d rotatV(3, 1);
	rotatV[0][0] = newZ.x; rotatV[1][0] = newZ.y; rotatV[2][0] = newZ.z;
	cv::Rodrigues(rotatV, rotatM);
	for (unsigned int i= 0; i< retContour.size(); i++)
	{
		cv::Vec3d curPoint = retContour[i];
		double t = A[0] * curPoint[0] + A[1] * curPoint[1] + A[2] * curPoint[2] + A[3];
		double x = rotatM[0][0] * curPoint[0] + rotatM[0][1] * curPoint[1] + rotatM[0][2] * curPoint[2];
		double y = rotatM[1][0] * curPoint[0] + rotatM[1][1] * curPoint[1] + rotatM[1][2] * curPoint[2];
		double z = rotatM[2][0] * curPoint[0] + rotatM[2][1] * curPoint[1] + rotatM[2][2] * curPoint[2];
		printf("%10f %10f %10f %10f\n", t, x, y, z);
	}

	return contourLength;
}
//
//void findFloorContour(MyMesh &inputMesh, vector<int> &floorContour)
//{
//	MyMesh::VertexIter  v_it, v_end(inputMesh.vertices_end());
//	MyMesh::VertexIter	v_it_s[3];
//	MyMesh::Point curPoint;
//	double minX, maxX;
//	cout << "Finding Floor Contour..." << endl;
//	v_it = inputMesh.vertices_begin();
//	minX = maxX = inputMesh.point(*v_it)[0];
//	cv::Mat1d tmpMat(3, inputMesh.n_vertices());
//	for (; v_it != v_end; v_it++)
//	{
//		curPoint = inputMesh.point(*v_it);
//		tmpMat[0][v_it->idx()] = curPoint[0];
//		tmpMat[1][v_it->idx()] = curPoint[1];
//		tmpMat[2][v_it->idx()] = curPoint[2];
//	}
//	cv::Mat1d pointMat(3, inputMesh.n_vertices());
//
//	cv::Mat1i dst;
//	cv::sortIdx(tmpMat.row(0), dst, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);
//
//	for (int i = 0; i < pointMat.cols; i++)
//	{
//		pointMat[0][i] = tmpMat[0][dst[0][i]];
//		pointMat[1][i] = tmpMat[1][dst[0][i]];
//		pointMat[2][i] = tmpMat[2][dst[0][i]];
//	}
//
//	double xSampleS = 0.01;
//	int xSampleL = int((maxX - minX) / xSampleS + 1);
//	vector<vector<int>> xSamplePtV(xSampleL);
//	v_it = inputMesh.vertices_begin();
//	for (; v_it != v_end; v_it++)
//	{
//		curPoint = inputMesh.point(*v_it);
//		int xIndex = int((curPoint[0] - minX) / xSampleS);
//		xSamplePtV[xIndex].push_back(v_it->idx());
//	}
//	floorContour.clear();
//	int sum = 0;
//	for (int i = 0; i < xSampleL; i++)
//	{
//		sum += xSamplePtV[i].size();
//		if (xSamplePtV[i].size() == 0)
//			continue;
//		double minZ = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][0]))[2];
//		for (unsigned int j = 1; j < xSamplePtV[i].size(); j++)
//		{
//			curPoint = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][j]));
//			if (curPoint[2] < minZ)
//				minZ = curPoint[2];
//		}
//		int left = 0;
//		int right = 0;
//		for (unsigned int j = 0; j < xSamplePtV[i].size(); j++)
//		{
//			curPoint = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][j]));
//
//			if ((curPoint[2] - minZ) > 5)
//				continue;
//			floorContour.push_back(xSamplePtV[i][j]);			
//		}
//	}
//	//int n = inputMesh.n_vertices();
//	//vector<cv::Point2d> cloudMat;
//	//vector<int> edgeL;
//	//for (unsigned int i = 0; i < floorContour.size(); i++)
//	//{
//	//	curPoint = inputMesh.point(inputMesh.vertex_handle(floorContour[i]));
//	//	cloudMat.push_back(cv::Point2d(curPoint[0], curPoint[1]));
//	//}
//	//findHull(cloudMat, edgeL);
//
//	writeXYZN(inputMesh, floorContour);
//}

void findFloorContour(MyMesh &inputMesh, vector<int> &floorContour)
{
	MyMesh::VertexIter  v_it, v_end(inputMesh.vertices_end());
	MyMesh::VertexIter	v_it_s[3];
	MyMesh::Point curPoint;
	double minX, maxX;
	cout << "Finding Floor Contour..." << endl;
	v_it = inputMesh.vertices_begin();
	minX = maxX = inputMesh.point(*v_it)[0];
	for (; v_it != v_end; v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		if (curPoint[0] > maxX)
			maxX = curPoint[0];
		if (curPoint[0] < minX)
			minX = curPoint[0];
	}
	double xSampleS = 1;
	int xSampleL = int((maxX - minX) / xSampleS + 1);
	vector<vector<int>> xSamplePtV(xSampleL);
	v_it = inputMesh.vertices_begin();
	for (; v_it != v_end; v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		int xIndex = int((curPoint[0] - minX) / xSampleS);
		xSamplePtV[xIndex].push_back(v_it->idx());
	}
	floorContour.clear();
	int sum = 0;
	for (int i = 0; i < xSampleL; i++)
	{
		sum += xSamplePtV[i].size();
		if (xSamplePtV[i].size() == 0)
			continue;
		double minZ = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][0]))[2];
		for (unsigned int j = 1; j < xSamplePtV[i].size(); j++)
		{
			curPoint = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][j]));
			if (curPoint[2] < minZ)
				minZ = curPoint[2];
		}
		int left = 0;
		int right = 0;
		for (unsigned int j = 0; j < xSamplePtV[i].size(); j++)
		{
			curPoint = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][j]));

			if ((curPoint[2] - minZ) > 5)
				continue;
			floorContour.push_back(xSamplePtV[i][j]);			
		}
	}
	//int n = inputMesh.n_vertices();
	//vector<cv::Point2d> cloudMat;
	//vector<int> edgeL;
	//for (unsigned int i = 0; i < floorContour.size(); i++)
	//{
	//	curPoint = inputMesh.point(inputMesh.vertex_handle(floorContour[i]));
	//	cloudMat.push_back(cv::Point2d(curPoint[0], curPoint[1]));
	//}
	//findHull(cloudMat, edgeL);

	writeXYZN(inputMesh, floorContour);
}

//void findFloorContour(MyMesh &inputMesh, vector<int> &floorContour)
//{
//	MyMesh::VertexIter  v_it, v_end(inputMesh.vertices_end());
//	MyMesh::VertexIter	v_it_s[3];
//	MyMesh::Point curPoint;
//	double minX, maxX, minY, maxY;
//	cout << "Finding Floor Contour..." << endl;
//	v_it = inputMesh.vertices_begin();
//	minX = maxX = inputMesh.point(*v_it)[0];
//	minY = maxY = inputMesh.point(*v_it)[1];
//	for (; v_it != v_end; v_it++)
//	{
//		curPoint = inputMesh.point(*v_it);
//		if (curPoint[0] > maxX)
//			maxX = curPoint[0];
//		if (curPoint[0] < minX)
//			minX = curPoint[0];
//		if (curPoint[1] > maxY)
//			maxY = curPoint[1];
//		if (curPoint[1] < minY)
//			minY = curPoint[1];
//		inputMesh.set_point(*v_it, curPoint);
//		v_it->idx();
//	}
//	double xSampleS = 0.05;
//	int xSampleL = int((maxX - minX) / xSampleS + 1);
//	vector<vector<int>> xSamplePtV(xSampleL);
//	v_it = inputMesh.vertices_begin();
//	for (; v_it != v_end; v_it++)
//	{
//		curPoint = inputMesh.point(*v_it);
//		int xIndex = int((curPoint[0] - minX) / xSampleS);
//		xSamplePtV[xIndex].push_back(v_it->idx());
//	}
//	floorContour.clear();
//	for (int i = 0; i < xSampleL; i++)
//	{
//		if (xSamplePtV[i].size() == 0)
//			continue;
//		double minZ = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][0]))[2];
//		for (unsigned int j = 1; j < xSamplePtV[i].size(); j++)
//		{
//			curPoint = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][j]));
//			if (curPoint[2] < minZ)
//				minZ = curPoint[2];
//		}
//		minY = maxY = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][0]))[1];
//		int left = 0;
//		int right = 0;
//		for (unsigned int j = 0; j < xSamplePtV[i].size(); j++)
//		{
//			curPoint = inputMesh.point(inputMesh.vertex_handle(xSamplePtV[i][j]));
//			floorContour.push_back(curPoint);
//
//			if ((curPoint[2] - minZ) > 5)
//				continue;
//			if (curPoint[1] < minY)
//			{
//				minY = curPoint[1];
//				left = j;
//			}
//			if (curPoint[1] > maxY)
//			{
//				maxY = curPoint[1];
//				right = j;
//			}
//		}
//		floorContour.push_back(xSamplePtV[i][left]);
//		floorContour.push_back(xSamplePtV[i][right]);
//	}
//	vector<cv::Point2d> cloudMat;
//	vector<int> edgeL;
//	for (unsigned int i = 0; i < floorContour.size(); i++)
//	{
//		curPoint = inputMesh.point(inputMesh.vertex_handle(floorContour[i]));
//		cloudMat.push_back(cv::Point2d(curPoint[0], curPoint[1]));
//	}
//	//findHull(cloudMat, edgeL);
//
//	FILE * fp = fopen("zzz.xyzn", "w");
//	for (unsigned int i = 0; i < edgeL.size(); i++)
//	{
//		fprintf(fp, "%f %f 10 0 0 1\n", cloudMat[edgeL[i]].x, cloudMat[edgeL[i]].y);
//	}
//	fclose(fp);
//
//	writeXYZN(inputMesh, floorContour);
//}
void rotateMesh(MyMesh &inputMesh)
{
	cv::Mat1d rotatM = cv::Mat1d::eye(3, 3);
	cv::Mat1d rotatV(3, 1);
	rotatV[0][0] = 0; rotatV[1][0] = 1; rotatV[2][0] = 0;
	cv::Rodrigues(rotatV, rotatM);
	for (int i = 0; i < 3; i++)
	{
		printf("%f\t%f\t%f\n", rotatM[i][0], rotatM[i][1], rotatM[i][2]);
	}
	MyMesh::VertexIter  v_it, v_end(inputMesh.vertices_end());
	MyMesh::VertexIter	v_it_s[3];
	MyMesh::Point curPoint;
	cout << "Finding Floor Contour..." << endl;
	v_it = inputMesh.vertices_begin();
	for (; v_it != v_end; v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		double x = rotatM[0][0] * curPoint[0] + rotatM[0][1] * curPoint[1] + rotatM[0][2] * curPoint[2];
		double y = rotatM[1][0] * curPoint[0] + rotatM[1][1] * curPoint[1] + rotatM[1][2] * curPoint[2];
		double z = rotatM[2][0] * curPoint[0] + rotatM[2][1] * curPoint[1] + rotatM[2][2] * curPoint[2];
		curPoint[0] = x;
		curPoint[1] = y;
		curPoint[2] = z;
		inputMesh.set_point(*v_it, curPoint);
	}
}

void rotateModel(MyMesh &inputMesh)
{
	cv::Mat1d rotatM = cv::Mat1d::eye(3, 3);
	cv::Mat1d rotatV(3, 1);
	rotatV[0][0] = 0; rotatV[1][0] = 1; rotatV[2][0] = 0;
	cv::Rodrigues(rotatV, rotatM);
	for (int i = 0; i < 3; i++)
		printf("%f\t%f\t%f\n", rotatM[i][0], rotatM[i][1], rotatM[i][2]);

	MyMesh::VertexIter  v_it, v_end(inputMesh.vertices_end());
	MyMesh::Point curPoint;
	cout << "Finding Floor Contour..." << endl;
	v_it = inputMesh.vertices_begin();
	for (; v_it != v_end; v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		double x = rotatM[0][0] * curPoint[0] + rotatM[0][1] * curPoint[1] + rotatM[0][2] * curPoint[2];
		double y = rotatM[1][0] * curPoint[0] + rotatM[1][1] * curPoint[1] + rotatM[1][2] * curPoint[2];
		double z = rotatM[2][0] * curPoint[0] + rotatM[2][1] * curPoint[1] + rotatM[2][2] * curPoint[2];
		curPoint[0] = x;
		curPoint[1] = y;
		curPoint[2] = z;
		inputMesh.set_point(*v_it, curPoint);
	}
}

//PCA寻找主轴
void findMainO(MyMesh &inputMesh)
{
	MyMesh::VertexIter  v_it  = inputMesh.vertices_begin(); 
	MyMesh::VertexIter  v_end = inputMesh.vertices_end();
	MyMesh::Point curPoint;
	cv::Mat1d cloudMat(3, inputMesh.n_vertices());
	double sumX = 0, sumY = 0, sumZ = 0;
	for (; v_it != v_end; v_it++)
	{
		curPoint = inputMesh.point(*v_it);
		cloudMat[0][v_it->idx()] = curPoint[0];
		cloudMat[1][v_it->idx()] = curPoint[1];
		cloudMat[2][v_it->idx()] = curPoint[2];
	}
	cv::PCA cloudMatMainO(cloudMat, cv::Mat(), CV_PCA_DATA_AS_COL);
	cv::Mat1d rotM3;
	cloudMatMainO.eigenvectors.copyTo(rotM3);
	outputMat(rotM3);
	cloudMat = rotM3 * cloudMat;

	cv::Mat1d rotM4 = cv::Mat1d::eye(4, 4);
	rotM3.copyTo(rotM4(cv::Rect(0, 0, 3, 3)));

	bool isversX = false;
	double minX, maxX;
	cv::minMaxIdx(cloudMat.row(0), &minX, &maxX);
	cv::Mat1d rotA = cv::Mat1d::zeros(1, 4);

	rotA[0][0] = 1; rotA[0][3] = -(minX + (maxX - minX) / 5);
	rotA = rotA * rotM4;
	cv::Vec4d planeA(rotA[0][0], rotA[0][1], rotA[0][2], rotA[0][3]);
	vector<cv::Point3d> aContour, bContour;
	double aLength = findContour(inputMesh, planeA, aContour);
	writeXYZN(aContour, "data/aContour.xyzn");

	rotA = cv::Mat1d::zeros(1, 4);
	rotA[0][0] = 1; rotA[0][3] = -(maxX - (maxX - minX) / 5);
	rotA = rotA * rotM4;
	planeA= cv::Vec4d (rotA[0][0], rotA[0][1], rotA[0][2], rotA[0][3]);
	double bLength = findContour(inputMesh, planeA, bContour);
	vector<cv::Point3d> *paC = &bContour;
	if (aLength< bLength)
	{
		paC = &aContour;
		isversX = true;
	}
	int i0, j0;
	double maxDis = -1;
	for (unsigned int i = 0; i < paC->size() - 1; i++)
	{
		for (unsigned int j = i + 1; j < paC->size(); j++)
		{
			double dx = paC->at(j).x - paC->at(i).x;
			double dy = paC->at(j).y - paC->at(i).y;
			double dis = sqrt(dx * dx + dy * dy);
			if (dis > maxDis)
			{
				maxDis = dis;
				i0 = i; j0 = j;
			}
		}
	}
	writeXYZN(*paC, "data/bContour.xyzn");

	cv::Mat1d xRotV(1, 3), xRotM(3, 3);
	xRotV[0][0] = paC->at(i0).x - paC->at(j0).x;
	xRotV[0][1] = paC->at(i0).y - paC->at(j0).y;
	xRotV[0][2] = paC->at(i0).z - paC->at(j0).z;
	xRotV = xRotV * rotM3.t();
	outputMat(xRotV);
	double dtheta = atan2(xRotV[0][2], xRotV[0][1]);
	xRotV[0][0] = dtheta;
	xRotV[0][1] = 0;
	xRotV[0][2] = 0;
	//cv::Rodrigues(xRotV, xRotM);
	//rotM3 = xRotM * rotM3;

	/*cv::Vec3d xDir = cv::Vec3d(rotM3[0][0], rotM3[0][1], rotM3[0][2]);

	if (isversX)
	{
		xDir = -xDir;
	}	cv::Vec3d zDir = xDir.cross(yDir);*/
	//memcpy(rotM3[0], &xDir[0], sizeof(double) * 3);
	//memcpy(rotM3[1], &yDir[0], sizeof(double) * 3);
	//memcpy(rotM3[2], &zDir[0], sizeof(double) * 3);
	//outputMat(rotM3);
	MyMesh::Point newPoint;

	for (v_it = inputMesh.vertices_begin(); v_it != v_end; v_it++)
	{
		curPoint = inputMesh.point(*v_it);

		newPoint[0] = rotM3[0][0] * curPoint[0] + rotM3[0][1] * curPoint[1] + rotM3[0][2] * curPoint[2];
		newPoint[1] = rotM3[1][0] * curPoint[0] + rotM3[1][1] * curPoint[1] + rotM3[1][2] * curPoint[2];
		newPoint[2] = rotM3[2][0] * curPoint[0] + rotM3[2][1] * curPoint[1] + rotM3[2][2] * curPoint[2];
		//newPoint[0] = cloudMat[0][v_it->idx()];
		//newPoint[1] = cloudMat[1][v_it->idx()];
		//newPoint[2] = cloudMat[2][v_it->idx()];

		inputMesh.set_point(*v_it, newPoint);
	}

	/*double minX, minY, minZ, maxX, maxY, maxZ;
	cv::minMaxIdx(cloudMat.row(0), &minX, &maxX);
	cv::minMaxIdx(cloudMat.row(1), &minY, &maxY);
	cv::minMaxIdx(cloudMat.row(2), &minZ, &maxZ);
	cv::Mat1d xrotM = cv::Mat1d::eye(3, 3);
	if (fabs(minX)> fabs(maxX))
	{ 
		cv::Mat1d tmprotV = cv::Mat1d::zeros(3, 1);
		tmprotV[0][2] = CV_PI / 2;
		cv::Rodrigues(tmprotV, xrotM);
	}
	cloudMat = xrotM * cloudMat;*/

//	if ((maxY - minY)> (maxZ - minZ))
//	{
//		cv::Mat1d tmprotV = cv::Mat1d::zeros(3, 1);
//		tmprotV[0][0] = CV_PI / 2;
//		cv::Mat1d tmprotM(3, 3);
//		cv::Rodrigues(tmprotV, tmprotM);
//		outputMat(tmprotM);
////		xrotM = tmprotM * xrotM;
//	}

}

int find3rdPoint(vector<cv::Point2d> &cloudMat, vector<int> validIndex, cv::Mat1d &disMat, vector<int> &edgeL)
{
	if (edgeL.size() < 2)
		return -1;
	//outputMat(disMat);
	int startP, endP;
	for (int k= 0; k< edgeL.size() - 1; k++)
	{
		startP = edgeL[k];
		endP   = edgeL[k + 1];
		cv::Point2d halfLine = cloudMat[endP] - cloudMat[startP];
		double maxDis = -1;
		int    retIndex = -1;
		for (unsigned int i = 0; i < validIndex.size(); i++)
		{
			if (validIndex[i] == startP || validIndex[i] == endP)
				continue;
			cv::Point2d newLine = cloudMat[validIndex[i]] - cloudMat[endP];
			if (halfLine.cross(newLine) <= 0)
				continue;
			double p = (disMat[startP][endP] + disMat[validIndex[i]][startP] + disMat[validIndex[i]][endP]) / 2;
			double dis = sqrt(p * (p - disMat[startP][endP]) * (p - disMat[validIndex[i]][startP]) * (p - disMat[validIndex[i]][endP])) / disMat[startP][endP];
			if (dis > maxDis)
			{
				maxDis = dis;
				retIndex = validIndex[i];
			}
		}
		if (retIndex> -1)
		{
			edgeL.insert(edgeL.begin() + k + 1, retIndex);
			k = k - 1;
		}
	}
	return 0;
}
//寻找二维点云的凸包
//返回有序的索引点
void findHull(vector<cv::Point2d> cloudMat, vector<int> &edgeL)
{
	if (cloudMat.size() < 3)
		return;

	int i0, j0;
	double maxDis = -1;
	cv::Mat1d disMat = cv::Mat1d::zeros(cloudMat.size(), cloudMat.size());
	for (unsigned int i= 0; i< cloudMat.size() - 1; i++)
	{
		for (unsigned int j = i + 1; j < cloudMat.size(); j++)
		{
			double dx = cloudMat[j].x - cloudMat[i].x;
			double dy = cloudMat[j].y - cloudMat[i].y;
			double dis = sqrt(dx * dx + dy * dy);
			disMat[i][j] = disMat[j][i] = dis;
			if (dis> maxDis)
			{
				maxDis = dis;
				i0 = i; j0 = j;
			}
		}
	}
	vector<int>  hullIndex;
	hullIndex.push_back(i0);
	hullIndex.push_back(j0);

	vector<int>  edgeR;
	edgeL.push_back(i0); edgeL.push_back(j0);
	edgeR.push_back(j0); edgeR.push_back(i0);

	vector<int> validIndexL, validIndexR;
	for (int i= 0; i< cloudMat.size(); i++)
	{
		if (i == i0 || i == j0)
			continue;

		cv::Point2d halfLine = cloudMat[j0] - cloudMat[i0];
		cv::Point2d newLine  = cloudMat[i]  - cloudMat[j0];
		if (halfLine.cross(newLine) > 0)
			validIndexL.push_back(i);
		else
			validIndexR.push_back(i);
	}

	find3rdPoint(cloudMat, validIndexL, disMat, edgeL);
	find3rdPoint(cloudMat, validIndexR, disMat, edgeR);
	if (edgeR.size() > 2)
	{
		for (unsigned int j = 1; j < edgeR.size() - 1; j++)
			edgeL.push_back(edgeR[j]);
	}
	FILE * fp = fopen("zhull.xyzn", "w");
	for (unsigned int i = 0; i < edgeL.size(); i++)
	{
		fprintf(fp, "%f %f 10 0 0 1\n", cloudMat[edgeL[i]].x, cloudMat[edgeL[i]].y);
	}
	fclose(fp);
}