#include "stdafx.h"
#include "ManageObj.h"

ManageObj::ManageObj(char * name,vector<Vector3f>*p,vector<Vector3i>*t):
	mrfilename(name),
	mvpoint(p),
	mptri(t)
{
}
/*read xyz file*/
ManageObj::ManageObj(char * name, vector<Vector3f>*p) :
	mrfilename(name),
	mvpoint(p)
{
}
ManageObj::ManageObj(char * name) :
	mrfilename(name)
{
	OpenObj();
}


ManageObj::~ManageObj()
{
}

int ManageObj::OpenObj() {
	FILE * pFile;
	long lSize;
	
	size_t result;
	fopen_s(&pFile, mrfilename, "r");
	if (pFile == NULL)
	{
		cout << "File error" << endl;
		return 1;
	}
	fseek(pFile, 0, SEEK_END);
	lSize = ftell(pFile);
	rewind(pFile);

	mbuffer = (char*)malloc(sizeof(char)*lSize);
	if (mbuffer == NULL)
	{
		cout << "Memory error" << endl;
		return 2;
	}

	result = fread(mbuffer, 1, lSize, pFile);
	if (result != lSize)
	{
		fputs("Reading error", stderr);
		return 0;
	}
	fclose(pFile);
	return 0;
}

void ManageObj::ReadObj() {
	if (OpenObj()) {
		return;
	}
	stringstream ss(mbuffer);
	string s;
	ios::sync_with_stdio(false);
	printf("Now is reading OBJ file %s ...", mrfilename);
	while (getline(ss, s)) {
		if (s[0] == 'v'){
			mvpoint->push_back(StringToFloatV(s));
		}
		if (s[0] == 'f') {
			mptri->push_back(StringToIntV(s));
		}
	}
	ios::sync_with_stdio(true);
	FreeBuf();
}

void ManageObj::ReadStl() {
	const char* p = mbuffer;
	char name[80];
	int i, j;
	memcpy(name, p, 80);
	p += 80;
	uint unTriangles = cpyint(p);
	for (i = 0; i < unTriangles; i++)  //mMeshTri
	{
		VertexTriMesh a;
		/*三角面片法向量*/
		a.a=Vector3f(cpyfloat(p),cpyfloat(p),cpyfloat(p));
		a.i = i;
		for (j = 0; j < 3; j++)//读取三顶点
		{
			Vector3f mp(cpyfloat(p), cpyfloat(p), cpyfloat(p));
			if (mVMap.find(mp) == mVMap.end()) {
				mVMap[mp] = mvpoint->size();
				mvpoint->push_back(mp);
			}
			a.p[j] = mVMap[mp];
		}
		p += 2;//跳过尾部标志
	}
}

Vector3f ManageObj::StringToFloatV(string s) {
	float a[3];
	char c = ' ';
	sscanf(s.c_str(),"%c %f %f %f\n",&c,&a[0],&a[1],&a[2]);
	return Vector3f(a[0],a[1],a[2]);
}

Vector3i ManageObj::StringToIntV(string s) {
	int a[3];
	char c = ' ';
	sscanf(s.c_str(), "%c %d %d %d\n", &c, &a[0], &a[1], &a[2]);
	return Vector3i(a[0], a[1], a[2]);
}

void ManageObj::ReadXyz() {
	stringstream ss(mbuffer);
	string s;
	ios::sync_with_stdio(false);
	cout << "Now is reading Xyz file....." << endl;
	while (getline(ss, s)) {
		mvpoint->push_back(StringFloat2(s));
	}
	ios::sync_with_stdio(true);
	FreeBuf();
}

Vector3f ManageObj::StringFloat2(string s) {
	float  d[3];
	sscanf(s.c_str(), "%f %f %f\n",d,d+1,d+2);
	return Vector3f(d[0], d[1], d[2]);
}

void ManageObj::ReadMeshPoints(vector<MyMesh::Point>& ab) {
	stringstream ss(mbuffer);
	string s;
	ios::sync_with_stdio(false);
	cout << "Now is reading ReadMeshPoints file....." << endl;
	while (getline(ss, s)) {
		ab.push_back(StringFloat3(s));
	}
	ios::sync_with_stdio(true);
	FreeBuf();
}

MyMesh::Point ManageObj::StringFloat3(string s) {
	float  d[3];
	sscanf(s.c_str(), "%f,%f,%f\n", d, d + 1, d + 2);
	return MyMesh::Point(d[0], d[1], d[2]);
}


void ManageObj::OutFileOutlinePointXyz(vector<Vector3f>* vp,char * outfilename) {
	FILE *fp;
	fopen_s(&fp, outfilename, "w");
	for (unsigned int i = 0; i < (*vp).size(); i++) {
		fprintf(fp, "%f %f %f\n", (*vp)[i][0], (*vp)[i][1], (*vp)[i][2]);
	}
	fclose(fp);
}

/*重构！*/
void ManageObj::OutFilePointObj(vector<Vector3f>* vp, const char * outfilenam) {
	FILE *fp;
	fopen_s(&fp, outfilenam, "w");
	fprintf(fp,"# Studio\n");
	fprintf(fp, "g Point_Model_1\n");
	for (unsigned int i = 0; i < (*vp).size(); i++) {
		fprintf(fp, "v %f %f %f\n", (*vp)[i][0], (*vp)[i][1], (*vp)[i][2]);
		fprintf(fp, "p %d\n", (i+1));
	}
	fprintf(fp, "# end of file\n");
	fclose(fp);
}

void ManageObj::OutFilePointObj(vector<float> *vp, const char * outfilenam) {
	FILE *fp;
	fopen_s(&fp, outfilenam, "w");
	for (unsigned int i = 0; i < (*vp).size(); i++) {
		fprintf(fp,"%d %f\n", i, (*vp)[i]);
	}
	fprintf(fp, "# end of file\n");
	fclose(fp);
}