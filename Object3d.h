/*
	Object3d - Library for manipulating 3d objects in realtime
	Created by Noel Bundy, August 31, 2014
	Released into the public domain.
*/
#ifndef Object3d_h
#define Object3d_h

typedef struct { 
	int ID;
	float depth;
} DepthMap;

typedef struct { 
	int a;
	int b;
	int c;
} FaceList;

class Point3d
{
	public:
	long lx,ly,lz,lt;
	long wx,wy,wz,wt;
	long ax,ay,az,at;
	long sx,sy,sz,st;
	Point3d();
	int operator == (Point3d &V);
	int operator != (Point3d &V);
	Point3d operator -  (Point3d &p);
	Point3d operator +  (Point3d &p);
	Point3d operator *  (Point3d &p);
	Point3d operator /  (Point3d &p);
	Point3d &operator -= (Point3d &p);
	Point3d &operator += (Point3d &p);
	Point3d &operator *= (Point3d &p);
	Point3d &operator /= (Point3d &p);
	int operator == (double p);
	int operator != (double p);
	Point3d operator -  (double p);
	Point3d operator +  (double p);
	Point3d operator *  (double p);
	Point3d operator /  (double p);
	Point3d &operator -= (double p);
	Point3d &operator += (double p);
	Point3d &operator *= (double p);
	Point3d &operator /= (double p);
	void Print();
};

class Matrix3d
{
	public:
	double Matrix[4][4];
	Matrix3d();
	void MatrixReset();
	void MatrixIdentity();
	void MatrixCopy(Matrix3d &M);
	void MatrixMult(Matrix3d &M1, Matrix3d &M2);

};

class Object3d
{
public:
	Point3d * mesh;
	int nVertices;
	int nFaces;
	FaceList * faces;
	DepthMap * depthMap;
	float minDepth;
	float maxDepth;
	Object3d(int numVerts, int numFaces);
	void init();
	void Translate(float x, float y, float z);
	void Rotate(float x, float y, float z);
	void Scale(float scale);
	Point3d ChangeLocalObject(Point3d &p);
	Point3d ChangeObjectPoint(Point3d &p);
	Point3d getNormal(Point3d a, Point3d b, Point3d c);
	Point3d Normalize(Point3d V);
	int Render(float Center, float FOV,int faceID, char *dir);
	void sortDepthMap();
	void ApplyTransforms();
	void loadMesh(float verts[][3], int faceArray[][3]);
	Matrix3d matrix, Rmat, rmatrix, objectmatrix;
	char local;
};
#endif