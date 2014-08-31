#include "Object3d.h"
#include "math.h"


Point3d::Point3d()
 { lx=ly=lz=ax=ay=az=at=wx=wy=wz=sx=sy=sz=0;
   lt=wt=st=1;
 }

 int Point3d::operator == (Point3d &p)
{int rvalue=0;
 if(p.lx==lx)
  if(p.ly==ly)
   if(p.lz==lz)
     rvalue=1;
 return rvalue;
}

int Point3d::operator != (Point3d &p)
{ int rvalue=0;
  if((p.lx!=lx)||(p.ly != ly)||(p.lz != lz))
     rvalue=1;
  return rvalue;
}

Point3d Point3d::operator - ( Point3d &p )
  {Point3d Temp;
   Temp.lx = lx - p.lx;
   Temp.ly = ly - p.ly;
   Temp.lz = lz - p.lz;

   return Temp;
  }

Point3d Point3d::operator + ( Point3d &p )
  {Point3d Temp;
   Temp.lx = lx + p.lx;
   Temp.ly = ly + p.ly;
   Temp.lz = lz + p.lz;

   return Temp;
  }
Point3d Point3d::operator * ( Point3d &p )
  {Point3d Temp;
   Temp.lx = lx * p.lx;
   Temp.ly = ly * p.ly;
   Temp.lz = lz * p.lz;


   return Temp;
  }
Point3d Point3d::operator / ( Point3d &p )
  {Point3d Temp;
   Temp.lx = lx / p.lx;
   Temp.ly = ly / p.ly;
   Temp.lz = lz / p.lz;

   return Temp;
  }

Point3d &Point3d::operator -= ( Point3d &p )
  {lx -= p.lx;
   ly -= p.ly;
   lz -= p.lz;


   return *this;
  }
Point3d &Point3d::operator += ( Point3d &p )
  {lx += p.lx;
   ly += p.ly;
   lz += p.lz;

   return *this;
  }
Point3d &Point3d::operator *= ( Point3d &p )
  {lx *= p.lx;
   ly *= p.ly;
   lz *= p.lz;


   return *this;
  }
Point3d &Point3d::operator /= ( Point3d &p )
  {lx /= p.lx;
   ly /= p.ly;
   lz /= p.lz;


   return *this;
  }

Point3d Point3d::operator - ( double Value )
  {Point3d Temp;
   Temp.lx = (long)(lx- Value);
   Temp.ly = (long)(ly - Value);
   Temp.lz = (long)(lz - Value);


   return Temp;
  }
Point3d Point3d::operator + ( double Value )
  {Point3d Temp;
   Temp.lx = (long)(lx + Value);
   Temp.ly = (long)(ly + Value);
   Temp.lz = (long)(lz + Value);

   return Temp;
  }
Point3d Point3d::operator * ( double Value )
  {Point3d Temp;
   Temp.lx = (long)(lx * Value);
   Temp.ly = (long)(ly * Value);
   Temp.lz = (long)(lz * Value);

   return Temp;
  }
Point3d Point3d::operator / ( double Value )
  {Point3d Temp;
   Temp.lx = (long)(lx / Value);
   Temp.ly = (long)(ly / Value);
   Temp.lz = (long)(lz / Value);

   return Temp;
  }
  
Point3d &Point3d::operator -= ( double Value )
  {lx -= (long)Value;
   ly -= (long)Value;
   lz -= (long)Value;

   return *this;
  }
Point3d &Point3d::operator += ( double V )
  {lx += (long)V;
   ly += (long)V;
   lz += (long)V;

   return *this;
  }
Point3d &Point3d::operator *= ( double Value )
  {lx *= (long)Value;
   ly *= (long)Value;
   lz *= (long)Value;
   return *this;
  }
Point3d &Point3d::operator /= ( double V )
  {lx /= (long)V;
   ly /= (long)V;
   lz /= (long)V;
   return *this;
  }



Point3d Object3d::getNormal(Point3d a, Point3d b, Point3d c) {
  Point3d n;
  Point3d U = (b - a);
  Point3d V = (c - a);
  n.lx = (U.ly * V.lz) - (U.lz * V.ly);
  n.ly = (U.lz * V.lx) - (U.lx * V.lz);
  n.lz = (U.lx * V.ly) - (U.ly * V.lx);

 // return Normalize(n);
  return n;
}

Point3d Object3d::Normalize(Point3d V) { 
  Point3d N;
 float length = sqrt((V.lx * V.lx) + (V.ly * V.ly) + (V.lx * V.lz));
  N.lx = V.lx/length;
  N.ly = V.ly/length;
  N.lz = V.lz/length;
  return N;
}

void Object3d::sortDepthMap() {
  int swapped;
  int i;
  int tempID;
  float tempDepth;
  for(int i=0; i<nFaces; i++) { 
    depthMap[i].depth = (mesh[faces[i].a].az + mesh[faces[i].b].az  + mesh[faces[i].c].az )/3;
    depthMap[i].ID = i;
  }
  do {
    swapped = 0;
    for (i = 1; i < this->nFaces; i++) {
      if (depthMap[i-1].depth > depthMap[i].depth) {
        tempID = depthMap[i].ID;
        tempDepth = depthMap[i].depth;
        depthMap[i].ID = depthMap[i-1].ID;
        depthMap[i].depth = depthMap[i-1].depth;
        depthMap[i-1].ID = tempID;
        depthMap[i-1].depth = tempDepth;
        swapped = 1;
      }
    }
  } while(swapped != 0);

  minDepth = depthMap[0].depth;
  maxDepth = depthMap[nFaces -1].depth;
}

void Object3d::loadMesh(float verts[][3], int faceArray[][3]) { 
  // load vertices into internal array
  for(int i = 0; i<this->nVertices;i++) { 
    mesh[i].lx = verts[i][0];
    mesh[i].ly = verts[i][1];
    mesh[i].lz = verts[i][2];
  }
  // load faces into internal array
  for(int i=0; i<this->nFaces; i++) { 
    faces[i].a = faceArray[i][0];
    faces[i].b = faceArray[i][1];
    faces[i].c = faceArray[i][2];
  }
}

Matrix3d::Matrix3d() 
{
  MatrixIdentity();
}

void Matrix3d::MatrixReset() { 
  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
      Matrix[i][j]=0;
    }
  }
}

void Matrix3d::MatrixIdentity() { 
  MatrixReset();
  Matrix[0][0] = Matrix[1][1] = Matrix [2][2] = Matrix[3][3] = 1;
}

void Matrix3d::MatrixCopy(Matrix3d &NewM){
  Matrix3d temp;
  int i,j;
  for(i = 0;i<4;i++){
    for(j =0; j<4; j++) {
      temp.Matrix[i][j] = (Matrix[i][0]*NewM.Matrix[0][j]) + (Matrix[i][1]*NewM.Matrix[1][j]) + (Matrix[i][2]*NewM.Matrix[2][j]) + (Matrix[i][3]*NewM.Matrix[3][j]);
    }
  }
  for(i=0;i<4;i++) { 
    for(j=0;j<4;j++) {
      Matrix[i][j] = temp.Matrix[i][j];
    }
  }
}

void Matrix3d::MatrixMult(Matrix3d &M1, Matrix3d &M2){
  Matrix3d temp;
  int i,j;
  for(i = 0;i<4;i++){
    for(j =0; j<4; j++) {
      temp.Matrix[i][j] = (M2.Matrix[i][0]*M1.Matrix[0][j]) + (M2.Matrix[i][1]*M1.Matrix[1][j]) + (M2.Matrix[i][2]*M1.Matrix[2][j]) + (M2.Matrix[i][3]*M1.Matrix[3][j]);
    }
  }
  for(i=0;i<4;i++) { 
    for(j=0;j<4;j++) {
      M1.Matrix[i][j] = temp.Matrix[i][j];
    }
  }
}

Object3d::Object3d(int numVerts, int numFaces) { 
  mesh = new Point3d[numVerts];
  depthMap = new DepthMap[numFaces];
  faces = new FaceList[numFaces];
  this->nVertices = numVerts;
  this->nFaces = numFaces;
  init();
  local = 1;
}

void Object3d::init() { 
  matrix.MatrixIdentity();
  objectmatrix.MatrixIdentity();
}

void Object3d::Translate(float x, float y, float z) {
  Rmat.MatrixIdentity();
  Rmat.Matrix[3][0]=x;
  Rmat.Matrix[3][1]=y;
  Rmat.Matrix[3][2]=z;
  if(local) { 
    objectmatrix.MatrixCopy(Rmat);
  } else {
    matrix.MatrixCopy(Rmat);
  }
}

void Object3d::Rotate(float x, float y, float z) { 
  rmatrix.MatrixIdentity();
  Rmat.MatrixIdentity();
  Rmat.Matrix[1][1]=cos(x); Rmat.Matrix[1][2]=sin(x);
  Rmat.Matrix[2][1]=-(sin(x)); Rmat.Matrix[2][2]=cos(x);
  rmatrix.MatrixMult(rmatrix,Rmat);
  Rmat.MatrixIdentity();
  Rmat.Matrix[0][0]=cos(y);Rmat.Matrix[0][2]=-(sin(y));
  Rmat.Matrix[2][0]=sin(y);Rmat.Matrix[2][2]=cos(y);
  Rmat.MatrixMult(rmatrix,Rmat);
  Rmat.MatrixIdentity();
  Rmat.Matrix[0][0]=cos(z); Rmat.Matrix[0][1]=sin(z);
  Rmat.Matrix[1][0]=-(sin(z)); Rmat.Matrix[1][1]=cos(z);
  Rmat.MatrixMult(rmatrix,Rmat);

  if(local)
   {
    objectmatrix.MatrixIdentity();
    objectmatrix.MatrixCopy(rmatrix);
   }
  else
   {
    matrix.MatrixCopy(rmatrix);
   } 
}

void Object3d::Scale(float scale) {
  Rmat.MatrixIdentity();
  Rmat.Matrix[0][0] = scale;
  Rmat.Matrix[1][1] = scale;
  Rmat.Matrix[2][2] = scale;
  if(local) { 
    objectmatrix.MatrixCopy(Rmat);
  } else { 
    matrix.MatrixCopy(Rmat);
  }
}

Point3d Object3d::ChangeLocalObject(Point3d &p)
{ 
  p.wx=(long)(p.ax*matrix.Matrix[0][0]+p.ay*matrix.Matrix[1][0]+p.az*matrix.Matrix[2][0]+matrix.Matrix[3][0]);
  p.wy=(long)(p.ax*matrix.Matrix[0][1]+p.ay*matrix.Matrix[1][1]+p.az*matrix.Matrix[2][1]+matrix.Matrix[3][1]);
  p.wz=(long)(p.ax*matrix.Matrix[0][2]+p.ay*matrix.Matrix[1][2]+p.az*matrix.Matrix[2][2]+matrix.Matrix[3][2]);
  return p;
}

Point3d Object3d::ChangeObjectPoint(Point3d &p)
{
  p.ax=(long)(p.lx*objectmatrix.Matrix[0][0]+p.ly*objectmatrix.Matrix[1][0]+(long)p.lz*objectmatrix.Matrix[2][0]+objectmatrix.Matrix[3][0]);
  p.ay=(long)(p.lx*objectmatrix.Matrix[0][1]+p.ly*objectmatrix.Matrix[1][1]+(long)p.lz*objectmatrix.Matrix[2][1]+objectmatrix.Matrix[3][1]);
  p.az=(long)(p.lx*objectmatrix.Matrix[0][2]+p.ly*objectmatrix.Matrix[1][2]+(long)p.lz*objectmatrix.Matrix[2][2]+objectmatrix.Matrix[3][2]);
  return p;
}
int Object3d::Render(float Center, float FOV, int faceID, char *dir) {
  int fid = depthMap[faceID].ID;
  if(dir == "ax") {
    return Center + FOV * mesh[faces[fid].a].wx / mesh[faces[fid].a].wz;
  } else if (dir == "ay"){ 
    return Center - FOV * mesh[faces[fid].a].wy / mesh[faces[fid].a].wz;
  } else if (dir == "bx") {
    return Center + FOV * mesh[faces[fid].b].wx / mesh[faces[fid].b].wz;
  } else if (dir == "by") {
    return Center - FOV * mesh[faces[fid].b].wy / mesh[faces[fid].b].wz;
  } else if (dir == "cx") {
    return Center + FOV * mesh[faces[fid].c].wx / mesh[faces[fid].c].wz;
  } else { 
    return Center - FOV * mesh[faces[fid].c].wy / mesh[faces[fid].c].wz;
  }
}

void Object3d::ApplyTransforms() { 
  for(int i=0; i<nVertices;i++) {
    ChangeObjectPoint(mesh[i]);
    ChangeLocalObject(mesh[i]);
  }
}