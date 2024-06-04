//#include "StdAfx.h"
//#include ".\simplexvol.h"
#include "math.h"
#include "float.h"

#include <vtkTriangle.h>
#include <vtkPointData.h> 
#include <vtkCellData.h> 
#include <vtkPolygon.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkPointSet.h> 
#include <vtkPolyDataWriter.h> 
#include <vtkPolyDataNormals.h> 
#include <vtkPolyDataReader.h>
#include <vtkJPEGWriter.h>
#include <vtkSTLWriter.h>
#include <vtkImageCast.h>

#include <vtkSmartPointer.h>

#include "deformation.h"


const bool USE_ISOTROPICDAMPING=false;


CSimplexMesh::CSimplexMesh(void)
{
TimeCounter=0;
Monitor=false; Monitor_edgelenght=false; Monitor_density=false;  Monitor_displacement=false;
gamma=0;timestep=0.5;
mass_modified=false;    mass_scalar=1; 
Flipped=false;
}

CSimplexMesh::~CSimplexMesh(void)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void CSimplexMesh::GetPrincipalAxis(P_float pm[3],P_float axis[3][3])	// get OBB centred on center of gravity
{
GetPrincipalDirections(pm,axis,nb_points,points,false);
}



// TO MODIFY !!!!
void CSimplexMesh::GetPrincipalAxis(const P_float p1[3],const P_float p2[3],const P_float pm[3],P_float axis[3][3]) // get OBB with 2 constraint points
{
int i;
P_float u[3],v[3],w[3],x,y,dp;

// principal axis
	//1st axis
axis[0][0]=p2[0]-p1[0]; axis[0][1]=p2[1]-p1[1]; axis[0][2]=p2[2]-p1[2];
dp=norm(axis[0]); axis[0][0]=axis[0][0]/dp; axis[0][1]=axis[0][1]/dp; axis[0][2]=axis[0][2]/dp;

	//2d axis
v[0]=1-p1[0];v[1]=-p1[1];v[2]=-p1[2];
dp=dotproduct(v,axis[0]); v[0]=v[0]-dp*axis[0][0]; v[1]=v[1]-dp*axis[0][1]; v[2]=v[2]-dp*axis[0][2];
dp=norm(v); v[0]=v[0]/dp; v[1]=v[1]/dp; v[2]=v[2]/dp;
crossproduct(w,axis[0],v);

x=0; y=0;
for(i=0;i<nb_points;i++)  
	{
	u[0]=GetPoint(i)[0]-p1[0];u[1]=GetPoint(i)[1]-p1[1];u[2]=GetPoint(i)[2]-p1[2];
	dp=dotproduct(u,axis[0]); u[0]=u[0]-dp*axis[0][0]; u[1]=u[1]-dp*axis[0][1]; u[2]=u[2]-dp*axis[0][2];
	dp=dotproduct(u,v); x+=dp*dp;
	dp=dotproduct(u,w); y+=dp*dp;
	}
axis[1][0]=x*v[0]+y*w[0]; axis[1][1]=x*v[1]+y*w[1]; axis[1][2]=x*v[2]+y*w[2];
dp=norm(axis[1]); axis[1][0]=axis[1][0]/dp; axis[1][1]=axis[1][1]/dp; axis[1][2]=axis[1][2]/dp;

	//3rd axis
crossproduct(axis[2],axis[0],axis[1]);

// distances max along axis
P_float dmax=0;
for(i=0;i<nb_points;i++)  {u[0]=GetPoint(i)[0]-pm[0];u[1]=GetPoint(i)[1]-pm[1];u[2]=GetPoint(i)[2]-pm[2];	dp=abs(dotproduct(u,axis[0])); if(dp>dmax) dmax=dp; }
axis[0][0]=2*axis[0][0]*dmax; axis[0][1]=2*axis[0][1]*dmax; axis[0][2]=2*axis[0][2]*dmax;
dmax=0;
for(i=0;i<nb_points;i++)  {u[0]=GetPoint(i)[0]-pm[0];u[1]=GetPoint(i)[1]-pm[1];u[2]=GetPoint(i)[2]-pm[2];	dp=abs(dotproduct(u,axis[1])); if(dp>dmax) dmax=dp; }
axis[1][0]=2*axis[1][0]*dmax; axis[1][1]=2*axis[1][1]*dmax; axis[1][2]=2*axis[1][2]*dmax;
dmax=0;
for(i=0;i<nb_points;i++)  {u[0]=GetPoint(i)[0]-pm[0];u[1]=GetPoint(i)[1]-pm[1];u[2]=GetPoint(i)[2]-pm[2];	dp=abs(dotproduct(u,axis[2])); if(dp>dmax) dmax=dp; }
axis[2][0]=2*axis[2][0]*dmax; axis[2][1]=2*axis[2][1]*dmax; axis[2][2]=2*axis[2][2]*dmax;
}

// TO MODIFY !!!!
void CSimplexMesh::GetPrincipalAxis(const P_float p1[3],const P_float p2[3],const P_float p3[3],const P_float pm[3],P_float axis[3][3]) // get OBB with 3 constraint points
{
int i;
P_float u[3],dp;

// principal axis
axis[0][0]=p2[0]-p1[0]; axis[0][1]=p2[1]-p1[1]; axis[0][2]=p2[2]-p1[2];
dp=norm(axis[0]); axis[0][0]=axis[0][0]/dp; axis[0][1]=axis[0][1]/dp; axis[0][2]=axis[0][2]/dp;

axis[1][0]=0; axis[1][1]=0; axis[1][2]=0;
axis[1][0]=p3[0]-p1[0];axis[1][1]=p3[1]-p1[1];axis[1][2]=p3[2]-p1[2];
dp=dotproduct(axis[1],axis[0]); axis[1][0]=axis[1][0]-dp*axis[0][0]; axis[1][1]=axis[1][1]-dp*axis[0][1]; axis[1][2]=axis[1][2]-dp*axis[0][2];
dp=norm(axis[1]); axis[1][0]=axis[1][0]/dp; axis[1][1]=axis[1][1]/dp; axis[1][2]=axis[1][2]/dp;
crossproduct(axis[2],axis[0],axis[1]);

// distances max along axis
P_float dmax=0;
for(i=0;i<nb_points;i++)  {u[0]=GetPoint(i)[0]-pm[0];u[1]=GetPoint(i)[1]-pm[1];u[2]=GetPoint(i)[2]-pm[2];	dp=abs(dotproduct(u,axis[0])); if(dp>dmax) dmax=dp; }
axis[0][0]=2*axis[0][0]*dmax; axis[0][1]=2*axis[0][1]*dmax; axis[0][2]=2*axis[0][2]*dmax;
dmax=0;
for(i=0;i<nb_points;i++)  {u[0]=GetPoint(i)[0]-pm[0];u[1]=GetPoint(i)[1]-pm[1];u[2]=GetPoint(i)[2]-pm[2];	dp=abs(dotproduct(u,axis[1])); if(dp>dmax) dmax=dp; }
axis[1][0]=2*axis[1][0]*dmax; axis[1][1]=2*axis[1][1]*dmax; axis[1][2]=2*axis[1][2]*dmax;
dmax=0;
for(i=0;i<nb_points;i++)  {u[0]=GetPoint(i)[0]-pm[0];u[1]=GetPoint(i)[1]-pm[1];u[2]=GetPoint(i)[2]-pm[2];	dp=abs(dotproduct(u,axis[2])); if(dp>dmax) dmax=dp; }
axis[2][0]=2*axis[2][0]*dmax; axis[2][1]=2*axis[2][1]*dmax; axis[2][2]=2*axis[2][2]*dmax;
}

void CSimplexMesh::RigidTransform(P_float M[16])
{
P_float p[3];
for(int i=0;i<nb_points;i++) {GetPoint(i,p); Transform(p,points+3*i,M);}
Equilibrium();
}


void CSimplexMesh::GetBoundingBox(P_float bounds[3][2])
{
P_float p[3];
bounds[0][0]=points[0]; bounds[0][1]=points[0]; bounds[1][0]=points[1]; bounds[1][1]=points[1]; bounds[2][0]=points[2]; bounds[2][1]=points[2];
for(int i=1;i<nb_points;i++)
    {
    GetPoint(i,p);
    if(p[0]<bounds[0][0]) bounds[0][0]=p[0];    if(p[0]>bounds[0][1]) bounds[0][1]=p[0];
    if(p[1]<bounds[1][0]) bounds[1][0]=p[1];    if(p[1]>bounds[1][1]) bounds[1][1]=p[1];
    if(p[2]<bounds[2][0]) bounds[2][0]=p[2];    if(p[2]>bounds[2][1]) bounds[2][1]=p[2];
    }
}

void CSimplexMesh::GetBoundingDOP(P_float bounds[9][2])
{
P_float p[3],val;
bounds[0][0]=points[0]; bounds[0][1]=points[0];  // +-(1,0,0)
bounds[1][0]=points[1]; bounds[1][1]=points[1];  // +-(0,1,0)
bounds[2][0]=points[2]; bounds[2][1]=points[2]; // +-(0,0,1)
bounds[3][0]=points[0]+points[1]; bounds[3][1]=points[0]+points[1]; // +-(1,1,0)
bounds[4][0]=points[0]-points[1]; bounds[4][1]=points[0]-points[1]; // +-(1,-1,0)
bounds[5][0]=points[0]+points[2]; bounds[5][1]=points[0]+points[2]; // +-(1,0,1)
bounds[6][0]=points[0]-points[1]; bounds[6][1]=points[0]-points[2]; // +-(1,0,-1)
bounds[7][0]=points[1]+points[2]; bounds[7][1]=points[1]+points[2]; // +-(0,1,1)
bounds[8][0]=points[1]-points[2]; bounds[8][1]=points[1]-points[2]; // +-(0,1,-1)

for(int i=1;i<nb_points;i++)
    {
    GetPoint(i,p);
    if(p[0]<bounds[0][0]) bounds[0][0]=p[0];    if(p[0]>bounds[0][1]) bounds[0][1]=p[0];
    if(p[1]<bounds[1][0]) bounds[1][0]=p[1];    if(p[1]>bounds[1][1]) bounds[1][1]=p[1];
    if(p[2]<bounds[2][0]) bounds[2][0]=p[2];    if(p[2]>bounds[2][1]) bounds[2][1]=p[2];

    val=p[0]+p[1]; if(val<bounds[3][0]) bounds[3][0]=val;   if(p[0]+p[1]>bounds[3][1]) bounds[3][1]=val;
    val=p[0]-p[1]; if(val<bounds[4][0]) bounds[4][0]=val;   if(p[0]+p[1]>bounds[4][1]) bounds[4][1]=val;
    val=p[0]+p[2]; if(val<bounds[5][0]) bounds[5][0]=val;   if(p[0]+p[1]>bounds[5][1]) bounds[5][1]=val;
    val=p[0]-p[2]; if(val<bounds[6][0]) bounds[6][0]=val;   if(p[0]+p[1]>bounds[6][1]) bounds[6][1]=val;
    val=p[1]+p[2]; if(val<bounds[7][0]) bounds[7][0]=val;   if(p[0]+p[1]>bounds[7][1]) bounds[7][1]=val;
    val=p[1]-p[2]; if(val<bounds[8][0]) bounds[8][0]=val;   if(p[0]+p[1]>bounds[8][1]) bounds[8][1]=val;
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::Free() {}
int CSimplexMesh::GetNumberOfPoints() {return nb_points; }
int CSimplexMesh::GetNumberOfCells() {return nb_cells;}
void CSimplexMesh::Copy(CSimplexMesh* mesh) {}

/////////////////////////////////////////////////////////////////////////////////////////////////
CSimplexMesh* CSimplexMesh::IncreaseResolution() {return NULL;} // virtual
CSimplexMesh* CSimplexMesh::DecreaseResolution() {return NULL;} // virtual
void CSimplexMesh::UpdatePointsFromLowerRes(CSimplexMesh* mesh) {} // virtual
void CSimplexMesh::UpdatePointsFromHigherRes(CSimplexMesh* mesh) {} // virtual
void CSimplexMesh::UpdateValsFromHigherRes(CSimplexMesh* mesh) {} // virtual
void CSimplexMesh::UpdateValsFromLowerRes(CSimplexMesh* mesh) {} // virtual
void CSimplexMesh::UpdateForcesFromLowerRes(CSimplexMesh* mesh,bool updatederivatives) {} // virtual
void CSimplexMesh::UpdateVal_c() {}// virtual
/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::SetPoint(int index,P_float x,P_float y, P_float z) {points[3*index]=x;points[3*index+1]=y;points[3*index+2]=z;}
void CSimplexMesh::SetPoint(int index,P_float p[3]) {memcpy(points+3*index,p,3*sizeof(P_float));}
void CSimplexMesh::SetPoint_tm1(int index,P_float x,P_float y, P_float z) {points_tm1[3*index]=x;points_tm1[3*index+1]=y;points_tm1[3*index+2]=z;}
void CSimplexMesh::SetPoint_tm1(int index,P_float p[3]) {memcpy(points_tm1+3*index,p,3*sizeof(P_float));}
void CSimplexMesh::GetPoint(int index,P_float p[3]) {memcpy(p,points+3*index,3*sizeof(P_float));}
P_float* CSimplexMesh::GetPoints() {return points;}
P_float* CSimplexMesh::GetPoint(int index) {return points+3*index;}
void CSimplexMesh::GetPoint_tm1(int index,P_float p[3]) {memcpy(p,points_tm1+3*index,3*sizeof(P_float));}
P_float* CSimplexMesh::GetPoints_tm1() {return points_tm1;}
P_float* CSimplexMesh::GetPoint_tm1(int index) {return points_tm1+3*index;}
void CSimplexMesh::SetSpeed(int index,P_float x,P_float y, P_float z) {speeds[3*index]=x;speeds[3*index+1]=y;speeds[3*index+2]=z;}
void CSimplexMesh::SetSpeed(int index,P_float p[3]) {memcpy(speeds+3*index,p,3*sizeof(P_float));}
void CSimplexMesh::SetSpeed_tm1(int index,P_float x,P_float y, P_float z) {speeds_tm1[3*index]=x;speeds_tm1[3*index+1]=y;speeds_tm1[3*index+2]=z;}
void CSimplexMesh::SetSpeed_tm1(int index,P_float p[3]) {memcpy(speeds_tm1+3*index,p,3*sizeof(P_float));}
void CSimplexMesh::GetSpeed(int index,P_float p[3]) {memcpy(p,speeds+3*index,3*sizeof(P_float));}
P_float* CSimplexMesh::GetSpeeds() {return speeds;}
P_float* CSimplexMesh::GetSpeed(int index) {return speeds+3*index;}
void CSimplexMesh::GetSpeed_tm1(int index,P_float p[3]) {memcpy(p,speeds_tm1+3*index,3*sizeof(P_float));}
P_float* CSimplexMesh::GetSpeeds_tm1() {return speeds_tm1;}
P_float* CSimplexMesh::GetSpeed_tm1(int index) {return speeds_tm1+3*index;}
void CSimplexMesh::GetForce(int index,P_float f[3]) {memcpy(f,forces+3*index,3*sizeof(P_float));}
void CSimplexMesh::SetForce(int index,P_float f[3]) {memcpy(forces+3*index,f,3*sizeof(P_float));}
void CSimplexMesh::SetForce(int index,P_float x,P_float y, P_float z) {forces[3*index]=x;forces[3*index+1]=y;forces[3*index+2]=z;}
void CSimplexMesh::AddForce(int index,P_float f[3]) {forces[3*index]+=f[0];forces[3*index+1]+=f[1];forces[3*index+2]+=f[2];}
P_float* CSimplexMesh::GetForces()  {return forces;} 
P_float* CSimplexMesh::GetForce(int index)  {return forces+3*index;}
P_float* CSimplexMesh::GetDfs() {return df_s;}
P_float* CSimplexMesh::GetDfs(int index) {return df_s+6*index;}
void CSimplexMesh::GetDfs(int index,P_float dfs[6]) {memcpy(dfs,df_s+6*index,6*sizeof(P_float));}
void CSimplexMesh::AddDfs(int index,P_float dfs[6]) {for(int i=0;i<6;i++) df_s[6*index+i]+=dfs[i];}
P_float* CSimplexMesh::GetDfp() {return df_p;}
P_float* CSimplexMesh::GetDfp(int index) {return df_p+6*index;}
void CSimplexMesh::GetDfp(int index,P_float dfp[6]) {memcpy(dfp,df_p+6*index,6*sizeof(P_float));}
void CSimplexMesh::AddDfp(int index,P_float dfp[6]) {for(int i=0;i<6;i++) df_p[6*index+i]+=dfp[i];}

void CSimplexMesh::SetVal_p(const int index,const P_float v[3]) {memcpy(Val_p+3*index,v,3*sizeof(P_float));}
void CSimplexMesh::SetVal_c(const int index,const P_float v[3]) {memcpy(Val_c+3*index,v,3*sizeof(P_float));}
P_float* CSimplexMesh::GetVal_p(const int index) {return Val_p+3*index;}
P_float* CSimplexMesh::GetVal_c(const int index) {return Val_c+3*index;}
void CSimplexMesh::FillVal0() {for(int i=0;i<3*nb_points;i++) *(Val_p+i)=0; for(int i=0;i<3*nb_cells;i++) *(Val_c+i)=0;}

P_float* CSimplexMesh::GetMass_inv(int index) {return mass_inv+3*index;}
P_float* CSimplexMesh::GetMass(int index) {return mass+3*index;}
P_float* CSimplexMesh::GetMass_constraint(int index) {return mass_constraint[index];}
bool CSimplexMesh::GetMass_Modified() {return mass_modified;}

void CSimplexMesh::SetMass(P_float m) {if(m==-1) SwicthOffAllForces(); for(int i=0;i<nb_points;i++) SetMass(i,m);}
void CSimplexMesh::SetMass() {for(int i=0;i<nb_points;i++) SetMass(i);}

void CSimplexMesh::SetMass(int index,P_float m) 
{
// Mi^-1=1/mi.I
if(m==-1) {mass_inv[3*index]=0; mass_inv[3*index+1]=0; mass_inv[3*index+2]=0; mass[3*index]=-1; mass[3*index+1]=-1; mass[3*index+2]=-1;}
else {mass_inv[3*index]=1/m; mass_inv[3*index+1]=1/m; mass_inv[3*index+2]=1/m; mass[3*index]=m; mass[3*index+1]=m; mass[3*index+2]=m;}
SetSpeed(index,0,0,0);
SetSpeed_tm1(index,0,0,0);
SetPoint_tm1(index,GetPoint(index));
}

void CSimplexMesh::SetMass(int index,P_float mx,P_float my,P_float mz) 
{
// Mi^-1=1/mi.I
if(mx==-1) {mass_inv[3*index]=0; mass[3*index]=-1;} else if(mx!=0) {mass_inv[3*index]=1/mx; mass[3*index]=mx;}
if(my==-1) {mass_inv[3*index+1]=0; mass[3*index+1]=-1;}  else if(my!=0) {mass_inv[3*index+1]=1/my; mass[3*index+1]=my;}
if(mz==-1) {mass_inv[3*index+2]=0; mass[3*index+2]=-1;}  else if(mz!=0) {mass_inv[3*index+2]=1/mz; mass[3*index+2]=mz;}
SetSpeed(index,0,0,0);
SetSpeed_tm1(index,0,0,0);
SetPoint_tm1(index,GetPoint(index));
mass_modified=true;
}


void CSimplexMesh::SetMass(int index) 
{
// Mi^-1=1/mi.I
mass_inv[3*index]=1/mass_scalar;            mass_inv[3*index+1]=1/mass_scalar;          mass_inv[3*index+2]=1/mass_scalar; 
mass[3*index]=mass_scalar;					mass[3*index+1]=mass_scalar;				mass[3*index+2]=mass_scalar;
SetSpeed(index,0,0,0);
SetSpeed_tm1(index,0,0,0);
SetPoint_tm1(index,GetPoint(index));
}


void CSimplexMesh::SetMass_planeconstraint(int index,P_float n[3])
{
// Mi^-1=1/mi.(I-N.NT)
P_float nrm=norm(n);
P_float u[3]={n[0]/nrm,n[1]/nrm,n[2]/nrm};
if(mass_constraint[index]==NULL) free(mass_constraint[index]);  mass_constraint[index]=new P_float[6];
mass_constraint[index][0]=(1-u[0]*u[0]);    
mass_constraint[index][1]=(-u[0]*u[1]);   mass_constraint[index][2]=(1-u[1]*u[1]);  
mass_constraint[index][3]=(-u[0]*u[2]);   mass_constraint[index][4]=(-u[1]*u[2]);   mass_constraint[index][5]=(1-u[2]*u[2]);
mass_modified=true;

P_float dp=dotproduct(n,speeds+3*index);
speeds[3*index]-=dp*n[0]; speeds[3*index+1]-=dp*n[1]; speeds[3*index+2]-=dp*n[2];
dp=dotproduct(n,speeds_tm1+3*index);
speeds_tm1[3*index]-=dp*n[0]; speeds_tm1[3*index+1]-=dp*n[1]; speeds_tm1[3*index+2]-=dp*n[2];
SetPoint_tm1(index,GetPoint(index));
}


void CSimplexMesh::SetMass_vectorconstraint(int index,P_float v[3])
{
// Mi^-1=1/mi.N.NT
P_float nrm=norm(v);
P_float u[3]={v[0]/nrm,v[1]/nrm,v[2]/nrm};
if(mass_constraint[index]==NULL) free(mass_constraint[index]);  mass_constraint[index]=new P_float[6];
mass_constraint[index][0]=(u[0]*u[0]);  
mass_constraint[index][1]=(u[0]*u[1]);    mass_constraint[index][2]=(u[1]*u[1]);    
mass_constraint[index][3]=(u[0]*u[2]);    mass_constraint[index][4]=(u[1]*u[2]);    mass_constraint[index][5]=(u[2]*u[2]);
mass_modified=true;

P_float dp=dotproduct(v,speeds+3*index);
speeds[3*index]=dp*v[0]; speeds[3*index+1]=dp*v[1]; speeds[3*index+2]=dp*v[2];
dp=dotproduct(v,speeds_tm1+3*index);
speeds_tm1[3*index]=dp*v[0]; speeds_tm1[3*index+1]=dp*v[1]; speeds_tm1[3*index+2]=dp*v[2];
SetPoint_tm1(index,GetPoint(index));
}

void CSimplexMesh::UpdateMass() {} // virtual
void CSimplexMesh::SwicthOffAllForces() {} // virtual



/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::SetCell(int index,int* f) {if(cells[index]!=NULL) free(cells[index]); cells[index]=new int[f[0]+1]; memcpy(cells[index],f,(f[0]+1)*sizeof(int));}
void CSimplexMesh::SetCell(int index,int nb_p, int* p)
{
if(cells[index]!=NULL) free(cells[index]);
cells[index]=new int[nb_p+1];
for(int i=0;i<nb_p;i++) cells[index][i+1]=p[i];
cells[index][0]=nb_p;
}


int CSimplexMesh::GetCell(int index,int* f) {memcpy(f,cells[index]+1,cells[index][0]*sizeof(int)); return cells[index][0];}
int* CSimplexMesh::GetCell(int index) {return cells[index];}


int** CSimplexMesh::GetCells() {return cells;}
void CSimplexMesh::DeleteCell(int index) {} // virtual
void CSimplexMesh::InsertNextCell(int* f) {}  // virtual
/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::SaveMesh(const char* filename) {}
void CSimplexMesh::SaveParams(const char* filename) {}
void CSimplexMesh::LoadParams(const char* filename) {}
void CSimplexMesh::LoadMesh(const char* filename) {}

void CSimplexMesh::SaveVTKMesh(const char* filename,int lift)
{
vtkPolyData* polydata=GetPolyData(lift);

vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
	normals->SetInputData(polydata);
	normals->SplittingOff();
	normals->ConsistencyOn();
	normals->ComputePointNormalsOn();
	if(Flipped) normals->FlipNormalsOn();

	vtkPolyDataWriter* writer=vtkPolyDataWriter::New();
	writer->SetInputData(normals->GetOutput());
    writer->SetFileName(filename);
    writer->SetFileType(1);
    writer->Write();

writer->Delete();
polydata->Delete();
normals->Delete();
}

void CSimplexMesh::SaveSTLMesh(const char* filename,int lift)
{
vtkPolyData* polydata=GetPolyData(lift);

vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
	normals->SetInputData(polydata);
	normals->SplittingOff();
	normals->ConsistencyOn();
	normals->ComputePointNormalsOn();
	if(Flipped) normals->FlipNormalsOn();

vtkSTLWriter* writer=vtkSTLWriter::New();
	writer->SetInputData(normals->GetOutput());
    writer->SetFileName(filename);
    writer->SetFileType(1);
    writer->Write();

writer->Delete();
polydata->Delete();
normals->Delete();
}


vtkUnstructuredGrid* CSimplexMesh::GetOBB()
{
P_float pm[3],axis[3][3],l[3]; GetPrincipalDirections(pm,axis,GetNumberOfPoints(),GetPoints(),false);
l[0]=norm(axis[0]); l[1]=norm(axis[1]); l[2]=norm(axis[2]); 

vtkUnstructuredGrid* UGrid=vtkUnstructuredGrid::New();  UGrid->Allocate();

vtkPoints* pts=vtkPoints::New();
    pts->SetNumberOfPoints(8);
	pts->SetPoint(0,pm[0]-axis[0][0]/2.-axis[1][0]/2.-axis[2][0]/2.,pm[1]-axis[0][1]/2.-axis[1][1]/2.-axis[2][1]/2.,pm[2]-axis[0][2]/2.-axis[1][2]/2.-axis[2][2]/2.);
	pts->SetPoint(1,pm[0]-axis[0][0]/2.+axis[1][0]/2.-axis[2][0]/2.,pm[1]-axis[0][1]/2.+axis[1][1]/2.-axis[2][1]/2.,pm[2]-axis[0][2]/2.+axis[1][2]/2.-axis[2][2]/2.);
	pts->SetPoint(2,pm[0]-axis[0][0]/2.-axis[1][0]/2.+axis[2][0]/2.,pm[1]-axis[0][1]/2.-axis[1][1]/2.+axis[2][1]/2.,pm[2]-axis[0][2]/2.-axis[1][2]/2.+axis[2][2]/2.);
	pts->SetPoint(3,pm[0]-axis[0][0]/2.+axis[1][0]/2.+axis[2][0]/2.,pm[1]-axis[0][1]/2.+axis[1][1]/2.+axis[2][1]/2.,pm[2]-axis[0][2]/2.+axis[1][2]/2.+axis[2][2]/2.);
	pts->SetPoint(4,pm[0]+axis[0][0]/2.-axis[1][0]/2.-axis[2][0]/2.,pm[1]+axis[0][1]/2.-axis[1][1]/2.-axis[2][1]/2.,pm[2]+axis[0][2]/2.-axis[1][2]/2.-axis[2][2]/2.);
	pts->SetPoint(5,pm[0]+axis[0][0]/2.+axis[1][0]/2.-axis[2][0]/2.,pm[1]+axis[0][1]/2.+axis[1][1]/2.-axis[2][1]/2.,pm[2]+axis[0][2]/2.+axis[1][2]/2.-axis[2][2]/2.);
	pts->SetPoint(6,pm[0]+axis[0][0]/2.-axis[1][0]/2.+axis[2][0]/2.,pm[1]+axis[0][1]/2.-axis[1][1]/2.+axis[2][1]/2.,pm[2]+axis[0][2]/2.-axis[1][2]/2.+axis[2][2]/2.);
	pts->SetPoint(7,pm[0]+axis[0][0]/2.+axis[1][0]/2.+axis[2][0]/2.,pm[1]+axis[0][1]/2.+axis[1][1]/2.+axis[2][1]/2.,pm[2]+axis[0][2]/2.+axis[1][2]/2.+axis[2][2]/2.);
UGrid->SetPoints(pts);

vtkPolygon* polygon=vtkPolygon::New(); polygon->GetPointIds()->SetNumberOfIds(4);
polygon->GetPointIds()->SetId(0,0); polygon->GetPointIds()->SetId(1,1); polygon->GetPointIds()->SetId(2,3); polygon->GetPointIds()->SetId(3,2); UGrid->InsertNextCell(polygon->GetCellType(),polygon->GetPointIds()); 
polygon->GetPointIds()->SetId(0,4); polygon->GetPointIds()->SetId(1,5); polygon->GetPointIds()->SetId(2,7); polygon->GetPointIds()->SetId(3,6); UGrid->InsertNextCell(polygon->GetCellType(),polygon->GetPointIds()); 
polygon->GetPointIds()->SetId(0,0); polygon->GetPointIds()->SetId(1,1); polygon->GetPointIds()->SetId(2,5); polygon->GetPointIds()->SetId(3,4); UGrid->InsertNextCell(polygon->GetCellType(),polygon->GetPointIds()); 
polygon->GetPointIds()->SetId(0,1); polygon->GetPointIds()->SetId(1,3); polygon->GetPointIds()->SetId(2,7); polygon->GetPointIds()->SetId(3,5); UGrid->InsertNextCell(polygon->GetCellType(),polygon->GetPointIds()); 
polygon->GetPointIds()->SetId(0,3); polygon->GetPointIds()->SetId(1,2); polygon->GetPointIds()->SetId(2,6); polygon->GetPointIds()->SetId(3,7); UGrid->InsertNextCell(polygon->GetCellType(),polygon->GetPointIds()); 
polygon->GetPointIds()->SetId(0,0); polygon->GetPointIds()->SetId(1,2); polygon->GetPointIds()->SetId(2,6); polygon->GetPointIds()->SetId(3,4); UGrid->InsertNextCell(polygon->GetCellType(),polygon->GetPointIds()); 
polygon->Delete();

UGrid->SetPoints(pts);
UGrid->Update();
return(UGrid);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::UpdateUnstructuredGrid() {} // virtual
void CSimplexMesh::UpdateUnstructuredGridCells() {} // virtual
vtkUnstructuredGrid* CSimplexMesh::GetUnstructuredGrid() {return UGrid;}
void CSimplexMesh::UpdateUnstructuredGridPoints()
{
vtkPoints* pts=vtkPoints::New();
    pts->SetNumberOfPoints(nb_points);
for(int i=0;i<nb_points;i++)  pts->SetPoint(i,GetPoint(i)); // copy points
UGrid->SetPoints(pts);
pts->Delete();
}

vtkPolyData* CSimplexMesh::GetPolyData(int lift) {return NULL;} // virtual

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::UpdateAll() {} // virtual
void CSimplexMesh::UpdateCellCenters() {} // virtual
P_float* CSimplexMesh::GetCellCenters() {return cellcenters;}
P_float* CSimplexMesh::GetCellCenter(int cell_index) {return cellcenters+3*cell_index;}
void CSimplexMesh::GetCellCenter(int cell_index,P_float c[3]) {memcpy(c,cellcenters+3*cell_index,3*sizeof(P_float));}
void CSimplexMesh::SetCellCenter(int cell_index,P_float c[3]) {memcpy(cellcenters+3*cell_index,c,3*sizeof(P_float));}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::UpdateParams() {} // virtual
void CSimplexMesh::UpdateParams(int index) {} // virtual
P_float* CSimplexMesh::GetParams(int index) {return NULL;} // virtual

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::SetGamma(P_float gam) {gamma=gam;}
void CSimplexMesh::SetTimeStep(P_float tstep) {timestep=tstep; timestep_inv2=1/(tstep*tstep);}
void CSimplexMesh::SetUseForces(bool useforces) {USEFORCES=useforces;} 
void CSimplexMesh::UpdateExternalForces() {} // virtual
void CSimplexMesh::UpdateInternalForces() {} // virtual
void CSimplexMesh::UpdateAlphas()  {} // virtual
void CSimplexMesh::TimeStep() {} // virtual
void CSimplexMesh::Equilibrium() {} // virtual

void CSimplexMesh::UpdateForcesDamping() {
	
	int i, k;
	P_float f[3], dp;

	if (USE_ISOTROPICDAMPING) {
		// iostropic damping
		for (i = 0; i < nb_points; i++) {
			f[0] = -gamma * speeds[3 * i]; 
			f[1] = -gamma * speeds[3 * i + 1]; 
			f[2] = -gamma * speeds[3 * i + 2];
			
			if (USEFORCES) {
				AddForce(i, f); 
				df_s[6 * i] = -gamma; 
				df_s[6 * i + 1] = 0; 
				df_s[6 * i + 2] = -gamma; 
				df_s[6 * i + 3] = 0; 
				df_s[6 * i + 4] = 0; 
				df_s[6 * i + 5] = -gamma;
			}
			else {
				AddDisp(i, f);
			}
		}
	}
	else {
		// anisotropic damping
		for (i = 0; i < nb_points; i++) {
			GetForce(i, f); 
			if (f[0] != 0 || f[1] != 0 || f[2] != 0) {
				dp = -gamma * dotproduct(speeds + 3 * i, f) / dotproduct(f, f); 
				f[0] = f[0] * dp; 
				f[1] = f[1] * dp; 
				f[2] = f[2] * dp;
			}

			if (USEFORCES) {
				AddForce(i, f); 
				GetDerivatives(1, f, df_s + 6 * i); 
				
				for (k = 0; k < 6; k++) {
					df_s[6 * i + k] = gamma * df_s[6 * i + k];
				}
			}
			else {
				AddDisp(i, f);
			}
		}
	}
}

void CSimplexMesh::UpdatePointsEnergyMinimisation()
{
P_float *p,*pm1;
for(int i=0; i<nb_points;i++)
 {
 SetSpeed_tm1(i,GetSpeed(i));
 pm1=GetPoint_tm1(i); p=GetPoint(i); 
 SetSpeed(i,(p[0]-pm1[0])/timestep,(p[1]-pm1[1])/timestep,(p[2]-pm1[2])/timestep);
 SetPoint_tm1(i,GetPoint(i));
 }
}

void CSimplexMesh::Add_indep(int index,P_float *PPref,P_float alpha,P_float extrapolation) {
	if(USEFORCES) {
		AddForce_indep(index,PPref,alpha,extrapolation); 
	}
	else {
		AddDisp_indep(index,PPref,alpha);
	}
}

void CSimplexMesh::AddForce_indep(int index,P_float *PPref,P_float alpha,P_float extrapolation)
{
int k;
P_float F[3]={alpha*PPref[0],alpha*PPref[1],alpha*PPref[2]};
P_float DFP[6]; GetDerivatives(extrapolation,PPref,DFP); for(k=0;k<6;k++) DFP[k]=alpha*DFP[k];

if(GetMass_constraint(index)!=NULL)	{} // to complete
else if(mass_modified) {} // to complete
else
	{
	for(k=0;k<3;k++) F[k]=mass[3*index]*F[k];
	for(k=0;k<6;k++) DFP[k]=mass[3*index]*DFP[k];
	}

for(k=0;k<3;k++) F[k]=timestep_inv2*F[k];
for(k=0;k<6;k++) DFP[k]=timestep_inv2*DFP[k];

AddForce(index,F);
AddDfp(index,DFP);
}

void CSimplexMesh::AddDisp_indep(int index,P_float *PPref,P_float alpha)
{
if(mass_inv[3*index]==0)  return;
P_float F[3]={alpha*PPref[0],alpha*PPref[1],alpha*PPref[2]};
AddDisp(index,F); 
}

void CSimplexMesh::AddDisp(int index,P_float *F) {} // virtual
/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexMesh::SetMonitor_edgelenght(bool val,const char* filename) {Monitor_edgelenght=val; Monitor_edgelenght_filename=filename; Update_Monitor();}
void CSimplexMesh::SetMonitor_density(bool val,const char* filename) {Monitor_density=val; Monitor_density_filename=filename; Update_Monitor();}
void CSimplexMesh::SetMonitor_displacement(bool val,const char* filename) {Monitor_displacement=val; Monitor_displacement_filename=filename; Update_Monitor();}
void CSimplexMesh::Update_Monitor() {} // virtual
void CSimplexMesh::Update_Monitor_edgelenght() {} // virtual
void CSimplexMesh::Update_Monitor_density() {} // virtual
void CSimplexMesh::Update_Monitor_displacement() // virtual
{
FILE* f=fopen(Monitor_displacement_filename.c_str(),"rt+");
if(f==NULL) f=fopen(Monitor_displacement_filename.c_str(),"wt"); 
else fseek(f,0,SEEK_END);
// MEAN / std dev
P_float mean=0,stddev=0,d;
int i;
for(i=0;i<nb_points;i++) mean+=dist3D(points+3*i,points_tm1+3*i);
mean=mean/(P_float)(nb_points);
for(i=0;i<nb_points;i++)  {d=dist3D(points+3*i,points_tm1+3*i); stddev+=(d-mean)*(d-mean);}
stddev=sqrt(stddev/(P_float)(nb_points));
fprintf(f,"%lf %lf\n",mean,stddev);
fclose(f);
}


/////////////////////////////////////////////////////////////////////////////////////////////////

void CSimplexMesh::SetNeighbors(int index,int i, int j, int k) {neighbors[3*index]=i; neighbors[3*index+1]=j; neighbors[3*index+2]=k;}

void CSimplexMesh::GetNeighbors(int index, int n[3]) {
	//memcpy(n,neighbors+3*index,3*sizeof(int));
		
	if (neighbors[3 * index + 0] != -2 && neighbors[3 * index + 1] != -2 && neighbors[3 * index + 2] != -2) {
		n[0] = neighbors[3 * index + 0];
		n[1] = neighbors[3 * index + 1];
		n[2] = neighbors[3 * index + 2];
	}
	else {
		//cout << "\n\nWARNING Point " << index << " is a multimaterial point and has two sets of neighbors.\nReturning only one set of neighbors.\nPlease consider using GetMultiMaterialNeighbors() with appropriate modifications.\n\n" << endl;
		int *mat = GetPointMaterialIndices(index);
		int *mmn = GetMultiMaterialNeighbors(mat[1], index);
		
		n[0] = mmn[0];
		n[1] = mmn[1];
		n[2] = mmn[2];
	}
}

int* CSimplexMesh::GetNeighbors(int index) {
	//return neighbors+3*index;
	if (neighbors[3 * index + 0] != -2 && neighbors[3 * index + 1] != -2 && neighbors[3 * index + 2] != -2) {
		int *ret = new int[3];
		ret[0] = neighbors[3 * index + 0];
		ret[1] = neighbors[3 * index + 1];
		ret[2] = neighbors[3 * index + 2];

		return ret;
	}
	else {
		//cout << "\n\nWARNING Point " << index << " is a multimaterial point and has two sets of neighbors.\nReturning only one set of neighbors.\nPlease consider using GetMultiMaterialNeighbors() with appropriate modifications.\n\n" << endl;
		int *mat = GetPointMaterialIndices(index);
		int *mmn = GetMultiMaterialNeighbors(mat[1], index);
		
		int *ret = new int[3];
		ret[0] = mmn[0];
		ret[1] = mmn[1];
		ret[2] = mmn[2];

		return ret;
	}
}



void CSimplexMesh::SetNeighbors(int index,int i,int j)
{
if(index<0 || i<0 || j<0) return;

if(neighbors[3*index]!=i && neighbors[3*index+1]!=i && neighbors[3*index+2]!=i)
    {
    if(neighbors[3*index]!=j && neighbors[3*index+1]!=j && neighbors[3*index+2]!=j)
        {
        // add both
        if(neighbors[3*index]==-1 && neighbors[3*index+1]==-1) {neighbors[3*index]=i; neighbors[3*index+1]=j;}
        if(neighbors[3*index+1]==-1 && neighbors[3*index+2]==-1) {neighbors[3*index+1]=i; neighbors[3*index+2]=j;}
        if(neighbors[3*index+2]==-1 && neighbors[3*index]==-1) {neighbors[3*index+2]=i; neighbors[3*index]=j;}
        }
    else
        {
        // add i
        if(neighbors[3*index]==j) neighbors[3*index+2]=i;
        if(neighbors[3*index+1]==j) neighbors[3*index]=i;
        if(neighbors[3*index+2]==j) neighbors[3*index+1]=i;
        }
    }
else
    {
    if(neighbors[3*index]!=j && neighbors[3*index+1]!=j && neighbors[3*index+2]!=j)
        {
        // add j
        if(neighbors[3*index]==i) neighbors[3*index+1]=j;
        if(neighbors[3*index+1]==i) neighbors[3*index+2]=j;
        if(neighbors[3*index+2]==i) neighbors[3*index]=j;
        }
    }
}


void CSimplexMesh::AllocatePointMaterialIndices() {
	//materialIndices = new int*[GetNumberOfPoints()];

	//// First two entries are material indices, and the last one is a scalar value. 
	//// Initialize with default material indices 0 and 2, and scalar value 5
	//// 0 and 2 are used for legacy purposes. 
	//// MM2M DC generates triangular meshes with material indices 0 and 2 for each triangular cell as default (single material mesh) and scalar value of 5. 
	//for (int i = 0; i < GetNumberOfPoints(); i++) {
	//	materialIndices[i] = new int[3];
	//	materialIndices[i][0] = 0;
	//	materialIndices[i][1] = 2;
	//	materialIndices[i][2] = 5;
	//}

	//hasMultiMaterialIndices = true;

	materialIndices = new int[3 * GetNumberOfPoints()]; 
	for (int i = 0; i < GetNumberOfPoints(); i++) {
		materialIndices[3 * i + 0] = 0;
		materialIndices[3 * i + 1] = 2;
		materialIndices[3 * i + 2] = 5;
	}
	hasMultiMaterialIndices = true;
}

int* CSimplexMesh::GetPointMaterialIndices(int index) {
	//int *ret = new int[3];
	//ret[0] = materialIndices[index][0];
	//ret[1] = materialIndices[index][1];
	//ret[2] = materialIndices[index][2];
	//return ret;

	int *ret = new int[3];
	ret[0] = materialIndices[3 * index + 0];
	ret[1] = materialIndices[3 * index + 1];
	ret[2] = materialIndices[3 * index + 2];
	return ret;
}

void CSimplexMesh::SetPointMaterialIndices(int index, int mat0, int mat1, int scalar) {
	materialIndices[3 * index + 0] = mat0;
	materialIndices[3 * index + 1] = mat1;
	materialIndices[3 * index + 2] = scalar;
}

bool CSimplexMesh::GetHasMultiMaterialIndices() {
	return hasMultiMaterialIndices;
}

void CSimplexMesh::SetHasMultiMaterialIndices(bool b) {
	hasMultiMaterialIndices = b;
}

/**
 * Initializing multiMaterial neighbors. 
 * 
 * 
 */
void CSimplexMesh::AllocateMultiMaterialNeighbors(int matLowerRange, int matUpperRange) {
	int numberOfMaterials = matUpperRange - matLowerRange + 1;
	multiMaterialNeighbors = new int*[numberOfMaterials];

	for (int i = 0; i < numberOfMaterials; i++) {
		multiMaterialNeighbors[i] = new int[GetNumberOfPoints() * 3];

		for (int j = 0; j < GetNumberOfPoints(); j++) {
			multiMaterialNeighbors[i][3 * j + 0] = -1;
			multiMaterialNeighbors[i][3 * j + 1] = -1;
			multiMaterialNeighbors[i][3 * j + 2] = -1;
		}
	}
}

void CSimplexMesh::SetMultiMaterialNeighbors(int matIndex, int ptIndex, int n0, int n1, int n2) {
	// Assign -2 in the regular neighbors array.
	// The -2 value indicates that the point is a multimaterial point 
	// and has neighbors in the multiMaterialNeighbors array. 
	SetNeighbors(ptIndex, -2, -2, -2); 

	// Set the points neighbors wrt to material indices.
	// Using [matIndex - 2] so that we can directly use the material
	// indices when calling this function. Since material index values
	// all start from 2, 3, 4, ...,
	// For example, when calling this function as
	// SetMultiMaterialNeighbors(2, x, n0, n1, n2)
	// it means setting the neighbors for point x for material 2. 
	// In the multimaterialNeighbors array, [matIndex - 2] index means
	// the neighbors for material index 2 will be stored in multiMaterialNeighbors[0][3 * x + 0], 
	// multiMaterialNeighbors[0][3 * x + 1] and multiMaterialNeighbors[0][3 * x + 2]
	multiMaterialNeighbors[matIndex - 2][3 * ptIndex + 0] = n0;
	multiMaterialNeighbors[matIndex - 2][3 * ptIndex + 1] = n1;
	multiMaterialNeighbors[matIndex - 2][3 * ptIndex + 2] = n2;
}

int* CSimplexMesh::GetMultiMaterialNeighbors(int matIndex, int ptIndex) {
	int *ret = new int[3];
	
	// Using [matIndex - 2] so that we can directly use the material
	// indices when calling this function. Since material index values
	// all start from 2, 3, 4, ...,
	// Same as in SetMultiMaterialNeighbors()
	ret[0] = multiMaterialNeighbors[matIndex - 2][3 * ptIndex + 0];
	ret[1] = multiMaterialNeighbors[matIndex - 2][3 * ptIndex + 1];
	ret[2] = multiMaterialNeighbors[matIndex - 2][3 * ptIndex + 2];

	return ret;
}

bool CSimplexMesh::IsMultiMaterialPoint(int index) {
	bool ret = isMultiMaterialPoints[index];
	return ret;
}

void CSimplexMesh::SetIsMultiMaterialPoint(int index, bool mm) {
	isMultiMaterialPoints[index] = mm;
}

void CSimplexMesh::AllocateIsMultiMaterialPoints() {
	isMultiMaterialPoints = new bool[GetNumberOfPoints()]; 
	for (int i = 0; i < GetNumberOfPoints(); i++) 
		isMultiMaterialPoints[i] = false;
}

void CSimplexMesh::SetSplitMesh(bool b) {
	isSplitMesh = b;

	// If the CSimplexMesh is a split mesh, then initialize 
	// the originalPointIndices array.
	if (isSplitMesh == true) {
		originalPointIndices = new int[GetNumberOfPoints()];

		for (int i = 0; i < GetNumberOfPoints(); i++) {
			originalPointIndices[i] = -1;
		}
	}
}

bool CSimplexMesh::GetSplitMesh() {
	return isSplitMesh;
}

void CSimplexMesh::SetOriginalPointIndex(int pointIndex, int originalIndex) {
	originalPointIndices[pointIndex] = originalIndex;
}

int CSimplexMesh::GetOriginalPointIndex(int pointIndex) {
	if (originalPointIndices[pointIndex] >= 0) {
		return originalPointIndices[pointIndex];
	}
	else {
		std::cerr << "Original point index is negetive. Not possible. " << std::endl;
		return -1;
	}
}
