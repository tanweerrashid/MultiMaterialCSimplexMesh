//#include "StdAfx.h"
#include "StdFunctions.h"
#include "float.h"
#include <iostream>
#include <vtkMatrix4x4.h>
#include <vtkMatrixToHomogeneousTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkImageSeedConnectivity.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkLandmarkTransform.h> 
#include <vtkTriangle.h> 
#include <vtkPolyDataNormals.h> 
#include <vtkImageGaussianSmooth.h> 
#include <vtkIdList.h>
#include <vtkDataArray.h>
#include <vtkImageConstantPad.h> 
#include <vtkImageAnisotropicDiffusion3D.h> 
#include <vtkImageAnisotropicDiffusion2D.h> 
#include <vtkImageGradientMagnitude.h> 
#include <vtkImageMirrorPad.h> 
#include <vtkDoubleArray.h>
#include <vtkImageEuclideanDistance.h> 
#include <vtkImageReslice.h>

#include <vtkCommand.h>

//#include "vtkPowerCrustSurfaceReconstruction.cxx"

using namespace std;

inline P_float dotproduct(P_float u[3], P_float v[3]) {return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);}
inline void crossproduct( P_float cp[3], P_float u[3], P_float v[3]) {cp[0]=u[1]*v[2]-u[2]*v[1]; cp[1]=u[2]*v[0]-u[0]*v[2]; cp[2]=u[0]*v[1]-u[1]*v[0];}
inline void crossproduct( P_float cp[3], const P_float u[3], const P_float v[3]) {cp[0]=u[1]*v[2]-u[2]*v[1]; cp[1]=u[2]*v[0]-u[0]*v[2]; cp[2]=u[0]*v[1]-u[1]*v[0];}
inline P_float tripleproduct( P_float u[3], P_float v[3], P_float w[3]) {P_float x[3]; crossproduct(x,u,v); return dotproduct(x,w);}
inline P_float norm(P_float u[3]) {return sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);}
inline P_float norm(const P_float u[3]) {return sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);}
inline P_float dist3D(P_float p1[3],P_float p2[3]) {P_float u[3]={p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]}; return norm(u);}

inline P_float GetTetrahedronVolume(P_float P0[3],P_float P1[3],P_float P2[3],P_float P3[3])
{
P_float u[3]={P0[0]-P1[0],P0[1]-P1[1],P0[2]-P1[2]};
P_float v[3]={P0[0]-P2[0],P0[1]-P2[1],P0[2]-P2[2]};
P_float w[3]={P0[0]-P3[0],P0[1]-P3[1],P0[2]-P3[2]};
P_float vol=(u[0]*(v[1]*w[2]-v[2]*w[1])+v[0]*(w[1]*u[2]-w[2]*u[1])+w[0]*(u[1]*v[2]-u[2]*v[1]))/6.;
return vol;
}

P_float GetCircumscribedRadius(P_float P0[3],P_float P1[3],P_float P2[3],P_float P3[3])
{
P_float vol=abs(GetTetrahedronVolume(P0,P1,P2,P3)); if(vol==0) return 1E10;
P_float P1P0[3]={P0[0]-P1[0],P0[1]-P1[1],P0[2]-P1[2]}; P_float p1p0=norm(P1P0);
P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]}; P_float p1p2=norm(P1P2);
P_float P1P3[3]={P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]}; P_float p1p3=norm(P1P3);
P_float P3P0[3]={P0[0]-P3[0],P0[1]-P3[1],P0[2]-P3[2]}; P_float p3p0=norm(P3P0);
P_float P2P3[3]={P3[0]-P2[0],P3[1]-P2[1],P3[2]-P2[2]}; P_float p2p3=norm(P2P3);
P_float P2P0[3]={P0[0]-P2[0],P0[1]-P2[1],P0[2]-P2[2]}; P_float p2p0=norm(P2P0);
P_float a=p3p0*p1p2,b=p1p0*p2p3,c=p1p3*p2p0;
P_float r=sqrt((a+b+c)*(a+b-c)*(b+c-a)*(a-b+c))/(24.*vol);

ofstream ff;
ff.open("C:\\Users\\ssult003\\Desktop\\CircumscribedRadius.txt");
ff << "P0" << P0[0] << ", " << P0[1] << ", " << P0[2] << endl;
ff << "P1" << P1[0] << ", " << P1[1] << ", " << P1[2] << endl;
ff << "P2" << P2[0] << ", " << P2[1] << ", " << P2[2] << endl;
ff << "P3" << P3[0] << ", " << P3[1] << ", " << P3[2] << endl;
ff << "radius: " << r << endl;

return r;
}


P_float GetInscribedRadius(P_float P0[3],P_float P1[3],P_float P2[3],P_float P3[3])
{
P_float vol=abs(GetTetrahedronVolume(P0,P1,P2,P3)); if(vol==0) return 0;
P_float P1P0[3]={P0[0]-P1[0],P0[1]-P1[1],P0[2]-P1[2]};
P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]};
P_float P1P3[3]={P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]};
P_float P3P0[3]={P0[0]-P3[0],P0[1]-P3[1],P0[2]-P3[2]};
P_float P3P2[3]={P2[0]-P3[0],P2[1]-P3[1],P2[2]-P3[2]};
P_float u[3],s=0; 
crossproduct(u,P1P0,P1P2); s+=norm(u)/2.;
crossproduct(u,P1P0,P1P3); s+=norm(u)/2.;
crossproduct(u,P1P2,P1P3); s+=norm(u)/2.;
crossproduct(u,P3P0,P3P2); s+=norm(u)/2.;
P_float r=3.*vol/s;
return r;
}

P_float GetTriangleQuality(P_float P1[3],P_float P2[3],P_float P3[3])
{
// qual=2r/R=16s2/(abc(a+b+c))
P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]}; P_float a=norm(P1P2);
P_float P1P3[3]={P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]}; P_float b=norm(P1P3);
P_float P3P2[3]={P2[0]-P3[0],P2[1]-P3[1],P2[2]-P3[2]}; P_float c=norm(P3P2);
P_float u[3]; crossproduct(u,P1P2,P1P3); 
P_float seize_s2=4.*dotproduct(u,u);
P_float qual=seize_s2/((a+b+c)*a*b*c);
return qual;
}

void Decompose(P_float vtg[3],P_float vn[3],P_float v[3],P_float n[3]) 
{
P_float cp=dotproduct(v,n);
vn[0]=cp*n[0];vn[1]=cp*n[1];vn[2]=cp*n[2];
vtg[0]=v[0]-vn[0]; vtg[1]=v[1]-vn[1]; vtg[2]=v[2]-vn[2];
}

void GetDerivatives(P_float extrapolation,P_float f[3],P_float df[6])  //0: isotropic ; 1: anisotropic
{
if(extrapolation==0) {
	df[0]=-1; df[1]=0; df[2]=-1; df[3]=0; df[4]=0; df[5]=-1; 
	return;
}

df[0]=f[0]*f[0];
df[1]=f[0]*f[1]; df[2]=f[1]*f[1];
df[3]=f[0]*f[2]; df[4]=f[1]*f[2]; df[5]=f[2]*f[2];
P_float l=df[0]+df[2]+df[5];
if(l==0) {
	for(int i=0;i<6;i++) df[i]=0;
}
else {
	for(int i=0;i<6;i++) df[i]=-extrapolation*df[i]/l;
}

df[0]-=1-extrapolation;df[2]-=1-extrapolation;df[5]-=1-extrapolation;
}

void GetForces_Shape(const P_float *P,const P_float *P1,const P_float *P2,const P_float *P3,const P_float *n,const P_float *pref,const P_float h,P_float F[3],P_float F1[3],P_float F2[3],P_float F3[3],P_float lambda[4])
{
F[0]=0;	 F[1]=0;  F[2]=0; F1[0]=0; F1[1]=0; F1[2]=0; F2[0]=0; F2[1]=0; F2[2]=0; F3[0]=0; F3[1]=0; F3[2]=0;

P_float N[3];
if(n==NULL) // normal is recomputed (for the energy method, we need to update the normal)
	{
	P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]},P1P3[3]={P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]},P2P3[3]={P3[0]-P2[0],P3[1]-P2[1],P3[2]-P2[2]};
	P_float N[3]; crossproduct(N,P2P3,P1P3);  P_float s=norm(N); N[0]=N[0]/s; N[1]=N[1]/s; N[2]=N[2]/s;
	}
else {N[0]=n[0]; N[1]=n[1]; N[2]=n[2];}

// compute forces for 4 particles based on linear combination
P_float PPi[9]={P1[0]-P[0],P2[0]-P[0],P3[0]-P[0],P1[1]-P[1],P2[1]-P[1],P3[1]-P[1],P1[2]-P[2],P2[2]-P[2],P3[2]-P[2]};
// invert PPi
P_float PPi_inv[9]; Invert_M33(PPi,PPi_inv);
F[0]=pref[0]*P1[0]+pref[1]*P2[0]+pref[2]*P3[0]-h*N[0]-P[0]; F[1]=pref[0]*P1[1]+pref[1]*P2[1]+pref[2]*P3[1]-h*N[1]-P[1]; F[2]=pref[0]*P1[2]+pref[1]*P2[2]+pref[2]*P3[2]-h*N[2]-P[2];
// compute weights
lambda[1]=PPi_inv[0]*F[0]+PPi_inv[1]*F[1]+PPi_inv[2]*F[2];
lambda[2]=PPi_inv[3]*F[0]+PPi_inv[4]*F[1]+PPi_inv[5]*F[2];
lambda[3]=PPi_inv[6]*F[0]+PPi_inv[7]*F[1]+PPi_inv[8]*F[2];
lambda[0]=lambda[1]+lambda[2]+lambda[3];
if(abs(lambda[0])<0.2) {lambda[0]=1; lambda[1]=0; lambda[2]=0; lambda[3]=0; return;}
// compute forces
lambda[1]=lambda[1]/lambda[0]; lambda[2]=lambda[2]/lambda[0]; lambda[3]=lambda[3]/lambda[0]; lambda[0]=1; 
F1[0]=-F[0]*lambda[1]; F1[1]=-F[1]*lambda[1]; F1[2]=-F[2]*lambda[1];
F2[0]=-F[0]*lambda[2]; F2[1]=-F[1]*lambda[2]; F2[2]=-F[2]*lambda[2];
F3[0]=-F[0]*lambda[3]; F3[1]=-F[1]*lambda[3]; F3[2]=-F[2]*lambda[3];


/*
// compute the gradients of the energy |P~ - P|= |pref[0]P1+pref[1]P2+(1-pref[0]-pref[2])P3+hn - P|
lambda[1]=pref[0]; lambda[2]=pref[1]; lambda[3]=pref[2]; lambda[0]=1; 
P_float DC[3],DC1[3],DC2[3],DC3[3];
P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]},P1P3[3]={P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]},P2P3[3]={P3[0]-P2[0],P3[1]-P2[1],P3[2]-P2[2]},C,df[9];
P_float N[3];crossproduct(N,P2P3,P1P3);  P_float s=norm(N); N[0]=N[0]/s; N[1]=N[1]/s; N[2]=N[2]/s;
P_float Pt[3]={pref[0]*P1[0]+pref[1]*P2[0]+pref[2]*P3[0]-h*N[0],pref[0]*P1[1]+pref[1]*P2[1]+pref[2]*P3[1]-h*N[1],pref[0]*P1[2]+pref[1]*P2[2]+pref[2]*P3[2]-h*N[2]};
// DC=PtP/|PtP|
DC[0]=P[0]-Pt[0]; DC[1]=P[1]-Pt[1]; DC[2]=P[2]-Pt[2]; C=norm(DC); 
if(C<1E-6) return;
DC[0]=DC[0]/C; DC[1]=DC[1]/C; DC[2]=DC[2]/C;
P_float dc=dotproduct(DC,DC);
// DC1=(pref[0]I + h dn/dP1)PPt/|PtP| = (pref[0]I + h d(P2P1 n P2P3/s)/d(P2P1))PPt/|PtP| 
GetDerivativesOfCrossProduct(P2P3,N,s,df); for(int i=0;i<9;i++) df[i]=h*df[i]; df[0]+=pref[0]; df[4]+=pref[0]; df[8]+=pref[0];
DC1[0]=-df[0]*DC[0]-df[1]*DC[1]-df[2]*DC[2]; DC1[1]=-df[3]*DC[0]-df[4]*DC[1]-df[5]*DC[2]; DC1[2]=-df[6]*DC[0]-df[7]*DC[1]-df[8]*DC[2];
P_float dc1=dotproduct(DC1,DC1);
// DC2=(pref[1]I + h dn/dP1)PPt/|PtP| = (pref[1]I + h d(P1P3 n P1P2/s)/d(P1P2))PPt/|PtP| 
GetDerivativesOfCrossProduct(P1P3,N,s,df); for(int i=0;i<9;i++) df[i]=-h*df[i]; df[0]+=pref[1]; df[4]+=pref[1]; df[8]+=pref[1];
DC2[0]=-df[0]*DC[0]-df[1]*DC[1]-df[2]*DC[2]; DC2[1]=-df[3]*DC[0]-df[4]*DC[1]-df[5]*DC[2]; DC2[2]=-df[6]*DC[0]-df[7]*DC[1]-df[8]*DC[2];
P_float dc2=dotproduct(DC2,DC2);
// DC3=(pref[2]I + h dn/dP1)PPt/|PtP| = (pref[2]I + h d(P1P3 n P1P2/s)/d(P1P3))PPt/|PtP| 
GetDerivativesOfCrossProduct(P1P2,N,s,df); for(int i=0;i<9;i++) df[i]=h*df[i]; df[0]+=pref[2]; df[4]+=pref[2]; df[8]+=pref[2];
DC3[0]=-df[0]*DC[0]-df[1]*DC[1]-df[2]*DC[2]; DC3[1]=-df[3]*DC[0]-df[4]*DC[1]-df[5]*DC[2]; DC3[2]=-df[6]*DC[0]-df[7]*DC[1]-df[8]*DC[2];
P_float dc3=dotproduct(DC3,DC3);
// compute the displacment Fi= -s * DCi = -C/sum(DC^2) * DCi 
P_float scaling=C/(dc+dc1+dc2+dc3);
F[0]=-scaling*DC[0]; F[1]=-scaling*DC[1]; F[2]=-scaling*DC[2];
F1[0]=-scaling*DC1[0]; F1[1]=-scaling*DC1[1]; F1[2]=-scaling*DC1[2];
F2[0]=-scaling*DC2[0]; F2[1]=-scaling*DC2[1]; F2[2]=-scaling*DC2[2];
F3[0]=-scaling*DC3[0]; F3[1]=-scaling*DC3[1]; F3[2]=-scaling*DC3[2];
*/

//test linear momentum
//P_float temp[3]={F[0]+F1[0]+F2[0]+F3[0],F[1]+F1[1]+F2[1]+F3[1],F[2]+F1[2]+F2[2]+F3[2]};
//test angular momentum
//P_float temp2[3],temp3[3]={0,0,0};
//crossproduct(temp2,P,F); temp3[0]+=temp2[0];temp3[1]+=temp2[1];temp3[2]+=temp2[2];
//crossproduct(temp2,P1,F1); temp3[0]+=temp2[0];temp3[1]+=temp2[1];temp3[2]+=temp2[2];
//crossproduct(temp2,P2,F2); temp3[0]+=temp2[0];temp3[1]+=temp2[1];temp3[2]+=temp2[2];
//crossproduct(temp2,P3,F3); temp3[0]+=temp2[0];temp3[1]+=temp2[1];temp3[2]+=temp2[2];
}


void GetDerivativesOfCrossProduct(P_float v2[3],P_float n[3],P_float normcross,P_float df[9])  
{
// return the dn/dv1=d(v1 n v2/|v1 n v2|)/d(v1)
P_float cp[3]; crossproduct(cp,n,v2);
df[0]=cp[0]*n[0];		df[1]=cp[0]*n[1]-v2[2];	df[2]=cp[0]*n[2]+v2[1];
df[3]=cp[1]*n[0]+v2[2];	df[4]=cp[1]*n[1];		df[5]=cp[1]*n[2]-v2[0];
df[6]=cp[2]*n[0]-v2[1];	df[7]=cp[2]*n[1]+v2[0];	df[8]=cp[2]*n[2];
for(int i=0;i<9;i++) df[i]=df[i]/normcross;
}


void GetPolyDataNormal(vtkPolyData* model,const int indexp,P_float n[3])
{
n[0]=model->GetPointData()->GetNormals()->GetComponent(indexp,0);
n[1]=model->GetPointData()->GetNormals()->GetComponent(indexp,1);
n[2]=model->GetPointData()->GetNormals()->GetComponent(indexp,2);
}

P_float* GetPolyDataNormals(vtkPolyData* model)
{
P_float *n=new P_float[3*model->GetNumberOfPoints()];
for(int i=0;i<model->GetNumberOfPoints();i++) 
	{
	n[3*i]=model->GetPointData()->GetNormals()->GetComponent(i,0);
	n[3*i+1]=model->GetPointData()->GetNormals()->GetComponent(i,1);
	n[3*i+2]=model->GetPointData()->GetNormals()->GetComponent(i,2);
	}
return n;
}

P_float* GetPolyDataPoints(vtkPolyData* model)
{
P_float *pts=new P_float[3*model->GetNumberOfPoints()];
for(int i=0;i<model->GetNumberOfPoints();i++) model->GetPoint(i,pts+3*i);
return pts;
}

void Identity(P_float M[16]) {M[0]=1; M[1]=0; M[2]=0; M[3]=0; M[4]=0; M[5]=1; M[6]=0; M[7]=0; M[8]=0; M[9]=0; M[10]=1; M[11]=0; M[12]=0; M[13]=0; M[14]=0; M[15]=1;}
void Multi_M(P_float M[16],P_float M2[16],P_float M3[16]) {M[0]= M2[0]*M3[0]+M2[1]*M3[4]+M2[2]*M3[8]+M2[3]*M3[12]; M[1]= M2[0]*M3[1]+M2[1]*M3[5]+M2[2]*M3[9]+M2[3]*M3[13]; M[2]= M2[0]*M3[2]+M2[1]*M3[6]+M2[2]*M3[10]+M2[3]*M3[14]; M[3]= M2[0]*M3[3]+M2[1]*M3[7]+M2[2]*M3[11]+M2[3]*M3[15]; M[4]= M2[4]*M3[0]+M2[5]*M3[4]+M2[6]*M3[8]+M2[7]*M3[12]; M[5]= M2[4]*M3[1]+M2[5]*M3[5]+M2[6]*M3[9]+M2[7]*M3[13]; M[6]= M2[4]*M3[2]+M2[5]*M3[6]+M2[6]*M3[10]+M2[7]*M3[14]; M[7]= M2[4]*M3[3]+M2[5]*M3[7]+M2[6]*M3[11]+M2[7]*M3[15]; M[8]= M2[8]*M3[0]+M2[9]*M3[4]+M2[10]*M3[8]+M2[11]*M3[12]; M[9]= M2[8]*M3[1]+M2[9]*M3[5]+M2[10]*M3[9]+M2[11]*M3[13]; M[10]=M2[8]*M3[2]+M2[9]*M3[6]+M2[10]*M3[10]+M2[11]*M3[14]; M[11]=M2[8]*M3[3]+M2[9]*M3[7]+M2[10]*M3[11]+M2[11]*M3[15]; M[12]=M2[12]*M3[0]+M2[13]*M3[4]+M2[14]*M3[8]+M2[15]*M3[12]; M[13]=M2[12]*M3[1]+M2[13]*M3[5]+M2[14]*M3[9]+M2[15]*M3[13]; M[14]=M2[12]*M3[2]+M2[13]*M3[6]+M2[14]*M3[10]+M2[15]*M3[14]; M[15]=M2[12]*M3[3]+M2[13]*M3[7]+M2[14]*M3[11]+M2[15]*M3[15];}
void Multi_MsT(P_float U[3],P_float M[6],P_float T[3]) {U[0]= M[0]*T[0]+M[1]*T[1]+M[3]*T[2]; U[1]= M[1]*T[0]+M[2]*T[1]+M[4]*T[2]; U[2]= M[3]*T[0]+M[4]*T[1]+M[5]*T[2];}
void Invert_M(P_float M[16],P_float M_inv[16])
{
int i,j;
double** M_d=new double*[4];  double** M_inv_d=new double*[4];
for(i=0;i<4;i++) { M_d[i]=new double[4]; M_inv_d[i]=new double[4]; for(j=0;j<4;j++) M_d[i][j]=M[j+i*4]; }
//vtkMath* math=vtkMath::New();
vtkMath::InvertMatrix(M_d,M_inv_d,4);
//	math->Delete();
for(i=0;i<4;i++) { for(j=0;j<4;j++) M_inv[j+i*4]=M_inv_d[i][j]; free(M_d[i]); free(M_inv_d[i]); }
free(M_d); free(M_inv_d);
}

void Invert_M33(P_float M[9],P_float M_inv[9])
{
P_float det = M[0] * ( M[4]*M[8] - M[7]*M[5] )- M[1] * ( M[3]*M[8] - M[6]*M[5] ) + M[2] * ( M[3]*M[7] - M[6]*M[4] );
if(det==0) {M_inv[0]=1;	M_inv[1]=0;	M_inv[2]=0; M_inv[3]=0;	M_inv[4]=1;	M_inv[5]=0; M_inv[6]=0;	M_inv[7]=0;	M_inv[8]=1; return;}
M_inv[0] =  ( M[4]*M[8] - M[5]*M[7] ) / det;	M_inv[1] = -( M[1]*M[8] - M[7]*M[2] ) / det;	M_inv[2] =  ( M[1]*M[5] - M[4]*M[2] ) / det;
M_inv[3] = -( M[3]*M[8] - M[5]*M[6] ) / det;	M_inv[4] =  ( M[0]*M[8] - M[6]*M[2] ) / det;	M_inv[5] = -( M[0]*M[5] - M[3]*M[2] ) / det;
M_inv[6] =  ( M[3]*M[7] - M[6]*M[4] ) / det;	M_inv[7] = -( M[0]*M[7] - M[6]*M[1] ) / det;	M_inv[8] =  ( M[0]*M[4] - M[1]*M[3] ) / det;
}

void Invert_MSymmetric(P_float M[6],P_float M_inv[6])
{
P_float det=M[0]*M[2]*M[5]+2*M[1]*M[3]*M[4]-M[0]*M[4]*M[4]-M[2]*M[3]*M[3]-M[5]*M[1]*M[1];
if(det!=0)
	{
	det=1/det;  
	M_inv[0]=det*(M[2]*M[5]-M[4]*M[4]);
	M_inv[1]=det*(M[4]*M[3]-M[1]*M[5]); M_inv[2]=det*(M[0]*M[5]-M[3]*M[3]);
	M_inv[3]=det*(M[1]*M[4]-M[2]*M[3]); M_inv[4]=det*(M[1]*M[3]-M[0]*M[4]); M_inv[5]=det*(M[0]*M[2]-M[1]*M[1]);
	}
else
	{
	M_inv[0]=1;
	M_inv[1]=0; M_inv[2]=1;
	M_inv[3]=0; M_inv[4]=0; M_inv[5]=1;
	}
}

void Transform(P_float pin[3], P_float pout[3], P_float M[16]) { 
	pout[0] = M[0] * pin[0] + M[1] * pin[1] + M[2] * pin[2] + M[3]; 
	pout[1] = M[4] * pin[0] + M[5] * pin[1] + M[6] * pin[2] + M[7]; 
	pout[2] = M[8] * pin[0] + M[9] * pin[1] + M[10] * pin[2] + M[11]; 
}

void Transform_R(P_float pin[3],P_float pout[3], P_float M[16]) { pout[0]=M[0]*pin[0]+M[1]*pin[1]+M[2]*pin[2]; pout[1]=M[4]*pin[0]+M[5]*pin[1]+M[6]*pin[2]; pout[2]=M[8]*pin[0]+M[9]*pin[1]+M[10]*pin[2]; }
void Transform(vtkPolyData* model,vtkPolyData* model_out, P_float rt[3],P_float tr[3])
{
P_float M[16]; NtoM(tr,rt,M);
Transform(model,model_out,M);
}

void Transform(vtkPolyData* model,vtkPolyData* model_out, P_float M[16])
{
vtkMatrix4x4* M44=vtkMatrix4x4::New();
for(int i=0;i<4;i++) for(int j=0;j<4;j++) M44->SetElement(i,j,(double)M[4*i+j]);

vtkMatrixToHomogeneousTransform* transform=vtkMatrixToHomogeneousTransform::New();
	transform->SetInput(M44);
	transform->Update();
vtkTransformPolyDataFilter* transf=vtkTransformPolyDataFilter::New();
	transf->SetTransform(transform);
	transf->SetInput(model);
	transf->Update();
model_out->DeepCopy(transf->GetOutput());
transform->Delete();
transf->Delete();
M44->Delete();
}
/*
void Transform(CSimplexMesh* model, P_float M[16])
{
P_float p[3],po[3];
for(int i=0;i<model->GetNumberOfPoints();i++)
	{
	model->GetPoint(i,p);
	Transform(p,po,M);
	model->SetPoint(i,po);
	}
model->Equilibrium();
}
*/


/*void PowerCrust(vtkPolyData* model)
{
vtkPowerCrustSurfaceReconstruction* filter=vtkPowerCrustSurfaceReconstruction::New();
filter->SetInput(model);
filter->Update();
vtkPolyData* ma=vtkPolyData::New(); ma->DeepCopy(filter->GetMedialSurface());
filter->Delete();

// mapping
int i;
vtkPointLocator* PointsLocator=vtkPointLocator::New(); PointsLocator->SetDataSet(model); PointsLocator->BuildLocator();
int nbp=ma->GetNumberOfPoints();
double *uca=new double[nbp]; for(i=0;i<nbp;i++) uca[i]=0;
for(i=0;i<nbp;i++) PointsLocator->FindClosestPointWithinRadius(100000,ma->GetPoint(i),uca[i]);
for(i=0;i<nbp;i++) uca[i]=sqrt(uca[i])/15.406418144828262;

model->DeepCopy(ma); ma->Delete();
vtkDoubleArray* da=vtkDoubleArray::New();	da->SetArray(uca,nbp,0);	model->GetPointData()->SetScalars(da);
}*/


void PerturbSurf(int nb_points,P_float* pts,P_float* normals,P_float* mass,P_float maxlength)
{
P_float d,*p,*n;
if(mass==NULL)
for(int i=0; i<nb_points;i++)
	{
	p=pts+3*i; n=normals+3*i;
	d=rand()*2*maxlength/RAND_MAX-maxlength; 
	p[0]+=d*n[0]; p[1]+=d*n[1]; p[2]+=d*n[2]; 
	}
else 
for(int i=0; i<nb_points;i++)
if(mass[3*i]!=0)
	{
	p=pts+3*i; n=normals+3*i;
	d=rand()*2*maxlength/RAND_MAX-maxlength; 
	p[0]+=d*n[0]; p[1]+=d*n[1]; p[2]+=d*n[2]; 
	}
}


void QtoM(P_float tr[3],P_float q[4],P_float M[16],bool post_rotation)
{
P_float xs = q[1]*2., ys = q[2]*2., zs = q[3]*2.;
P_float wx = q[0]*xs, wy = q[0]*ys, wz = q[0]*zs;
P_float xx = q[1]*xs, xy = q[1]*ys, xz = q[1]*zs;
P_float yy = q[2]*ys, yz = q[2]*zs, zz = q[3]*zs;

M[0] = 1.0 - (yy + zz); M[1]= xy - wz; M[2] = xz + wy;
M[4] = xy + wz; M[5] = 1.0 - (xx + zz); M[6] = yz - wx;
M[8] = xz - wy; M[9] = yz + wx; M[10] = 1.0 - (xx + yy);
M[12] = M[13] = M[14] = 0.0; M[15] = 1.0;

if(!post_rotation) { M[3]=tr[0]; M[7]=tr[1]; M[11]=tr[2]; }
else 
	{
	M[3]=M[0]*tr[0]+M[1]*tr[1]+M[2]*tr[2];
	M[7]=M[4]*tr[0]+M[5]*tr[1]+M[6]*tr[2];
	M[11]=M[8]*tr[0]+M[9]*tr[1]+M[10]*tr[2];
	} 
}
void MtoQ(P_float tr[3],P_float q[4],P_float M[16],bool post_rotation)
{
q[0]=0.5*sqrt(M[0]+M[5]+M[10]+M[15]);
q[1]=0.25*(M[9]-M[6])/q[0];
q[2]=0.25*(M[2]-M[8])/q[0];
q[3]=0.25*(M[4]-M[1])/q[0];

if(!post_rotation) {tr[0]=M[3];tr[1]=M[7];tr[2]=M[11];}
else 
	{
	P_float R_inv[16]; P_float tr2[3]={M[3],M[7],M[11]}; M[3]=0; M[7]=0; M[11]=0;
	Invert_M(M,R_inv);
	tr[0]=R_inv[0]*tr2[0]+R_inv[1]*tr2[1]+R_inv[2]*tr2[2];
	tr[1]=R_inv[4]*tr2[0]+R_inv[5]*tr2[1]+R_inv[6]*tr2[2];
	tr[2]=R_inv[8]*tr2[0]+R_inv[9]*tr2[1]+R_inv[10]*tr2[2];
	}
}
/*
void QtoAN(P_float v1[3],P_float v2[3],P_float q[4],P_float *theta)
{
// V1 = normal of Qvect with x=0 (or y=0 if infinite number of sol) 
v1[0]=0; v1[1]=0; v1[2]=0;
P_float temp=sqrt(q[2]*q[2]+q[3]*q[3]);
if(temp==0) { temp=sqrt(q[1]*q[1]+q[3]*q[3]); v1[0]=-q[3]/temp; v1[2]=q[1]/temp; }
else { v1[1]=-q[3]/temp; v1[2]=q[2]/temp; }

// V2 = cross-product(Q,v1)
v2[0]=q[2]*v1[2] - q[3]*v1[1];
v2[1]=q[3]*v1[0] - q[1]*v1[2];
v2[2]=q[1]*v1[1] - q[2]*v1[0];

P_float n=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
v2[0]=v2[0]/n;
v2[1]=v2[1]/n;
v2[2]=v2[2]/n;

*theta=2*(P_float)acos(q[0]); // Qreal=cos(theta/2)
}
*/
/*
void ANtoQ(P_float v1[3],P_float p1,P_float v2[3],P_float p2,P_float theta,P_float q[4])
{
q[0]=cos(theta/2); // Qreal=cos(theta/2)

q[1]=q[1]+p1*v1[0]+p2*v2[0];
q[2]=q[2]+p1*v1[1]+p2*v2[1];
q[3]=q[3]+p1*v1[2]+p2*v2[2];
QNorm(q);
}
*/

void QtoAN(P_float q[4],P_float *alpha,P_float *beta,P_float *theta)
{
P_float r=sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
*alpha=acos(q[3]/r);
*beta=atan2(q[2],q[1]); 
//if(q[1]<0 && q[2]<0) *beta+=PI;
*theta=2*(P_float)acos(q[0]); // Qreal=cos(theta/2)
}


void ANtoQ(P_float alpha,P_float beta,P_float theta,P_float q[4])
{
q[0]=cos(theta/2); // Qreal=cos(theta/2)
q[1]=sin(alpha)*cos(beta);
q[2]=sin(alpha)*sin(beta);
q[3]=cos(alpha);
QNorm(q);
}

void MNorm(P_float M[16])
{
P_float n;
n=sqrt(M[0]*M[0]+M[4]*M[4]+M[8]*M[8]);
M[0]=M[0]/n;M[4]=M[4]/n;M[8]=M[8]/n;
n=sqrt(M[1]*M[1]+M[5]*M[5]+M[9]*M[9]);
M[1]=M[1]/n;M[5]=M[5]/n;M[9]=M[9]/n;
M[2]=M[4]*M[9]-M[8]*M[5];
M[6]=M[8]*M[1]-M[9]*M[0];
M[10]=M[0]*M[5]-M[1]*M[4];
n=sqrt(M[2]*M[2]+M[6]*M[6]+M[10]*M[10]);
M[2]=M[2]/n;M[6]=M[6]/n;M[10]=M[10]/n;
M[12]=0; M[13]=0; M[14]=0; M[15]=1;
}


void QNorm(P_float q[4])
{
// n=1/a=sqrt(1-r2)/norm

P_float n=sqrt(1-q[0]*q[0])/sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
q[1]=q[1]*n;
q[2]=q[2]*n;
q[3]=q[3]*n;
}

void QComplete(P_float qout[4],P_float q[3])
{
qout[0]=q[0];
qout[1]=q[1];
qout[2]=q[2];
qout[3]=sqrt(1-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
}


void QSlerp(P_float q[4],P_float t,P_float qs[4],P_float qe[4])
{
int i,flag=0;
P_float dotp = qs[0]*qe[0] + qs[1]*qe[1] + qs[2]*qe[2] + qs[3]*qe[3];

if ( dotp < 0 ) { dotp *= -1; for (i=0;i<4;i++) qs[i] *= -1; }
if(dotp>=1) { memcpy(q,qs,4*sizeof(P_float)); return; }

P_float omega=acos(dotp);
P_float sinOmega=sin(omega);

if(t>1) {t=1/t; flag=1;}
if(t<0) {t=t/(t-1); flag=2;}
P_float a=sin((1.0 - t)*omega) / sinOmega;
P_float b=sin(t*omega) / sinOmega;

if(flag==0) for (i=0;i<4;i++) q[i] = a*qs[i]+b*qe[i];
else if(flag==1) for (i=0;i<4;i++) q[i] = (qe[i]-a*qs[i])/b;
else if(flag==2) for (i=0;i<4;i++) q[i] = (qs[i]-b*qe[i])/a;
}


void NtoM(P_float tr[3],P_float rt[3],P_float M[16])
{
P_float ca=cos(PI*rt[0]/180.),cb=cos(PI*rt[1]/180.),cc=cos(PI*rt[2]/180.);
P_float sa=sin(PI*rt[0]/180.),sb=sin(PI*rt[1]/180.),sc=sin(PI*rt[2]/180.);

M[0]=(ca*cc)-(sa*sb*sc); M[1]=-sa*cb; M[2]=(sc*ca)+(sa*sb*cc); 
M[4]=(sa*cc)+(ca*sb*sc); M[5]=ca*cb; M[6]=(sa*sc)-(ca*cc*sb); 
M[8]=-sc*cb; M[9]=sb; M[10]=cb*cc; 
M[12]=0; M[13]=0; M[14]=0; M[15]=1;


P_float m5=M[5],m9=M[9];
if(m5==0) m5=1E-10; if(m9==1) m9=1+1E-10; // singularities
M[11]=tr[0];
M[7]=(tr[2]-tr[0]*m9-M[1]*tr[1]/m5)*m5/(1-m9*m9); 
M[3]=(tr[1]+M[7]*M[1])/m5; 
}


void MtoN(P_float tr[3],P_float rt[3],P_float M[16])
{
P_float F[3]={M[5],-M[1],0};
P_float n=sqrt(F[0]*F[0]+F[1]*F[1]); F[0]=F[0]/n; F[1]=F[1]/n; 

rt[0]=(P_float)180*acos(F[0])/PI; if(F[1]<0) rt[0]=-rt[0]; //flexion=angle(X,F)=rotation around Z
rt[1]=(P_float)180*acos(sqrt(M[1]*M[1]+M[5]*M[5]))/PI; if(M[9]<0) rt[1]=-rt[1]; //abduction=angle(y,P(XY))=rotation around F
n=sqrt(M[8]*M[8]+M[10]*M[10]);
rt[2]=(P_float)180*acos(M[10]/n)/PI; 	if((M[8]/n)>0) rt[2]=-rt[2]; //rotation=angle(x,F)=rotation around y

P_float u[3]={M[3],M[7],M[11]}; tr[0]=u[2]; // tr0=u.e1=u.[0,0,1]
P_float e3[3]={M[1],M[5],M[9]}; tr[2]=dotproduct(u,e3); // tr2=u.e3
P_float e2[3]={M[5],-M[1],0}; tr[1]=dotproduct(u,e2); // tr1=u.(e3 n e1)
}


void GetConic(P_float M[16],P_float alpha,P_float beta) // alpha: elevation ; beta: revolution
{
P_float q[4],tr[3]={0,0,0};
q[0]=cos(PI*alpha/(2.*180.)); q[1]=cos(PI*beta/180.); q[2]=0; q[3]=sin(PI*beta/180.); QNorm(q);
QtoM(tr,q,M,false);
}




P_float* GetDistMap_pointdistances(vtkPolyData* model1,vtkPolyData* model2,vtkIdList** PointCells2,P_float* CellNormals2,vtkPointLocator* PointsLocator2,P_float M[16])
{
// returns model1(M.model2)
P_float* distances=new P_float[model1->GetNumberOfPoints()];
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1];

P_float M_inv[16];  Invert_M(M,M_inv);

for(int i=0;i<model1->GetNumberOfPoints();i++)
    {
    distances[i]=1E10;
	closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10;
    GetClosest(model2,M,M_inv,PointCells2,CellNormals2,PointsLocator2,closest,model1->GetPoint(i),NULL,NULL);
	if(closest->dist2!=1E10)  {distances[i]=sqrt(closest->dist2); if(closest->side==-1) distances[i]=-distances[i];}
	free(closest->pts); free(closest->weights);
	}
free(closest); 
return distances;
}





void GetDistMap_collisions(P_float *distances,P_float *distances_diff,vtkPolyData* model1,vtkIdList** PointCells1,P_float* CellNormals1,vtkPointLocator* PointsLocator1,vtkPolyData* model2,bool *selectP,P_float *distances_ref,P_float M[16])
{
// returns distance_ref(model2) - distmap1(M.model2)
P_float M_inv[16];  Invert_M(M,M_inv);
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1]; 

for(int i=0;i<model2->GetNumberOfPoints();i++)
    {
    distances_diff[i]=0; distances[i]=1E10;
    if(selectP[i])
        {
		closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10;
		GetClosest(model1,M_inv,M,PointCells1,CellNormals1,PointsLocator1,closest,model2->GetPoint(i),NULL,NULL);
        if(closest->dist2!=1E10)  
            {
			distances[i]=sqrt(closest->dist2); if(closest->side==-1) distances[i]=-distances[i];
            distances_diff[i]=distances_ref[i] - sqrt(closest->dist2);
            }
		free(closest->pts); free(closest->weights);
        }
    }
free(closest); 
}



P_float* GetDistMap_pointdistances(vtkPolyData* model1,vtkStructuredPoints* distmap2,P_float M[16])
{
// returns model1(M.distmap2)
P_float* distances=new P_float[model1->GetNumberOfPoints()];
double* ptr2=(double*)distmap2->GetScalarPointer(); 
P_float *p,pt[3];
P_float M_inv[16];	Invert_M(M,M_inv);
P_float bounds[3][2]; GetBoundingBox(distmap2,M,bounds);
int ID;
for(int i=0;i<model1->GetNumberOfPoints();i++)
	{
	distances[i]=1E10;
	p=model1->GetPoint(i);
	if(p[0]>bounds[0][0] && p[0]<bounds[0][1] && p[1]>bounds[1][0] && p[1]<bounds[1][1] && p[2]>bounds[2][0] && p[2]<bounds[2][1])
		{
		Transform(p,pt,M_inv);
		ID=distmap2->FindPoint(pt);
		if(ID!=-1)	distances[i]=*(ptr2+ID); 
		}
	}
return distances;
}

vtkStructuredPoints* GetDistMap_thickness(vtkStructuredPoints* distmap1,vtkStructuredPoints* distmap2,P_float M[16])
{
// returns distmap1+M.distmap2
vtkStructuredPoints* Thickness=vtkStructuredPoints::New();
	Thickness->DeepCopy(distmap2);
double* ptr_thickness=(double*)Thickness->GetScalarPointer(); 
double* ptr1=(double*)distmap1->GetScalarPointer(); 
P_float *p,pt[3];
P_float M_inv[16];	Invert_M(M,M_inv);
P_float bounds[3][2]; GetBoundingBox(distmap1,M_inv,bounds);
bool ok;
int ID;
for(int i=0;i<Thickness->GetNumberOfPoints();i++)
	{
	ok=false;
	p=Thickness->GetPoint(i);
	if(p[0]>bounds[0][0] && p[0]<bounds[0][1] && p[1]>bounds[1][0] && p[1]<bounds[1][1] && p[2]>bounds[2][0] && p[2]<bounds[2][1])
		{
		Transform(p,pt,M);
		ID=distmap1->FindPoint(pt);
		if(ID!=-1)	{*(ptr_thickness+i)+=*(ptr1+ID); ok=true;}
		}
	if(!ok) *(ptr_thickness+i)=(double)1E10;
	}
return Thickness;
}

void GetDistMap_collisions(P_float *distances,P_float *distances_diff,vtkStructuredPoints* distmap1,vtkPolyData* model2,bool *selectP,P_float *distances_ref,vtkStructuredPoints* Thickness,P_float M[16])
{
// returns [thickness(model2)+distance_ref(model2)]/2 - distmap1(M.model2)
P_float M_inv[16];	Invert_M(M,M_inv);
P_float bounds[3][2]; GetBoundingBox(distmap1,M_inv,bounds);
P_float *p,pt[3];
int ID;
double* ptr1=(double*)distmap1->GetScalarPointer(); 
double* ptr_thickness;
if(Thickness!=NULL) ptr_thickness=(double*)Thickness->GetScalarPointer(); 
for(int i=0;i<model2->GetNumberOfPoints();i++)
	{
	distances_diff[i]=0; distances[i]=1E10;
	if(selectP[i])
		{
		p=model2->GetPoint(i);
		if(p[0]>bounds[0][0] && p[0]<bounds[0][1] && p[1]>bounds[1][0] && p[1]<bounds[1][1] && p[2]>bounds[2][0] && p[2]<bounds[2][1])
			{
			Transform(p,pt,M);
			ID=distmap1->FindPoint(pt);
			if(ID!=-1)	
				if(Thickness!=NULL) 
					{if(*(ptr_thickness+ID)!=(double)1E10) {distances[i]=*(ptr1+ID); distances_diff[i]=(*(ptr_thickness+ID)+distances_ref[i])/2. - *(ptr1+ID);}}
				else {distances[i]=*(ptr1+ID); distances_diff[i]=distances_ref[i] - *(ptr1+ID);}
			}
		}
	}
}



int* GetDistMap_collisions(P_float force[3],vtkStructuredPoints* distmap1,vtkStructuredPoints* distmap2,P_float dmin,P_float M[16])
{
P_float f[3]={0,0,0};
int* ret=new int[1]; ret[0]=0; // collisions: [nb,col1ID_1,col1ID_2,col2ID_1,col2ID_2,...]
int* ptr1_vol=(int*)distmap1->GetScalarPointer(); 
int* ptr2_vol=(int*)distmap2->GetScalarPointer(); 
double p[3];
P_float d,pt[3],g1[3],g2[3],nrm,u[3];
double *spac1=distmap1->GetSpacing(),*org1=distmap1->GetOrigin(); 
double *spac2=distmap2->GetSpacing(),*org2=distmap2->GetOrigin();
int ID,count=0;

for(int i=0;i<distmap2->GetNumberOfPoints();i++)
	{
	d=*(ptr2_vol+i);
	if(d<dmin)
		{
		distmap2->GetPoint(i,p);
		Transform(p,pt,M);
		ID=distmap1->FindPoint(pt);
		if(ID!=-1)
			{
			d+=*(ptr1_vol+ID);
			if(d<dmin)
				{
				ret[0]++; ret=(int*)realloc(ret,(2*ret[0]+1)*sizeof(int)); 
				ret[2*ret[0]-1]=ID; ret[2*ret[0]]=i;

				if(Gradient_Int(g1,distmap1,pt[0],pt[1],pt[2],1,1)!=-1E10)
					if(Gradient_Int(g2,distmap2,p[0],p[1],p[2],1,1)!=-1E10)
						{
						nrm=norm(g1); if(nrm!=0) {g1[0]=g1[0]/nrm; g1[1]=g1[1]/nrm; g1[2]=g1[2]/nrm;}
						nrm=norm(g2); if(nrm!=0) {g2[0]=g2[0]/nrm; g2[1]=g2[1]/nrm; g2[2]=g2[2]/nrm;}
						u[0]=g2[0]-g1[0]; u[1]=g2[1]-g1[1]; u[2]=g2[2]-g1[2];  
						nrm=norm(u); if(nrm!=0) {u[0]=u[0]/nrm; u[1]=u[1]/nrm; u[2]=u[2]/nrm;}
						f[0]+=(dmin-d)*u[0]; f[1]+=(dmin-d)*u[1]; f[2]+=(dmin-d)*u[2];
						count++;
						}
				}
			}
		}
	}
f[0]=f[0]/count; f[1]=f[1]/count; f[2]=f[2]/count;
P_float Minv[16]; Invert_M(M,Minv); Transform_R(f,force,Minv);
return ret;
}


P_float FitSphere(P_float C[3],int nb_pts, P_float *pts,P_float thresh)
{
int i,count=0;
P_float d=1E10,nb_pts2=nb_pts*nb_pts,li,L[3],Lm,m[3]={0,0,0},C2[3]; 
for(i=0;i<nb_pts;i++) {m[0]+=pts[3*i]; m[1]+=pts[3*i+1]; m[2]+=pts[3*i+2];} m[0]=m[0]/nb_pts; m[1]=m[1]/nb_pts; m[2]=m[2]/nb_pts;
memcpy(C,m,3*sizeof(P_float));

while(d>thresh)
	{
	L[0]=0; L[1]=0; L[2]=0; Lm=0;
	for(i=0;i<nb_pts;i++)
		{
		li=dist3D(C,pts+3*i);
		if(li!=0) {L[0]+=(C[0]-pts[3*i])/li; L[1]+=(C[1]-pts[3*i+1])/li;	L[2]+=(C[2]-pts[3*i+2])/li; 	Lm+=li; }
		}
	C2[0]=m[0] + ((Lm*L[0])/nb_pts2); C2[1]=m[1] + ((Lm*L[1])/nb_pts2); C2[2]=m[2] + ((Lm*L[2])/nb_pts2);
	d=dist3D(C,C2);
	memcpy(C,C2,3*sizeof(P_float));
	count++; if(count==10000) d=0;
	}
return Lm/(P_float)nb_pts;
}


P_float FitDoubleSphere(P_float C[3],int nb_pts1, P_float *pts1,int nb_pts2, P_float *pts2,P_float thresh)
{
int i,count=0,nb_pts=nb_pts1+nb_pts2;;
P_float d=1E10,li,L1[3],L2[3],Lm1,Lm2,m[3]={0,0,0},C2[3]; 
for(i=0;i<nb_pts1;i++) {m[0]+=pts1[3*i]; m[1]+=pts1[3*i+1]; m[2]+=pts1[3*i+2];} for(i=0;i<nb_pts2;i++) {m[0]+=pts2[3*i]; m[1]+=pts2[3*i+1]; m[2]+=pts2[3*i+2];} m[0]=m[0]/nb_pts; m[1]=m[1]/nb_pts; m[2]=m[2]/nb_pts;
memcpy(C,m,3*sizeof(P_float));
while(d>thresh)
    {
    L1[0]=0; L1[1]=0; L1[2]=0; Lm1=0; L2[0]=0; L2[1]=0; L2[2]=0; Lm2=0;
    for(i=0;i<nb_pts1;i++) {li=dist3D(C,pts1+3*i); if(li!=0) {L1[0]+=(C[0]-pts1[3*i])/li; L1[1]+=(C[1]-pts1[3*i+1])/li;    L1[2]+=(C[2]-pts1[3*i+2])/li;     Lm1+=li; }}
    for(i=0;i<nb_pts2;i++) {li=dist3D(C,pts2+3*i); if(li!=0) {L2[0]+=(C[0]-pts2[3*i])/li; L2[1]+=(C[1]-pts2[3*i+1])/li;    L2[2]+=(C[2]-pts2[3*i+2])/li;     Lm2+=li; }}
    C2[0]=m[0] + Lm1*L1[0]/(nb_pts1*nb_pts) + Lm2*L2[0]/(nb_pts2*nb_pts); C2[1]=m[1] + Lm1*L1[1]/(nb_pts1*nb_pts) + Lm2*L2[1]/(nb_pts2*nb_pts); C2[2]=m[2] + Lm1*L1[2]/(nb_pts1*nb_pts) + Lm2*L2[2]/(nb_pts2*nb_pts);
    d=dist3D(C,C2);
    memcpy(C,C2,3*sizeof(P_float));
    count++; if(count==10000) d=0;
    }
return abs(Lm1/(P_float)nb_pts1 - Lm2/(P_float)nb_pts2);
}





void Barycenter(P_float e[3],P_float Pproj[3], P_float P1[3],P_float P2[3], P_float P3[3]) //barycenter of 3 points
{
P_float M[3][3];
M[0][0]=P1[0]; M[0][1]=P2[0]; M[0][2]=P3[0];
M[1][0]=P1[1]; M[1][1]=P2[1]; M[1][2]=P3[1];
M[2][0]=P1[2]; M[2][1]=P2[2]; M[2][2]=P3[2];
vtkMath::LinearSolve3x3(M,Pproj,e);
}

void Barycenter(P_float e[4],P_float P[4], P_float P1[3],P_float P2[3], P_float P3[3],P_float P4[3]) //barycenter of 4 points
{
P_float PP1[3]={P1[0]-P[0],P1[1]-P[1],P1[2]-P[2]},PP4[3]={P4[0]-P[0],P4[1]-P[1],P4[2]-P[2]},P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]},P1P3[3]={P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]};
P_float n[3]; crossproduct(n,P1P2,P1P3);
P_float proj[3];
P_float dp=dotproduct(PP4,n);
if(dp!=0) {P_float alpha=dotproduct(PP1,n)/dp; proj[0]=P[0]+alpha*PP4[0]; proj[1]=P[1]+alpha*PP4[1]; proj[2]=P[2]+alpha*PP4[2]; e[3]=alpha/(alpha-1.);}
else {memcpy(proj,P,3*sizeof(P_float)); e[3]=0;}
P_float eproj[3]; Barycenter(eproj,proj,P1,P2,P3);
e[0]=eproj[0]*(1-e[3]);
e[1]=eproj[1]*(1-e[3]);
e[2]=eproj[2]*(1-e[3]);
}

void GetCurvatureTriangle(P_float Kg_Kh_S[3],P_float n[3],P_float p[3],int nbneighb,P_float *pn,P_float Ncheck[3])
{
// return the normal, the 1-ring voronoi surface, and the curvatures [Desbrun00]
int i;
P_float *p1,*pi,*p2,u[3],v1[3],w1[3],v2[3],w2[3],UU,cotab,theta=2*PI,temp[3],x,y,s=0;
n[0]=0; n[1]=0; n[2]=0;
for(i=0;i<nbneighb;i++)
	{
	p1=(i==0)?(pn+3*nbneighb-3):(pn+3*i-3); pi=pn+3*i; p2=(i==nbneighb-1)?(pn):(pn+3*i+3);
	u[0]=pi[0]-p[0]; u[1]=pi[1]-p[1]; u[2]=pi[2]-p[2]; UU=dotproduct(u,u);
	v1[0]=p1[0]-p[0]; v1[1]=p1[1]-p[1]; v1[2]=p1[2]-p[2]; 	w1[0]=pi[0]-p1[0]; w1[1]=pi[1]-p1[1]; w1[2]=pi[2]-p1[2]; 
	v2[0]=p2[0]-p[0]; v2[1]=p2[1]-p[1]; v2[2]=p2[2]-p[2]; 	w2[0]=pi[0]-p2[0]; w2[1]=pi[1]-p2[1]; w2[2]=pi[2]-p2[2]; 
	// theta= atan(||u n v2||/abs(u.v2))
	y=abs(dotproduct(u,v2)); if(y==0) theta-=PI/2.;	else {crossproduct(temp,u,v2); x=norm(temp); theta-=atan2(x,y);}
	if(theta<0)
		theta=0;
	// cota= abs(v1.w1)/||v1 n w1||  cotb= abs(v2.w2)/||v2 n w2||
	x=abs(dotproduct(v1,w1)); crossproduct(temp,v1,w1); y=norm(temp); if(y==0) y=1E-10; cotab=x/y;
	x=abs(dotproduct(v2,w2)); crossproduct(temp,v2,w2); y=norm(temp); if(y==0) y=1E-10; cotab+=x/y;
	// results
	s+=cotab*UU;
	n[0]-=cotab*u[0]; n[1]-=cotab*u[1];	n[2]-=cotab*u[2];
	}
s=s/8.;
x=norm(n); if(x==0) memcpy(n,Ncheck,3*sizeof(P_float)); else {if(dotproduct(n,Ncheck)<0) {x=-x; theta=-theta;} n[0]=n[0]/x;n[1]=n[1]/x;n[2]=n[2]/x;}

Kg_Kh_S[0]=theta/s;
Kg_Kh_S[1]=x/(4.*s);
Kg_Kh_S[2]=s;
}


void GetPrincipalDirections(P_float pm[3],P_float axis[3][3],int nb_points,P_float *points,bool normalize)	// get OBB 
{
int i,j,p01,p02,p11,p12,p21,p22;
P_float dmax,d,u[3],dp,l[3];

// principal axis
	//1st axis
dmax=0;
for(i=0;i<nb_points;i++) 
	for(j=i;j<nb_points;j++)
		{
		u[0]=points[3*i]-points[3*j]; 
		u[1]=points[3*i+1]-points[3*j+1];
		u[2]=points[3*i+2]-points[3*j+2];
		d=dotproduct(u,u);

		if(d>dmax) 
			{
				dmax=d; 
				p01=i; 
				p02=j;
			}
		}

axis[0][0]=points[3*p01]-points[3*p02]; axis[0][1]=points[3*p01+1]-points[3*p02+1]; axis[0][2]=points[3*p01+2]-points[3*p02+2]; l[0]=dotproduct(axis[0],axis[0]);
pm[0]=points[3*p02]+axis[0][0]/2.; pm[1]=points[3*p02+1]+axis[0][1]/2.; pm[2]=points[3*p02+2]+axis[0][2]/2.;

	//2d axis
// projection on plane
dmax=0;
for(i=0;i<nb_points;i++) 
	for(j=i+1;j<nb_points;j++)
		{
		u[0]=points[3*i]-points[3*j]; u[1]=points[3*i+1]-points[3*j+1]; u[2]=points[3*i+2]-points[3*j+2];
		dp=dotproduct(u,axis[0]); u[0]=u[0]-dp*axis[0][0]/l[0]; u[1]=u[1]-dp*axis[0][1]/l[0]; u[2]=u[2]-dp*axis[0][2]/l[0]; 
		d=dotproduct(u,u); if(d>dmax) {dmax=d; p11=i; p12=j;}
		}
axis[1][0]=points[3*p11]-points[3*p12]; axis[1][1]=points[3*p11+1]-points[3*p12+1]; axis[1][2]=points[3*p11+2]-points[3*p12+2];
dp=dotproduct(axis[0],axis[1]); axis[1][0]=axis[1][0]-dp*axis[0][0]/l[0]; axis[1][1]=axis[1][1]-dp*axis[0][1]/l[0]; axis[1][2]=axis[1][2]-dp*axis[0][2]/l[0]; l[1]=dotproduct(axis[1],axis[1]);
u[0]=points[3*p12]-pm[0];u[1]=points[3*p12+1]-pm[1];u[2]=points[3*p12+2]-pm[2]; dp=dotproduct(u,axis[0]); u[0]=u[0]-dp*axis[0][0]/l[0]; u[1]=u[1]-dp*axis[0][1]/l[0]; u[2]=u[2]-dp*axis[0][2]/l[0];
pm[0]+=u[0]+axis[1][0]/2.; pm[1]+=u[1]+axis[1][1]/2.; pm[2]+=u[2]+axis[1][2]/2.;

	//3rd axis
// projection on line
crossproduct(axis[2],axis[0],axis[1]); l[2]=norm(axis[2]); axis[2][0]=axis[2][0]/l[2]; axis[2][1]=axis[2][1]/l[2]; axis[2][2]=axis[2][2]/l[2];

dmax=0; 
for(i=0;i<nb_points;i++) 
	for(j=i+1;j<nb_points;j++)
		{
		u[0]=points[3*i]-points[3*j]; u[1]=points[3*i+1]-points[3*j+1]; u[2]=points[3*i+2]-points[3*j+2];
		d=dotproduct(u,axis[2]);
		if(d>dmax) {dmax=d; p21=i; p22=j;}
		if(-d>dmax) {dmax=-d; p21=j; p22=i;}
		}
axis[2][0]=axis[2][0]*dmax; axis[2][1]=axis[2][1]*dmax; axis[2][2]=axis[2][2]*dmax;

u[0]=points[3*p22]-pm[0];
u[1]=points[3*p22+1]-pm[1];
u[2]=points[3*p22+2]-pm[2]; 

dp=dotproduct(u,axis[2]); 
u[0]=dp*axis[2][0]; 
u[1]=dp*axis[2][1]; 
u[2]=dp*axis[2][2]; 

l[2]=dotproduct(axis[2],axis[2]); 
u[0]=u[0]/l[2]; 
u[1]=u[1]/l[2]; 
u[2]=u[2]/l[2];

pm[0]+=u[0]+axis[2][0]/2.; pm[1]+=u[1]+axis[2][1]/2.; pm[2]+=u[2]+axis[2][2]/2.;

l[0]=sqrt(l[0]); l[1]=sqrt(l[1]); l[2]=sqrt(l[2]);

if(normalize)
	{
	axis[0][0]=axis[0][0]/l[0]; axis[0][1]=axis[0][1]/l[0]; axis[0][2]=axis[0][2]/l[0];
	axis[1][0]=axis[1][0]/l[1]; axis[1][1]=axis[1][1]/l[1]; axis[1][2]=axis[1][2]/l[1];
	axis[2][0]=axis[2][0]/l[2]; axis[2][1]=axis[2][1]/l[2]; axis[2][2]=axis[2][2]/l[2];
	}


/*
int i;
P_float u[3],v[3],w[3],x,y,dp;

// center of gravity
pm[0]=0; pm[1]=0; pm[2]=0;
for(i=0;i<nb_points;i++){pm[0]+=points[3*i];pm[1]+=points[3*i+1];pm[2]+=points[3*i+2];}
pm[0]=pm[0]/nb_points;pm[1]=pm[1]/nb_points;pm[2]=pm[2]/nb_points;

// principal axis
	//1st axis
axis[0][0]=0; axis[0][1]=0; axis[0][2]=0;
for(i=0;i<nb_points;i++)  {axis[0][0]+=(points[3*i]-pm[0])*(points[3*i]-pm[0]); axis[0][1]+=(points[3*i+1]-pm[1])*(points[3*i+1]-pm[1]); axis[0][2]+=(points[3*i+2]-pm[2])*(points[3*i+2]-pm[2]);}
dp=norm(axis[0]); axis[0][0]=axis[0][0]/dp; axis[0][1]=axis[0][1]/dp; axis[0][2]=axis[0][2]/dp;

	//2d axis
v[0]=1-pm[0];v[1]=-pm[1];v[2]=-pm[2];
dp=dotproduct(v,axis[0]); v[0]=v[0]-dp*axis[0][0]; v[1]=v[1]-dp*axis[0][1]; v[2]=v[2]-dp*axis[0][2];
dp=norm(v); v[0]=v[0]/dp; v[1]=v[1]/dp; v[2]=v[2]/dp;
crossproduct(w,axis[0],v);

x=0; y=0;
for(i=0;i<nb_points;i++)  
	{
	u[0]=points[3*i]-pm[0];u[1]=points[3*i+1]-pm[1];u[2]=points[3*i+2]-pm[2];
	dp=dotproduct(u,axis[0]); u[0]=u[0]-dp*axis[0][0]; u[1]=u[1]-dp*axis[0][1]; u[2]=u[2]-dp*axis[0][2];
	dp=dotproduct(u,v); x+=dp*dp;
	dp=dotproduct(u,w); y+=dp*dp;
	}
axis[1][0]=x*v[0]+y*w[0]; axis[1][1]=x*v[1]+y*w[1]; axis[1][2]=x*v[2]+y*w[2];
dp=norm(axis[1]); axis[1][0]=axis[1][0]/dp; axis[1][1]=axis[1][1]/dp; axis[1][2]=axis[1][2]/dp;

	//3rd axis
crossproduct(axis[2],axis[0],axis[1]);

// distances max along axis
if(!normalize)
	{
	P_float dmax=0;
	for(i=0;i<nb_points;i++)  {u[0]=points[3*i]-pm[0];u[1]=points[3*i+1]-pm[1];u[2]=points[3*i+2]-pm[2];	dp=abs(dotproduct(u,axis[0])); if(dp>dmax) dmax=dp; }
	axis[0][0]=2*axis[0][0]*dmax; axis[0][1]=2*axis[0][1]*dmax; axis[0][2]=2*axis[0][2]*dmax;

	dmax=0;
	for(i=0;i<nb_points;i++)  {u[0]=points[3*i]-pm[0];u[1]=points[3*i+1]-pm[1];u[2]=points[3*i+2]-pm[2];	dp=abs(dotproduct(u,axis[1])); if(dp>dmax) dmax=dp; }
	axis[1][0]=2*axis[1][0]*dmax; axis[1][1]=2*axis[1][1]*dmax; axis[1][2]=2*axis[1][2]*dmax;

	dmax=0;
	for(i=0;i<nb_points;i++)  {u[0]=points[3*i]-pm[0];u[1]=points[3*i+1]-pm[1];u[2]=points[3*i+2]-pm[2];	dp=abs(dotproduct(u,axis[2])); if(dp>dmax) dmax=dp; }
	axis[2][0]=2*axis[2][0]*dmax; axis[2][1]=2*axis[2][1]*dmax; axis[2][2]=2*axis[2][2]*dmax;
	}*/
}


void Reject_Histo(int nb, double *data, P_float reject_percentage)
{
if(reject_percentage==0) return; 
int i,nb_bins=100000;
int *histo=new int[nb_bins]; for(i=0;i<nb_bins;i++) histo[i]=0;
for(i=0;i<nb;i++)	histo[(int)(((P_float)nb_bins-1.)*data[i])]++;
int m=0; for(i=0;i<nb_bins;i++) if(histo[i]>histo[m]) m=i;
int index=0,count=histo[m]; while(100.-100.*(P_float)count/(P_float)nb>reject_percentage) {index++; count+=histo[m-index]; count+=histo[m+index];}
double max=(P_float)(m+index+1)/(P_float)nb_bins; double min=(P_float)(m-index-1)/(P_float)nb_bins;
for(i=0;i<nb;i++) {if(data[i]>max) data[i]=max; if(data[i]<min) data[i]=min; data[i]=(data[i]-min)/(max-min);}
}

bool GetFaceOrdering(int* l,int p1,int p2)
{
int count=1;
while(l[count]!=p1) count++;
if(count==l[0]) {if(p2==l[1]) return true; else return false;}
else {if(p2==l[count+1]) return true; else return false;}
}


int* Concatenate(int* l1,int* l2,int p1,int p2,bool checkorder) //concatenate 2 faces with ordering
{
bool reverse=false; if(checkorder)	{bool o1=GetFaceOrdering(l1,p1,p2),o2=GetFaceOrdering(l2,p1,p2); if(o1==o2) reverse=true;}

int* l=new int[l1[0]+l2[0]-3]; l[0]=l1[0]+l2[0]-4;
int count=1,index,flag;
index=0; flag=0; while(count!=l1[0]-1) {index++; if(index==l1[0]+1) index=1; if(l1[index]==p1 || l1[index]==p2) flag++; else if(flag>=2) {l[count]=l1[index]; count++;}} 
if(reverse) {index=l2[0]+1; flag=0; while(count!=l[0]+1) {index--; if(index==0) index=l2[0]; if(l2[index]==p1 || l2[index]==p2) flag++; else if(flag>=2) {l[count]=l2[index]; count++;}} }
else  {index=0; flag=0; while(count!=l[0]+1) {index++; if(index==l2[0]+1) index=1; if(l2[index]==p1 || l2[index]==p2) flag++; else if(flag>=2) {l[count]=l2[index]; count++;}} }
return l;
}


int* Concatenate(int* l1,int* l2) //concatenate 2 cells without ordering
{
int *l=new int[l1[0]+l2[0]-1]; l[0]=l1[0]+l2[0]-2;
int i,j,count=1;
int p; for(i=0;i<l1[0];i++) for(j=0;j<l2[0];j++) if(l1[i+1]==l2[j+1] || l1[i+1]==-l2[j+1]) p=l1[i+1];
for(i=0;i<l1[0];i++) if(l1[i+1]!=p && l1[i+1]!=-p) {l[count]=l1[i+1]; count++;}
for(i=0;i<l2[0];i++) if(l2[i+1]!=p && l2[i+1]!=-p) {l[count]=l2[i+1]; count++;}
return l;
}


void Separate(int* l,int p1,int a, int p2,int p3, int b, int p4,int* l1,int* l2)
{
l1[0]=3;l1[1]=p1;l1[2]=a;l1[3]=b; if(p1!=p4) {l1[4]=p4; l1[0]++;}
l2[0]=3;l2[1]=p3;l2[2]=b;l2[3]=a; if(p2!=p3) {l2[4]=p2; l2[0]++;}

int i=0,flag=-1;

for(i=0;i<l[0];i++) 
    {
    if(l[i+1]!=p1 && l[i+1]!=p2 && l[i+1]!=p3 && l[i+1]!=p4)
        {
        if(flag==p4) {l1[l1[0]+1]=l[i+1]; l1[0]++;}
        if(flag==p2) {l2[l2[0]+1]=l[i+1]; l2[0]++;}
        }
    if(l[i+1]==p4) flag=p4;
    if(l[i+1]==p2) flag=p2;
    }

i=0;
while(l[i+1]!=p1 && l[i+1]!=p2 && l[i+1]!=p3 && l[i+1]!=p4)
    {
    if(flag==p4) {l1[l1[0]+1]=l[i+1]; l1[0]++;}
    if(flag==p2) {l2[l2[0]+1]=l[i+1]; l2[0]++;}
    i++;
    }
}



void Line3D_UC(vtkStructuredPoints* vol,P_float iX0, P_float iY0, P_float iZ0, P_float iX1, P_float iY1, P_float iZ1)
{
	unsigned char* ptr=(unsigned char*)vol->GetScalarPointer();
	// starting point of line
    P_float iX = iX0, iY = iY0, iZ = iZ0;

    // direction of line
    P_float iDx = iX1-iX0, iDy = iY1-iY0, iDz = iZ1-iZ0;

    // increment or decrement depending on direction of line
    P_float iSx = (iDx > 0 ? 0.05 : (iDx < 0 ? -0.05 : 0));
    P_float iSy = (iDy > 0 ? 0.05 : (iDy < 0 ? -0.05 : 0));
    P_float iSz = (iDz > 0 ? 0.05 : (iDz < 0 ? -0.05 : 0));

    // decision parameters for voxel selection
    if ( iDx < 0 ) iDx = -iDx;
    if ( iDy < 0 ) iDy = -iDy;
    if ( iDz < 0 ) iDz = -iDz;
    P_float iAx = 2*iDx, iAy = 2*iDy, iAz = 2*iDz;
    P_float iDecX, iDecY, iDecZ;

    // determine largest direction component, single-step related variable
    P_float iMax = iDx; 
	int iVar = 0;
    if ( iDy > iMax ) { iMax = iDy; iVar = 1; }
    if ( iDz > iMax ) { iVar = 2; }

    // traverse Bresenham line
    switch ( iVar )
    {
    case 0:  // single-step in iX-direction
        iDecY = iAy - iDx;
        iDecZ = iAz - iDx;
        for (/**/; /**/; iX += iSx, iDecY += iAy, iDecZ += iAz)
        {
            // process voxel
            int ptId=vol->FindPoint(iX,iY,iZ);
			if(ptId!=-1) *(ptr+ptId)=255; 

            // take Bresenham step
            if ( iX+0.05>iX1 && iX-0.05<iX1 ) break;
            if ( iDecY >= 0 ) { iDecY -= iAx; iY += iSy; }
            if ( iDecZ >= 0 ) { iDecZ -= iAx; iZ += iSz; }
        }
        break;
    case 1:  // single-step in iY-direction
        iDecX = iAx - iDy;
        iDecZ = iAz - iDy;
        for (/**/; /**/; iY += iSy, iDecX += iAx, iDecZ += iAz)
        {
            // process voxel
            int ptId=vol->FindPoint(iX,iY,iZ);
			if(ptId!=-1) *(ptr+ptId)=255; 

            // take Bresenham step
            if ( iY+0.05>iY1 && iY-0.05<iY1 ) break;
            if ( iDecX >= 0 ) { iDecX -= iAy; iX += iSx; }
            if ( iDecZ >= 0 ) { iDecZ -= iAy; iZ += iSz; }
        }
        break;
    case 2:  // single-step in iZ-direction
        iDecX = iAx - iDz;
        iDecY = iAy - iDz;
        for (/**/; /**/; iZ += iSz, iDecX += iAx, iDecY += iAy)
        {
            // process voxel
            int ptId=vol->FindPoint(iX,iY,iZ);
			if(ptId!=-1) *(ptr+ptId)=255; 

            // take Bresenham step
            if ( iZ+0.05>iZ1 && iZ-0.05<iZ1 ) break;
            if ( iDecX >= 0 ) { iDecX -= iAz; iX += iSx; }
            if ( iDecY >= 0 ) { iDecY -= iAz; iY += iSy; }
        }
        break;
    }
}

void Line3D_Int(vtkStructuredPoints* vol,P_float iX0, P_float iY0, P_float iZ0, P_float iX1, P_float iY1, P_float iZ1)
{
	int* ptr=(int*)vol->GetScalarPointer();
	// starting point of line
    P_float iX = iX0, iY = iY0, iZ = iZ0;

    // direction of line
    P_float iDx = iX1-iX0, iDy = iY1-iY0, iDz = iZ1-iZ0;

    // increment or decrement depending on direction of line
    P_float iSx = (iDx > 0 ? 0.05 : (iDx < 0 ? -0.05 : 0));
    P_float iSy = (iDy > 0 ? 0.05 : (iDy < 0 ? -0.05 : 0));
    P_float iSz = (iDz > 0 ? 0.05 : (iDz < 0 ? -0.05 : 0));

    // decision parameters for voxel selection
    if ( iDx < 0 ) iDx = -iDx;
    if ( iDy < 0 ) iDy = -iDy;
    if ( iDz < 0 ) iDz = -iDz;
    P_float iAx = 2*iDx, iAy = 2*iDy, iAz = 2*iDz;
    P_float iDecX, iDecY, iDecZ;

    // determine largest direction component, single-step related variable
    P_float iMax = iDx; 
	int iVar = 0;
    if ( iDy > iMax ) { iMax = iDy; iVar = 1; }
    if ( iDz > iMax ) { iVar = 2; }

    // traverse Bresenham line
    switch ( iVar )
    {
    case 0:  // single-step in iX-direction
        iDecY = iAy - iDx;
        iDecZ = iAz - iDx;
        for (/**/; /**/; iX += iSx, iDecY += iAy, iDecZ += iAz)
        {
            // process voxel
            int ptId=vol->FindPoint(iX,iY,iZ);
			if(ptId!=-1) *(ptr+ptId)=255; 

            // take Bresenham step
            if ( iX+0.05>iX1 && iX-0.05<iX1 ) break;
            if ( iDecY >= 0 ) { iDecY -= iAx; iY += iSy; }
            if ( iDecZ >= 0 ) { iDecZ -= iAx; iZ += iSz; }
        }
        break;
    case 1:  // single-step in iY-direction
        iDecX = iAx - iDy;
        iDecZ = iAz - iDy;
        for (/**/; /**/; iY += iSy, iDecX += iAx, iDecZ += iAz)
        {
            // process voxel
            int ptId=vol->FindPoint(iX,iY,iZ);
			if(ptId!=-1) *(ptr+ptId)=255; 

            // take Bresenham step
            if ( iY+0.05>iY1 && iY-0.05<iY1 ) break;
            if ( iDecX >= 0 ) { iDecX -= iAy; iX += iSx; }
            if ( iDecZ >= 0 ) { iDecZ -= iAy; iZ += iSz; }
        }
        break;
    case 2:  // single-step in iZ-direction
        iDecX = iAx - iDz;
        iDecY = iAy - iDz;
        for (/**/; /**/; iZ += iSz, iDecX += iAx, iDecY += iAy)
        {
            // process voxel
            int ptId=vol->FindPoint(iX,iY,iZ);
			if(ptId!=-1) *(ptr+ptId)=255; 

            // take Bresenham step
            if ( iZ+0.05>iZ1 && iZ-0.05<iZ1 ) break;
            if ( iDecX >= 0 ) { iDecX -= iAz; iX += iSx; }
            if ( iDecY >= 0 ) { iDecY -= iAz; iY += iSy; }
        }
        break;
    }
}

void Line2D(vtkStructuredPoints* im,int iX0, int iY0,int iX1, int iY1,unsigned char* color)
{
	unsigned char* ptr=(unsigned char*)im->GetScalarPointer();
	int i,n_components=im->GetNumberOfScalarComponents();
	int x,y;
	// starting point of line
    P_float iX = iX0, iY = iY0;

	// direction of line
    P_float iDx = iX1-iX0, iDy = iY1-iY0;

    // increment or decrement depending on direction of line
    P_float iSx = (iDx > 0 ? 0.05 : (iDx < 0 ? -0.05 : 0));
    P_float iSy = (iDy > 0 ? 0.05 : (iDy < 0 ? -0.05 : 0));

    // decision parameters for voxel selection
    if ( iDx < 0 ) iDx = -iDx;
    if ( iDy < 0 ) iDy = -iDy;
    P_float iAx = 2*iDx, iAy = 2*iDy;
    P_float iDecX, iDecY;

    // determine largest direction component, single-step related variable
    P_float iMax = iDx; 
	int iVar = 0;
    if ( iDy > iMax ) { iMax = iDy; iVar = 1; }

    // traverse Bresenham line
    switch ( iVar )
    {
    case 0:  // single-step in iX-direction
        iDecY = iAy - iDx;
        for (/**/; /**/; iX += iSx, iDecY += iAy)
        {
            // process voxel
			x=((iX-floor(iX))<0.5)?floor(iX):ceil(iX); y=((iY-floor(iY))<0.5)?floor(iY):ceil(iY);
			for(i=0;i<n_components;i++)	*(ptr+n_components*((int)x+(int)y*im->GetDimensions()[0])+i)=*(color+i);

			// take Bresenham step
            if ( iX+0.05>iX1 && iX-0.05<iX1 ) break;
            if ( iDecY >= 0 ) { iDecY -= iAx; iY += iSy; }
        }
        break;
    case 1:  // single-step in iY-direction
        iDecX = iAx - iDy;
        for (/**/; /**/; iY += iSy, iDecX += iAx)
        {
            // process voxel
			x=((iX-floor(iX))<0.5)?floor(iX):ceil(iX); y=((iY-floor(iY))<0.5)?floor(iY):ceil(iY);
			for(i=0;i<n_components;i++)	*(ptr+n_components*((int)x+(int)y*im->GetDimensions()[0])+i)=*(color+i);

            // take Bresenham step
            if ( iY+0.05>iY1 && iY-0.05<iY1 ) break;
            if ( iDecX >= 0 ) { iDecX -= iAy; iX += iSx; }
        }
        break;
    }
}


void Line2D(vtkStructuredPoints* im,int iX0, int iY0,int iX1, int iY1,unsigned short* color)
{
	unsigned short* ptr=(unsigned short*)im->GetScalarPointer();
	int i,n_components=im->GetNumberOfScalarComponents();
	int x,y;
	// starting point of line
    P_float iX = iX0, iY = iY0;

	// direction of line
    P_float iDx = iX1-iX0, iDy = iY1-iY0;

    // increment or decrement depending on direction of line
    P_float iSx = (iDx > 0 ? 0.05 : (iDx < 0 ? -0.05 : 0));
    P_float iSy = (iDy > 0 ? 0.05 : (iDy < 0 ? -0.05 : 0));

    // decision parameters for voxel selection
    if ( iDx < 0 ) iDx = -iDx;
    if ( iDy < 0 ) iDy = -iDy;
    P_float iAx = 2*iDx, iAy = 2*iDy;
    P_float iDecX, iDecY;

    // determine largest direction component, single-step related variable
    P_float iMax = iDx; 
	int iVar = 0;
    if ( iDy > iMax ) { iMax = iDy; iVar = 1; }

    // traverse Bresenham line
    switch ( iVar )
    {
    case 0:  // single-step in iX-direction
        iDecY = iAy - iDx;
        for (/**/; /**/; iX += iSx, iDecY += iAy)
        {
            // process voxel
			x=((iX-floor(iX))<0.5)?floor(iX):ceil(iX); y=((iY-floor(iY))<0.5)?floor(iY):ceil(iY);
			for(i=0;i<n_components;i++)	*(ptr+n_components*((int)x+(int)y*im->GetDimensions()[0])+i)=*(color+i);

			// take Bresenham step
            if ( iX+0.05>iX1 && iX-0.05<iX1 ) break;
            if ( iDecY >= 0 ) { iDecY -= iAx; iY += iSy; }
        }
        break;
    case 1:  // single-step in iY-direction
        iDecX = iAx - iDy;
        for (/**/; /**/; iY += iSy, iDecX += iAx)
        {
            // process voxel
			x=((iX-floor(iX))<0.5)?floor(iX):ceil(iX); y=((iY-floor(iY))<0.5)?floor(iY):ceil(iY);
			for(i=0;i<n_components;i++)	*(ptr+n_components*((int)x+(int)y*im->GetDimensions()[0])+i)=*(color+i);

            // take Bresenham step
            if ( iY+0.05>iY1 && iY-0.05<iY1 ) break;
            if ( iDecX >= 0 ) { iDecX -= iAy; iX += iSx; }
        }
        break;
    }
}


vtkStructuredPoints*  distancemap(vtkStructuredPoints* input, bool converttoint)
{
double maxdist=1E10;

int i,*dim=input->GetDimensions();
unsigned char* ptr=(unsigned char*)input->GetScalarPointer(); for(i=0;i<input->GetNumberOfPoints();i++) if(*(ptr+i)==0) *(ptr+i)=1; else *(ptr+i)=0;

// compute distmap
vtkImageEuclideanDistance* dist=vtkImageEuclideanDistance::New();
dist->InitializeOn();
dist->ConsiderAnisotropyOn();
dist->SetMaximumDistance(maxdist);
dist->SetAlgorithmToSaito ();
dist->SetInput(input);
dist->SetDimensionality(3);
//dist->setNumberOfThreads(2);
dist->Update();
double* ptr_dist=(double*)dist->GetOutput()->GetScalarPointer(); 

vtkStructuredPoints* volume=vtkStructuredPoints::New();
	volume->SetOrigin(input->GetOrigin());
	volume->SetDimensions(input->GetDimensions());
	volume->SetSpacing(input->GetSpacing());
int* ptr_vol; double* ptrf_vol;
if(converttoint)
	{
	volume->SetScalarTypeToInt();
	volume->AllocateScalars();
    ptr_vol=(int*)volume->GetScalarPointer(); 
	for(i=0;i<input->GetNumberOfPoints();i++) *(ptr_vol+i)=(int)sqrt(*(ptr_dist+i));
	}
else
	{
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
    ptrf_vol=(double*)volume->GetScalarPointer(); 
	for(i=0;i<input->GetNumberOfPoints();i++) *(ptrf_vol+i)=sqrt(*(ptr_dist+i));
	}
dist->Delete();

// sign the distance
vtkImageSeedConnectivity* connect=vtkImageSeedConnectivity::New();
	connect->SetInput(input);
	connect->SetDimensionality(3);
	connect->SetInputConnectValue(1); 
	connect->SetOutputConnectedValue(1); 
	connect->SetOutputUnconnectedValue(0);  
//	connect->AddSeed(0,0,0); connect->AddSeed(dim[0]-1,0,0); connect->AddSeed(0,dim[1]-1,0); connect->AddSeed(dim[0]-1,dim[1]-1,0); 
	connect->AddSeed(0,0,dim[2]-1);	connect->AddSeed(0,dim[1]-1,dim[2]-1); connect->AddSeed(dim[0]-1,0,dim[2]-1); connect->AddSeed(dim[0]-1,dim[1]-1,dim[2]-1);
	connect->Update();
ptr=(unsigned char*)connect->GetOutput()->GetScalarPointer();

if(converttoint) {for(i=0;i<input->GetNumberOfPoints();i++) if(*(ptr+i)==0) *(ptr_vol+i)=- *(ptr_vol+i);}
else {for(i=0;i<input->GetNumberOfPoints();i++) if(*(ptr+i)==0) *(ptrf_vol+i)=- *(ptrf_vol+i);}
connect->Delete();

return volume;
}


vtkStructuredPoints* voxelize(P_float* bounds,P_float* spac, vtkPolyData* model)
{
int i,min[3]={1E10,1E10,1E10},max[3]={-1E10,-1E10,-1E10};
double dblp[3];
double spacing[3];
double tol=10;

if(spac!=NULL) memcpy(spacing,spac,3*sizeof(P_float));
else {spacing[0]=1; spacing[1]=1; spacing[2]=1;}

if(bounds==NULL)
	{
	for(i=0; i<model->GetNumberOfPoints();i++)
		{
		model->GetPoint(i,dblp);
		if(dblp[0]<min[0]) min[0]=(int)floor(dblp[0]); 	if(dblp[1]<min[1]) min[1]=(int)floor(dblp[1]);	if(dblp[2]<min[2]) min[2]=(int)floor(dblp[2]);
		if(dblp[0]>max[0]) max[0]=(int)ceil(dblp[0]); 	if(dblp[1]>max[1]) max[1]=(int)ceil(dblp[1]);	if(dblp[2]>max[2]) max[2]=(int)ceil(dblp[2]);
		}
	min[0]-=tol;min[1]-=tol;min[2]-=tol;
	max[0]+=tol;max[1]+=tol;max[2]+=tol;
	}
else
	{
	min[0]=bounds[0]; 	max[0]=bounds[1];
	min[1]=bounds[2]; 	max[1]=bounds[3];
	min[2]=bounds[4]; 	max[2]=bounds[5];
	}

vtkStructuredPoints* volume=vtkStructuredPoints::New();
	volume->SetOrigin(min[0],min[1],min[2]);
	volume->SetDimensions((max[0]-min[0])/spacing[0],(max[1]-min[1])/spacing[1],(max[2]-min[2])/spacing[2]);
	volume->SetSpacing(spacing);
	volume->SetScalarTypeToUnsignedChar();
	volume->AllocateScalars();

unsigned char* ptr_uc=(unsigned char*)volume->GetScalarPointer(); for(i=0;i<volume->GetNumberOfPoints();i++) *(ptr_uc+i)=0;

int ids[3];
P_float p1[3], p2[3], p3[3], p1p2[3], p1p3[3], p2p1[3], p2p3[3], p3p1[3], p3p2[3], ps[3], pt[3];

//		FILE* f=fopen("C:\\temp\\test.txt","rt+"); ///// debug
//		if(f==NULL) f=fopen("C:\\temp\\test.txt","wt");  ///// debug
//		else fseek(f,0,SEEK_END); ///// debug
//		CPerfTimer start; start.Start();

//voxelize triangles using 3D bresenham
for(i=0; i<model->GetNumberOfCells(); i++)
	{
	ids[0]=model->GetCell(i)->GetPointIds()->GetId(0); ids[1]=model->GetCell(i)->GetPointIds()->GetId(1); ids[2]=model->GetCell(i)->GetPointIds()->GetId(2);
	p1[0]=model->GetPoint(ids[0])[0]; p1[1]=model->GetPoint(ids[0])[1]; p1[2]=model->GetPoint(ids[0])[2]; p2[0]=model->GetPoint(ids[1])[0]; p2[1]=model->GetPoint(ids[1])[1]; p2[2]=model->GetPoint(ids[1])[2]; p3[0]=model->GetPoint(ids[2])[0]; p3[1]=model->GetPoint(ids[2])[1]; p3[2]=model->GetPoint(ids[2])[2];
	p1p2[0]=p2[0]-p1[0]; p1p2[1]=p2[1]-p1[1]; p1p2[2]=p2[2]-p1[2]; p1p3[0]=p3[0]-p1[0]; p1p3[1]=p3[1]-p1[1]; p1p3[2]=p3[2]-p1[2];
	p2p1[0]=p1[0]-p2[0]; p2p1[1]=p1[1]-p2[1]; p2p1[2]=p1[2]-p2[2]; p2p3[0]=p3[0]-p2[0]; p2p3[1]=p3[1]-p2[1]; p2p3[2]=p3[2]-p2[2];
	p3p1[0]=p1[0]-p3[0]; p3p1[1]=p1[1]-p3[1]; p3p1[2]=p1[2]-p3[2]; p3p2[0]=p2[0]-p3[0]; p3p2[1]=p2[1]-p3[1]; p3p2[2]=p2[2]-p3[2];
	for(P_float s=0.0; s<1.0; s+=0.1)
		{
		ps[0]=p1[0]+s*p1p2[0]; ps[1]=p1[1]+s*p1p2[1]; ps[2]=p1[2]+s*p1p2[2];
		pt[0]=p1[0]+s*p1p3[0]; pt[1]=p1[1]+s*p1p3[1]; pt[2]=p1[2]+s*p1p3[2];
		Line3D_UC(volume,ps[0],ps[1],ps[2],pt[0],pt[1],pt[2]); 
		ps[0]=p2[0]+s*p2p1[0]; ps[1]=p2[1]+s*p2p1[1]; ps[2]=p2[2]+s*p2p1[2];
		pt[0]=p2[0]+s*p2p3[0]; pt[1]=p2[1]+s*p2p3[1]; pt[2]=p2[2]+s*p2p3[2];
		Line3D_UC(volume,ps[0],ps[1],ps[2],pt[0],pt[1],pt[2]); 
		ps[0]=p3[0]+s*p3p1[0]; ps[1]=p3[1]+s*p3p1[1]; ps[2]=p3[2]+s*p3p1[2];
		pt[0]=p3[0]+s*p3p2[0]; pt[1]=p3[1]+s*p3p2[1]; pt[2]=p3[2]+s*p3p2[2];
		Line3D_UC(volume,ps[0],ps[1],ps[2],pt[0],pt[1],pt[2]);
		}
	}

//		start.Stop(); double duration=start.Elapsedms(); ///// debug
//		fprintf(f,"BRESENHAM: %lf (%d points)\n",duration,model->GetNumberOfPoints()); ///// debug
//		fclose(f); ///// debug

return volume;
}



vtkStructuredPoints* MergeMRI(vtkStructuredPoints** vols,int nbMRI,bool fillholes,bool normalize,bool overwrite)
{
int i,j;

// calcul volume englobant et resolution minimale
double minx=vols[0]->GetOrigin()[0],miny=vols[0]->GetOrigin()[1],minz=vols[0]->GetOrigin()[2],spacx=vols[0]->GetSpacing()[0],spacy=vols[0]->GetSpacing()[1],spacz=vols[0]->GetSpacing()[2],maxx=vols[0]->GetDimensions()[0]*spacx+minx,maxy=vols[0]->GetDimensions()[1]*spacy+miny,maxz=vols[0]->GetDimensions()[2]*spacz+minz;
for(i=1;i<nbMRI;i++) {if(vols[i]->GetOrigin()[0]<minx) minx=vols[i]->GetOrigin()[0]; if(vols[i]->GetOrigin()[1]<miny) miny=vols[i]->GetOrigin()[1]; if(vols[i]->GetOrigin()[2]<minz) minz=vols[i]->GetOrigin()[2]; if(vols[i]->GetSpacing()[0]<spacx) spacx=vols[i]->GetSpacing()[0]; 	if(vols[i]->GetSpacing()[1]<spacy) spacy=vols[i]->GetSpacing()[1]; if(vols[i]->GetSpacing()[2]<spacz) spacz=vols[i]->GetSpacing()[2]; if((vols[i]->GetDimensions()[0]*vols[i]->GetSpacing()[0]+vols[i]->GetOrigin()[0])>maxx) maxx=(vols[i]->GetDimensions()[0]*vols[i]->GetSpacing()[0]+vols[i]->GetOrigin()[0]); 	if((vols[i]->GetDimensions()[1]*vols[i]->GetSpacing()[1]+vols[i]->GetOrigin()[1])>maxy) maxy=(vols[i]->GetDimensions()[1]*vols[i]->GetSpacing()[1]+vols[i]->GetOrigin()[1]); if((vols[i]->GetDimensions()[2]*vols[i]->GetSpacing()[2]+vols[i]->GetOrigin()[2])>maxz) maxz=(vols[i]->GetDimensions()[2]*vols[i]->GetSpacing()[2]+vols[i]->GetOrigin()[2]);	}

int dimx=ceil((maxx-minx)/spacx)+1,dimy=ceil((maxy-miny)/spacy)+1,dimz=ceil((maxz-minz)/spacz)+1;

// creation du nouveau volume
vtkStructuredPoints* ret=vtkStructuredPoints::New();
	ret->SetDimensions(dimx,dimy,dimz);
	ret->SetSpacing(spacx,spacy,spacz); 
	ret->SetScalarType(VTK_UNSIGNED_SHORT);
	ret->SetNumberOfScalarComponents(1);
	ret->AllocateScalars();
	ret->SetOrigin(minx,miny,minz);
	ret->Update();
	unsigned short* vol_ptr = (unsigned short *) ret->GetScalarPointer();
for(i=0;i<dimx*dimy*dimz;i++) *(vol_ptr+i)=0; 

if(normalize)
	{
	unsigned short min,max,Min=30000,Max=0;
	for(i=0;i<nbMRI;i++) {min=GetMin(vols[i]); max=GetMax(vols[i]); if(min<Min) Min=min; if(max>Max) Max=max;}
	for(i=0;i<nbMRI;i++) Normalize(vols[i],Max,Min);
	}

// meme resolution pour chaque volume et Remplissage
vtkImageReslice** Reslice=(vtkImageReslice**)malloc(nbMRI*sizeof(vtkImageReslice*));
unsigned short* Reslice_ptr;
for(i=0;i<nbMRI;i++)
	{
	Reslice[i]=vtkImageReslice::New();
		Reslice[i]->SetInput(vols[i]);
		Reslice[i]->SetInterpolationModeToLinear();
		Reslice[i]->SetOutputSpacing(spacx,spacy,spacz);
		Reslice[i]->SetOutputExtent(0,dimx-1,0,dimy-1,0,dimz-1);
		Reslice[i]->SetOutputOrigin(minx,miny,minz);
		Reslice[i]->Update();
	}	

if(overwrite)
	for(i=0;i<nbMRI;i++)
		{
		Reslice_ptr= (unsigned short *) (Reslice[nbMRI-i-1]->GetOutput())->GetScalarPointer();
		for(j=0;j<Reslice[nbMRI-i-1]->GetOutput()->GetNumberOfPoints();j++) 
			if(*(Reslice_ptr+j)!=0) *(vol_ptr+j)=*(Reslice_ptr+j);	
		}	
else
	{
	int k,*zlim=new int[nbMRI+1]; zlim[0]=0; zlim[nbMRI]=dimz;
	P_float zmin,zmax,zmid;

	for(i=0;i<nbMRI-1;i++) {zmin=(vols[i]->GetOrigin()[2]>vols[i+1]->GetOrigin()[2])?vols[i]->GetOrigin()[2]:vols[i+1]->GetOrigin()[2]; zmax=(vols[i]->GetOrigin()[2]+vols[i]->GetSpacing()[2]*(P_float)(vols[i]->GetDimensions()[2]-1)<vols[i+1]->GetOrigin()[2]+vols[i+1]->GetSpacing()[2]*(P_float)(vols[i+1]->GetDimensions()[2]-1))?vols[i]->GetOrigin()[2]+vols[i]->GetSpacing()[2]*(P_float)(vols[i]->GetDimensions()[2]-1):vols[i+1]->GetOrigin()[2]+vols[i+1]->GetSpacing()[2]*(P_float)(vols[i+1]->GetDimensions()[2]-1); zmid=((zmin+zmax)/2. - ret->GetOrigin()[2])/ret->GetSpacing()[2]; 
	zlim[i+1]=((zmid-floor(zmid))<(ceil(zmid)-zmid))?floor(zmid):ceil(zmid);
	}

	for(i=0;i<nbMRI;i++)
		{
		Reslice_ptr= (unsigned short *) (Reslice[nbMRI-i-1]->GetOutput())->GetScalarPointer();
		for(j=zlim[nbMRI-i-1];j<zlim[nbMRI-i];j++) 
			for(k=0;k<dimx*dimy;k++) 
				if(*(Reslice_ptr+j*dimx*dimy+k)!=0) *(vol_ptr+j*dimx*dimy+k)=*(Reslice_ptr+j*dimx*dimy+k);	
		}
	free(zlim);
	}
for(i=1;i<nbMRI;i++) {if(vols[i]->GetOrigin()[0]<minx) minx=vols[i]->GetOrigin()[0]; if(vols[i]->GetOrigin()[1]<miny) miny=vols[i]->GetOrigin()[1]; if(vols[i]->GetOrigin()[2]<minz) minz=vols[i]->GetOrigin()[2]; if(vols[i]->GetSpacing()[0]<spacx) spacx=vols[i]->GetSpacing()[0]; 	if(vols[i]->GetSpacing()[1]<spacy) spacy=vols[i]->GetSpacing()[1]; if(vols[i]->GetSpacing()[2]<spacz) spacz=vols[i]->GetSpacing()[2]; if((vols[i]->GetDimensions()[0]*vols[i]->GetSpacing()[0]+vols[i]->GetOrigin()[0])>maxx) maxx=(vols[i]->GetDimensions()[0]*vols[i]->GetSpacing()[0]+vols[i]->GetOrigin()[0]); 	if((vols[i]->GetDimensions()[1]*vols[i]->GetSpacing()[1]+vols[i]->GetOrigin()[1])>maxy) maxy=(vols[i]->GetDimensions()[1]*vols[i]->GetSpacing()[1]+vols[i]->GetOrigin()[1]); if((vols[i]->GetDimensions()[2]*vols[i]->GetSpacing()[2]+vols[i]->GetOrigin()[2])>maxz) maxz=(vols[i]->GetDimensions()[2]*vols[i]->GetSpacing()[2]+vols[i]->GetOrigin()[2]);	}
		
for(i=0;i<nbMRI;i++) Reslice[i]->Delete();

// fill holes
bool flag; int x,y,z;
if(fillholes)
	for(z=1;z<dimz;z++) 
		{ 
		flag=false; 
		for(y=0;y<dimy;y++) for(x=0;x<dimx;x++) { if(*(vol_ptr+z*dimx*dimy+y*dimx+x)!=0) flag=true;} 
		if(!flag) memcpy(vol_ptr+z*dimx*dimy,vol_ptr+(z-1)*dimx*dimy,dimx*dimy*sizeof(unsigned short));
		}
return ret;
}



bool RotateVolume(vtkStructuredPoints* vol,int x,int y)
{
unsigned short* pt=(unsigned short*)vol->GetScalarPointer();

int z;
if(x==1)  {if(y==2) z=3;  else if(y==-2) z=-3; else if(y==3) z=-2; else if(y==-3) z=2; else return false;}
else if(x==-1)  {if(y==2) z=-3;  else if(y==-2) z=3; else if(y==3) z=2; else if(y==-3) z=-2; else return false;}
else if(x==2)  {if(y==1) z=-3;  else if(y==-1) z=3; else if(y==3) z=1; else if(y==-3) z=-1; else return false;}
else if(x==-2)  {if(y==1) z=3;  else if(y==-1) z=-3; else if(y==3) z=-1; else if(y==-3) z=1; else return false;}
else if(x==3)  {if(y==1) z=2;  else if(y==-1) z=-2; else if(y==2) z=-1; else if(y==-2) z=1; else return false;}
else if(x==-3)  {if(y==1) z=-2;  else if(y==-1) z=2; else if(y==2) z=1; else if(y==-2) z=-1; else return false;}
else return false;

// rotate
int *dim=vol->GetDimensions(),dim2[3];
P_float *spac=vol->GetSpacing(),spac2[3];
P_float *org=vol->GetOrigin(),org2[3];

if(x==1) {dim2[0]=dim[0]; org2[0]=org[0]; spac2[0]=spac[0];} else if(x==-1) {dim2[0]=dim[0]; org2[0]=org[0]+(P_float)dim[0]*spac[0]; spac2[0]=spac[0];}
else if(x==2) {dim2[1]=dim[0]; org2[1]=org[1]; spac2[1]=spac[0];} else if(x==-2) {dim2[1]=dim[0]; org2[1]=org[1]; spac2[1]=spac[0];}
else if(x==3) {dim2[2]=dim[0]; org2[2]=org[2]; spac2[2]=spac[0];} else if(x==-3) {dim2[2]=dim[0]; org2[2]=org[2]; spac2[2]=spac[0];}
if(y==1) {dim2[0]=dim[1]; org2[0]=org[0]; spac2[0]=spac[1];} else if(y==-1) {dim2[0]=dim[1]; org2[0]=org[0]; spac2[0]=spac[1];}
else if(y==2) {dim2[1]=dim[1]; org2[1]=org[1]; spac2[1]=spac[1];} else if(y==-2) {dim2[1]=dim[1]; org2[1]=org[1]+(P_float)dim[1]*spac[1]; spac2[1]=spac[1];}
else if(y==3) {dim2[2]=dim[1]; org2[2]=org[2]; spac2[2]=spac[1];} else if(y==-3) {dim2[2]=dim[1]; org2[2]=org[2]; spac2[2]=spac[1];}
if(z==1) {dim2[0]=dim[2]; org2[0]=org[0]; spac2[0]=spac[2];} else if(z==-1) {dim2[0]=dim[2]; org2[0]=org[0]; spac2[0]=spac[2];}
else if(z==2) {dim2[1]=dim[2]; org2[1]=org[1]; spac2[1]=spac[2];} else if(z==-2) {dim2[1]=dim[2]; org2[1]=org[1]; spac2[1]=spac[2];}
else if(z==3) {dim2[2]=dim[2]; org2[2]=org[2]; spac2[2]=spac[2];} else if(z==-3) {dim2[2]=dim[2]; org2[2]=org[2]+(P_float)dim[2]*spac[2]; spac2[2]=spac[2];}

vtkStructuredPoints* vol2=vtkStructuredPoints::New();
	vol2->SetDimensions(dim2);
	vol2->SetOrigin(org2);
	vol2->SetSpacing(spac2);
	vol2->SetScalarTypeToUnsignedShort();
	vol2->SetNumberOfScalarComponents(1);
	vol2->AllocateScalars();
unsigned short* pt2=(unsigned short*)vol2->GetScalarPointer();

int i,j,k;
for(int k2=0;k2<dim2[2];k2++)
	for(int j2=0;j2<dim2[1];j2++)
		for(int i2=0;i2<dim2[0];i2++)
			{
			if(x==1) i=i2; else if(x==-1) i=dim2[0]-1-i2;
			else if(x==2) i=j2; else if(x==-2) i=dim2[1]-1-j2;
			else if(x==3) i=k2; else if(x==-3) i=dim2[2]-1-k2;
			if(y==1) j=i2; else if(y==-1) j=dim2[0]-1-i2;
			else if(y==2) j=j2; else if(y==-2) j=dim2[1]-1-j2;
			else if(y==3) j=k2; else if(y==-3) j=dim2[2]-1-k2;
			if(z==1) k=i2; else if(z==-1) k=dim2[0]-1-i2;
			else if(z==2) k=j2; else if(z==-2) k=dim2[1]-1-j2;
			else if(z==3) k=k2; else if(z==-3) k=dim2[2]-1-k2;
			*(pt2+k2*dim2[0]*dim2[1]+j2*dim2[0]+i2)=*(pt+k*dim[0]*dim[1]+j*dim[0]+i);
			}
vol->DeepCopy(vol2);
vol2->Delete();
return true;
}

P_float* GetGlobalMeanStddev(vtkStructuredPoints* input)
{
unsigned short val,*pt=(unsigned short*)input->GetScalarPointer();
P_float mean=0,stddev=0,*res=new P_float[2],nb=input->GetNumberOfPoints();
for(int i=0;i<nb;i++)
	{
	val=*(pt+i);
	mean+=(P_float)val; stddev+=(P_float)val*(P_float)val;
	}
stddev=sqrt((stddev-mean*mean/nb)/(nb-1));  mean=mean/nb;
res[0]=mean; res[1]=stddev;
return res;
}

P_float* GetGlobalMeanStddev(vtkStructuredPoints* input,int nbtrials)
{
unsigned short val,*pt=(unsigned short*)input->GetScalarPointer();
P_float mean=0,stddev=0,*res=new P_float[2],nb=input->GetNumberOfPoints();
int j;
for(int i=0;i<nbtrials;i++)
	{
	j=rand()*(nb-1)/RAND_MAX; val=*(pt+j);
	mean+=(P_float)val; stddev+=(P_float)val*(P_float)val;
	}
stddev=sqrt((stddev-mean*mean/nbtrials)/(nbtrials-1));  mean=mean/nbtrials;
res[0]=mean; res[1]=stddev;
return res;
}

P_float* GetLocalMeanStddev(vtkStructuredPoints* input,int nbdim,int windowsize)
{
unsigned short val,*pt=(unsigned short*)input->GetScalarPointer();
P_float mean,stddev,*res=new P_float[2*input->GetNumberOfPoints()],nb;
int pos[3],x,y,z,X,Y,Z,*dim=input->GetDimensions();

if(nbdim==3)
	{
	nb=(2*windowsize+1)*(2*windowsize+1)*(2*windowsize+1);
	for(Z=0;Z<dim[2];Z++) for(Y=0;Y<dim[1];Y++) for(X=0;X<dim[0];X++)
		{
		mean=0; stddev=0;
		for(z=-windowsize;z<=windowsize;z++) for(y=-windowsize;y<=windowsize;y++) for(x=-windowsize;x<=windowsize;x++)
			{
			pos[0]=X+x; if(pos[0]<0) pos[0]=0; if(pos[0]>dim[0]-1) pos[0]=dim[0]-1; pos[1]=Y+y; if(pos[1]<0) pos[1]=0; if(pos[1]>dim[1]-1) pos[1]=dim[1]-1; pos[2]=Z+z; if(pos[2]<0) pos[2]=0; if(pos[2]>dim[2]-1) pos[2]=dim[2]-1;
			val=*(pt+pos[0]+pos[1]*dim[0]+pos[2]*dim[0]*dim[1]); mean+=(P_float)val; stddev+=(P_float)val*(P_float)val;
			}
		stddev=sqrt((stddev-mean*mean/nb)/(nb-1));  mean=mean/nb;
		res[2*(X+Y*dim[0]+Z*dim[0]*dim[1])]=mean;	res[2*(X+Y*dim[0]+Z*dim[0]*dim[1])+1]=stddev;
		}
	}
else if(nbdim==2)
	{
	nb=(2*windowsize+1)*(2*windowsize+1);
	for(Z=0;Z<dim[2];Z++) for(Y=0;Y<dim[1];Y++) for(X=0;X<dim[0];X++)
		{
		mean=0; stddev=0;
		for(y=-windowsize;y<=windowsize;y++) for(x=-windowsize;x<=windowsize;x++)
			{
			pos[0]=X+x; if(pos[0]<0) pos[0]=0; if(pos[0]>dim[0]-1) pos[0]=dim[0]-1; pos[1]=Y+y; if(pos[1]<0) pos[1]=0; if(pos[1]>dim[1]-1) pos[1]=dim[1]-1; 
			val=*(pt+pos[0]+pos[1]*dim[0]+Z*dim[0]*dim[1]); mean+=(P_float)val; stddev+=(P_float)val*(P_float)val;
			}
		stddev=sqrt((stddev-mean*mean/nb)/(nb-1));  mean=mean/nb;
		res[2*(X+Y*dim[0]+Z*dim[0]*dim[1])]=mean;	res[2*(X+Y*dim[0]+Z*dim[0]*dim[1])+1]=stddev;
		}
	}
return res;
}



P_float GetNoiseStddev(vtkStructuredPoints* input,P_float *MeanStddev_Vol,P_float treshold) // here MeanStddev_Vol is precomputed
{
// compute noise deviation (sigma=1/2(mean(v^2))) from background points (mean<treshold)
unsigned short val,*pt=(unsigned short*)input->GetScalarPointer();
P_float stddev=0;
int count=0;
for(int i=0;i<input->GetNumberOfPoints();i++)	if(MeanStddev_Vol[2*i]<treshold && MeanStddev_Vol[2*i]!=0) {val=*(pt+i); stddev+=val*val; count++;}
stddev=sqrt(stddev/(2.*(P_float)count));
return stddev;
}

P_float GetNoiseStddev(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold)
{
// compute noise deviation (sigma=1/2(mean(v^2))) from background points (mean<treshold)
unsigned short val,*pt=(unsigned short*)input->GetScalarPointer();
P_float mean,stddev=0,nb;
int count=0,x,y,z,X,Y,Z,*dim=input->GetDimensions();
if(nbdim==3)
	{
	nb=(2*windowsize+1)*(2*windowsize+1)*(2*windowsize+1);
	P_float sumtreshold=nb*treshold;
	for(Z=windowsize;Z<dim[2]-windowsize;Z++)	for(Y=windowsize;Y<dim[1]-windowsize;Y++)		for(X=windowsize;X<dim[0]-windowsize;X++)
		{
		mean=0; 
		for(z=-windowsize;z<=windowsize;z++)	for(y=-windowsize;y<=windowsize;y++)		for(x=-windowsize;x<=windowsize;x++)	if(mean<sumtreshold) mean+=(P_float)*(pt+(x+X)+(y+Y)*dim[0]+(z+Z)*dim[0]*dim[1]); 
		if(mean<sumtreshold && mean!=0) {val=*(pt+X+Y*dim[0]+Z*dim[0]*dim[1]); stddev+=val*val; count++;}
		}
	}
else if(nbdim==2)
	{
	nb=(2*windowsize+1)*(2*windowsize+1);
	P_float sumtreshold=nb*treshold;
	for(Z=0;Z<dim[2]-1;Z++)	for(Y=windowsize;Y<dim[1]-windowsize;Y++)		for(X=windowsize;X<dim[0]-windowsize;X++)
		{
		mean=0; 
		for(y=-windowsize;y<=windowsize;y++)		for(x=-windowsize;x<=windowsize;x++)	if(mean<sumtreshold) mean+=(P_float)*(pt+(x+X)+(y+Y)*dim[0]+Z*dim[0]*dim[1]); 
		if(mean<sumtreshold && mean!=0) {val=*(pt+X+Y*dim[0]+Z*dim[0]*dim[1]); stddev+=val*val; count++;}
		}
	}
stddev=sqrt(stddev/(2.*(P_float)count));
return stddev;
}


P_float GetNoiseStddev(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold,int nbtrials)
{
// compute noise deviation (sigma=1/2(mean(v^2))) from background points (mean<treshold)
unsigned short val,*pt=(unsigned short*)input->GetScalarPointer();
P_float mean,stddev=0,nb;
int count=0,x,y,z,X,Y,Z,*dim=input->GetDimensions();
if(nbdim==3)
	{
	nb=(2*windowsize+1)*(2*windowsize+1)*(2*windowsize+1);
	P_float sumtreshold=nb*treshold;
	while(count<nbtrials)
		{
		Z=rand()*(dim[2]-2*windowsize-1)/RAND_MAX+windowsize; 	Y=rand()*(dim[1]-2*windowsize-1)/RAND_MAX+windowsize; 	X=rand()*(dim[0]-2*windowsize-1)/RAND_MAX+windowsize; 
		mean=0; 
		for(z=-windowsize;z<=windowsize;z++)	for(y=-windowsize;y<=windowsize;y++)		for(x=-windowsize;x<=windowsize;x++)	if(mean<sumtreshold) mean+=(P_float)*(pt+(x+X)+(y+Y)*dim[0]+(z+Z)*dim[0]*dim[1]); 
		if(mean<sumtreshold && mean!=0) {val=*(pt+X+Y*dim[0]+Z*dim[0]*dim[1]); stddev+=val*val; count++;}
		}
	}
else if(nbdim==2)
	{
	nb=(2*windowsize+1)*(2*windowsize+1);
	P_float sumtreshold=nb*treshold;
	while(count<nbtrials)
		{
		Z=rand()*(dim[2]-1)/RAND_MAX; 	Y=rand()*(dim[1]-2*windowsize-1)/RAND_MAX+windowsize; 	X=rand()*(dim[0]-2*windowsize-1)/RAND_MAX+windowsize; 
		mean=0; 
		for(y=-windowsize;y<=windowsize;y++)		for(x=-windowsize;x<=windowsize;x++)	if(mean<sumtreshold) mean+=(P_float)*(pt+(x+X)+(y+Y)*dim[0]+Z*dim[0]*dim[1]); 
		if(mean<sumtreshold && mean!=0) {val=*(pt+X+Y*dim[0]+Z*dim[0]*dim[1]); stddev+=val*val; count++;}
		}
	}
stddev=sqrt(stddev/(2.*(P_float)count));
return stddev;
}

P_float GetSNR(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold)
{
P_float *meanstd=GetGlobalMeanStddev(input);
P_float noisestd=GetNoiseStddev(input,nbdim,windowsize,treshold);
P_float SNR=meanstd[1]/noisestd;
return SNR;
}


P_float GetSNR(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold,int nbtrials)
{
P_float *meanstd=GetGlobalMeanStddev(input,nbtrials);
P_float noisestd=GetNoiseStddev(input,nbdim,windowsize,treshold,nbtrials);
P_float SNR=meanstd[1]/noisestd;
return SNR;
}


P_float NLDenoising(vtkStructuredPoints* input,int nbdim,int neighborsize,int maxdist,int nbtrials,P_float h)
{
P_float* MeanStddev=GetLocalMeanStddev(input,nbdim,neighborsize);
P_float noisestd=GetNoiseStddev(input,MeanStddev,300);
P_float *globalmeanstd=GetGlobalMeanStddev(input); P_float SNR1=globalmeanstd[1]/noisestd; free(globalmeanstd); // input SNR

int X,Y,Z,x,y,z,i,j,i2,j2,k,*dim=input->GetDimensions();
P_float H=h*noisestd;
P_float dist,weight,sum,val,h2=H*H,testm,testsdev,m_tresh[2]={.8,1.2},sdev_tresh[2]={0.5,1.5};
//P_float dist,weight,sum,val,h2=H*H,testm,testsdev,m_tresh[2]={.95,1.05},sdev_tresh[2]={0.7,1.22};

vtkImageMirrorPad* pad=vtkImageMirrorPad::New();
	pad->SetInput(input); 	
	pad->SetOutputWholeExtent(-neighborsize,dim[0]-1+neighborsize,-neighborsize,dim[1]-1+neighborsize,-neighborsize,dim[2]-1+neighborsize);	
	pad->Update(); 	
vtkStructuredPoints* vol=vtkStructuredPoints::New();
	vol->DeepCopy(pad->GetOutput());
pad->Delete();
unsigned short *ptout=(unsigned short*)input->GetScalarPointer(),*pt=(unsigned short*)vol->GetScalarPointer();
int *dimvol=vol->GetDimensions();

if(nbdim==3)
	{
	for(Z=neighborsize;Z<dimvol[2]-neighborsize;Z++)	for(Y=neighborsize;Y<dimvol[1]-neighborsize;Y++)		for(X=neighborsize;X<dimvol[0]-neighborsize;X++)
		{
		i=X+Y*dimvol[0]+Z*dimvol[0]*dimvol[1]; i2=(X-neighborsize)+(Y-neighborsize)*dim[0]+(Z-neighborsize)*dim[0]*dim[1]; val=*(pt+i); sum=1;
		for(k=0;k<nbtrials;k++)
			{
			x=0; while(x<neighborsize || x>dimvol[0]-neighborsize-1) {x=X+rand()*(2*maxdist)/RAND_MAX-maxdist;}		y=0; while(y<neighborsize || y>dimvol[1]-neighborsize-1) {y=Y+rand()*(2*maxdist)/RAND_MAX-maxdist;}		z=0; while(z<neighborsize || z>dimvol[2]-neighborsize-1) {z=Z+rand()*(2*maxdist)/RAND_MAX-maxdist;}
			j=x+y*dimvol[0]+z*dimvol[0]*dimvol[1]; j2=(x-neighborsize)+(y-neighborsize)*dim[0]+(z-neighborsize)*dim[0]*dim[1];
			testm=*(MeanStddev+2*i2)/ *(MeanStddev+2*j2); testsdev=*(MeanStddev+2*i2+1)/ *(MeanStddev+2*j2+1);
			if(testm>m_tresh[0] && testm<m_tresh[1])
				if(testsdev>sdev_tresh[0] && testsdev<sdev_tresh[1])
					{		
					dist=EuclideanDistance_Vol(pt,dimvol,i,j,neighborsize);
					weight=exp(-dist/h2); sum+=weight;
					val+=*(pt+j)*weight;
					}
			}
		val=val/sum;
		*(ptout+i2)=val;
		}
	}
else if(nbdim==2)
	{
	for(Z=neighborsize;Z<dimvol[2]-neighborsize;Z++)	for(Y=neighborsize;Y<dimvol[1]-neighborsize;Y++)		for(X=neighborsize;X<dimvol[0]-neighborsize;X++)
		{
		i=X+Y*dimvol[0]+Z*dimvol[0]*dimvol[1]; i2=(X-neighborsize)+(Y-neighborsize)*dim[0]+(Z-neighborsize)*dim[0]*dim[1]; val=*(pt+i); sum=1;
		for(k=0;k<nbtrials;k++)
			{
			x=0; while(x<neighborsize || x>dimvol[0]-neighborsize-1) {x=X+rand()*(2*maxdist)/RAND_MAX-maxdist;}		y=0; while(y<neighborsize || y>dimvol[1]-neighborsize-1) {y=Y+rand()*(2*maxdist)/RAND_MAX-maxdist;}		
			j=x+y*dimvol[0]+Z*dimvol[0]*dimvol[1]; j2=(x-neighborsize)+(y-neighborsize)*dim[0]+(Z-neighborsize)*dim[0]*dim[1];
			testm=*(MeanStddev+2*i2)/ *(MeanStddev+2*j2); testsdev=*(MeanStddev+2*i2+1)/ *(MeanStddev+2*j2+1);
			if(testm>m_tresh[0] && testm<m_tresh[1])
				if(testsdev>sdev_tresh[0] && testsdev<sdev_tresh[1])
					{		
					dist=EuclideanDistance_Im(pt,dimvol,i,j,neighborsize);
					weight=exp(-dist/h2); sum+=weight;
					val+=*(pt+j)*weight;
					}
			}
		val=val/sum;
		*(ptout+i2)=val;
		}
	}

vol->Delete();
noisestd=GetNoiseStddev(input,MeanStddev,300); globalmeanstd=GetGlobalMeanStddev(input); P_float SNR2=globalmeanstd[1]/noisestd; free(globalmeanstd); // output SNR
free(MeanStddev);
return SNR2-SNR1;
}

inline P_float EuclideanDistance_Vol(unsigned short* ptvol,int* dim,int i,int j,int neighborsize)
{
int offset; P_float val1,val2,ret=0,card=(2*neighborsize+1)*(2*neighborsize+1)*(2*neighborsize+1); 
for(int z=-neighborsize;z<=neighborsize;z++)	for(int y=-neighborsize;y<=neighborsize;y++)		for(int x=-neighborsize;x<=neighborsize;x++)	{offset=x+y*dim[0]+z*dim[0]*dim[1]; 	val1=(P_float)*(ptvol+i+offset); 	val2=(P_float)*(ptvol+j+offset);	ret+=(val1-val2)*(val1-val2);}
return ret;///card;
}

inline P_float EuclideanDistance_Im(unsigned short* ptvol,int* dim,int i,int j,int neighborsize)
{
int offset; P_float val1,val2,ret=0; 
for(int y=-neighborsize;y<=neighborsize;y++)		for(int x=-neighborsize;x<=neighborsize;x++)	{offset=x+y*dim[0]; 	val1=(P_float)*(ptvol+i+offset); 	val2=(P_float)*(ptvol+j+offset);	ret+=(val1-val2)*(val1-val2);}
return ret;
}

void ComputeDifference_vol(vtkStructuredPoints* input1,vtkStructuredPoints* input2) 
{
// put in input2 the RMS difference
unsigned short *pt1=(unsigned short*)input1->GetScalarPointer(),*pt2=(unsigned short*)input2->GetScalarPointer();
for(int i=0;i<input2->GetNumberOfPoints();i++)	*(pt2+i)=(unsigned short)sqrt(((double)*(pt2+i)- (double)*(pt1+i))*((double)*(pt2+i)- (double)*(pt1+i)));
}


void SmoothVol(vtkStructuredPoints* input,P_float stddev,int dim)
{
P_float Stddev[3]={stddev/input->GetSpacing()[0],stddev/input->GetSpacing()[1],stddev/input->GetSpacing()[2]};
vtkImageGaussianSmooth *gaussian=vtkImageGaussianSmooth::New();
     gaussian->SetDimensionality(dim);
     gaussian->SetStandardDeviations(Stddev[0],Stddev[1],Stddev[2]);
     gaussian->SetRadiusFactors(3*Stddev[0],3*Stddev[1],3*Stddev[2]);
     gaussian->SetInput(input);
     gaussian->Update();
 input->DeepCopy(gaussian->GetOutput());
     gaussian->Delete();
}



void SmoothVolAnisotropic(vtkStructuredPoints* input,P_float difffactor,P_float difftresh,int nbit,int dim)
{
if(dim==2)
	{
	vtkImageAnisotropicDiffusion2D *diff=vtkImageAnisotropicDiffusion2D::New();
		diff->FacesOn(); diff->EdgesOn(); diff->CornersOff();
		diff->SetInput(input);
	diff->GradientMagnitudeThresholdOff(); ///// TO TEST
		diff->SetDiffusionFactor(difffactor);
		diff->SetDiffusionThreshold(difftresh);
		diff->SetNumberOfIterations(nbit);
		diff->Update();
	input->DeepCopy(diff->GetOutput());
		diff->Delete();
	}
else
	{
	vtkImageAnisotropicDiffusion3D *diff=vtkImageAnisotropicDiffusion3D::New();
		diff->FacesOn(); diff->EdgesOn(); diff->CornersOff();
		diff->SetInput(input);
	diff->GradientMagnitudeThresholdOff(); ///// TO TEST
		diff->SetDiffusionFactor(difffactor);
		diff->SetDiffusionThreshold(difftresh);
		diff->SetNumberOfIterations(nbit);
	
		try {
			diff->Update();
		}
		catch (std::exception &e) {
			cout << e.what();
		}

	input->DeepCopy(diff->GetOutput());
		diff->Delete();
	}
}


vtkStructuredPoints* GradientMagnitude(vtkStructuredPoints* input,int dim)
{
if(input==NULL) return NULL;
vtkImageGradientMagnitude* grad=vtkImageGradientMagnitude::New();
	grad->HandleBoundariesOn ();
	grad->SetDimensionality (dim);
	grad->SetInput(input);
    grad->Update();
vtkStructuredPoints* output=vtkStructuredPoints::New();
output->DeepCopy(grad->GetOutput());
    grad->Delete();
	return output;
}

vtkStructuredPoints* Gradient(vtkStructuredPoints* vol, int nbdim, P_float** Transform_RTMRI) {
	if (vol == NULL) return NULL;

	vtkStructuredPoints* grad = vtkStructuredPoints::New();
	grad->SetDimensions(vol->GetDimensions());
	grad->SetSpacing(vol->GetSpacing());
	grad->SetOrigin(vol->GetOrigin());
	grad->SetNumberOfScalarComponents(3);
	grad->SetScalarTypeToFloat();
	grad->AllocateScalars();
	
	float* ptr = (float*)grad->GetScalarPointer();
	unsigned short* pvol = (unsigned short*)vol->GetScalarPointer();
	int* dim = vol->GetDimensions();
	int X, Y, Z; 
	if (nbdim == 3)	{
		for (int z = 0; z < dim[2]; z++)
			for (int y = 0; y < dim[1]; y++)
				for (int x = 0; x < dim[0]; x++) {
					if (x == 0) X = 1; else if (x == dim[0] - 1) X = dim[0] - 2; else X = x;	
					if (y == 0) Y = 1; else if (y == dim[1] - 1) Y = dim[1] - 2; else Y = y;	
					if (z == 0) Z = 1; else if (z == dim[2] - 1) Z = dim[2] - 2; else Z = z;	
					
					*(ptr+3*(x+y*dim[0]+z*dim[0]*dim[1]))=((float)*(pvol+(X+1)+(Y)*dim[0]+(Z)*dim[0]*dim[1])-(float)*(pvol+(X-1)+(Y)*dim[0]+(Z)*dim[0]*dim[1]))/(2.*vol->GetSpacing()[0]);
					*(ptr+3*(x+y*dim[0]+z*dim[0]*dim[1])+1)=((float)*(pvol+(X)+(Y+1)*dim[0]+(Z)*dim[0]*dim[1])-(float)*(pvol+(X)+(Y-1)*dim[0]+(Z)*dim[0]*dim[1]))/(2.*vol->GetSpacing()[1]);
					*(ptr+3*(x+y*dim[0]+z*dim[0]*dim[1])+2)=((float)*(pvol+(X)+(Y)*dim[0]+(Z+1)*dim[0]*dim[1])-(float)*(pvol+(X)+(Y)*dim[0]+(Z-1)*dim[0]*dim[1]))/(2.*vol->GetSpacing()[2]);
				}
	}

	else if (nbdim == 2) {
		P_float v[3] = {0, 0, 0}, v2[3];
		for (int z = 0; z < dim[2]; z++)
			for (int y = 0; y < dim[1]; y++)
				for (int x = 0; x < dim[0]; x++) {
					if (x == 0) X = 1; else if (x == dim[0] - 1) X = dim[0] - 2; else X = x;	
					if (y == 0) Y = 1; else if (y == dim[1] - 1) Y = dim[1] - 2; else Y = y;	
					
					Z = z;
					v[0] = ((float)*(pvol+(X+1)+(Y)*dim[0]+(Z)*dim[0]*dim[1])-(float)*(pvol+(X-1)+(Y)*dim[0]+(Z)*dim[0]*dim[1]))/(2.*vol->GetSpacing()[0]);
					v[1] = ((float)*(pvol+(X)+(Y+1)*dim[0]+(Z)*dim[0]*dim[1])-(float)*(pvol+(X)+(Y-1)*dim[0]+(Z)*dim[0]*dim[1]))/(2.*vol->GetSpacing()[1]);
					
					Transform_R(v, v2, Transform_RTMRI[z]);
					ptr[3*(x+y*dim[0]+z*dim[0]*dim[1])] = v2[0];
					ptr[3*(x+y*dim[0]+z*dim[0]*dim[1])+1] = v2[1];
					ptr[3*(x+y*dim[0]+z*dim[0]*dim[1])+2] = v2[2];
				}
	}
	return grad;
}



int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float pos[3],int interpolationmode,P_float step) {return Gradient(grad,vol,pos[0],pos[1],pos[2],interpolationmode,step);}
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step)
{
P_float i1,i2;
i1=Interpolation(vol,posx+step,posy,posz,interpolationmode); if(i1==-1) return -1;	i2=Interpolation(vol,posx-step,posy,posz,interpolationmode); if(i2==-1) return -1;
grad[0]=i1-i2;
i1=Interpolation(vol,posx,posy+step,posz,interpolationmode); if(i1==-1) return -1;	i2=Interpolation(vol,posx,posy-step,posz,interpolationmode); if(i2==-1) return -1;
grad[1]=i1-i2;
i1=Interpolation(vol,posx,posy,posz+step,interpolationmode); if(i1==-1) return -1;	i2=Interpolation(vol,posx,posy,posz-step,interpolationmode); if(i2==-1) return -1;
grad[2]=i1-i2;
return 1;
}

int Gradient_Int(P_float grad[3],vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step)
{
P_float i1,i2;
i1=Interpolation_Int(vol,posx+step,posy,posz,interpolationmode); if(i1==-1E10) return -1E10;	i2=Interpolation_Int(vol,posx-step,posy,posz,interpolationmode); if(i2==-1) return -1E10;
grad[0]=i1-i2;
i1=Interpolation_Int(vol,posx,posy+step,posz,interpolationmode); if(i1==-1E10) return -1E10;	i2=Interpolation_Int(vol,posx,posy-step,posz,interpolationmode); if(i2==-1) return -1E10;
grad[1]=i1-i2;
i1=Interpolation_Int(vol,posx,posy,posz+step,interpolationmode); if(i1==-1E10) return -1E10;	i2=Interpolation_Int(vol,posx,posy,posz-step,interpolationmode); if(i2==-1) return -1E10;
grad[2]=i1-i2;
return 1;
}

int Gradient_Dbl(P_float grad[3],vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step)
{
P_float i1,i2;
i1=Interpolation_Dbl(vol,posx+step,posy,posz,interpolationmode); if(i1==-1E10) return -1E10;	i2=Interpolation_Dbl(vol,posx-step,posy,posz,interpolationmode); if(i2==-1) return -1E10;
grad[0]=i1-i2;
i1=Interpolation_Dbl(vol,posx,posy+step,posz,interpolationmode); if(i1==-1E10) return -1E10;	i2=Interpolation_Dbl(vol,posx,posy-step,posz,interpolationmode); if(i2==-1) return -1E10;
grad[1]=i1-i2;
i1=Interpolation_Dbl(vol,posx,posy,posz+step,interpolationmode); if(i1==-1E10) return -1E10;	i2=Interpolation_Dbl(vol,posx,posy,posz-step,interpolationmode); if(i2==-1) return -1E10;
grad[2]=i1-i2;
return 1;
}
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float *transform_P,P_float *spacing_P,P_float pos[3],int interpolationmode,P_float step) {return Gradient(grad,vol,transform_P,spacing_P,pos[0],pos[1],pos[2],interpolationmode,step);}
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float *transform_P,P_float *spacing_P,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step)
{
P_float i1,i2;
i1=Interpolation(vol,posx+step,posy,posz,transform_P,spacing_P,interpolationmode); if(i1==-1) return -1;	i2=Interpolation(vol,posx-step,posy,posz,transform_P,spacing_P,interpolationmode); if(i2==-1) return -1;
grad[0]=i1-i2;
i1=Interpolation(vol,posx,posy+step,posz,transform_P,spacing_P,interpolationmode); if(i1==-1) return -1;	i2=Interpolation(vol,posx,posy-step,posz,transform_P,spacing_P,interpolationmode); if(i2==-1) return -1;
grad[1]=i1-i2;
i1=Interpolation(vol,posx,posy,posz+step,transform_P,spacing_P,interpolationmode); if(i1==-1) return -1;	i2=Interpolation(vol,posx,posy,posz-step,transform_P,spacing_P,interpolationmode); if(i2==-1) return -1;
grad[2]=i1-i2;
return 1;
}


P_float Interpolation(vtkStructuredPoints* vol,P_float pos[3],int interpolationmode) {return Interpolation(vol,pos[0],pos[1],pos[2],interpolationmode);}
P_float Interpolation(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode)
{
//	interpolationmode=0: trilinear
//	interpolationmode=1: nearest neighbor 
unsigned short* pt= (unsigned short*) vol->GetScalarPointer();

P_float pos[3]={(double)posx,(double)posy,(double)posz};

double* spac=vol->GetSpacing();
double* org=vol->GetOrigin();
int* dim=vol->GetDimensions();

pos[0]=(pos[0]-org[0])/spac[0];pos[1]=(pos[1]-org[1])/spac[1];pos[2]=(pos[2]-org[2])/spac[2];
if(pos[0]<1 || pos[1]<1 || pos[2]<1 || pos[0]>dim[0]-2 || pos[1]>dim[1]-2 || pos[2]>dim[2]-2) return -1;

P_float w[3][2];
w[0][0]=(floor(pos[0])-pos[0]); w[1][0]=(floor(pos[1])-pos[1]); w[2][0]=(floor(pos[2])-pos[2]);
if(w[0][0]<-0.5) {w[0][1]=w[0][0]; w[0][0]=1+w[0][1];} else w[0][1]=1+w[0][0];
if(w[1][0]<-0.5) {w[1][1]=w[1][0]; w[1][0]=1+w[1][1];} else w[1][1]=1+w[1][0];
if(w[2][0]<-0.5) {w[2][1]=w[2][0]; w[2][0]=1+w[2][1];} else w[2][1]=1+w[2][0];

unsigned short ng[8];
ng[0]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);

if(interpolationmode==1) return ng[0];

ng[1]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[2]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[3]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[4]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[5]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[6]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[7]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);

w[0][0]=abs(w[0][0]);w[0][1]=abs(w[0][1]);
w[1][0]=abs(w[1][0]);w[1][1]=abs(w[1][1]);
w[2][0]=abs(w[2][0]);w[2][1]=abs(w[2][1]);

P_float result=w[0][1]*w[1][1]*w[2][1]*ng[0]
			+w[0][0]*w[1][1]*w[2][1]*ng[1]
			+w[0][0]*w[1][0]*w[2][1]*ng[2]
			+w[0][1]*w[1][0]*w[2][1]*ng[3]
			+w[0][1]*w[1][1]*w[2][0]*ng[4]
			+w[0][0]*w[1][1]*w[2][0]*ng[5]
			+w[0][0]*w[1][0]*w[2][0]*ng[6]
			+w[0][1]*w[1][0]*w[2][0]*ng[7];

return result;
}


// rt mri
P_float Interpolation(vtkStructuredPoints* vol,P_float pos[3],int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode) {return Interpolation(vol,pos[0],pos[1],pos[2],pselect,transform_P,spacing_P,tol,interpolationmode);}
P_float Interpolation(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode)
{
//	interpolationmode=0: trilinear
//	interpolationmode=1: nearest neighbor 
if(pselect==-1) return -1;
unsigned short* pt= (unsigned short*) vol->GetScalarPointer();

P_float x,y,dp,OP[3],pos[3]={(double)posx,(double)posy,(double)posz};

int* dim=vol->GetDimensions();

OP[0]=pos[0]-transform_P[pselect][3]; OP[1]=pos[1]-transform_P[pselect][7];	OP[2]=pos[2]-transform_P[pselect][11];
dp=OP[0]*transform_P[pselect][2]+OP[1]*transform_P[pselect][6]+OP[2]*transform_P[pselect][10];
if(abs(dp)>tol) return -1;

OP[0]-=dp*transform_P[pselect][2]; OP[1]-=dp*transform_P[pselect][6]; OP[2]-=dp*transform_P[pselect][10];
x=(OP[0]*transform_P[pselect][0]+OP[1]*transform_P[pselect][4]+OP[2]*transform_P[pselect][8])/spacing_P[pselect][0];
y=(OP[0]*transform_P[pselect][1]+OP[1]*transform_P[pselect][5]+OP[2]*transform_P[pselect][9])/spacing_P[pselect][1];
if(x<1 || x>dim[0]-2 || y<1 || y>dim[1]-2) return -1;

P_float w[2][2];
w[0][0]=(floor(x)-x); w[1][0]=(floor(y)-y); 
if(w[0][0]<-0.5) {w[0][1]=w[0][0]; w[0][0]=1+w[0][1];} else w[0][1]=1+w[0][0];
if(w[1][0]<-0.5) {w[1][1]=w[1][0]; w[1][0]=1+w[1][1];} else w[1][1]=1+w[1][0];

P_float ng[4];
ng[0]=*(pt+(int)(x+w[0][0])+(int)(y+w[1][0])*dim[0]+pselect*dim[0]*dim[1]);

if(interpolationmode==1) return ng[0];

ng[1]=*(pt+(int)(x+w[0][1])+(int)(y+w[1][0])*dim[0]+pselect*dim[0]*dim[1]);
ng[2]=*(pt+(int)(x+w[0][1])+(int)(y+w[1][1])*dim[0]+pselect*dim[0]*dim[1]);
ng[3]=*(pt+(int)(x+w[0][0])+(int)(y+w[1][1])*dim[0]+pselect*dim[0]*dim[1]);

P_float result=abs(w[0][1])*abs(w[1][1])*ng[0]
			+abs(w[0][0])*abs(w[1][1])*ng[1]
			+abs(w[0][0])*abs(w[1][0])*ng[2]
			+abs(w[0][1])*abs(w[1][0])*ng[3];

return result;
}

void vInterpolation(P_float* v,vtkStructuredPoints* vol,P_float pos[3],int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode) {vInterpolation(v,vol,pos[0],pos[1],pos[2],pselect,transform_P,spacing_P,tol,interpolationmode);}
void vInterpolation(P_float* v,vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode)
{
//	interpolationmode=0: trilinear
//	interpolationmode=1: nearest neighbor 
int i;
int n_component=vol->GetNumberOfScalarComponents(); 
float* pt= (float*) vol->GetScalarPointer();

P_float x,y,dp,OP[3],pos[3]={(double)posx,(double)posy,(double)posz};

int* dim=vol->GetDimensions();

OP[0]=pos[0]-transform_P[pselect][3]; OP[1]=pos[1]-transform_P[pselect][7];	OP[2]=pos[2]-transform_P[pselect][11];
dp=OP[0]*transform_P[pselect][2]+OP[1]*transform_P[pselect][6]+OP[2]*transform_P[pselect][10];
if(abs(dp)>tol) {for(i=0;i<n_component;i++) v[i]=1E5; return;}

OP[0]-=dp*transform_P[pselect][2]; OP[1]-=dp*transform_P[pselect][6]; OP[2]-=dp*transform_P[pselect][10];
x=(OP[0]*transform_P[pselect][0]+OP[1]*transform_P[pselect][4]+OP[2]*transform_P[pselect][8])/spacing_P[pselect][0];
y=(OP[0]*transform_P[pselect][1]+OP[1]*transform_P[pselect][5]+OP[2]*transform_P[pselect][9])/spacing_P[pselect][1];
if(x<1 || x>dim[0]-2 || y<1 || y>dim[1]-2) {for(i=0;i<n_component;i++) v[i]=1E5; return;}

P_float w[2][2];
w[0][0]=(floor(x)-x); w[1][0]=(floor(y)-y);
if(w[0][0]<-0.5) {w[0][1]=w[0][0]; w[0][0]=1+w[0][1];} else w[0][1]=1+w[0][0];
if(w[1][0]<-0.5) {w[1][1]=w[1][0]; w[1][0]=1+w[1][1];} else w[1][1]=1+w[1][0];

P_float ng[4];
for(i=0;i<n_component;i++)
{
ng[0]=*(pt+n_component*((int)(x+w[0][0])+(int)(y+w[1][0])*dim[0]+pselect*dim[0]*dim[1])+i);

if(interpolationmode==1) v[i]=ng[0];
else
	{
	ng[1]=*(pt+n_component*((int)(x+w[0][1])+(int)(y+w[1][0])*dim[0]+pselect*dim[0]*dim[1])+i);
	ng[2]=*(pt+n_component*((int)(x+w[0][1])+(int)(y+w[1][1])*dim[0]+pselect*dim[0]*dim[1])+i);
	ng[3]=*(pt+n_component*((int)(x+w[0][0])+(int)(y+w[1][1])*dim[0]+pselect*dim[0]*dim[1])+i);

	v[i]=		abs(w[0][1])*abs(w[1][1])*ng[0]
				+abs(w[0][0])*abs(w[1][1])*ng[1]
				+abs(w[0][0])*abs(w[1][0])*ng[2]
				+abs(w[0][1])*abs(w[1][0])*ng[3];
	}
}
}


P_float Interpolation_Int(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode)
{
//	interpolationmode=0: trilinear
//	interpolationmode=1: nearest neighbor 
int* pt= (int*) vol->GetScalarPointer();

P_float pos[3]={(double)posx,(double)posy,(double)posz};

double* spac=vol->GetSpacing();
double* org=vol->GetOrigin();
int* dim=vol->GetDimensions();

pos[0]=(pos[0]-org[0])/spac[0];pos[1]=(pos[1]-org[1])/spac[1];pos[2]=(pos[2]-org[2])/spac[2];
if(pos[0]<1 || pos[1]<1 || pos[2]<1 || pos[0]>dim[0]-2 || pos[1]>dim[1]-2 || pos[2]>dim[2]-2) return -1E10;

P_float w[3][2];
w[0][0]=(floor(pos[0])-pos[0]); w[1][0]=(floor(pos[1])-pos[1]); w[2][0]=(floor(pos[2])-pos[2]);
if(w[0][0]<-0.5) {w[0][1]=w[0][0]; w[0][0]=1+w[0][1];} else w[0][1]=1+w[0][0];
if(w[1][0]<-0.5) {w[1][1]=w[1][0]; w[1][0]=1+w[1][1];} else w[1][1]=1+w[1][0];
if(w[2][0]<-0.5) {w[2][1]=w[2][0]; w[2][0]=1+w[2][1];} else w[2][1]=1+w[2][0];

P_float ng[8];
ng[0]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);

if(interpolationmode==1) return ng[0];

ng[1]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[2]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[3]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[4]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[5]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[6]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[7]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);

w[0][0]=abs(w[0][0]);w[0][1]=abs(w[0][1]);
w[1][0]=abs(w[1][0]);w[1][1]=abs(w[1][1]);
w[2][0]=abs(w[2][0]);w[2][1]=abs(w[2][1]);

P_float result=w[0][1]*w[1][1]*w[2][1]*ng[0]
			+w[0][0]*w[1][1]*w[2][1]*ng[1]
			+w[0][0]*w[1][0]*w[2][1]*ng[2]
			+w[0][1]*w[1][0]*w[2][1]*ng[3]
			+w[0][1]*w[1][1]*w[2][0]*ng[4]
			+w[0][0]*w[1][1]*w[2][0]*ng[5]
			+w[0][0]*w[1][0]*w[2][0]*ng[6]
			+w[0][1]*w[1][0]*w[2][0]*ng[7];

return result;
}


P_float Interpolation_Dbl(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode)
{
//	interpolationmode=0: trilinear
//	interpolationmode=1: nearest neighbor 
double* pt= (double*) vol->GetScalarPointer();

P_float pos[3]={(double)posx,(double)posy,(double)posz};

double* spac=vol->GetSpacing();
double* org=vol->GetOrigin();
int* dim=vol->GetDimensions();

pos[0]=(pos[0]-org[0])/spac[0];pos[1]=(pos[1]-org[1])/spac[1];pos[2]=(pos[2]-org[2])/spac[2];
if(pos[0]<1 || pos[1]<1 || pos[2]<1 || pos[0]>dim[0]-2 || pos[1]>dim[1]-2 || pos[2]>dim[2]-2) return -1E10;

P_float w[3][2];
w[0][0]=(floor(pos[0])-pos[0]); w[1][0]=(floor(pos[1])-pos[1]); w[2][0]=(floor(pos[2])-pos[2]);
if(w[0][0]<-0.5) {w[0][1]=w[0][0]; w[0][0]=1+w[0][1];} else w[0][1]=1+w[0][0];
if(w[1][0]<-0.5) {w[1][1]=w[1][0]; w[1][0]=1+w[1][1];} else w[1][1]=1+w[1][0];
if(w[2][0]<-0.5) {w[2][1]=w[2][0]; w[2][0]=1+w[2][1];} else w[2][1]=1+w[2][0];

P_float ng[8];
ng[0]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);

if(interpolationmode==1) return ng[0];

ng[1]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[2]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[3]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[4]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[5]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[6]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[7]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);

w[0][0]=abs(w[0][0]);w[0][1]=abs(w[0][1]);
w[1][0]=abs(w[1][0]);w[1][1]=abs(w[1][1]);
w[2][0]=abs(w[2][0]);w[2][1]=abs(w[2][1]);

P_float result=w[0][1]*w[1][1]*w[2][1]*ng[0]
			+w[0][0]*w[1][1]*w[2][1]*ng[1]
			+w[0][0]*w[1][0]*w[2][1]*ng[2]
			+w[0][1]*w[1][0]*w[2][1]*ng[3]
			+w[0][1]*w[1][1]*w[2][0]*ng[4]
			+w[0][0]*w[1][1]*w[2][0]*ng[5]
			+w[0][0]*w[1][0]*w[2][0]*ng[6]
			+w[0][1]*w[1][0]*w[2][0]*ng[7];

return result;
}

int GetSelectP(P_float pos[3],vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol)
{
P_float x,y,dp,OP[3],dmin=1E10;
int pselect=-1;
int* dim=vol->GetDimensions();

for(int i=0;i<dim[2];i++)
	{
	OP[0]=pos[0]-transform_P[i][3]; OP[1]=pos[1]-transform_P[i][7];	OP[2]=pos[2]-transform_P[i][11];
	dp=OP[0]*transform_P[i][2]+OP[1]*transform_P[i][6]+OP[2]*transform_P[i][10];
	if(abs(dp)<tol && abs(dp)<dmin)
		{
		OP[0]-=dp*transform_P[i][2]; OP[1]-=dp*transform_P[i][6]; OP[2]-=dp*transform_P[i][10];
		x=(OP[0]*transform_P[i][0]+OP[1]*transform_P[i][4]+OP[2]*transform_P[i][8])/spacing_P[i][0];
		y=(OP[0]*transform_P[i][1]+OP[1]*transform_P[i][5]+OP[2]*transform_P[i][9])/spacing_P[i][1];
		if(x>=1 && x<dim[0]-1 && y>=1 && y<dim[1]-1) {pselect=i; dmin=abs(dp);}
		}
	}
return pselect;
}

int GetSelectP(P_float n2[3],P_float pos[3],P_float n[3],vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol)
{
P_float x,y,dp,OP[3],dmin=1E10;
int pselect=-1;
int* dim=vol->GetDimensions();
n2[0]=0; n2[1]=0; n2[2]=0;
for(int i=0;i<dim[2];i++)
	{
	OP[0]=pos[0]-transform_P[i][3]; OP[1]=pos[1]-transform_P[i][7];	OP[2]=pos[2]-transform_P[i][11];
	dp=OP[0]*transform_P[i][2]+OP[1]*transform_P[i][6]+OP[2]*transform_P[i][10];
	if(abs(dp)<tol && abs(dp)<dmin) 
		{
		OP[0]-=dp*transform_P[i][2]; OP[1]-=dp*transform_P[i][6]; OP[2]-=dp*transform_P[i][10];
		x=(OP[0]*transform_P[i][0]+OP[1]*transform_P[i][4]+OP[2]*transform_P[i][8])/spacing_P[i][0];
		y=(OP[0]*transform_P[i][1]+OP[1]*transform_P[i][5]+OP[2]*transform_P[i][9])/spacing_P[i][1];
		if(x>=1 && x<dim[0]-1 && y>=1 && y<dim[1]-1) 
			{
			dp=n[0]*transform_P[i][2]+n[1]*transform_P[i][6]+n[2]*transform_P[i][10];
			n2[0]=n[0]-dp*transform_P[i][2]; n2[1]=n[1]-dp*transform_P[i][6];	n2[2]=n[2]-dp*transform_P[i][10];
			pselect=i; dmin=abs(dp);
			}
		}
	}
return pselect;
}


void vInterpolation(P_float *v, vtkStructuredPoints* vol, P_float pos[3], int interpolationmode) {
	vInterpolation(v, vol, pos[0], pos[1], pos[2], interpolationmode);
}

void vInterpolation(P_float *v, vtkStructuredPoints* vol, P_float posx, P_float posy, P_float posz, int interpolationmode) {
	//	interpolationmode=0: trilinear
	//	interpolationmode=1: nearest neighbor 
	int i;
	int n_component=vol->GetNumberOfScalarComponents(); 
	float* pt = (float*) vol->GetScalarPointer();

	P_float pos[3] = {
		(double)posx, (double)posy, (double)posz
	};

	double* spac = vol->GetSpacing();
	double* org = vol->GetOrigin();
	int* dim = vol->GetDimensions();

	pos[0] = (pos[0] - org[0]) / spac[0]; 
	pos[1] = (pos[1] - org[1]) / spac[1]; 
	pos[2] = (pos[2] - org[2]) / spac[2];
	
	if (pos[0] < 1 || pos[1] < 1 || pos[2] < 1 || pos[0] > dim[0] - 2 || pos[1] > dim[1] - 2 || pos[2] > dim[2] - 2) {
		for (i = 0; i < n_component; i++) 
			v[i] = 1E5; 
	
		return;
	}

	P_float w[3][2];
	w[0][0] = (floor(pos[0]) - pos[0]); 
	w[1][0] = (floor(pos[1]) - pos[1]); 
	w[2][0] = (floor(pos[2]) - pos[2]);

	if (w[0][0] < -0.5) {
		w[0][1] = w[0][0]; 
		w[0][0] = 1 + w[0][1];
	} 
	else {
		w[0][1] = 1 + w[0][0];
	}
	
	if (w[1][0] < -0.5) { 
		w[1][1] = w[1][0]; 
		w[1][0] = 1 + w[1][1];
	} 
	else {
		w[1][1] = 1 + w[1][0];
	}
	if (w[2][0] < -0.5) {
		w[2][1] = w[2][0]; 
		w[2][0] = 1 + w[2][1];
	} 
	else {
		w[2][1] = 1 + w[2][0];
	}

	P_float ng[8];
	for (i = 0; i < n_component; i++) {
		ng[0] = *(pt + n_component * ((int)(pos[0] + w[0][0]) + (int)(pos[1] + w[1][0]) * dim[0] + (int)(pos[2] + w[2][0]) * dim[0] * dim[1]) + i);

		if (interpolationmode == 1) v[i] = ng[0];
		else {
			ng[1] = *(pt + n_component * ((int)(pos[0] + w[0][1]) + (int)(pos[1] + w[1][0]) * dim[0] + (int)(pos[2] + w[2][0]) * dim[0] * dim[1]) + i);
			ng[2] = *(pt + n_component * ((int)(pos[0] + w[0][1]) + (int)(pos[1] + w[1][1]) * dim[0] + (int)(pos[2] + w[2][0]) * dim[0] * dim[1]) + i);
			ng[3] = *(pt + n_component * ((int)(pos[0] + w[0][0]) + (int)(pos[1] + w[1][1]) * dim[0] + (int)(pos[2] + w[2][0]) * dim[0] * dim[1]) + i);
			ng[4] = *(pt + n_component * ((int)(pos[0] + w[0][0]) + (int)(pos[1] + w[1][0]) * dim[0] + (int)(pos[2] + w[2][1]) * dim[0] * dim[1]) + i);
			ng[5] = *(pt + n_component * ((int)(pos[0] + w[0][1]) + (int)(pos[1] + w[1][0]) * dim[0] + (int)(pos[2] + w[2][1]) * dim[0] * dim[1]) + i);
			ng[6] = *(pt + n_component * ((int)(pos[0] + w[0][1]) + (int)(pos[1] + w[1][1]) * dim[0] + (int)(pos[2] + w[2][1]) * dim[0] * dim[1]) + i);
			ng[7] = *(pt + n_component * ((int)(pos[0] + w[0][0]) + (int)(pos[1] + w[1][1]) * dim[0] + (int)(pos[2] + w[2][1]) * dim[0] * dim[1]) + i);

			v[i] =		abs(w[0][1]) * abs(w[1][1]) * abs(w[2][1]) * ng[0]
						+ abs(w[0][0]) * abs(w[1][1]) * abs(w[2][1]) * ng[1]
						+ abs(w[0][0]) * abs(w[1][0]) * abs(w[2][1]) * ng[2]
						+ abs(w[0][1]) * abs(w[1][0]) * abs(w[2][1]) * ng[3]
						+ abs(w[0][1]) * abs(w[1][1]) * abs(w[2][0]) * ng[4]
						+ abs(w[0][0]) * abs(w[1][1]) * abs(w[2][0]) * ng[5]
						+ abs(w[0][0]) * abs(w[1][0]) * abs(w[2][0]) * ng[6]
						+ abs(w[0][1]) * abs(w[1][0]) * abs(w[2][0]) * ng[7];
		}
	}
}

// radial
P_float Interpolation(vtkStructuredPoints* vol,P_float pos[3],const P_float M[16],const P_float spacing[3],int interpolationmode) {return Interpolation(vol,pos[0],pos[1],pos[2],M,spacing,interpolationmode);}
P_float Interpolation(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,const P_float M[16],const P_float spacing[3],int interpolationmode)
{
//	interpolationmode=0: trilinear
//	interpolationmode=1: nearest neighbor 
unsigned short* pt= (unsigned short*) vol->GetScalarPointer();

P_float cc[3],pos[3]={(double)posx,(double)posy,(double)posz}; GetCylindricCoordinates(cc,pos,M);

int* dim=vol->GetDimensions();

pos[0]=cc[0]/spacing[0]; pos[1]=cc[1]/spacing[1]; pos[2]=-cc[2]/spacing[2];
if(pos[0]<1 || pos[0]>dim[0]-2 || pos[1]>dim[1]/2-2) return -1;

if(pos[2]<0) {pos[2]+=dim[2]-1; pos[1]=(P_float)(dim[1]-1)/2.-pos[1];}
else {pos[1]=(P_float)(dim[1]-1)/2.+pos[1];}

P_float w[3][2];
w[0][0]=(floor(pos[0])-pos[0]); w[1][0]=(floor(pos[1])-pos[1]); w[2][0]=(floor(pos[2])-pos[2]);
if(w[0][0]<-0.5) {w[0][1]=w[0][0]; w[0][0]=1+w[0][1];} else w[0][1]=1+w[0][0];
if(w[1][0]<-0.5) {w[1][1]=w[1][0]; w[1][0]=1+w[1][1];} else w[1][1]=1+w[1][0];
if(w[2][0]<-0.5) {w[2][1]=w[2][0]; w[2][0]=1+w[2][1];} else w[2][1]=1+w[2][0];

unsigned short ng[8];
ng[0]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);

if(interpolationmode==1) return ng[0];

ng[1]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[2]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[3]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][0])*dim[0]*dim[1]);
ng[4]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[5]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][0])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[6]=*(pt+(int)(pos[0]+w[0][1])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);
ng[7]=*(pt+(int)(pos[0]+w[0][0])+(int)(pos[1]+w[1][1])*dim[0]+(int)(pos[2]+w[2][1])*dim[0]*dim[1]);

w[0][0]=abs(w[0][0]);w[0][1]=abs(w[0][1]);
w[1][0]=abs(w[1][0]);w[1][1]=abs(w[1][1]);
w[2][0]=abs(w[2][0]);w[2][1]=abs(w[2][1]);

P_float result=w[0][1]*w[1][1]*w[2][1]*ng[0]
			+w[0][0]*w[1][1]*w[2][1]*ng[1]
			+w[0][0]*w[1][0]*w[2][1]*ng[2]
			+w[0][1]*w[1][0]*w[2][1]*ng[3]
			+w[0][1]*w[1][1]*w[2][0]*ng[4]
			+w[0][0]*w[1][1]*w[2][0]*ng[5]
			+w[0][0]*w[1][0]*w[2][0]*ng[6]
			+w[0][1]*w[1][0]*w[2][0]*ng[7];

return result;

}


void GetCylindricCoordinates(P_float cc[3],const P_float p[3],const P_float M[16])
{
// axis along x direction
P_float op[3]={p[0]-M[3],p[1]-M[7],p[2]-M[11]};
cc[0]=op[0]*M[0]+op[1]*M[4]+op[2]*M[8];
P_float projp[3]={op[0]-cc[0]*M[0],op[1]-cc[0]*M[4],op[2]-cc[0]*M[8]};
cc[1]=norm(projp);
P_float projpv=projp[0]*M[1]+projp[1]*M[5]+projp[2]*M[9];
P_float projpw=projp[0]*M[2]+projp[1]*M[6]+projp[2]*M[10];
cc[2]=atan2(projpw,projpv);
}

unsigned short InterpolationSmooth(vtkStructuredPoints* vol,P_float p[3],P_float t1[3],P_float t2[3],P_float **GaussianWeights,int GaussianR,int interpolationmode)
{
unsigned short val=GaussianWeights[0][0]*Interpolation(vol,p[0],p[1],p[2],interpolationmode);
if(val==-1) return -1;
for(int x=1;x<GaussianR;x++)
	{
	val+=GaussianWeights[0][x]*Interpolation(vol,p[0]+x*t1[0],p[1]+x*t1[1],p[2]+x*t1[2],interpolationmode); 
	val+=GaussianWeights[0][x]*Interpolation(vol,p[0]-x*t1[0],p[1]-x*t1[1],p[2]-x*t1[2],interpolationmode); 
	}
for(int y=1;y<GaussianR;y++)
	{
	val+=GaussianWeights[y][0]*Interpolation(vol,p[0]+y*t2[0],p[1]+y*t2[1],p[2]+y*t2[2],interpolationmode); 
	val+=GaussianWeights[y][0]*Interpolation(vol,p[0]-y*t2[0],p[1]-y*t2[1],p[2]-y*t2[2],interpolationmode); 
	}
for(int y=1;y<GaussianR;y++)
	for(int x=1;x<GaussianR;x++)
		{
		val+=GaussianWeights[y][x]*Interpolation(vol,p[0]+x*t1[0]+y*t2[0],p[1]+x*t1[1]+y*t2[1],p[2]+x*t1[2]+y*t2[2],interpolationmode); 
		val+=GaussianWeights[y][x]*Interpolation(vol,p[0]-x*t1[0]+y*t2[0],p[1]-x*t1[1]+y*t2[1],p[2]-x*t1[2]+y*t2[2],interpolationmode); 
		val+=GaussianWeights[y][x]*Interpolation(vol,p[0]+x*t1[0]-y*t2[0],p[1]+x*t1[1]-y*t2[1],p[2]+x*t1[2]-y*t2[2],interpolationmode); 
		val+=GaussianWeights[y][x]*Interpolation(vol,p[0]-x*t1[0]-y*t2[0],p[1]-x*t1[1]-y*t2[1],p[2]-x*t1[2]-y*t2[2],interpolationmode); 
		}
return val;
}


void Normalize(vtkStructuredPoints* vol,const unsigned short Max,const unsigned short Min)
{
unsigned short *ptr = (unsigned short*) vol->GetScalarPointer(),max=GetMax(vol),min=GetMin(vol);
for(int i=0;i<vol->GetNumberOfPoints();i++) *(ptr+i)=(unsigned short)(Min+ ((Max-Min)*(*(ptr+i)-min)) / (max-min) );
}

unsigned short GetMax(vtkStructuredPoints* vol)
{
unsigned short max=0;
unsigned short* ptr = (unsigned short*) vol->GetScalarPointer();
for(int i=0;i<vol->GetNumberOfPoints();i++) 	if(*(ptr+i)>max) max=*(ptr+i); 
return max;
}

unsigned short GetMin(vtkStructuredPoints* vol)
{
unsigned short min=30000;
unsigned short* ptr = (unsigned short*) vol->GetScalarPointer();
for(int i=0;i<vol->GetNumberOfPoints();i++) 	if(*(ptr+i)<min) min=*(ptr+i); 
return min;
}


vtkPolyData* Concatenate(vtkPolyData* model1,vtkPolyData* model2)
{
vtkPolyData* ret=vtkPolyData::New();
int nb1=model1->GetNumberOfPoints();
int nb2=model2->GetNumberOfPoints();
int i;

model1->Update();
model2->Update();

vtkPoints* pts=vtkPoints::New();
pts->SetNumberOfPoints(nb1+nb2);
ret->SetPoints(pts);

ret->Allocate(model1->GetNumberOfCells()+model2->GetNumberOfCells(),model1->GetNumberOfCells()+model2->GetNumberOfCells());

for(i=0;i<nb2;i++) pts->SetPoint(i,model2->GetPoint(i));
for(i=0;i<nb1;i++) pts->SetPoint(i+nb2,model1->GetPoint(i));

vtkTriangle* triangle=vtkTriangle::New();

for(i=0;i<model2->GetNumberOfCells();i++) 
	{ 
	triangle->GetPointIds()->SetId(0,model2->GetCell(i)->GetPointId(0)); 	
	triangle->GetPointIds()->SetId(1,model2->GetCell(i)->GetPointId(1)); 	
	triangle->GetPointIds()->SetId(2,model2->GetCell(i)->GetPointId(2));
	ret->InsertNextCell(triangle->GetCellType(),triangle->GetPointIds()); // Insert triangle
	}
for(i=0;i<model1->GetNumberOfCells();i++) 
	{ 
	triangle->GetPointIds()->SetId(0,model1->GetCell(i)->GetPointId(0)+nb2); 	
	triangle->GetPointIds()->SetId(1,model1->GetCell(i)->GetPointId(1)+nb2); 	
	triangle->GetPointIds()->SetId(2,model1->GetCell(i)->GetPointId(2)+nb2);
	ret->InsertNextCell(triangle->GetCellType(),triangle->GetPointIds()); // Insert triangle
	}

ret->Update();
triangle->Delete();

vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
	normals->SetInput(ret);
	normals->SplittingOff();
	normals->ConsistencyOn();
	normals->ComputePointNormalsOn();
	normals->Update();

ret->DeepCopy(normals->GetOutput());
normals->Delete();

return ret;
}


void Pad(vtkStructuredPoints* vol,unsigned short val, P_float dist)
{
int offset[3]={dist/vol->GetSpacing()[0],dist/vol->GetSpacing()[1],dist/vol->GetSpacing()[2]};
vtkImageConstantPad* pad=vtkImageConstantPad::New();
	pad->SetInput(vol);
	pad->SetOutputWholeExtent(-offset[0],vol->GetDimensions()[0]-1+offset[0],-offset[1],vol->GetDimensions()[1]-1+offset[1],-offset[2],vol->GetDimensions()[2]-1+offset[2]);
	pad->SetConstant((double)val);
	pad->Update();
	vol->DeepCopy(pad->GetOutput());
	vol->UpdateData();
pad->Delete();
}


void Differences(P_float diff[3],P_float M1[16],P_float M2[16])
{
// returns the difference (error) betwen the two transforms
P_float M[16],M1_inv[16];
Invert_M(M1,M1_inv);
Multi_M(M,M2,M1_inv);
P_float rt[3],tr[3],q[4]; 
MtoQ(tr,q,M,false); diff[0]=norm(tr); diff[1]=2.*180.*(P_float)acos(q[0])/PI; 
MtoN(tr,rt,M); diff[2]=norm(rt);
}


void Differences(P_float distdev[2],vtkPolyData* model1,vtkPolyData* model2,bool pointtopoint)
{
int i;
P_float c=0,d=0,std_dev=0,p[3],p2[3];

if(pointtopoint)
	{
	for(i=0;i<model1->GetNumberOfPoints();i++)
		{
		model1->GetPoint(i,p);	model2->GetPoint(i,p2); 
		d+=dist3D(p,p2); 
		}
	d=d/model1->GetNumberOfPoints();
	for(i=0;i<model1->GetNumberOfPoints();i++)
		{
		model1->GetPoint(i,p);	model2->GetPoint(i,p2); 
		std_dev+=(dist3D(p,p2)-d)*(dist3D(p,p2)-d); 
		}
	std_dev=sqrt(std_dev/(model1->GetNumberOfPoints()-1));
	}
else
	{
/// TO COMPLETE
	}

distdev[0]=d; distdev[1]=std_dev;
}


P_float GaussianFilter(P_float std_dev2, int nb_val, P_float* val, P_float* dist2)
{
P_float v=0,w,s=0;
for(int i=0;i<nb_val;i++)
	{
	w=exp(-dist2[i]/std_dev2);
	v+=w*val[i];
	s+=w;
	}
v=v/s;
return v;
}

void LdmRigidTransform(const char* filein,const char* fileout,P_float M[16])
{
int i;
// load points
vtkPoints* ldm_in=vtkPoints::New();
vtkPoints* ldm_out=vtkPoints::New();
int nb_p; P_float p[3]; double dp[3];
FILE *f,*f2; 		f=fopen(filein,"rt");		fscanf(f,"%d \n",&nb_p);		f2=fopen(fileout,"rt");		fscanf(f2,"%d \n",&nb_p);

ldm_in->SetDataTypeToFloat(); ldm_in->SetNumberOfPoints(nb_p);		ldm_out->SetDataTypeToFloat(); ldm_out->SetNumberOfPoints(nb_p);
for(i=0;i<nb_p;i++) 
	{ 
	fscanf(f,"%f %f %f \n",&p[0],&p[1],&p[2]); dp[0]=(double)p[0];dp[1]=(double)p[1];dp[2]=(double)p[2];			ldm_in->SetPoint(i,dp);
	fscanf(f2,"%f %f %f \n",&p[0],&p[1],&p[2]); dp[0]=(double)p[0];dp[1]=(double)p[1];dp[2]=(double)p[2];			ldm_out->SetPoint(i,dp);
	}
fclose(f); fclose(f2);
LdmRigidTransform(ldm_in,ldm_out,M);
// clean
ldm_in->Delete();
ldm_out->Delete();
}

void LdmRigidTransform(const char* filein,vtkPoints* ldm_out,P_float M[16])
{
int i;
// load points
vtkPoints* ldm_in=vtkPoints::New();
int nb_p; P_float p[3]; double dp[3];
FILE *f; 		f=fopen(filein,"rt");		fscanf(f,"%d \n",&nb_p);		

ldm_in->SetDataTypeToFloat(); ldm_in->SetNumberOfPoints(nb_p);		
for(i=0;i<nb_p;i++) 
	{ 
	fscanf(f,"%lf %lf %lf \n",&p[0],&p[1],&p[2]); dp[0]=(double)p[0];dp[1]=(double)p[1];dp[2]=(double)p[2];			ldm_in->SetPoint(i,dp);
	}
fclose(f);
LdmRigidTransform(ldm_in,ldm_out,M);
// clean
ldm_in->Delete();
}

void LdmRigidTransform(vtkPoints* ldm_in,vtkPoints* ldm_out,P_float M[16])
{
/*
int nb=ldm_in->GetNumberOfPoints();
P_float* pts=new P_float[3*nb];
P_float* f=new P_float[3*nb];
bool* ign=new bool[nb];
for(int i=0;i<nb;i++) {ldm_in->GetPoint(i,pts+3*i); ldm_out->GetPoint(i,f+3*i); f[3*i]-=pts[3*i]; f[3*i+1]-=pts[3*i+1]; f[3*i+2]-=pts[3*i+2]; ign[i]=false;}
Regularize_Rigid(M,nb,pts,f,ign);
free(pts); free(ign); free(f); return;
*/
vtkLandmarkTransform* T=vtkLandmarkTransform::New();
T->SetSourceLandmarks(ldm_in);
T->SetTargetLandmarks(ldm_out);
T->SetModeToRigidBody ();
T->Update();
for(int i=0;i<4;i++) for(int j=0;j<4;j++) M[j+i*4]=T->GetMatrix ()->GetElement(i,j);
}

void FlipCSystem_X(P_float M[16])
{
M[0]=-M[0]; M[1]=-M[1]; M[3]=-M[3];
P_float x[3]={M[0],M[4],M[8]};
P_float y[3]={M[1],M[5],M[9]};
P_float z[3]; crossproduct(z,x,y); M[2]=z[0]; M[6]=z[1]; M[10]=z[2];
}

void FlipModel_X(vtkPolyData* model)
{
if(model->GetNumberOfPoints()==0) return;
for(int i=0;i<model->GetNumberOfPoints();i++) model->GetPoints()->SetPoint(i,-model->GetPoint(i)[0],model->GetPoint(i)[1],model->GetPoint(i)[2]);
vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
	normals->SetInput(model);
	normals->SplittingOff();
	normals->ConsistencyOn();
	normals->ComputePointNormalsOn();
	normals->FlipNormalsOn();
	normals->Update();
model->DeepCopy(normals->GetOutput());
normals->Delete();
}

void FlipVolume_X(vtkStructuredPoints* vol)
{
vtkStructuredPoints* vol2=vtkStructuredPoints::New();
vol2->DeepCopy(vol);
vol2->SetOrigin(-vol->GetSpacing()[0]*(double)vol->GetDimensions()[0]-vol->GetOrigin()[0],vol->GetOrigin()[1],vol->GetOrigin()[2]);
int* p_vol= (int*) vol->GetScalarPointer();
int* p_vol2= (int*) vol2->GetScalarPointer();
for(int z=0;z<vol->GetDimensions()[2];z++)
for(int y=0;y<vol->GetDimensions()[1];y++)
for(int x=0;x<vol->GetDimensions()[0];x++)
*(p_vol2+x+y*vol->GetDimensions()[0]+z*vol->GetDimensions()[0]*vol->GetDimensions()[1])=*(p_vol+vol->GetDimensions()[0]-x-1+y*vol->GetDimensions()[0]+z*vol->GetDimensions()[0]*vol->GetDimensions()[1]);
vol->DeepCopy(vol2);
}

P_float GetEquidistant2(P_float F[4],P_float Fp[4],P_float P[3],P_float P1[3],P_float P2[3],P_float d12,P_float d22)
{
P_float P2P1[3]={P1[0]-P2[0],P1[1]-P2[1],P1[2]-P2[2]};
P_float P2P12=P2P1[0]*P2P1[0]+P2P1[1]*P2P1[1]+P2P1[2]*P2P1[2];
P_float f=(d22-d12)/(2*P2P12);	F[0]=P2P1[0]*f; F[1]=P2P1[1]*f; F[2]=P2P1[2]*f;	F[3]=norm(F);
Fp[0]=(P[0]-(P1[0]+P2[0])/2)+F[0]; Fp[1]=(P[1]-(P1[1]+P2[1])/2)+F[1]; Fp[2]=(P[2]-(P1[2]+P2[2])/2)+F[2]; Fp[3]=norm(Fp);
P_float R[3]; R[0]=P[0]-P1[0]+F[0]; R[1]=P[1]-P1[1]+F[1]; R[2]=P[2]-P1[2]+F[2]; return norm(R);
}

P_float GetEquidistant3(P_float F[4],P_float P[3],P_float P1[3],P_float P2[3],P_float P3[3])
{
P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]};	P_float P1pP2[3]={P2[0]+P1[0],P2[1]+P1[1],P2[2]+P1[2]};	P_float DP1=dotproduct(P1P2,P1pP2);
P_float P2P3[3]={P3[0]-P2[0],P3[1]-P2[1],P3[2]-P2[2]};	P_float P2pP3[3]={P3[0]+P2[0],P3[1]+P2[1],P3[2]+P2[2]};	P_float DP2=dotproduct(P2P3,P2pP3);
P_float u[3]; crossproduct(u,P1P2,P2P3);	P_float u2=u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
P_float v[3]; v[0]=DP2*P1P2[0]-DP1*P2P3[0];	v[1]=DP2*P1P2[1]-DP1*P2P3[1];	v[2]=DP2*P1P2[2]-DP1*P2P3[2];
P_float QP[3]; crossproduct(QP,u,v); QP[0]=P[0]-QP[0]/(2*u2); QP[1]=P[1]-QP[1]/(2*u2); QP[2]=P[2]-QP[2]/(2*u2);
P_float dp=dotproduct(QP,u); dp=dp/u2;
F[0]=dp*u[0]-QP[0];	F[1]=dp*u[1]-QP[1];	F[2]=dp*u[2]-QP[2];	F[3]=norm(F);
P_float R[3]; R[0]=P[0]-P1[0]+F[0]; R[1]=P[1]-P1[1]+F[1]; R[2]=P[2]-P1[2]+F[2]; return norm(R);
}

P_float Closest(P_float fa,P_float fb,P_float fc,P_float fd,P_float fe,P_float ff,P_float fout[2]) 
{
P_float fdet=abs(fa*fc-fb*fb);
P_float fs=fb*fe-fc*fd;
P_float ft=fb*fd-fa*fe;
P_float dist;

if ( fs + ft <= fdet )
    {
        if ( fs < 0 )
        {
            if ( ft < 0 )  // region 4
            {
                if ( fd < 0 )
                {
                    ft = 0;
                    if ( -fd >= fa )
                    {
                        fs = 1;
                        dist = fa+2*fd+ff;
                    }
                    else
                    {
                        fs = -fd/fa;
                        dist = fd*fs+ff;
                    }
                }
                else
                {
                    fs = 0;
                    if ( fe >= 0 )
                    {
                        ft = 0;
                        dist = ff;
                    }
                    else if ( -fe >= fc )
                    {
                        ft = 1;
                        dist = fc+2*fe+ff;
                    }
                    else
                    {
                        ft = -fe/fc;
                        dist = fe*ft+ff;
                    }
                }
            }
            else  // region 3
            {
                fs = 0;
                if ( fe >= 0 )
                {
                    ft = 0;
                    dist = ff;
                }
                else if ( -fe >= fc )
                {
                    ft = 1;
                    dist = fc+2*fe+ff;
                }
                else
                {
                    ft = -fe/fc;
                    dist = fe*ft+ff;
                }
            }
        }
        else if ( ft < 0 )  // region 5
        {
            ft = 0;
            if ( fd >= 0 )
            {
                fs = 0;
                dist = ff;
            }
            else if ( -fd >= fa )
            {
                fs = 1;
                dist = fa+2*fd+ff;
            }
            else
            {
                fs = -fd/fa;
                dist = fd*fs+ff;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            P_float fInvDet = 1/fdet;
            fs *= fInvDet;
            ft *= fInvDet;
            dist = fs*(fa*fs+fb*ft+2*fd) +  ft*(fb*fs+fc*ft+2*fe)+ff;
        }
    }
    else
    {
        P_float ftmp0, ftmp1, fNumer, fDenom;

        if ( fs < 0 )  // region 2
        {
            ftmp0 = fb + fd;
            ftmp1 = fc + fe;
            if ( ftmp1 > ftmp0 )
            {
                fNumer = ftmp1 - ftmp0;
                fDenom = fa-2*fb+fc;
                if ( fNumer >= fDenom )
                {
                    fs = 1;
                    ft = 0;
                    dist = fa+2*fd+ff;
                }
                else
                {
                    fs = fNumer/fDenom;
                    ft = 1 - fs;
                    dist = fs*(fa*fs+fb*ft+2*fd) +   ft*(fb*fs+fc*ft+2*fe)+ff;
                }
            }
            else
            {
                fs = 0;
                if ( ftmp1 <= 0 )
                {
                    ft = 1;
                    dist = fc+2*fe+ff;
                }
                else if ( fe >= 0 )
                {
                    ft = 0;
                    dist = ff;
                }
                else
                {
                    ft = -fe/fc;
                    dist = fe*ft+ff;
                }
            }
        }
        else if ( ft < 0 )  // region 6
        {
            ftmp0 = fb + fe;
            ftmp1 = fa + fd;
            if ( ftmp1 > ftmp0 )
            {
                fNumer = ftmp1 - ftmp0;
                fDenom = fa-2*fb+fc;
                if ( fNumer >= fDenom )
                {
                    ft = 1;
                    fs = 0;
                    dist = fc+2*fe+ff;
                }
                else
                {
                    ft = fNumer/fDenom;
                    fs = 1 - ft;
                    dist = fs*(fa*fs+fb*ft+2*fd) +   ft*(fb*fs+fc*ft+2*fe)+ff;
                }
            }
            else
            {
                ft = 0;
                if ( ftmp1 <= 0 )
                {
                    fs = 1;
                    dist = fa+2*fd+ff;
                }
                else if ( fd >= 0 )
                {
                    fs = 0;
                    dist = ff;
                }
                else
                {
                    fs = -fd/fa;
                    dist = fd*fs+ff;
                }
            }
        }
        else  // region 1
        {
            fNumer = fc + fe - fb - fd;
            if ( fNumer <= 0 )
            {
                fs = 0;
                ft = 1;
                dist = fc+2*fe+ff;
            }
            else
            {
                fDenom = fa-2*fb+fc;
                if ( fNumer >= fDenom )
                {
                    fs = 1;
                    ft = 0;
                    dist = fa+2*fd+ff;
                }
                else
                {
                    fs = fNumer/fDenom;
                    ft = 1 - fs;
                    dist = fs*(fa*fs+fb*ft+2*fd) +         ft*(fb*fs+fc*ft+2*fe)+ff;
                }
            }
        }
    }

fout[0]=fs; fout[1]=ft;
return dist;
}


void Regularize_Rigid(P_float T[16],int nb_points,P_float *points,P_float* f,bool* ign)
{
int i,j,k;
P_float tr[3]={0,0,0},e=0;

P_float c[3],C[3],cm[3]={0,0,0},Cm[3]={0,0,0},M[16],M2[16],w[4],**v=new P_float*[4],**A=new P_float*[4];
for(k=0;k<4;k++) {v[k]=new P_float[4]; A[k]=new P_float[4]; for(j=0;j<4;j++) A[k][j]=0;}
//barycenter
int count=0;
for(i=0;i<nb_points;i++)
	if(!ign[i])
		{
		cm[0]+=points[3*i]; cm[1]+=points[3*i+1]; cm[2]+=points[3*i+2]; 
		Cm[0]+=points[3*i]+f[3*i]; Cm[1]+=points[3*i+1]+f[3*i+1]; Cm[2]+=points[3*i+2]+f[3*i+2]; 
		count++;
		}
cm[0]=cm[0]/count; cm[1]=cm[1]/count; cm[2]=cm[2]/count;
Cm[0]=Cm[0]/count; Cm[1]=Cm[1]/count; Cm[2]=Cm[2]/count;
//rotation
for(i=0;i<nb_points;i++)
	if(!ign[i])
		{
		c[0]=points[3*i]-cm[0]; c[1]=points[3*i+1]-cm[1];	c[2]=points[3*i+2]-cm[2];
		C[0]=points[3*i]+f[3*i]-Cm[0]; C[1]=points[3*i+1]+f[3*i+1]-Cm[1];	C[2]=points[3*i+2]+f[3*i+2]-Cm[2];
		M[0]=0;				M[1]=c[0]-C[0];		M[2]=c[1]-C[1];		M[3]=c[2]-C[2];
		M[4]=C[0]-c[0];		M[5]=0;				M[6]=-c[2]-C[2];	M[7]=c[1]+C[1];
		M[8]=C[1]-c[1];		M[9]=c[2]+C[2];		M[10]=0;			M[11]=-c[0]-C[0];
		M[12]=C[2]-c[2];	M[13]=-c[1]-C[1];	M[14]=c[0]+C[0];	M[15]=0;
		Multi_M(M2,M,M);
		for(k=0;k<4;k++) for(j=0;j<4;j++) A[k][j]-=M2[j+k*4]; 
		}
vtkMath::JacobiN(A,4,w,v);
P_float V[4]={v[0][3],v[1][3],v[2][3],v[3][3]};
QtoM(tr,V,T,0);

//translation
Transform(cm,c,T); T[3]=Cm[0]-c[0]; T[7]=Cm[1]-c[1]; T[11]=Cm[2]-c[2]; 
for(k=0;k<4;k++) {free(v[k]); free(A[k]);} free(v); free(A);
}

void Regularize_Simi(P_float T[16],int nb_points,P_float *points,P_float* f,bool* ign)
{
int i,j,k;
P_float tr[3]={0,0,0},e=0;
P_float sx=0,sxy[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},c[3],C[3],cm[3]={0,0,0},Cm[3]={0,0,0},M[16],M2[16],w[4],**v=new P_float*[4],**A=new P_float*[4];
for(k=0;k<4;k++) {v[k]=new P_float[4]; A[k]=new P_float[4]; for(j=0;j<4;j++) A[k][j]=0;}
//barycenter
int count=0;
for(i=0;i<nb_points;i++)
	if(!ign[i])
		{
		cm[0]+=points[3*i]; cm[1]+=points[3*i+1]; cm[2]+=points[3*i+2]; 
		Cm[0]+=points[3*i]+f[3*i]; Cm[1]+=points[3*i+1]+f[3*i+1]; Cm[2]+=points[3*i+2]+f[3*i+2]; 
		count++;
		}
cm[0]=cm[0]/count; cm[1]=cm[1]/count; cm[2]=cm[2]/count;
Cm[0]=Cm[0]/count; Cm[1]=Cm[1]/count; Cm[2]=Cm[2]/count;
//rotation
for(i=0;i<nb_points;i++)
	if(!ign[i])
		{
		c[0]=points[3*i]-cm[0]; c[1]=points[3*i+1]-cm[1];	c[2]=points[3*i+2]-cm[2];
		C[0]=points[3*i]+f[3*i]-Cm[0]; C[1]=points[3*i+1]+f[3*i+1]-Cm[1];	C[2]=points[3*i+2]+f[3*i+2]-Cm[2];
		sx+=dotproduct(c,c); 
		sxy[0]+=C[0]*c[0];    sxy[1]+=C[1]*c[0];  sxy[2]+=C[2]*c[0]; 
		sxy[4]+=C[0]*c[1];    sxy[5]+=C[1]*c[1];  sxy[6]+=C[2]*c[1]; 
		sxy[8]+=C[0]*c[2];    sxy[9]+=C[1]*c[2];  sxy[10]+=C[2]*c[2]; 
		M[0]=0;				M[1]=c[0]-C[0];		M[2]=c[1]-C[1];		M[3]=c[2]-C[2];
		M[4]=C[0]-c[0];		M[5]=0;				M[6]=-c[2]-C[2];	M[7]=c[1]+C[1];
		M[8]=C[1]-c[1];		M[9]=c[2]+C[2];		M[10]=0;			M[11]=-c[0]-C[0];
		M[12]=C[2]-c[2];	M[13]=-c[1]-C[1];	M[14]=c[0]+C[0];	M[15]=0;
		Multi_M(M2,M,M);
		for(k=0;k<4;k++) for(j=0;j<4;j++) A[k][j]-=M2[j+k*4]; 
		}
vtkMath::JacobiN(A,4,w,v);
P_float V[4]={v[0][3],v[1][3],v[2][3],v[3][3]};
QtoM(tr,V,T,0);
// homothetie
Multi_M(M,T,sxy); e=(M[0]+M[5]+M[10])/sx;
T[0]=e*T[0]; T[1]=e*T[1]; T[2]=e*T[2]; T[4]=e*T[4]; T[5]=e*T[5]; T[6]=e*T[6]; T[8]=e*T[8]; T[9]=e*T[9]; T[10]=e*T[10]; 
//translation
Transform(cm,c,T); T[3]=Cm[0]-c[0]; T[7]=Cm[1]-c[1]; T[11]=Cm[2]-c[2]; 
for(k=0;k<4;k++) {free(v[k]); free(A[k]);} free(v); free(A);
}





void Regularize_Affine(P_float T[16], int nb_points, P_float *points, P_float* f, bool* ign) {
	int i;
	P_float tr[3] = {0, 0, 0}, e = 0;
	P_float sxx[6] = {0, 0, 0, 0, 0, 0}, sxx_i[6], syx[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, c[3], C[3], cm[3] = {0, 0, 0}, Cm[3] = {0, 0, 0};
	
	//barycenter
	int count = 0;
	for (i = 0; i < nb_points; i++) {
		if (!ign[i]) {
			cm[0] += points[3 * i]; 
			cm[1] += points[3 * i + 1]; 
			cm[2] += points[3 * i + 2]; 
			
			Cm[0] += points[3 * i] + f[3 * i]; 
			Cm[1] += points[3 * i + 1] + f[3 * i + 1]; 
			Cm[2] += points[3 * i + 2] + f[3 * i + 2]; 
			
			count++;
		}
	}

	cm[0] = cm[0] / count; 
	cm[1] = cm[1] / count; 
	cm[2] = cm[2] / count;

	Cm[0] = Cm[0] / count; 
	Cm[1] = Cm[1] / count; 
	Cm[2] = Cm[2] / count;

	//affine
	for (i = 0; i < nb_points; i++) {
		if (!ign[i]) {
			c[0] = points[3 * i] - cm[0]; 
			c[1] = points[3 * i + 1] - cm[1];	
			c[2] = points[3 * i + 2] - cm[2];
			
			C[0] = points[3 * i] + f[3 * i] - Cm[0]; 
			C[1] = points[3 * i + 1] + f[3 * i + 1] - Cm[1];
			C[2] = points[3 * i + 2] + f[3 * i + 2] - Cm[2];

			syx[0] += C[0] * c[0];		
			syx[1] += C[0] * c[1];  
			syx[2] += C[0] * c[2];
			syx[3] += C[1] * c[0];
			syx[4] += C[1] * c[1];
			syx[5] += C[1] * c[2];
			syx[6] += C[2] * c[0];		
			syx[7] += C[2] * c[1]; 
			syx[8] += C[2] * c[2]; 

			sxx[0] += c[0] * c[0];   
			sxx[1] += c[1] * c[0];
			sxx[2] += c[1] * c[1];  
			sxx[3] += c[2] * c[0];
			sxx[4] += c[2] * c[1];  
			sxx[5] += c[2] * c[2]; 
		}
	}
	
	Invert_MSymmetric(sxx, sxx_i);
	T[0] = syx[0] * sxx_i[0] + syx[1] * sxx_i[1] + syx[2] * sxx_i[3]; 
	T[1] = syx[0] * sxx_i[1] + syx[1] * sxx_i[2] + syx[2] * sxx_i[4]; 
	T[2] = syx[0] * sxx_i[3] + syx[1] * sxx_i[4] + syx[2] * sxx_i[5];
	T[3] = 0;
	T[4] = syx[3] * sxx_i[0] + syx[4] * sxx_i[1] + syx[5] * sxx_i[3];
	T[5] = syx[3] * sxx_i[1] + syx[4] * sxx_i[2] + syx[5] * sxx_i[4];
	T[6] = syx[3] * sxx_i[3] + syx[4] * sxx_i[4] + syx[5] * sxx_i[5];
	T[7] = 0;
	T[8] = syx[6] * sxx_i[0] + syx[7] * sxx_i[1] + syx[8] * sxx_i[3]; 
	T[9] = syx[6] * sxx_i[1] + syx[7] * sxx_i[2] + syx[8] * sxx_i[4];
	T[10] = syx[6] * sxx_i[3] + syx[7] * sxx_i[4] + syx[8] * sxx_i[5];
	T[11] = 0;
	T[12] = 0; 
	T[13] = 0; 
	T[14] = 0; 
	T[15] = 0; 
	
	//translation
	Transform(cm, c, T); 
	T[3] = Cm[0] - c[0]; 
	T[7] = Cm[1] - c[1]; 
	T[11] = Cm[2] - c[2]; 
}




void Regularize_Centered(P_float T[16],int nb_points,P_float *points,P_float* f,P_float RC[3],bool TR,bool* ign)
{
int i,j,k;
P_float tr[3]={0,0,0},e=0;

P_float c[3],C[3],cm[3]={0,0,0},Cm[3]={0,0,0},M[16],M2[16],w[4],**v=new P_float*[4],**A=new P_float*[4];
for(k=0;k<4;k++) {v[k]=new P_float[4]; A[k]=new P_float[4]; for(j=0;j<4;j++) A[k][j]=0;}
//barycenter
if(TR)
	{
	int count=0;
	for(i=0;i<nb_points;i++)
		if(!ign[i])
			{
			cm[0]+=points[3*i]-RC[0]; cm[1]+=points[3*i+1]-RC[1]; cm[2]+=points[3*i+2]-RC[2]; 
			Cm[0]+=points[3*i]+f[3*i]-RC[0]; Cm[1]+=points[3*i+1]+f[3*i+1]-RC[1]; Cm[2]+=points[3*i+2]+f[3*i+2]-RC[2]; 
			count++;
			}
	cm[0]=cm[0]/count; cm[1]=cm[1]/count; cm[2]=cm[2]/count;
	Cm[0]=Cm[0]/count; Cm[1]=Cm[1]/count; Cm[2]=Cm[2]/count;
	}
//rotation
for(i=0;i<nb_points;i++)
	if(!ign[i])
		{
		c[0]=points[3*i]-cm[0]-RC[0]; c[1]=points[3*i+1]-cm[1]-RC[1];	c[2]=points[3*i+2]-cm[2]-RC[2];
		C[0]=points[3*i]+f[3*i]-Cm[0]-RC[0]; C[1]=points[3*i+1]+f[3*i+1]-Cm[1]-RC[1];	C[2]=points[3*i+2]+f[3*i+2]-Cm[2]-RC[2];
		M[0]=0;				M[1]=c[0]-C[0];		M[2]=c[1]-C[1];		M[3]=c[2]-C[2];
		M[4]=C[0]-c[0];		M[5]=0;				M[6]=-c[2]-C[2];	M[7]=c[1]+C[1];
		M[8]=C[1]-c[1];		M[9]=c[2]+C[2];		M[10]=0;			M[11]=-c[0]-C[0];
		M[12]=C[2]-c[2];	M[13]=-c[1]-C[1];	M[14]=c[0]+C[0];	M[15]=0;
		Multi_M(M2,M,M);
		for(k=0;k<4;k++) for(j=0;j<4;j++) A[k][j]-=M2[j+k*4]; 
		}
vtkMath::JacobiN(A,4,w,v);
P_float V[4]={v[0][3],v[1][3],v[2][3],v[3][3]};
QtoM(tr,V,T,0);
//translation
if(TR)
	{
	Transform(cm,c,T); 
	tr[0]=Cm[0]-c[0]; tr[1]=Cm[1]-c[1]; tr[2]=Cm[2]-c[2]; 
	}
// rotation recentration
T[3]=RC[0]-T[0]*RC[0]-T[1]*RC[1]-T[2]*RC[2]+tr[0];
T[7]=RC[1]-T[4]*RC[0]-T[5]*RC[1]-T[6]*RC[2]+tr[1];
T[11]=RC[2]-T[8]*RC[0]-T[9]*RC[1]-T[10]*RC[2]+tr[2];
for(k=0;k<4;k++) {free(v[k]); free(A[k]);} free(v); free(A);
}


void Regularize_Translation(P_float T[16],int nbpoints,P_float* f,bool* ign)
{
Identity(T);
int count=0; for(int i=0;i<nbpoints;i++) if(!ign[i]) {T[3]+=f[3*i]; T[7]+=f[3*i+1]; T[11]+=f[3*i+2]; count++;}
if(count!=0) {T[3]=T[3]/count; T[7]=T[7]/count; T[11]=T[11]/count; }
}



void AnalyseEnergy(P_float ACC_DO_NOM_RON_CR[5],P_float *E,int nbpoints,int depth,P_float s,P_float thresh)
{
int i,j,L=2*depth+1,kdo=2;
P_float ACC=0,DO=0,NOM=0,RON=0,CR=0;

P_float* dE=new P_float[nbpoints*L];
for(i=0;i<nbpoints;i++) {dE[i*L]=0;	for(j=1;j<L;j++)  dE[i*L+j]=(E[i*L+j]-E[i*L+j-1]);}

int jmin,jcr;
P_float Emin,E1;
for(i=0;i<nbpoints;i++)
	{
	jcr=1E10; jmin=0; Emin=thresh;
	for(j=0;j<L;j++) if(E[i*L+j]<Emin) {Emin=E[i*L+j]; jmin=j-depth;}
	if(Emin!=thresh) 
		{
		ACC+=(abs(jmin)*s);
		DO-=2*Emin; if(jmin+depth-kdo>=0) E1=E[i*L+jmin+depth-kdo]; else E1=E[i*L]; if(E1!=thresh) DO+=E1; if(jmin+depth+kdo<L) E1=E[i*L+jmin+depth+kdo]; else E1=E[i*L+L-1]; if(E1!=thresh) DO+=E1;
		for(j=1;j<L;j++) if(dE[i*L+j-1]<0 && dE[i*L+j]>0) {NOM++; if(abs(jmin-j+depth)<jcr) jcr=abs(jmin-j+depth);} //local maxima
		for(j=0;j<L;j++) {if(j<jmin+depth && dE[i*L+j]>0) RON+=dE[i*L+j]; else if(j>jmin+depth && dE[i*L+j]<0) RON-=dE[i*L+j];}//false positive
		if(jcr!=(int)1E10) CR+=s*jcr;
		}
	}
free(E); free(dE);

ACC=ACC/(P_float)nbpoints;  DO=DO/((P_float)nbpoints*2*kdo*s); NOM=NOM/(P_float)nbpoints; RON=RON/((P_float)nbpoints*2*L*s); CR=CR/(P_float)nbpoints;

ACC_DO_NOM_RON_CR[0]=ACC;ACC_DO_NOM_RON_CR[1]=DO;ACC_DO_NOM_RON_CR[2]=NOM;ACC_DO_NOM_RON_CR[3]=RON;ACC_DO_NOM_RON_CR[4]=CR;
}






void GetIPForces(P_float* f,bool *ign,P_float *E,int nbpoints,int mn,int mnref,P_float s,P_float *normals,P_float alpha,P_float thresh)
{
int depth=mn-mnref;
GetIPForces(f,ign,E,nbpoints,depth,s,normals,alpha,thresh);
}

void GetIPForces(P_float* f, bool *ign, P_float *E, int nbpoints, int depth, P_float s, P_float *normals, P_float alpha, P_float thresh) {
	int i, j, jmin;
	P_float Emin;
	
	ProcessE(E, nbpoints, 2 * depth + 1);
	
	for (i = 0; i < nbpoints; i++) {
		jmin = 0; 
		Emin = thresh;
		


		

		// Line files were written out in VTK PolyData format in GetGradEnergies()
		// Here, we append those files with color values from E[] array. 
		//std::stringstream ss;
		//ss << "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterial_Deformation_TwoBoxes\\Lines\\Line_Point_" << i << ".vtk";
		//ofstream file(ss.str().c_str(), ios::app);	
		//std::vector<double> Es;

		//cout << "Point " << i << ": ";
		for (j = -depth; j <= depth; j++)   {
			//Es.push_back(E[i * (2 * depth + 1) + j + depth]);
			//cout << E[i * (2 * depth + 1) + j + depth] << ", ";
			
			if (E[i * (2 * depth + 1) + j + depth] < Emin) {
				Emin = E[i * (2 * depth + 1) + j + depth]; 
				jmin = j;
			}
		}
		//cout << endl;
		//if (
		//	i == 112 || i == 113 || i == 136 || i == 137 || i == 244 || i == 245 || i == 174 || i == 157 || 
		//	i == 156 || i == 225 || i == 224 || i == 249 || i == 248 || i == 213 || i == 236 || i == 237
		//	) {
		//		cout << "Point " << i << ": ";
		//		for (j = -depth; j <= depth; j++)   {
		//			cout << E[i * (2 * depth + 1) + j + depth] << "  ";		
		//		}
		//		cout << "  jmin=" << jmin << endl;
		//}		
		//if (			
		//	i == 246 || i == 247 || i == 184 || i == 185 || i == 190 || i == 181 || i == 180 || i == 235 || 
		//	i == 234 || i == 231 || i == 230 || i == 251 || i == 238 || i == 239 || i == 165 || i == 164
		//	) {
		//		cout << "  Point " << i << ": ";
		//		for (j = -depth; j <= depth; j++)   {
		//			cout << E[i * (2 * depth + 1) + j + depth] << "  ";		
		//		}
		//		cout << "  jmin=" << jmin << endl;
		//}
		
		
		
		//file << "\n\nVERTICES 1 2\n1 " << (jmin + depth) << "\n";		
		//file << "\nPOINT_DATA " << (2 * depth + 1) << "\nSCALARS E_Array_Values float 1\nLOOKUP_TABLE default\n";
		//for (int q = 0; q < (2 * depth + 1); q++) {
		//	file << Es.at(q) << "  ";
		//}
		//file.close();



		if (Emin != thresh) ign[i] = false;
	
		f[3 * i] += alpha * (jmin * s * normals[3 * i]);
		f[3 * i + 1] += alpha * (jmin * s * normals[3 * i + 1]); 
		f[3 * i + 2] += alpha * (jmin * s * normals[3 * i + 2]);

		//cout << "\nF = {" << f[3*i] << ", " << f[3*i+1] << ", " << f[3*i+2] << "}" << endl << endl;
	}
	free(E);
}



void GetDemonForces(P_float* f,bool *ign,int nbpoints,P_float* normals,int **IP,int **IPref,P_float step,int mn,int pn,int mnref,int pnref,P_float alpha)
{
int i,j,size=pnref+mnref+1,count,offset=mn-mnref;
P_float l=1./(2.*step*(P_float)(mn-mnref)),d,diff,grad;
P_float l2=l*l;
P_float epsilon=1E-10;

for(i=0;i<nbpoints;i++)
	{
	d=0; count=0;
	for(j=1;j<size-1;j++)  if(IP[i][j+offset]!=-1E5 && IPref[i][j]!=-1E5 && IPref[i][j-1]!=-1E5 && IPref[i][j+1]!=-1E5) {grad=(P_float)IPref[i][j+1]-(P_float)IPref[i][j-1]; diff=(P_float)IPref[i][j]-(P_float)IP[i][j+offset]; d+=(diff*grad)/(grad*grad+l2*diff*diff+epsilon); count++;}
	if(count!=0) {d=d/(P_float)count; ign[i]=false;} 
	f[3*i]+=alpha*(d*normals[3*i]); f[3*i+1]+=alpha*(d*normals[3*i+1]); f[3*i+2]+=alpha*(d*normals[3*i+2]);
	}
}



void ProcessE(P_float* E,int nb_points,int depth)
{/*
int i,j;
P_float* tmp=new P_float[depth];
for(i=0;i<nb_points;i++)
	{
	memcpy(tmp,E+i*depth,depth*sizeof(P_float));
	for(j=1;j<depth-1;j++)	
		if(tmp[j-1]==1E10 || tmp[j]==1E10 || tmp[j+1]==1E10) E[j+i*depth]=1E10;	else E[j+i*depth]=(tmp[i-1]+2*tmp[i]+tmp[i+1])/4.;
	if(tmp[0]==1E10 || tmp[1]==1E10) E[i*depth]=1E10;	else E[i*depth]=(tmp[0]+2*tmp[0]+tmp[1])/4.;
	if(tmp[depth-2]==1E10 || tmp[depth-1]==1E10) E[depth-1+i*depth]=1E10; else E[depth-1+i*depth]=(tmp[depth-2]+2*tmp[depth-1]+tmp[depth-1])/4.;
	}
free(tmp);
*/}


void GetGradEnergies(P_float* E, int nb_points, P_float* points, P_float* normals, vtkStructuredPoints* MRI_grad, P_float step, int depth, int interpolationmode, P_float* ROI, bool oppositedirection) {
	if (MRI_grad == NULL) return;

	int i, j, index;
	P_float p[3], n[3];
	P_float vals[3], val;
	if (MRI_grad->GetNumberOfScalarComponents() == 1) {// NG metric
		for (i = 0; i < nb_points; i++) {
			n[0] = normals[3 * i] * step; 
			n[1] = normals[3 * i + 1] * step;  
			n[2] = normals[3 * i + 2] * step; 
			
			p[0] = points[3 * i] - (P_float)depth * n[0];
			p[1] = points[3 * i + 1] - (P_float)depth * n[1];		
			p[2] = points[3 * i + 2] - (P_float)depth * n[2];
			
			index = i * (2 * depth + 1);
			if (p[0] > ROI[0] && p[0] < ROI[1] && p[1] > ROI[2] && p[1] < ROI[3] && p[2] > ROI[4] && p[2] < ROI[5]) {
				for (j = 0; j < 2 * depth + 1; j++) {
					val = Interpolation(MRI_grad, p, interpolationmode);
					if (val != -1) 
						if (E[index] == 1E10) 
							E[index] = -(P_float)val;

					p[0] += n[0];	
					p[1] += n[1];
					p[2] += n[2];

					index++;
				}
			}
		}
	}
	else { // G metric
	
		for (i = 0; i < nb_points; i++) {
			n[0] = normals[3 * i] * step; 
			n[1] = normals[3 * i + 1] * step;  
			n[2] = normals[3 * i + 2] * step; 

			p[0] = points[3 * i] - (P_float)depth * n[0];		
			p[1] = points[3 * i + 1] - (P_float)depth * n[1];		
			p[2] = points[3 * i + 2] - (P_float)depth * n[2];

			index = i * (2 * depth + 1);

			//std::vector<double*> line;
			
			if (p[0] > ROI[0] && p[0] < ROI[1] && p[1] > ROI[2] && p[1] < ROI[3] && p[2] > ROI[4] && p[2] < ROI[5]) {
				for (j = 0; j < 2 * depth + 1; j++)	{
					vInterpolation(vals, MRI_grad, p, interpolationmode);
					
					if (vals[0] != 1E5) {
						if (E[index] == 1E10) 
							if (oppositedirection)
								E[index] = (dotproduct(vals, normals + 3 * i)); 
							else 
								E[index] = -(dotproduct(vals, normals + 3 * i)); 
					}
					
					//double *temp2 = new double[3]; 
					//temp2[0] = p[0]; temp2[1] = p[1]; temp2[2] = p[2];
					//line.push_back(temp2);

					p[0] += n[0];	
					p[1] += n[1];	
					p[2] += n[2];

					index++;
				}

				//std::stringstream ss;
				//ss << "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterial_Deformation_TwoBoxes\\Lines\\Line_Point_" << i << ".vtk";

				//ofstream f;
				//f.open(ss.str().c_str());

				//f << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS " << line.size() << " float\n";
				//for (int q = 0; q < line.size(); q++) {
				//	double *qp = line.at(q);
				//	f << qp[0] << " " << qp[1] << " " << qp[2] << "\n";
				//}
				//f << "\nLINES 1 " << (line.size() + 1) << "\n";
				//f << line.size() << " ";
				//for (int q = 0; q < line.size(); q++) 
				//	f << q << " ";
				
				//f << "\n";
				//f.close();


			}
		}
	}
}

// gradient interpolation near planes
P_float* GetGradEnergies(P_float* E,int nb_points,P_float* points,P_float* normals,vtkStructuredPoints* MRI_grad,P_float **transform_P,P_float **spacing_P,P_float tol,P_float step,int depth,int interpolationmode,P_float* ROI,bool oppositedirection)
{
if(MRI_grad==NULL) return NULL;
int i,j,index,pselect;
P_float p[3],n[3],n2[3];
P_float vals[3],val;
P_float *Normals=new P_float[3*nb_points];

if(MRI_grad->GetNumberOfScalarComponents()==1) // NG metric
	{
	for(i=0;i<nb_points;i++)
		{
		pselect=GetSelectP(n2,points+3*i,normals+3*i,MRI_grad,transform_P,spacing_P,tol);
		n[0]=n2[0]*step; n[1]=n2[1]*step; n[2]=n2[2]*step; 
		memcpy(Normals+3*i,n2,3*sizeof(P_float));
		p[0]=points[3*i]-(P_float)depth*n[0];		p[1]=points[3*i+1]-(P_float)depth*n[1];		p[2]=points[3*i+2]-(P_float)depth*n[2];
		index=i*(2*depth+1);
		if(pselect!=-1)
			for(j=0;j<2*depth+1;j++)
				{
				if(p[0]>ROI[0] && p[0]<ROI[1] && p[1]>ROI[2] && p[1]<ROI[3] && p[2]>ROI[4] && p[2]<ROI[5]) 
					{
					val=Interpolation(MRI_grad,p,pselect,transform_P,spacing_P,tol,interpolationmode);
					if(val!=-1) if(E[index]==1E10) E[index]=-(P_float)val; 
					}
				p[0]+=n[0];	p[1]+=n[1];	p[2]+=n[2];
				index++;
				}
		}
	}
else	// G metric
	{
	for(i=0;i<nb_points;i++)	
		{
		pselect=GetSelectP(n2,points+3*i,normals+3*i,MRI_grad,transform_P,spacing_P,tol);
		n[0]=n2[0]*step; n[1]=n2[1]*step; n[2]=n2[2]*step; 
		memcpy(Normals+3*i,n2,3*sizeof(P_float));
		p[0]=points[3*i]-(P_float)depth*n[0];		p[1]=points[3*i+1]-(P_float)depth*n[1];		p[2]=points[3*i+2]-(P_float)depth*n[2];
		index=i*(2*depth+1);
		if(pselect!=-1)
			for(j=0;j<2*depth+1;j++)
				{
				if(p[0]>ROI[0] && p[0]<ROI[1] && p[1]>ROI[2] && p[1]<ROI[3] && p[2]>ROI[4] && p[2]<ROI[5]) 
					{
					vInterpolation(vals,MRI_grad,p,pselect,transform_P,spacing_P,tol,interpolationmode);
					if(vals[0]!=1E5) {if(E[index]==1E10) if(oppositedirection) E[index]=(dotproduct(vals,normals+3*i)); else E[index]=-(dotproduct(vals,normals+3*i)); }
					}
				p[0]+=n[0];	p[1]+=n[1];	p[2]+=n[2];
				index++;
				}
		}
	}
return Normals;
}


void GetIPEnergies(P_float* E,int metric,int nbpoints,int **IP,int **IPref,int mn,int pn,int mnref,int pnref)
{
int depth=mn-mnref;
int i,j;
bool outside;
P_float val;
for(i=0;i<nbpoints;i++)
	{
	outside=false;
	for(j=0;j<(mn+pn+1);j++) 	if(IP[i][j]==-1E5) outside=true;
//	for(j=0;j<(mnref+pnref+1);j++) if(IPref[i][j]==-1E5) outside=true;
	if(!outside)
		for(j=0;j<2*depth+1;j++) 
			{
			val=GetIPSimilarity(metric,IP[i],IPref[i],mn,pn,mnref,pnref,j);
			if(val!=1E10) {if(E[i*(2*depth+1)+j]==1E10) E[i*(2*depth+1)+j]=val; else E[i*(2*depth+1)+j]+=val;}
			}
	}
}


void GetIPEnergies(P_float* E,int metric,int nbpoints,P_float **IP,P_float **IPref,int mn,int pn,int mnref,int pnref)
{
int depth=mn-mnref;
int i,j;
bool outside;
P_float val;
for(i=0;i<nbpoints;i++)
	{
	outside=false;
	for(j=0;j<3*(mn+pn+1);j++) 	if(IP[i][j]==-1E5) outside=true;
//	for(j=0;j<3*(mnref+pnref+1);j++) if(IPref[i][j]==-1E5) outside=true;
	if(!outside)
		for(j=0;j<2*depth+1;j++) 
			{
			val=GetIPSimilarity(metric,IP[i],IPref[i],mn,pn,mnref,pnref,j);
			if(val!=1E10) {if(E[i*(2*depth+1)+j]==1E10) E[i*(2*depth+1)+j]=val; else E[i*(2*depth+1)+j]+=val;}
			}
	}
}


P_float GetIPSimilarity(int metric,int *IP,int *IPref,int mn,int pn,int mnref,int pnref,int offset)
{
int j,size=pnref+mnref+1;
for(j=0;j<size;j++)  if(IP[j+offset]==-1E5 || IPref[j]==-1E5) return 1E10;
P_float S=0;

if(metric==0) // AD 
    {
    for(j=0;j<size;j++)  S+=fabs( (P_float)IP[j+offset] - (P_float)IPref[j]);
    S=S/(P_float)size;
    }
else if(metric==1) // NCC 
    {
	P_float m1=0,m2=0,norm1=0,norm2=0,diff1,diff2;
    for(j=0;j<size;j++)  {  m1+=IP[j+offset]; m2+=IPref[j];}   
    m1=m1/(P_float)size; m2=m2/(P_float)size;
    for(j=0;j<size;j++) 
        {
        diff1=IP[j+offset]-m1; diff2=IPref[j]-m2;
        S+=diff1*diff2; norm1+=diff1*diff1; norm2+=diff2*diff2;
        }  
    if(norm1!=0 && norm2!=0) S=-S/(sqrt(norm1*norm2)); else S=-1;
    }
else if(metric==2) // I 
    {
    S=(P_float)IP[j+offset];
    }
return S;
}


P_float GetIPSimilarity(int metric,P_float *IP,P_float *IPref,int mn,int pn,int mnref,int pnref,int offset)
{
int j,size=pnref+mnref+1;
for(j=0;j<size;j++)  if(IP[3*(j+offset)]==-1E5 || IPref[3*j]==-1E5) return 1E10;
P_float v[3],S=0;

if(metric==0) // AD 
    {
	for(j=0;j<size;j++) { v[0]=IP[3*(j+offset)] - IPref[3*j]; v[1]=IP[3*(j+offset)+1] - IPref[3*j+1]; v[2]=IP[3*(j+offset)+2] - IPref[3*j+2]; S+=dotproduct(v,v); }
    S=S/((P_float)size+0.0001);
    }
else if(metric==1 ) // NCC 
    {
	P_float m1[3]={0,0,0},m2[3]={0,0,0},norm1=0,norm2=0,diff1[3],diff2[3];
    for(j=0;j<size;j++)  {  m1[0]+=IP[3*(j+offset)]; m1[1]+=IP[3*(j+offset)+1]; m1[2]+=IP[3*(j+offset)+2]; m2[0]+=IPref[3*j]; m2[1]+=IPref[3*j+1]; m2[2]+=IPref[3*j+2];}   
    m1[0]=m1[0]/(P_float)size; m1[1]=m1[1]/(P_float)size; m1[2]=m1[2]/(P_float)size; m2[0]=m2[0]/(P_float)size; m2[1]=m2[1]/(P_float)size; m2[2]=m2[2]/(P_float)size;
    for(j=0;j<size;j++) 
        {
        diff1[0]=IP[3*(j+offset)]-m1[0];   diff1[1]=IP[3*(j+offset)+1]-m1[1];     diff1[2]=IP[3*(j+offset)+2]-m1[2]; 
		diff2[0]=IPref[3*j]-m2[0]; diff2[1]=IPref[3*j+1]-m2[1]; diff2[2]=IPref[3*j+2]-m2[2];
        S+=dotproduct(diff1,diff2); norm1+=dotproduct(diff1,diff1); norm2+=dotproduct(diff2,diff2);
        }  
    if(norm1!=0 && norm2!=0) S=-S/(sqrt(norm1*norm2)); else S=-1;
    }
return S;
}



// global similarity (rigid registration)
P_float GetIPSimilarity(int metric,int nbpoints,int **IP,int **IPref,int mn,int pn)
{
P_float S=0;
int i,j,size=pn+mn+1,count=0;
P_float m1=0,m2=0,dev1=0,norm1=0,norm2=0,norm=0,maxval1=0,maxval2=0; 
int nb_samples=10000,nbbins=256; 
P_float *j_buff=new P_float[nbbins], *i_buff=new P_float[nbbins], *buff=new P_float[nbbins*nbbins]; 

switch(metric)
	{
	case 0: // AD 
		for(i=0;i<nbpoints;i++)	for(j=0;j<size;j++)	if(IP[i][j]!=-1E5 && IPref[i][j]!=-1E5) {S+=fabs( (P_float)IP[i][j] - (P_float)IPref[i][j]); count++;}
		if(count!=0) S=S/(P_float)(size*count); else S=1E10;
	break;

	case 1: // NCC 
		for(i=0;i<nbpoints;i++) for(j=0;j<size;j++) if(IP[i][j]!=-1E5 && IPref[i][j]!=-1E5) {m1+=IP[i][j]; dev1+=IP[i][j]*IP[i][j]; m2+=IPref[i][j]; count++;}
		if(count!=0) {m1=m1/(P_float)count; m2=m2/(P_float)count; dev1=sqrt(dev1-(P_float)count*m1*m1);} else return 1E10;
		for(i=0;i<nbpoints;i++) for(j=0;j<size;j++) if(IP[i][j]!=-1E5 && IPref[i][j]!=-1E5)	S+=IP[i][j]*(IPref[i][j]-m2); 
		S=-S/dev1;
	break;

	case 2: // MI
		for(i=0;i<nbbins;i++) { *(i_buff+i)=0; *(j_buff+i)=0; }		for(i=0;i<nbbins*nbbins;i++) *(buff+i)=0; 
		for(i=0;i<nbpoints;i++) for(j=0;j<size;j++) if(IP[i][j]!=-1E5 && IPref[i][j]!=-1E5) {if(IP[i][j]>maxval1) maxval1=IP[i][j];  if(IPref[i][j]>maxval2) maxval2=IPref[i][j]; }
		if(maxval1==0 || maxval2==0) return 1E10;
		maxval1=(P_float)nbbins/(maxval1+1); maxval2=(P_float)nbbins/(maxval2+1);
		for(i=0;i<nbpoints;i++) for(j=0;j<size;j++)	if(IP[i][j]!=-1E5 && IPref[i][j]!=-1E5)
				{
				*(buff+(int)floor((P_float)IPref[i][j]*maxval2)+ (int)floor((P_float)IP[i][j]*maxval1)*(int)nbbins)+=1;
				*(i_buff+(int)floor((P_float)IP[i][j]*maxval1))+=1;
				*(j_buff+(int)floor((P_float)IPref[i][j]*maxval2))+=1;
				}
		for(i=0;i<nbbins;i++) { if(*(i_buff+i)!=0) norm1+=*(i_buff+i)*log10(*(i_buff+i)); if(*(j_buff+i)!=0) norm2+=*(j_buff+i)*log10(*(j_buff+i));}
		for(i=0;i<nbbins;i++) for(j=0;j<nbbins;j++) if(*(buff+j+i*nbbins)!=0) norm+=*(buff+j+i*nbbins)*log10(*(buff+j+i*nbbins));
		S=(norm1+norm2)/norm;
	break;
	}
free(j_buff); free(i_buff); free(buff); 
return S;
}


P_float GetGradEnergy(int nb_points,P_float* points,P_float* normals,vtkStructuredPoints* MRI_grad, int interpolationmode,P_float *OffsetTransform)
{
if(MRI_grad==NULL) return 1E10;
int i,count=0;
P_float p[3],n[3],S=0;
P_float vals[3],val;
if(MRI_grad->GetNumberOfScalarComponents()==1) // NG metric
	{
	for(i=0;i<nb_points;i++)
		{
		if(OffsetTransform!=NULL)	{Transform(points+3*i,p,OffsetTransform); }	else  {memcpy(p,points+3*i,3*sizeof(P_float)); }
		val=Interpolation(MRI_grad,p,interpolationmode);	if(val!=-1) {S-=(P_float)val;  count++;}
		}
	}
else	// G metric
	{
	for(i=0;i<nb_points;i++)	
		{
		if(OffsetTransform!=NULL)	{Transform(points+3*i,p,OffsetTransform); Transform_R(normals+3*i,n,OffsetTransform);	}	else  {memcpy(p,points+3*i,3*sizeof(P_float)); memcpy(n,normals+3*i,3*sizeof(P_float));  }
		vInterpolation(vals,MRI_grad,p,interpolationmode);
		if(vals[0]!=1E5) {S-=dotproduct(vals,n); count++;}
		}
	}
S=S/(P_float)count;
return S;
}



//gl/gradient interpolation in cartesian volume
	//scalar
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,int** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI,P_float *OffsetTransform)
{
int i,j;
double p[3],n[3],interp; 

P_float pf[3]; 

for(i=0;i<nbpoints;i++)  
    {
	if(OffsetTransform!=NULL)	{Transform(pts+3*i,p,OffsetTransform); Transform_R(normals+3*i,n,OffsetTransform);	}	else  {memcpy(p,pts+3*i,3*sizeof(P_float)); memcpy(n,normals+3*i,3*sizeof(P_float));  }
	n[0]=n[0]*s; n[1]=n[1]*s; n[2]=n[2]*s; 
	for(j=-mn;j<=pn;j++)  
		{
		pf[0]=p[0]+j*n[0];	pf[1]=p[1]+j*n[1];	pf[2]=p[2]+j*n[2];
		if(pf[0]>ROI[0] && pf[0]<ROI[1] && pf[1]>ROI[2] && pf[1]<ROI[3] && pf[2]>ROI[4] && pf[2]<ROI[5])
			{
			interp=Interpolation(vol,pf,interpolationmode); 
			if(interp!=-1) IP[i][j+mn]=(int)interp;
			}
		}
	}
}
	

	//vector
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI,P_float *OffsetTransform)
{
int i,j;
P_float p[3],n[3]; 
P_float vals[3];

P_float pf[3]; 

for(i=0;i<nbpoints;i++)  
	{
	if(OffsetTransform!=NULL)	{Transform(pts+3*i,p,OffsetTransform); Transform_R(normals+3*i,n,OffsetTransform);	}	else  {memcpy(p,pts+3*i,3*sizeof(P_float)); memcpy(n,normals+3*i,3*sizeof(P_float));  }
	n[0]=n[0]*s; n[1]=n[1]*s; n[2]=n[2]*s; 
	for(j=-mn;j<=pn;j++)  
		{
		pf[0]=p[0]+j*n[0];	pf[1]=p[1]+j*n[1];	pf[2]=p[2]+j*n[2];
		if(pf[0]>ROI[0] && pf[0]<ROI[1] && pf[1]>ROI[2] && pf[1]<ROI[3] && pf[2]>ROI[4] && pf[2]<ROI[5])
			{
			vInterpolation(vals,vol,pf,interpolationmode); 
			if(vals[0]!=1E5) {memcpy(IP[i]+3*(j+mn),vals,3*sizeof(P_float)); }
			}
		}
	}
}


	//+ smoothing
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,P_float *tangents1,P_float *tangents2,vtkStructuredPoints* vol,int** IP,int mn,int pn,P_float s,P_float** IP_GaussianWeights,int IP_GaussianR,int interpolationmode,P_float *ROI,P_float *OffsetTransform)
{
int i,j;
double t1[3],t2[3],p[3],n[3],interp; 

P_float pf[3]; 

for(i=0;i<nbpoints;i++)  
    {
	if(OffsetTransform!=NULL)	{Transform(pts+3*i,p,OffsetTransform); Transform_R(normals+3*i,n,OffsetTransform);	}	else  {memcpy(p,pts+3*i,3*sizeof(P_float)); memcpy(n,normals+3*i,3*sizeof(P_float));  }
	n[0]=n[0]*s; n[1]=n[1]*s; n[2]=n[2]*s; 
	for(j=-mn;j<=pn;j++)  
		{
		pf[0]=p[0]+j*n[0];	pf[1]=p[1]+j*n[1];	pf[2]=p[2]+j*n[2];
		if(pf[0]>ROI[0] && pf[0]<ROI[1] && pf[1]>ROI[2] && pf[1]<ROI[3] && pf[2]>ROI[4] && pf[2]<ROI[5])
			{
			interp=InterpolationSmooth(vol,pf,t1,t2,IP_GaussianWeights,IP_GaussianR,interpolationmode); 
			if(interp!=-1) IP[i][j+mn]=(int)interp; 
			}
		}
	}
}


//gl interpolation near planes
P_float* ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol,int** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI)
{
int i,j,pselect;
double n[3],n2[3],interp; 

P_float pf[3]; 
P_float *Normals=new P_float[3*nbpoints];

for(i=0;i<nbpoints;i++)  
    {
	pselect=GetSelectP(n2,pts+3*i,normals+3*i,vol,transform_P,spacing_P,tol);
	n[0]=n2[0]*s; n[1]=n2[1]*s; n[2]=n2[2]*s; 
	memcpy(Normals+3*i,n2,3*sizeof(P_float));
	if(pselect!=-1)
		for(j=-mn;j<=pn;j++)  
			{
			pf[0]=pts[3*i]+j*n[0];	pf[1]=pts[3*i+1]+j*n[1];	pf[2]=pts[3*i+2]+j*n[2];
			if(pf[0]>ROI[0] && pf[0]<ROI[1] && pf[1]>ROI[2] && pf[1]<ROI[3] && pf[2]>ROI[4] && pf[2]<ROI[5])
				{
				interp=Interpolation(vol,pf,pselect,transform_P,spacing_P,tol,interpolationmode); 
				if(interp!=-1) IP[i][j+mn]=(int)interp; 
				}
			}
	}
return Normals;
}

P_float* ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol,P_float** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI)
{
int i,j,pselect;
P_float n[3],n2[3],vals[3]; 
P_float pf[3]; 
P_float *Normals=new P_float[3*nbpoints];

for(i=0;i<nbpoints;i++)  
    {
	pselect=GetSelectP(n2,pts+3*i,normals+3*i,vol,transform_P,spacing_P,tol);
	n[0]=n2[0]*s; n[1]=n2[1]*s; n[2]=n2[2]*s; 
	memcpy(Normals+3*i,n2,3*sizeof(P_float));
	if(pselect!=-1)
//	if(pselect==1 || pselect==4 || pselect==5) // hack for testing
		for(j=-mn;j<=pn;j++)  
			{
			pf[0]=pts[3*i]+j*n[0];	pf[1]=pts[3*i+1]+j*n[1];	pf[2]=pts[3*i+2]+j*n[2];
			if(pf[0]>ROI[0] && pf[0]<ROI[1] && pf[1]>ROI[2] && pf[1]<ROI[3] && pf[2]>ROI[4] && pf[2]<ROI[5])
				{
				vInterpolation(vals,vol,pf,pselect,transform_P,spacing_P,tol,interpolationmode); 
				if(vals[0]!=1E5)  {memcpy(IP[i]+3*(j+mn),vals,3*sizeof(P_float)); }
				}
			}
	}
return Normals;
}

//gl/gradient interpolation in cylindric volume
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float *transform_P,P_float *spacing_P,int** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI)
{
int i,j;
double *p,n[3],interp; 
P_float pf[3];

for(i=0;i<nbpoints;i++)  
    {
	p=pts+3*i;	n[0]=normals[3*i]*s; n[1]=normals[3*i+1]*s; n[2]=normals[3*i+2]*s; 
    for(j=-mn;j<=pn;j++)  
		{
		pf[0]=p[0]+j*n[0];	pf[1]=p[1]+j*n[1];	pf[2]=p[2]+j*n[2];
		if(pf[0]>ROI[0] && pf[0]<ROI[1] && pf[1]>ROI[2] && pf[1]<ROI[3] && pf[2]>ROI[4] && pf[2]<ROI[5])
			{
			interp=Interpolation(vol,pf,transform_P,spacing_P,interpolationmode); 
			if(interp!=-1) IP[i][j+mn]=(int)interp; 
			}
		}
	}
}

void ProcessIP(int nb_points,int** IP,int mn,int pn,int* IPmask,int IP_masksize) {ProcessIP(IP,nb_points,IP,mn,pn,IPmask,IP_masksize);}
void ProcessIP(int** IPproc,int nb_points,int** IP,int mn,int pn,int* IPmask,int IP_masksize)
{
if(IPmask==NULL) return;

int i,j,k,IPsize=mn+pn+1;
int offset_mask=IP_masksize/2,ind;
int* proftemp=new int[IPsize];
bool flag;

for(i=0;i<nb_points;i++)
    {
    for(j=0;j<IPsize;j++)
        {
		flag=false;
		proftemp[j]=0;
		for(k=0;k<IP_masksize;k++)
			{
			ind=k-offset_mask+j;
			if(ind<0) {proftemp[j]+=IPmask[k] * IP[i][0]; if(IP[i][0]==-1E5) flag=true;}
			else if(ind>=IPsize) {proftemp[j]+=IPmask[k] * IP[i][IPsize-1]; if(IP[i][IPsize-1]==-1E5) flag=true;}
			else {proftemp[j]+=IPmask[k] * IP[i][ind]; if(IP[i][ind]==-1E5) flag=true;}
			}
		if(flag) proftemp[j]=-1E5;
		}
    memcpy(IPproc[i],proftemp,IPsize*sizeof(int));
    }
free(proftemp);
}




void ProcessIP(int nb_points,P_float** IP,int mn,int pn,int* IPmask,int IP_masksize) {ProcessIP(IP,nb_points,IP,mn,pn,IPmask,IP_masksize);}
void ProcessIP(P_float** IPproc,int nb_points,P_float** IP,int mn,int pn,int* IPmask,int IP_masksize)
{
if(IPmask==NULL) return;

int i,j,k,IPsize=mn+pn+1;
int offset_mask=IP_masksize/2,ind;
P_float* proftemp=new P_float[3*IPsize];
bool flag;

for(i=0;i<nb_points;i++)
    {
    for(j=0;j<IPsize;j++)
        {
		flag=false;
		proftemp[j]=0;
		for(k=0;k<IP_masksize;k++)
			{
			ind=k-offset_mask+j;
			if(ind<0) {proftemp[3*j]+=(P_float)IPmask[k] * IP[i][0];	proftemp[3*j+1]+=(P_float)IPmask[k] * IP[i][1];	proftemp[3*j+2]+=(P_float)IPmask[k] * IP[i][2]; if(IP[i][0]==-1E5) flag=true;}
			else if(ind>=IPsize) {proftemp[3*j]+=(P_float)IPmask[k] * IP[i][3*(IPsize-1)]; proftemp[3*j+1]+=(P_float)IPmask[k] * IP[i][3*(IPsize-1)+1]; proftemp[3*j+2]+=(P_float)IPmask[k] * IP[i][3*(IPsize-1)+2]; if(IP[i][3*(IPsize-1)]==-1E5) flag=true;}
			else {proftemp[3*j]+=(P_float)IPmask[k] * IP[i][3*ind];	proftemp[3*j+1]+=(P_float)IPmask[k] * IP[i][3*ind+1];	proftemp[3*j+2]+=(P_float)IPmask[k] * IP[i][3*ind+2];	if(IP[i][3*ind]==-1E5) flag=true;}
			}
		if(flag) {proftemp[3*j]=-1E5; proftemp[3*j+1]=-1E5; proftemp[3*j+2]=-1E5;}
		}
    memcpy(IPproc[i],proftemp,3*IPsize*sizeof(P_float));
    }
free(proftemp);
}



inline P_float Get3x3Determinant(const P_float M[3][3]) {return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])+M[0][1]*(M[1][2]*M[2][0]-M[1][0]*M[2][2])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);}

P_float Get3x3Determinant(const P_float M[4][4],const int x,const int y)
{
P_float M2[3][3]; int ci=0,cj=0;
for(int i=0;i<4;i++) if(i!=x) for(int j=0;j<4;j++) if(j!=y) {M2[ci][cj]=M[i][j]; ci++;cj++;}
return Get3x3Determinant(M2)*pow((double)-1,x+y);
}

P_float Get4x4Determinant(const P_float M[4][4])
{
P_float ret=0;
for(int j=0;j<4;j++) ret+=M[0][j]*Get3x3Determinant(M,0,j);
return ret;
}

void GetThreadCoefficientsFromLocalCoordinates(P_float R[12],const P_float P0[3],const P_float P1[3],const P_float P2[3],const P_float P3[3])
{
P_float M[4][4]={{1,P0[0],P0[1],P0[2]},{1,P1[0],P1[1],P1[2]},{1,P2[0],P2[1],P2[2]},{1,P3[0],P3[1],P3[2]}};
P_float d=Get4x4Determinant(M);
for(int i=0;i<3;i++) for(int j=0;j<4;j++) R[i+3*j]=Get3x3Determinant(M,i,j)/d;
}

void GetThreadCoefficients(P_float R[12],const P_float Pproj[3],const P_float P1[3],const P_float P2[3],const P_float P3[3],const P_float n[3],const P_float h,const P_float vol)
{
P_float u[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]};
P_float nu=norm(u); u[0]=u[0]/nu; u[1]=u[1]/nu; u[2]=u[2]/nu;
P_float v[3]; crossproduct(v,u,n);
P_float det=u[0]*v[1]-u[1]*v[0];
P_float M[4]={v[1]/det,-v[0]/det,-u[1]/det,u[0]/det};
P_float p1[2]={(P1[0]-Pproj[0])*M[0]+(P1[1]-Pproj[1])*M[1],(P1[0]-Pproj[0])*M[2]+(P1[1]-Pproj[1])*M[3]};
P_float p2[2]={(P2[0]-Pproj[0])*M[0]+(P2[1]-Pproj[1])*M[1],(P2[0]-Pproj[0])*M[2]+(P2[1]-Pproj[1])*M[3]};
P_float p3[2]={(P3[0]-Pproj[0])*M[0]+(P3[1]-Pproj[1])*M[1],(P3[0]-Pproj[0])*M[2]+(P3[1]-Pproj[1])*M[3]};
R[3]=h*(p3[1]-p2[1])/vol;	R[4]=h*(p2[0]-p3[0])/vol;	R[5]=(p3[1]*p2[0]-p3[0]*p2[1])/vol;
R[6]=h*(p1[1]-p3[1])/vol;	R[7]=h*(p3[0]-p1[0])/vol;	R[8]=(p1[1]*p3[0]-p1[0]*p3[1])/vol;
R[9]=h*(p2[1]-p1[1])/vol;	R[10]=h*(p1[0]-p2[0])/vol;	R[11]=(p2[1]*p1[0]-p2[0]*p1[1])/vol;
R[0]=0;						R[1]=0;						R[2]=-R[5]-R[8]-R[11];
}


void GetThreadVectors(P_float U[3],P_float V[3],P_float W[3],const P_float R[12],const P_float P0[3],const P_float P1[3],const P_float P2[3],const P_float P3[3])
{
U[0]=(R[0]*P0[0]+R[3]*P1[0]+R[6]*P2[0]+R[9]*P3[0]);
U[1]=(R[0]*P0[1]+R[3]*P1[1]+R[6]*P2[1]+R[9]*P3[1]);
U[2]=(R[0]*P0[2]+R[3]*P1[2]+R[6]*P2[2]+R[9]*P3[2]);
V[0]=(R[1]*P0[0]+R[4]*P1[0]+R[7]*P2[0]+R[10]*P3[0]);
V[1]=(R[1]*P0[1]+R[4]*P1[1]+R[7]*P2[1]+R[10]*P3[1]);
V[2]=(R[1]*P0[2]+R[4]*P1[2]+R[7]*P2[2]+R[10]*P3[2]);
W[0]=(R[2]*P0[0]+R[5]*P1[0]+R[8]*P2[0]+R[11]*P3[0]);
W[1]=(R[2]*P0[1]+R[5]*P1[1]+R[8]*P2[1]+R[11]*P3[1]);
W[2]=(R[2]*P0[2]+R[5]*P1[2]+R[8]*P2[2]+R[11]*P3[2]);
}

P_float GetModePatternUU(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12])
{
P_float nrm=norm(U);
P_float u[3]={U[0]/nrm,U[1]/nrm,U[2]/nrm};

H[0][0]=R[0]*u[0]; H[0][1]=R[0]*u[1]; H[0][2]=R[0]*u[2];
H[1][0]=R[3]*u[0]; H[1][1]=R[3]*u[1]; H[1][2]=R[3]*u[2];
H[2][0]=R[6]*u[0]; H[2][1]=R[6]*u[1]; H[2][2]=R[6]*u[2];
H[3][0]=R[9]*u[0]; H[3][1]=R[9]*u[1]; H[3][2]=R[9]*u[2];

return nrm-1;
}

P_float GetModePatternVV(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12])
{
P_float nrm=norm(V);
P_float u[3]={V[0]/nrm,V[1]/nrm,V[2]/nrm};

H[0][0]=R[1]*u[0]; H[0][1]=R[1]*u[1]; H[0][2]=R[1]*u[2];
H[1][0]=R[4]*u[0]; H[1][1]=R[4]*u[1]; H[1][2]=R[4]*u[2];
H[2][0]=R[7]*u[0]; H[2][1]=R[7]*u[1]; H[2][2]=R[7]*u[2];
H[3][0]=R[10]*u[0]; H[3][1]=R[10]*u[1]; H[3][2]=R[10]*u[2];

return nrm-1;
}

P_float GetModePatternWW(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12])
{
P_float nrm=norm(W);
P_float u[3]={W[0]/nrm,W[1]/nrm,W[2]/nrm};

H[0][0]=R[2]*u[0]; H[0][1]=R[2]*u[1]; H[0][2]=R[2]*u[2];
H[1][0]=R[5]*u[0]; H[1][1]=R[5]*u[1]; H[1][2]=R[5]*u[2];
H[2][0]=R[8]*u[0]; H[2][1]=R[8]*u[1]; H[2][2]=R[8]*u[2];
H[3][0]=R[11]*u[0]; H[3][1]=R[11]*u[1]; H[3][2]=R[11]*u[2];

return nrm-1;
}

P_float GetModePatternUV(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12])
{
P_float UpV[3]={U[0]+V[0],U[1]+V[1],U[2]+V[2]};
P_float UmV[3]={U[0]-V[0],U[1]-V[1],U[2]-V[2]};
P_float nrmUpV=norm(UpV),nrmUmV=norm(UmV);
P_float upv[3]={UpV[0]/nrmUpV,UpV[1]/nrmUpV,UpV[2]/nrmUpV},umv[3]={UmV[0]/nrmUmV,UmV[1]/nrmUmV,UmV[2]/nrmUmV};

H[0][0]=DSQRT2*((R[0]+R[1])*upv[0]-(R[0]-R[1])*umv[0]); H[0][1]=DSQRT2*((R[0]+R[1])*upv[1]-(R[0]-R[1])*umv[1]); H[0][2]=DSQRT2*((R[0]+R[1])*upv[2]-(R[0]-R[1])*umv[2]);
H[1][0]=DSQRT2*((R[3]+R[4])*upv[0]-(R[3]-R[4])*umv[0]); H[1][1]=DSQRT2*((R[3]+R[4])*upv[1]-(R[3]-R[4])*umv[1]); H[1][2]=DSQRT2*((R[3]+R[4])*upv[2]-(R[3]-R[4])*umv[2]);
H[2][0]=DSQRT2*((R[6]+R[7])*upv[0]-(R[6]-R[7])*umv[0]); H[2][1]=DSQRT2*((R[6]+R[7])*upv[1]-(R[6]-R[7])*umv[1]); H[2][2]=DSQRT2*((R[6]+R[7])*upv[2]-(R[6]-R[7])*umv[2]);
H[3][0]=DSQRT2*((R[9]+R[10])*upv[0]-(R[9]-R[10])*umv[0]); H[3][1]=DSQRT2*((R[9]+R[10])*upv[1]-(R[9]-R[10])*umv[1]); H[3][2]=DSQRT2*((R[9]+R[10])*upv[2]-(R[9]-R[10])*umv[2]);

return DSQRT2*(nrmUpV-nrmUmV);
}

P_float GetModePatternUW(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12])
{
P_float UpV[3]={U[0]+W[0],U[1]+W[1],U[2]+W[2]};
P_float UmV[3]={U[0]-W[0],U[1]-W[1],U[2]-W[2]};
P_float nrmUpV=norm(UpV),nrmUmV=norm(UmV);
P_float upv[3]={UpV[0]/nrmUpV,UpV[1]/nrmUpV,UpV[2]/nrmUpV},umv[3]={UmV[0]/nrmUmV,UmV[1]/nrmUmV,UmV[2]/nrmUmV};

H[0][0]=DSQRT2*((R[0]+R[2])*upv[0]-(R[0]-R[2])*umv[0]); H[0][1]=DSQRT2*((R[0]+R[2])*upv[1]-(R[0]-R[2])*umv[1]); H[0][2]=DSQRT2*((R[0]+R[2])*upv[2]-(R[0]-R[2])*umv[2]);
H[1][0]=DSQRT2*((R[3]+R[5])*upv[0]-(R[3]-R[5])*umv[0]); H[1][1]=DSQRT2*((R[3]+R[5])*upv[1]-(R[3]-R[5])*umv[1]); H[1][2]=DSQRT2*((R[3]+R[5])*upv[2]-(R[3]-R[5])*umv[2]);
H[2][0]=DSQRT2*((R[6]+R[8])*upv[0]-(R[6]-R[8])*umv[0]); H[2][1]=DSQRT2*((R[6]+R[8])*upv[1]-(R[6]-R[8])*umv[1]); H[2][2]=DSQRT2*((R[6]+R[8])*upv[2]-(R[6]-R[8])*umv[2]);
H[3][0]=DSQRT2*((R[9]+R[11])*upv[0]-(R[9]-R[11])*umv[0]); H[3][1]=DSQRT2*((R[9]+R[11])*upv[1]-(R[9]-R[11])*umv[1]); H[3][2]=DSQRT2*((R[9]+R[11])*upv[2]-(R[9]-R[11])*umv[2]);

return DSQRT2*(nrmUpV-nrmUmV);
}

P_float GetModePatternVW(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12])
{
P_float UpV[3]={V[0]+W[0],V[1]+W[1],V[2]+W[2]};
P_float UmV[3]={V[0]-W[0],V[1]-W[1],V[2]-W[2]};
P_float nrmUpV=norm(UpV),nrmUmV=norm(UmV);
P_float upv[3]={UpV[0]/nrmUpV,UpV[1]/nrmUpV,UpV[2]/nrmUpV},umv[3]={UmV[0]/nrmUmV,UmV[1]/nrmUmV,UmV[2]/nrmUmV};

H[0][0]=DSQRT2*((R[1]+R[2])*upv[0]-(R[1]-R[2])*umv[0]); H[0][1]=DSQRT2*((R[1]+R[2])*upv[1]-(R[1]-R[2])*umv[1]); H[0][2]=DSQRT2*((R[1]+R[2])*upv[2]-(R[1]-R[2])*umv[2]);
H[1][0]=DSQRT2*((R[4]+R[5])*upv[0]-(R[4]-R[5])*umv[0]); H[1][1]=DSQRT2*((R[4]+R[5])*upv[1]-(R[4]-R[5])*umv[1]); H[1][2]=DSQRT2*((R[4]+R[5])*upv[2]-(R[4]-R[5])*umv[2]);
H[2][0]=DSQRT2*((R[7]+R[8])*upv[0]-(R[7]-R[8])*umv[0]); H[2][1]=DSQRT2*((R[7]+R[8])*upv[1]-(R[7]-R[8])*umv[1]); H[2][2]=DSQRT2*((R[7]+R[8])*upv[2]-(R[7]-R[8])*umv[2]);
H[3][0]=DSQRT2*((R[10]+R[11])*upv[0]-(R[10]-R[11])*umv[0]); H[3][1]=DSQRT2*((R[10]+R[11])*upv[1]-(R[10]-R[11])*umv[1]); H[3][2]=DSQRT2*((R[10]+R[11])*upv[2]-(R[10]-R[11])*umv[2]);

return DSQRT2*(nrmUpV-nrmUmV);
}

P_float GetStress(P_float elongation,P_float a) {return a*elongation;}
P_float GetStressDerivative(P_float elongation,P_float a) {return a;}


void GetBoundingBox(vtkStructuredPoints* vol, P_float *offsettransform, P_float bounds[3][2]) {
	if(offsettransform==NULL) 	{
		bounds[0][0] = vol->GetOrigin()[0]; 	
		bounds[1][0] = vol->GetOrigin()[1];	
		bounds[2][0] = vol->GetOrigin()[2];

		bounds[0][1] = vol->GetOrigin()[0] + (P_float)(vol->GetDimensions()[0] - 1) * vol->GetSpacing()[0]; 
		bounds[1][1] = vol->GetOrigin()[1] + (P_float)(vol->GetDimensions()[1] - 1) * vol->GetSpacing()[1]; 
		bounds[2][1] = vol->GetOrigin()[2] + (P_float)(vol->GetDimensions()[2] - 1) * vol->GetSpacing()[2]; 	
	}
	else
	{
		P_float pt[3];
		bounds[0][0]=1E10; 	bounds[1][0]=1E10;	bounds[2][0]=1E10;	bounds[0][1]=-1E10; bounds[1][1]=-1E10; bounds[2][1]=-1E10; 	
		pt[0]=vol->GetOrigin()[0];																pt[1]=vol->GetOrigin()[1];																pt[2]=vol->GetOrigin()[2]; 																if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0];																pt[1]=vol->GetOrigin()[1]+(P_float)(vol->GetDimensions()[1]-1)*vol->GetSpacing()[1];	pt[2]=vol->GetOrigin()[2]; 																if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0];																pt[1]=vol->GetOrigin()[1];																pt[2]=vol->GetOrigin()[2]+(P_float)(vol->GetDimensions()[2]-1)*vol->GetSpacing()[2]; 	if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0];																pt[1]=vol->GetOrigin()[1]+(P_float)(vol->GetDimensions()[1]-1)*vol->GetSpacing()[1];	pt[2]=vol->GetOrigin()[2]+(P_float)(vol->GetDimensions()[2]-1)*vol->GetSpacing()[2]; 	if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0]+(P_float)(vol->GetDimensions()[0]-1)*vol->GetSpacing()[0]; 	pt[1]=vol->GetOrigin()[1];																pt[2]=vol->GetOrigin()[2]; 																if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0]+(P_float)(vol->GetDimensions()[0]-1)*vol->GetSpacing()[0]; 	pt[1]=vol->GetOrigin()[1]+(P_float)(vol->GetDimensions()[1]-1)*vol->GetSpacing()[1];	pt[2]=vol->GetOrigin()[2]; 																if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0]+(P_float)(vol->GetDimensions()[0]-1)*vol->GetSpacing()[0]; 	pt[1]=vol->GetOrigin()[1];																pt[2]=vol->GetOrigin()[2]+(P_float)(vol->GetDimensions()[2]-1)*vol->GetSpacing()[2]; 	if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
		pt[0]=vol->GetOrigin()[0]+(P_float)(vol->GetDimensions()[0]-1)*vol->GetSpacing()[0]; 	pt[1]=vol->GetOrigin()[1]+(P_float)(vol->GetDimensions()[1]-1)*vol->GetSpacing()[1];	pt[2]=vol->GetOrigin()[2]+(P_float)(vol->GetDimensions()[2]-1)*vol->GetSpacing()[2]; 	if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
	}
}

//rt
void GetBoundingBox(vtkStructuredPoints* vol,P_float **Transform,P_float **Spacing,P_float bounds[3][2])
{
P_float pt[3],u[3],v[3];
bounds[0][0]=1E10; 	bounds[1][0]=1E10;	bounds[2][0]=1E10;	bounds[0][1]=-1E10; bounds[1][1]=-1E10; bounds[2][1]=-1E10; 	
for(int i=0;i<vol->GetDimensions()[2];i++)
	{
	u[0]=(P_float)(vol->GetDimensions()[0]-1)*Spacing[i][0]*Transform[i][0];	u[1]=(P_float)(vol->GetDimensions()[0]-1)*Spacing[i][0]*Transform[i][4];	u[2]=(P_float)(vol->GetDimensions()[0]-1)*Spacing[i][0]*Transform[i][8];
	v[0]=(P_float)(vol->GetDimensions()[1]-1)*Spacing[i][1]*Transform[i][1]; 	v[1]=(P_float)(vol->GetDimensions()[1]-1)*Spacing[i][1]*Transform[i][5];	v[2]=(P_float)(vol->GetDimensions()[1]-1)*Spacing[i][1]*Transform[i][9];
	pt[0]=Transform[i][3];																pt[1]=Transform[i][7];																pt[2]=Transform[i][11]; 																if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
	pt[0]=Transform[i][3]+u[0];														pt[1]=Transform[i][7]+u[1];														pt[2]=Transform[i][11]+u[2]; 															if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
	pt[0]=Transform[i][3]+v[0];														pt[1]=Transform[i][7]+v[1];														pt[2]=Transform[i][11]+v[2]; 															if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
	pt[0]=Transform[i][3]+u[0]+v[0];													pt[1]=Transform[i][7]+u[1]+v[1];													pt[2]=Transform[i][11]+u[2]+v[2]; 														if(pt[0]<bounds[0][0]) bounds[0][0]=pt[0]; 	if(pt[1]<bounds[1][0]) bounds[1][0]=pt[1];	if(pt[2]<bounds[2][0]) bounds[2][0]=pt[2];	if(pt[0]>bounds[0][1]) bounds[0][1]=pt[0]; 	if(pt[1]>bounds[1][1]) bounds[1][1]=pt[1];	if(pt[2]>bounds[2][1]) bounds[2][1]=pt[2];
	}
}


void GetBoundingBox(vtkPolyData* model,P_float bounds[3][2])
{
double p[3];
model->GetPoint(0,p);
bounds[0][0]=p[0]; bounds[0][1]=p[0]; bounds[1][0]=p[1]; bounds[1][1]=p[1]; bounds[2][0]=p[2]; bounds[2][1]=p[2];
for(int i=1;i<model->GetNumberOfPoints();i++)
	{
	model->GetPoint(i,p);
	if(p[0]<bounds[0][0]) bounds[0][0]=p[0]; 	if(p[0]>bounds[0][1]) bounds[0][1]=p[0];
	if(p[1]<bounds[1][0]) bounds[1][0]=p[1]; 	if(p[1]>bounds[1][1]) bounds[1][1]=p[1];
	if(p[2]<bounds[2][0]) bounds[2][0]=p[2]; 	if(p[2]>bounds[2][1]) bounds[2][1]=p[2];
	}
}


int GetClosestPoint(vtkPolyData* model,vtkPointLocator* PointsLocator,P_float p[3],bool* selectcells,P_float n[3])
{
int i,ID=-1;
P_float d;
if(selectcells==NULL && n==NULL) {ID=PointsLocator->FindClosestPointWithinRadius(DMAX_ICP,p,d); }
else if(selectcells==NULL)	// closest pt according to a normal direction
	{
	P_float dmin=1E10,dp;
	for(i=0;i<model->GetNumberOfPoints();i++)
		{
		dp=n[0]*model->GetPointData()->GetNormals()->GetComponent(i,0)+n[1]*model->GetPointData()->GetNormals()->GetComponent(i,1)+n[2]*model->GetPointData()->GetNormals()->GetComponent(i,2);
		if(dp<0) {d=dist3D(p,model->GetPoint(i)); if(d<dmin) {dmin=d; ID=i;}}
		}
	}
else	// closest pt according to a list of slected cells
	{
	P_float dmin=1E10; int pt;
	for(i=0;i<model->GetNumberOfCells();i++)
		if(selectcells[i])	
			{
			pt=model->GetCell(i)->GetPointIds()->GetId(0); d=dist3D(p,model->GetPoint(pt)); if(d<dmin) {dmin=d; ID=pt;}
			pt=model->GetCell(i)->GetPointIds()->GetId(1); d=dist3D(p,model->GetPoint(pt)); if(d<dmin) {dmin=d; ID=pt;}
			pt=model->GetCell(i)->GetPointIds()->GetId(2); d=dist3D(p,model->GetPoint(pt)); if(d<dmin) {dmin=d; ID=pt;}
			}
	}
return ID;
}

void GetClosest(vtkPolyData* model,P_float M[16],P_float M_inv[16],vtkIdList** PointCells,P_float* CellNormals,vtkPointLocator* PointsLocator,CLOSESTSTRUCT* closest,P_float p[3],bool* selectcells,P_float n[3])
{
if(model->GetNumberOfPoints()==0) { closest->nb=0; closest->side=0; return;}

P_float tp[3]; Transform(p,tp,M_inv);
P_float tn[3]; if(n!=NULL) Transform_R(n,tn,M_inv); 
int i,index_cell,ids[3],side;
P_float fa,fb,fc,fd,fe,ff,p1[3],p2[3],p3[3],p1p2[3],p1p3[3],pp1[3],c2[3],*norml;

int ID;
if(n!=NULL) ID=GetClosestPoint(model,PointsLocator,tp,selectcells,tn); 
else ID=GetClosestPoint(model,PointsLocator,tp,selectcells,NULL); 
if(ID==-1) 	{ closest->nb=0; closest->side=0; return;}

for(i=0;i<PointCells[ID]->GetNumberOfIds();i++)
	{
	index_cell=PointCells[ID]->GetId(i);
	norml=CellNormals+3*index_cell;

	ids[0]=model->GetCell(index_cell)->GetPointIds()->GetId(0);		ids[1]=model->GetCell(index_cell)->GetPointIds()->GetId(1);			ids[2]=model->GetCell(index_cell)->GetPointIds()->GetId(2);
	p1[0]=model->GetPoint(ids[0])[0]; p1[1]=model->GetPoint(ids[0])[1]; p1[2]=model->GetPoint(ids[0])[2];	p2[0]=model->GetPoint(ids[1])[0]; p2[1]=model->GetPoint(ids[1])[1]; p2[2]=model->GetPoint(ids[1])[2];	p3[0]=model->GetPoint(ids[2])[0]; p3[1]=model->GetPoint(ids[2])[1]; p3[2]=model->GetPoint(ids[2])[2];

	p1p2[0]=p2[0]-p1[0];p1p2[1]=p2[1]-p1[1];p1p2[2]=p2[2]-p1[2];
	p1p3[0]=p3[0]-p1[0];p1p3[1]=p3[1]-p1[1];p1p3[2]=p3[2]-p1[2];
	pp1[0]=p1[0]-tp[0];pp1[1]=p1[1]-tp[1];pp1[2]=p1[2]-tp[2];

	if(dotproduct(pp1,norml)<=0) side=1; else side=-1;
	fa=dotproduct(p1p2,p1p2);fb=dotproduct(p1p2,p1p3);fc=dotproduct(p1p3,p1p3);fd=dotproduct(p1p2,pp1);fe=dotproduct(p1p3,pp1);ff=dotproduct(pp1,pp1);
	P_float fout[2],fdist; fdist=Closest(fa,fb,fc,fd,fe,ff,fout);

	if(fdist<closest->dist2)
		{
		P_float fs=fout[0],ft=fout[1],fu=1-fs-ft;
		closest->side=side;
		c2[0]=p1[0]+fs*p1p2[0]+ft*p1p3[0];c2[1]=p1[1]+fs*p1p2[1]+ft*p1p3[1];c2[2]=p1[2]+fs*p1p2[2]+ft*p1p3[2];
		Transform(c2,closest->c,M);
		closest->dist2=fdist;
		closest->cell=index_cell;
		closest->nb=(fu==0)?(0):(1);closest->nb+=(fs==0)?(0):(1);closest->nb+=(ft==0)?(0):(1);
		free(closest->weights); free(closest->pts); 
		closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb]; 
		int count=0;
		if(fu!=0) {closest->weights[0]=fu; closest->pts[0]=ids[0]; count++;}
		if(fs!=0) {closest->weights[count]=fs; closest->pts[count]=ids[1]; count++;}
		if(ft!=0) {closest->weights[count]=ft; closest->pts[count]=ids[2]; count++;}
		}
	}
return;
}

void IniPolyData(vtkPolyData* model,vtkIdList** PointCells,P_float* CellNormals,vtkPointLocator* PointsLocator)
{
model->BuildCells(); 
model->Modified();

PointsLocator->SetDataSet(model);
PointsLocator->BuildLocator();

int i,ids[3]; double cref[3]={0,0,0},c[3]={0,0,0}; P_float p1[3],p2[3],p3[3],p1p2[3],p1p3[3],n[3],nm;
for(i=0;i<model->GetNumberOfPoints();i++) {PointCells[i]=vtkIdList::New(); model->GetPointCells(i,PointCells[i]); }	
for(i=0;i<model->GetNumberOfCells();i++) 
	{
	ids[0]=model->GetCell(i)->GetPointIds()->GetId(0);		ids[1]=model->GetCell(i)->GetPointIds()->GetId(1);			ids[2]=model->GetCell(i)->GetPointIds()->GetId(2);
	p1[0]=model->GetPoint(ids[0])[0]; p1[1]=model->GetPoint(ids[0])[1]; p1[2]=model->GetPoint(ids[0])[2];p2[0]=model->GetPoint(ids[1])[0]; p2[1]=model->GetPoint(ids[1])[1]; p2[2]=model->GetPoint(ids[1])[2]; p3[0]=model->GetPoint(ids[2])[0]; p3[1]=model->GetPoint(ids[2])[1]; p3[2]=model->GetPoint(ids[2])[2];
	p1p2[0]=p2[0]-p1[0];	p1p2[1]=p2[1]-p1[1];	p1p2[2]=p2[2]-p1[2];	p1p3[0]=p3[0]-p1[0];	p1p3[1]=p3[1]-p1[1];	p1p3[2]=p3[2]-p1[2];
	crossproduct(n,p1p2,p1p3); 	nm=norm(n); CellNormals[3*i]=n[0]/nm;	CellNormals[3*i+1]=n[1]/nm; 	CellNormals[3*i+2]=n[2]/nm;
	}	
}


bool* SelectCellsWithSpline(vtkPolyData* model,P_float M[16],P_float M_inv[16],vtkIdList** PointCells,P_float* CellNormals,vtkPointLocator* PointsLocator,P_float* spline,int splineres)
{
int nb_points=model->GetNumberOfPoints(),nb_cells=model->GetNumberOfCells();

int i,j,selected,count,seed_cell=0,*flagcelli=new int[nb_cells]; 
bool ok=false,flag,stop=false,*flagcell=new bool[nb_cells],*fc=new bool[nb_cells];

for(i=0;i<nb_cells;i++) {flagcell[i]=false; flagcelli[i]=-1;}

//select contour cells
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1]; 
for(i=1;i<splineres;i++)
       {
		closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10;
        GetClosest(model,M,M_inv,PointCells,CellNormals,PointsLocator,closest,spline+3*i,NULL,NULL);
		if(closest->nb!=0) flagcelli[closest->cell]=1;
		free(closest->pts); free(closest->weights);
	   }
free(closest);
// inflate selection
for(i=0;i<nb_cells;i++) fc[i]=false;
for(i=0;i<nb_cells;i++) 
	if(flagcelli[i]==1)
		{
		selected=model->GetCell(i)->GetPointIds()->GetId(0);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) fc[PointCells[selected]->GetId(j)]=true;
		selected=model->GetCell(i)->GetPointIds()->GetId(1);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) fc[PointCells[selected]->GetId(j)]=true;
		selected=model->GetCell(i)->GetPointIds()->GetId(2);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) fc[PointCells[selected]->GetId(j)]=true;
		}
for(i=0;i<nb_cells;i++) if(fc[i]) flagcelli[i]=1;

// fill interior
while(!ok)
	{
	for(i=0;i<nb_cells;i++) if(flagcelli[i]!=1) flagcelli[i]=-1;
	// select seed cell
	for(i=seed_cell+1;i<nb_cells;i++) if(flagcelli[i]!=1) {seed_cell=i; break;} // seed=first cell not selected
	if(seed_cell==nb_cells-1) {free(flagcelli); free(fc); return NULL;} // all selected
	flagcelli[seed_cell]=0;

    // propagate from seed cell
    while(!stop)
        {
        flag=true;
        for(i=0;i<nb_cells;i++)
            if(flagcelli[i]==0)
                {
				selected=model->GetCell(i)->GetPointIds()->GetId(0);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) if(flagcelli[PointCells[selected]->GetId(j)]==-1) {flagcelli[PointCells[selected]->GetId(j)]=0; flag=false;}
				selected=model->GetCell(i)->GetPointIds()->GetId(1);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) if(flagcelli[PointCells[selected]->GetId(j)]==-1) {flagcelli[PointCells[selected]->GetId(j)]=0; flag=false;}
				selected=model->GetCell(i)->GetPointIds()->GetId(2);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) if(flagcelli[PointCells[selected]->GetId(j)]==-1) {flagcelli[PointCells[selected]->GetId(j)]=0; flag=false;}
                }
        stop=flag;
        }

	// count nb of selected cells
	count=0; for(i=0;i<nb_cells;i++) if(flagcelli[i]!=0) count++; 
	if(count<nb_cells-count)		ok=true;
	}

for(i=0;i<nb_cells;i++) if(flagcelli[i]!=0) flagcell[i]=true;
// deinflate selection
for(i=0;i<nb_cells;i++) fc[i]=false;
for(i=0;i<nb_cells;i++) 
	if(flagcell[i])
		{
		selected=model->GetCell(i)->GetPointIds()->GetId(0);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) if(!flagcell[PointCells[selected]->GetId(j)]) fc[i]=true;
		selected=model->GetCell(i)->GetPointIds()->GetId(1);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) if(!flagcell[PointCells[selected]->GetId(j)]) fc[i]=true;
		selected=model->GetCell(i)->GetPointIds()->GetId(2);	for(j=0;j<PointCells[selected]->GetNumberOfIds();j++) if(!flagcell[PointCells[selected]->GetId(j)]) fc[i]=true;
		}
for(i=0;i<nb_cells;i++) if(fc[i]) flagcell[i]=false;


free(flagcelli); free(fc);
return flagcell;
}


vtkStructuredPoints* RBF(int nbp,P_float* coords,P_float* vals,int dim[2])
{
P_float **A=new P_float*[nbp],*w=new P_float[nbp],p[2];
int i,j,k;
for(i=0;i<nbp;i++)
	{
	A[i]=new P_float[nbp];
	for(j=0;j<nbp;j++) A[i][j]=RBF_kernel(coords+2*i,coords+2*j);
	w[i]=vals[i];
	}

vtkMath::SolveLinearSystem(A,w,nbp);

for(i=0;i<nbp;i++)	free(A[i]);
free(A);

vtkStructuredPoints* im=vtkStructuredPoints::New();
	im->SetOrigin(0,0,0);
	im->SetDimensions(dim[0],dim[1],1);
	im->SetSpacing(1./(double)dim[0],1./(double)dim[1],.1);
	im->SetScalarTypeToDouble();
	im->AllocateScalars();
double* ptr=(double*)im->GetScalarPointer();

for(j=0;j<dim[1];j++)
	for(i=0;i<dim[0];i++)
		{
		*(ptr+i+j*dim[0])=0;
		for(k=0;k<nbp;k++)
			{
			p[0]=(P_float)i/(double)dim[0]; p[1]=(P_float)j/(double)dim[1];
			*(ptr+i+j*dim[0])+=w[k]*RBF_kernel(p,coords+2*k);
			}
		}

free(w);
return im;
}


P_float RBF_kernel(P_float* p1,P_float* p2)
{
//P_float sigma=.1;
P_float d2=(p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]);
//return exp(-d2/(sigma*sigma*2));
return sqrt(d2);
}


void DrawPolyData(vtkStructuredPoints* image,P_float ori[3],P_float op1[3],P_float op2[3],P_float n[3],vtkPolyData* poly,unsigned char color[4])
{
vtkCell* cl;
P_float fp[3],fp2[3],alph,d=ori[0]*n[0]+ori[1]*n[1]+ori[2]*n[2],v,pc[10][3],OP1=norm(op1),OP2=norm(op2);
int *dim=image->GetDimensions();
double *spac=image->GetSpacing();
int cross,x,y,x2,y2,k;
bool sgn;

for(int i=0; i<poly->GetNumberOfCells();i++)
	{
	cross=0; 
	cl=poly->GetCell(i);  
	poly->GetPoint(cl->GetPointId(cl->GetNumberOfPoints()-1),fp2);  sgn=((fp2[0]*n[0]+fp2[1]*n[1]+fp2[2]*n[2])>d)?true:false;
	for(k=0; k<cl->GetNumberOfPoints();k++)
		{
		poly->GetPoint(cl->GetPointId(k),fp);  v=fp[0]*n[0]+fp[1]*n[1]+fp[2]*n[2]; 
		if((sgn && v<d)|| (!sgn && v>d)) 
			{
			alph=fp2[0]*n[0]+fp2[1]*n[1]+fp2[2]*n[2]-v;
			if(alph!=0)	{alph=(d-v)/alph; pc[cross][0]=alph*fp2[0]+(1-alph)*fp[0]; pc[cross][1]=alph*fp2[1]+(1-alph)*fp[1]; pc[cross][2]=alph*fp2[2]+(1-alph)*fp[2];}
			else memcpy(pc[cross],fp,3*sizeof(P_float));
			cross++; 
			}
		sgn=(v>d)?true:false;
		memcpy(fp2,fp,3*sizeof(P_float));
		}
	if(cross>1)
		{
		x=(int)((op1[0]*(pc[0][0]-ori[0])+op1[1]*(pc[0][1]-ori[1])+op1[2]*(pc[0][2]-ori[2]))/(OP1*spac[0]));
		y=(int)((op2[0]*(pc[0][0]-ori[0])+op2[1]*(pc[0][1]-ori[1])+op2[2]*(pc[0][2]-ori[2]))/(OP2*spac[1]));
		x2=(int)((op1[0]*(pc[1][0]-ori[0])+op1[1]*(pc[1][1]-ori[1])+op1[2]*(pc[1][2]-ori[2]))/(OP1*spac[0]));
		y2=(int)((op2[0]*(pc[1][0]-ori[0])+op2[1]*(pc[1][1]-ori[1])+op2[2]*(pc[1][2]-ori[2]))/(OP2*spac[1]));
		if(x>=0 && y>=0 && x<dim[0] && y<dim[1] && x2>=0 && y2>=0 && x2<dim[0] && y2<dim[1]) 
		Line2D(image,x,y,x2,y2,color);
		}
	}
}


void DrawMesh(vtkStructuredPoints* image,P_float ori[3],P_float op1[3],P_float op2[3],P_float n[3],P_float* pts,int nb_cells,int** cells,unsigned char color[4])
{
int* cl;
P_float *fp,*fp2,alph,d=ori[0]*n[0]+ori[1]*n[1]+ori[2]*n[2],v,pc[10][3],OP1=norm(op1),OP2=norm(op2);
int *dim=image->GetDimensions();
double *spac=image->GetSpacing();
int cross,x,y,x2,y2,k;
bool sgn;

for(int i=0; i<nb_cells;i++)
	{
	cross=0; 
	cl=cells[i];
	fp2=pts+3*cl[cl[0]];  sgn=((fp2[0]*n[0]+fp2[1]*n[1]+fp2[2]*n[2])>d)?true:false;
	for(k=0; k<cl[0];k++)
		{
		fp=pts+3*cl[k+1];  v=fp[0]*n[0]+fp[1]*n[1]+fp[2]*n[2]; 
		if((sgn && v<d)|| (!sgn && v>d)) 
			{
			alph=fp2[0]*n[0]+fp2[1]*n[1]+fp2[2]*n[2]-v;
			if(alph!=0)	{alph=(d-v)/alph; pc[cross][0]=alph*fp2[0]+(1-alph)*fp[0]; pc[cross][1]=alph*fp2[1]+(1-alph)*fp[1]; pc[cross][2]=alph*fp2[2]+(1-alph)*fp[2];}
			else memcpy(pc[cross],fp,3*sizeof(P_float));
			cross++; 
			}
		sgn=(v>d)?true:false;
		fp2=fp;
		}
	if(cross>1)
		{
		x=(int)((op1[0]*(pc[0][0]-ori[0])+op1[1]*(pc[0][1]-ori[1])+op1[2]*(pc[0][2]-ori[2]))/(OP1*spac[0]));
		y=(int)((op2[0]*(pc[0][0]-ori[0])+op2[1]*(pc[0][1]-ori[1])+op2[2]*(pc[0][2]-ori[2]))/(OP2*spac[1]));
		x2=(int)((op1[0]*(pc[1][0]-ori[0])+op1[1]*(pc[1][1]-ori[1])+op1[2]*(pc[1][2]-ori[2]))/(OP1*spac[0]));
		y2=(int)((op2[0]*(pc[1][0]-ori[0])+op2[1]*(pc[1][1]-ori[1])+op2[2]*(pc[1][2]-ori[2]))/(OP2*spac[1]));
		if(x>=0 && y>=0 && x<dim[0] && y<dim[1] && x2>=0 && y2>=0 && x2<dim[0] && y2<dim[1]) 
		Line2D(image,x,y,x2,y2,color);
		}
	}
}

void SetScalars(vtkDataSet* mesh,double *uca)
{
if(mesh->GetPointData()->GetScalars()!=NULL) 
	{
	vtkDoubleArray* da=dynamic_cast<vtkDoubleArray*>(mesh->GetPointData()->GetScalars());
	for(int j=0;j<mesh->GetNumberOfPoints();j++) da->SetValue(j,uca[j]);
	delete [] uca;
	}
else 
	{
	vtkDoubleArray* da=vtkDoubleArray::New();
	da->SetArray(uca,mesh->GetNumberOfPoints(),0);
	mesh->GetPointData()->SetScalars(da);
	}
}