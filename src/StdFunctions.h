#pragma once

#include <vtkPolyData.h>
#include <vtkStructuredPoints.h>
#include <sstream>

#include <string> 
#include <vtkSmartPointer.h>
#include <vtkCommand.h>

#define P_float double

const P_float PI=3.14159265359;
const P_float DSQRT2=1./sqrt(2.);
const int NB_LDM_MAX=100;
const int NB_TIME_MAX=50;
const int NB_SMRI_MAX=5;
const int DEFAULTWINDOW=2000;
const int DEFAULTLEVEL=800;
const P_float DMAX_ICP=1000;

typedef struct { int nb; int* pts; P_float* weights; int side; P_float dist2; int cell; int model; P_float c[3];P_float n[3]; P_float t; bool crossing;} CLOSESTSTRUCT;

inline P_float dotproduct(P_float u[3], P_float v[3]);
inline void crossproduct( P_float cp[3], P_float u[3], P_float v[3]);
inline void crossproduct( P_float cp[3], const P_float u[3], const P_float v[3]);
inline P_float tripleproduct( P_float u[3], P_float v[3], P_float w[3]);
inline P_float norm(P_float u[3]);
inline P_float norm(const P_float u[3]);
inline P_float dist3D(P_float p1[3],P_float p2[3]);
inline P_float GetTetrahedronVolume(P_float P0[3],P_float P1[3],P_float P2[3],P_float P3[3]);
P_float GetCircumscribedRadius(P_float P0[3],P_float P1[3],P_float P2[3],P_float P3[3]);
P_float GetInscribedRadius(P_float P0[3],P_float P1[3],P_float P2[3],P_float P3[3]);
P_float GetTriangleQuality(P_float P1[3],P_float P2[3],P_float P3[3]);
void GetCurvatureTriangle(P_float Kg_Kh_S[3],P_float n[3],P_float p[3],int nbneighb,P_float *pn,P_float Ncheck[3]);
void Decompose(P_float vtg[3],P_float vn[3],P_float v[3],P_float n[3]);
void GetDerivatives(P_float extrapolation,P_float f[3],P_float df[6]);
void GetForces_Shape(const P_float *P,const P_float *P1,const P_float *P2,const P_float *P3,const P_float *n,const P_float *pref,const P_float h,P_float F[3],P_float F1[3],P_float F2[3],P_float F3[3],P_float lambda[4]);
void GetDerivativesOfCrossProduct(P_float v2[3],P_float n[3],P_float normcross,P_float df[9]);

void Identity(P_float M[16]);
void Invert_M(P_float M[16],P_float M_inv[16]);
void Invert_M33(P_float M[9],P_float M_inv[9]);
void Invert_MSymmetric(P_float M[6],P_float M_inv[6]);
void Multi_M(P_float M[16],P_float M2[16],P_float M3[16]);
void Multi_MsT(P_float U[3],P_float M[6],P_float T[3]);

void GetPolyDataNormal(vtkPolyData* model,const int indexp,P_float n[3]);
P_float* GetPolyDataNormals(vtkPolyData* model);
P_float* GetPolyDataPoints(vtkPolyData* model);

void Transform(P_float pin[3],P_float pout[3], P_float M[16]);
void Transform_R(P_float pin[3],P_float pout[3], P_float M[16]);
void Transform(vtkPolyData* model,vtkPolyData* model_out, P_float rt[3],P_float tr[3]);
void Transform(vtkPolyData* model,vtkPolyData* model_out, P_float M[16]);
//void Transform(CSimplexMesh* model, P_float M[16]);
void GetExtP(P_float p[3],P_float s[3], P_float extp[18]);
void PowerCrust(vtkPolyData* model);
void PerturbSurf(int nb_points,P_float* pts,P_float* normals,P_float* mass,P_float maxlength);

void QtoM(P_float tr[3],P_float q[4],P_float M[16],bool post_rotation);
void MtoQ(P_float tr[3],P_float q[4],P_float M[16],bool post_rotation);
void QtoAN(P_float q[4],P_float *alpha,P_float *beta,P_float *theta);
void ANtoQ(P_float alpha,P_float beta,P_float theta,P_float q[4]);
void QNorm(P_float q[4]);
void MNorm(P_float M[16]);
void QComplete(P_float qout[4],P_float q[3]);
void QSlerp(P_float q[4],P_float t,P_float qs[4],P_float qe[4]);
void GetConic(P_float M[16],P_float alpha,P_float beta); // alpha: elevation ; beta:revolution
void NtoM(P_float tr[3],P_float rt[3],P_float M[16]);
void MtoN(P_float tr[3],P_float rt[3],P_float M[16]);
void Differences(P_float diff[2],P_float M1[16],P_float M2[16]);

vtkStructuredPoints* GetDistMap_thickness(vtkStructuredPoints* distmap1,vtkStructuredPoints* distmap2,P_float M[16]);
P_float* GetDistMap_pointdistances(vtkPolyData* model1,vtkPolyData* model2,vtkIdList** PointCells2,P_float* CellNormals2,vtkPointLocator* PointsLocator2,P_float M[16]);
P_float* GetDistMap_pointdistances(vtkPolyData* model1,vtkStructuredPoints* distmap2,P_float M[16]);
void GetDistMap_collisions(P_float *distances,P_float *distances_diff,vtkPolyData* model1,vtkIdList** PointCells1,P_float* CellNormals1,vtkPointLocator* PointsLocator1,vtkPolyData* model2,bool *selectP,P_float *distances_ref,P_float M[16]);
void GetDistMap_collisions(P_float *distances,P_float *distances_diff,vtkStructuredPoints* distmap1,vtkPolyData* model2,bool *selectP,P_float *distances_ref,vtkStructuredPoints* Thickness,P_float M[16]);
int* GetDistMap_collisions(P_float force[3],vtkStructuredPoints* distmap1,vtkStructuredPoints* distmap2,P_float dmin,P_float M[16]);
P_float FitSphere(P_float C[3],int nb_pts, P_float *pts,P_float thresh);
P_float FitDoubleSphere(P_float C[3],int nb_pts1, P_float *pts1,int nb_pts2, P_float *pts2,P_float thresh);

void GetPrincipalDirections(P_float pm[3],P_float axis[3][3],int nb_points,P_float *points,bool normalize);	// get OBB centred on center of gravity
void Barycenter(P_float e[3],P_float Pproj[3], P_float P1[3],P_float P2[3], P_float P3[3]); //barycenter of 3 points
void Barycenter(P_float e[4],P_float P[4], P_float P1[3],P_float P2[3], P_float P3[3],P_float P4[3]); //barycenter of 4 points
void Line3D_UC(vtkStructuredPoints* vol,P_float iX0, P_float iY0, P_float iZ0, P_float iX1, P_float iY1, P_float iZ1);
void Line2D(vtkStructuredPoints* im,int iX0, int iY0,int iX1, int iY1,unsigned short* color);
void Line2D(vtkStructuredPoints* im,int iX0, int iY0,int iX1, int iY1,unsigned char* color);
void DrawPolyData(vtkStructuredPoints* image,P_float ori[3],P_float op1[3],P_float op2[3],P_float n[3],vtkPolyData* poly,unsigned char color[4]);
void DrawMesh(vtkStructuredPoints* image,P_float ori[3],P_float op1[3],P_float op2[3],P_float n[3],P_float* pts,int nb_cells,int** cells,unsigned char color[4]);
vtkStructuredPoints* voxelize(P_float* bounds,P_float* spac,vtkPolyData* model);
vtkStructuredPoints* distancemap(vtkStructuredPoints* input, bool converttoint=true);
vtkPolyData* Concatenate(vtkPolyData* model1,vtkPolyData* model2);
void Differences(P_float distdev[2],vtkPolyData* model1,vtkPolyData* model2,bool pointtopoint); // returns dist and std dev
void Pad(vtkStructuredPoints* vol,unsigned short val, P_float dist);
bool RotateVolume(vtkStructuredPoints* vol,int x,int y);

void Reject_Histo(int nb, double *data, P_float reject_percentage);

int* Concatenate(int* l1,int* l2,int p1,int p2,bool checkorder=false);
int* Concatenate(int* l1,int* l2);
void Separate(int* l,int p1,int a, int p2,int p3, int b, int p4,int* l1,int* l2);

void IniPolyData(vtkPolyData* model,vtkIdList** PointCells,P_float* CellNormals,vtkPointLocator* PointsLocator);
void SetScalars(vtkDataSet* mesh,double *uca);

void GetBoundingBox(vtkStructuredPoints* vol,P_float *offsettransform,P_float bounds[3][2]);
void GetBoundingBox(vtkStructuredPoints* vol,P_float **Transform,P_float **Spacing,P_float bounds[3][2]);
void GetBoundingBox(vtkPolyData* model,P_float bounds[3][2]);
int GetClosestPoint(vtkPolyData* model,vtkPointLocator* PointsLocator,P_float p[3],bool* selectcells,P_float n[3]);
void GetClosest(vtkPolyData* model,P_float M[16],P_float M_inv[16],vtkIdList** PointCells,P_float* CellNormals,vtkPointLocator* PointsLocator,CLOSESTSTRUCT* closest,P_float p[3],bool* selectcells,P_float n[3]);
bool* SelectCellsWithSpline(vtkPolyData* model,P_float M[16],P_float M_inv[16],vtkIdList** PointCells,P_float* CellNormals,vtkPointLocator* PointsLocator,P_float* spline,int splineres);

P_float* GetGlobalMeanStddev(vtkStructuredPoints* input); // global mean/stddev
P_float* GetGlobalMeanStddev(vtkStructuredPoints* input,int nbtrials);  // global mean/stddev with nbtrials points
P_float* GetLocalMeanStddev(vtkStructuredPoints* input,int nbdim,int windowsize); // local mean/stddev
P_float GetNoiseStddev(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold); // noise stddev 
P_float GetNoiseStddev(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold,int nbtrials); // noise stddev with nbtrials points
P_float GetNoiseStddev(vtkStructuredPoints* input,P_float *MeanStddev_Vol,P_float treshold); // noise stddev using precalculated local mean/stddev
P_float GetSNR(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold); // SNR
P_float GetSNR(vtkStructuredPoints* input,int nbdim,int windowsize,P_float treshold,int nbtrials); // SNR with nbtrials points
void SmoothVol(vtkStructuredPoints* input,P_float stddev,int dim);
void SmoothVolAnisotropic(vtkStructuredPoints* input,P_float difffactor,P_float difftresh,int nbit,int dim);
P_float NLDenoising(vtkStructuredPoints* input,int nbdim,int neighborsize,int maxdist,int nbtrials,P_float h);
inline P_float EuclideanDistance_Vol(unsigned short* ptvol,int* dim,int i,int j,int neighborsize);
inline P_float EuclideanDistance_Im(unsigned short* ptvol,int* dim,int i,int j,int neighborsize);
void ComputeDifference_vol(vtkStructuredPoints* input1,vtkStructuredPoints* input2);

vtkStructuredPoints* MergeMRI(vtkStructuredPoints** vols,int nbMRI,bool fillholes,bool normalize,bool overwrite);
vtkStructuredPoints* GradientMagnitude(vtkStructuredPoints* input,int dim);
P_float GaussianFilter(P_float std_dev2, int nb_val, P_float* val, P_float* dist2);
vtkStructuredPoints* Gradient(vtkStructuredPoints* vol,int nbdim,P_float** Transform_RTMRI=NULL);
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float pos[3],int interpolationmode,P_float step);
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step);
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float *transform_P,P_float *spacing_P,P_float pos[3],int interpolationmode,P_float step);
int Gradient(P_float grad[3],vtkStructuredPoints* vol,P_float *transform_P,P_float *spacing_P,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step);
int Gradient_Int(P_float grad[3],vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step);
int Gradient_Dbl(P_float grad[3],vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode,P_float step);

//gl/gradient interpolation in cartesian volume
P_float Interpolation(vtkStructuredPoints* vol,P_float pos[3],int interpolationmode);
P_float Interpolation(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode);
P_float Interpolation_Int(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode);
P_float Interpolation_Dbl(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode);
void vInterpolation(P_float *v,vtkStructuredPoints* vol,P_float pos[3],int interpolationmode);
void vInterpolation(P_float *v,vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int interpolationmode);
//gl interpolation near planes
P_float Interpolation(vtkStructuredPoints* vol,P_float pos[3],int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode);
P_float Interpolation(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode);
void vInterpolation(P_float *v,vtkStructuredPoints* vol,P_float pos[3],int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode);
void vInterpolation(P_float *v,vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,int pselect,P_float **transform_P,P_float **spacing_P,P_float tol,int interpolationmode);
//gl interpolation in cylindric volume
P_float Interpolation(vtkStructuredPoints* vol,P_float posx,P_float posy,P_float posz,const P_float M[16],const P_float spacing[3],int interpolationmode);
P_float Interpolation(vtkStructuredPoints* vol,P_float pos[3],const P_float M[16],const P_float spacing[3],int interpolationmode);
void GetCylindricCoordinates(P_float cc[3],const P_float p[3],const P_float M[16]);

unsigned short InterpolationSmooth(vtkStructuredPoints* vol,P_float p[3],P_float t1[3],P_float t2[3],P_float **GaussianWeights,int GaussianR,int interpolationmode);
void Normalize(vtkStructuredPoints* vol,const unsigned short Max,const unsigned short Min);
unsigned short GetMax(vtkStructuredPoints* vol);
unsigned short GetMin(vtkStructuredPoints* vol);

int GetSelectP(P_float pos[3],vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol);
int GetSelectP(P_float n2[3],P_float pos[3],P_float n[3],vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol);

void LdmRigidTransform(const char* filein,const char* fileout,P_float M[16]);
void LdmRigidTransform(const char* filein,vtkPoints* ldm_out,P_float M[16]);
void LdmRigidTransform(vtkPoints* ldm_in,vtkPoints* ldm_out,P_float M[16]);

void FlipCSystem_X(P_float M[16]);
void FlipModel_X(vtkPolyData* model);
void FlipVolume_X(vtkStructuredPoints* vol);

P_float Closest(P_float fa,P_float fb,P_float fc,P_float fd,P_float fe,P_float ff,P_float fout[2]);
P_float GetEquidistant2(P_float F[4],P_float Fp[4],P_float P[3],P_float P1[3],P_float P2[3],P_float d12,P_float d22);
P_float GetEquidistant3(P_float F[4],P_float P[3],P_float P1[3],P_float P2[3],P_float P3[3]);

void Regularize_Rigid(P_float T[16],int nb_points,P_float *points,P_float* f,bool* ign);
void Regularize_Simi(P_float T[16],int nb_points,P_float *points,P_float* f,bool* ign);
void Regularize_Affine(P_float T[16],int nb_points,P_float *points,P_float* f,bool* ign);
void Regularize_Centered(P_float T[16],int nb_points,P_float *points,P_float* f,P_float RC[3],bool TR,bool* ign);
void Regularize_Translation(P_float T[16],int nbpoints,P_float* f,bool* ign);

void AnalyseEnergy(P_float ACC_DO_NOM_RON_CR[5],P_float *E,int nbpoints,int depth,P_float s,P_float thresh);
void GetGradEnergies(P_float* E,int nb_points,P_float* points,P_float* normals,vtkStructuredPoints* MRI_grad,P_float step,int depth,int interpolationmode,P_float* ROI,bool oppositedirection=false);
P_float* GetGradEnergies(P_float* E,int nb_points,P_float* points,P_float* normals,vtkStructuredPoints* MRI_grad,P_float **transform_P,P_float **spacing_P,P_float tol,P_float step,int depth,int interpolationmode,P_float* ROI,bool oppositedirection);
P_float GetGradEnergy(int nb_points,P_float* points,P_float* normals,vtkStructuredPoints* MRI_grad, int interpolationmode,P_float *OffsetTransform);
void GetIPEnergies(P_float* E,int metric,int nbpoints,int **IP,int **IPref,int mn,int pn,int mnref,int pnref);
void GetIPEnergies(P_float* E,int metric,int nbpoints,P_float **IP,P_float **IPref,int mn,int pn,int mnref,int pnref);
P_float GetIPSimilarity(int metric,int *IP,int *IPref,int mn,int pn,int mnref,int pnref,int offset);
P_float GetIPSimilarity(int metric,P_float *IP,P_float *IPref,int mn,int pn,int mnref,int pnref,int offset);
P_float GetIPSimilarity(int metric,int nbpoints,int **IP,int **IPref,int mn,int pn);
void ProcessE(P_float* E,int nb_points,int depth);
void GetIPForces(P_float* f,bool *ign,P_float *E,int nbpoints,int depth,P_float s,P_float *normals,P_float alpha,P_float thresh=1E10);
void GetIPForces(P_float* f,bool *ign,P_float *E,int nbpoints,int mn,int mnref,P_float s,P_float *normals,P_float alpha,P_float tresh=1E10);
void GetDemonForces(P_float* f,bool *ign,int nbpoints,P_float* normals,int **IP,int **IPref,P_float step,int mn,int pn,int mnref,int pnref,P_float alpha);

//gl/gradient interpolation in cartesian volume
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,int** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI,P_float *OffsetTransform);
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI,P_float *OffsetTransform);
	//+ smoothing
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,P_float *tangents1,P_float *tangents2,vtkStructuredPoints* vol,int** IP,int mn,int pn,P_float s,P_float** IP_GaussianWeights,int IP_GaussianR,int interpolationmode,P_float *ROI,P_float *OffsetTransform);

//gl interpolation near planes
P_float* ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol,int** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI);
P_float* ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float **transform_P,P_float **spacing_P,P_float tol,P_float** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI);

//gl/gradient interpolation in cylindric volume
void ComputeIP(int nbpoints,P_float* pts,P_float* normals,vtkStructuredPoints* vol,P_float *transform_P,P_float *spacing_P,int** IP,int mn,int pn,P_float s,int interpolationmode,P_float *ROI);

void ProcessIP(int nb_points,int** IP,int mn,int pn,int* IPmask,int IP_masksize);
void ProcessIP(int nb_points,P_float** IP,int mn,int pn,int* IPmask,int IP_masksize);
void ProcessIP(int** IPproc,int nb_points,int** IP,int mn,int pn,int* IPmask,int IP_masksize);
void ProcessIP(P_float** IPproc,int nb_points,P_float** IP,int mn,int pn,int* IPmask,int IP_masksize);

// thread-based method functions
inline P_float Get3x3Determinant(const P_float M[3][3]);
P_float Get3x3Determinant(const P_float M[4][4],const int x,const int y);
P_float Get4x4Determinant(const P_float M[4][4]);
void GetThreadCoefficientsFromLocalCoordinates(P_float R[12],const P_float P0[3],const P_float P1[3],const P_float P2[3],const P_float P3[3]);
void GetThreadCoefficients(P_float R[12],const P_float Pproj[3],const P_float P1[3],const P_float P2[3],const P_float P3[3],const P_float n[3],const P_float h,const P_float vol);
void GetThreadVectors(P_float U[3],P_float V[3],P_float W[3],const P_float R[12],const P_float P0[3],const P_float P1[3],const P_float P2[3],const P_float P3[3]);
P_float GetModePatternUU(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12]);
P_float GetModePatternVV(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12]);
P_float GetModePatternWW(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12]);
P_float GetModePatternUV(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12]);
P_float GetModePatternUW(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12]);
P_float GetModePatternVW(P_float H[4][3],const P_float U[3],const P_float V[3],const P_float W[3],const P_float R[12]);

P_float GetStress(P_float elongation,P_float a);
P_float GetStressDerivative(P_float elongation,P_float a);

vtkStructuredPoints* RBF(int nbp,P_float* coords,P_float* vals,int dim[2]);
P_float RBF_kernel(P_float* p1,P_float* p2);

class ErrorObserver : public vtkCommand
{
public:
  ErrorObserver():
    Error(false),
    Warning(false),
    ErrorMessage(""),
    WarningMessage("") {}
  static ErrorObserver *New()
  {
  return new ErrorObserver;
  }
  bool GetError() const
  {
  return this->Error;
  }
  bool GetWarning() const
  {
  return this->Warning;
  }
  void Clear()
  {
  this->Error = false;
  this->Warning = false;
  this->ErrorMessage = "";
  this->WarningMessage = "";
  }
  virtual void Execute(vtkObject *vtkNotUsed(caller),
                       unsigned long event,
                       void *calldata)
  {
  switch(event)
    {
    case vtkCommand::ErrorEvent:
      ErrorMessage = static_cast<char *>(calldata);
      this->Error = true;
      break;
    case vtkCommand::WarningEvent:
      WarningMessage = static_cast<char *>(calldata);
      this->Warning = true;
      break;
    }
  }
  std::string GetErrorMessage()
  {
  return ErrorMessage;
  }
std::string GetWarningMessage()
  {
  return WarningMessage;
  }
private:
  bool        Error;
  bool        Warning;
  std::string ErrorMessage;
  std::string WarningMessage;
};