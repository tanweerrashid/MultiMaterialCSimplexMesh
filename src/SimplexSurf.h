#pragma once

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointLocator.h>
#include <vtkStructuredPoints.h>
#include <vtkPoints.h> 
#include "StdFunctions.h"
#include "SimplexMesh.h"

#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

typedef struct DOPCELL{P_float ext[18]; P_float iext[18]; P_float text[18]; int* cells; int nb_cells; int nb_points; P_float** dp;struct DOPCELL* dopm;struct DOPCELL* dopn; struct DOPCELL* dopo; struct DOPCELL* dopp;} DOPCELLSTRUCT;
typedef struct DOPPOINT{P_float ext[18]; P_float iext[18]; P_float text[18]; int* points; int nb_points; P_float** dp;struct DOPPOINT* dopm; struct DOPPOINT* dopn;struct DOPPOINT* dopo; struct DOPPOINT* dopp;} DOPPOINTSTRUCT;
typedef struct AXIALLINK{int nb; int* pts; P_float* w;  P_float f[3]; P_float Rref; P_float dR; P_float W; P_float dist; bool isin; P_float err; bool side;} AXIALLINKSTRUCT; 

class CSimplexSurf  : public CSimplexMesh 
{

public:
	CSimplexSurf(void);
	~CSimplexSurf(void);

	void Free();
	void Allocate(int nb_p,int nb_c);
	void Copy(CSimplexSurf* mesh);

	void SaveMesh(const char* filename,bool writeneighbors,bool writemergedpoints);
	void SaveOFFMesh(const char* filename,int lift);
	void LoadMesh(const char* filename);
	
	bool GetIsAxis(); 	void SetIsAxis(bool isaxis);
	bool GetIsSkin(); 	void SetIsSkin(bool isskin);


	void SetMRI_GradientMRI(vtkStructuredPoints* mri);
	void AddMRI_GradientMRI(vtkStructuredPoints* mri);
	void FreeMRI(int MRI_index);
	void SetMRI(vtkStructuredPoints* mri,vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]);
	void SetMRI(vtkStructuredPoints* mri,P_float offset[16],vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]); 
	void SetMRI(vtkStructuredPoints* mri,P_float **transform_RTMRI,P_float **spacing_RTMRI,vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]); 
	void SetMRI(vtkStructuredPoints* mri,P_float *transform_RadialMRI,P_float *spacing_RadialMRI,vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]); 
	void AddMRI();
	void AddMRI(vtkStructuredPoints* mri,vtkStructuredPoints* mri_grad,P_float bounds[3][2]); 
	void AddMRI(vtkStructuredPoints* mri,P_float offset[16],vtkStructuredPoints* mri_grad,P_float bounds[3][2]); 
	void AddMRI(vtkStructuredPoints* mri,P_float **transform_RTMRI,P_float **spacing_RTMRI,vtkStructuredPoints* mri_grad,P_float bounds[3][2]); 
	void AddMRI(vtkStructuredPoints* mri,P_float *transform_RadialMRI,P_float *spacing_RadialMRI,vtkStructuredPoints* mri_grad,P_float bounds[3][2]); 
	void SetMRI_ROI(int MRI_index,P_float bounds[3][2]);
	//	void SetMRI2(vtkStructuredPoints* mri,P_float transform[16],bool* flagpt,vtkStructuredPoints* mri_grad=NULL);

	void ComputeIntensityProfile(int mn,int pn,P_float s,bool grad);
		void ComputeIntensityProfileFromHigherRes(int mn,int pn,P_float s,CSimplexSurf* mesh);
		void SimplifyIntensityProfile(unsigned short* levels,int nb_levels);
		unsigned short* ComputeLevelsIntensityProfile(int nb_levels);
		void ComputeIntensityProfile_GaussianWeights(P_float r,P_float stddev,P_float s);
		void ResizeIntensityProfile(int mn,int pn);
		void ResizeIntensityProfileRef(int mn,int pn);
	void ComputeSimilarityImage(int depth,const char* filename);
	void SaveIntensityProfile(const char* filename,bool grad);
	void LoadIntensityProfileRef(const char* filename,bool grad);
	void LoadIntensityProfileRefFromHigherRes(CSimplexSurf* mesh,bool grad);
	void GetIntensityProfile(int indexpt,int* prof);	int* GetIntensityProfile(int indexpt);
	void* GetIntensityProfileRef(int indexpt,bool grad,P_float mnpns[2]); 
	void SetIntensityProfileRef(int mn,int pn,P_float l,int* prof);
	void SetIntensityProfileProcessingMask(int* mask,int masksize);
	void ProcessIntensityProfile();
	void ProcessIntensityProfileRef();
	void RegisterIntensityProfile(int mn,int pn,P_float s,P_float searchdepth, P_float searchstep);
	void Test_IntensityProfile(const char* filename);
	void Analyse_IntensityProfile(const char* filename);

	void SetMultiresMethod(int val);
	void SetHandleborder(bool hb);
	CSimplexSurf* IncreaseResolution();
	CSimplexSurf* DecreaseResolution();
	void UpdatePointsFromLowerRes(CSimplexSurf* mesh);
	void UpdatePointsFromLowerRes_Disp(CSimplexSurf* mesh);
	void UpdatePointsFromHigherRes(CSimplexSurf* mesh);
	void UpdateValsFromHigherRes(CSimplexSurf* mesh);
	void UpdateValsFromLowerRes(CSimplexSurf* mesh);
	void UpdateForcesFromLowerRes(CSimplexSurf* mesh,bool updatederivatives);
	void UpdateVal_c();
	void UpdateMergedPointsFromHigherRes(CSimplexSurf* mesh);

	void SetTopo_Increase_L(int frequ,int nb,bool ignoreborder=false);    
	void SetTopo_Decrease_L(int frequ,int nb,bool ignoreborder=false);    
	void SetTopo_Elongation(int frequ,int nb,bool ignoreborder=false); 
	void SetTopo_Exchange(int frequ,int nb,bool ignoreborder=false);   
	void SetTopo_Lref(P_float lref);
	void SetTopo_Tol(P_float tol);
	bool GetTopo();
	void Subdivide_Merge(P_float p1[3],P_float p2[3]);
	void Subdivide_test(P_float p[3]);

	void InsertNextPoint(P_float p[3]); void DeletePoint(int index);

	void UpdateMass();

	void UpdateNeighbors(); void UpdateNeighbors(int cell_index); void UpdateNeighbors2(int pt_index);
	void UpdateNeighbors_c(int pt_index);
	//void SetNeighbors(int index,int i, int j, int k); void SetNeighbors(int index,int i,int j);
	//int* GetNeighbors(int index); void GetNeighbors(int index,int n[3]); 
	int* GetNeighbors_c(int index); void GetNeighbors_c(int index,int nc[3]);
	
	void SetNeighbors_c(int index, int i, int j, int k);
	
	void OrderNeighbors();
	void UpdateNeighbors_vtkPolyData();
	void UpdateNeighbors_c_vtkPolyData();
	bool IsSharedBoundaryCell(vtkSmartPointer<vtkCell> cell);
	vtkSmartPointer<vtkIdList> GetNeighboringPoints(int index, vtkSmartPointer<vtkPolyData> mesh);

	int* GetBorder();
	void Flip(bool reversepts);	
	P_float GetShortestDist(const int istart,const int iend);
	P_float GetShortestPath(const int istart,const int iend,int* path);
	P_float* GetShortestDistGraph(const int index,bool returnpredecessors,int* predecessors);
	P_float* GetNeighborsDistGraph();
	P_float* GetDistCoordinates();
	void LoadRadiiRef(const char* filename);
	void SaveRadii(const char* filename);
	void SaveDistCoordinates(const char* filename,int dim[2],bool normalizeVal,bool isdR);
	void SaveDistCurve(const char* filename,bool isdR);
	void SetElongationRef(const char* filename);
	P_float GetElongation(int index);
	P_float GetElongationPercentage(int index);

	void DeleteCell(int index); void InsertNextCell(int* f);    
	void SearchCell(int p,int c[3]); int SearchCell(int p1,int p2,int c[2][2]); 	int SearchCell(int p1, int p2,int p3);
	void InsertPointInCell(int c,int p,int p1);
	void DeletePointInCell(int c,int p);
	int TO_DeleteCell(int index);

	void UpdateAll();
	
	void UpdateCellCenters();

	void SaveParams(const char* filename);	void LoadParams(const char* filename);
	void UpdateParams(); void UpdateParams(int index); 	P_float* GetParams(int index); 	void GetParams(int index,P_float p[3]); void SetParams(int index,P_float p[3]);	void SetParams(int index,P_float p1,P_float p2,P_float p3);

	void UpdateNormals();
	void SetNormal(int index,P_float x,P_float y, P_float z); 	void SetNormal(int index,P_float p[3]);
	void GetNormal(int index,P_float p[3]);
	P_float* GetNormal(int index);
	void GetTangentVectors(int index,P_float t1[3],P_float t2[3]); 
	void GetTangentVectors(P_float *t1,P_float *t2); 
	void ComputeVolume();
	P_float GetVolume();
	P_float GetVolume_ref();
	void SetVolume(P_float vol);
	void SetVolume_ref(P_float val,int mode);
	P_float GetDv();
	P_float GetHdiff(int index);
	P_float GetCurvature(int index);

	void SetSurface_ref(P_float val,int mode);

	void UpdateUnstructuredGrid();
	void UpdateUnstructuredGridCells();
	vtkPolyData* GetPolyData(int lift);

	void SetRefModel(vtkPolyData* refmodel,int inimode);
	void InsertInternalConstraintPoint(P_float p[3]);	void InsertInternalConstraintPoint(P_float x, P_float y, P_float z);
	void InsertExternalConstraintPoint(P_float p[3]);	void InsertExternalConstraintPoint(P_float x, P_float y, P_float z);
	void InsertFrontierConstraintPoint(P_float p[3]);	void InsertFrontierConstraintPoint(P_float x, P_float y, P_float z);
	void DeleteConstraintPoint(P_float x, P_float y, P_float z,P_float radius);
	void SetInternalConstraintPoints(vtkPoints* pts);	
	void SetExternalConstraintPoints(vtkPoints* pts);	
	void SetFrontierConstraintPoints(vtkPoints* pts);	

	void RemoveAllInternalConstraintPoints();
	void RemoveAllExternalConstraintPoints();
	void RemoveAllFrontierConstraintPoints();
	void SetInternalConstraintPointsFromRefModel();
	void InsertConstraintPointForce(int index,P_float f[3]);
	void InsertConstraintPointForce(int index,P_float x,P_float y,P_float z);
	vtkPoints* GetInternalConstraintPoints();
	vtkPoints* GetExternalConstraintPoints();
	vtkPoints* GetFrontierConstraintPoints();
	void SaveInternalConstraintPoints(const char* filename); void LoadInternalConstraintPoints(const char* filename);
	void SaveExternalConstraintPoints(const char* filename); void LoadExternalConstraintPoints(const char* filename); 
	void SaveFrontierConstraintPoints(const char* filename); void LoadFrontierConstraintPoints(const char* filename);

	void RemoveAllAttachPoints();
	void AttachBorder();
	void SetAttachPoint(int index,P_float pa[3]); 	void SetAttachPoint(int index,P_float x, P_float y, P_float z);	void SetAttachPoint(int index,P_float pa[3],P_float n[3]);
	void SetAttachedPoint(int index,int val);
	int GetAttachPoint(int index,P_float pa[3]);
	int GetAttachedPoint(int index);
	void SetAttachCell(int index,int val);
	int GetAttachCell(int index);
	void SelectCellsWithSpline(P_float* spline,int splineres,int attachment_index,int influence_dist);  
	int* SelectAttachCellsCentroid(int attachment_index);
	int* UpdateAttachedPoints(int attachment_index,P_float* seed_point=NULL,P_float* seed_point2=NULL); // return contour[contlenght]
	void UpdateAttachedFromHigherRes(CSimplexSurf* mesh);
	void SetSlideAttach(bool val);
	void SetMergedPoints(int* mp);
	int* GetMergedPoints();

	void Equilibrium();
	void TimeStep();
	void UpdateInternalForces();
	void UpdateExternalForces(); 

	void SwicthOffAllForces();
	void SetInternalForceAxialConstraint(bool i,int mode,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceLaplacianFlexion(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceLaplacian(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceSurfEquFlexion(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceSurfEqu(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceRefShape(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceVolPres(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetInternalForceSurfPres(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
//	void SetInternalForceGaussianFilterDisplacement(bool i,P_float rigidity);
	void SetExternalForceICP(bool i,double dmax,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetExternalForceGradientMRI(bool i,double size,int depth,P_float alpha,P_float alpha_lim,P_float alpha_incr,vtkStructuredPoints* gradvol,bool oppositedirection);
	void SetExternalForceDemon(bool i,int depth,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetExternalForceIntensityProfile(bool i,int metric,P_float depth,P_float alpha,P_float alpha_lim,P_float alpha_incr,P_float depth_lim=0,P_float depth_incr=0);
	void SetExternalForceConstraintPoints(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr);
	void SetExternalForceConstraintMesh(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr,P_float tresh=0,CSimplexSurf* mesh=NULL);
	void SetExternalForceRegularization(bool en,bool smooth,int transform,P_float lambda,P_float lambda_lim=0,P_float lambda_incr=0);

	int GetDepth_GradientMRI() {return depth_GradientMRI;}

	void SetAxis_model(CSimplexSurf* model);	CSimplexSurf* GetAxis_model();

	AXIALLINK* GetAxialLinks(); void AllocateAxialLinks();
	bool GetInternalForce_AxialConstraint();
	void SetAxialConstraint_mode(int mode);
	void SetAxialLinks_dR(P_float dr);
	P_float* ComputeAxialLinks_error();
	void SetAxialLinks_smoothR(int nb_it);
	void SetAxialLinks_smoothR(CSimplexSurf* axis2);
	void SetAxialLinks_blockside(int side);
	void UpdateAxialLinks_R();
	void UpdateAxialLinks_Rmean();	
	void UpdateAxialLinks_Ropt(int nb_it,P_float epsilon);
	void UpdateAxialLinks_links();
	void UpdateAxialLinks_forces();
	void UpdateAxialLinks_SetUniformR(P_float R);
	P_float GetAxialLinks_RaxisFromBaryCoord(const int nbpt,const int *pts,const P_float* w);
	int* ExternalForce_ConstraintMesh_pts; 

	void GetCrossing_P(CLOSESTSTRUCT* closest,P_float p[3],P_float s[3],int pointindex);
	void GetCrossing_C(CLOSESTSTRUCT* closest,P_float p[3],P_float s[3],int cellindex);
	void GetCrossing_C_ini(CLOSESTSTRUCT* closest,P_float p[3],int cellindex);

	int GetClosestPoint(P_float p[3]);
	int GetClosestBorderPoint(P_float p[3]);
	int GetClosestAttachPoint(P_float p[3]);
	int GetClosestAttachPoint(P_float p[3],int attachindex);
	void GetClosestBorder2Points(int closests[2],P_float dists2[2],P_float p[3]);
	void GetClosest3Points(int closests[3],P_float dists2[3],P_float p[3]);  //returns the 3 closests points
	void GetClosest(CLOSESTSTRUCT* closest,P_float p[3]);
	void GetClosest_P(CLOSESTSTRUCT* closest,P_float p[3],int pointindex); 
	void GetClosest_C(CLOSESTSTRUCT* closest,P_float p[3],int cellindex); 
	void RegisterToSimplexMesh(const char * filename);

	DOPCELLSTRUCT* GetDOPCellTree();
	void BuildDOPCellTree(int depth,P_float tolerance);
	void UpdateDOPCellTree(P_float tolerance);
	void PrintDOPCellTree(const char* filename);
	void DOPCellIniInflation();

	DOPPOINTSTRUCT* GetDOPPointTree();
	void BuildDOPPointTree(int depth,P_float tolerance);
	void UpdateDOPPointTree(P_float tolerance);
	void DOPPointIniInflation();

    void SetMonitor_surfaces(bool val,const char* filename);
    void SetMonitor_axiallinks_error(bool val,const char* filename);

	void verifyPointsAndParams();
	void writeCSimplexMeshAsVTKPolyData(std::string filename);
	vtkSmartPointer<vtkPolyData> getAsVTKPolyData();
	void verify3PointConnectivity();

	// For decimation //////////////////////////////
	void decimateByRemovingTriangularCells();
	void decimateByMergingSmallSurfaceAreaCells(double factor);
	////////////////////////////////////////////////



	
	void Subdivide_Exchange();
	void Subdivide_Increase_L(); 
	void Subdivide_Decrease_L(); 
	void Subdivide_Elongation();
	void Subdivide_Update();

	void TO1(int p1,int p2); bool TO1_test(int p1,int p2);	
	void TO2(int c,int p1,int p2,int p3,int p4);	bool TO2_test(const int c,const P_float tresh,int pt[4],P_float l[2]);

	void TOech(int p1,int p2);	bool TOech_test(int p1,int p2);	void TOech2(int p);	bool TOech2_test(int p);
	void TOcell_i(int c); bool TOcell_i_test(int c);	void TOcell_d(int c); bool TOcell_d_test(int c);
	void TOvertex_i(int p); bool TOvertex_i_test(int p);	void TOvertex_d(int c,bool first); bool TOvertex_d_test(int c,bool* first);
	void TOmerge(int c1,int c2,int offset);

	P_float getSurface();
	P_float getSurfaces_cell(int index);
	
private:	

	/*int* flagmultires;
	void SetFlagMultires(int index,int val);
	int GetFlagMultires(int index);*/

// variables
	int nb_points_border;
	P_float* normals;
	P_float* normals2;
	P_float* surfaces; // Surface area around each point??
	P_float* surfaces_cell; // Surface area of cells
	P_float* h; 
	P_float GetHmean(int index);

	AXIALLINK* axiallinks;

	void UpdateFlipped();

	int SlideAttach;

	P_float volume;
	P_float surface;
	P_float dv;
	P_float volume_ref;
	P_float surface_ref;

	int TimeCounter;

	int MultiresMethod;

	vtkPoints* InternalConstraintPoints;
	vtkPoints* ExternalConstraintPoints;
	vtkPoints* FrontierConstraintPoints;

	int* attached_point;
	int* attached_cell;

	P_float* externalforces;	bool* externalforces_ign;
	P_float* constraintpointsforces;  

	CSimplexSurf* Axis_model;
	bool AxialLinks_UniformR;
	void UpdateAxialLinks_SmoothR(); int Axial_smoothR; CSimplexSurf* Axial_smoothR_axis2; int Axial_blockside;
	void UpdateAxialLinks_SmoothR(CSimplexSurf* axis2);
	void UpdateAxialLinks_FilterR();
	void UpdateAxialLinks_RegularizeSide();
	
	bool IsAxis;
	bool IsSkin;

	int* MergedPoints;

	vtkPolyData* RefModel;
	vtkIdList** RefModelPointCells;
	P_float* RefModelCellNormals;
	vtkPointLocator* RefPointsLocator;
	P_float RefModelVolume;
	P_float RefModelSurface;
	P_float Threshold_volsurfpres;

// DOP
	DOPCELLSTRUCT *DOPCellTree;
	DOPPOINTSTRUCT *DOPPointTree;
	void DOPCellSubdivision(DOPCELLSTRUCT* dop,int* flag_c,int level,int depth);
	void DOPCellBinaryToQuadtree(DOPCELLSTRUCT* dop);
	void DOPCellBinaryToQuadtree2(DOPCELLSTRUCT* dop);
	void DOPCellFree(DOPCELLSTRUCT* dop);
	void DOP2CellFree(DOPCELLSTRUCT* dop);
	void DOPCellExtUpdate(DOPCELLSTRUCT* dop);
	void DOPCellExtUpdate_Axis(DOPCELLSTRUCT* dop);
	void DOPCellInflation(DOPCELLSTRUCT* dop,P_float tolerance,P_float tolerance_r2);
	void PrintDOPCellTree(FILE* f, DOPCELLSTRUCT* dop,int depth);

	void DOPPointSubdivision(DOPPOINTSTRUCT* dop,int* flag_p,int level,int depth);
	void DOPPointBinaryToQuadtree(DOPPOINTSTRUCT* dop);
	void DOPPointFree(DOPPOINTSTRUCT* dop);
	void DOP2PointFree(DOPPOINTSTRUCT* dop);
	void DOPPointInflation(DOPPOINTSTRUCT* dop,P_float tolerance,P_float tolerance_r2);
	void DOPPointExtUpdate(DOPPOINTSTRUCT* dop);

// IP

	vtkStructuredPoints** MRI; vtkStructuredPoints** MRI_grad; 
	P_float ***MRI_Transform; P_float ***MRI_Spacing; P_float** MRI_Normals; P_float **MRI_OffsetTransform;	P_float **MRI_ROI;
	int* MRI_Type; // 0:cartesian	1:rt	 2:radial
	int MRI_nb;
	vtkStructuredPoints* MRI_GradientMRI; P_float MRI_ROI_GradientMRI[6];
	vtkStructuredPoints* MRI_GradientMRI2; P_float MRI_ROI_GradientMRI2[6];

	void IniIntensityProfile(int mn,int pn,P_float s,bool grad);
	void IniIntensityProfileRef(int mn,int pn,P_float s,bool grad);

	int** IntensityProfileRef; P_float** IntensityProfileRef_grad; 
	int IntensityProfile_mn_ref; int IntensityProfile_pn_ref;  

	int* IntensityProfile_Proc_Mask; int IntensityProfile_Proc_MaskSize; 	
	int** IntensityProfile; P_float** IntensityProfile_grad; int IntensityProfile_mn; int IntensityProfile_pn; P_float IntensityProfile_s;
	int IntensityProfile_GaussianR; P_float IntensityProfile_GaussianS;	P_float **IntensityProfile_GaussianWeights;
	unsigned short GetSmoothIntensityProfile(vtkStructuredPoints* vol,P_float px,P_float py,P_float pz,P_float t1[3],P_float t2[3]);

	int Interpolationmode;

	void Update_Monitor();
    
	bool Monitor_surfaces; 
	std::string Monitor_surfaces_filename; 
	virtual void Update_Monitor_surfaces();
    
	bool Monitor_axiallinks_error; 
	std::string Monitor_axiallinks_error_filename; 
	virtual void Update_Monitor_axiallinks_error();
    
	void Update_Monitor_edgelenght();
    void Update_Monitor_density();

// forces variables
	bool InternalForce_AxialConstraint; int AxialConstraint_mode; P_float Alpha_AxialConstraint; P_float Alpha_AxialConstraint_lim; P_float Alpha_AxialConstraint_incr; 
	bool InternalForce_LaplacianFlexion; P_float Alpha_LaplacianFlexion; P_float Alpha_LaplacianFlexion_lim; P_float Alpha_LaplacianFlexion_incr; 
	bool InternalForce_SurfEquFlexion; P_float Alpha_SurfEquFlexion; P_float Alpha_SurfEquFlexion_lim; P_float Alpha_SurfEquFlexion_incr;
	bool InternalForce_Laplacian; P_float Alpha_Laplacian; P_float Alpha_Laplacian_lim; P_float Alpha_Laplacian_incr; 
	bool InternalForce_SurfEqu; P_float Alpha_SurfEqu; P_float Alpha_SurfEqu_lim; P_float Alpha_SurfEqu_incr;
	bool InternalForce_RefShape; P_float Alpha_RefShape;  P_float Alpha_RefShape_lim;  P_float Alpha_RefShape_incr;
	bool InternalForce_VolPres; P_float Alpha_VolPres; P_float Alpha_VolPres_lim; P_float Alpha_VolPres_incr;
	bool InternalForce_SurfPres; P_float Alpha_SurfPres; P_float Alpha_SurfPres_lim; P_float Alpha_SurfPres_incr;
//	bool InternalForce_GaussianFilterDisplacement;
	bool ExternalForce_ICP; double DMax_ICP; P_float Alpha_ICP; P_float Alpha_ICP_lim; P_float Alpha_ICP_incr; 
	bool ExternalForce_GradientMRI; P_float Size_GradientMRI; int depth_GradientMRI; P_float Alpha_GradientMRI;  P_float Alpha_GradientMRI_lim; bool oppositedirection_GradientMRI; P_float Alpha_GradientMRI_incr; 
	bool ExternalForce_Demon; int depth_Demon; P_float Alpha_Demon;  P_float Alpha_Demon_lim;; P_float Alpha_Demon_incr; 
	bool ExternalForce_IntensityProfile; int Metric_IntensityProfile; int Depth_IntensityProfile; P_float Depth_IntensityProfile_float; P_float Depth_IntensityProfile_incr; P_float Depth_IntensityProfile_lim; int Edgedivide_IntensityProfile; P_float Alpha_IntensityProfile; P_float Alpha_IntensityProfile_lim; P_float Alpha_IntensityProfile_incr; 
	bool ExternalForce_ConstraintPoints; P_float Alpha_ConstraintPoints;  P_float Alpha_ConstraintPoints_lim;  P_float Alpha_ConstraintPoints_incr;
	bool ExternalForce_ConstraintMesh;	P_float ExternalForce_ConstraintMesh_tresh; P_float Alpha_ConstraintMesh;	P_float Alpha_ConstraintMesh_lim;	P_float Alpha_ConstraintMesh_incr;	CSimplexSurf* ExternalForce_ConstraintMesh_mesh;

// forces functions
    void UpdateForcesLaplacian();			    
    void UpdateForcesLaplacianFlexion();	    
    void UpdateForcesSurfEqu();				    
    void UpdateForcesSurfEquFlexion();		    
    void UpdateForcesRefShape();			    

	void UpdateForcesAxialConstraint();			
	void UpdateForcesVolPres();					
	void UpdateForcesSurfPres();				
	void UpdateForcesICP();						void GetClosestRefModel(P_float p[3],P_float n[3],P_float u[3]);
	void UpdateForcesExternal();				
	void UpdateForcesConstraintPoints();		
	void UpdateForcesConstraintMesh();			void Init_ConstraintMeshForces();

	void RegularizeExternalForces();
	void RegularizeEnergies(P_float *E);
	void UpdateIntensityProfileForces(int metric); void UpdateGradientMRIForces(); void UpdateDemonForces();
	bool Regularization_En; bool Regularization_Smooth; bool Regularization_Rigid; 	bool Regularization_Simi; bool Regularization_Affine; 
	P_float Regularization_Lambda; P_float Regularization_Lambda_incr; P_float Regularization_Lambda_lim;

	void Add_multi(int index,P_float *pref,P_float h,P_float alpha,P_float extrapolation);

	void AddDisp_indep(int index,P_float *PPref,P_float alpha);
	void AddDisp_multi(int index,P_float *pref,P_float h,P_float alpha);
	void AddDisp(int index,P_float *F);

	void AddForce_indep(int index,P_float *PPref,P_float alpha,P_float extrapolation);
	void AddForce_multi(int index,P_float *pref,P_float h,P_float alpha,P_float extrapolation);

	bool Forcehandleborder;
	bool IsInternalPointsHandlingMethodCollision;
	bool IsExternalPointsHandlingMethodCollision;
	bool IsFrontierPointsHandlingMethodCollision;

	void UpdateAlphas();

// misc
	void IniRefModel(int mode);
	P_float LiftCellcenters(int nb_it);
	P_float GetTriCurvature_Point(int index,P_float Kg_Kh_S[3],P_float n[3]);
	P_float GetTriCurvature_Cell(int index,P_float Kg_Kh_S[3],P_float n[3],P_float nm[3]);
	void MapCurvatureOnPolyData(vtkPolyData* model,int reject_percentage);
	void MapRadiusOnPolyData(vtkPolyData* model);
	void DOPInclusion(P_float p[3],CLOSESTSTRUCT *crossing2, DOPCELLSTRUCT* dop2);
	P_float* Elongations_ref; int* Elongations_refindex;

// axial links
	void GetAxialLinks_PaxisFromPmodel(int index,P_float p[3]);
	void GetAxialLinks_NaxisFromPmodel(int index,P_float n[3]);

// topological operations
	int Topo; void Update_Topo();
	bool Topo_ignoreborder;
	int Topo_Exchange[2];
	int Topo_Decrease_L[2];
	int Topo_Increase_L[2];
	int Topo_Elongation[2];
	P_float Topo_sref;
	P_float Topo_scellref;
	int Topo_nref;
	P_float Topo_lref;
	P_float Topo_Tol;

	//void Subdivide_Exchange();
	//void Subdivide_Increase_L(); 
	//void Subdivide_Decrease_L(); 
	//void Subdivide_Elongation();
	//void Subdivide_Update();

	/*void TO1(int p1,int p2); bool TO1_test(int p1,int p2);	void TO2(int c,int p1,int p2,int p3,int p4);	bool TO2_test(const int c,const P_float tresh,int pt[4],P_float l[2]);
	void TOech(int p1,int p2);	bool TOech_test(int p1,int p2);	void TOech2(int p);	bool TOech2_test(int p);
	void TOcell_i(int c); bool TOcell_i_test(int c);	void TOcell_d(int c); bool TOcell_d_test(int c);
	void TOvertex_i(int p); bool TOvertex_i_test(int p);	void TOvertex_d(int c,bool first); bool TOvertex_d_test(int c,bool* first);
	void TOmerge(int c1,int c2,int offset);*/

	void Monito_TO(const char* filename);
	int Monitor_TO2; int Monitor_TO2_test;	int Monitor_TO1; int Monitor_TO1_test;
	int Monitor_TOech; int Monitor_TOech_test;	int Monitor_TOech2; int Monitor_TOech2_test;
	int Monitor_TOcell_i; int Monitor_TOcell_i_test;	int Monitor_TOcell_d; int Monitor_TOcell_d_test;
	int Monitor_TOvertex_i; int Monitor_TOvertex_i_test;	int Monitor_TOvertex_d; int Monitor_TOvertex_d_test;
	P_float Monitor_TOduration;

	void TO_test();
//	void TOstar(int p);	bool TOstar_test(int p);	void TOtri(int c);	bool TOtri_test(int c);

// basic functions
	void GetCrossing_T(CLOSESTSTRUCT* closest,P_float p[3],P_float s[3],int pindex1,int pindex2,int cindex,P_float sc[3]);
	void GetCrossing_T2(CLOSESTSTRUCT* closest,P_float sl,P_float n[3],P_float p[3],P_float s[3],int pindex1,int pindex2,int cindex,P_float p1[3],P_float s1[3],P_float p2[3],P_float s2[3],P_float p3[3], P_float s3[3]);
	void decompose(int index,P_float f[3],P_float ftg[3],P_float fn[3]); 
	bool GetCutEdges_PrincipalDirection(const int index_cell,P_float tresh,int pt[4],P_float l[2]);
};
	int qsortintmaxtomin (const void * a, const void * b); 
	int qsortP_floatmaxtomin (const void * a, const void * b); 

