#pragma once
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointLocator.h>
#include <vtkStructuredPoints.h>
#include <vtkPoints.h> 
#include "StdFunctions.h"
#include "SimplexMesh.h"


const int NB_NEIGHBORHOOD=2;
const P_float INTERNALFORCEEXTRAPOLATION=0; // fully isotropic
const P_float ISOTROPICEXTRAPOLATION=0.7; // isotropic
const P_float ANISOTROPICEXTRAPOLATION=0.9; // anisotropic

class CSimplexMesh
{

public:
    CSimplexMesh(void);
    ~CSimplexMesh(void);

    virtual void Free();
    virtual void Copy(CSimplexMesh* mesh);

    virtual void LoadMesh(const char* filename);
    virtual void SaveMesh(const char* filename);
	void SaveVTKMesh(const char* filename,int lift);
	void SaveSTLMesh(const char* filename,int lift);

    int GetNumberOfPoints();
    int GetNumberOfCells();

	void GetBoundingBox(P_float bounds[3][2]);
	void GetBoundingDOP(P_float bounds[9][2]);
	void GetPrincipalAxis(P_float pm[3],P_float axis[3][3]);
	void GetPrincipalAxis(const P_float p1[3],const P_float p2[3],const P_float pm[3],P_float axis[3][3]);
	void GetPrincipalAxis(const P_float p1[3],const P_float p2[3],const P_float p3[3],const P_float pm[3],P_float axis[3][3]);
	void RigidTransform(P_float M[16]);

	virtual CSimplexMesh* IncreaseResolution();
    virtual CSimplexMesh* DecreaseResolution();
    virtual void UpdatePointsFromLowerRes(CSimplexMesh* mesh);
    virtual void UpdatePointsFromHigherRes(CSimplexMesh* mesh);

    void SetPoint(int index,P_float x,P_float y, P_float z);    void SetPoint(int index,P_float p[3]);  void SetPoint_tm1(int index,P_float x,P_float y, P_float z);    void SetPoint_tm1(int index,P_float p[3]);
	void GetPoint(int index,P_float p[3]); P_float* GetPoints(); P_float* GetPoint(int index);  void GetPoint_tm1(int index,P_float p[3]); P_float* GetPoints_tm1(); P_float* GetPoint_tm1(int index);
	void SetSpeed(int index,P_float x,P_float y, P_float z);    void SetSpeed(int index,P_float p[3]); void SetSpeed_tm1(int index,P_float x,P_float y, P_float z);    void SetSpeed_tm1(int index,P_float p[3]);
    void GetSpeed(int index,P_float p[3]);  P_float* GetSpeeds(); P_float* GetSpeed(int index);  void GetSpeed_tm1(int index,P_float p[3]);  P_float* GetSpeeds_tm1(); P_float* GetSpeed_tm1(int index);
	void GetForce(int index,P_float f[3]);  P_float* GetForces(); P_float* GetForce(int index); void SetForce(int index,P_float x,P_float y, P_float z);    void SetForce(int index,P_float f[3]);  void AddForce(int index,P_float f[3]);
    P_float* GetDfs();  P_float* GetDfs(int index); void GetDfs(int index,P_float dfs[6]); void AddDfs(int index,P_float dfs[6]);
    P_float* GetDfp();  P_float* GetDfp(int index); void GetDfp(int index,P_float dfp[6]); void AddDfp(int index,P_float dfp[6]);

	void SetVal_p(const int index,const P_float v[3]);
	void SetVal_c(const int index,const P_float v[3]);
	P_float* GetVal_p(const int index);
	P_float* GetVal_c(const int index);
	virtual void UpdateVal_c();
	void FillVal0();
	virtual void UpdateValsFromHigherRes(CSimplexMesh* mesh);
	virtual void UpdateValsFromLowerRes(CSimplexMesh* mesh);
	virtual void UpdateForcesFromLowerRes(CSimplexMesh* mesh,bool updatederivatives);

    P_float* GetMass_inv(int index);        P_float* GetMass(int index);     P_float* GetMass_constraint(int index);
    void SetMass(P_float m);    void SetMass();    void SetMass(int index,P_float m);    void SetMass(int index,P_float mx,P_float my,P_float mz);    void SetMass(int index);
    void SetMass_planeconstraint(int index,P_float n[3]);    void SetMass_vectorconstraint(int index,P_float v[3]);
    bool GetMass_Modified();
    virtual void UpdateMass();
	virtual void SwicthOffAllForces();

    void SetCell(int index,int* f); void SetCell(int index,int nb_p, int* p); int GetCell(int index,int* f); int* GetCell(int index); int** GetCells();	virtual void InsertNextCell(int* f);    virtual void DeleteCell(int index);

    virtual void UpdateAll();
    
	P_float* GetCellCenters();    P_float* GetCellCenter(int cell_index); void GetCellCenter(int cell_index,P_float c[3]); void SetCellCenter(int cell_index,P_float c[3]); void UpdateCellCenters();  

    virtual void SaveParams(const char* filename); virtual void LoadParams(const char* filename);
    virtual void UpdateParams(); virtual void UpdateParams(int index); virtual P_float* GetParams(int index); 

    vtkUnstructuredGrid* GetUnstructuredGrid();
    void UpdateUnstructuredGridPoints(); virtual void UpdateUnstructuredGrid(); virtual void UpdateUnstructuredGridCells();
    virtual vtkPolyData* GetPolyData(int lift);
	vtkUnstructuredGrid* CSimplexMesh::GetOBB();

    void SetGamma(P_float gam);    void SetTimeStep(P_float tstep);
    virtual void Equilibrium(); virtual void TimeStep();  
	virtual void UpdateInternalForces();  virtual void UpdateExternalForces(); virtual void UpdateAlphas();
	void SetUseForces(bool useforces);
	void UpdateForcesDamping();
	void UpdatePointsEnergyMinimisation();

    void SetMonitor_edgelenght(bool val,const char* filename);
    void SetMonitor_density(bool val,const char* filename);
    void SetMonitor_displacement(bool val,const char* filename);
   
	void AllocatePointMaterialIndices();
	int* GetPointMaterialIndices(int index);
	void SetPointMaterialIndices(int index, int mat0, int mat1, int scalar);
		
	void AllocateMultiMaterialNeighbors(int matLowerRange, int matUpperRange);
	void SetMultiMaterialNeighbors(int matIndex, int ptIndex, int n0, int n1, int n2);
	int* GetMultiMaterialNeighbors(int matIndex, int ptIndex);

	bool IsMultiMaterialPoint(int index);
	void SetIsMultiMaterialPoint(int index, bool mm);
	void AllocateIsMultiMaterialPoints();

	void SetSplitMesh(bool b);
	bool GetSplitMesh();

	void SetOriginalPointIndex(int pointIndex, int originalIndex); 
	int GetOriginalPointIndex(int pointIndex);
	
	bool GetHasMultiMaterialIndices();
	void SetHasMultiMaterialIndices(bool b);
	
	
private:
	//int* neighbors;   // 1-neighbors
	
	int **multiMaterialNeighbors; //
	bool *isMultiMaterialPoints; // An array to indicate whether a point is multimaterial or not. 

	// For multi-material mesh
	// and then splitting the meshes and deforming each split mesh separately
	// and then combining the split deformed meshes into a whole multi-material mesh.
	// If true, then this CSimplex mesh represents a material mesh of a larger whole multi-material CSimplexMesh.
	bool isSplitMesh; 
	int *originalPointIndices; // Original point indices of the split mesh's points. 
	
	// Indicates whether the materialIndices** array is populated or not. 
	bool hasMultiMaterialIndices;


protected:
	//int **materialIndices; // Size is number of points. Array is of the format [mat0][mat1][scalar]
	int *materialIndices; // Size is number of points. Format is the same as P_float* points array. [P0_mat0][P0_mat1][P0_scalar][P1_mat0][P1_mat1][P1_scalar][P2_mat0][P2_mat1][P2_scalar]...

// variables
	bool USEFORCES;
    int nb_points;
	
    P_float* points; P_float* points_tm1;
    P_float* speeds; P_float* speeds_tm1;
    P_float* params;        
    
    int nb_cells;   int** cells;       P_float* cellcenters;

    P_float* df_p;
    P_float* df_s;
    P_float* forces;

	P_float* Val_p; // used for interpolation across resolutions
	P_float* Val_c; // used for interpolation across resolutions

	bool Flipped; 

	int* neighbors;   // 1-neighbors
	int* neighbors_c; // cell neighbors, i.e. which cells are sharing the ith point. 
	int*** neighbors2;  // n-neighbors 1<=n<=NB_NEIGHBORHOOD. in other words, second order neighbors. Neighbor of neighbors.

	int TimeCounter;
    P_float gamma;
    P_float timestep;
    P_float timestep_inv2;

    P_float* mass_inv;
    P_float* mass;
    P_float** mass_constraint;
    P_float mass_scalar;
    bool mass_modified;

    vtkUnstructuredGrid* UGrid;


// forces
	void Add_indep(int index,P_float *PPref,P_float alpha,P_float extrapolation);
	void AddForce_indep(int index,P_float *PPref,P_float alpha,P_float extrapolation);
	void AddDisp_indep(int index,P_float *PPref,P_float alpha);
	virtual void AddDisp(int index,P_float *F);

// monitor variables
    bool Monitor;
    virtual void Update_Monitor();
    
	bool Monitor_edgelenght; 
	std::string Monitor_edgelenght_filename;  
	virtual void Update_Monitor_edgelenght();
    
	bool Monitor_density; 
	std::string Monitor_density_filename; 
	virtual void Update_Monitor_density();
    
	bool Monitor_displacement; 
	std::string Monitor_displacement_filename; 
	void Update_Monitor_displacement();

	/// Get/Set Neighbors copied from CSimplexSurf
	void SetNeighbors(int index,int i, int j, int k);
	void GetNeighbors(int index,int n[3]);
	int* GetNeighbors(int index);
	void SetNeighbors(int index,int i,int j);
	
	
};
