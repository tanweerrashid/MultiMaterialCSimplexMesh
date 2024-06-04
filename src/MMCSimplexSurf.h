
#include <iostream>
#include <string>
#include <cstdlib>

#include "SimplexMesh.h"
#include "SimplexSurf.h"

#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkGenericCell.h>
#include <vtkPolyData.h>

class MMCSimplexSurf {

private:
	CSimplexSurf **meshes;
	int numberOfMeshes;

	bool IsPointOn2ManifoldEdge(int pointIndex, vtkSmartPointer<vtkPolyData> pdata);
	bool HasPointsOnNonManifoldEdge(vtkSmartPointer<vtkIdList> pointList, vtkSmartPointer<vtkPolyData> pdata);

public:
	MMCSimplexSurf(CSimplexSurf** m, int num);
	~MMCSimplexSurf();

	CSimplexSurf* GetMesh(int index);
	void SetMesh(int index, CSimplexSurf *m);

	int GetNumberOfMeshes();
	void SetNumberOfMeshes(int num);

	CSimplexSurf** GetMeshes();

	void UpdateMeshArray();

	void WriteAsVTKPolyData(std::string filename, int index);

	bool CheckMeshIntegrity();

	bool CheckMeshIntegrity(int meshIndex);

	// By default, apply this operator on the multi-material mesh 		
	bool ApplyTO1_Test(int p1, int p2, int meshIndex);
	void ApplyTO1(int p1, int p2, int meshIndex = 0);

	void DecimateByRemovingTriangularCells(int meshIndex = 0);

	void DecimateByMergingSmallSurfaceAreaCells(int meshIndex = 0, bool basedOnSmallSurfaceAreaCells = true);

	bool CellPointsAreSameMaterial(int *simplexCell, int meshIndex);

	void DecimateBySubMesh(int matIndex0, int matIndex1, int meshIndex = 0);
};