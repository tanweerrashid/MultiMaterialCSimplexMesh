
#include <float.h>

#include "MMCSimplexSurf.h"

MMCSimplexSurf::MMCSimplexSurf(CSimplexSurf** m, int num) {
	meshes = m;
	numberOfMeshes = num;

	//meshes[0]->UpdatePointsOn2ManifoldEdgeArray();
}

MMCSimplexSurf::~MMCSimplexSurf() {
	for (int i = 0; i < GetNumberOfMeshes(); i++) {
		meshes[i]->Free();
	}

	delete [] meshes;
}

CSimplexSurf* MMCSimplexSurf::GetMesh(int index) {
	if (index >=0 && index <= numberOfMeshes) 
		return meshes[index];
	else {
		std::cerr << "ERROR\n  MMCSimplexSurf::GetMesh() - index out of bounds" << std::endl;
		return NULL;
	}
}

void MMCSimplexSurf::SetMesh(int index, CSimplexSurf *m) {
	if (index >=0 && index <= numberOfMeshes) {
		meshes[index]->Free();
		meshes[index] = m;
	}
	else {
		std::cerr << "ERROR\n  MMCSimplexSurf::SetMesh() - index out of bounds" << std::endl;
	}
}

int MMCSimplexSurf::GetNumberOfMeshes() {
	return numberOfMeshes;
}

void MMCSimplexSurf::SetNumberOfMeshes(int num) {
	numberOfMeshes = num;
}

CSimplexSurf** MMCSimplexSurf::GetMeshes() {
	return meshes;
}

void MMCSimplexSurf::WriteAsVTKPolyData(std::string filename, int index) {
	if (index >= 0 && index <= numberOfMeshes) {
		meshes[index]->writeCSimplexMeshAsVTKPolyData(filename);
	}
	else {
		std::cerr << "ERROR\n  MMCSimplexSurf::WriteAsVTKPolyData() - index out of bounds" << std::endl;
	}
}

//void MMCSimplexSurf::UpdateMeshArray() {
//
//}
//
//bool MMCSimplexSurf::CheckMeshIntegrity() {
//
//}

//bool MMCSimplexSurf::CheckMeshIntegrity(int meshIndex) {
	//// For the sub-meshes, which should all be pure 2-simplex meshes. 
	//if ((meshIndex != 0) && (meshIndex > 0)) {
	//	bool ret = true;

	//	vtkSmartPointer<vtkPolyData> pdata = meshes[meshIndex]->getAsVTKPolyData();
	//	for (int i = 0; i < pdata->GetNumberOfPoints(); i++) {
	//		vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
	//		pdata->GetPointCells(i, cellList);
	//		if (cellList->GetNumberOfIds() != 3) {
	//			ret = false;
	//			break;
	//		}
	//	}
	//	return ret;
	//}
	//// For the first mesh in the array, i.e. the whole multi-material mesh
	//else {
	//	
	//}
//}


// A multi-material 2-Simplex mesh will have some points having more than 3 neighboring points.
// These points are always found on the non-manifold edge of the shared boundary. 
// This function checks if any of the points in a cell are on the non-manifold edge. 
bool MMCSimplexSurf::HasPointsOnNonManifoldEdge(vtkSmartPointer<vtkIdList> pointList, vtkSmartPointer<vtkPolyData> pdata) {
	bool ret = false;
	for (int i = 0; i < pointList->GetNumberOfIds(); i++) {
		if (IsPointOn2ManifoldEdge(pointList->GetId(i), pdata) != true) {
			ret = true;
			break;
		}
	}
	return ret;
}

// In a closed, watertight 2-simplex mesh, a point will always be shared by exactly 3 cells. 
// In a multi-material 2-simplex mesh, there will be non-manifold edges, where points of the
// non-manifold edge will be shared by more than 3 points. 
bool MMCSimplexSurf::IsPointOn2ManifoldEdge(int pointIndex, vtkSmartPointer<vtkPolyData> pdata) {
	vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
	pdata->GetPointCells(pointIndex, cellList);
	
	if (cellList->GetNumberOfIds() == 3) {
		return true;
	}
	else {
		return false;
	}
}

bool MMCSimplexSurf::ApplyTO1_Test(int p1, int p2, int meshIndex) {
	/*bool ret = false;
	// Ensure that the two points are on 2 manifold edges. 
	if (meshes[meshIndex]->IsPointOn2ManifoldEdge(p1) == true && meshes[meshIndex]->IsPointOn2ManifoldEdge(p2) == true) {
		ret = meshes[meshIndex]->TO1_test(p1, p2);
	}
	return ret;*/
}
//
//bool MMCSimplexSurf::ApplyTO1_Test(int p1, int p2, int meshIndex) {
//	bool res = false;
//	res = meshes[meshIndex]->TO1_test(p1, p2);
//	return res;
//}

void MMCSimplexSurf::ApplyTO1(int p1, int p2, int meshIndex) {	
	meshes[meshIndex]->TO1(p1, p2);	
}

void MMCSimplexSurf::DecimateByRemovingTriangularCells(int meshIndex) {
	int *cell;
	
	std::cout << "Starting decimation by removing triangular cells..." << std::endl;

	int counter = 0;

	int i = 0;
	while (i < meshes[meshIndex]->GetNumberOfCells()) {
		cell = meshes[meshIndex]->GetCell(i);
		
		if (cell[0] == 3) {			
			if (ApplyTO1_Test(cell[1], cell[2], meshIndex) == true) {
				ApplyTO1(cell[1], cell[2], meshIndex);
				counter++;
			}
			else if (ApplyTO1_Test(cell[2], cell[3], meshIndex) == true) {
				ApplyTO1(cell[2], cell[3], meshIndex);
				counter++;
			}
			else if (ApplyTO1_Test(cell[3], cell[1], meshIndex) == true) {
				ApplyTO1(cell[3], cell[1], meshIndex);
				counter++;
			}
			else {
				std::cout << "  WARNING: TO1() operator cannot be safely applied to any points in triangular cell " << i << std::endl;	
			}			
		}
		i = i + 1;
	}

	//meshes[meshIndex]->writeCSimplexMeshAsVTKPolyData("data\\triTemp.vtk");
	std::cout << "  Decimation done...\n  Removed " << counter << " cells \n  Running UpdateAll()" << std::endl;
	meshes[meshIndex]->UpdateAll();
}

void MMCSimplexSurf::DecimateByMergingSmallSurfaceAreaCells(int meshIndex, bool basedOnSmallSurfaceAreaCells) {	
	vtkSmartPointer<vtkPolyData> pdata = meshes[meshIndex]->getAsVTKPolyData();

	int selectedCell = -1;
	if (basedOnSmallSurfaceAreaCells == true) {
		// First locate the cell that has the smallest surface area. 
		double surfaceArea = DBL_MAX;
		
		for (int i = 0; i < meshes[meshIndex]->GetNumberOfCells(); i++) {
			if ((meshes[meshIndex]->getSurfaces_cell(i) < surfaceArea) && (HasPointsOnNonManifoldEdge(pdata->GetCell(i)->GetPointIds(), pdata) == false)) {
				surfaceArea = meshes[meshIndex]->getSurfaces_cell(i);
				selectedCell = i;
			}
		}
	}
	else {
		// Randomly select a cell
		int max = meshes[meshIndex]->GetNumberOfCells() - 1;
		int min = 0;
		selectedCell = rand() % (max - min + 1) + min;
	}


	// Next identify all the cells that are neighbors to the current cell using the current cell's points,
	// and store all the neighbors in a list
	vtkSmartPointer<vtkIdList> neighborsList = vtkSmartPointer<vtkIdList>::New();

	std::cout << "Selected cell is : " << selectedCell << std::endl;
	if (selectedCell >= 0) { // Make sure this cell has a valid index. 
		vtkSmartPointer<vtkIdList> cellPoints = pdata->GetCell(selectedCell)->GetPointIds();

		for (int ptid = 0; ptid < cellPoints->GetNumberOfIds(); ptid++) {
			vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New(); 
			pdata->GetCellEdgeNeighbors(selectedCell, cellPoints->GetId(ptid), cellPoints->GetId((ptid + 1) % cellPoints->GetNumberOfIds()), neighbors);

			// Insert neighbor into neighborsList, using InsertUniqueId() to prevent duplicate entries. 
			for (int k = 0; k < neighbors->GetNumberOfIds(); k++) {
				neighborsList->InsertUniqueId(neighbors->GetId(k));
			}
		}
	}

	std::cout << "  Cell " << selectedCell << " has neighboring cells: ";
	for (int j = 0; j < neighborsList->GetNumberOfIds(); j++) {
		std::cout << neighborsList->GetId(j) << ", ";
	}
	std::cout << std::endl;




	bool done = false;
	int cellCount = 0;	
	while (done != true) {
		//std::cout << "   Trying out TO1() with neighboring cell " << neighborsList->GetId(cellCount) << std::endl;

		vtkSmartPointer<vtkGenericCell> cell1 = vtkSmartPointer<vtkGenericCell>::New();
		vtkSmartPointer<vtkGenericCell> cell2 = vtkSmartPointer<vtkGenericCell>::New();
		pdata->GetCell(selectedCell, cell1);
		pdata->GetCell(neighborsList->GetId(cellCount), cell2);
		
		vtkSmartPointer<vtkIdList> cell1pts = cell1->GetPointIds();
		vtkSmartPointer<vtkIdList> cell2pts = cell2->GetPointIds();


		// Locate the two points that make up the edge shared by the two cells. 
		int edge[2] = {-1, -1};
		int edgeCount = 0;
		for (int ci1 = 0; ci1 < cell1pts->GetNumberOfIds(); ci1++) {
			for (int ci2 = 0; ci2 < cell2pts->GetNumberOfIds(); ci2++) {
				//std::cout << "   " << cell1pts->GetId(ci1) << "   " << cell2pts->GetId(ci2) << std::endl;
				if (cell1pts->GetId(ci1) == cell2pts->GetId(ci2)) {
					if (edgeCount < 2) {
						edge[edgeCount] = cell1pts->GetId(ci1);
						edgeCount = edgeCount + 1;
						break;
					}
				}
			}
		}

		// Make sure that the two point ids are valid. 
		// Then apply TO1() test. If test is true, then apply TO1() operator. 
		// If test is false, then move on to the next cell in neighborsList. 
		if (edge[0] >= 0 && edge[1] >= 0) {
			//std::cout << "   Edge selected is made of points: [" << edge[0] << ", " << edge[1] << "]" << std::endl;

			int pt0 = edge[0]; 
			int pt1 = edge[1];

			if (ApplyTO1_Test(pt0, pt1, 0) == true) {
				ApplyTO1(pt0, pt1, 0);
				done = true;

				//std::cout << "   TO1() successful with edge [" << edge[0] << ", " << edge[1] << "]" << std::endl;
				meshes[meshIndex]->UpdateCellCenters();
				meshes[meshIndex]->ComputeVolume();
			}
			else {
				//std::cout << "   Edge [" << edge[0] << ", " << edge[1] << "] is not suitable for TO1() operation" << std::endl;
				cellCount = cellCount + 1;
			}
		}
		
		// If TO1() cannot be safely applied with any neighboring cells, then terminate loop. 
		if (cellCount == (neighborsList->GetNumberOfIds())) {
			done = true;
			std::cout << "  WARNING: TO1() operator cannot be safely applied with any neighboring cells for cell " << selectedCell << std::endl;
		}
	}
	//meshes[meshIndex]->writeCSimplexMeshAsVTKPolyData("data/mesh[0]_AfterSurfaceDecimation.vtk");
}

/**
 * Check to see if all the points of a cell are of the same material.
 * Uses point scalars to make comparison. 
 *
 * Returns true of all points belong to the same material. False otherwise.
 */
bool MMCSimplexSurf::CellPointsAreSameMaterial(int *simplexCell, int meshIndex) {
	bool ret = true;

	int scalar0 = meshes[meshIndex]->GetPointMaterialIndices(simplexCell[1])[2];

	//std::cout << "Point " << simplexCell[1] << " has material scalar: " << scalar0 << std::endl;

	for (int i = 2; i <= simplexCell[0]; i++) {
		int newScalar = meshes[meshIndex]->GetPointMaterialIndices(simplexCell[i])[2];
		
		//std::cout << "Point " << simplexCell[i] << " has material scalar: " << newScalar << std::endl;

		if (newScalar != scalar0) {

			//std::cout << "\nNOT SAME MATERIAL\n" << std::endl;

			ret = false;
			break;
		}
	}
	return ret;
}

void MMCSimplexSurf::DecimateBySubMesh(int matIndex0, int matIndex1, int meshIndex) {
		
	meshes[meshIndex]->ComputeVolume();
	meshes[meshIndex]->UpdateNeighbors();

	// First locate the cell that has the smallest surface area. 
	double surfaceArea = DBL_MAX;
	int selectedCell = -1;
	
	for (int i = 0; i < meshes[meshIndex]->GetNumberOfCells(); i++) {
		int *icell = meshes[meshIndex]->GetCell(i);

		if ((meshes[meshIndex]->getSurfaces_cell(i) < surfaceArea) && (CellPointsAreSameMaterial(icell, meshIndex) == true)) {
			surfaceArea = meshes[meshIndex]->getSurfaces_cell(i);
			selectedCell = i;
		}
	}
	
	std::cout << "  selectedCell is: " << selectedCell << std::endl;

	// Identify the cells that are neighbors to selectedCell
	// Also identify the neighbor with the smallest surface area, and the points that make up the shared edge
	if (selectedCell >= 0) {
		int *scell = meshes[meshIndex]->GetCell(selectedCell);
		
		double area = DBL_MAX;
		int p0 = -1, p1 = -1;

		for (int i = 1; i <= scell[0]; i++) {
			
			int pt0 = scell[i];
			int pt1 = scell[(i % scell[0]) + 1];

			//std::cout << "   pt0: " << pt0 << "\n   pt1: " << pt1 << std::endl << std::endl;

			int neighbors[2][2];
			meshes[meshIndex]->SearchCell(pt0, pt1, neighbors);

			// In neighbors[][], the elements [1][0] and [1][1] will contain the cells shared by pt0 and pt1. According to Gilles' code
			if (neighbors[1][0] != selectedCell) {
				if (meshes[meshIndex]->getSurfaces_cell(neighbors[1][0]) < area) {
					area = meshes[meshIndex]->getSurfaces_cell(neighbors[1][0]);
					p0 = pt0;
					p1 = pt1;
				}
			}
			else {
				if (meshes[meshIndex]->getSurfaces_cell(neighbors[1][1]) < area) {
					area = meshes[meshIndex]->getSurfaces_cell(neighbors[1][1]);
					p0 = pt0;
					p1 = pt1;
				}
			}
		}

		std::cout << "  Selected edge is: " << p0 << ", " << p1 << std::endl;

		if (p0 >=0 && p1 >= 0) {
			if (meshes[meshIndex]->TO1_test(p0, p1) == true) {
				std::cout << "  Applying TO()" << std::endl;
				meshes[meshIndex]->TO1(p0, p1);

				//meshes[meshIndex]->ComputeVolume();
				//meshes[meshIndex]->writeCSimplexMeshAsVTKPolyData("data/mesh[0]_tempDecimation.vtk");
			}
			else {
				std::cout << "  WARNING: Cannot apply TO() on edge " << p0 << " and " << p1 << "for cell " << selectedCell << std::endl;
			}
		}
	}
	
}