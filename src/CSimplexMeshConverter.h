#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkIntArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFieldData.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>

#include <vtkInformation.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "SimplexMesh.h"
#include "SimplexSurf.h"

vtkSmartPointer<vtkPolyData> reverseCellOrder(vtkSmartPointer<vtkPolyData> oldpdata, int cellIndex);

vtkSmartPointer<vtkPolyData> orientSimplexSurfMesh(vtkSmartPointer<vtkPolyData> pdata);

void writeOrVisualizePolyData(vtkSmartPointer<vtkPolyData> simplex, bool write, bool visualize, std::string output);

bool compareFunction (int i, int j);

bool isDuplicateCell(vtkSmartPointer<vtkIdList> iList, vtkSmartPointer<vtkIdList> jList);

bool isDuplicateCell(vtkSmartPointer<vtkCell> iCell, vtkSmartPointer<vtkCell> jCell);

bool isNeighbor(int list1[], int list2[]);

vtkSmartPointer<vtkIdList> orderCellList(vtkSmartPointer<vtkIdList> list, vtkSmartPointer<vtkPolyData> objectData);

vtkSmartPointer<vtkPolyData> convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh(std::string input);

CSimplexSurf** Split_MultiMaterialSimplexMesh_Into_Separate_CSimplexMesh(vtkSmartPointer<vtkPolyData> multimaterialSimplex, std::vector<bool> flipNormals);

vtkSmartPointer<vtkPolyData> LoadMMSimplexMeshFromVTKPolyData(std::string filename);