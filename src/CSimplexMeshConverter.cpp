
#include "CSimplexMeshConverter.h"

/** 
 * Takes a polydata along with a cell index as input. Returns a new copy of the polydata
 * with everything same, except for the cell having index cellIndex. This cell's point Ids are
 * reversed. 
*/
vtkSmartPointer<vtkPolyData> reverseCellOrder(vtkSmartPointer<vtkPolyData> oldpdata, int cellIndex) {
	vtkSmartPointer<vtkPolyData> newpdata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pts = oldpdata->GetPoints();

	vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

	
	for (int i = 0; i < oldpdata->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkCell> cell = oldpdata->GetCell(i);

		if (i == cellIndex) {
			vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
			for (int j = cell->GetNumberOfPoints() - 1; j >= 0; j--) {
				list->InsertNextId(cell->GetPointId(j));
			}

			cellarray->InsertNextCell(list);
		}
		else {
			vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
			for (int j = 0; j < cell->GetNumberOfPoints(); j++) {
				list->InsertNextId(cell->GetPointId(j));
			}

			cellarray->InsertNextCell(list);
		}
	}

	newpdata->SetPoints(pts);
	newpdata->SetPolys(cellarray);

	newpdata->BuildCells();
	newpdata->BuildLinks();

	/*vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInput(newpdata);
	writer->SetFileName("E:\\Gilles_MRI_Segmentation\\MRISegmentation\\data\\changedcells.vtk");
	writer->SetFileTypeToASCII();
	writer->Write();*/

	return newpdata;
}

/**
 * Takes a CSimplex mesh in vtkPolyData format as input, and then orients all face in one direction.
 * Returns a CSimplex mesh in vtkPolyData format whose faces are all oriented in either clockwise or counter-clockwise direction. 
 */
vtkSmartPointer<vtkPolyData> orientSimplexSurfMesh(vtkSmartPointer<vtkPolyData> pdata) {
	//vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	//reader->SetFileName("C:\\Users\\Tanweer Rashid\\Desktop\\Object26_part2_ConvertedToSimplex.vtk");
	//reader->Update();

	//vtkSmartPointer<vtkPolyData> pdata = reader->GetOutput();
	pdata->BuildCells();
	pdata->BuildLinks();

	
	std::vector<int> cellList; // List of indices of all cells that need to be processed for orientation determination

	bool *orientedCells = new bool[pdata->GetNumberOfCells()]; // Track which cells have been oriented.
	for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
		orientedCells[i] = false;
	}
	

	// Insert the first cell index to initialize the process
	// All other cells will follow this cell's orientation
	cellList.push_back(0); 
	orientedCells[0] = true;

	while (cellList.size() > 0) {
		// Retrieve the top-most cell number in the list, then remove the index from the list
		// Then retrieve the cell from the simplex polydata
		int icell = cellList.at(0); 

		// Check if the selected cell is oriented. If not yet oriented, then remove the cell and use another cell. 
		if (orientedCells[icell] == true) {
			vtkSmartPointer<vtkCell> ithcell = pdata->GetCell(icell);

			for (int ci = 0; ci < ithcell->GetNumberOfPoints(); ci++) {
				ithcell = pdata->GetCell(icell);

				//cout << "Working on cell " << icell << endl;
				//cout << "Cell " << icell << " has points ";
				//for (int q = 0; q < ithcell->GetNumberOfPoints(); q++) {
				//	cout << ithcell->GetPointId(q) << ", ";
				//}cout << endl;


				// Take two points in the cell that form an edge
				int size =  ithcell->GetNumberOfPoints();
				int pi = ithcell->GetPointId(ci % size);
				int pj = ithcell->GetPointId((ci + 1) % size);

				//cout << "\tUsing points " << pi << " and " << pj << endl;

				// Using the two points, retrieve the cell that shares the edge
				vtkSmartPointer<vtkIdList> neighboringcell = vtkSmartPointer<vtkIdList>::New();
				pdata->GetCellEdgeNeighbors(icell, pi, pj, neighboringcell);

				int jcell = neighboringcell->GetId(0); // The cell Id of the neighboring cell. 

				if (orientedCells[jcell] == false) {
					vtkSmartPointer<vtkCell> jthcell = pdata->GetCell(jcell);

					//cout << "\tPoints " << pi << " and " << pj << " has neighboring cell " << jcell << endl;

					// Check and see if the point ordering of the two cells are the same or not. 
					for (int cj = 0; cj < jthcell->GetNumberOfPoints(); cj++) {
						int currentPt = jthcell->GetPointId(cj % jthcell->GetNumberOfPoints());
				
						if (pi == currentPt) {
							int nextPt = jthcell->GetPointId((cj + 1) % jthcell->GetNumberOfPoints());
							int prevPt;
					
							if (cj == 0) {
								prevPt = jthcell->GetPointId((jthcell->GetNumberOfPoints() - 1));
							}
							else {
								prevPt = jthcell->GetPointId(cj - 1);
							}
					
							// If the jth cell has opposite orientation, then reverse the cell's point ordering
							if (pj == nextPt) {
								cout << "\tCell " << jcell << " has opposite orientation. " << endl;
								//vtkSmartPointer<vtkPolyData> copy = vtkSmartPointer<vtkPolyData>::New();
								//copy->DeepCopy(pdata);
								pdata->ReverseCell(jcell);
								//pdata = reverseCellOrder(copy, jcell);


							}
							// Otherwise do nothing. 
							else if (pj == prevPt) {
								cout << "\tCell " << jcell << " has same orientation. " << endl;
							}

							break;
						}
					}
				
					cellList.insert(cellList.end(), jcell); // Insert the jth cell into cellList. 
					orientedCells[jcell] = true; // Mark the jth cell as having been processed. 
					ithcell = pdata->GetCell(icell);
				}

			}

			cellList.erase(cellList.begin()); // Remove the ith cell from the list. 


		}
		else {
			cellList.erase(cellList.begin()); // Removing cells that are not oriented. 
		}
		
		//cellList.erase(cellList.begin());		
		vtkSmartPointer<vtkCell> ithcell = pdata->GetCell(icell);

		// In situations where the mesh is composed of separate pieces
		// such as the left STN and right STN, one part of the mesh will get
		// oriented, and the other parts may not be in cellList.
		// Therefore, search for unoriented cells, and add then to cellList.
		if (cellList.size() == 0) {
			for (int q = 0; q < pdata->GetNumberOfCells(); q++) {
				if (orientedCells[q] == false) {
					cout << "More unoriented cells found... cell no. " << q << endl;
					cellList.push_back(q);
					break;
				}
			}
		}
	}

	//vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	//writer->SetInput(pdata);
	//writer->SetFileName("C:\\Users\\Tanweer Rashid\\Desktop\\Object26_part2_ConvertedToSimplex_Oriented.vtk");
	//writer->SetFileTypeToASCII();
	//writer->Write();

	return pdata;
}

void writeOrVisualizePolyData(vtkSmartPointer<vtkPolyData> simplex, bool write, bool visualize, std::string output) {
    if (write) {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(output.c_str());
        writer->SetInputData(simplex);
        writer->Update();
    }
    if (visualize) {
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(simplex);

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();
        vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

        ren->AddActor(actor);
        renwin->AddRenderer(ren);
        iren->SetRenderWindow(renwin);
        iren->Start();
    }
}

/** 
 * Used for performing sort operation. 
 */
bool compareFunction (int i, int j) { 
	return (i < j); 
}

/**
 * Checks whether the two vtkIdLists iList and jList are the same by comparing their respective points Ids. 
 * 
 * @param iList - List containing the point ids of one vtkCell. 
 * @param jList - List containing the point ids of another vtkCell. 
 * @return - True if the two lists have the same point ids. False otherwise. 
 */
bool isDuplicateCell(vtkSmartPointer<vtkIdList> iList, vtkSmartPointer<vtkIdList> jList) {
    // If the two cells have different number of points, then they cannot be duplicate cells.     
    if (iList->GetNumberOfIds() != jList->GetNumberOfIds()) {
        return false;
    }
    else {
        // Copy point Ids into vectors. 
        std::vector<int> ic;
        for (int i = 0; i < iList->GetNumberOfIds(); i++) {
            ic.push_back(iList->GetId(i));
        }
                
        std::vector<int> jc;
        for (int j = 0; j < jList->GetNumberOfIds(); j++) {
            jc.push_back((jList->GetId(j)));
        }
        
        // Sort the vectors. 
        std::sort(ic.begin(), ic.end(), compareFunction);
        std::sort(jc.begin(), jc.end(), compareFunction);
        
        // Assuming that the cells are duplicate, compare each respective cell point Ids. 
        // If the Ids are not the same, then the cells are not duplicate. 
        bool same = true;
        for (int i = 0; i < ic.size(); i++) {            
            if (ic.at(i) != jc.at(i)) {
                same = false;
                break;
            }
        }
        
        //if (same) cout << "Duplicates detected. " << endl;
        
        return same;
    }    
}

/**
 * Checks whether the two vtkCells iCell and jCell are the same by comparing their respective point Ids. 
 * 
 * @param iCell - One vtkCell
 * @param jCell - Another vtkCell
 * @return - True if the two cells have the same point ids. False otherwise. 
 */
bool isDuplicateCell(vtkSmartPointer<vtkCell> iCell, vtkSmartPointer<vtkCell> jCell) {
    // If the two cells have different number of points, then they cannot be duplicate cells. 
    if (iCell->GetNumberOfPoints() != jCell->GetNumberOfPoints()) {
        return false;
    }
    else {
        // Copy point Ids into vectors. 
        std::vector<int> ic;
        for (int i = 0; i < iCell->GetNumberOfPoints(); i++) {
            ic.push_back(iCell->GetPointId(i));
        }
        
        std::vector<int> jc;
        for (int j = 0; j < jCell->GetNumberOfPoints(); j++) {
            jc.push_back((jCell->GetPointId(j)));
        }
        
        // Sort the vectors. 
        std::sort(ic.begin(), ic.end(), compareFunction);
        std::sort(jc.begin(), jc.end(), compareFunction);
        
        // Assuming that the cells are duplicate, compare each respective cell point Ids. 
        // If the Ids are not the same, then the cells are not duplicate. 
        bool same = true;
        for (int i = 0; i < ic.size(); i++) {            
            if (ic.at(i) != jc.at(i)) {
                same = false;
                break;
            }            
        }
        return same;
    }    
}

/**
 * Determines whether two triangular cells are neighbors. list1[] and list2[] contain the points of the two cells.  
 * Assumes that the cells have only 3 points. 
 * 
 * @param list1 - List of point ids for one cell. 
 * @param list2 - List of point ids for another cell. 
 * @return - True if the two cells are neighbors. False otherwise. 
 */
bool isNeighbor(int list1[], int list2[]) {    
    bool b = false;
    int count = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (list1[i] == list2[j]) {
                count = count + 1;
            }
        }
    }

    if (count == 2) b = true;
    else b = false;

    return b;
}

/**
 * Rearranges the cell Ids in *list according to their spatial occurance. 
 * 
 * @param list - List of cell Ids that surround a point in the triangular mesh.  
 * @param objectData - The triangular mesh. 
 * @return - A list of cell Ids that are arranged in their order of their spatial occurance. 
 */
vtkSmartPointer<vtkIdList> orderCellList(vtkSmartPointer<vtkIdList> list, vtkSmartPointer<vtkPolyData> objectData) {   
    vtkSmartPointer<vtkIdList> sortedCellList = list;   // Contains the list of cell numbers in ascending/descending order

    int initialSize_sortedCellList = sortedCellList->GetNumberOfIds();

    int *unsortedCellList = new int[initialSize_sortedCellList]; // Contains the list of cell numbers in order of occurance. 
    int unsorted_p = 0; // Pointer to the unsorted list. 

    // Remove the first cell number in the sorted list and place it in the unsorted list. 
    unsortedCellList[unsorted_p] = sortedCellList->GetId(0);		
    unsorted_p = unsorted_p + 1;
    sortedCellList->DeleteId(sortedCellList->GetId(0));

    // While the sorted list is not empty. 
    int count = 0;
    int j = 0;
    while (sortedCellList->GetNumberOfIds() > 0) {
        // Take the last cell number from the unsorted list. Use this number to retrieve the vtkCell from the polydata object. 
        // Then get the list of point numbers of that vtkCell. Then assign those point numbers to an int array. 
        vtkSmartPointer<vtkCell> iCell = objectData->GetCell(unsortedCellList[unsorted_p - 1]);			
        vtkSmartPointer<vtkIdList> iList = vtkSmartPointer<vtkIdList>::New(); 
        iList = iCell->GetPointIds();

        int *array1 = new int[iList->GetNumberOfIds()];
        for (int a = 0; a < iList->GetNumberOfIds(); a++) {
            array1[a] = iList->GetId(a);
        }

        // Take the jth cell number from the sorted list. Use this number to retrieve the vtkCell from the polydata object.
        // Then get the list of point numbers of that vtkCell. Then assign those point numbers to an int array. 
        vtkSmartPointer<vtkCell> jCell = objectData->GetCell(sortedCellList->GetId(j));
        vtkSmartPointer<vtkIdList> jList = vtkSmartPointer<vtkIdList>::New();
        jList = jCell->GetPointIds();

        int *array2 = new int[jList->GetNumberOfIds()];
        for (int b = 0; b < jList->GetNumberOfIds(); b++) {
            array2[b] = jList->GetId(b);
        }

        // If the two cells iCell and jCell are neighboring cells, then remove the jth cell number
        // from the sorted list and append it to the unsorted list. Reset the while-loop. 
        if (isNeighbor(array1, array2) == true) {
            unsortedCellList[unsorted_p] = sortedCellList->GetId(j);
            sortedCellList->DeleteId(sortedCellList->GetId(j));
            unsorted_p = unsorted_p + 1;

            j = 0;
        }

        if (j == sortedCellList->GetNumberOfIds() - 1) {
            j = 0;
        }
        else {
            j = j + 1;
        }
        
        count = count + 1;
        
        if (count >= 100) {
            vtkSmartPointer<vtkIdList> ret = vtkSmartPointer<vtkIdList>::New();
            ret->InsertNextId(-1);
            return ret;
        }
    }

    vtkSmartPointer<vtkIdList> ret = vtkSmartPointer<vtkIdList>::New();
    for (int i = 0; i < initialSize_sortedCellList; i++) {
        ret->InsertNextId(unsortedCellList[i]);
    }
    
    return ret;
}

/**
 * Converts a multimaterial triangular mesh into a simplex mesh, while preserving 
 * material information. The input triangular mesh is assumed to be in vtkPolyData format. 
 * The centroids of the triangles of the input mesh are the points of the simplex mesh. 
 * So, the triangular cell's material info becomes the material info of the simplex point. 
 * 
 * @param input - The address of the input vtkPolyData file. 
 * @param output - Output vtkPolyData file location. 
 */
vtkSmartPointer<vtkPolyData> convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh(std::string input) {
	double p0[3], p1[3], p2[3];
	double cellCentroid[3];

    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(input.c_str());
    reader->Update();
    
    vtkSmartPointer<vtkPolyData> pdata = reader->GetOutput();
    
    vtkSmartPointer<vtkPoints> points = pdata->GetPoints();    
    
    int numberOfCells = pdata->GetNumberOfCells();
    int numberOfPoints = pdata->GetNumberOfPoints();
    
    cout << "Number of cells: " << numberOfCells << endl;
    cout << "Number of points: " << numberOfPoints << endl;        
    cout << "Scalar range: [" << pdata->GetScalarRange()[0] << ", " << pdata->GetScalarRange()[1] << "]" << endl;
        
    vtkSmartPointer<vtkDataArray> scalarArray = pdata->GetCellData()->GetArray(0); // Array 0 contains scalar info, 
    vtkSmartPointer<vtkDataArray> dataArray1 = pdata->GetCellData()->GetArray(1); // Array 1 contains material index, each row having two indices. 0 indicates background. 

    double *range = dataArray1->GetRange();
    cout << "Range: [" << range[0] << ", " << range[1] << "]" << endl;
    int maxMaterialRange = dataArray1->GetRange()[1];
    
    // Material range always start from 2. 
    cout << "Materials range: [2, " << maxMaterialRange << "]" << endl;
    
       
    vtkSmartPointer<vtkPolyData> simplex = vtkSmartPointer<vtkPolyData>::New(); // Simplex mesh of triangulated mesh. 
    vtkSmartPointer<vtkPoints> simplexPoints = vtkSmartPointer<vtkPoints>::New(); // Points of the simplex mesh. 
    vtkSmartPointer<vtkCellArray> simplexVertices = vtkSmartPointer<vtkCellArray>::New(); // To set the vertices of the simplex mesh. 
    vtkSmartPointer<vtkCellArray> simplexCells = vtkSmartPointer<vtkCellArray>::New(); // To set the cells of the simplex mesh. 
    
    typedef vtkSmartPointer<vtkIdList> idList;
    std::vector<idList> candidateCells;
    
    
    // Compute the centroid of each cell in the triangular mesh and set as the vertices of the simplex mesh. 
    vtkSmartPointer<vtkIdList> verticesList = vtkSmartPointer<vtkIdList>::New();
    for (int i = 0; i < numberOfCells; i++) {        
        vtkSmartPointer<vtkCell> cell = pdata->GetCell(i);
        vtkSmartPointer<vtkIdList> cellPoints = cell->GetPointIds();

        //double p0[3], p1[3], p2[3];
        points->GetPoint(cellPoints->GetId(0), p0);
        points->GetPoint(cellPoints->GetId(1), p1);
        points->GetPoint(cellPoints->GetId(2), p2);

        //double cellCentroid[3];
        cellCentroid[0] = (p0[0] + p1[0] + p2[0]) / 3;
        cellCentroid[1] = (p0[1] + p1[1] + p2[1]) / 3;
        cellCentroid[2] = (p0[2] + p1[2] + p2[2]) / 3;

        verticesList->InsertNextId(simplexPoints->InsertNextPoint(cellCentroid));        
    }    
    simplexVertices->InsertNextCell(verticesList);
    

    // For each primary material indices in the triangular mesh. 
    for (int m = 2; m <= maxMaterialRange; m++) {
        int currentMaterial = m; 
        
        // For each point of the triangular mesh. 
        for (int i = 0; i < numberOfPoints; i++) {
            //cout << currentMaterial << ": " << i << endl;
            // Get the cells that surround a point i. 
            vtkSmartPointer<vtkIdList> vertex = vtkSmartPointer<vtkIdList>::New();
            vertex->InsertNextId(i);

            vtkSmartPointer<vtkIdList> cellsSharingPoint = vtkSmartPointer<vtkIdList>::New();
            pdata->GetPointCells(i, cellsSharingPoint);
            
            // Use only current material cells. 
            double *materials = new double[2];
            vtkSmartPointer<vtkIdList> unorderedCellList = vtkSmartPointer<vtkIdList>::New();
            
            for (int j = 0; j < cellsSharingPoint->GetNumberOfIds(); j++) {
                materials = dataArray1->GetTuple(cellsSharingPoint->GetId(j));
                
                if ((int)materials[0] == currentMaterial || (int)materials[1] == currentMaterial) {
                    unorderedCellList->InsertNextId(cellsSharingPoint->GetId(j));
                }                
            }
            
            // Spatially order the cell ids, and add the list of cell ids as a simplex cell. 
            if (unorderedCellList->GetNumberOfIds() > 0) {
                vtkSmartPointer<vtkIdList> organizedCellList = orderCellList(unorderedCellList, pdata);
                
                // For each material's iteration, the shared boundaries will be encountered multiple times. 
                // For shared boundaries, the materials[] array will contain non-zero values. 
                // If the cell is a shared boundary, then store it in candidateCells array 
                // to compare with at later times. 
                if ((int)materials[0] != 0 && (int)materials[1] != 0) {
                    if (organizedCellList->GetId(0) != -1) {
                        // Check if the cell is already created and in candidateCells. 
                        bool alreadyCreated = false;
                        for (int z = 0; z < candidateCells.size(); z++) {
                            vtkSmartPointer<vtkIdList> aCell = candidateCells.at(z);
                            
                            if (isDuplicateCell(aCell, organizedCellList) == true) {
                                alreadyCreated = true;
                                break;
                            }
                        }

                        // If the cell is not created, then add it into candidateCells. 
                        if (alreadyCreated == false) {
                            candidateCells.push_back(organizedCellList);
                        }
                    }
                }
                
                // If the materials[] array contains 0, then that means the cell 
                // is not a shared boundary. It will be encountered once only. 
                // No need for comparing. 
                else {
                    simplexCells->InsertNextCell(organizedCellList);
                }
            }
        }
    }
    
    // Copy all the cells from candidateCells into simplexCells. 
    for (int a = 0; a < candidateCells.size(); a++) {
        simplexCells->InsertNextCell(candidateCells.at(a));
    }
    
	

    // Set simplex points, cells, and material info. 
    simplex->SetPoints(simplexPoints);
    simplex->SetPolys(simplexCells);    
	//simplex->SetFieldData(pointData);
    //simplex->SetVerts(simplexVertices);
    
	    
    // Set the material indices of the simplex points
    // Since the centroid of a triangle in the input mesh is a point of the simplex mesh
    // the triangle's material info is the material info of the simplex point. 
	// In the triangular mesh, scalar info is the first array, and the 2 material indices
	// are the second array. 
	// For the Simplex mesh, the 2 material indices are separated into 2 arrays, and the
	// scalar info is the third array. 

	// Set up two arrays to contain the two material indices
	vtkSmartPointer<vtkIntArray> matvals1 = vtkSmartPointer<vtkIntArray>::New();
	vtkSmartPointer<vtkIntArray> matvals2 = vtkSmartPointer<vtkIntArray>::New();
	
	matvals1->SetNumberOfValues(simplex->GetNumberOfPoints());
	matvals1->SetNumberOfComponents(1);
	matvals1->SetName("PointMaterialIndices1");

	matvals2->SetNumberOfValues(simplex->GetNumberOfPoints());
	matvals2->SetNumberOfComponents(1);
	matvals2->SetName("PointMaterialIndices2");

	for (int a = 0; a < simplex->GetNumberOfPoints(); a++) {
		double *d = dataArray1->GetTuple2(a);
		matvals1->SetValue(a, d[0]);
		matvals2->SetValue(a, d[1]);
	}

	// Add all three arrays to the simplex mesh as vtkPointData
	simplex->GetPointData()->AddArray(matvals1); // First material indices
	simplex->GetPointData()->AddArray(matvals2); // Second material indices
	simplex->GetPointData()->AddArray(scalarArray); // Scalars

    cout << "Unoriented simplex mesh" << endl;
    cout << "\tNumber of points: " << simplex->GetNumberOfPoints() << endl;
    cout << "\tNumber of cells: " << simplex->GetNumberOfCells() << endl;
    
	//vtkSmartPointer<vtkPolyData> oriented = orientSimplexSurfMesh(simplex);
	vtkSmartPointer<vtkPolyDataNormals> orientPolyData = vtkSmartPointer<vtkPolyDataNormals>::New();
	orientPolyData->SplittingOff();
	orientPolyData->ConsistencyOn();
	orientPolyData->FlipNormalsOff();
	orientPolyData->NonManifoldTraversalOff();
	orientPolyData->ComputeCellNormalsOff();
	
	orientPolyData->SetInputData(simplex);
	orientPolyData->Update();

	vtkSmartPointer<vtkPolyData> oriented = orientPolyData->GetOutput();

	if (oriented->GetNumberOfCells() != simplex->GetNumberOfCells() && oriented->GetNumberOfPoints() != simplex->GetNumberOfPoints()) {
		std::cerr << "ERROR: Mesh orientation using vtkPolyDataNormals resulted in a new mesh with different number of cells and/or points." << endl;
	}

    //writeOrVisualizePolyData(oriented, true, false, output);

	return oriented;
}

/**
 * Takes a multi-material Simplex mesh in vtkPolyData format
 * and returns a vector of vtkSmartPointer<vtkPolyData>
 * where each row of the vector contains a separate material mesh. 
 * Each new vtkPolyData contains their three arrays (Mat1 Index, Mat2 Index and Scalar)
 * plus an additional array which contains the original point indices. 
 */
std::vector<vtkSmartPointer<vtkPolyData> > splitMesh(vtkSmartPointer<vtkPolyData> pdata) {
	std::vector<vtkSmartPointer<vtkPolyData> > ret;
	
	pdata->BuildCells();
	pdata->BuildLinks();

	// Read material indices and determine the range of material indices. 
	// MMSimplexMesh's material indices are stored as a vtkPointData of 3 arrays, 
	// each array containing 1 component. Array 0 and Array 1 are the two material 
	// indices, and Array 2 contains the scalar values
	vtkSmartPointer<vtkDataArray> matsArray0 = pdata->GetPointData()->GetArray(0);
	vtkSmartPointer<vtkDataArray> matsArray1 = pdata->GetPointData()->GetArray(1);

	double *range0 = matsArray0->GetRange();
	double *range1 = matsArray1->GetRange();
	
	int maxRange = -1;
	if (range0[1] >= range1[1]) {
		maxRange = (int)range0[1];
	}
	else {
		maxRange = (int)range1[1];
	}

	// For each material starting from 2 to maxRange
	// Since this splitting process relies on using points to get cells containing the points,
	// we first select all the points ids that are not the current material 	
	std::cout << "Starting splitting process..." << std::endl;
	for (int currentMaterial = 2; currentMaterial <= maxRange; currentMaterial++) {
		std::cout << " Working on Material " << currentMaterial << std::endl;
		vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
		ids->SetNumberOfComponents(1);

		for (vtkIdType i = 0; i < pdata->GetNumberOfPoints(); i++) {
			double mat0 = matsArray0->GetTuple1(i);
			double mat1 = matsArray1->GetTuple1(i);

			if (mat0 != currentMaterial && mat1 != currentMaterial) {
				ids->InsertNextValue(i);
			}
		}
		
		// Using the point ids that do not belong to the current material,
		// we use vtkSelectionNode::CONTAINING_CELLS() to extract all cells that
		// contain the selected point ids. 
		// Then doing an inversion selection with vtkSelectionNode::INVERSE()
		// gives only the points and cells that belong to the current material. 
		vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
		selectionNode->SetFieldType(vtkSelectionNode::POINT);
		selectionNode->SetContentType(vtkSelectionNode::INDICES);
		selectionNode->SetSelectionList(ids);
		selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);
		selectionNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1);

		vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
		selection->AddNode(selectionNode);

		vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
		extractSelection->SetInputConnection(0, pdata->GetProducerPort());
		extractSelection->SetInputData(1, selection);
		extractSelection->Update();

		vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
		selected->ShallowCopy(extractSelection->GetOutput());
		
		std::cout << " Done splitting for Material " << currentMaterial << std::endl;
		std::cout << " Converting vtkUnstructuredGrid to vtkPolyData..." << std::endl;

		// Convert vtkUnstructuredGrid into vtkPolyData
		vtkSmartPointer<vtkPolyData> iMatMesh = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> meshPts = selected->GetPoints();
		vtkSmartPointer<vtkCellArray> meshCells = vtkSmartPointer<vtkCellArray>::New();

		// Copy unstructuredgrid cells into polydata cells. 
		for (int j = 0; j < selected->GetNumberOfCells(); j++) {
			vtkSmartPointer<vtkIdList> list = selected->GetCell(j)->GetPointIds();
			meshCells->InsertNextCell(list);
		}

		//std::cout << "Number of PointData arrays: " << selected->GetPointData()->GetNumberOfArrays() << std::endl;
		//std::cout << "Number of CellData arrays: " << selected->GetCellData()->GetNumberOfArrays() << std::endl;
		//std::cout << "Number of FieldData arrays: " << selected->GetFieldData()->GetNumberOfArrays() << std::endl;
		//
		//std::cout << "Names of PointData arrays: " << std::endl;
		//for (int q = 0; q < selected->GetPointData()->GetNumberOfArrays(); q++) {
		//	std::cout << "  " << selected->GetPointData()->GetArray(q)->GetName() << std::endl;
		//}

		// Use array names to locate the required arrays. 
		vtkSmartPointer<vtkDataArray> dataArray0 = selected->GetPointData()->GetArray("PointMaterialIndices1"); // Material index 1
		vtkSmartPointer<vtkDataArray> dataArray1 = selected->GetPointData()->GetArray("PointMaterialIndices2"); // Material index 2
		vtkSmartPointer<vtkDataArray> dataArray2 = selected->GetPointData()->GetArray("scalars"); // Scalars
		vtkSmartPointer<vtkDataArray> dataArray3 = selected->GetPointData()->GetArray("vtkOriginalPointIds"); // Original point indices
		
		// Set points, cells and data arrays. 
		iMatMesh->SetPoints(meshPts);
		iMatMesh->SetPolys(meshCells);
		iMatMesh->GetPointData()->AddArray(dataArray0);
		iMatMesh->GetPointData()->AddArray(dataArray1);
		iMatMesh->GetPointData()->AddArray(dataArray2);
		iMatMesh->GetPointData()->AddArray(dataArray3);

		//std::stringstream ss;
		//std::string MMSimplexMesh = "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterial_Deformation_TwoBoxes\\TwoBoxesUnequalSideBySide_MM2M_Transformed.vtk";
		//ss << MMSimplexMesh.substr(0, MMSimplexMesh.size() - 4) << "_Simplex_Split_Mat_" << currentMaterial << ".vtk";
		//
		//vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		//writer->SetInput(iMatMesh);
		//writer->SetFileName(ss.str().c_str());
		//writer->Write();


		// Make sure that the cells of the simplex mesh are all oriented the same way. 
		vtkSmartPointer<vtkPolyDataNormals> orientPolyData = vtkSmartPointer<vtkPolyDataNormals>::New();
		orientPolyData->SplittingOff();
		orientPolyData->ConsistencyOn();
		orientPolyData->FlipNormalsOff();
		orientPolyData->NonManifoldTraversalOff();
		orientPolyData->ComputeCellNormalsOff();
	
		orientPolyData->SetInputData(iMatMesh);
		orientPolyData->Update();

		vtkSmartPointer<vtkPolyData> oriented = orientPolyData->GetOutput();
		for (int ar = 0; ar < oriented->GetPointData()->GetNumberOfArrays(); ar++) 	{
			oriented->GetPointData()->RemoveArray(ar);
		}
		for (int ar = 0; ar < oriented->GetCellData()->GetNumberOfArrays(); ar++) 	{
			oriented->GetCellData()->RemoveArray(ar);
		}


		oriented->GetPointData()->AddArray(dataArray0);
		oriented->GetPointData()->AddArray(dataArray1);
		oriented->GetPointData()->AddArray(dataArray2);
		oriented->GetPointData()->AddArray(dataArray3);


		if (oriented->GetNumberOfCells() != iMatMesh->GetNumberOfCells() && oriented->GetNumberOfPoints() != iMatMesh->GetNumberOfPoints()) {
			std::cerr << "ERROR: Mesh orientation using vtkPolyDataNormals resulted in a new mesh with different number of cells and/or points." << endl;
		}

		ret.push_back(oriented);
	}

	return ret;
}

/**
 * The input is a multi-material simplex mesh in vtkPolyData format. This polydata should
 * come from convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh()
 */
CSimplexSurf** Split_MultiMaterialSimplexMesh_Into_Separate_CSimplexMesh(vtkSmartPointer<vtkPolyData> multimaterialSimplex) {
	multimaterialSimplex->BuildCells();
	multimaterialSimplex->BuildLinks();

	// Retrieve vtkPointData arrays and determine max range of material indices. 
	vtkSmartPointer<vtkDataArray> matsArray0 = multimaterialSimplex->GetPointData()->GetArray(0);
	vtkSmartPointer<vtkDataArray> matsArray1 = multimaterialSimplex->GetPointData()->GetArray(1);
	vtkSmartPointer<vtkDataArray> scalarArray = multimaterialSimplex->GetPointData()->GetArray(2);
	
	double *range0 = matsArray0->GetRange();
	double *range1 = matsArray1->GetRange();
	
	int maxRange = -1;
	if (range0[1] >= range1[1]) {
		maxRange = (int)range0[1];
	}
	else {
		maxRange = (int)range1[1];
	}

	int numberOfMeshes = maxRange; // Total number of meshes: numberOfMeshes = originalMMSimplexMesh + (No of separate material meshes) = maxRange

	CSimplexSurf **ret = new CSimplexSurf*[numberOfMeshes];
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// First create the whole Multi Material Simplex Mesh and store it as the first mesh in ret
	/////////////////////////////////////////////////////////////////////////////////////////////////
	ret[0] = new(CSimplexSurf);
	ret[0]->Free();
	ret[0]->Allocate(multimaterialSimplex->GetNumberOfPoints(), multimaterialSimplex->GetNumberOfCells());
	ret[0]->SetSplitMesh(false); // Basically saying that this mesh is the whole multi-material simplex mesh.

	if ((matsArray0->GetNumberOfTuples() > 0) && (matsArray1->GetNumberOfTuples() > 0) && (scalarArray->GetNumberOfTuples() > 0)) { // If the dataArray has material data
		for (int i = 0; i < multimaterialSimplex->GetNumberOfPoints(); i++) {
			double *d = multimaterialSimplex->GetPoint(i);
			ret[0]->SetPoint(i, d);

			double m0 = matsArray0->GetTuple1(i);
			double m1 = matsArray1->GetTuple1(i);
			double scalar = scalarArray->GetTuple1(i);

			ret[0]->SetPointMaterialIndices(i, (int)m0, (int)m1, (int)scalar);
		}
		ret[0]->AllocateMultiMaterialNeighbors(2, maxRange); // Material indices always start from 2
	}
	else { // If the dataArray does not have any material data. 
		std::cerr << "ERROR: Mesh does not appear to contain multi-material elements" << endl;
	}

	for (int i = 0; i < multimaterialSimplex->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkCell> cell = multimaterialSimplex->GetCell(i);

		int n = cell->GetNumberOfPoints();
		int *c = new int[n];
		
		for (int j = 0; j < n; j++) {
			c[j] = cell->GetPointId(j);

			ret[0]->SetCell(i, n, c);
		}
	}

	ret[0]->UpdateNeighbors();
	ret[0]->ComputeVolume();
	ret[0]->Equilibrium();
	ret[0]->UpdateParams();
	ret[0]->UpdateMass();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Now split each material mesh as a regular 2-simplex mesh
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	std::vector<vtkSmartPointer<vtkPolyData> > otherMeshes = splitMesh(multimaterialSimplex);
	std::cout << "Number of split meshes: " << otherMeshes.size() << std::endl;

	int count = 1;
	for (int i = 0; i < otherMeshes.size(); i++) {
		vtkSmartPointer<vtkPolyData> iMesh = otherMeshes.at(i);
		iMesh->BuildCells();
		iMesh->BuildLinks();

		//vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
		//writer->SetFileName("data\\tempvtk.vtk");
		//writer->SetInput(iMesh);
		//writer->Write();

		// Retrieve material and point data
		vtkSmartPointer<vtkDataArray> m0Array = iMesh->GetPointData()->GetArray(0);
		vtkSmartPointer<vtkDataArray> m1Array = iMesh->GetPointData()->GetArray(1);
		vtkSmartPointer<vtkDataArray> sArray = iMesh->GetPointData()->GetArray(2);
		vtkSmartPointer<vtkDataArray> originalPointIndicesArray = iMesh->GetPointData()->GetArray(3);
		
		ret[count] = new(CSimplexSurf);
		ret[count]->Free();
		ret[count]->Allocate(iMesh->GetNumberOfPoints(), iMesh->GetNumberOfCells());
		ret[count]->SetSplitMesh(true); // Basically saying that this mesh is a split material mesh of the whole simplex mesh

		// Set simplex points and material data
		for (int j = 0; j < iMesh->GetNumberOfPoints(); j++) {
			double *d = iMesh->GetPoint(j);
			ret[count]->SetPoint(j, d);
			
			double m0 = m0Array->GetTuple1(j);
			double m1 = m1Array->GetTuple1(j);
			double scalar = sArray->GetTuple1(j);
			double ind = originalPointIndicesArray->GetTuple1(j);

			ret[count]->SetPointMaterialIndices(j, (int)m0, (int)m1, (int)scalar);;
			ret[count]->SetOriginalPointIndex(j, (int)ind);
		}

		// Set cells
		for (int k = 0; k < iMesh->GetNumberOfCells(); k++) {
			vtkSmartPointer<vtkCell> cell = iMesh->GetCell(k);

			int n = cell->GetNumberOfPoints();
			int *c = new int[n];
		
			for (int kj = 0; kj < n; kj++) {
				c[kj] = cell->GetPointId(kj);

				ret[count]->SetCell(k, n, c);
			}
		}

		ret[count]->UpdateNeighbors();
		ret[count]->ComputeVolume();
		ret[count]->Equilibrium();
		ret[count]->UpdateParams();
		ret[count]->UpdateMass();

		count = count + 1;
	}

	// Free up the memory used by otherMeshes vector.
	// Using clear() is not going to free up the memory, therefore use swap() with an empty vector. 
	// According to stackoverflow. 
	std::vector<vtkSmartPointer<vtkPolyData> > emptyV;
	otherMeshes.swap(emptyV);

	return ret;
}

CSimplexSurf** Load_MultiMaterialCSimplexMeshFromVTKPolyData(std::string filename) {
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();

	vtkSmartPointer<vtkPolyData> MMSimplexMesh = reader->GetOutput();

	return Split_MultiMaterialSimplexMesh_Into_Separate_CSimplexMesh(MMSimplexMesh);
}