
#include <stdlib.h>
#include <iostream>
#include <sstream>

#include "SimplexMesh.h"
#include "SimplexSurf.h"
#include "Deformation.h"
#include "MMCSimplexSurf.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <vtkFeatureEdges.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageShiftScale.h>
#include <vtkImageCast.h>
#include <vtkImageAnisotropicDiffusion3D.h>

#include <vtkAppendPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkTriangleFilter.h>
#include <vtkMassProperties.h>

#include "CSimplexMeshConverter.h"


#include <vtkDiskSource.h>

using namespace std;

P_float dotproduct(P_float u[3], P_float v[3]) {return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);}
void crossproduct( P_float cp[3], P_float u[3], P_float v[3]) {cp[0]=u[1]*v[2]-u[2]*v[1]; cp[1]=u[2]*v[0]-u[0]*v[2]; cp[2]=u[0]*v[1]-u[1]*v[0];}
void crossproduct( P_float cp[3], const P_float u[3], const P_float v[3]) {cp[0]=u[1]*v[2]-u[2]*v[1]; cp[1]=u[2]*v[0]-u[0]*v[2]; cp[2]=u[0]*v[1]-u[1]*v[0];}


vtkStructuredPoints* Load_dMRI(const char* filename) {
	vtkStructuredPoints *dMRI = vtkStructuredPoints::New();

	vtkSmartPointer<vtkStructuredPointsReader> dMRI_reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	dMRI_reader->SetFileName(filename);
	dMRI_reader->Update();

	vtkSmartPointer<vtkImageCast> Cast = vtkSmartPointer<vtkImageCast>::New(); 
	Cast->SetInputData(dMRI_reader->GetOutput()); 
	Cast->SetOutputScalarTypeToUnsignedShort(); 
	Cast->Update(); 
		
	dMRI->DeepCopy(Cast->GetOutput());
		
	return dMRI;
}


	
CSimplexSurf* generateSphere(double R) {
	int i, j;
	double PI = 3.14159265359;
	
	// icosaedre
	double al = 2 * asin(cos(PI / 5) / sin(PI / 3));	// angle d'inclinaison
	double l = (2 * R) / (tan(PI / 5) * tan(al / 2));	// longueur d'une arête
	double k = (2 * R * R - l * l) / (2 * R);						// plan d'intersection
	double r = sqrt(R * R - k * k);					// rayon du cercle d'intersection
	double da = 2 * PI / 5;	
	double pts[12][3] = {
		{r, 0, k}, {r * cos(da), r * sin(da), k}, {r * cos(2 * da), r * sin(2 * da), k}, {r * cos(3 * da), r * sin(3 * da), k}, {r * cos(4 * da), r * sin(4 * da), k}, {0, 0, R}, {0, 0, -R}, {-r, 0, -k}, {-r * cos(da), -r * sin(da), -k}, {-r * cos(2 * da), -r * sin(2 * da), -k}, {-r * cos(3 * da), -r * sin(3 * da), -k}, {-r * cos(4 * da), -r * sin(4 * da), -k}
	};
	int edges[12][5] = {
		{1, 5, 4, 9, 10}, {2, 5, 0, 10, 11}, {1, 11, 7, 3, 5}, {5, 2, 7, 8, 4}, {5, 3, 8, 9, 0}, {0, 1, 2, 3, 4}, {9, 8, 7, 11, 10}, {11, 6, 8, 3, 2}, {7, 6, 9, 4, 3}, {10, 0, 4, 8, 6}, {0, 9, 6, 11, 1}, {6, 7, 2, 1, 10} 
	};

	CSimplexSurf *ret = new(CSimplexSurf); 

	ret->Free();
	ret->Allocate(60, 32);
	
	double p[3];
	for (i = 0; i < 12; i++) {
		for (j = 0; j < 5; j++) {
			p[0] = (2 * pts[i][0] / 3 + pts[edges[i][j]][0] / 3) + 50;
			p[1] = (2 * pts[i][1] / 3 + pts[edges[i][j]][1] / 3);
			p[2] = (2 * pts[i][2] / 3 + pts[edges[i][j]][2] / 3);
			ret->SetPoint(5 * i + j, p);
		}
	}
	int pentas[12][5] = {
		{0, 1, 2, 3, 4}, {5, 6, 7, 8, 9}, {10, 11, 12, 13, 14}, {15, 16, 17, 18, 19}, 
		{20, 21, 22, 23, 24}, {25, 26, 27, 28, 29}, {30, 31, 32, 33, 34}, {35, 36, 37, 38, 39}, 
		{40, 41, 42, 43, 44}, {45, 46, 47, 48, 49}, {50, 51, 52, 53, 54}, {55, 56, 57, 58, 59}
	};

	int hexas[20][6] = {
		{0, 7, 6, 26, 25, 1}, {1, 25, 29, 20, 24, 2}, {2, 24, 23, 47, 46, 3}, {3, 46, 45, 51, 50, 4}, 
		{4, 50, 54, 8, 7, 0}, {15, 28, 27, 14, 13, 16}, {16, 13, 12, 39, 38, 17}, {17, 38, 37, 40, 44, 18},
		{18, 44, 43, 22, 21, 19}, {19, 21, 20, 29, 28, 15}, {6, 5, 10, 14, 27, 26}, {5, 9, 58, 57, 11, 10}, 
		{9, 8, 54, 53, 59, 58}, {53, 52, 34, 33, 55, 59}, {52, 51, 45, 49, 30, 34}, {12, 11, 57, 56, 35, 39},
		{37, 36, 32, 31, 41, 40}, {23, 22, 43, 42, 48, 47}, {31, 30, 49, 48, 42, 41}, {33, 32, 36, 35, 56, 55}
	};

	for (i = 0; i < 12; i++) ret->SetCell(i, 5, pentas[i]);
	for (i = 0; i < 20; i++) ret->SetCell(i + 12, 6, hexas[i]);

	ret->UpdateNeighbors();	
	ret->ComputeVolume();
	ret->Equilibrium();
	ret->UpdateParams();
	ret->UpdateMass();

	return ret;
}

//
//         *___________#___________ @
//         /|          /|         /|
//        / |         / |        / |
//      */__|_______#/__|______@/  |
//       |  |        |  |       |  |
//       | *|________|__|_______|__/@
//       | /         | / #      | /
//       |/__________|/_________|/
//       *           #          @
//
// * -> Material 1 Points
// @ -> Material 2 Points
// # -> Points on the shared boundary
//
CSimplexSurf* generateSyntheticSimpleMultiMaterialBox() {
	double pts[12][3] = {
		{0, 0, 0},  // 0
		{0, 1, 0},  // 1 
		{0, 1, 1},  // 2
		{0, 0, 1},  // 3
		{1, 0, 0},  // 4
		{1, 1, 0},  // 5
		{1, 1, 1},  // 6
		{1, 0, 1},  // 7
		{2, 0, 0},  // 8
		{2, 1, 0},  // 9
		{2, 1, 1},  // 10
		{2, 0, 1}   // 11
	};
	
	int cells[11][4] = {
		{0, 1, 2, 3}, 
		{0, 4, 5, 1}, 
		{1, 5, 6, 2}, 
		{4, 7, 6, 5}, 
		{4, 8, 9, 5}, 
		{8, 11, 10, 9}, 
		{5, 9, 10, 6},
		{4, 7, 11, 8}, 
		{7, 6, 10, 11}, 
		{0, 3, 7, 4}, 
		{3, 2, 6, 7}
	};

	CSimplexSurf *cs = new(CSimplexSurf);
	cs->Free();
	cs->Allocate(12, 11);

	// Set points
	for (int i = 0; i < 12; i++) {
		cs->SetPoint(i, pts[i][0], pts[i][1], pts[i][2]);
	}

	// Set cells
	for (int i = 0; i < 11; i++) {
		cs->SetCell(i, 4, cells[i]);
	}

	// Initialize multimaterial neighbors array
	cs->AllocateMultiMaterialNeighbors(2, 3);

	// Set point material indices
	// For material 1
	for (int i = 0; i < 4; i++) {
		cs->SetPointMaterialIndices(i, 0, 2, 5);
	}
	// For shared boundary
	for (int i = 4; i < 8; i++) {
		cs->SetPointMaterialIndices(i, 2, 3, 12);
	}
	// For material 2
	for (int i = 8; i < 12; i++) {
		cs->SetPointMaterialIndices(i, 0, 3, 9);
	}

	cs->UpdateNeighbors();
	cs->ComputeVolume();
	cs->Equilibrium();
	cs->UpdateParams();
	cs->UpdateMass();
	cs->UpdateNormals();
	return cs;
}

CSimplexSurf* generateSyntheticSimpleSimplexBox() {
	double pts[12][3] = {
		{0, 0, 0},  // 0
		{0, 1, 0},  // 1 
		{0, 1, 1},  // 2
		{0, 0, 1},  // 3
		{1, 0, 0},  // 4
		{1, 1, 0},  // 5
		{1, 1, 1},  // 6
		{1, 0, 1},  // 7
	};

	int cells[6][4] = {
		{0, 1, 2, 3}, 
		{0, 4, 5, 1}, 
		{1, 5, 6, 2}, 
		{4, 7, 6, 5}, 		
		{0, 3, 7, 4}, 
		{3, 2, 6, 7}
	};

	CSimplexSurf *cs = new(CSimplexSurf);
	cs->Free();
	cs->Allocate(8, 6);

	// Set points
	for (int i = 0; i < 8; i++) {
		cs->SetPoint(i, pts[i][0], pts[i][1], pts[i][2]);
	}

	// Set cells
	for (int i = 0; i < 6; i++) {
		cs->SetCell(i, 4, cells[i]);
	}

	cs->UpdateNeighbors();
	cs->ComputeVolume();
	cs->Equilibrium();
	cs->UpdateParams();
	cs->UpdateMass();
	cs->UpdateNormals();
	return cs;
}

vtkStructuredPoints* loadVTKStructuredPointsVolume(const char *filename) {
	vtkStructuredPoints *ret = vtkStructuredPoints::New();

	vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName(filename);
	reader->Update();
	ret->DeepCopy(reader->GetOutput());

	return ret;
}

int BasicDeformation(std::string triangleMeshFilename, std::string volumeFilename, int numberOfDeformableMeshes) {

	//std::string triangleMeshFilename = "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Cropped_Labels_And_Colins27_t1\\Object_39-47_MM2M_World.vtk";
	//std::string volumeFilename = "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Cropped_Labels_And_Colins27_t1\\colins_t1_STN_SN_cropped.vtk";
	//std::string triangleMeshFilename = "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Object39_close_open\\Object39_close_open_2Manifold_Smoother_Transformed.vtk";
	//std::string volumeFilename = "E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Object39_close_open\\labels_on_colin_Nov2010_minc2_mincreshape_xyz_Object39.vtk";

	//int numberOfDeformableMeshes = 2;
	
	vtkSmartPointer<vtkPolyData> pdata = convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh(triangleMeshFilename);
	CSimplexSurf **meshes = Split_MultiMaterialSimplexMesh_Into_Separate_CSimplexMesh(pdata);
	Deformation **def = new Deformation*[numberOfDeformableMeshes + 1]; // def[0] is a dummy. Used for consistency for meshes[i]

	for (int m = 0; m <= numberOfDeformableMeshes; m++) {
		std::stringstream ss;
		ss << triangleMeshFilename.substr(0, triangleMeshFilename.size() - 4) << "_Simplex_Split[" << m << "].vtk";
		meshes[m]->writeCSimplexMeshAsVTKPolyData(ss.str());
	}


	for (int i = 1; i <= numberOfDeformableMeshes; i++) {

		meshes[i]->UpdateNeighbors();
		meshes[i]->ComputeVolume();
		meshes[i]->Equilibrium();
		meshes[i]->UpdateParams();
		meshes[i]->UpdateMass();

		meshes[i]->SetGamma(0.0);
		meshes[i]->SetTimeStep(0.1);

		// Load image file.
		vtkStructuredPoints *dMRI = Load_dMRI(volumeFilename.c_str());
	
		// Smooth image
		double difffactor = 0.2;
		double difftresh = 50.0;
		int nbit = 4;
		int dim = 3;
		SmoothVolAnisotropic(dMRI, difffactor, difftresh, nbit, dim);
		
	
		// Compute gradient image
		vtkStructuredPoints *dMRI_grad = GradientMagnitude(dMRI, dim);
		//vtkStructuredPoints *dMRI_grad = Gradient(dMRI, dim);

		//vtkSmartPointer<vtkStructuredPointsWriter> rw_writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
		//rw_writer->SetInput(dMRI);
		//std::stringstream rwss;
		//rwss << volumeFilename.substr(0, volumeFilename.size() - 4) << "_ReadWrite.vtk";
		//rw_writer->SetFileName(rwss.str().c_str());
		//rw_writer->Write();

		//vtkSmartPointer<vtkStructuredPointsWriter> sp_writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
		//sp_writer->SetInput(dMRI);
		//std::stringstream spss;
		//spss << volumeFilename.substr(0, volumeFilename.size() - 4) << "_SmoothAnisotropicVol.vtk";
		//sp_writer->SetFileName(spss.str().c_str());
		//sp_writer->Write();

		//vtkSmartPointer<vtkStructuredPointsWriter> g_writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
		//g_writer->SetInput(dMRI_grad);
		//std::stringstream grss;
		//grss << volumeFilename.substr(0, volumeFilename.size() - 4) << "_GradientMagnitude.vtk";
		////grss << volumeFilename.substr(0, volumeFilename.size() - 4) << "_Gradient.vtk";
		//g_writer->SetFileName(grss.str().c_str());
		//g_writer->Write();
	
		// Compute bounding box
		P_float bounds[3][2]; 
		GetBoundingBox(dMRI, NULL, bounds);

		// Set image and gradient image for mesh
		meshes[i]->SetMRI(dMRI, dMRI_grad, 0, bounds); 
		meshes[i]->SetMRI_GradientMRI(dMRI_grad); 

		//// Set force parameters
		//bool val_surfequ = true;
		//double alpha_surfequ = 0.5;
		//double alpha_lim_surfequ = 0.1;
		//double alpha_incr_surfequ = 0.0;
		//meshes[i]->SetInternalForceSurfEqu(val_surfequ, alpha_surfequ, alpha_lim_surfequ, alpha_incr_surfequ);

		//bool val_refshape = true;
		//double alpha_refshape = 0.5;
		//double alpha_lim_refshape = 0.3;
		//double alpha_incr_refshape = 0.0;
		//meshes[i]->SetInternalForceRefShape(val_refshape, alpha_refshape, alpha_lim_refshape, alpha_incr_refshape);
	
		bool val_laplacian = true;
		double alpha_laplacian = 0.3;
		double alpha_lim_alpha_laplacian = 0.1;
		double alpha_incr_alpha_laplacian = 0.0;
		meshes[i]->SetInternalForceLaplacian(val_laplacian, alpha_laplacian, alpha_lim_alpha_laplacian, alpha_incr_alpha_laplacian);
	
		

		double size_gradientmri;	
		double depth_gradientmri;// = 20;
		bool opposite;
		if (i == 1) {
			depth_gradientmri = 5;
			size_gradientmri = 0.1;
			opposite = false;
		}
		else if (i == 2) {
			depth_gradientmri = 5;
			size_gradientmri = 0.1;
			opposite = false;
		}
		
		bool val_gradientmri = true;
		double alpha_gradientmri = 0.8;
		double alpha_lim_gradientmri = 0.6;
		double alpha_incr_gradientmri = 0.0;
		vtkStructuredPoints *vol = dMRI_grad;
		
		meshes[i]->SetExternalForceGradientMRI(val_gradientmri, size_gradientmri, depth_gradientmri, alpha_gradientmri, alpha_lim_gradientmri, alpha_incr_gradientmri, vol, opposite);

		
		// Regularization
		bool en = false;
		bool smooth = true;
		int transform = 3;
		double lambda = 0.0;
		double lambda_lim = 0.7;

		double lambda_incr = 400.00;
		lambda_incr = (lambda_lim - (double)lambda) / lambda_incr;
		meshes[i]->SetExternalForceRegularization(en, smooth, transform, lambda, lambda_lim, lambda_incr);

		double timestep = 0.1;
		double alpha = 1.0;
		//Deformation *def = new Deformation(meshes[i], timestep, alpha);
		def[i] = new Deformation(meshes[i], timestep, alpha);
	
		//for (int iter = 0; iter <= 200; iter++) {
		//	cout << "Iteration: " << iter << endl;
		//	def->ImpEuler();

		//	if (iter % 10 == 0) {
		//		std::stringstream ss;
		//		ss << "DeformResults\\DeformedMesh_" << iter << ".vtk";
		//		def->getMesh()->writeCSimplexMeshAsVTKPolyData(ss.str().c_str());
		//	}
		//}
	}

	cout << "\nBeginning deformation..." << endl;
	for (int iter = 0; iter <= 200; iter++) {
		for (int j = 1; j <= numberOfDeformableMeshes; j++) { // Do the deformation for each split mesh. 
			def[j]->ImpEuler();
			
			//CSimplexSurf *jmesh = def[j]->getMesh();
			//std::vector<double*> sharedBoundaryPts;
			//for (int p = 0; p < jmesh->GetNumberOfPoints(); p++) {
			//	int *mats = jmesh->GetPointMaterialIndices(p);
			//	
			//	if (mats[0] != 0 && mats[1] != 0) {
			//		double *pt = jmesh->GetPoint(p);
			//		sharedBoundaryPts.push_back(pt);
			//	}
			//}

			//vtkSmartPointer<vtkPolyData> sb = vtkSmartPointer<vtkPolyData>::New();
			//vtkSmartPointer<vtkPoints> sbpts = vtkSmartPointer<vtkPoints>::New();
			//vtkSmartPointer<vtkCellArray> sbarray = vtkSmartPointer<vtkCellArray>::New();
			//
			//for (int p = 0; p < sharedBoundaryPts.size(); p++) {
			//	double *qp = sharedBoundaryPts.at(p);
			//	sbpts->InsertNextPoint(qp[0], qp[1], qp[2]);
			//	
			//	vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
			//	list->InsertNextId(p);
			//	sbarray->InsertNextCell(list);
			//}		

			//sb->SetPoints(sbpts);
			//sb->SetVerts(sbarray);


			//std::stringstream sbss;
			//sbss << "DeformResults\\DeformedMesh_Split[" << j << "]_iter_" << iter << "_PointsOnSharedBoundary.vtk";
			//vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
			//writer->SetInput(sb);
			//writer->SetFileTypeToASCII();
			//writer->SetFileName(sbss.str().c_str());
			//writer->Write();

		//	if (iter % 10 == 0) {
			std::stringstream ss;
			ss << "DeformResults\\DeformedMesh_Split[" << j << "]_iter_" << iter << ".vtk";
			def[j]->getMesh()->writeCSimplexMeshAsVTKPolyData(ss.str().c_str());
				
			cout << "  Iteration " << iter << " for mesh[" << j << "] complete." << endl;
		}

		// Update the parent multi-material mesh's points using the points from the split meshes. 
		// For simplicity, copy all points of all split meshes into a new array newPoints[][].
		// Create another array ptsCounter[] that will contain the number of times an addition 
		// occurs for a row in newPoints[][].
		// Then average all the values in newPoints[][] with ptsCounter[]. 
		// At this stage, the points between material and background should be the same for split meshes 
		// and the updated parent multi-material mesh. Only the points on the shared boundary might be different.
		// Feedback these different shared boundary points from newly updated parent mesh to the split meshes. 
		cout << "  Updating parent multi-material mesh..." << endl << endl;
		
		double **newPoints = new double*[meshes[0]->GetNumberOfPoints()]; // Add all new point to this array
		int *ptsCounter = new int[meshes[0]->GetNumberOfPoints()]; // Increase by 1 for each new point added to newPoints[][]. Used later for averaging

		// Initialize arrays.
		for (int t = 0; t < meshes[0]->GetNumberOfPoints(); t++) {
			newPoints[t] = new double[3];
			newPoints[t][0] = 0;
			newPoints[t][1] = 0;
			newPoints[t][2] = 0;

			ptsCounter[t] = 0;
		}

		// Add points from all split meshes
		for (int j = 1; j <= numberOfDeformableMeshes; j++) {
			CSimplexSurf *jmesh = def[j]->getMesh();

			for (int k = 0; k < jmesh->GetNumberOfPoints(); k++) {
				int originalId = jmesh->GetOriginalPointIndex(k);
				double *pt = jmesh->GetPoint(k);				

				newPoints[originalId][0] = newPoints[originalId][0] + pt[0];
				newPoints[originalId][1] = newPoints[originalId][1] + pt[1];
				newPoints[originalId][2] = newPoints[originalId][2] + pt[2];
				ptsCounter[originalId] = ptsCounter[originalId] + 1;
			}
		}

		// Averaging
		for (int a = 0; a < meshes[0]->GetNumberOfPoints(); a++) {
			newPoints[a][0] = newPoints[a][0] / ptsCounter[a];
			newPoints[a][1] = newPoints[a][1] / ptsCounter[a];
			newPoints[a][2] = newPoints[a][2] / ptsCounter[a];

			//cout << "  New computed value for point " << a << ": "<< newPoints[a][0] << ", " << newPoints[a][1] << ", " << newPoints[a][2] << endl;

			meshes[0]->SetPoint(a, newPoints[a][0], newPoints[a][1], newPoints[a][2]);
		}

		std::stringstream ss0;
	    ss0 << "DeformResults\\DeformedMesh_Split[0]_iter_" << iter << ".vtk";
		meshes[0]->writeCSimplexMeshAsVTKPolyData(ss0.str().c_str());

		//Feedback stage
		for (int a = 1; a <= numberOfDeformableMeshes; a++) {
			// Regular points, i.e. points between material and background
			// are essentially the same in split meshes and parent mesh.
			// Only points on shared boundary might differ.
			// So only considering shared boundary points for each split mesh
			CSimplexSurf *c = def[a]->getMesh();
			for (int b = 0; b < c->GetNumberOfPoints(); b++) {
				int *mats = c->GetPointMaterialIndices(b);
				if (mats[0] != 0 && mats[1] != 0) {
					int originalIndex = c->GetOriginalPointIndex(b);
					double *d = meshes[0]->GetPoint(originalIndex);
					c->SetPoint(b, d[0], d[1], d[2]);
				}
			}

			//std::stringstream ssc;
			//ssc << "DeformResults\\DeformedMesh_Split[" << a << "]_iter_" << iter << "_Feedback.vtk";
			//c->writeCSimplexMeshAsVTKPolyData(ssc.str().c_str());
		}

		// Free up memory, avoid memory leaks
		for (int w = 0; w < meshes[0]->GetNumberOfPoints(); w++) {
			delete [] newPoints[w];
		}

		delete [] newPoints;
		delete [] ptsCounter;
	}






	/*
	vtkSmartPointer<vtkPolyData> pdata = convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterial_Deformation_TwoBoxes\\TwoBoxesUnequalSideBySide_MM2M_Transformed.vtk");
	//vtkSmartPointer<vtkPolyData> pdata = convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Object39_close_open\\Object39_close_open_2Manifold_Smoother_Transformed.vtk");

	vtkSmartPointer<vtkDataArray> matsArray0 = pdata->GetPointData()->GetArray(0);
	vtkSmartPointer<vtkDataArray> matsArray1 = pdata->GetPointData()->GetArray(1);
	vtkSmartPointer<vtkDataArray> scalarArray = pdata->GetPointData()->GetArray(2);

	double *range0 = matsArray0->GetRange();
	double *range1 = matsArray1->GetRange();
	
	int maxRange = -1;
	if (range0[1] >= range1[1]) {
		maxRange = (int)range0[1];
	}
	else {
		maxRange = (int)range1[1];
	}
    
	pdata->BuildCells();
	pdata->BuildLinks();

	CSimplexSurf *cs = new(CSimplexSurf);
	cs->Free();
	cs->Allocate(pdata->GetNumberOfPoints(), pdata->GetNumberOfCells());
	
	if ((matsArray0->GetNumberOfTuples() > 0) && (matsArray1->GetNumberOfTuples() > 0) && (scalarArray->GetNumberOfTuples() > 0)) { // If the dataArray has material data
		for (int i = 0; i < pdata->GetNumberOfPoints(); i++) {
			double *d = pdata->GetPoint(i);
			cs->SetPoint(i, d);

			double m0 = matsArray0->GetTuple1(i);
			double m1 = matsArray1->GetTuple1(i);
			double scalar = scalarArray->GetTuple1(i);

			cs->SetPointMaterialIndices(i, (int)m0, (int)m1, (int)scalar);
		}
		cs->AllocateMultiMaterialNeighbors(2, maxRange); // Material indices always start from 2
	}
	else { // If the dataArray does not have any material data. 
		for (int i = 0; i < pdata->GetNumberOfPoints(); i++) {
			double *d = pdata->GetPoint(i);
			cs->SetPoint(i, d);
		}
	}

	for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkCell> cell = pdata->GetCell(i);

		int n = cell->GetNumberOfPoints();
		int *c = new int[n];
		
		for (int j = 0; j < n; j++) {
			c[j] = cell->GetPointId(j);

			cs->SetCell(i, n, c);
		}
	}
	
	cs->writeCSimplexMeshAsVTKPolyData("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterial_Deformation_TwoBoxes\\TwoBoxesUnequalSideBySide_MM2M_Transformed_Simplex.vtk");
	

	//vtkSmartPointer<vtkPolyData> pdata = convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterial_Deformation_TwoBoxes\\TwoBoxesUnequalSideBySide_MM2M_Transformed.vtk");
	//CSimplexSurf **meshes = SplitMultiMaterialSimplexMesh_Into_Separate_CSimplexMesh(pdata);



	cs->SaveMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\SyntheticSimpleMultiMaterialBox\\SyntheticSimpleMultiMaterialBox_Simplex.defm", true, false);
	cs->writeCSimplexMeshAsVTKPolyData("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\SyntheticSimpleMultiMaterialBox\\SyntheticSimpleMultiMaterialBox_Simplex.vtk");

	//cs->SaveMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Object39_close_open\\Object39_close_open_2Manifold_Smoother_Transformed_Simplex.defm", true, false);
	//cs->writeCSimplexMeshAsVTKPolyData("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Object39_close_open\\Object39_close_open_2Manifold_Smoother_Transformed_Simplex.vtk");

	cs->UpdateNeighbors();
	cs->ComputeVolume();
	cs->Equilibrium();
	cs->UpdateParams();
	cs->UpdateMass();
	
	cs->SaveMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\SyntheticSimpleMultiMaterialBox\\SyntheticSimpleMultiMaterialBox_Simplex.defm", true, false);
	cs->writeCSimplexMeshAsVTKPolyData("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\SyntheticSimpleMultiMaterialBox\\SyntheticSimpleMultiMaterialBox_Simplex.vtk");


	//cs = cs->IncreaseResolution();
	//cs->writeCSimplexMeshAsVTKPolyData("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\SyntheticSimpleMultiMaterialBox\\SyntheticSimpleMultiMaterialBox_Simplex_Res1.vtk");
	
	
	//CSimplexSurf *cs = new(CSimplexSurf);
	//cs->LoadMesh("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\Object39\\Object39_2Manifold_Smoother_Transformed_Half.defm");
	*/
	
	
	/*
	cs->SetGamma(0.0);
	cs->SetTimeStep(1.0);

	
	cs->UpdateNeighbors();
	cs->ComputeVolume();
	cs->Equilibrium();
	cs->UpdateParams();
	cs->UpdateMass();
	

	double *n = cs->GetNormal(4);
	cout << n[0] << ", " << n[1] << ", " << n[2] << endl;

	//double *mass = cs->GetMass(10);
	//cout << mass[0] << ", " << mass[1] << ", " << mass[2] << endl;

//	cs->SaveMesh("translatedmesh.defm", true, false);

	//cs->LoadMesh("Sphere_R100_Res1.defm");
	//cs->writeCSimplexMeshAsVTKPolyData("DeformResults\\OriginalMesh.vtk");

	//cs->SetTimeStep(1.0);


	
	// Load image file. 
	vtkStructuredPoints *dMRI = Load_dMRI("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\SyntheticSimpleMultiMaterialBox\\Vol_SyntheticSimpleMultiMaterialBox.vtk");
	
	// Smooth image
	double difffactor = 0.2;
	double difftresh = 50.0;
	int nbit = 4;
	int dim = 3;
	SmoothVolAnisotropic(dMRI, difffactor, difftresh, nbit, dim);
	
	// Compute gradient image
	vtkStructuredPoints *dMRI_grad = Gradient(dMRI, dim);
	
	// Compute bounding box
	P_float bounds[3][2]; 
	GetBoundingBox(dMRI, NULL, bounds);

	// Set image and gradient image for mesh
	cs->SetMRI(dMRI, dMRI_grad, 0, bounds); 
	cs->SetMRI_GradientMRI(dMRI_grad); 

	// Set force parameters
	bool val_surfequ = true;
	double alpha_surfequ = 0.9;
	double alpha_lim_surfequ = 0.1;
	double alpha_incr_surfequ = 0.0;
	cs->SetInternalForceSurfEqu(val_surfequ, alpha_surfequ, alpha_lim_surfequ, alpha_incr_surfequ);

	bool val_refshape = true;
	double alpha_refshape = 0.9;
	double alpha_lim_refshape = 0.3;
	double alpha_incr_refshape = 0.0;
	cs->SetInternalForceRefShape(val_refshape, alpha_refshape, alpha_lim_refshape, alpha_incr_refshape);
	
	bool val_laplacian = true;
	double alpha_laplacian = 0.9;
	double alpha_lim_alpha_laplacian = 0.1;
	double alpha_incr_alpha_laplacian = 0.0;
	cs->SetInternalForceLaplacian(val_laplacian, alpha_laplacian, alpha_lim_alpha_laplacian, alpha_incr_alpha_laplacian);
	
	bool val_gradientmri = true;
	double size_gradientmri = 0.1;
	double depth_gradientmri = 20;
	double alpha_gradientmri = 0.8;
	double alpha_lim_gradientmri = 0.6;
	double alpha_incr_gradientmri = 0.0;
	vtkStructuredPoints *vol = dMRI_grad;
	bool opposite = true;
	cs->SetExternalForceGradientMRI(val_gradientmri, size_gradientmri, depth_gradientmri, alpha_gradientmri, alpha_lim_gradientmri, alpha_incr_gradientmri, vol, opposite);

	// Regularization
	bool en = false;
	bool smooth = true;
	int transform = 3;
	double lambda = 0.0;
	double lambda_lim = 0.7;
	double lambda_incr = 400.00;
	lambda_incr = (lambda_lim - (double)lambda) / lambda_incr;
	cs->SetExternalForceRegularization(en, smooth, transform, lambda, lambda_lim, lambda_incr);

	double timestep = 0.1;
	double alpha = 1.0;
	Deformation *def = new Deformation(cs, timestep, alpha);
	
	for (int i = 0; i <= 200; i++) {
		cout << "Iteration: " << i << endl;
		def->ImpEuler();

		if (i % 10 == 0) {
			std::stringstream ss;
			ss << "DeformResults\\DeformedMesh_" << i << ".vtk";
			def->getMesh()->writeCSimplexMeshAsVTKPolyData(ss.str().c_str());
		}
	}
	*/

	return 0;
}

int main() {

	CSimplexSurf *a = new(CSimplexSurf);
	a->LoadMesh("C:\\Users\\Tanweer Rashid\\Downloads\\mesh1.defm");
	a->UpdateAll();

	a->writeCSimplexMeshAsVTKPolyData("C:\\Users\\Tanweer Rashid\\Desktop\\mesh1.vtk");


	//CSimplexSurf **a = Load_MultiMaterialCSimplexMeshFromVTKPolyData("data\\mesh[0]_Tmp_AfterSurfaceDecimation.vtk");
	//
	//MMCSimplexSurf *aa = new MMCSimplexSurf(a, 2);
	//aa->GetMesh(0)->writeCSimplexMeshAsVTKPolyData("read.vtk");

	std::string triangularFilename = "C:\\Users\\Tanweer Rashid\\Desktop\\DeformationResults_V2\\Left_SN_STN_Again\\Object_39-47_Left_padded_step0.9_L_39_47_V_0.7_0.7_m2b_0.4_m2m_0.3_MM2M_World_AffineRegistered2.vtk";
	vtkSmartPointer<vtkPolyData> mmSimplex = convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh(triangularFilename);
		
	CSimplexSurf **meshes = Split_MultiMaterialSimplexMesh_Into_Separate_CSimplexMesh(mmSimplex);

	MMCSimplexSurf *cm = new MMCSimplexSurf(meshes, 2);
	//cm->GetMesh(0)->writeCSimplexMeshAsVTKPolyData("data/mesh[0].vtk");
	//cm->GetMesh(1)->writeCSimplexMeshAsVTKPolyData("data/mesh[1].vtk");
	//cm->GetMesh(2)->writeCSimplexMeshAsVTKPolyData("data/mesh[2].vtk");
		

	//cm->GetMesh(0)->writeCSimplexMeshAsVTKPolyData("data/mesh[0].vtk");
	//cm->GetMesh(1)->writeCSimplexMeshAsVTKPolyData("data/mesh[1].vtk");
	//cm->GetMesh(2)->writeCSimplexMeshAsVTKPolyData("data/mesh[2].vtk");

	

	cm->DecimateByRemovingTriangularCells(0);
	cm->DecimateByRemovingTriangularCells(0);
	cm->GetMesh(0)->writeCSimplexMeshAsVTKPolyData("C:\\Users\\Tanweer Rashid\\Desktop\\DeformationResults_V2\\Left_SN_STN_Again\\Simplified.vtk");
	
	
		
	
	//cm->DecimateByMergingSmallSurfaceAreaCells(0, false);
	
	int i = 0;
	while (i < 300) {
		cout << "Iteration: " << i << " ";
		cm->DecimateBySubMesh(0, 0, 0);


		cm->GetMesh(0)->writeCSimplexMeshAsVTKPolyData("data\\mesh[0]_Tmp_AfterSurfaceDecimation.vtk");

		i = i + 1;
	}
	


	//cm->GetMesh(0)->TO1(21940, 21934);
	cm->GetMesh(0)->writeCSimplexMeshAsVTKPolyData("data\\mesh[0]_AfterSurfaceAreaDecimation.vtk");


	cout << "\n\nExecution completed. " << endl;
	getchar();
	return 0;
}




// Takes a CSimplexSurf object and visualizes it. 
void visualizeSimplex(CSimplexSurf *c) {	
	vtkSmartPointer<vtkPolyData> vp = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();	
	
	double *allPoints = c->GetPoints();
	
	for (int i = 0; i < (3 * c->GetNumberOfPoints()); i = i + 3) {
		points->InsertNextPoint(allPoints[i], allPoints[i + 1], allPoints[i + 2]);		
	}
	
	for (int j = 0; j < c->GetNumberOfCells(); j++) {
		int *q = c->GetCell(j);
		
		int df = 1;
		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		while (df <= q[0]) {
			cellPointIds->InsertNextId(q[df]);
			df = df + 1;
		}
		cells->InsertNextCell(cellPointIds);
	}

	vp->SetPoints(points);
	vp->SetPolys(cells);

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(vp);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();
	renwin->AddRenderer(ren);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renwin);

	ren->AddActor(actor);

	renwin->Render();
	iren->Start();
	
	return;
}


void boxSimplex() {
	CSimplexSurf *cs = new CSimplexSurf;
	cs->Free();
	cs->Allocate(8, 6);

	cs->SetPoint(0, 0, 0, 0);
	cs->SetPoint(1, 0, 0, 1);
	cs->SetPoint(2, 0, 1, 0);
	cs->SetPoint(3, 0, 1, 1);
	cs->SetPoint(4, 1, 0, 0);
	cs->SetPoint(5, 1, 0, 1);
	cs->SetPoint(6, 1, 1, 0);
	cs->SetPoint(7, 1, 1, 1);

	int f1[] = {0, 1, 5, 4};
	int f2[] = {2, 3, 7, 6};
	int f3[] = {7, 6, 4, 5};
	int f4[] = {3, 2, 0, 1};
	int f5[] = {1, 3, 7, 5};
	int f6[] = {0, 2, 6, 4};
	
	cs->SetCell(0, 4, f1);
	cs->SetCell(1, 4, f2);
	cs->SetCell(2, 4, f3);
	cs->SetCell(3, 4, f4);
	cs->SetCell(4, 4, f5);
	cs->SetCell(5, 4, f6);

	visualizeSimplex(cs);
}

void computeCellSurfaceAreaUsingVTKMassProperties(CSimplexSurf *cs) {
	ofstream f("C:\\Users\\Tanweer Rashid\\Desktop\\CSimplexCellSurfaceArea_comparison_against_vtkMassProperties.txt");
	f << "Cell surface area computed using vtkMassProperties; Cell surface area computed using Gilles' method\n";

	for (int i = 0; i < cs->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
			
		int *icell = cs->GetCell(i);

		double *center = cs->GetCellCenter(i);
		points->InsertNextPoint(center);
		
		for (int j = 1; j <= icell[0]; j++) {
			double *pt = cs->GetPoint(icell[j]);
			points->InsertNextPoint(pt);
		}

		for (int j = 1; j < icell[0]; j++) {
			vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
			list->InsertNextId(0);
			list->InsertNextId(j);
			list->InsertNextId(j + 1);

			cellArray->InsertNextCell(list);
		}

		vtkSmartPointer<vtkIdList> list1 = vtkSmartPointer<vtkIdList>::New();
		list1->InsertNextId(0);
		list1->InsertNextId(icell[0]);
		list1->InsertNextId(1);
		cellArray->InsertNextCell(list1);

		pdata->SetPoints(points);
		pdata->SetPolys(cellArray);

		vtkSmartPointer<vtkTriangleFilter> trifilter = vtkSmartPointer<vtkTriangleFilter>::New();
		trifilter->SetInputData(pdata);
		trifilter->Update();

		vtkSmartPointer<vtkMassProperties> massP = vtkSmartPointer<vtkMassProperties>::New();
		massP->SetInputData(trifilter->GetOutput());
		massP->Update();

		f << massP->GetSurfaceArea() << "; " << cs->getSurfaces_cell(i) << "\n";
	}
	
	f.close();
}

/*
void main() {
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName("C:\\Users\\Tanweer Rashid\\Desktop\\Object39_Decimation_Results\\CloseOpen_Smoothed\\Object39_close_open_2Manifold_Smoother_SimplexOriented.vtk");
	reader->Update();
	
	vtkSmartPointer<vtkPolyData> pdata = reader->GetOutput();
	pdata->BuildCells();
	pdata->BuildLinks();

	CSimplexSurf *cs = new(CSimplexSurf); 

	cs->Free();
	cs->Allocate(pdata->GetNumberOfPoints(), pdata->GetNumberOfCells());
	
	for (int i = 0; i < pdata->GetNumberOfPoints(); i++) {
		double *d = pdata->GetPoint(i);
		cs->SetPoint(i, d);
	}

	for (int i = 0; i < pdata->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkCell> cell = pdata->GetCell(i);

		int n = cell->GetNumberOfPoints();
		int *c = new int[n];
		
		//c[0] = n;
		for (int j = 0; j < n; j++) {
			c[j] = cell->GetPointId(j);

			cs->SetCell(i, n, c);
		}
	}

	
	//CSimplexSurf *cs = new(CSimplexSurf);
	//cs->LoadMesh("C:\\Users\\Tanweer Rashid\\Desktop\\Object26_part2_ConvertedToSimplex_Oriented.defm");
	cs->UpdateNeighbors();
	cs->UpdateAll();

	int numPts_original = cs->GetNumberOfPoints();
	int numCells_original = cs->GetNumberOfCells();

	cs->decimateByRemovingTriangularCells();
	//cs->writeCSimplexMeshAsVTKPolyData("C:\\Users\\Tanweer Rashid\\Desktop\\Object26_part2_ConvertedToSimplex_Oriented_TRI_Decimated.vtk");
	//cs->SaveMesh("C:\\Users\\Tanweer Rashid\\Desktop\\Object26_part2_ConvertedToSimplex_Oriented_TRI_Decimated.defm", true, false);

	//cs->decimateByMergingSmallSurfaceAreaCells(0.1);
	//cs->decimateByMergingSmallSurfaceAreaCells(0.5);
	//cs->decimateByMergingSmallSurfaceAreaCells(0.25);
	
	cs->writeCSimplexMeshAsVTKPolyData("C:\\Users\\Tanweer Rashid\\Desktop\\Object39_Decimation_Results\\CloseOpen_Smoothed\\Object39_close_open_2Manifold_Smoother_SimplexOriented_Decimated.vtk");
	cs->SaveMesh("C:\\Users\\Tanweer Rashid\\Desktop\\Object39_Decimation_Results\\CloseOpen_Smoothed\\Object39_close_open_2Manifold_Smoother_SimplexOriented_Decimated.defm", true, false);

	//visualizeSimplex(cs);
	
	int numPts_decimated = cs->GetNumberOfPoints();
	int numCells_decimated = cs->GetNumberOfCells();

	cout << "Points decimation % " << (100 * (numPts_original - numPts_decimated) / numPts_original) << endl;
	cout << "Cells decimation % " << (100 * (numCells_original - numCells_decimated) / numCells_original) << endl;


	//vtkSmartPointer<vtkPolyData> pdata = convertMultiMaterialTriangularMeshToMultiMaterialSimplexMesh("C:\\Users\\Tanweer Rashid\\Desktop\\Object39_Decimation_Results\\CloseOpen_Smoothed\\Object39_close_open_2Manifold_Smoother.vtk");
	//writeOrVisualizePolyData(pdata, true, true, "C:\\Users\\Tanweer Rashid\\Desktop\\Object39_Decimation_Results\\CloseOpen_Smoothed\\Object39_close_open_2Manifold_Smoother_SimplexOriented.vtk");

	cout << "Execution completed." << endl;
	getchar();
	
}
*/