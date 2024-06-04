//#include "StdAfx.h"
#include "SimplexMesh.h"
#include "math.h"
#include "float.h"
#include <ostream>
#include <sstream>

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
#include <vtkDoubleArray.h>
#include <vtkCharArray.h>
//#include "perftimer.h"

#include "Deformation.h"

const P_float MIN_RADIUS=1.5;
const bool LIFTING=true;
const bool USEMULTIPLEPARTICLES=false;
//const P_float SCALEINVARIANT_REFSHAPE=0.5;  // scale invariant
const P_float SCALEINVARIANT_REFSHAPE=0.5;  // 
//const P_float SCALEINVARIANT_REFSHAPE=0.25;  // 

CSimplexSurf::CSimplexSurf(void)
{ 
Allocate(1,1);
InternalForce_AxialConstraint=false;InternalForce_SurfEqu=false;InternalForce_SurfEquFlexion=false;InternalForce_Laplacian=false;InternalForce_LaplacianFlexion=false;InternalForce_RefShape=false;InternalForce_VolPres=false;InternalForce_SurfPres=false;ExternalForce_ICP=false;ExternalForce_GradientMRI=false;ExternalForce_Demon=false;ExternalForce_IntensityProfile=false;ExternalForce_ConstraintPoints=false;ExternalForce_ConstraintMesh=false;
SlideAttach=false;
RefModel=NULL;
RefModelPointCells=NULL;
RefModelCellNormals=NULL;
RefPointsLocator=NULL;
normals2=NULL;
MultiresMethod=0;
Topo=0; Topo_Decrease_L[0]=0; Topo_Increase_L[0]=0; Topo_Elongation[0]=0; Topo_Exchange[0]=0; Topo_ignoreborder=false;
//Topo_Decrease_R[0]=0; Topo_Increase_R[0]=0; Topo_Decrease_Err[0]=0; Topo_Increase_Err[0]=0; 
IntensityProfile_mn=0; IntensityProfile_pn=0; IntensityProfile_mn_ref=0; IntensityProfile_pn_ref=0;  IntensityProfile_s=1; 
Interpolationmode=0; //trilinear
Edgedivide_IntensityProfile=-1;
IsInternalPointsHandlingMethodCollision=true;
IsExternalPointsHandlingMethodCollision=true;
IsFrontierPointsHandlingMethodCollision=true;
IsAxis=false; IsSkin=false;
MergedPoints=NULL;
axiallinks=NULL;
Axis_model=NULL;
AxialConstraint_mode=0;
AxialLinks_UniformR=false;
Topo_Tol=1;
Forcehandleborder=false; 
Threshold_volsurfpres=10;
Metric_IntensityProfile=0;
Regularization_En=false; Regularization_Smooth=false; Regularization_Rigid=false; Regularization_Simi=false; Regularization_Affine=false; 
Regularization_Lambda=1; Regularization_Lambda_incr=0; Regularization_Lambda_lim=0;
IntensityProfile=NULL;			IntensityProfile_grad=NULL;
IntensityProfileRef=NULL;		IntensityProfileRef_grad=NULL;		IntensityProfile_Proc_Mask=NULL;
IntensityProfile_GaussianWeights=NULL; IntensityProfile_GaussianR=1;  IntensityProfile_GaussianS=0;	
Axial_smoothR=-1; Axial_smoothR_axis2=NULL; Axial_blockside=-1;
ExternalForce_ConstraintMesh_pts=NULL;
Monitor_axiallinks_error=false;
Monitor_surfaces=false;

Monitor_TO2=0; Monitor_TO2_test=0;	Monitor_TO1=0; Monitor_TO1_test=0;
Monitor_TOech=0; Monitor_TOech_test=0;	Monitor_TOech2=0; Monitor_TOech2_test=0;
Monitor_TOcell_i=0; Monitor_TOcell_i_test=0;	Monitor_TOcell_d=0; Monitor_TOcell_d_test=0;
Monitor_TOvertex_i=0; Monitor_TOvertex_i_test=0;	Monitor_TOvertex_d=0; Monitor_TOvertex_d_test=0;
Monitor_TOduration=0;

MRI_nb=0; MRI=NULL; MRI_grad=NULL; MRI_Transform=NULL; MRI_Spacing=NULL; MRI_Normals=NULL; MRI_OffsetTransform=NULL;	MRI_ROI=NULL; MRI_Type=NULL;
MRI_GradientMRI=NULL;  MRI_GradientMRI2=NULL;  oppositedirection_GradientMRI=false;

Elongations_ref=NULL;  Elongations_refindex=NULL;
}

CSimplexSurf::~CSimplexSurf(void)
{
Free();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::Free()
{
int i,j;
for(i=0;i<nb_cells;i++) free(cells[i]);
for(j=0;j<nb_points;j++)  {for(i=0;i<NB_NEIGHBORHOOD;i++) free(neighbors2[j][i]); free(neighbors2[j]); } free(neighbors2);
free(points); free(points_tm1);free(params);free(speeds);free(speeds_tm1);free(mass_inv);free(mass);free(forces);free(df_s);free(df_p);free(normals);free(neighbors);free(neighbors_c);free(attached_point);free(cellcenters);free(surfaces_cell);free(attached_cell);free(externalforces);free(externalforces_ign);free(constraintpointsforces);free(cells); free(surfaces); free(h); 
//free(volumes);
if(normals2!=NULL) free(normals2);
//free(flagmultires); 
free(Val_p); free(Val_c);
//free(params_multires);
for(j=0;j<nb_points;j++)  if(mass_constraint[j]!=NULL) free(mass_constraint[j]);
if(IntensityProfileRef!=NULL)			{for(j=0;j<nb_points;j++)  free(IntensityProfileRef[j]);			free(IntensityProfileRef);} 
if(IntensityProfileRef_grad!=NULL)		{for(j=0;j<nb_points;j++)  free(IntensityProfileRef_grad[j]);			free(IntensityProfileRef_grad);} 
if(IntensityProfile!=NULL)				{for(j=0;j<nb_points;j++)  free(IntensityProfile[j]);				free(IntensityProfile);}
if(IntensityProfile_grad!=NULL)			{for(j=0;j<nb_points;j++)  free(IntensityProfile_grad[j]);			free(IntensityProfile_grad);}
if(IntensityProfile_GaussianWeights!=NULL) {for(j=0;j<IntensityProfile_GaussianR;j++) free(IntensityProfile_GaussianWeights[j]); free(IntensityProfile_GaussianWeights);}
if(IntensityProfile_Proc_Mask!=NULL) free(IntensityProfile_Proc_Mask);
free(mass_constraint);
if(RefModelPointCells!=NULL) {for(i=0;i<RefModel->GetNumberOfPoints();i++) RefModelPointCells[i]->Delete(); free(RefModelPointCells);}
if(RefModelCellNormals!=NULL) free(RefModelCellNormals);
if(RefModel!=NULL) RefModel->Delete();
if(RefPointsLocator!=NULL) RefPointsLocator->Delete();
if(InternalConstraintPoints!=NULL) {InternalConstraintPoints->Delete(); InternalConstraintPoints=NULL;} 
if(ExternalConstraintPoints!=NULL) {ExternalConstraintPoints->Delete(); ExternalConstraintPoints=NULL;} 
if(FrontierConstraintPoints!=NULL) {FrontierConstraintPoints->Delete(); FrontierConstraintPoints=NULL;} 
if(ExternalForce_ConstraintMesh_pts!=NULL) free(ExternalForce_ConstraintMesh_pts);
if(axiallinks!=NULL) {for(i=0;i<nb_points;i++) {if(axiallinks[i].w!=NULL) free(axiallinks[i].w); if(axiallinks[i].pts!=NULL) free(axiallinks[i].pts);} free(axiallinks);}

DOPCellFree(DOPCellTree);
free(DOPCellTree);
DOPPointFree(DOPPointTree);
free(DOPPointTree);

for(i=0;i<MRI_nb;i++) FreeMRI(i);
if(MRI_nb!=0) {free(MRI); free(MRI_grad); free(MRI_Transform); free(MRI_Spacing); free(MRI_Normals); free(MRI_OffsetTransform);	free(MRI_ROI); free(MRI_Type);}
if(Elongations_ref!=NULL) free(Elongations_ref); if(Elongations_refindex!=NULL) free(Elongations_refindex);
}

void CSimplexSurf::Allocate(int nb_p,int nb_c)
{
int i,j; nb_points=nb_p; nb_cells=nb_c;
points=new P_float[3*nb_points]; 
speeds=new P_float[3*nb_points];
points_tm1=new P_float[3*nb_points]; 
speeds_tm1=new P_float[3*nb_points];
df_p=new P_float[6*nb_points]; 
df_s=new P_float[6*nb_points]; 
forces=new P_float[3*nb_points]; 
Val_p=new P_float[3*nb_points]; Val_c=new P_float[3*nb_cells]; 
mass_inv=new P_float[3*nb_points];  for(i=0;i<3*nb_points;i++) mass_inv[i]=1;  
mass=new P_float[3*nb_points];  for(i=0;i<3*nb_points;i++) mass[i]=1;  
mass_constraint=new P_float*[nb_points];  for(i=0;i<nb_points;i++) mass_constraint[i]=NULL;  
params=new P_float[3*nb_points]; 
//params_multires=new P_float[3*nb_points];
normals=new P_float[3*nb_points]; 
surfaces=new P_float[nb_points]; 
//volumes=new P_float[nb_points]; 
h=new P_float[nb_points]; 
externalforces=new P_float[3*nb_points]; 
externalforces_ign=new bool[nb_points]; 
constraintpointsforces=new P_float[4*nb_points];
neighbors=new int[3*nb_points]; for(i=0;i<3*nb_points;i++) neighbors[i]=-1;
neighbors_c=new int[3*nb_points]; 

//// TR mods
//for(i=0;i<3*nb_points;i++) neighbors_c[i]=-1;
////

attached_point=new int[nb_points]; for(i=0;i<nb_points;i++) attached_point[i]=0;
attached_cell=new int[nb_cells]; for(i=0;i<nb_cells;i++) attached_cell[i]=0;
cellcenters=new P_float[3*nb_cells]; 
surfaces_cell=new P_float[nb_cells]; 
//flagmultires=new int[nb_cells]; for(i=0;i<nb_cells;i++) flagmultires[i]=0;
cells=(int**)malloc(nb_cells*sizeof(int*)); for(i=0;i<nb_cells;i++) {cells[i]=new int[1]; cells[i][0]=0;}
neighbors2=(int***)malloc(nb_points*sizeof(int**)); for(j=0;j<nb_points;j++) {neighbors2[j]=(int**)malloc(NB_NEIGHBORHOOD*sizeof(int*)); for(i=0;i<NB_NEIGHBORHOOD;i++) {neighbors2[j][i]=new int[1]; neighbors2[j][i][0]=0;}}
UGrid=vtkUnstructuredGrid::New();
UGrid->Allocate(nb_points,nb_cells);
InternalConstraintPoints=vtkPoints::New(); ExternalConstraintPoints=vtkPoints::New();  FrontierConstraintPoints=vtkPoints::New();  InternalConstraintPoints->SetDataTypeToFloat(); ExternalConstraintPoints->SetDataTypeToFloat(); FrontierConstraintPoints->SetDataTypeToFloat();

DOPCellTree=new DOPCELLSTRUCT;
DOPCellTree->nb_cells=1;
DOPCellTree->cells=new int[1]; *DOPCellTree->cells=0;
DOPCellTree->dopm=NULL; DOPCellTree->dopn=NULL; DOPCellTree->dopo=NULL; DOPCellTree->dopp=NULL;

DOPPointTree=new DOPPOINTSTRUCT;
DOPPointTree->nb_points=1;
DOPPointTree->points=new int[1]; *DOPPointTree->points=0;
DOPPointTree->dopm=NULL; DOPPointTree->dopn=NULL; DOPPointTree->dopo=NULL; DOPPointTree->dopp=NULL;

AllocatePointMaterialIndices();
AllocateIsMultiMaterialPoints();

}


void CSimplexSurf::Copy(CSimplexSurf* mesh)
{
Free();
int i,nb_p=mesh->GetNumberOfPoints(),nb_c=mesh->GetNumberOfCells();
P_float p[3];

Allocate(nb_p,nb_c);
for(i=0;i<nb_p;i++) {mesh->GetPoint(i,p); SetPoint(i,p);}
for(i=0;i<nb_c;i++) {SetCell(i,mesh->GetCell(i)); SetAttachCell(i,mesh->GetAttachCell(i));}
//int nb=mesh->GetFlagMultires(i); SetFlagMultires(i,nb); 
for(i=0;i<nb_p;i++) if(mesh->GetAttachPoint(i,p)!=0) {SetAttachPoint(i,p); attached_point[i]=mesh->GetAttachPoint(i,p);}
UpdateNeighbors();
Equilibrium();
}

bool CSimplexSurf::GetIsAxis(){return IsAxis;} 	
void CSimplexSurf::SetIsAxis(bool isaxis) 
{
IsAxis=isaxis;
if(IsAxis) {if(normals2!=NULL) free(normals2); normals2=new P_float[3*nb_points]; }
}
bool CSimplexSurf::GetIsSkin(){return IsSkin;} 	
void CSimplexSurf::SetIsSkin(bool isskin) {IsSkin=isskin;}



void CSimplexSurf::SetMultiresMethod(int val) {MultiresMethod=val;}

CSimplexSurf* CSimplexSurf::IncreaseResolution() {
	MultiresMethod = 1;
	if (MultiresMethod == 0) {
		/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2    
		int i, j;
		P_float *c, p[3], p1[3], alpha, nm[3], pm[3], n[3];
		int pt_index = 0, cell_index = 0, cell[200], pt1, pt2;

		int nb_p = 4 * nb_points;
		int nb_c = nb_cells + 3 * nb_points / 2;
		int* int_cells = new int[12 * nb_points]; 
		for (i = 0; i < 12 * nb_points; i++) 
			int_cells[i] = -1; 

		CSimplexSurf* mesh = new CSimplexSurf;
		mesh->Free();
		mesh->Allocate(nb_p, nb_c);
		mesh->SetIsAxis(IsAxis);

		for (i = 0; i < nb_points; i++) {
			mesh->SetPoint(i, points + 3 * i); 
			if (GetMass_inv(i)[0] == 0) 
				mesh->SetMass(i, -1);
		}
	
		pt_index = nb_points;

		for (i = 0; i < nb_cells; i++) { 
			// cell center
			c = cellcenters + 3 * i;

			// LIFTING
			if (LIFTING) {
				nm[0] = 0; nm[1] = 0; nm[2] = 0;
				for (j = 0; j < cells[i][0]; j++) {
					nm[0] += normals[3 * cells[i][j + 1]]; 
					nm[1] += normals[3 * cells[i][j + 1] + 1];
					nm[2] += normals[3 * cells[i][j + 1] + 2];
				}
				alpha = norm(nm);	
				nm[0] = nm[0] / alpha; 
				nm[1] = nm[1] / alpha;
				nm[2] = nm[2] / alpha;
			}

			// compute points (CPi centers)
			for (j = 0; j < cells[i][0]; j++) {
				pt1 = cells[i][j + 1]; 
				GetPoint(pt1, p1); 

				pt2 = cells[i][(j == cells[i][0] - 1) ? (1) : (j + 2)];
				
				p[0] = (p1[0] + c[0]) / 2.;    
				p[1] = (p1[1] + c[1]) / 2.;    
				p[2] = (p1[2] + c[2]) / 2.;
				
				mesh->SetPoint(pt_index, p); 
				cell[j] = pt_index; 
        
				// LIFTING
				if (LIFTING) {
					pm[0] = p1[0] + normals[3 * pt1] * h[pt1] / 4.; 
					pm[1] = p1[1] + normals[3 * pt1 + 1] * h[pt1] / 4.;
					pm[2] = p1[2] + normals[3 * pt1 + 2] * h[pt1] / 4.;

					n[0] = (nm[0] + normals[3 * pt1]) / 2.;
					n[1] = (nm[1] + normals[3 * pt1 + 1]) / 2.; 
					n[2] = (nm[2] + normals[3 * pt1 + 2]) / 2.; 

					alpha = norm(n); n[0] = n[0] / alpha; 
					n[1] = n[1] / alpha;
					n[2] = n[2] / alpha; 
					//		GetNormal(pt1,n);
					alpha = (dotproduct(pm, normals + 3 * pt1) - dotproduct(p, normals + 3 * pt1)) / dotproduct(n, normals + 3 * pt1); 
					//		alpha=-h[pt1]/2.;
					p[0] += alpha * n[0];	p[1] += alpha * n[1]; p[2] += alpha * n[2];
					mesh->SetPoint(pt_index, p);
				}
		
				// int cells[edgeindex] :  pt2 i2 i1 pt1
				int_cells[4 * (pt_index - nb_points)] = pt2;
				int_cells[4 * (pt_index - nb_points) + 1] = (j == cells[i][0] - 1) ? (pt_index - cells[i][0] + 1) : (pt_index + 1);
				int_cells[4 * (pt_index - nb_points) + 2] = pt_index;
				int_cells[4 * (pt_index - nb_points) + 3] = pt1;

				pt_index++;
			}
			
			//mesh->SetFlagMultires(cell_index,GetFlagMultires(i)+1);
			mesh->SetCell(cell_index, cells[i][0], cell); 
			cell_index++;
		}



		for (i = 0; i < 3 * nb_points; i++) {
			for (j = 0; j < 3 * nb_points; j++)
				if (int_cells[4 * i + 3] == int_cells[4 * j] && int_cells[4 * i] == int_cells[4 * j + 3] && int_cells[4 * j] != -1 && int_cells[4 * i] != -1) {
					memcpy(cell, int_cells + 4 * i, 4 * sizeof(int));
					memcpy(cell + 4, int_cells + 4 * j + 1, 2 * sizeof(int));
					mesh->SetCell(cell_index, 6, cell); cell_index++;
					int_cells[4 * i + 3] = -1;
					int_cells[4 * j] = -1;  
				}
		}
		mesh->SetMultiresMethod(0);

		double dp[3];
		for (i = 0; i < InternalConstraintPoints->GetNumberOfPoints(); i++) {
			InternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertInternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < ExternalConstraintPoints->GetNumberOfPoints(); i++) {
			ExternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertExternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < FrontierConstraintPoints->GetNumberOfPoints(); i++) {
			FrontierConstraintPoints->GetPoint(i, dp);
			mesh->InsertFrontierConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}

		// update
		mesh->UpdateNeighbors();
		mesh->Equilibrium();
		free(int_cells);
		return mesh;
	}

	else {
		/////////// METHODE 1 : nb_p=3*nb_points nb_c=nb_cells+nb_points
		int i, j;
		P_float *c, p[3], p1[3], p2[3];
		int pt_index = 0, cell_index = 0, int_cells_index, cell[20], pt1, pt2;

		int nb_p = 0; 
		for (i = 0; i < nb_cells; i++) 
			nb_p += cells[i][0];
		
		int nb_c = nb_cells + nb_points;
		int* int_cells = new int[6 * nb_points]; 

		CSimplexSurf* mesh = new CSimplexSurf;
		mesh->Free();
		mesh->Allocate(nb_p, nb_c);
		mesh->SetIsAxis(IsAxis); 
		mesh->SetIsSkin(IsSkin);

		for (i = 0; i < nb_cells; i++) { 
			// cell center
			c = cellcenters + 3 * i;

			// compute points (triangle centers)
			for (j = 0; j < cells[i][0]; j++) {
				pt1 = cells[i][j + 1];
				pt2 = cells[i][(j == cells[i][0] - 1) ? (1) : (j + 2)];
				
				GetPoint(pt1, p1); 
				GetPoint(pt2, p2);
				
				p[0] = (p1[0] + p2[0] + c[0]) / 3; 
				p[1] = (p1[1] + p2[1] + c[1]) / 3;  
				p[2] = (p1[2] + p2[2] + c[2]) / 3;
				
				mesh->SetPoint(pt_index, p);
				cell[j] = pt_index; 
        
				// int cells[i] :  neighbor3(pt2,pt1),  neighbor2(pt2,pt1),  neighbor1(pt2,pt1)
				
				
				//int_cells_index = (neighbors[3 * pt1] == pt2) ? 5 : 0 + (neighbors[3 * pt1 + 1] == pt2) ? 3 : 0 + (neighbors[3 * pt1 + 2] == pt2) ? 1 : 0; 
				int *n1 = GetNeighbors(pt1);
				int_cells_index = (n1[0] == pt2) ? 5 : 0 + (n1[1] == pt2) ? 3 : 0 + (n1[2] == pt2) ? 1 : 0; 
				int_cells[6 * pt1 + int_cells_index] = pt_index;

				//int_cells_index = (neighbors[3 * pt2] == pt1) ? 4 : 0 + (neighbors[3 * pt2 + 1] == pt1) ? 2 : 0; 
				int *n2 = GetNeighbors(pt2);				
				int_cells_index = (n2[0] == pt1) ? 4 : 0 + (n2[1] == pt1) ? 2 : 0; 
				int_cells[6 * pt2 + int_cells_index] = pt_index;
				
				pt_index++;
			}
			//mesh->SetFlagMultires(cell_index,GetFlagMultires(i)+1);
			mesh->SetCell(cell_index, cells[i][0], cell); 
			cell_index++;
		}

		for (i = 0; i < nb_points; i++) {
			memcpy(cell, int_cells + 6 * i, 6 * sizeof(int));
			mesh->SetCell(cell_index, 6, cell); cell_index++;
		}
		mesh->SetMultiresMethod(1);

		double dp[3];
		for (i = 0; i < InternalConstraintPoints->GetNumberOfPoints(); i++) {
			InternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertInternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < ExternalConstraintPoints->GetNumberOfPoints(); i++) {
			ExternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertExternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < FrontierConstraintPoints->GetNumberOfPoints(); i++) {
			FrontierConstraintPoints->GetPoint(i, dp);
			mesh->InsertFrontierConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}

		// update
		mesh->UpdateNeighbors();
		mesh->Equilibrium();
		free(int_cells);
		return mesh;
	}
}




CSimplexSurf* CSimplexSurf::DecreaseResolution() {
	if (MultiresMethod == 0) {
		/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
		int i, j;
		int pt1, pt2, pt3, cell[100], n[3];

		int nb_c = 0, nb_p = nb_points / 4;  nb_c = nb_cells - 3 * nb_p / 2; 
		
		//for(i=0;i<nb_cells;i++) if(flagmultires[i]!=0) nb_c++;
		// test decreasability
		if (nb_c == 0) return NULL;

		CSimplexSurf* mesh = new CSimplexSurf;
		mesh->Free();
		mesh->Allocate(nb_p, nb_c);
		mesh->SetIsAxis(IsAxis); 
		mesh->SetIsSkin(IsSkin);

		for (i = 0; i < nb_c; i++) {
			for (j = 0; j < cells[i][0]; j++) {
				pt1 = cells[i][j + 1];  
				pt2 = cells[i][(j == cells[i][0] - 1) ? (1) : (j + 2)]; 
				pt3 = cells[i][(j == 0) ? (cells[i][0]) : j];
				
				GetNeighbors(pt1, n);
				if (n[0] != pt2 && n[0] != pt3)  cell[j] = n[0];
				if (n[1] != pt2 && n[1] != pt3)  cell[j] = n[1];
				if (n[2] != pt2 && n[2] != pt3)  cell[j] = n[2];
			}
			//mesh->SetFlagMultires(i,GetFlagMultires(i)-1);
			mesh->SetCell(i, cells[i][0], cell);
		}

		for (i = 0; i < nb_p; i++) 
			mesh->SetPoint(i, points + 3 * i);

		mesh->SetMultiresMethod(0);

		double dp[3];
		for (i = 0; i < InternalConstraintPoints->GetNumberOfPoints(); i++) {
			InternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertInternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < ExternalConstraintPoints->GetNumberOfPoints(); i++) {
			ExternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertExternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < FrontierConstraintPoints->GetNumberOfPoints(); i++){
			FrontierConstraintPoints->GetPoint(i, dp);
			mesh->InsertFrontierConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}

		// update
		mesh->UpdateNeighbors();
		mesh->Equilibrium();
		return mesh;
	}

	else {
		/////////// METHODE 1 : nb_p=3*nb_points nb_c=nb_cells+nb_points
		int i, j, i1, i2;
		int pt_index = 0, cell_index = 0, index_int_cells, cell[10];

		int nb_c = 0; int nb_p = 0; 
		//for(i=0;i<nb_cells;i++) if(flagmultires[i]!=0) nb_c++; else nb_p++;

		// test decreasability
		if (nb_c == 0) return NULL;

		CSimplexSurf* mesh = new CSimplexSurf;
		mesh->Free();
		mesh->Allocate(nb_p, nb_c);
		mesh->SetIsAxis(IsAxis); 
		mesh->SetIsSkin(IsSkin);

		int* int_cells = new int[2 * nb_points]; 
		for (i = 0; i < 2 * nb_points; i++) 
			int_cells[i] = -1;

		for (i = 0; i < nb_cells; i++) { 

			//if(flagmultires[i]==0) // the last created cells -> cell center=point
			for (j = 0; j < cells[i][0]; j++)  {
		
				// index last created cells number into corresponding points
				if (int_cells[2 * cells[i][j + 1]] == -1) 
					index_int_cells = 2 * cells[i][j + 1]; 
				else 
					index_int_cells = 2 * cells[i][j + 1] + 1; 
				
				int_cells[index_int_cells] = pt_index;
			}
			// set point
			mesh->SetPoint(pt_index, cellcenters + 3 * i); 
			pt_index++;
		}

		// create cells
		for (i = 0; i < nb_cells; i++) {
			//if(flagmultires[i]!=0) 
		
			for (j = 0; j < cells[i][0]; j++) {
				i1 = 2 * cells[i][j + 1];
				i2 = 2 * cells[i][(j == cells[i][0] - 1) ? 1 : (j + 2)];
				
				if (int_cells[i1] == int_cells[i2] || int_cells[i1] == int_cells[i2 + 1]) 
					cell[j] = int_cells[i1 + 1];
				else 
					cell[j] = int_cells[i1];
			}
			//mesh->SetFlagMultires(cell_index,GetFlagMultires(i)-1);
			mesh->SetCell(cell_index, cells[i][0], cell); 
			cell_index++;
		}
		mesh->SetMultiresMethod(1);

		double dp[3];
		for (i = 0; i < InternalConstraintPoints->GetNumberOfPoints(); i++) {
			InternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertInternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < ExternalConstraintPoints->GetNumberOfPoints(); i++) {
			ExternalConstraintPoints->GetPoint(i, dp);
			mesh->InsertExternalConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}
		for (i = 0; i < FrontierConstraintPoints->GetNumberOfPoints(); i++) {
			FrontierConstraintPoints->GetPoint(i, dp);
			mesh->InsertFrontierConstraintPoint((P_float)dp[0], (P_float)dp[1], (P_float)dp[2]);
		}

		// update
		mesh->UpdateNeighbors();
		mesh->Equilibrium();
		free(int_cells);
		return mesh;
	}
}


void CSimplexSurf::UpdateValsFromLowerRes(CSimplexSurf* mesh)
{
if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
int *c,i,j,pt_index=0;
P_float *p1,*ct;
pt_index=mesh->GetNumberOfPoints();

memcpy(Val_p,mesh->GetVal_p(0),pt_index*3*sizeof(P_float));

if(pt_index!=nb_points)
	{
	for(i=0;i<mesh->GetNumberOfCells();i++)
		{
		c=mesh->GetCell(i);    ct=mesh->GetVal_c(i);  
		for(j=0;j<c[0];j++)
			{
			p1=mesh->GetVal_p(c[j+1]); 
			Val_p[3*pt_index]=(p1[0]+ct[0])/2;  Val_p[3*pt_index+1]=(p1[1]+ct[1])/2;  Val_p[3*pt_index+2]=(p1[2]+ct[2])/2;
			pt_index++;
			}
		SetVal_c(i,ct);
		}
	for(i=mesh->GetNumberOfCells();i<GetNumberOfCells();i++)
		{
		c=GetCell(i); Val_c[3*i]=0;Val_c[3*i+1]=0;Val_c[3*i+2]=0;
		for(j=0;j<c[0];j++)
			{
			p1=GetVal_p(c[j+1]); 
			Val_c[3*i]+=p1[0]; Val_c[3*i+1]+=p1[1]; Val_c[3*i+2]+=p1[2];
			}
		Val_c[3*i]=Val_c[3*i]/c[0]; Val_c[3*i+1]=Val_c[3*i+1]/c[0]; Val_c[3*i+2]=Val_c[3*i+2]/c[0];
		}
	}

}
}





void CSimplexSurf::UpdateForcesFromLowerRes(CSimplexSurf* mesh,bool updatederivatives)
{
if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
int *c,i,j,pt_index=0;
P_float *f,fm[3],tmp[6],*dfp,dfpm[6];
pt_index=mesh->GetNumberOfPoints();

if(updatederivatives)
	{
	for(j=0;j<pt_index;j++)     {AddForce(j,mesh->GetForce(j)); AddDfp(j,mesh->GetDfp(j));}
	for(i=0;i<mesh->GetNumberOfCells();i++)
		{
		c=mesh->GetCell(i);    
		fm[0]=0;    fm[1]=0;    fm[2]=0;	dfpm[0]=0;  dfpm[1]=0;  dfpm[2]=0;  dfpm[3]=0;  dfpm[4]=0;  dfpm[5]=0; for(j=0;j<c[0];j++) {f=mesh->GetForce(c[j+1]); fm[0]+=f[0];   fm[1]+=f[1];    fm[2]+=f[2]; dfp=mesh->GetDfp(c[j+1]); dfpm[0]+=dfp[0]; dfpm[1]+=dfp[1];  dfpm[2]+=dfp[2];    dfpm[3]+=dfp[3];    dfpm[4]+=dfp[4];    dfpm[5]+=dfp[5];  }	fm[0]=fm[0]/c[0]; fm[1]=fm[1]/c[0]; fm[2]=fm[2]/c[0]; dfpm[0]=dfpm[0]/c[0]; dfpm[1]=dfpm[1]/c[0]; dfpm[2]=dfpm[2]/c[0];  dfpm[3]=dfpm[3]/c[0];    dfpm[4]=dfpm[4]/c[0]; dfpm[5]=dfpm[5]/c[0];
		for(j=0;j<c[0];j++) {mesh->GetForce(c[j+1],tmp); tmp[0]=(tmp[0]+fm[0])/2;    tmp[1]=(tmp[1]+fm[1])/2;    tmp[2]=(tmp[2]+fm[2])/2; AddForce(pt_index,tmp); mesh->GetDfp(c[j+1],tmp); tmp[0]=(tmp[0]+dfpm[0])/2; tmp[1]=(tmp[1]+dfpm[1])/2;  tmp[2]=(tmp[2]+dfpm[2])/2; tmp[3]=(tmp[3]+dfpm[3])/2;   tmp[4]=(tmp[4]+dfpm[4])/2;  tmp[5]=(tmp[5]+dfpm[5])/2;  AddDfp(pt_index,tmp); pt_index++; }
		}
	}
else
	{
	for(j=0;j<pt_index;j++)    AddForce(j,mesh->GetForce(j));
	for(i=0;i<mesh->GetNumberOfCells();i++)
		{
		c=mesh->GetCell(i);    
		fm[0]=0;    fm[1]=0;    fm[2]=0;	for(j=0;j<c[0];j++) {f=mesh->GetForce(c[j+1]); fm[0]+=f[0];   fm[1]+=f[1];    fm[2]+=f[2]; }	fm[0]=fm[0]/c[0]; fm[1]=fm[1]/c[0]; fm[2]=fm[2]/c[0];
		for(j=0;j<c[0];j++) {mesh->GetForce(c[j+1],tmp); tmp[0]=(tmp[0]+fm[0])/2;    tmp[1]=(tmp[1]+fm[1])/2;    tmp[2]=(tmp[2]+fm[2])/2; AddForce(pt_index,tmp); pt_index++; }
		}
	}
}
}
		
void CSimplexSurf::UpdatePointsFromLowerRes(CSimplexSurf* mesh)
{

if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
int *c,i,j,pt_index=0;
P_float pt[3],p1[3],ct[3],*pts;
pt_index=mesh->GetNumberOfPoints();
pts=mesh->GetPoints(); 
for(i=0;i<pt_index;i++)	if(mass_inv[3*i]!=0) memcpy(points+3*i,pts+3*i,3*sizeof(P_float));

if(pt_index!=nb_points)
	for(i=0;i<mesh->GetNumberOfCells();i++)
		{
		c=mesh->GetCell(i);
		mesh->GetCellCenter(i,ct);  
		for(j=0;j<c[0];j++)
			{
			mesh->GetPoint(c[j+1],p1); 
			pt[0]=(p1[0]+ct[0])/2;  pt[1]=(p1[1]+ct[1])/2;  pt[2]=(p1[2]+ct[2])/2;
			if(mass_inv[3*pt_index]!=0) SetPoint(pt_index,pt); 
			pt_index++; 
			}
		}
UpdateAll();
}

else {
/////////// METHODE 1 : nb_p=3*nb_points nb_c=nb_cells+nb_points
int c[200],i,j,pt_index=0,nb;
P_float pt[3],p1[3],p2[3],ct[3];
for(i=0;i<mesh->GetNumberOfCells();i++)
    {
    nb=mesh->GetCell(i,c);
    // cell center
    mesh->GetCellCenter(i,ct);  
    // compute points (triangle centers) 
    for(j=0;j<nb;j++)
        {
        mesh->GetPoint(c[j],p1); if(j==(nb-1)) mesh->GetPoint(c[0],p2); else mesh->GetPoint(c[j+1],p2);
        pt[0]=(p1[0]+p2[0]+ct[0])/3;    pt[1]=(p1[1]+p2[1]+ct[1])/3;    pt[2]=(p1[2]+p2[2]+ct[2])/3;
        SetPoint(pt_index,pt); pt_index++; 
        }
    }
UpdateAll();
}

}


void CSimplexSurf::UpdatePointsFromLowerRes_Disp(CSimplexSurf* mesh) // interpolate the displacement
{
if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
int *c,i,j,pt_index=0;
P_float dp[3],dm[3];
pt_index=mesh->GetNumberOfPoints();
for(i=0;i<mesh->GetNumberOfCells();i++)
	{
	c=mesh->GetCell(i);    
	dm[0]=0;    dm[1]=0;    dm[2]=0;	for(j=0;j<c[0];j++) {dp[0]=mesh->GetPoint(c[j+1])[0]-GetPoint(c[j+1])[0]; dp[1]=mesh->GetPoint(c[j+1])[1]-GetPoint(c[j+1])[1]; dp[2]=mesh->GetPoint(c[j+1])[2]-GetPoint(c[j+1])[2]; dm[0]+=dp[0];   dm[1]+=dp[1];    dm[2]+=dp[2]; }	dm[0]=dm[0]/c[0]; dm[1]=dm[1]/c[0]; dm[2]=dm[2]/c[0];
	for(j=0;j<c[0];j++) 
		{
		dp[0]=mesh->GetPoint(c[j+1])[0]-GetPoint(c[j+1])[0]; dp[1]=mesh->GetPoint(c[j+1])[1]-GetPoint(c[j+1])[1]; dp[2]=mesh->GetPoint(c[j+1])[2]-GetPoint(c[j+1])[2]; dp[0]=(dp[0]+dm[0])/2;    dp[1]=(dp[1]+dm[1])/2;    dp[2]=(dp[2]+dm[2])/2; 
		if(GetAttachedPoint(pt_index)==0) 
			{
			GetPoint(pt_index)[0]+=dp[0]; GetPoint(pt_index)[1]+=dp[1]; GetPoint(pt_index)[2]+=dp[2];
			if(axiallinks!=NULL) {axiallinks[pt_index].f[0]-=dp[0]; axiallinks[pt_index].f[1]-=dp[1]; axiallinks[pt_index].f[2]-=dp[2];}
			externalforces[3*pt_index]-=dp[0]; externalforces[3*pt_index+1]-=dp[1]; externalforces[3*pt_index+2]-=dp[2];
			constraintpointsforces[4*pt_index+1]-=dp[0]; constraintpointsforces[4*pt_index+2]-=dp[1]; constraintpointsforces[4*pt_index+3]-=dp[2];
			} 
		pt_index++; 
		}
	}
for(i=0;i<mesh->GetNumberOfPoints();i++)  
	if(GetAttachedPoint(i)==0) 
		{
		dp[0]=mesh->GetPoint(i)[0]-GetPoint(i)[0]; dp[1]=mesh->GetPoint(i)[1]-GetPoint(i)[1]; dp[2]=mesh->GetPoint(i)[2]-GetPoint(i)[2];
		GetPoint(i)[0]+=dp[0]; GetPoint(i)[1]+=dp[1]; GetPoint(i)[2]+=dp[2];
		if(axiallinks!=NULL) {axiallinks[i].f[0]-=dp[0]; axiallinks[i].f[1]-=dp[1]; axiallinks[i].f[2]-=dp[2];}
		externalforces[3*i]-=dp[0]; externalforces[3*i+1]-=dp[1]; externalforces[3*i+2]-=dp[2];
		constraintpointsforces[4*i+1]-=dp[0]; constraintpointsforces[4*i+2]-=dp[1]; constraintpointsforces[4*i+3]-=dp[2];
		}
UpdateAll();
}
}


void CSimplexSurf::UpdateValsFromHigherRes(CSimplexSurf* mesh)
{
if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
P_float* pt;
pt=mesh->GetVal_p(0); memcpy(Val_p,pt,nb_points*3*sizeof(P_float));
pt=mesh->GetVal_c(0); memcpy(Val_c,pt,nb_cells*3*sizeof(P_float));
}
}

void CSimplexSurf::UpdatePointsFromHigherRes(CSimplexSurf* mesh)
{
/*if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
P_float* pt;
pt=mesh->GetPoints(); memcpy(points,pt,nb_points*3*sizeof(P_float));
pt=mesh->GetPoints_tm1(); memcpy(points_tm1,pt,nb_points*3*sizeof(P_float));
pt=mesh->GetSpeeds(); memcpy(speeds,pt,nb_points*3*sizeof(P_float));
pt=mesh->GetSpeeds_tm1(); memcpy(speeds_tm1,pt,nb_points*3*sizeof(P_float));
UpdateAll();
}

else {
/////////// METHODE 1 : nb_p=3*nb_points nb_c=nb_cells+nb_points
int i,pt_index=0;
P_float pt2[3];
for(i=0;i<mesh->GetNumberOfCells();i++)
    if(mesh->GetFlagMultires(i)==0)
        {
        mesh->GetCellCenter(i,pt2);     
        SetPoint(pt_index,pt2); 
        SetSpeed(pt_index,0,0,0);
        SetForce(pt_index,0,0,0);
        pt_index++;
        }
UpdateAll();
}*/
}

void CSimplexSurf::UpdateVal_c() 
{
int *c,i,j; P_float *p1;
for(i=0;i<GetNumberOfCells();i++)
	{
	c=GetCell(i); Val_c[3*i]=0;Val_c[3*i+1]=0;Val_c[3*i+2]=0;
	for(j=0;j<c[0];j++)
        {
		p1=GetVal_p(c[j+1]); 
		Val_c[3*i]+=p1[0]; Val_c[3*i+1]+=p1[1]; Val_c[3*i+2]+=p1[2];
		}
	Val_c[3*i]=Val_c[3*i]/c[0]; Val_c[3*i+1]=Val_c[3*i+1]/c[0]; Val_c[3*i+2]=Val_c[3*i+2]/c[0];
	}
}


void CSimplexSurf::UpdateMergedPointsFromHigherRes(CSimplexSurf* mesh)
{
int i,nb=0,*mp=mesh->GetMergedPoints();
if(MultiresMethod==0)
	{
	for(i=0;i<mp[0];i++) if(mp[2*i+1]<nb_points && mp[2*i+2]<nb_points) nb++;
	if(nb!=0)
		{
		MergedPoints=new int[2*nb+1]; MergedPoints[0]=nb; nb=0;
		for(i=0;i<mp[0];i++) if(mp[2*i+1]<nb_points && mp[2*i+2]<nb_points) {MergedPoints[2*nb+1]=mp[2*i+1]; MergedPoints[2*nb+2]=mp[2*i+2]; nb++;}
		}
	}
}

/*
void CSimplexSurf::UpdatePointsFromHigherRes_Bound(CSimplexSurf* mesh)
{

if(MultiresMethod==0)
{
/////////// METHODE 0 : nb_p=4*nb_points nb_c=nb_cells+3*nb_points/2
int i,j,k,l,pt_index,i0,i1;
P_float W,ColResponse[3],dep[3],c0[3],c1[3],p[3],p1[3],p0[3];
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1];

for(i=0;i<nb_points;i++) {mesh->GetPoint(i,p0); SetPoint(i,p0);}
UpdateAll();

int flag=1,count=0,count2=2;
while(flag!=0 && count2!=0)
    {
    count=0;
    pt_index=nb_points;
    for(i=0;i<nb_cells;i++)
        {
        mesh->GetCellCenter(i,c1); 
        GetCellCenter(i,c0);

        for(j=0;j<cells[i][0];j++)
            {
            mesh->GetPoint(pt_index,p); 

            for(l=0;l<cells[i][0];l++)
                {
                i0=cells[i][l+1];   
                i1=cells[i][(l==cells[i][0]-1)?1:(l+2)];    
                GetPoint(i0,p0); GetPoint(i1,p1); 
                GetClosest(closest,p,c0,p0,p1);

                if(closest->nb==3 && closest->side==1 && closest->dist2>0.01)
                    {
                    count++;
                    W=0; for(k=0;k<closest->nb;k++) W+=closest->weights[k]* closest->weights[k]; W=1/W;
                    ColResponse[0]=p[0]-closest->c[0]; ColResponse[1]=p[1]-closest->c[1]; ColResponse[2]=p[2]-closest->c[2];
                    for(k=0;k<closest->nb;k++)
                        {
                        dep[0]=closest->weights[k]* W * ColResponse[0]; dep[1]=closest->weights[k]* W * ColResponse[1]; dep[2]=closest->weights[k]* W * ColResponse[2]; 
                        GetClosest(closest,p,c0,p0,p1);
                        if(closest->pts[k]==0) 
                            {
                            c0[0]+=dep[0]; c0[1]+=dep[1]; c0[2]+=dep[2];
                            SetCellCenter(i,c0);
                            }
                        else if (closest->pts[k]==1) 
                            {
                            p0[0]+=dep[0]; p0[1]+=dep[1]; p0[2]+=dep[2];
                            SetPoint(i0,p0);
                            }
                        else if (closest->pts[k]==2) 
                            {
                            p1[0]+=dep[0]; p1[1]+=dep[1]; p1[2]+=dep[2];
                            SetPoint(i1,p1);
                            }
                        }               
                    }
                if(closest->nb!=0) {    free(closest->weights); free(closest->pts); free(closest);}
                }           
            pt_index++; 
            }

        }
    flag=count;
    count2--;
    }
}

}
*/


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::DeletePointInCell(int c,int p)
{
if(c==-1) return;
int* cell=new int[cells[c][0]]; cell[0]=cells[c][0]-1;
for(int i=0;i<cells[c][0];i++) 
	if(cells[c][i+1]==p) 
		{memcpy(cell+1,cells[c]+1,i*sizeof(int)); memcpy(cell+i+1,cells[c]+i+2,(cells[c][0]-i-1)*sizeof(int));}
SetCell(c,cell); free(cell);
}

void CSimplexSurf::InsertPointInCell(int c,int p,int p1)
{
if(c==-1) return;
// insert p before p1 
int* cell=new int[cells[c][0]+2]; cell[0]=cells[c][0]+1; //bool flag=false;
for(int i=0;i<cells[c][0];i++) 
    if(cells[c][i+1]==p1) 
        {
//		flag=true;
        memcpy(cell+1,cells[c]+1,i*sizeof(int));
        memcpy(cell+i+2,cells[c]+i+1,(cells[c][0]-i)*sizeof(int));
        cell[i+1]=p;
        }
//if(!flag)  return; // debug
SetCell(c,cell); free(cell);
}

//void CSimplexSurf::SetFlagMultires(int index,int val) {flagmultires[index]=val;}
//int CSimplexSurf::GetFlagMultires(int index) {return flagmultires[index];}


DOPCELLSTRUCT* CSimplexSurf::GetDOPCellTree() {return DOPCellTree;}
DOPPOINTSTRUCT* CSimplexSurf::GetDOPPointTree() {return DOPPointTree;}

void CSimplexSurf::BuildDOPCellTree(int depth,P_float tolerance)
{
int i,j,*c;
int *flag_c=new int[nb_cells]; 
P_float** dp=new P_float*[9]; for(i=0;i<9;i++) {dp[i]=new P_float[2*nb_cells]; for(j=0;j<nb_cells;j++) {dp[i][2*j+1]=j; dp[i][2*j]=0;} }
// +-(1,0,0)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[0][2*i]+=points[3*c[j+1]]; dp[0][2*i]=dp[0][2*i]/c[0]; }
// +-(0,1,0)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[1][2*i]+=points[3*c[j+1]+1]; dp[1][2*i]=dp[1][2*i]/c[0]; }
// +-(0,0,1)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[2][2*i]+=points[3*c[j+1]+2]; dp[2][2*i]=dp[2][2*i]/c[0]; }
// +-(1,1,0)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[3][2*i]+=points[3*c[j+1]]+points[3*c[j+1]+1]; dp[3][2*i]=dp[3][2*i]/c[0]; }
// +-(1,-1,0)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[4][2*i]+=points[3*c[j+1]]-points[3*c[j+1]+1]; dp[4][2*i]=dp[4][2*i]/c[0]; }
// +-(1,0,1)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[5][2*i]+=points[3*c[j+1]]+points[3*c[j+1]+2]; dp[5][2*i]=dp[5][2*i]/c[0]; }
// +-(1,0,-1)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[6][2*i]+=points[3*c[j+1]]-points[3*c[j+1]+2]; dp[6][2*i]=dp[6][2*i]/c[0]; }
// +-(0,1,1)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[7][2*i]+=points[3*c[j+1]+1]+points[3*c[j+1]+2]; dp[7][2*i]=dp[7][2*i]/c[0]; }
// +-(0,1,-1)
for(i=0;i<nb_cells;i++) {c=GetCell(i); for(j=0;j<c[0];j++) dp[8][2*i]+=points[3*c[j+1]+1]-points[3*c[j+1]+2]; dp[8][2*i]=dp[8][2*i]/c[0]; }

for(i=0;i<9;i++) qsort(dp[i], nb_cells, 2*sizeof(P_float),qsortP_floatmaxtomin);

DOPCellFree(DOPCellTree);

c=new int[nb_cells]; for(i=0;i<nb_cells;i++) c[i]=i;
for(i=0;i<9;i++) {DOPCellTree->ext[2*i]=dp[i][0]; DOPCellTree->ext[2*i+1]=dp[i][2*(nb_cells-1)];} 
//for(i=0;i<18;i++) DOPCellTree->border[i]=true;
DOPCellTree->nb_cells=nb_cells;
DOPCellTree->nb_points=nb_points;
DOPCellTree->cells=c;
DOPCellTree->dp=dp;

DOPCellSubdivision(DOPCellTree,flag_c,0,depth);
DOPCellBinaryToQuadtree(DOPCellTree);
DOPCellBinaryToQuadtree2(DOPCellTree);
DOPCellExtUpdate(DOPCellTree);
DOPCellInflation(DOPCellTree,tolerance,tolerance*sqrt((P_float)2));

free(flag_c);
}


void CSimplexSurf::PrintDOPCellTree(const char* filename)
{
FILE *f=fopen(filename,"wt");
PrintDOPCellTree(f,DOPCellTree,0);
fclose(f);
}

void CSimplexSurf::PrintDOPCellTree(FILE* f, DOPCELLSTRUCT* dop,int depth)
{
    int i;
    for(i=0;i<depth;i++) fprintf(f,"     ");
    fprintf(f,"%d ",dop->nb_cells);
    fprintf(f,"\n");
    for(i=0;i<depth;i++) fprintf(f,"     ");
    for(i=0;i<dop->nb_cells;i++) fprintf(f,"%d ",dop->cells[i]);
    fprintf(f,"\n");
    for(i=0;i<depth;i++) fprintf(f,"     ");
    for(i=0;i<18;i++) fprintf(f,"%lf ",dop->ext[i]);
    fprintf(f,"\n");
    fprintf(f,"\n");


if(dop->dopp!=NULL) PrintDOPCellTree(f,dop->dopp,depth+1);
if(dop->dopm!=NULL) PrintDOPCellTree(f,dop->dopm,depth+1);
if(dop->dopn!=NULL) PrintDOPCellTree(f,dop->dopn,depth+1);
if(dop->dopo!=NULL) PrintDOPCellTree(f,dop->dopo,depth+1);
}



void CSimplexSurf::BuildDOPPointTree(int depth,P_float tolerance)
{
int i;
int *flag_p=new int[nb_points];
P_float** dp=new P_float*[9]; for(i=0;i<9;i++) dp[i]=new P_float[2*nb_points];
// +-(1,0,0)
for(i=0;i<nb_points;i++) {dp[0][2*i]=points[3*i]; dp[0][2*i+1]=i;}
// +-(0,1,0)
for(i=0;i<nb_points;i++) {dp[1][2*i]=points[3*i+1]; dp[1][2*i+1]=i;}
// +-(0,0,1)
for(i=0;i<nb_points;i++) {dp[2][2*i]=points[3*i+2]; dp[2][2*i+1]=i;}
// +-(1,1,0)
for(i=0;i<nb_points;i++) {dp[3][2*i]=points[3*i]+points[3*i+1]; dp[3][2*i+1]=i;}
// +-(1,-1,0)
for(i=0;i<nb_points;i++) {dp[4][2*i]=points[3*i]-points[3*i+1]; dp[4][2*i+1]=i;}
// +-(1,0,1)
for(i=0;i<nb_points;i++) {dp[5][2*i]=points[3*i]+points[3*i+2]; dp[5][2*i+1]=i;}
// +-(1,0,-1)
for(i=0;i<nb_points;i++) {dp[6][2*i]=points[3*i]-points[3*i+2]; dp[6][2*i+1]=i;}
// +-(0,1,1)
for(i=0;i<nb_points;i++) {dp[7][2*i]=points[3*i+1]+points[3*i+2]; dp[7][2*i+1]=i;}
// +-(0,1,-1)
for(i=0;i<nb_points;i++) {dp[8][2*i]=points[3*i+1]-points[3*i+2]; dp[8][2*i+1]=i;}

for(i=0;i<9;i++) qsort(dp[i], nb_points, 2*sizeof(P_float),qsortP_floatmaxtomin);

DOPPointFree(DOPPointTree);

for(i=0;i<9;i++) {DOPPointTree->ext[2*i]=dp[i][0]; DOPPointTree->ext[2*i+1]=dp[i][2*(nb_points-1)];} //memcpy(DOPPointTree->ext_max,DOPPointTree->ext,18*sizeof(P_float));
//for(i=0;i<18;i++) DOPPointTree->border[i]=true;
DOPPointTree->nb_points=nb_points;
DOPPointTree->points=NULL;
DOPPointTree->dp=dp;

DOPPointSubdivision(DOPPointTree,flag_p,0,depth);
DOPPointBinaryToQuadtree(DOPPointTree);
DOPPointInflation(DOPPointTree,tolerance,tolerance*sqrt((P_float)2));

free(flag_p);
}


void CSimplexSurf::DOPPointBinaryToQuadtree(DOPPOINTSTRUCT* dop)
{
if(dop->dopp==NULL) {dop->dopn=NULL; dop->dopo=NULL; return;}
if(dop->dopp->dopp!=NULL && dop->dopm->dopp!=NULL)
    {
    DOPPOINTSTRUCT* dopm=dop->dopm;     DOPPOINTSTRUCT* dopp=dop->dopp;
    dop->dopn=dop->dopm->dopp;  dop->dopo=dop->dopp->dopm; dop->dopm=dop->dopm->dopm; dop->dopp=dop->dopp->dopp;
    free(dopp->points); free(dopm->points);
    free(dopm); free(dopp);
    DOPPointBinaryToQuadtree(dop->dopm); DOPPointBinaryToQuadtree(dop->dopn); DOPPointBinaryToQuadtree(dop->dopo); DOPPointBinaryToQuadtree(dop->dopp);
    }   
else
    {
    DOP2PointFree(dop->dopp); DOP2PointFree(dop->dopm);
    free(dop->dopp); free(dop->dopm);
    dop->dopm=NULL; dop->dopn=NULL; dop->dopo=NULL; dop->dopp=NULL;
    }
}

void CSimplexSurf::DOPCellBinaryToQuadtree(DOPCELLSTRUCT* dop)
{
if(dop->dopp==NULL) {dop->dopn=NULL; dop->dopo=NULL; return;}
if(dop->dopp->dopp!=NULL && dop->dopm->dopp!=NULL)
    {
    DOPCELLSTRUCT* dopm=dop->dopm;  DOPCELLSTRUCT* dopp=dop->dopp;
    dop->dopn=dop->dopm->dopp;  dop->dopo=dop->dopp->dopm; dop->dopm=dop->dopm->dopm; dop->dopp=dop->dopp->dopp;
    free(dopp->cells); free(dopm->cells);
    free(dopm); free(dopp);
    DOPCellBinaryToQuadtree(dop->dopm); DOPCellBinaryToQuadtree(dop->dopn); DOPCellBinaryToQuadtree(dop->dopo); DOPCellBinaryToQuadtree(dop->dopp);
    }   
else
    {
    DOP2CellFree(dop->dopp); DOP2CellFree(dop->dopm);
    free(dop->dopp); free(dop->dopm);
    dop->dopm=NULL; dop->dopn=NULL; dop->dopo=NULL; dop->dopp=NULL;
    }
}

void CSimplexSurf::DOPCellBinaryToQuadtree2(DOPCELLSTRUCT* dop)
{
if(dop->dopp==NULL)
    {
    if(dop->nb_cells!=1)
        {
        dop->dopm=new DOPCELLSTRUCT;    dop->dopp=new DOPCELLSTRUCT;
        dop->dopm->nb_cells=1; dop->dopp->nb_cells=1;
        dop->dopm->dopm=NULL; dop->dopm->dopn=NULL; dop->dopm->dopo=NULL; dop->dopm->dopp=NULL;     dop->dopp->dopm=NULL; dop->dopp->dopn=NULL; dop->dopp->dopo=NULL; dop->dopp->dopp=NULL; 
        dop->dopm->cells=new int[1]; dop->dopm->cells[0]=dop->cells[0]; dop->dopp->cells=new int[1]; dop->dopp->cells[0]=dop->cells[1];
        if(dop->nb_cells==3)
            {
            dop->dopo=new DOPCELLSTRUCT;    
            dop->dopo->nb_cells=1;
            dop->dopo->dopm=NULL; dop->dopo->dopn=NULL; dop->dopo->dopo=NULL; dop->dopo->dopp=NULL;     
            dop->dopo->cells=new int[1]; dop->dopo->cells[0]=dop->cells[2]; 
            }
        if(dop->nb_cells>3)
            {
            int i=0;
            }
        }
    return;
    }
else {DOPCellBinaryToQuadtree2(dop->dopm); DOPCellBinaryToQuadtree2(dop->dopn); DOPCellBinaryToQuadtree2(dop->dopo); DOPCellBinaryToQuadtree2(dop->dopp);}
}

void CSimplexSurf::UpdateDOPCellTree(P_float tolerance)
{
if(IsAxis) DOPCellExtUpdate_Axis(DOPCellTree); else DOPCellExtUpdate(DOPCellTree);
DOPCellInflation(DOPCellTree,tolerance,tolerance*sqrt((P_float)2));
}

void CSimplexSurf::UpdateDOPPointTree(P_float tolerance)
{
DOPPointExtUpdate(DOPPointTree);
DOPPointInflation(DOPPointTree,tolerance,tolerance*sqrt((P_float)2));
}

void CSimplexSurf::DOPCellExtUpdate(DOPCELLSTRUCT* dop)
{
int i;
if(dop->dopm==NULL && dop->dopn==NULL && dop->dopo==NULL && dop->dopp==NULL)
    {
    P_float val,max,min;
    int j,i2;
    // +-(1,0,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]; if(val>max) max=val; if(val<min) min=val; }  dop->iext[0]=std::max(max,dop->ext[0]); dop->iext[1]=std::min(min,dop->ext[1]); dop->ext[0]=max;    dop->ext[1]=min; 
    // +-(0,1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+1]; if(val>max) max=val; if(val<min) min=val; }  dop->iext[2]=std::max(max,dop->ext[2]); dop->iext[3]=std::min(min,dop->ext[3]);   dop->ext[2]=max;    dop->ext[3]=min; 
    // +-(0,0,1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }  dop->iext[4]=std::max(max,dop->ext[4]); dop->iext[5]=std::min(min,dop->ext[5]);   dop->ext[4]=max;    dop->ext[5]=min; 
    // +-(1,1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]+points[3*i2+1]; if(val>max) max=val; if(val<min) min=val; }     dop->iext[6]=std::max(max,dop->ext[6]); dop->iext[7]=std::min(min,dop->ext[7]);           dop->ext[6]=max;    dop->ext[7]=min; 
    // +-(1,-1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]-points[3*i2+1]; if(val>max) max=val; if(val<min) min=val; }     dop->iext[8]=std::max(max,dop->ext[8]); dop->iext[9]=std::min(min,dop->ext[9]);           dop->ext[8]=max;    dop->ext[9]=min; 
    // +-(1,0,1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]+points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }     dop->iext[10]=std::max(max,dop->ext[10]); dop->iext[11]=std::min(min,dop->ext[11]);       dop->ext[10]=max;   dop->ext[11]=min; 
    // +-(1,0,-1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]-points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }     dop->iext[12]=std::max(max,dop->ext[12]); dop->iext[13]=std::min(min,dop->ext[13]);       dop->ext[12]=max;   dop->ext[13]=min; 
    // +-(0,1,1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+1]+points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }   dop->iext[14]=std::max(max,dop->ext[14]); dop->iext[15]=std::min(min,dop->ext[15]);   dop->ext[14]=max;   dop->ext[15]=min; 
    // +-(0,1,-1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+1]-points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }   dop->iext[16]=std::max(max,dop->ext[16]); dop->iext[17]=std::min(min,dop->ext[17]);   dop->ext[16]=max;   dop->ext[17]=min; 
    }
else
    {
    if(dop->dopm!=NULL) DOPCellExtUpdate(dop->dopm);
    if(dop->dopn!=NULL) DOPCellExtUpdate(dop->dopn);
    if(dop->dopo!=NULL) DOPCellExtUpdate(dop->dopo);
    if(dop->dopp!=NULL) DOPCellExtUpdate(dop->dopp);
    
    for(i=0;i<9;i++) 
        {
        dop->ext[2*i]=-1E10; dop->iext[2*i]=-1E10; 
        if(dop->dopm!=NULL) {if(dop->dopm->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopm->ext[2*i];  if(dop->dopm->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopm->iext[2*i];}
        if(dop->dopn!=NULL) {if(dop->dopn->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopn->ext[2*i];  if(dop->dopn->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopn->iext[2*i];}
        if(dop->dopo!=NULL) {if(dop->dopo->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopo->ext[2*i];  if(dop->dopo->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopo->iext[2*i];}
        if(dop->dopp!=NULL) {if(dop->dopp->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopp->ext[2*i];  if(dop->dopp->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopp->iext[2*i];}

        dop->ext[2*i+1]=1E10; dop->ext[2*i+1]=1E10; 
        if(dop->dopm!=NULL) {if(dop->dopm->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopm->ext[2*i+1];  if(dop->dopm->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopm->iext[2*i+1];}
        if(dop->dopn!=NULL) {if(dop->dopn->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopn->ext[2*i+1];  if(dop->dopn->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopn->iext[2*i+1];}
        if(dop->dopo!=NULL) {if(dop->dopo->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopo->ext[2*i+1];  if(dop->dopo->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopo->iext[2*i+1];}
        if(dop->dopp!=NULL) {if(dop->dopp->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopp->ext[2*i+1];  if(dop->dopp->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopp->iext[2*i+1];}
        }
    }
}



void CSimplexSurf::DOPCellExtUpdate_Axis(DOPCELLSTRUCT* dop)
{
int i;
if(dop->dopm==NULL && dop->dopn==NULL && dop->dopo==NULL && dop->dopp==NULL)
    {
    P_float val,max,min,R;
    int j,i2;
    // +-(1,0,0)
	max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]; R=axiallinks[i2].Rref; if(val+R>max) max=val+R; if(val-R<min) min=val-R; }  dop->iext[0]=std::max(max,dop->ext[0]); dop->iext[1]=std::min(min,dop->ext[1]); dop->ext[0]=max;    dop->ext[1]=min; 
    // +-(0,1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+1]; R=axiallinks[i2].Rref; if(val+R>max) max=val+R; if(val-R<min) min=val-R; }  dop->iext[2]=std::max(max,dop->ext[2]); dop->iext[3]=std::min(min,dop->ext[3]);   dop->ext[2]=max;    dop->ext[3]=min; 
    // +-(0,0,1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+2]; R=axiallinks[i2].Rref; if(val+R>max) max=val+R; if(val-R<min) min=val-R; }  dop->iext[4]=std::max(max,dop->ext[4]); dop->iext[5]=std::min(min,dop->ext[5]);   dop->ext[4]=max;    dop->ext[5]=min; 
    // +-(1,1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]+points[3*i2+1]; R=axiallinks[i2].Rref*sqrt((P_float)2); if(val+R>max) max=val+R; if(val-R<min) min=val-R; }     dop->iext[6]=std::max(max,dop->ext[6]); dop->iext[7]=std::min(min,dop->ext[7]);           dop->ext[6]=max;    dop->ext[7]=min; 
    // +-(1,-1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]-points[3*i2+1]; R=axiallinks[i2].Rref*sqrt((P_float)2); if(val+R>max) max=val+R; if(val-R<min) min=val-R; }     dop->iext[8]=std::max(max,dop->ext[8]); dop->iext[9]=std::min(min,dop->ext[9]);           dop->ext[8]=max;    dop->ext[9]=min; 
    // +-(1,0,1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]+points[3*i2+2]; R=axiallinks[i2].Rref*sqrt((P_float)2); if(val+R>max) max=val+R; if(val-R<min) min=val-R; }     dop->iext[10]=std::max(max,dop->ext[10]); dop->iext[11]=std::min(min,dop->ext[11]);       dop->ext[10]=max;   dop->ext[11]=min; 
    // +-(1,0,-1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2]-points[3*i2+2]; R=axiallinks[i2].Rref*sqrt((P_float)2); if(val+R>max) max=val+R; if(val-R<min) min=val-R; }     dop->iext[12]=std::max(max,dop->ext[12]); dop->iext[13]=std::min(min,dop->ext[13]);       dop->ext[12]=max;   dop->ext[13]=min; 
    // +-(0,1,1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+1]+points[3*i2+2]; R=axiallinks[i2].Rref*sqrt((P_float)2); if(val+R>max) max=val+R; if(val-R<min) min=val-R; }  dop->iext[14]=std::max(max,dop->ext[14]); dop->iext[15]=std::min(min,dop->ext[15]);   dop->ext[14]=max;   dop->ext[15]=min; 
    // +-(0,1,-1)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_cells;i++) for(j=0;j<cells[dop->cells[i]][0];j++) {i2=cells[dop->cells[i]][j+1]; val=points[3*i2+1]-points[3*i2+2]; R=axiallinks[i2].Rref*sqrt((P_float)2); if(val+R>max) max=val+R; if(val-R<min) min=val-R; }  dop->iext[16]=std::max(max,dop->ext[16]); dop->iext[17]=std::min(min,dop->ext[17]);   dop->ext[16]=max;   dop->ext[17]=min; 
    }
else
    {
    if(dop->dopm!=NULL) DOPCellExtUpdate_Axis(dop->dopm);
    if(dop->dopn!=NULL) DOPCellExtUpdate_Axis(dop->dopn);
    if(dop->dopo!=NULL) DOPCellExtUpdate_Axis(dop->dopo);
    if(dop->dopp!=NULL) DOPCellExtUpdate_Axis(dop->dopp);
    
    for(i=0;i<9;i++) 
        {
        dop->ext[2*i]=-1E10; dop->iext[2*i]=-1E10; 
        if(dop->dopm!=NULL) {if(dop->dopm->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopm->ext[2*i];  if(dop->dopm->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopm->iext[2*i];}
        if(dop->dopn!=NULL) {if(dop->dopn->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopn->ext[2*i];  if(dop->dopn->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopn->iext[2*i];}
        if(dop->dopo!=NULL) {if(dop->dopo->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopo->ext[2*i];  if(dop->dopo->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopo->iext[2*i];}
        if(dop->dopp!=NULL) {if(dop->dopp->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopp->ext[2*i];  if(dop->dopp->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopp->iext[2*i];}

        dop->ext[2*i+1]=1E10; dop->ext[2*i+1]=1E10; 
        if(dop->dopm!=NULL) {if(dop->dopm->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopm->ext[2*i+1];  if(dop->dopm->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopm->iext[2*i+1];}
        if(dop->dopn!=NULL) {if(dop->dopn->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopn->ext[2*i+1];  if(dop->dopn->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopn->iext[2*i+1];}
        if(dop->dopo!=NULL) {if(dop->dopo->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopo->ext[2*i+1];  if(dop->dopo->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopo->iext[2*i+1];}
        if(dop->dopp!=NULL) {if(dop->dopp->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopp->ext[2*i+1];  if(dop->dopp->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopp->iext[2*i+1];}
        }
    }
}

void CSimplexSurf::DOPPointExtUpdate(DOPPOINTSTRUCT* dop)
{
int i;
if(dop->dopm==NULL && dop->dopn==NULL && dop->dopo==NULL && dop->dopp==NULL)
    {
    P_float val,max,min;
    int i2;
    // +-(1,0,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=dop->points[i]; val=points[3*i2]; if(val>max) max=val; if(val<min) min=val; } dop->iext[0]=std::max(max,dop->ext[0]); dop->iext[1]=std::min(min,dop->ext[1]);   dop->ext[0]=max;    dop->ext[0]=max;    dop->ext[1]=min; 
    // +-(0,1,0)                                                                                                                        
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2+1]; if(val>max) max=val; if(val<min) min=val; } dop->iext[2]=std::max(max,dop->ext[2]); dop->iext[3]=std::min(min,dop->ext[3]);  dop->ext[2]=max;  dop->ext[2]=max;  dop->ext[3]=min; 
    // +-(0,0,1)                                                                                                                        
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; } dop->iext[4]=std::max(max,dop->ext[4]); dop->iext[5]=std::min(min,dop->ext[5]);  dop->ext[4]=max;  dop->ext[4]=max;  dop->ext[5]=min; 
    // +-(1,1,0)
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2]+points[3*i2+1]; if(val>max) max=val; if(val<min) min=val; }   dop->iext[6]=std::max(max,dop->ext[6]); dop->iext[7]=std::min(min,dop->ext[7]);       dop->ext[6]=max;    dop->ext[7]=min; 
    // +-(1,-1,0)                                                                                                                                   
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2]-points[3*i2+1]; if(val>max) max=val; if(val<min) min=val; }   dop->iext[8]=std::max(max,dop->ext[8]); dop->iext[9]=std::min(min,dop->ext[9]);       dop->ext[8]=max;    dop->ext[9]=min; 
    // +-(1,0,1)                                                                                                                                    
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2]+points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }   dop->iext[10]=std::max(max,dop->ext[10]); dop->iext[11]=std::min(min,dop->ext[11]);   dop->ext[10]=max;   dop->ext[11]=min; 
    // +-(1,0,-1)                                                                                                                                   
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2]-points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; }   dop->iext[12]=std::max(max,dop->ext[12]); dop->iext[13]=std::min(min,dop->ext[13]);   dop->ext[12]=max;   dop->ext[13]=min; 
    // +-(0,1,1)                                                                                                                                    
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2+1]+points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; } dop->iext[14]=std::max(max,dop->ext[14]); dop->iext[15]=std::min(min,dop->ext[15]);   dop->ext[14]=max;   dop->ext[15]=min; 
    // +-(0,1,-1)                                                                                                                                   
    max=-1E10;min=1E10; for(i=0;i<dop->nb_points;i++) {i2=points[i]; val=points[3*i2+1]-points[3*i2+2]; if(val>max) max=val; if(val<min) min=val; } dop->iext[16]=std::max(max,dop->ext[16]); dop->iext[17]=std::min(min,dop->ext[17]);   dop->ext[16]=max;   dop->ext[17]=min; 
    }
else
    {
    if(dop->dopm!=NULL) DOPPointExtUpdate(dop->dopm);
    if(dop->dopn!=NULL) DOPPointExtUpdate(dop->dopn);
    if(dop->dopo!=NULL) DOPPointExtUpdate(dop->dopo);
    if(dop->dopp!=NULL) DOPPointExtUpdate(dop->dopp);
    
    for(i=0;i<9;i++) 
        {
        dop->ext[2*i]=-1E10; dop->iext[2*i]=-1E10; 
        if(dop->dopm!=NULL) {if(dop->dopm->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopm->ext[2*i];  if(dop->dopm->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopm->iext[2*i];}
        if(dop->dopn!=NULL) {if(dop->dopn->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopn->ext[2*i];  if(dop->dopn->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopn->iext[2*i];}
        if(dop->dopo!=NULL) {if(dop->dopo->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopo->ext[2*i];  if(dop->dopo->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopo->iext[2*i];}
        if(dop->dopp!=NULL) {if(dop->dopp->ext[2*i]>dop->ext[2*i]) dop->ext[2*i]=dop->dopp->ext[2*i];  if(dop->dopp->iext[2*i]>dop->iext[2*i]) dop->iext[2*i]=dop->dopp->iext[2*i];}

        dop->ext[2*i+1]=1E10; dop->ext[2*i+1]=1E10; 
        if(dop->dopm!=NULL) {if(dop->dopm->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopm->ext[2*i+1];  if(dop->dopm->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopm->iext[2*i+1];}
        if(dop->dopn!=NULL) {if(dop->dopn->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopn->ext[2*i+1];  if(dop->dopn->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopn->iext[2*i+1];}
        if(dop->dopo!=NULL) {if(dop->dopo->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopo->ext[2*i+1];  if(dop->dopo->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopo->iext[2*i+1];}
        if(dop->dopp!=NULL) {if(dop->dopp->ext[2*i+1]<dop->ext[2*i+1]) dop->ext[2*i+1]=dop->dopp->ext[2*i+1];  if(dop->dopp->iext[2*i+1]<dop->iext[2*i+1]) dop->iext[2*i+1]=dop->dopp->iext[2*i+1];}
        }
    }
}


void CSimplexSurf::DOPCellInflation(DOPCELLSTRUCT* dop,P_float tolerance,P_float tolerance_r2)
{
for(int i=0;i<3;i++) {dop->text[2*i]=dop->ext[2*i]+tolerance; dop->text[2*i+1]=dop->ext[2*i+1]-tolerance;  dop->iext[2*i]+=tolerance; dop->iext[2*i+1]-=tolerance;}
for(int i=3;i<9;i++) {dop->text[2*i]=dop->ext[2*i]+tolerance_r2; dop->text[2*i+1]=dop->ext[2*i+1]-tolerance_r2;  dop->iext[2*i]+=tolerance_r2; dop->iext[2*i+1]-=tolerance_r2;}
if(dop->dopp!=NULL) DOPCellInflation(dop->dopp,tolerance,tolerance_r2);
if(dop->dopm!=NULL) DOPCellInflation(dop->dopm,tolerance,tolerance_r2);
if(dop->dopn!=NULL) DOPCellInflation(dop->dopn,tolerance,tolerance_r2);
if(dop->dopo!=NULL) DOPCellInflation(dop->dopo,tolerance,tolerance_r2);
}


void CSimplexSurf::DOPPointInflation(DOPPOINTSTRUCT* dop,P_float tolerance,P_float tolerance_r2)
{
for(int i=0;i<3;i++) {dop->text[2*i]=dop->ext[2*i]+tolerance; dop->text[2*i+1]=dop->ext[2*i+1]-tolerance;  dop->iext[2*i]+=tolerance; dop->iext[2*i+1]-=tolerance;}
for(int i=3;i<9;i++) {dop->text[2*i]=dop->ext[2*i]+tolerance_r2; dop->text[2*i+1]=dop->ext[2*i+1]-tolerance_r2;  dop->iext[2*i]+=tolerance_r2; dop->iext[2*i+1]-=tolerance_r2;}
if(dop->dopp!=NULL) DOPPointInflation(dop->dopp,tolerance,tolerance_r2);
if(dop->dopm!=NULL) DOPPointInflation(dop->dopm,tolerance,tolerance_r2);
if(dop->dopn!=NULL) DOPPointInflation(dop->dopn,tolerance,tolerance_r2);
if(dop->dopo!=NULL) DOPPointInflation(dop->dopo,tolerance,tolerance_r2);
}


void CSimplexSurf::DOPCellFree(DOPCELLSTRUCT* dop)
{
free(dop->cells);
if(dop->dopp!=NULL) {DOPCellFree(dop->dopp); free(dop->dopp);}
if(dop->dopm!=NULL) {DOPCellFree(dop->dopm); free(dop->dopm);}
if(dop->dopn!=NULL) {DOPCellFree(dop->dopn); free(dop->dopn);}
if(dop->dopo!=NULL) {DOPCellFree(dop->dopo); free(dop->dopo);}
}

void CSimplexSurf::DOP2CellFree(DOPCELLSTRUCT* dop)
{
free(dop->cells);
if(dop->dopp!=NULL) {DOP2CellFree(dop->dopp); free(dop->dopp);}
if(dop->dopm!=NULL) {DOP2CellFree(dop->dopm); free(dop->dopm);}
}


void CSimplexSurf::DOPPointFree(DOPPOINTSTRUCT* dop)
{
free(dop->points);
if(dop->dopp!=NULL) {DOPPointFree(dop->dopp); free(dop->dopp);}
if(dop->dopm!=NULL) {DOPPointFree(dop->dopm); free(dop->dopm);}
if(dop->dopn!=NULL) {DOPPointFree(dop->dopn); free(dop->dopn);}
if(dop->dopo!=NULL) {DOPPointFree(dop->dopo); free(dop->dopo);}
}

void CSimplexSurf::DOP2PointFree(DOPPOINTSTRUCT* dop)
{
free(dop->points);
if(dop->dopp!=NULL) {DOP2PointFree(dop->dopp); free(dop->dopp);}
if(dop->dopm!=NULL) {DOP2PointFree(dop->dopm); free(dop->dopm);}
}

void CSimplexSurf::DOPCellSubdivision(DOPCELLSTRUCT* dop,int* flag_c,int level,int depth)
{
int i,j;

if(level==depth || dop->nb_cells==1)  // end of subdivision
    {
    for(i=0;i<9;i++) free(dop->dp[i]); free(dop->dp);
    dop->dopm=NULL; dop->dopp=NULL;
    }
else    // subdivide
    {
    //ini
    P_float max=0; int ip,im,mode=0;
    dop->dopp=new DOPCELLSTRUCT;    dop->dopm=new DOPCELLSTRUCT;
    dop->dopp->dp=new P_float*[9]; dop->dopm->dp=new P_float*[9];
    for(i=0;i<nb_cells;i++) flag_c[i]=-1;

    // calculate mode
    for(i=0;i<9;i++) if(dop->ext[2*i]-dop->ext[2*i+1]>max) {max=dop->ext[2*i]-dop->ext[2*i+1]; mode=i;}

    // select cells
    for(i=0;i<dop->nb_cells/2;i++) flag_c[(int)dop->dp[mode][2*i+1]]=1;
    int* cp=new int[dop->nb_cells/2];   int* cm=new int[dop->nb_cells-dop->nb_cells/2];
    ip=0;im=0; for(i=0;i<dop->nb_cells;i++) if(flag_c[(int)dop->dp[mode][2*i+1]]==1) {cp[ip]=dop->dp[mode][2*i+1]; ip++;} else {cm[im]=dop->dp[mode][2*i+1]; im++;}
    dop->dopp->nb_cells=ip; dop->dopm->nb_cells=dop->nb_cells-ip;
    dop->dopp->cells=cp; dop->dopm->cells=cm;

    // ini children
    for(i=0;i<9;i++) {dop->dopp->dp[i]=new P_float[2*dop->dopp->nb_cells]; dop->dopm->dp[i]=new P_float[2*dop->dopm->nb_cells];}

    // split dp to children using flag
    for(i=0;i<9;i++) 
        {
        ip=0;im=0;
        for(j=0;j<dop->nb_cells;j++) {if(flag_c[(int)dop->dp[i][2*j+1]]==1) {memcpy(dop->dopp->dp[i]+2*ip,dop->dp[i]+2*j,2*sizeof(P_float)); ip++;} else  {memcpy(dop->dopm->dp[i]+2*im,dop->dp[i]+2*j,2*sizeof(P_float)); im++;} }
        }

    // update children
//  memcpy(dop->dopp->border,dop->border,18*sizeof(bool)); dop->dopp->border[2*mode+1]=false;
//  memcpy(dop->dopm->border,dop->border,18*sizeof(bool)); dop->dopm->border[2*mode]=false;
        // remove interdependencies | n1.n2>0
//      if(mode==0)     {dop->dopm->border[6]=false;dop->dopm->border[8]=false;dop->dopm->border[10]=false;dop->dopm->border[12]=false; dop->dopp->border[7]=false;dop->dopp->border[9]=false;dop->dopp->border[11]=false;dop->dopp->border[13]=false;}
//      if(mode==1)     {dop->dopm->border[6]=false;dop->dopm->border[9]=false;dop->dopm->border[14]=false;dop->dopm->border[16]=false; dop->dopp->border[7]=false;dop->dopp->border[8]=false;dop->dopp->border[15]=false;dop->dopp->border[17]=false;}
//      if(mode==2)     {dop->dopm->border[10]=false;dop->dopm->border[13]=false;dop->dopm->border[14]=false;dop->dopm->border[17]=false;   dop->dopp->border[11]=false;dop->dopp->border[12]=false;dop->dopp->border[15]=false;dop->dopp->border[16]=false;}
//      if(mode==3)     {dop->dopm->border[0]=false;dop->dopm->border[2]=false;dop->dopm->border[10]=false;dop->dopm->border[12]=false;dop->dopm->border[14]=false;dop->dopm->border[16]=false; dop->dopp->border[1]=false;dop->dopp->border[3]=false;dop->dopp->border[11]=false;dop->dopp->border[13]=false;dop->dopp->border[15]=false;dop->dopp->border[17]=false;}
//      if(mode==4)     {dop->dopm->border[0]=false;dop->dopm->border[3]=false;dop->dopm->border[10]=false;dop->dopm->border[12]=false;dop->dopm->border[15]=false;dop->dopm->border[17]=false; dop->dopp->border[1]=false;dop->dopp->border[2]=false;dop->dopp->border[11]=false;dop->dopp->border[13]=false;dop->dopp->border[14]=false;dop->dopp->border[16]=false;}
//      if(mode==5)     {dop->dopm->border[0]=false;dop->dopm->border[4]=false;dop->dopm->border[6]=false;dop->dopm->border[8]=false;dop->dopm->border[14]=false;dop->dopm->border[17]=false;   dop->dopp->border[1]=false;dop->dopp->border[5]=false;dop->dopp->border[7]=false;dop->dopp->border[9]=false;dop->dopp->border[15]=false;dop->dopp->border[16]=false;}
//      if(mode==6)     {dop->dopm->border[0]=false;dop->dopm->border[5]=false;dop->dopm->border[6]=false;dop->dopm->border[8]=false;dop->dopm->border[15]=false;dop->dopm->border[16]=false;   dop->dopp->border[1]=false;dop->dopp->border[4]=false;dop->dopp->border[7]=false;dop->dopp->border[9]=false;dop->dopp->border[14]=false;dop->dopp->border[17]=false;}
//      if(mode==7)     {dop->dopm->border[2]=false;dop->dopm->border[4]=false;dop->dopm->border[6]=false;dop->dopm->border[9]=false;dop->dopm->border[10]=false;dop->dopm->border[13]=false;   dop->dopp->border[3]=false;dop->dopp->border[5]=false;dop->dopp->border[7]=false;dop->dopp->border[8]=false;dop->dopp->border[11]=false;dop->dopp->border[12]=false;}
//      if(mode==8)     {dop->dopm->border[2]=false;dop->dopm->border[5]=false;dop->dopm->border[6]=false;dop->dopm->border[9]=false;dop->dopm->border[11]=false;dop->dopm->border[12]=false;   dop->dopp->border[3]=false;dop->dopp->border[4]=false;dop->dopp->border[7]=false;dop->dopp->border[8]=false;dop->dopp->border[10]=false;dop->dopp->border[13]=false;}

//  for(i=0;i<18;i++) if(dop->dopp->border[i]==false) dop->dopp->ext_max[i]=dop->ext_max[i]; else dop->dopp->ext_max[i]=dop->dopp->ext_min[i];
//  for(i=0;i<18;i++) if(dop->dopm->border[i]==false) dop->dopm->ext_max[i]=dop->ext_max[i]; else dop->dopm->ext_max[i]=dop->dopm->ext_min[i];
    
    for(i=0;i<9;i++) {dop->dopp->ext[2*i]=dop->dopp->dp[i][0]; dop->dopp->ext[2*i+1]=dop->dopp->dp[i][2*(dop->dopp->nb_cells-1)];  dop->dopm->ext[2*i]=dop->dopm->dp[i][0]; dop->dopm->ext[2*i+1]=dop->dopm->dp[i][2*(dop->dopm->nb_cells-1)];}
    for(i=0;i<9;i++) free(dop->dp[i]);  free(dop->dp);

    // recursive subdivision
    DOPCellSubdivision(dop->dopp,flag_c,level+1,depth);
    DOPCellSubdivision(dop->dopm,flag_c,level+1,depth);
    }
}



void CSimplexSurf::DOPPointSubdivision(DOPPOINTSTRUCT* dop,int* flag_p,int level,int depth)
{
int i,j;

if(level==depth || dop->nb_points==1)  // end of subdivision
    {
    for(i=0;i<9;i++) free(dop->dp[i]); free(dop->dp);
    dop->dopm=NULL; dop->dopp=NULL;
    }
else    // subdivide
    {
    //ini
    P_float max=0,tresh; int ip,im,mode;
    dop->dopp=new DOPPOINTSTRUCT;   dop->dopm=new DOPPOINTSTRUCT;
    dop->dopp->dp=new P_float*[9]; dop->dopm->dp=new P_float*[9];
    for(i=0;i<nb_points;i++) flag_p[i]=-1;

    // calculate treshold value
    for(i=0;i<9;i++) if(dop->ext[2*i]-dop->ext[2*i+1]>max) {max=dop->ext[2*i]-dop->ext[2*i+1]; mode=i;}
    tresh=dop->dp[mode][dop->nb_points/2];

    // select points
    ip=0; while(dop->dp[mode][2*ip]>=tresh) {flag_p[(int)dop->dp[mode][2*ip+1]]=1; ip++;}
    dop->dopp->nb_points=ip; dop->dopm->nb_points=dop->nb_points-ip;

    // ini children
    for(i=0;i<9;i++) {dop->dopp->dp[i]=new P_float[2*dop->dopp->nb_points]; dop->dopm->dp[i]=new P_float[2*dop->dopm->nb_points];}
    dop->dopp->points=NULL; dop->dopm->points=NULL;

    // split dp to children using flag
    for(i=0;i<9;i++) 
        {
        ip=0;im=0;
        for(j=0;j<dop->nb_points;j++) {if(flag_p[(int)dop->dp[i][2*j+1]]==1) {memcpy(dop->dopp->dp[i]+2*ip,dop->dp[i]+2*j,2*sizeof(P_float)); ip++;} else  {memcpy(dop->dopm->dp[i]+2*im,dop->dp[i]+2*j,2*sizeof(P_float)); im++;} }
        }
/*
    // update children
    memcpy(dop->dopp->border,dop->border,18*sizeof(bool)); dop->dopp->border[2*mode+1]=false;
    memcpy(dop->dopm->border,dop->border,18*sizeof(bool)); dop->dopm->border[2*mode]=false;
        // remove interdependencies | n1.n2>0
        if(mode==0)     {dop->dopm->border[6]=false;dop->dopm->border[8]=false;dop->dopm->border[10]=false;dop->dopm->border[12]=false; dop->dopp->border[7]=false;dop->dopp->border[9]=false;dop->dopp->border[11]=false;dop->dopp->border[13]=false;}
        if(mode==1)     {dop->dopm->border[6]=false;dop->dopm->border[9]=false;dop->dopm->border[14]=false;dop->dopm->border[16]=false; dop->dopp->border[7]=false;dop->dopp->border[8]=false;dop->dopp->border[15]=false;dop->dopp->border[17]=false;}
        if(mode==2)     {dop->dopm->border[10]=false;dop->dopm->border[13]=false;dop->dopm->border[14]=false;dop->dopm->border[17]=false;   dop->dopp->border[11]=false;dop->dopp->border[12]=false;dop->dopp->border[15]=false;dop->dopp->border[16]=false;}
        if(mode==3)     {dop->dopm->border[0]=false;dop->dopm->border[2]=false;dop->dopm->border[10]=false;dop->dopm->border[12]=false;dop->dopm->border[14]=false;dop->dopm->border[16]=false; dop->dopp->border[1]=false;dop->dopp->border[3]=false;dop->dopp->border[11]=false;dop->dopp->border[13]=false;dop->dopp->border[15]=false;dop->dopp->border[17]=false;}
        if(mode==4)     {dop->dopm->border[0]=false;dop->dopm->border[3]=false;dop->dopm->border[10]=false;dop->dopm->border[12]=false;dop->dopm->border[15]=false;dop->dopm->border[17]=false; dop->dopp->border[1]=false;dop->dopp->border[2]=false;dop->dopp->border[11]=false;dop->dopp->border[13]=false;dop->dopp->border[14]=false;dop->dopp->border[16]=false;}
        if(mode==5)     {dop->dopm->border[0]=false;dop->dopm->border[4]=false;dop->dopm->border[6]=false;dop->dopm->border[8]=false;dop->dopm->border[14]=false;dop->dopm->border[17]=false;   dop->dopp->border[1]=false;dop->dopp->border[5]=false;dop->dopp->border[7]=false;dop->dopp->border[9]=false;dop->dopp->border[15]=false;dop->dopp->border[16]=false;}
        if(mode==6)     {dop->dopm->border[0]=false;dop->dopm->border[5]=false;dop->dopm->border[6]=false;dop->dopm->border[8]=false;dop->dopm->border[15]=false;dop->dopm->border[16]=false;   dop->dopp->border[1]=false;dop->dopp->border[4]=false;dop->dopp->border[7]=false;dop->dopp->border[9]=false;dop->dopp->border[14]=false;dop->dopp->border[17]=false;}
        if(mode==7)     {dop->dopm->border[2]=false;dop->dopm->border[4]=false;dop->dopm->border[6]=false;dop->dopm->border[9]=false;dop->dopm->border[10]=false;dop->dopm->border[13]=false;   dop->dopp->border[3]=false;dop->dopp->border[5]=false;dop->dopp->border[7]=false;dop->dopp->border[8]=false;dop->dopp->border[11]=false;dop->dopp->border[12]=false;}
        if(mode==8)     {dop->dopm->border[2]=false;dop->dopm->border[5]=false;dop->dopm->border[6]=false;dop->dopm->border[9]=false;dop->dopm->border[11]=false;dop->dopm->border[12]=false;   dop->dopp->border[3]=false;dop->dopp->border[4]=false;dop->dopp->border[7]=false;dop->dopp->border[8]=false;dop->dopp->border[10]=false;dop->dopp->border[13]=false;}
*/
    for(i=0;i<9;i++) {dop->dopp->ext[2*i]=dop->dopp->dp[i][0]; dop->dopp->ext[2*i+1]=dop->dopp->dp[i][2*(dop->dopp->nb_points-1)];  dop->dopm->ext[2*i]=dop->dopm->dp[i][0]; dop->dopm->ext[2*i+1]=dop->dopm->dp[i][2*(dop->dopm->nb_points-1)];}

    for(i=0;i<9;i++) free(dop->dp[i]);  free(dop->dp);

    dop->dopp->points=new int[dop->dopp->nb_points];
    for(j=0;j<dop->dopp->nb_points;j++) dop->dopp->points[j]=(int)dop->dopp->dp[0][2*j+1]; 
    dop->dopm->points=new int[dop->dopm->nb_points];
    for(j=0;j<dop->dopm->nb_points;j++) dop->dopm->points[j]=(int)dop->dopm->dp[0][2*j+1]; 

    // recursive subdivision
    DOPPointSubdivision(dop->dopp,flag_p,level+1,depth);
    DOPPointSubdivision(dop->dopm,flag_p,level+1,depth);
    }
}

void CSimplexSurf::DOPInclusion(P_float p[3],CLOSESTSTRUCT *crossing2, DOPCELLSTRUCT* dop2)
{
if(dop2->dopm==NULL && dop2->dopn==NULL && dop2->dopo==NULL && dop2->dopp==NULL) 
    {
    for(int j=0;j<dop2->nb_cells;j++) GetCrossing_C_ini(crossing2,p,dop2->cells[j]);
    }
else 
	{
    if(dop2->dopp!=NULL) if(p[0]>=dop2->dopp->ext[1] && p[1]<=dop2->dopp->ext[2] && p[1]>=dop2->dopp->ext[3] &&   p[2]<=dop2->dopp->ext[4] && p[2]>=dop2->dopp->ext[5] &&   (p[0]+p[1])>=dop2->dopp->ext[7] && (p[0]-p[1])>=dop2->dopp->ext[9] &&  (p[0]+p[2])>=dop2->dopp->ext[11] &&   (p[0]-p[2])>=dop2->dopp->ext[13] &&  (p[1]+p[2])<=dop2->dopp->ext[14] && (p[1]+p[2])>=dop2->dopp->ext[15] &&  (p[1]-p[2])<=dop2->dopp->ext[16] && (p[1]-p[2])>=dop2->dopp->ext[17] ) DOPInclusion(p,crossing2,dop2->dopp); 
    if(dop2->dopm!=NULL) if(p[0]>=dop2->dopm->ext[1] && p[1]<=dop2->dopm->ext[2] && p[1]>=dop2->dopm->ext[3] &&   p[2]<=dop2->dopm->ext[4] && p[2]>=dop2->dopm->ext[5] &&   (p[0]+p[1])>=dop2->dopm->ext[7] && (p[0]-p[1])>=dop2->dopm->ext[9] &&  (p[0]+p[2])>=dop2->dopm->ext[11] &&   (p[0]-p[2])>=dop2->dopm->ext[13] &&  (p[1]+p[2])<=dop2->dopm->ext[14] && (p[1]+p[2])>=dop2->dopm->ext[15] &&  (p[1]-p[2])<=dop2->dopm->ext[16] && (p[1]-p[2])>=dop2->dopm->ext[17] ) DOPInclusion(p,crossing2,dop2->dopm); 
    if(dop2->dopn!=NULL) if(p[0]>=dop2->dopn->ext[1] && p[1]<=dop2->dopn->ext[2] && p[1]>=dop2->dopn->ext[3] &&   p[2]<=dop2->dopn->ext[4] && p[2]>=dop2->dopn->ext[5] &&   (p[0]+p[1])>=dop2->dopn->ext[7] && (p[0]-p[1])>=dop2->dopn->ext[9] &&  (p[0]+p[2])>=dop2->dopn->ext[11] &&   (p[0]-p[2])>=dop2->dopn->ext[13] &&  (p[1]+p[2])<=dop2->dopn->ext[14] && (p[1]+p[2])>=dop2->dopn->ext[15] &&  (p[1]-p[2])<=dop2->dopn->ext[16] && (p[1]-p[2])>=dop2->dopn->ext[17] ) DOPInclusion(p,crossing2,dop2->dopn); 
    if(dop2->dopo!=NULL) if(p[0]>=dop2->dopo->ext[1] && p[1]<=dop2->dopo->ext[2] && p[1]>=dop2->dopo->ext[3] &&   p[2]<=dop2->dopo->ext[4] && p[2]>=dop2->dopo->ext[5] &&   (p[0]+p[1])>=dop2->dopo->ext[7] && (p[0]-p[1])>=dop2->dopo->ext[9] &&  (p[0]+p[2])>=dop2->dopo->ext[11] &&   (p[0]-p[2])>=dop2->dopo->ext[13] &&  (p[1]+p[2])<=dop2->dopo->ext[14] && (p[1]+p[2])>=dop2->dopo->ext[15] &&  (p[1]-p[2])<=dop2->dopo->ext[16] && (p[1]-p[2])>=dop2->dopo->ext[17] ) DOPInclusion(p,crossing2,dop2->dopo); 
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateMass()
{
int i;
if(nb_points==1) return;
for(i=0;i<3*nb_points;i++) {mass_inv[i]=mass_scalar*mass_inv[i]; mass[i]=mass[i]/mass_scalar;}

// compute mean radius
/*P_float mean=0;
int n[3]; P_float p[3],n0[3],n1[3],n2[3],e[3];
for(i=0;i<nb_points;i++) 
    {
    GetPoint(i,p); GetNeighbors(i,n); GetPoint(n[0],n0);GetPoint(n[1],n1);GetPoint(n[2],n2);
    e[0]=n0[0]-p[0];e[1]=n0[1]-p[1];e[2]=n0[2]-p[2]; mean+=norm(e);
    e[0]=n1[0]-p[0];e[1]=n1[1]-p[1];e[2]=n1[2]-p[2]; mean+=norm(e);
    e[0]=n2[0]-p[0];e[1]=n2[1]-p[1];e[2]=n2[2]-p[2]; mean+=norm(e);
    }
mean=mean/(P_float)(3*nb_points);
mass_scalar=0.015*mean*mean;
*/

mass_scalar=0.1*surface/(P_float)nb_points;

for(i=0;i<3*nb_points;i++) {mass_inv[i]=mass_inv[i]/mass_scalar; mass[i]=mass[i]*mass_scalar;}
mass_modified=false;
}

void CSimplexSurf::InsertNextPoint(P_float p[3])
{
nb_points++;
points=(P_float*)realloc(points,3*nb_points*sizeof(P_float)); memcpy(points+3*(nb_points-1),p,3*sizeof(P_float));
points_tm1=(P_float*)realloc(points_tm1,3*nb_points*sizeof(P_float)); memcpy(points_tm1+3*(nb_points-1),p,3*sizeof(P_float));

neighbors=(int*)realloc(neighbors,3*nb_points*sizeof(int)); 
neighbors[3*(nb_points-1)]=-1; 
neighbors[3*(nb_points-1)+1]=-1; 
neighbors[3*(nb_points-1)+2]=-1;

neighbors_c=(int*)realloc(neighbors_c,3*nb_points*sizeof(int)); 
params=(P_float*)realloc(params,3*nb_points*sizeof(P_float)); 
//params_multires=(P_float*)realloc(params_multires,3*nb_points*sizeof(P_float)); 
speeds=(P_float*)realloc(speeds,3*nb_points*sizeof(P_float)); SetSpeed(nb_points-1,0,0,0);
speeds_tm1=(P_float*)realloc(speeds_tm1,3*nb_points*sizeof(P_float)); SetSpeed_tm1(nb_points-1,0,0,0);
Val_p=(P_float*)realloc(Val_p,3*nb_points*sizeof(P_float)); P_float O[3]={0,0,0}; SetVal_p(nb_points-1,O);
df_p=(P_float*)realloc(df_p,6*nb_points*sizeof(P_float)); 
df_s=(P_float*)realloc(df_s,6*nb_points*sizeof(P_float)); 
forces=(P_float*)realloc(forces,3*nb_points*sizeof(P_float)); 
mass_inv=(P_float*)realloc(mass_inv,3*nb_points*sizeof(P_float)); mass=(P_float*)realloc(mass,3*nb_points*sizeof(P_float)); SetMass(nb_points-1);
mass_constraint=(P_float**)realloc(mass_constraint,nb_points*sizeof(P_float*)); mass_constraint[nb_points-1]=NULL;
normals=(P_float*)realloc(normals,3*nb_points*sizeof(P_float)); 
if(normals2!=NULL) normals2=(P_float*)realloc(normals2,3*nb_points*sizeof(P_float)); 
surfaces=(P_float*)realloc(surfaces,nb_points*sizeof(P_float)); 
//volumes=(P_float*)realloc(volumes,nb_points*sizeof(P_float)); 
h=(P_float*)realloc(h,nb_points*sizeof(P_float)); 
if(axiallinks!=NULL) {axiallinks=(AXIALLINK*)realloc(axiallinks,nb_points*sizeof(AXIALLINK)); axiallinks[nb_points-1].pts=NULL; axiallinks[nb_points-1].w=NULL; axiallinks[nb_points-1].f[0]=0; axiallinks[nb_points-1].f[1]=0; axiallinks[nb_points-1].f[2]=0;} 
externalforces=(P_float*)realloc(externalforces,3*nb_points*sizeof(P_float)); 
externalforces_ign=(bool*)realloc(externalforces_ign,nb_points*sizeof(bool)); 
constraintpointsforces=(P_float*)realloc(constraintpointsforces,4*nb_points*sizeof(P_float)); 
attached_point=(int*)realloc(attached_point,nb_points*sizeof(int)); attached_point[nb_points-1]=0; 
neighbors2=(int***)realloc(neighbors2,nb_points*sizeof(int**)); neighbors2[nb_points-1]=(int**)malloc(NB_NEIGHBORHOOD*sizeof(int*)); for(int i=0;i<NB_NEIGHBORHOOD;i++) {neighbors2[nb_points-1][i]=new int[1]; neighbors2[nb_points-1][i][0]=0;}
}

void CSimplexSurf::DeletePoint(int index)
{
int i,j;
for(i=0;i<nb_cells;i++) for(j=0;j<cells[i][0];j++) if(cells[i][j+1]>index) cells[i][j+1]=cells[i][j+1]-1;
for(i=0;i<nb_points;i++) for(j=0;j<3;j++) if(neighbors[3*i+j]>index) neighbors[3*i+j]=neighbors[3*i+j]-1;

nb_points--;

P_float* buff=(P_float*)malloc(3*(nb_points-index)*sizeof(P_float));

memcpy(buff,points+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(points+3*index,buff,3*(nb_points-index)*sizeof(P_float)); points=(P_float*)realloc(points,3*nb_points*sizeof(P_float));
memcpy(buff,points_tm1+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(points_tm1+3*index,buff,3*(nb_points-index)*sizeof(P_float)); points_tm1=(P_float*)realloc(points_tm1,3*nb_points*sizeof(P_float));
memcpy(buff,speeds+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(speeds+3*index,buff,3*(nb_points-index)*sizeof(P_float)); speeds=(P_float*)realloc(speeds,3*nb_points*sizeof(P_float));
memcpy(buff,speeds_tm1+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(speeds_tm1+3*index,buff,3*(nb_points-index)*sizeof(P_float)); speeds_tm1=(P_float*)realloc(speeds_tm1,3*nb_points*sizeof(P_float));
memcpy(buff,externalforces+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(externalforces+3*index,buff,3*(nb_points-index)*sizeof(P_float)); externalforces=(P_float*)realloc(externalforces,3*nb_points*sizeof(P_float)); 
memcpy(buff,mass_inv+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(mass_inv+3*index,buff,3*(nb_points-index)*sizeof(P_float)); mass_inv=(P_float*)realloc(mass_inv,3*nb_points*sizeof(P_float));
memcpy(buff,mass+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(mass+3*index,buff,3*(nb_points-index)*sizeof(P_float)); mass=(P_float*)realloc(mass,3*nb_points*sizeof(P_float));
memcpy(buff,Val_p+3*(index+1),3*(nb_points-index)*sizeof(P_float)); memcpy(Val_p+3*index,buff,3*(nb_points-index)*sizeof(P_float)); Val_p=(P_float*)realloc(Val_p,3*nb_points*sizeof(P_float));
free(buff);

P_float** buffp=(P_float**)malloc((nb_points-index)*sizeof(P_float*));
memcpy(buffp,mass_constraint+(index+1),(nb_points-index)*sizeof(P_float*)); memcpy(mass_constraint+index,buffp,(nb_points-index)*sizeof(P_float*)); mass_constraint=(P_float**)realloc(mass_constraint,nb_points*sizeof(P_float*));
free(buffp);

int* buffi=(int*)malloc(3*(nb_points-index)*sizeof(int));
memcpy(buffi,neighbors+3*(index+1),3*(nb_points-index)*sizeof(int)); memcpy(neighbors+3*index,buffi,3*(nb_points-index)*sizeof(int)); neighbors=(int*)realloc(neighbors,3*nb_points*sizeof(int)); 
memcpy(buffi,neighbors_c+3*(index+1),3*(nb_points-index)*sizeof(int)); memcpy(neighbors_c+3*index,buffi,3*(nb_points-index)*sizeof(int)); neighbors_c=(int*)realloc(neighbors_c,3*nb_points*sizeof(int)); 
memcpy(buffi,attached_point+(index+1),(nb_points-index)*sizeof(int)); memcpy(attached_point+index,buffi,(nb_points-index)*sizeof(int)); attached_point=(int*)realloc(attached_point,nb_points*sizeof(int));  
free(buffi);

if(axiallinks!=NULL) axiallinks=(AXIALLINK*)realloc(axiallinks,nb_points*sizeof(AXIALLINK)); 
surfaces=(P_float*)realloc(surfaces,nb_points*sizeof(P_float)); 
//volumes=(P_float*)realloc(volumes,nb_points*sizeof(P_float)); 
h=(P_float*)realloc(h,nb_points*sizeof(P_float)); 
params=(P_float*)realloc(params,3*nb_points*sizeof(P_float)); 
//params_multires=(P_float*)realloc(params_multires,3*nb_points*sizeof(P_float)); 
forces=(P_float*)realloc(forces,3*nb_points*sizeof(P_float)); 
df_p=(P_float*)realloc(df_p,6*nb_points*sizeof(P_float)); 
df_s=(P_float*)realloc(df_s,6*nb_points*sizeof(P_float)); 
normals=(P_float*)realloc(normals,3*nb_points*sizeof(P_float)); 
if(normals2!=NULL) normals2=(P_float*)realloc(normals2,3*nb_points*sizeof(P_float)); 
neighbors2=(int***)realloc(neighbors2,nb_points*sizeof(int**)); 
constraintpointsforces=(P_float*)realloc(constraintpointsforces,4*nb_points*sizeof(P_float)); 
externalforces_ign=(bool*)realloc(externalforces_ign,nb_points*sizeof(bool)); 
}


/*
int* CSimplexSurf::GetDOPPointClosestPoints(CSimplexSurf* mesh)
{
nb_DOPIntersect=0; DOPPointsIntersect=new int*[1]; nb_DOPPointsIntersect=new int[1];
DOPPointIntersect(DOPPointTree,mesh->GetDOPPointTree());
if(nb_DOPIntersect==0) return NULL;

int i,j,k,l,m,n,pt_index;
P_float *p,*p2,d;
P_float **dist1,**dist2;

P_float *dmin1=new P_float[nb_points]; P_float *dmin2=new P_float[mesh->GetNumberOfPoints()];
int *i1=new int[nb_points]; int *i2=new int[mesh->GetNumberOfPoints()];
for(i=0;i<nb_points;i++) {dmin1[i]=1E10; i1[i]=-1;}
for(i=0;i<mesh->GetNumberOfPoints();i++) {dmin2[i]=1E10; i2[i]=-1;}

for(i=0;i<nb_DOPIntersect;i++)
    {
    // precompute distances
    dist1=new P_float*[nb_DOPPointsIntersect[2*i]]; for(j=0;j<nb_DOPPointsIntersect[2*i];j++) dist1[j]=new P_float[nb_DOPPointsIntersect[2*i+1]];
    dist2=new P_float*[nb_DOPPointsIntersect[2*i+1]]; for(j=0;j<nb_DOPPointsIntersect[2*i+1];j++) dist2[j]=new P_float[nb_DOPPointsIntersect[2*i]];
    for(j=0;j<nb_DOPPointsIntersect[2*i];j++)
        {
        p=GetPoint(DOPPointsIntersect[2*i][j]);
        for(k=0;k<nb_DOPPointsIntersect[2*i+1];k++) { p2=mesh->GetPoint(DOPPointsIntersect[2*i+1][k]); d=(p[0]-p2[0])*(p[0]-p2[0]) + (p[1]-p2[1])*(p[1]-p2[1]) + (p[2]-p2[2])*(p[2]-p2[2]); dist1[j][k]=d; dist2[k][j]=d;   }
        }
    // search for minimum distances
    for(j=0;j<nb_DOPPointsIntersect[2*i];j++) for(k=0;k<nb_DOPPointsIntersect[2*i+1];k++) if(dist1[j][k]<dmin1[DOPPointsIntersect[2*i][j]]) {i1[DOPPointsIntersect[2*i][j]]=DOPPointsIntersect[2*i+1][k]; dmin1[DOPPointsIntersect[2*i][j]]=dist1[j][k];} 
    for(j=0;j<nb_DOPPointsIntersect[2*i+1];j++) for(k=0;k<nb_DOPPointsIntersect[2*i];k++) if(dist2[j][k]<dmin2[DOPPointsIntersect[2*i+1][j]]) {i2[DOPPointsIntersect[2*i+1][j]]=DOPPointsIntersect[2*i][k]; dmin2[DOPPointsIntersect[2*i+1][j]]=dist2[j][k];} 
    // free
    for(j=0;j<nb_DOPPointsIntersect[2*i];j++) free(dist1[j]); free(dist1);
    for(j=0;j<nb_DOPPointsIntersect[2*i+1];j++) free(dist2[j]); free(dist2);
    }

// return result
n=0;
for(i=0;i<nb_points;i++) if(i1[i]!=-1) n++;
for(i=0;i<mesh->GetNumberOfPoints();i++) if(i2[i]!=-1) n++;

int* closestpoints= new int[2*n+1];
closestpoints[0]=n;

n=1;
for(i=0;i<nb_points;i++) if(i1[i]!=-1) {closestpoints[n]=-i; n++; closestpoints[n]=i1[i]; n++;}
for(i=0;i<mesh->GetNumberOfPoints();i++) if(i2[i]!=-1) {closestpoints[n]=i2[i]; n++; closestpoints[n]=-i; n++;}

free(i1); free(i2); free(dmin1); free(dmin2);
free(nb_DOPPointsIntersect);free(DOPPointsIntersect);
return closestpoints;
}

int* CSimplexSurf::GetDOPCellClosestPoints(CSimplexSurf* mesh)
{
nb_DOPIntersect=0; DOPPointsIntersect=new int*[1]; nb_DOPPointsIntersect=new int[1];
DOPCellIntersect(DOPCellTree,mesh->GetDOPCellTree(),this,mesh);
if(nb_DOPIntersect==0) return NULL;

int i,j,k,l,m,n,pt_index;
P_float *p,*p2,d;
P_float **dist1,**dist2;

P_float *dmin1=new P_float[nb_points]; P_float *dmin2=new P_float[mesh->GetNumberOfPoints()];
int *i1=new int[nb_points]; int *i2=new int[mesh->GetNumberOfPoints()];
for(i=0;i<nb_points;i++) {dmin1[i]=1E10; i1[i]=-1;}
for(i=0;i<mesh->GetNumberOfPoints();i++) {dmin2[i]=1E10; i2[i]=-1;}

for(i=0;i<nb_DOPIntersect;i++)
    {
    // precompute distances
    dist1=new P_float*[nb_DOPPointsIntersect[2*i]]; for(j=0;j<nb_DOPPointsIntersect[2*i];j++) dist1[j]=new P_float[nb_DOPPointsIntersect[2*i+1]];
    dist2=new P_float*[nb_DOPPointsIntersect[2*i+1]]; for(j=0;j<nb_DOPPointsIntersect[2*i+1];j++) dist2[j]=new P_float[nb_DOPPointsIntersect[2*i]];
    for(j=0;j<nb_DOPPointsIntersect[2*i];j++)
        {
        p=GetPoint(DOPPointsIntersect[2*i][j]);
        for(k=0;k<nb_DOPPointsIntersect[2*i+1];k++) { p2=mesh->GetPoint(DOPPointsIntersect[2*i+1][k]); d=(p[0]-p2[0])*(p[0]-p2[0]) + (p[1]-p2[1])*(p[1]-p2[1]) + (p[2]-p2[2])*(p[2]-p2[2]); dist1[j][k]=d; dist2[k][j]=d;   }
        }
    // search for minimum distances
    for(j=0;j<nb_DOPPointsIntersect[2*i];j++) for(k=0;k<nb_DOPPointsIntersect[2*i+1];k++) if(dist1[j][k]<dmin1[DOPPointsIntersect[2*i][j]]) {i1[DOPPointsIntersect[2*i][j]]=DOPPointsIntersect[2*i+1][k]; dmin1[DOPPointsIntersect[2*i][j]]=dist1[j][k];} 
    for(j=0;j<nb_DOPPointsIntersect[2*i+1];j++) for(k=0;k<nb_DOPPointsIntersect[2*i];k++) if(dist2[j][k]<dmin2[DOPPointsIntersect[2*i+1][j]]) {i2[DOPPointsIntersect[2*i+1][j]]=DOPPointsIntersect[2*i][k]; dmin2[DOPPointsIntersect[2*i+1][j]]=dist2[j][k];} 
    // free
    for(j=0;j<nb_DOPPointsIntersect[2*i];j++) free(dist1[j]); free(dist1);
    for(j=0;j<nb_DOPPointsIntersect[2*i+1];j++) free(dist2[j]); free(dist2);
    }

// return result
n=0;
for(i=0;i<nb_points;i++) if(i1[i]!=-1) n++;
for(i=0;i<mesh->GetNumberOfPoints();i++) if(i2[i]!=-1) n++;

int* closestpoints= new int[2*n+1];
closestpoints[0]=n;

n=1;
for(i=0;i<nb_points;i++) if(i1[i]!=-1) {closestpoints[n]=-i; n++; closestpoints[n]=i1[i]; n++;}
for(i=0;i<mesh->GetNumberOfPoints();i++) if(i2[i]!=-1) {closestpoints[n]=i2[i]; n++; closestpoints[n]=-i; n++;}

free(i1); free(i2); free(dmin1); free(dmin2);
free(nb_DOPPointsIntersect);
for(i=0;i<nb_DOPIntersect;i++) free(DOPPointsIntersect[i]);
free(DOPPointsIntersect);
return closestpoints;
}

*/

int CSimplexSurf::GetClosestPoint(P_float p[3])
{
P_float d,dmin=1E10; int index;
for(int i=0;i<nb_points;i++)
    {
    d=(points[3*i]-p[0])*(points[3*i]-p[0])+(points[3*i+1]-p[1])*(points[3*i+1]-p[1])+(points[3*i+2]-p[2])*(points[3*i+2]-p[2]);
    if(d<dmin) {index=i; dmin=d;}
    }
return index;
}

int CSimplexSurf::GetClosestBorderPoint(P_float p[3])
{
P_float d,dmin=1E10; int index=-1;
for(int i=0;i<nb_points;i++)
	if(GetNeighbors_c(i)[2]==-1)
		{
		d=(points[3*i]-p[0])*(points[3*i]-p[0])+(points[3*i+1]-p[1])*(points[3*i+1]-p[1])+(points[3*i+2]-p[2])*(points[3*i+2]-p[2]);
		if(d<dmin) {index=i; dmin=d;}
		}
return index;
}


int CSimplexSurf::GetClosestAttachPoint(P_float p[3])
{
P_float d,dmin=1E10; int index=-1;
for(int i=0;i<nb_points;i++)
	if(attached_point[i]!=0)
		{
		d=(points[3*i]-p[0])*(points[3*i]-p[0])+(points[3*i+1]-p[1])*(points[3*i+1]-p[1])+(points[3*i+2]-p[2])*(points[3*i+2]-p[2]);
		if(d<dmin) {index=i; dmin=d;}
		}
return index;
}

int CSimplexSurf::GetClosestAttachPoint(P_float p[3],int attachindex)
{
P_float d,dmin=1E10; int index=-1;
for(int i=0;i<nb_points;i++)
	if((attachindex<0 && (attached_point[i]==attachindex || attached_point[i]==-attachindex)) || (attachindex>0 && attached_point[i]==attachindex))
		{
		d=(points[3*i]-p[0])*(points[3*i]-p[0])+(points[3*i+1]-p[1])*(points[3*i+1]-p[1])+(points[3*i+2]-p[2])*(points[3*i+2]-p[2]);
		if(d<dmin) {index=i; dmin=d;}
		}
return index;
}


void CSimplexSurf::GetClosestBorder2Points(int closests[2],P_float dists2[2],P_float p[3])
{
P_float d;
dists2[0]=1E10; dists2[1]=1E10; 
for(int i=0;i<nb_points;i++)
	if(GetNeighbors_c(i)[2]==-1)
    {
    d=(points[3*i]-p[0])*(points[3*i]-p[0])+(points[3*i+1]-p[1])*(points[3*i+1]-p[1])+(points[3*i+2]-p[2])*(points[3*i+2]-p[2]);
    if(d<dists2[0]) 
		{
		closests[1]=closests[0]; closests[0]=i; 
		dists2[1]=dists2[0]; dists2[0]=d; 
		}
    else if(d<dists2[1]) 
		{
		closests[1]=i; 
		dists2[1]=d; 
		}
    }
}



void CSimplexSurf::GetClosest3Points(int closests[3],P_float dists2[3],P_float p[3])
{
P_float d;
dists2[0]=1E10; dists2[1]=1E10; dists2[2]=1E10;
for(int i=0;i<nb_points;i++)
    {
    d=(points[3*i]-p[0])*(points[3*i]-p[0])+(points[3*i+1]-p[1])*(points[3*i+1]-p[1])+(points[3*i+2]-p[2])*(points[3*i+2]-p[2]);
    if(d<dists2[0]) 
		{
		closests[2]=closests[1]; closests[1]=closests[0]; closests[0]=i; 
		dists2[2]=dists2[1]; dists2[1]=dists2[0]; dists2[0]=d; 
		}
    else if(d<dists2[1]) 
		{
		closests[2]=closests[1]; closests[1]=i; 
		dists2[2]=dists2[1]; dists2[1]=d; 
		}
    else if(d<dists2[2]) 
		{
		closests[2]=i; 
		dists2[2]=d; 
		}
    }
}

void CSimplexSurf::GetClosest(CLOSESTSTRUCT* closest,P_float p[3]) 
{
int CP=GetClosestPoint(p);
if(CP==-1)   return;
GetClosest_P(closest,p,CP);
}


void CSimplexSurf::GetCrossing_P(CLOSESTSTRUCT* closest,P_float p[3],P_float s[3],int pointindex)
{
int i,*c,n[3],nc[3]; 
GetNeighbors(pointindex,n);
GetNeighbors_c(pointindex,nc); 

P_float sc[3][3]={0,0,0,0,0,0,0,0,0};
int ntoc[3][2]={-1,-1,-1,-1,-1,-1};
c=GetCell(nc[0]); for(i=0;i<c[0];i++) {sc[0][0]+=speeds[3*c[i+1]];sc[0][1]+=speeds[3*c[i+1]+1];sc[0][2]+=speeds[3*c[i+1]+2];if(c[i+1]==n[0]) ntoc[0][0]=0; if(c[i+1]==n[1]) ntoc[1][0]=0;  if(c[i+1]==n[2]) ntoc[2][0]=0;  } sc[0][0]=sc[0][0]/c[0]; sc[0][1]=sc[0][1]/c[0]; sc[0][2]=sc[0][2]/c[0];
c=GetCell(nc[1]); for(i=0;i<c[0];i++) {sc[1][0]+=speeds[3*c[i+1]];sc[1][1]+=speeds[3*c[i+1]+1];sc[1][2]+=speeds[3*c[i+1]+2]; if(c[i+1]==n[0]) ntoc[0][(ntoc[0][0]==-1)?0:1]=1; if(c[i+1]==n[1]) ntoc[1][(ntoc[1][0]==-1)?0:1]=1;   if(c[i+1]==n[2]) ntoc[2][(ntoc[2][0]==-1)?0:1]=1;} sc[1][0]=sc[1][0]/c[0]; sc[1][1]=sc[1][1]/c[0]; sc[1][2]=sc[1][2]/c[0];
if(nc[2]!=-1) {c=GetCell(nc[2]); for(i=0;i<c[0];i++) {sc[2][0]+=speeds[3*c[i+1]];sc[2][1]+=speeds[3*c[i+1]+1];sc[2][2]+=speeds[3*c[i+1]+2]; if(c[i+1]==n[0]) ntoc[0][1]=2; if(c[i+1]==n[1]) ntoc[1][1]=2; if(c[i+1]==n[2]) ntoc[2][1]=2;} sc[2][0]=sc[2][0]/c[0]; sc[2][1]=sc[2][1]/c[0]; sc[2][2]=sc[2][2]/c[0];}

closest->nb=0;
if(ntoc[0][0]==ntoc[1][0] || ntoc[0][0]==ntoc[1][1]) GetCrossing_T(closest,p,s,n[0],pointindex,nc[ntoc[0][0]],sc[ntoc[0][0]]); else GetCrossing_T(closest,p,s,pointindex,n[0],nc[ntoc[0][0]],sc[ntoc[0][0]]);  if(closest->nb!=0) return;
if(ntoc[0][1]!=-1 && (ntoc[0][1]==ntoc[1][0] || ntoc[0][1]==ntoc[1][1])) GetCrossing_T(closest,p,s,n[0],pointindex,nc[ntoc[0][1]],sc[ntoc[0][1]]); else GetCrossing_T(closest,p,s,pointindex,n[0],nc[ntoc[0][1]],sc[ntoc[0][1]]); if(closest->nb!=0) return;
if(ntoc[1][0]==ntoc[2][0] || ntoc[1][0]==ntoc[2][1]) GetCrossing_T(closest,p,s,n[1],pointindex,nc[ntoc[1][0]],sc[ntoc[1][0]]); else GetCrossing_T(closest,p,s,pointindex,n[1],nc[ntoc[1][0]],sc[ntoc[1][0]]); if(closest->nb!=0) return;
if(ntoc[1][1]!=-1 && (ntoc[1][1]==ntoc[2][0] || ntoc[1][1]==ntoc[2][1])) GetCrossing_T(closest,p,s,n[1],pointindex,nc[ntoc[1][1]],sc[ntoc[1][1]]); else GetCrossing_T(closest,p,s,pointindex,n[1],nc[ntoc[1][1]],sc[ntoc[1][1]]); if(closest->nb!=0) return;
if(ntoc[2][0]==ntoc[0][0] || ntoc[2][0]==ntoc[0][1]) GetCrossing_T(closest,p,s,n[2],pointindex,nc[ntoc[2][0]],sc[ntoc[2][0]]); else GetCrossing_T(closest,p,s,pointindex,n[2],nc[ntoc[2][0]],sc[ntoc[2][0]]); if(closest->nb!=0) return;
if(ntoc[2][1]!=-1 && (ntoc[2][1]==ntoc[0][0] || ntoc[2][1]==ntoc[0][1])) GetCrossing_T(closest,p,s,n[2],pointindex,nc[ntoc[2][1]],sc[ntoc[2][1]]); else GetCrossing_T(closest,p,s,pointindex,n[2],nc[ntoc[2][1]],sc[ntoc[2][1]]); if(closest->nb!=0) return;
closest->nb=0; 
}


void CSimplexSurf::GetCrossing_C(CLOSESTSTRUCT* closest,P_float p[3],P_float s[3],int cellindex)
{
if(cellindex==-1) return;
int i,p1,p2,*c=GetCell(cellindex);
P_float sc[3]={0,0,0};
for(i=0;i<c[0];i++) {sc[0]+=speeds[3*c[i+1]];sc[1]+=speeds[3*c[i+1]+1];sc[2]+=speeds[3*c[i+1]+2];}
sc[0]=sc[0]/c[0]; sc[1]=sc[1]/c[0]; sc[2]=sc[2]/c[0];
for(i=0;i<c[0];i++) 
    {
    p1=c[i+1]; p2=(i==c[0]-1)?c[1]:c[i+2]; 
    GetCrossing_T(closest,p,s,p1,p2,cellindex,sc); 
    }
}


void CSimplexSurf::GetCrossing_C_ini(CLOSESTSTRUCT* closest,P_float p[3],int cellindex)
{
if(cellindex==-1) return;
int *c=GetCell(cellindex);
for(int i=0;i<c[0];i++) 
    {
    P_float *p1=cellcenters+3*cellindex,*p2=GetPoint(c[i+1]),*p3=GetPoint((i==c[0]-1)?c[1]:c[i+2]);
    P_float AB[3]={p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]}; P_float AC[3]={p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2]}; P_float AP[3]={p[0]-p1[0],p[1]-p1[1],p[2]-p1[2]}; 
    P_float nrml[3]={AB[1]*AC[2] - AB[2]*AC[1] , AB[2]*AC[0] - AB[0]*AC[2] , AB[0]*AC[1] - AB[1]*AC[0]};
    P_float sl;
    if(nrml[0]!=0)  sl=(AP[0]*nrml[0]+ AP[1]*nrml[1] + AP[2]*nrml[2])/nrml[0]; else    sl=-1;
    if(sl>=0)  
        {
        P_float P[3]={p[0]-sl,p[1],p[2]};
        P_float e[3]; Barycenter(e,P,p1,p2,p3);       
        if(e[0]>0 && e[1]>0 &&  e[2]>0 && e[0]<1 && e[1]<1 && e[2]<1)       
            closest->side*=-1;
        }
    }
}

void CSimplexSurf::GetCrossing_T(CLOSESTSTRUCT* closest,P_float p[3],P_float s[3],int pindex1,int pindex2,int cindex,P_float sc[3])
{
int nb_crossing=0; 
P_float *p1=cellcenters+3*cindex,*p2=GetPoint(pindex1),*p3=GetPoint(pindex2),*s1=sc,*s2=GetSpeed(pindex1),*s3=GetSpeed(pindex2);

P_float sl[3]={1E10,1E10,1E10},nrml[3][3];

//provot's methods

// initial positions
double VAB[3]={s2[0]-s1[0],s2[1]-s1[1],s2[2]-s1[2]}; double VAC[3]={s3[0]-s1[0],s3[1]-s1[1],s3[2]-s1[2]}; double VAP[3]={s[0]-s1[0],s[1]-s1[1],s[2]-s1[2]};
double AB[3]={p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]}; double AC[3]={p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2]}; double AP[3]={p[0]-p1[0],p[1]-p1[1],p[2]-p1[2]}; 

if(VAB[0]==0 && VAB[1]==0 && VAB[2]==0 && VAC[0]==0 && VAC[1]==0 && VAC[2]==0)
    {
    if(VAP[0]==0 && VAP[1]==0 && VAP[2]==0) return;
    // n= AB n AC 
    nrml[0][0]=AB[1]*AC[2] - AB[2]*AC[1]; nrml[0][1]=AB[2]*AC[0] - AB[0]*AC[2]; nrml[0][2]=AB[0]*AC[1] - AB[1]*AC[0];
    // dp= AP.n = f(t) + cte
    P_float dp[2]={   AP[0]*nrml[0][0]+ AP[1]*nrml[0][1] + AP[2]*nrml[0][2], VAP[0]*nrml[0][0]+ VAP[1]*nrml[0][1]+ VAP[2]*nrml[0][2]};

    sl[nb_crossing]=-dp[0]/dp[1];
    if(sl[nb_crossing]<=0 && sl[nb_crossing]>=-timestep)  nb_crossing++;
    }

else if(VAP[0]==0 && VAP[1]==0 && VAP[2]==0)
    {
    // n= AB n AC = f(t,t^2) + cte
    double n[3][3]={ AB[1]*AC[2] - AB[2]*AC[1] , VAB[1]*AC[2] + AB[1]*VAC[2] - VAB[2]*AC[1] - AB[2]*VAC[1] ,    VAB[1]*VAC[2] - VAB[2]*VAC[1] , 
                    AB[2]*AC[0] - AB[0]*AC[2] , VAB[2]*AC[0] + AB[2]*VAC[0] - VAB[0]*AC[2] - AB[0]*VAC[2] , VAB[2]*VAC[0] - VAB[0]*VAC[2] ,
                    AB[0]*AC[1] - AB[1]*AC[0] , VAB[0]*AC[1] + AB[0]*VAC[1] - VAB[1]*AC[0] - AB[1]*VAC[0] , VAB[0]*VAC[1] - VAB[1]*VAC[0] };

    // dp= AP.n = f(t,t^2) + cte
    double dp[3]={  AP[0]*n[0][0] + AP[1]*n[1][0] + AP[2]*n[2][0] ,
                    AP[0]*n[0][1] + AP[1]*n[1][1] + AP[2]*n[2][1] ,
                    AP[0]*n[0][2] + AP[1]*n[1][2] + AP[2]*n[2][2] };

    // resolve dp=0 for tE[0,timestep]

    double sol[2]; int num_sol[1],retv; 
    retv=vtkMath::SolveQuadratic(dp[2],dp[1],dp[0],sol,sol+1,num_sol);

    switch(retv)
        {
        case 1:  //(1)-one distinct real root of multiplicity 3 (stored in r1); 
            if(sol[0]<=0 && sol[0]>=-timestep)  {sl[nb_crossing]=sol[0]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++; }
        break;

        case 2:  //(2)-two distinct real roots, one of multiplicity 2 (stored in r1 & r2); 
            if(sol[0]<=0 && sol[0]>=-timestep)  {sl[nb_crossing]=sol[0]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sl[nb_crossing]+n[0][2]*sl[nb_crossing]*sl[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sl[nb_crossing]+n[1][2]*sl[nb_crossing]*sl[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sl[nb_crossing]+n[2][2]*sl[nb_crossing]*sl[nb_crossing]; nb_crossing++;}
            if(sol[1]<=0 && sol[1]>=-timestep)  {sl[nb_crossing]=sol[1]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sl[nb_crossing]+n[0][2]*sl[nb_crossing]*sl[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sl[nb_crossing]+n[1][2]*sl[nb_crossing]*sl[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sl[nb_crossing]+n[2][2]*sl[nb_crossing]*sl[nb_crossing]; nb_crossing++;}
        break;

        default:
            return; // no collision
        break;
        }
    }
else
    {
    // n= AB n AC = f(t,t^2) + cte
    double n[3][3]={ AB[1]*AC[2] - AB[2]*AC[1] , VAB[1]*AC[2] + AB[1]*VAC[2] - VAB[2]*AC[1] - AB[2]*VAC[1] ,    VAB[1]*VAC[2] - VAB[2]*VAC[1] , 
                    AB[2]*AC[0] - AB[0]*AC[2] , VAB[2]*AC[0] + AB[2]*VAC[0] - VAB[0]*AC[2] - AB[0]*VAC[2] , VAB[2]*VAC[0] - VAB[0]*VAC[2] ,
                    AB[0]*AC[1] - AB[1]*AC[0] , VAB[0]*AC[1] + AB[0]*VAC[1] - VAB[1]*AC[0] - AB[1]*VAC[0] , VAB[0]*VAC[1] - VAB[1]*VAC[0] };

    // dp= AP.n = f(t,t^2,t^3) + cte
    double dp[4]={  AP[0]*n[0][0] + AP[1]*n[1][0] + AP[2]*n[2][0] ,
                    AP[0]*n[0][1] + AP[1]*n[1][1] + AP[2]*n[2][1] + VAP[0]*n[0][0] + VAP[1]*n[1][0] + VAP[2]*n[2][0] ,
                    AP[0]*n[0][2] + AP[1]*n[1][2] + AP[2]*n[2][2] + VAP[0]*n[0][1] + VAP[1]*n[1][1] + VAP[2]*n[2][1] ,
                    VAP[0]*n[0][2] + VAP[1]*n[1][2] + VAP[2]*n[2][2] };

    // resolve dp=0 for tE[0,timestep]

    double sol[3]; int num_sol[1],retv; 
    retv=vtkMath::SolveCubic(dp[3],dp[2],dp[1],dp[0],sol,sol+1,sol+2,num_sol);

    switch(retv)
        {
        case 1:  //(1)-one distinct real root of multiplicity 3 (stored in r1); 
            if(sol[0]<=0 && sol[0]>=-timestep)  {sl[nb_crossing]=sol[0]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++; }
        break;

        case 2:  //(2)-two distinct real roots, one of multiplicity 2 (stored in r1 & r2); 
            if(sol[0]<=0 && sol[0]>=-timestep)  {sl[nb_crossing]=sol[0]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++;}
            if(sol[1]<=0 && sol[1]>=-timestep)  {sl[nb_crossing]=sol[1]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++;}
        break;

        case 3:  //(3)-three distinct real roots;
            if(sol[0]<=0 && sol[0]>=-timestep)  {sl[nb_crossing]=sol[0]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++;}
            if(sol[1]<=0 && sol[1]>=-timestep)  {sl[nb_crossing]=sol[1]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++;}
            if(sol[2]<=0 && sol[2]>=-timestep)  {sl[nb_crossing]=sol[2]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++;}
        break;

        case -3:    //(-3)-one real root and a complex conjugate pair (real root in r1 and real part of pair in r2 and imaginary in r3). 
            if(sol[0]<=0 && sol[0]>=-timestep)  {sl[nb_crossing]=sol[0]; nrml[nb_crossing][0]=n[0][0]+n[0][1]*sol[nb_crossing]+n[0][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][1]=n[1][0]+n[1][1]*sol[nb_crossing]+n[1][2]*sol[nb_crossing]*sol[nb_crossing]; nrml[nb_crossing][2]=n[2][0]+n[2][1]*sol[nb_crossing]+n[2][2]*sol[nb_crossing]*sol[nb_crossing]; nb_crossing++; }
        break;

        default:
            return; // no collision
        break;
        }
    }

for(int i=0;i<nb_crossing;i++) GetCrossing_T2(closest,sl[i],nrml[i],p,s,pindex1,pindex2,cindex,p1,s1,p2,s2,p3,s3);
}

void CSimplexSurf::GetCrossing_T2(CLOSESTSTRUCT* closest,P_float sl,P_float n[3],P_float p[3],P_float s[3],int pindex1,int pindex2,int cindex,P_float p1[3],P_float s1[3],P_float p2[3],P_float s2[3],P_float p3[3], P_float s3[3])
{
P_float P[3]={p[0]+sl*s[0],p[1]+sl*s[1],p[2]+sl*s[2]};
P_float A[3]={p1[0]+sl*s1[0],p1[1]+sl*s1[1],p1[2]+sl*s1[2]};
P_float B[3]={p2[0]+sl*s2[0],p2[1]+sl*s2[1],p2[2]+sl*s2[2]};
P_float C[3]={p3[0]+sl*s3[0],p3[1]+sl*s3[1],p3[2]+sl*s3[2]};

P_float e[3]; Barycenter(e,P,A,B,C);
//e[2]=1-e[1]-e[0];

// find t,s / ap = s.ab + t.ac  p=sb+tc+ua

bool copy=false;
if(e[0]>0 && e[1]>0 &&  e[2]>0 && e[0]<1 && e[1]<1 && e[2]<1)
    {
    closest->side*=-1;
    closest->crossing=true;
    if(sl<closest->t-timestep)  copy=true;

    if(copy)
        {
        int *c=GetCell(cindex);
        free(closest->weights); free(closest->pts);
        closest->nb=c[0];       closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb]; 
        for(int i=0;i<c[0];i++) {closest->weights[i]=e[0]/c[0]; closest->pts[i]=c[i+1]; if(c[i+1]==pindex1) closest->weights[i]+=e[1]; if(c[i+1]==pindex2) closest->weights[i]+=e[2];}
        closest->c[0]=e[0]*p1[0]+e[1]*p2[0]+e[2]*p3[0];     closest->c[1]=e[0]*p1[1]+e[1]*p2[1]+e[2]*p3[1];         closest->c[2]=e[0]*p1[2]+e[1]*p2[2]+e[2]*p3[2];
        closest->dist2=abs((p[0]-closest->c[0])*(p[0]-closest->c[0])+(p[1]-closest->c[1])*(p[1]-closest->c[1])+(p[2]-closest->c[2])*(p[2]-closest->c[2]));
        closest->t=sl+timestep;
		closest->n[0]=0; closest->n[1]=0; closest->n[2]=0; 
		for(int j=0;j<closest->nb;j++) {closest->n[0]+=normals[3*closest->pts[j]]*closest->weights[j]; closest->n[1]+=normals[3*closest->pts[j]+1]*closest->weights[j]; closest->n[2]+=normals[3*closest->pts[j]+2]*closest->weights[j]; }
        }
    }
}


void CSimplexSurf::GetClosest_C(CLOSESTSTRUCT* closest,P_float p[3],int cellindex) 
{
if(cellindex==-1) return;
int i,j,*c;
P_float *pi,*cc=cellcenters+3*cellindex;

c=GetCell(cellindex);
P_float *ccpi=new P_float[3*c[0]]; 
P_float *ppi=new P_float[3*c[0]];
P_float *n=new P_float[3*c[0]];
P_float pcc[3]={cc[0]-p[0],cc[1]-p[1],cc[2]-p[2]};

for(i=0;i<c[0];i++) {pi=GetPoint(c[i+1]); ccpi[3*i]=pi[0]-cc[0]; ccpi[3*i+1]=pi[1]-cc[1]; ccpi[3*i+2]=pi[2]-cc[2]; ppi[3*i]=pi[0]-p[0]; ppi[3*i+1]=pi[1]-p[1]; ppi[3*i+2]=pi[2]-p[2];}

P_float *fa=new P_float[c[0]];
P_float *fb=new P_float[c[0]];
P_float *fd=new P_float[c[0]];
P_float ff=dotproduct(pcc,pcc);

for(i=0;i<c[0];i++) {fa[i]=dotproduct(ccpi+3*i,ccpi+3*i); if(i!=c[0]-1) {crossproduct(n+3*i,ccpi+3*i,ccpi+3*i+3); fb[i]=dotproduct(ccpi+3*i,ccpi+3*i+3);} else   {fb[i]=dotproduct(ccpi+3*i,ccpi); crossproduct(n+3*i,ccpi+3*i,ccpi);} fd[i]=dotproduct(ccpi+3*i,pcc); }

P_float fout[2],fdist,dist_sav=closest->dist2,t_sav,s_sav;
int sav=-1;

for(i=0;i<c[0];i++)
    {
    if(i!=c[0]-1) 	fdist=Closest(fa[i],fb[i],fa[i+1],fd[i],fd[i+1],ff,fout);
    else  fdist=Closest(fa[i],fb[i],fa[0],fd[i],fd[0],ff,fout);
	if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=i;    }
	}

//if(_finite(dist_sav)==0 || sav==-1) {free(ccpi); free(ppi); free(n); free(fa); free(fb); free(fd); return;}

// return values
free(closest->pts); free(closest->weights);

int sav2=(sav!=c[0]-1)?(sav+1):0;
int count=0;
P_float u_sav=1-s_sav-t_sav;

closest->cell=cellindex;
closest->dist2=dist_sav;
closest->c[0]=cc[0]+s_sav*ccpi[3*sav]+t_sav*ccpi[3*sav2];        closest->c[1]=cc[1]+s_sav*ccpi[3*sav+1]+t_sav*ccpi[3*sav2+1];          closest->c[2]=cc[2]+s_sav*ccpi[3*sav+2]+t_sav*ccpi[3*sav2+2];

if(u_sav<=0) 
    {
    closest->nb=(t_sav==0)?(0):(1);
    closest->nb+=(s_sav==0)?(0):(1); 
    closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb]; 
    if(s_sav!=0) {closest->weights[0]=s_sav; closest->pts[0]=c[sav+1]; count++;}
    if(t_sav!=0) {closest->weights[count]=t_sav; closest->pts[count]=c[sav2+1]; }
    }
else 
    {
    closest->nb=c[0]; closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb];
    for(j=0;j<closest->nb;j++) { closest->weights[j]=u_sav/closest->nb; closest->pts[j]=c[j+1];    }
    closest->weights[sav]+=s_sav;
    closest->weights[sav2]+=t_sav;
	}   

closest->n[0]=0; closest->n[1]=0; closest->n[2]=0; 
for(j=0;j<closest->nb;j++) {closest->n[0]+=normals[3*closest->pts[j]]*closest->weights[j]; closest->n[1]+=normals[3*closest->pts[j]+1]*closest->weights[j]; closest->n[2]+=normals[3*closest->pts[j]+2]*closest->weights[j]; }
if(dotproduct(ppi+3*sav,closest->n)<=0) closest->side=1; else closest->side=-1;

free(ccpi); free(ppi); free(n); free(fa); free(fb); free(fd);
return;
}


void CSimplexSurf::GetClosest_P(CLOSESTSTRUCT* closest,P_float p[3],int pointindex) 
{
int j;
P_float *cp,temp[3]; 
cp=GetPoint(pointindex);

// neighbors
int *n; n=GetNeighbors(pointindex);  
P_float *n1,*n2,*n3; n1=GetPoint(n[0]); n2=GetPoint(n[1]); n3=GetPoint(n[2]);

P_float cpn1[3]={n1[0]-cp[0],n1[1]-cp[1],n1[2]-cp[2]};P_float cpn2[3]={n2[0]-cp[0],n2[1]-cp[1],n2[2]-cp[2]};P_float cpn3[3]={n3[0]-cp[0],n3[1]-cp[1],n3[2]-cp[2]};
P_float pcp[3]={cp[0]-p[0],cp[1]-p[1],cp[2]-p[2]};
P_float FA[3]={dotproduct(cpn1,cpn1),dotproduct(cpn2,cpn2),dotproduct(cpn3,cpn3)};
P_float ff=dotproduct(pcp,pcp);
P_float FD[3]={dotproduct(cpn1,pcp),dotproduct(cpn2,pcp),dotproduct(cpn3,pcp)};

// cell centers
int *nc; nc=GetNeighbors_c(pointindex); 
int ntoc[3][2]={-1,-1,-1,-1,-1,-1};
for(j=0;j<cells[nc[0]][0];j++) { if(cells[nc[0]][j+1]==n[0]) ntoc[0][0]=0; if(cells[nc[0]][j+1]==n[1]) ntoc[1][0]=0;    if(cells[nc[0]][j+1]==n[2]) ntoc[2][0]=0;   }
for(j=0;j<cells[nc[1]][0];j++) { if(cells[nc[1]][j+1]==n[0]) ntoc[0][(ntoc[0][0]==-1)?0:1]=1;   if(cells[nc[1]][j+1]==n[1]) ntoc[1][(ntoc[1][0]==-1)?0:1]=1;    if(cells[nc[1]][j+1]==n[2]) ntoc[2][(ntoc[2][0]==-1)?0:1]=1;    }
if(nc[2]!=-1) for(j=0;j<cells[nc[2]][0];j++) { if(cells[nc[2]][j+1]==n[0]) ntoc[0][1]=2; if(cells[nc[2]][j+1]==n[1]) ntoc[1][1]=2;    if(cells[nc[2]][j+1]==n[2]) ntoc[2][1]=2;   GetPoint(cells[nc[2]][j+1],temp);       }

P_float cpcc[3][3],FB[3][2];
cpcc[0][0]=cellcenters[3*nc[0]]-cp[0];		cpcc[0][1]=cellcenters[3*nc[0]+1]-cp[1];		cpcc[0][2]=cellcenters[3*nc[0]+2]-cp[2];
cpcc[1][0]=cellcenters[3*nc[1]]-cp[0];		cpcc[1][1]=cellcenters[3*nc[1]+1]-cp[1];		cpcc[1][2]=cellcenters[3*nc[1]+2]-cp[2];
if(nc[2]!=-1) {cpcc[2][0]=cellcenters[3*nc[2]]-cp[0];		cpcc[2][1]=cellcenters[3*nc[2]+1]-cp[1];		cpcc[2][2]=cellcenters[3*nc[2]+2]-cp[2];}
else {cpcc[2][0]=0; cpcc[2][1]=0; cpcc[2][2]=0;}

P_float FC[3]={dotproduct(cpcc[0],cpcc[0]),dotproduct(cpcc[1],cpcc[1]),dotproduct(cpcc[2],cpcc[2])};
P_float FE[3]={dotproduct(cpcc[0],pcp),dotproduct(cpcc[1],pcp),dotproduct(cpcc[2],pcp)};
FB[0][0]=dotproduct(cpn1,cpcc[ntoc[0][0]]);		if(ntoc[0][1]!=-1) FB[0][1]=dotproduct(cpn1,cpcc[ntoc[0][1]]);
FB[1][0]=dotproduct(cpn2,cpcc[ntoc[1][0]]);		if(ntoc[1][1]!=-1) FB[1][1]=dotproduct(cpn2,cpcc[ntoc[1][1]]);
FB[2][0]=dotproduct(cpn3,cpcc[ntoc[2][0]]);		if(ntoc[2][1]!=-1) FB[2][1]=dotproduct(cpn3,cpcc[ntoc[2][1]]);

// test proj for each triangle
P_float fa,fb,fc,fd,fe,fdist,fout[2];
P_float s_sav,t_sav,dist_sav=closest->dist2;
int sav=-1;

fa=FA[0]; fb=FB[0][0]; fc=FC[ntoc[0][0]]; fd=FD[0]; fe=FE[ntoc[0][0]];
fdist=Closest(fa,fb,fc,fd,fe,ff,fout);
if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=0;	}

if(ntoc[0][1]!=-1) 
{
fa=FA[0]; fb=FB[0][1]; fc=FC[ntoc[0][1]]; fd=FD[0]; fe=FE[ntoc[0][1]];
fdist=Closest(fa,fb,fc,fd,fe,ff,fout);
if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=1;    }
}

fa=FA[1]; fb=FB[1][0]; fc=FC[ntoc[1][0]]; fd=FD[1]; fe=FE[ntoc[1][0]];
fdist=Closest(fa,fb,fc,fd,fe,ff,fout);
if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=2;    }

if(ntoc[1][1]!=-1) 
{
fa=FA[1]; fb=FB[1][1]; fc=FC[ntoc[1][1]]; fd=FD[1]; fe=FE[ntoc[1][1]];
fdist=Closest(fa,fb,fc,fd,fe,ff,fout);
if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=3;    }
}

fa=FA[2]; fb=FB[2][0]; fc=FC[ntoc[2][0]]; fd=FD[2]; fe=FE[ntoc[2][0]];
fdist=Closest(fa,fb,fc,fd,fe,ff,fout);
if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=4;    }

if(ntoc[2][1]!=-1) 
{
fa=FA[2]; fb=FB[2][1]; fc=FC[ntoc[2][1]]; fd=FD[2]; fe=FE[ntoc[2][1]];
fdist=Closest(fa,fb,fc,fd,fe,ff,fout);
if(fdist<dist_sav)  {   dist_sav=fdist; s_sav=fout[0]; t_sav=fout[1]; sav=5;    }
}

//if(_finite(dist_sav)==0 || sav==-1) {return;}

// return values
free(closest->pts); free(closest->weights);

int *cell,count=0,n_sav;
P_float u_sav=1-s_sav-t_sav;

switch(sav)
    {
    case 0:
    n_sav=n[0]; closest->cell=nc[ntoc[0][0]]; 
    closest->c[0]=cp[0]+s_sav*cpn1[0]+t_sav*cpcc[ntoc[0][0]][0]; closest->c[1]=cp[1]+s_sav*cpn1[1]+t_sav*cpcc[ntoc[0][0]][1]; closest->c[2]=cp[2]+s_sav*cpn1[2]+t_sav*cpcc[ntoc[0][0]][2];
    break;

    case 1:
    n_sav=n[0]; closest->cell=nc[ntoc[0][1]];
    closest->c[0]=cp[0]+s_sav*cpn1[0]+t_sav*cpcc[ntoc[0][1]][0]; closest->c[1]=cp[1]+s_sav*cpn1[1]+t_sav*cpcc[ntoc[0][1]][1]; closest->c[2]=cp[2]+s_sav*cpn1[2]+t_sav*cpcc[ntoc[0][1]][2];
    break;

    case 2:
    n_sav=n[1]; closest->cell=nc[ntoc[1][0]];
    closest->c[0]=cp[0]+s_sav*cpn2[0]+t_sav*cpcc[ntoc[1][0]][0]; closest->c[1]=cp[1]+s_sav*cpn2[1]+t_sav*cpcc[ntoc[1][0]][1]; closest->c[2]=cp[2]+s_sav*cpn2[2]+t_sav*cpcc[ntoc[1][0]][2];
    break;

    case 3:
    n_sav=n[1]; closest->cell=nc[ntoc[1][1]];
    closest->c[0]=cp[0]+s_sav*cpn2[0]+t_sav*cpcc[ntoc[1][1]][0]; closest->c[1]=cp[1]+s_sav*cpn2[1]+t_sav*cpcc[ntoc[1][1]][1]; closest->c[2]=cp[2]+s_sav*cpn2[2]+t_sav*cpcc[ntoc[1][1]][2];
    break;

    case 4:
    n_sav=n[2]; closest->cell=nc[ntoc[2][0]];
    closest->c[0]=cp[0]+s_sav*cpn3[0]+t_sav*cpcc[ntoc[2][0]][0]; closest->c[1]=cp[1]+s_sav*cpn3[1]+t_sav*cpcc[ntoc[2][0]][1]; closest->c[2]=cp[2]+s_sav*cpn3[2]+t_sav*cpcc[ntoc[2][0]][2];
    break;

    case 5:
    n_sav=n[2]; closest->cell=nc[ntoc[2][1]];
    closest->c[0]=cp[0]+s_sav*cpn3[0]+t_sav*cpcc[ntoc[2][1]][0]; closest->c[1]=cp[1]+s_sav*cpn3[1]+t_sav*cpcc[ntoc[2][1]][1]; closest->c[2]=cp[2]+s_sav*cpn3[2]+t_sav*cpcc[ntoc[2][1]][2];
    break;
    }

if(t_sav==0) 
    {
    closest->nb=(u_sav==0)?(0):(1);
    closest->nb+=(s_sav==0)?(0):(1); 
    closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb]; 
    if(s_sav!=0) {closest->weights[0]=s_sav; closest->pts[0]=n_sav; count++;}
    if(u_sav!=0) {closest->weights[count]=u_sav; closest->pts[count]=pointindex; }
    }
else 
    {
    cell=cells[closest->cell];
    closest->nb=cell[0];
    closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb];
    for(j=0;j<closest->nb;j++) { closest->weights[j]=t_sav/closest->nb; closest->pts[j]=cell[j+1];    }
    if(s_sav!=0) for(j=0;j<closest->nb;j++) if(cell[j+1]==n_sav) closest->weights[j]+=s_sav; 
    if(u_sav!=0) for(j=0;j<closest->nb;j++) if(cell[j+1]==pointindex) closest->weights[j]+=u_sav; 
    }   

closest->n[0]=0; closest->n[1]=0; closest->n[2]=0; 
for(j=0;j<closest->nb;j++) {closest->n[0]+=normals[3*closest->pts[j]]*closest->weights[j]; closest->n[1]+=normals[3*closest->pts[j]+1]*closest->weights[j]; closest->n[2]+=normals[3*closest->pts[j]+2]*closest->weights[j]; }
if(dotproduct(pcp,closest->n)<=0) closest->side=1; else closest->side=-1;
closest->dist2=dist_sav;
return;
}


void CSimplexSurf::RemoveAllAttachPoints()
{
for(int i=0;i<nb_points;i++) {attached_point[i]=0; SetMass(i);}
mass_modified=false;
}

void CSimplexSurf::AttachBorder()
{
for(int i=0;i<nb_points;i++) if(neighbors_c[3*i+2]==-1) SetAttachPoint(i,points+3*i);
}


void CSimplexSurf::SetAttachedPoint(int index,int val)
{
attached_point[index]=val;
attached_cell[neighbors_c[3*index]]=abs(val);
attached_cell[neighbors_c[3*index+1]]=abs(val);
if(neighbors_c[3*index+2]!=-1) attached_cell[neighbors_c[3*index+2]]=abs(val);
}

void CSimplexSurf::SetAttachPoint(int index,P_float pa[3])
{
if(attached_point[index]==0) attached_point[index]=1;
memcpy(points+3*index,pa,3*sizeof(P_float));
SetMass(index,-1);
}

void CSimplexSurf::SetAttachPoint(int index,P_float x, P_float y, P_float z)
{
if(attached_point[index]==0) attached_point[index]=1;
points[3*index]=x; points[3*index+1]=y; points[3*index+2]=z;
SetMass(index,-1);
}

void CSimplexSurf::SetAttachPoint(int index,P_float pa[3],P_float n[3]) // attach point to a surface
{
if(attached_point[index]==0) attached_point[index]=1;
memcpy(points+3*index,pa,3*sizeof(P_float));
//SetMass_planeconstraint(index,n);
SetMass(index);
}


int CSimplexSurf::GetAttachPoint(int index,P_float pa[3])
{
memcpy(pa,points+3*index,3*sizeof(P_float));
return attached_point[index];
}

int CSimplexSurf::GetAttachedPoint(int index)	{return attached_point[index];}

void CSimplexSurf::SetAttachCell(int index,int val) {attached_cell[index]=val;}
int CSimplexSurf::GetAttachCell(int index) {return attached_cell[index];}


int* CSimplexSurf::UpdateAttachedPoints(int attachment_index,P_float* seed_point,P_float* seed_point2)
{
int i,j,k,n[3];
int* ret=NULL;

// detect and order contour points according to attached cells
int count=0,*cont=new int[nb_points],*flag=new int[nb_points],*selected=new int[nb_points]; for(i=0;i<nb_points;i++) {cont[i]=-1; flag[i]=0;}

for(i=0;i<nb_points;i++)	if(attached_point[i]==attachment_index) attached_point[i]=0;

// flag: number of attached cell per point
for(i=0;i<nb_cells;i++)
    if(attached_cell[i]==attachment_index)
        for(j=0;j<cells[i][0];j++)
            flag[cells[i][j+1]]++;

// flag=1 or 2 -> contour point. targetcount number of contour point 
int targetcount=0;
for(i=0;i<nb_points;i++)    if(flag[i]==1 || flag[i]==2) targetcount++;
if(targetcount==0) return NULL;

// handle merged splines
	int seed_point_index1=-1,seed_point_index2=-1; P_float d,dmin;
	if(seed_point2!=NULL) 
		{
		dmin=1E10; for(i=0;i<nb_points;i++) if(flag[i]==1) {d=(seed_point[0]-points[3*i])*(seed_point[0]-points[3*i])+(seed_point[1]-points[3*i+1])*(seed_point[1]-points[3*i+1])+(seed_point[2]-points[3*i+2])*(seed_point[2]-points[3*i+2]); if(d<dmin) {dmin=d; seed_point_index1=i;}}
		dmin=1E10; for(i=0;i<nb_points;i++) if(flag[i]==1) {d=(seed_point2[0]-points[3*i])*(seed_point2[0]-points[3*i])+(seed_point2[1]-points[3*i+1])*(seed_point2[1]-points[3*i+1])+(seed_point2[2]-points[3*i+2])*(seed_point2[2]-points[3*i+2]); if(d<dmin) {dmin=d; seed_point_index2=i;}}
		}

i=0;
while(count==0) 
    {
    for(k=0;k<nb_points;k++) selected[k]=0; // selected controls that points are selected once

    // search seed cell and point
	if(seed_point2!=NULL) {cont[0]=seed_point_index1; selected[seed_point_index1]=1; count=1; }// handle merged splines
	else
		{
		while(count==0 || i>=nb_cells) 
			if(attached_cell[i]==attachment_index)
				{
				for(j=1;j<cells[i][0]+1;j++) if(flag[cells[i][j]]==1)   { cont[0]=cells[i][j]; selected[cont[0]]=1; count=1; break;}
				i++;
				if(i>=nb_cells) {count=0; break;}
				}
			else { i++; if(i>=nb_cells) {count=0; break;}}
		if(i>=nb_cells) break;
		}
    // second point
    memcpy(n,neighbors+3*cont[count-1],3*sizeof(int));
    if(flag[n[0]]==0) {cont[count]=n[2]; selected[cont[count]]=1; count++;}
    if(flag[n[1]]==0) {cont[count]=n[0]; selected[cont[count]]=1; count++;}
    if(flag[n[2]]==0) {cont[count]=n[1]; selected[cont[count]]=1; count++;}

    // follow the loop
    for(k=0;k<targetcount-2;k++)
        {
        memcpy(n,neighbors+3*cont[count-1],3*sizeof(int));

        if(n[0]==cont[count-2]) {if(selected[n[2]]==0 && (flag[n[2]]==1 || flag[n[2]]==2)) {cont[count]=n[2]; selected[cont[count]]=1; count++;} else if(selected[n[1]]==0 && (flag[n[1]]==1 || flag[n[1]]==2)) {cont[count]=n[1]; selected[cont[count]]=1; count++;} }
        else if(n[1]==cont[count-2]) {if(selected[n[0]]==0 && (flag[n[0]]==1 || flag[n[0]]==2)) {cont[count]=n[0]; selected[cont[count]]=1; count++;} else if(selected[n[2]]==0 && (flag[n[2]]==1 || flag[n[2]]==2)) {cont[count]=n[2]; selected[cont[count]]=1; count++;} }
        else if(n[2]==cont[count-2]) {if(selected[n[1]]==0 && (flag[n[1]]==1 || flag[n[1]]==2)) {cont[count]=n[1]; selected[cont[count]]=1; count++;} else if(selected[n[0]]==0 && (flag[n[0]]==1 || flag[n[0]]==2)) {cont[count]=n[0]; selected[cont[count]]=1; count++;}}
        }

	if(seed_point2!=NULL)// handle merged splines
		{
		cont[count]=seed_point_index2; selected[cont[count]]=1; count++;
		memcpy(n,neighbors+3*cont[count-1],3*sizeof(int));
		if(flag[n[0]]==0) {cont[count]=n[2]; selected[cont[count]]=1; count++;}
		if(flag[n[1]]==0) {cont[count]=n[0]; selected[cont[count]]=1; count++;}
		if(flag[n[2]]==0) {cont[count]=n[1]; selected[cont[count]]=1; count++;}
		for(k=0;k<targetcount-2;k++)
			{
			memcpy(n,neighbors+3*cont[count-1],3*sizeof(int));

			if(n[0]==cont[count-2]) {if(selected[n[2]]==0 && (flag[n[2]]==1 || flag[n[2]]==2)) {cont[count]=n[2]; selected[cont[count]]=1; count++;} else if(selected[n[1]]==0 && (flag[n[1]]==1 || flag[n[1]]==2)) {cont[count]=n[1]; selected[cont[count]]=1; count++;} }
			else if(n[1]==cont[count-2]) {if(selected[n[0]]==0 && (flag[n[0]]==1 || flag[n[0]]==2)) {cont[count]=n[0]; selected[cont[count]]=1; count++;} else if(selected[n[2]]==0 && (flag[n[2]]==1 || flag[n[2]]==2)) {cont[count]=n[2]; selected[cont[count]]=1; count++;} }
			else if(n[2]==cont[count-2]) {if(selected[n[1]]==0 && (flag[n[1]]==1 || flag[n[1]]==2)) {cont[count]=n[1]; selected[cont[count]]=1; count++;} else if(selected[n[0]]==0 && (flag[n[0]]==1 || flag[n[0]]==2)) {cont[count]=n[0]; selected[cont[count]]=1; count++;}}
			}
		}
    if((double)count>(double)targetcount/2) break; else count=0; // check if the main loop has been selected
    }
if(targetcount!=count) 
    targetcount=count;

ret=new int[count+1]; ret[0]=count;
//memcpy(ret+1,cont+i,(count-i)*sizeof(int)); memcpy(ret+1+(count-i),cont,i*sizeof(int)); // ??
memcpy(ret+1,cont,count*sizeof(int));

for(i=0;i<count;i++) attached_point[cont[i]]=attachment_index; // contour points
for(i=0;i<nb_cells;i++) if(attached_cell[i]==attachment_index) for(j=0;j<cells[i][0];j++) if(attached_point[cells[i][j+1]]!=attachment_index) attached_point[cells[i][j+1]]=-attachment_index; // interior points

free(cont);free(flag);free(selected);
/*
// update influenced points (tendons)
for(i=0;i<nb_cells;i++)
    if(attached_cell[i]==-1)
        for(j=0;j<cells[i][0];j++)
            if(attached_point[cells[i][j+1]]==0)
                attached_point[cells[i][j+1]]=-1;*/

return ret;
}



int* CSimplexSurf::SelectAttachCellsCentroid(int attachment_index)
{
int i,j,*flagcell=new int[nb_cells],*ret; 
bool stop=false,stop2,flag;

// initialize distance to 0
for(i=0;i<nb_cells;i++) if(attached_cell[i]==attachment_index) flagcell[i]=0; else flagcell[i]=-1; 

// distance map
while(!stop)
	{
	stop2=true;
	for(i=0;i<nb_cells;i++) 
		if(flagcell[i]!=-1)
			{
			flag=false;
			for(j=0;j<cells[i][0];j++)
				{
				if(flagcell[neighbors_c[3*cells[i][j+1]]]<flagcell[i]) flag=true;
				if(flagcell[neighbors_c[3*cells[i][j+1]+1]]<flagcell[i]) flag=true;
				if(flagcell[neighbors_c[3*cells[i][j+1]+2]]<flagcell[i]) flag=true;
				}
			if(!flag) {flagcell[i]++; 	stop2=false;}
			}
	stop=stop2;
	}

// local maxima
int count=0;
for(i=0;i<nb_cells;i++) 
	if(flagcell[i]!=-1)
		{
		flag=false;
		for(j=0;j<cells[i][0];j++)
			{
			if(flagcell[neighbors_c[3*cells[i][j+1]]]>flagcell[i]) flag=true;
			if(flagcell[neighbors_c[3*cells[i][j+1]+1]]>flagcell[i]) flag=true;
			if(flagcell[neighbors_c[3*cells[i][j+1]+2]]>flagcell[i]) flag=true;
			}
		if(!flag) count++;
		}
ret=new int[count+1];
ret[0]=count;
count=0;
for(i=0;i<nb_cells;i++) 
	if(flagcell[i]!=-1)
		{
		flag=false;
		for(j=0;j<cells[i][0];j++)
			{
			if(flagcell[neighbors_c[3*cells[i][j+1]]]>flagcell[i]) flag=true;
			if(flagcell[neighbors_c[3*cells[i][j+1]+1]]>flagcell[i]) flag=true;
			if(flagcell[neighbors_c[3*cells[i][j+1]+2]]>flagcell[i]) flag=true;
			}
		if(!flag) {ret[count+1]=i; count++;}
		}

free(flagcell);
return ret;
}

void CSimplexSurf::SelectCellsWithSpline(P_float* spline,int splineres,int attachment_index,int influence_dist)
{
int selected;
int i,j,flag,stop=0,*flagcell=new int[nb_cells]; 
int seed_cell=-1;

// unattach cells
for(i=0;i<nb_cells;i++)	if(attached_cell[i]==attachment_index) attached_cell[i]=0;

//select contour cells
    selected=GetClosestPoint(spline);

/*	
	int pt;
	for(i=1;i<splineres;i++)
        {
        pt=GetClosestPoint(spline+3*i);
		attached_cell[neighbors_c[3*pt]]=attachment_index;
		attached_cell[neighbors_c[3*pt+1]]=attachment_index;
		attached_cell[neighbors_c[3*pt+2]]=attachment_index;
		}
*/
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1]; 
for(i=1;i<splineres;i++)
       {
		closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10;
        GetClosest(closest,spline+3*i);
		attached_cell[closest->cell]=attachment_index;
		free(closest->pts); free(closest->weights);
	   }


/*
		if(selected!=pt)
            {
            if(SearchCell(selected,pt,c)==1) SetAttachCell(c[1][0],attachment_index);
            else
                {
                k=1; while(cells[c[1][0]][k]!=selected) k++;
                if(cells[c[1][0]][(k==cells[c[1][0]][0])?1:(k+1)]==pt)  
                    {
                    SetAttachCell(c[1][0],attachment_index); 
                    if(seed_cell==-1) 
                        seed_cell=c[1][1]; 
                    } 
                else 
                    {
                    SetAttachCell(c[1][1],attachment_index); 
                    if(seed_cell==-1) 
                        seed_cell=c[1][0]; 
                    }
                }
            selected=pt;
            }*/
        
bool ok=false; seed_cell=0;
int count1,count2;
while(!ok)
	{
	// select seed cell
	for(i=seed_cell+1;i<nb_cells;i++) {if(attached_cell[i]!=attachment_index) seed_cell=i; break;}
	if(seed_cell==nb_cells-1)
		return;

	//Fill Contour Cells
    for(i=0;i<nb_cells;i++) if(attached_cell[i]==attachment_index) flagcell[i]=1; else flagcell[i]=-1;
    flagcell[seed_cell]=0;

    // propagate from seed cell
    while(stop==0)
        {
        flag=1;
        for(i=0;i<nb_cells;i++)
            if(flagcell[i]==0)
                {
                for(j=0;j<cells[i][0];j++)
                    {
                    if(flagcell[neighbors_c[3*cells[i][j+1]]]==-1) {flagcell[neighbors_c[3*cells[i][j+1]]]=0; flag=0;}
                    if(flagcell[neighbors_c[3*cells[i][j+1]+1]]==-1) {flagcell[neighbors_c[3*cells[i][j+1]+1]]=0; flag=0;}
                    if(flagcell[neighbors_c[3*cells[i][j+1]+2]]==-1) {flagcell[neighbors_c[3*cells[i][j+1]+2]]=0; flag=0;}
                    }
                }
        stop=flag;
        }

	// count nb of selected cells
	count1=0; count2=0;
	for(i=0;i<nb_cells;i++) if(flagcell[i]==0) count1++; else if(flagcell[i]==-1) count2++;
	if(count1>count2) ok=true;
	}

    // set untouched cells to attachment_index
    for(i=0;i<nb_cells;i++) if(flagcell[i]==-1) attached_cell[i]=attachment_index;

/*
// set influenced cells to -1 (tendons)
for(i=0;i<nb_cells;i++) if(attached_cell[i]==attachment_index) flagcell[i]=1; else flagcell[i]=-1;

for(k=0;k<influence_dist;k++)
    {
    for(i=0;i<nb_cells;i++)
        if(flagcell[i]==1) 
            for(j=0;j<cells[i][0];j++)
                {
                GetNeighbors_c(cells[i][j+1],nc);
                if(attached_cell[nc[0]]==0) {attached_cell[nc[0]]=-1; flagcell[nc[0]]=2;}
                if(attached_cell[nc[1]]==0) {attached_cell[nc[1]]=-1; flagcell[nc[1]]=2;}
                if(attached_cell[nc[2]]==0) {attached_cell[nc[2]]=-1; flagcell[nc[2]]=2;}
                }
    for(i=0;i<nb_cells;i++) if(flagcell[i]==2) flagcell[i]=1;
    }*/
free(flagcell);
}

void CSimplexSurf::UpdateAttachedFromHigherRes(CSimplexSurf* mesh)
{
int val; P_float pa[3];
for(int i=0;i<nb_points;i++)    
    {
    val=mesh->GetAttachPoint(i,pa);
	attached_point[i]=val;
	if(val!=0)
        {
		if(*mesh->GetMass_inv(i)==0) 	SetMass(i,-1);
		else SetMass(i);
		memcpy(points+3*i,mesh->GetPoint(i),3*sizeof(P_float)); 
        }
    }
}

void CSimplexSurf::SetSlideAttach(bool val) 
{
SlideAttach=val;
if(val)
    {for(int i=0;i<nb_points;i++) if(attached_point[i]<0) SetMass(i);}
else 
    {for(int i=0;i<nb_points;i++) if(attached_point[i]<0) SetMass(i,-1);}
}

void CSimplexSurf::SetMergedPoints(int* mp) {MergedPoints=mp;}
int* CSimplexSurf::GetMergedPoints() {return MergedPoints;}

/////////////////////////////////////////////////////////////////////////////////////////////////
AXIALLINK* CSimplexSurf::GetAxialLinks() {return axiallinks;}
void CSimplexSurf::AllocateAxialLinks() 
{
if(axiallinks!=NULL) {for(int i=0;i<nb_points;i++) {if(axiallinks[i].w!=NULL) free(axiallinks[i].w); if(axiallinks[i].pts!=NULL) free(axiallinks[i].pts);} free(axiallinks);}
axiallinks=new AXIALLINK[nb_points];
for(int i=0;i<nb_points;i++) {axiallinks[i].pts=NULL; axiallinks[i].w=NULL; axiallinks[i].f[0]=0; axiallinks[i].f[1]=0; axiallinks[i].f[2]=0;}
}

void CSimplexSurf::SetAxis_model(CSimplexSurf* model) {Axis_model=model;}
CSimplexSurf* CSimplexSurf::GetAxis_model() {return Axis_model;}



void CSimplexSurf::GetAxialLinks_PaxisFromPmodel(int index,P_float p[3])
{
p[0]=0; p[1]=0; p[2]=0;
AXIALLINK *modellink=Axis_model->GetAxialLinks()+index;
for(int j=0;j<modellink->nb;j++) 	{ p[0]+=modellink->w[j]*points[3*modellink->pts[j]]; p[1]+=modellink->w[j]*points[3*modellink->pts[j]+1]; p[2]+=modellink->w[j]*points[3*modellink->pts[j]+2]; }
}

void CSimplexSurf::GetAxialLinks_NaxisFromPmodel(int index,P_float n[3])
{
n[0]=0; n[1]=0; n[2]=0;
AXIALLINK *modellink=Axis_model->GetAxialLinks()+index;
for(int j=0;j<modellink->nb;j++) 	{ n[0]+=modellink->w[j]*normals[3*modellink->pts[j]]; n[1]+=modellink->w[j]*normals[3*modellink->pts[j]+1]; n[2]+=modellink->w[j]*normals[3*modellink->pts[j]+2]; }
}

P_float CSimplexSurf::GetAxialLinks_RaxisFromBaryCoord(const int nbpt,const int *pts,const P_float* w)
{
P_float R=0; for(int j=0;j<nbpt;j++) 	R+=w[j]*axiallinks[pts[j]].Rref;
return R;
}

P_float* CSimplexSurf::ComputeAxialLinks_error()		// compute links and compare to radius
{
if(Axis_model==NULL || axiallinks==NULL) return NULL;
if(Axis_model->GetAxialLinks()==NULL) return NULL;

if(GetTopo() || Axis_model->GetTopo()) 	{UpdateAxialLinks_links();	UpdateAxialLinks_Rmean();}

int i,j,count=0;	P_float *err=new P_float[2],derr,p[3];
err[0]=0; err[1]=0;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
for(i=0;i<nb_points;i++) axiallinks[i].err=0;
for(i=0;i<Axis_model->GetNumberOfPoints();i++)
	if(Axis_model->GetAttachedPoint(i)==0)
		{
		GetAxialLinks_PaxisFromPmodel(i,p);
		derr=abs(modellinks[i].Rref-dist3D(Axis_model->GetPoint(i),p));
		err[0]+=derr; err[1]+=derr*derr;
		for(j=0;j<modellinks[i].nb;j++) axiallinks[modellinks[i].pts[j]].err+=derr*modellinks[i].w[j];
		count++;
		}
err[0]=err[0]/count;
err[1]=sqrt(err[1]/count-err[0]*err[0]);

return err;
}


void CSimplexSurf::SetAxialLinks_dR(P_float dr)
{
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
for(int i=0;i<nb_points;i++) axiallinks[i].dR=0;
for(int i=0;i<Axis_model->GetNumberOfPoints();i++)	modellinks[i].dR=0;
}



void CSimplexSurf::UpdateAxialLinks_R()		// compute axis point radius
{
if(AxialLinks_UniformR) return;
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;
int i,j; 
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
for(i=0;i<nb_points;i++) {axiallinks[i].dR-=axiallinks[i].Rref; axiallinks[i].Rref=dist3D(GetPoint(i),Axis_model->GetPoint(Axis_model->GetClosestPoint(GetPoint(i)))); axiallinks[i].dR+=axiallinks[i].Rref; }
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	
	{
	modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	
	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref; 
	}
UpdateAxialLinks_FilterR();
}

void CSimplexSurf::UpdateAxialLinks_Rmean()		// compute axis point radius based on mean weighted distances
{
if(AxialLinks_UniformR) return;
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;

int i,j,index;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
for(i=0;i<nb_points;i++)	{axiallinks[i].dR-=axiallinks[i].Rref; axiallinks[i].Rref=0;	axiallinks[i].W=0;}
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	
	if(Axis_model->GetAttachedPoint(i)==0)
		for(j=0;j<modellinks[i].nb;j++)
			{
			index=modellinks[i].pts[j];
			axiallinks[index].Rref+=modellinks[i].dist*modellinks[i].w[j];
			axiallinks[index].W+=modellinks[i].w[j];
			}
for(i=0;i<nb_points;i++) 
	{
	if(axiallinks[i].W!=0) axiallinks[i].Rref=axiallinks[i].Rref/axiallinks[i].W;
	else axiallinks[i].Rref=dist3D(GetPoint(i),Axis_model->GetPoint(Axis_model->GetClosestPoint(GetPoint(i)))); 
	axiallinks[i].dR+=axiallinks[i].Rref;
	}
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	{modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref;}
UpdateAxialLinks_FilterR();
}


void CSimplexSurf::UpdateAxialLinks_Ropt(int nb_it,P_float epsilon)		// iterate to optimze radius
{
if(AxialLinks_UniformR) return;
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;

int i,j,k,index;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();

UpdateAxialLinks_Rmean();

// update Ri=sum(wij.Rj)  and isin
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	
	{
	modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref;  modellinks[i].dR+=modellinks[i].Rref; 
	if(modellinks[i].Rref>modellinks[i].dist) modellinks[i].isin=true; else modellinks[i].isin=false;
	}

P_float derr,dR;
//P_float err=ComputeAxialLinks_error()*Axis_model->GetNumberOfPoints();	//debug
for(k=0;k<nb_it;k++)	
	{
	derr=0;
	j= rand()*nb_points/RAND_MAX; if(j==nb_points) j=nb_points-1;
	for(i=0;i<axiallinks[j].nb;i++)	if(modellinks[axiallinks[j].pts[i]].isin) derr+=axiallinks[j].w[i]; else derr-=axiallinks[j].w[i];
	if(derr<0) dR=epsilon; else dR=-epsilon;
	axiallinks[j].Rref+=dR; axiallinks[j].dR+=dR;
	for(i=0;i<axiallinks[j].nb;i++)	
		{
		index=axiallinks[j].pts[i];
		modellinks[index].Rref+=dR*axiallinks[j].w[i];
		modellinks[index].dR+=dR*axiallinks[j].w[i];
		if(modellinks[index].Rref>modellinks[index].dist) modellinks[index].isin=true; else modellinks[index].isin=false;
		}
//	err+=derr*epsilon;	//debug
	}
//err=err/Axis_model->GetNumberOfPoints();	// estimated error reduction //debug
//P_float err2=ComputeAxialLinks_error();	// real error reduction //debug
UpdateAxialLinks_FilterR();
}


void CSimplexSurf::UpdateAxialLinks_SetUniformR(P_float R)
{
int i,j;
if(R==-1) {AxialLinks_UniformR=false; return;}
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
for(i=0;i<nb_points;i++) {axiallinks[i].dR-=axiallinks[i].Rref; axiallinks[i].Rref=R; axiallinks[i].dR+=axiallinks[i].Rref;}
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	{modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref;}
AxialLinks_UniformR=true;
}

void CSimplexSurf::UpdateAxialLinks_FilterR()		
{
int i,j;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
for(i=0;i<nb_points;i++) if(axiallinks[i].Rref<MIN_RADIUS) {axiallinks[i].dR-=axiallinks[i].Rref; axiallinks[i].Rref=MIN_RADIUS; axiallinks[i].dR+=axiallinks[i].Rref;}
if(Axial_smoothR_axis2!=NULL) UpdateAxialLinks_SmoothR(Axial_smoothR_axis2);
for(i=0;i<Axial_smoothR;i++) UpdateAxialLinks_SmoothR();
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	{modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref;}
}


void CSimplexSurf::SetAxialLinks_smoothR(int nb_it) {Axial_smoothR=nb_it;}
void CSimplexSurf::SetAxialLinks_smoothR(CSimplexSurf* axis2) {Axial_smoothR_axis2=axis2;}

void CSimplexSurf::UpdateAxialLinks_SmoothR()
{
int i,j;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
P_float* r=new P_float[nb_points];
for(i=0;i<nb_points;i++) r[i]=(axiallinks[i].Rref+axiallinks[neighbors[3*i]].Rref +axiallinks[neighbors[3*i+1]].Rref+axiallinks[neighbors[3*i+2]].Rref)/4;
for(i=0;i<nb_points;i++) {axiallinks[i].dR-=axiallinks[i].Rref; axiallinks[i].Rref=r[i]; axiallinks[i].dR+=axiallinks[i].Rref;}
for(i=0;i<Axis_model->GetNumberOfPoints();i++)	{modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref;}
free(r);
}


void CSimplexSurf::UpdateAxialLinks_SmoothR(CSimplexSurf* axis2)
{
AXIALLINK *axiallinks2=axis2->GetAxialLinks();
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
AXIALLINK *modellinks2=axis2->GetAxis_model()->GetAxialLinks();
int i,j; //P_float d;
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1]; 
for(i=0;i<GetNumberOfPoints();i++)
	{
	closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10;
	axiallinks[i].dR-=axiallinks[i].Rref; 
	axis2->GetClosest(closest,GetPoint(i)); 

//	if(closest->dist2<10) axiallinks[i].Rref=(axiallinks[i].Rref+sqrt(closest->dist2)/2.)/2.;
//	if(closest->dist2<10) axiallinks[i].Rref=(axiallinks[i].Rref+sqrt(closest->dist2))/3.;
	if(closest->dist2<20) axiallinks[i].Rref=sqrt(closest->dist2)/2.;
//	j=axis2->GetClosestPoint(GetPoint(i));  d=dist3D(GetPoint(i),axis2->GetPoint(j)); if(d<10) axiallinks[i].Rref=(axiallinks[i].Rref+axiallinks2[j].Rref)/2.;

	axiallinks[i].dR+=axiallinks[i].Rref;
	free(closest->pts); free(closest->weights);
	}
free(closest); 

for(i=0;i<Axis_model->GetNumberOfPoints();i++)	{modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref;}
}

void CSimplexSurf::SetAxialLinks_blockside(int side)	{Axial_blockside=side;}

void CSimplexSurf::UpdateAxialLinks_RegularizeSide()	
{
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
P_float dp,p[3],u[3];
if(Axial_blockside!=-1)
	{
	if(Axial_blockside==0) for(int i=0;i<Axis_model->GetNumberOfPoints();i++) modellinks[i].side=false;
	else for(int i=0;i<Axis_model->GetNumberOfPoints();i++) modellinks[i].side=true;
	}
else
	for(int i=0;i<Axis_model->GetNumberOfPoints();i++)
		if(modellinks[i].nb>2)
			{
			GetAxialLinks_PaxisFromPmodel(i,p);
			Axis_model->GetPoint(i,u); u[0]-=p[0]; u[1]-=p[1]; u[2]-=p[2];
			dp=dotproduct(u,Axis_model->GetNormal(i));
			if(dp<0)	modellinks[i].side=(modellinks[i].side)?false:true;
			}
}

void CSimplexSurf::UpdateAxialLinks_links()		// compute links
{
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;

int i,j,index;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1]; 

// update model links
for(i=0;i<Axis_model->GetNumberOfPoints();i++)
	{
	closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10;
	GetClosest(closest,Axis_model->GetPoint(i));
	modellinks[i].nb=closest->nb; 
	if(modellinks[i].pts!=NULL) free(modellinks[i].pts);	if(modellinks[i].w!=NULL) free(modellinks[i].w);
	modellinks[i].pts=new int[modellinks[i].nb];	modellinks[i].w=new P_float[modellinks[i].nb];
	memcpy(modellinks[i].pts,closest->pts,modellinks[i].nb*sizeof(int));	memcpy(modellinks[i].w,closest->weights,modellinks[i].nb*sizeof(P_float));
	modellinks[i].dist=sqrt(closest->dist2);
	modellinks[i].side=(closest->side==1)?true:false;
	free(closest->pts); free(closest->weights);
	
	modellinks[i].W=0;
	for(j=0;j<modellinks[i].nb;j++) modellinks[i].W+=modellinks[i].w[j]*modellinks[i].w[j];
	modellinks[i].W=1/modellinks[i].W;
	}
free(closest);
UpdateAxialLinks_RegularizeSide();

// update axis links
for(i=0;i<nb_points;i++)
	{
	if(axiallinks[i].pts!=NULL) free(axiallinks[i].pts); axiallinks[i].pts=new int[1];
	if(axiallinks[i].w!=NULL) free(axiallinks[i].w); axiallinks[i].w=new P_float[1];
	axiallinks[i].nb=0;
	}
for(i=0;i<Axis_model->GetNumberOfPoints();i++)
	{
	for(j=0;j<modellinks[i].nb;j++) 
		{
		index=modellinks[i].pts[j];
		axiallinks[index].nb++;
		axiallinks[index].w=(P_float*)realloc(axiallinks[index].w,axiallinks[index].nb*sizeof(P_float));
		axiallinks[index].pts=(int*)realloc(axiallinks[index].pts,axiallinks[index].nb*sizeof(int));
		axiallinks[index].w[axiallinks[index].nb-1]=modellinks[i].W*modellinks[i].w[j];
		axiallinks[index].pts[axiallinks[index].nb-1]=i;
		}
	}
}



void CSimplexSurf::UpdateAxialLinks_forces()
{
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;

int i,j;
P_float l,p[3],dp,n[3],u[3];
AXIALLINK *modellinks=Axis_model->GetAxialLinks();
bool flipping=true;

if(AxialConstraint_mode==0) // repulsive force to axis
	{
	// Rref matching using closest
	UpdateAxialLinks_Rmean();
	if(Axis_model->GetInternalForce_AxialConstraint())
		for(i=0;i<Axis_model->GetNumberOfPoints();i++)
			{
			Axis_model->GetPoint(i,p);		j=GetClosestPoint(p);
			p[0]-=points[3*j]; p[1]-=points[3*j+1]; p[2]-=points[3*j+2];
			l=norm(p); if(l==0) l=1E-10; l=(axiallinks[j].Rref-l)/l;
			modellinks[i].f[0]=l*p[0]; modellinks[i].f[1]=l*p[1]; modellinks[i].f[2]=l*p[2];
			}
	
	if(GetInternalForce_AxialConstraint())
		for(i=0;i<nb_points;i++)
			{
			Axis_model->GetPoint(Axis_model->GetClosestPoint(points+3*i),p);
			p[0]-=points[3*i]; p[1]-=points[3*i+1]; p[2]-=points[3*i+2];
			l=norm(p); if(l==0) l=1E-10; l=(axiallinks[i].Rref-l)/l;
			axiallinks[i].f[0]=-l*p[0]; axiallinks[i].f[1]=-l*p[1];	axiallinks[i].f[2]=-l*p[2];
			}
/*
	P_float friction=0.5; // 1: complete attenuation of tangential component

	P_float epsilon=1;
	for(i=0;i<nb_points;i++)
		{
		Axis_model->GetPoint(Axis_model->GetClosestPoint(points+3*i),p);
		p[0]-=points[3*i]; p[1]-=points[3*i+1]; p[2]-=points[3*i+2];
		axiallinks[i].dR-=axiallinks[i].Rref; axiallinks[i].Rref=norm(p);	if(axiallinks[i].Rref==0) axiallinks[i].Rref=1E-10; axiallinks[i].dR+=axiallinks[i].Rref;
		axiallinks[i].f[0]=-epsilon*p[0]/axiallinks[i].Rref; axiallinks[i].f[1]=-epsilon*p[1]/axiallinks[i].Rref;	axiallinks[i].f[2]=-epsilon*p[2]/axiallinks[i].Rref;
		}
	// attenuate tangential component
	decompose(i,axiallinks[i].f,ftg,fn);
	axiallinks[i].f[0]-=friction*ftg[0];		axiallinks[i].f[1]-=friction*ftg[1];		axiallinks[i].f[2]-=friction*ftg[2];
*/
/*
	P_float tolerence=5,epsilon=3;
	P_float dists2[3],F1[4],F2[4],F2p[4],F3[3],R2,R3; int closests[3];
	for(i=0;i<nb_points;i++)
		{
		Axis_model->GetClosestPoints(closests,dists2,points+3*i);
		F1[0]=GetPoint(i)[0]-Axis_model->GetPoint(closests[0])[0];		F1[1]=GetPoint(i)[1]-Axis_model->GetPoint(closests[0])[1];		F1[2]=GetPoint(i)[2]-Axis_model->GetPoint(closests[0])[2];	F1[3]=norm(F1);
		R2=GetEquidistant2(F2,F2p,GetPoint(i),Axis_model->GetPoint(closests[0]),Axis_model->GetPoint(closests[1]),dists2[0],dists2[1]);
		R3=GetEquidistant3(F3,GetPoint(i),Axis_model->GetPoint(closests[0]),Axis_model->GetPoint(closests[1]),Axis_model->GetPoint(closests[2]));
		if(F3[3]<tolerence && R3>5)	
			{
			memcpy(axiallinks[i].f,F3,3*sizeof(P_float));
			axiallinks[i].Rref=R3;
			}
		else if(F2p[3]<tolerence)	
			{
			axiallinks[i].f[0]=F2p[0]*epsilon/F2p[3]; axiallinks[i].f[1]=F2p[1]*epsilon/F2p[3]; axiallinks[i].f[2]=F2p[2]*epsilon/F2p[3];
			axiallinks[i].Rref=R2;			
			}
		else 
			{
			axiallinks[i].f[0]=F1[0]*epsilon/F1[3]; axiallinks[i].f[1]=F1[1]*epsilon/F1[3]; axiallinks[i].f[2]=F1[2]*epsilon/F1[3];
			axiallinks[i].Rref=F1[3];			
			}
		if(axiallinks[i].Rref==0) axiallinks[i].Rref=1E-10;	
		}
*/

	if(Axis_model->GetInternalForce_AxialConstraint())
		for(i=0;i<Axis_model->GetNumberOfPoints();i++)
			{
			Axis_model->GetPoint(i,p);		j=GetClosestPoint(p);
			p[0]-=points[3*j]; p[1]-=points[3*j+1]; p[2]-=points[3*j+2];
			l=norm(p); if(l==0) l=1E-10; l=(axiallinks[j].Rref-l)/l; 
			modellinks[i].f[0]=l*p[0]; modellinks[i].f[1]=l*p[1]; modellinks[i].f[2]=l*p[2];
			}
	}
else if(AxialConstraint_mode==1)	
	{
	for(i=0;i<Axis_model->GetNumberOfPoints();i++)
		{
		GetAxialLinks_PaxisFromPmodel(i,p); u[0]=p[0]-Axis_model->GetPoint(i)[0]; u[1]=p[1]-Axis_model->GetPoint(i)[1]; u[2]=p[2]-Axis_model->GetPoint(i)[2];
		modellinks[i].dist=norm(u); if(modellinks[i].dist!=0) {u[0]=u[0]/modellinks[i].dist; u[1]=u[1]/modellinks[i].dist; u[2]=u[2]/modellinks[i].dist;}
		p[0]-=u[0]*modellinks[i].Rref; p[1]-=u[1]*modellinks[i].Rref; p[2]-=u[2]*modellinks[i].Rref;

		if(flipping)
			{
			GetAxialLinks_NaxisFromPmodel(i,n); 	dp=dotproduct(n,u); 
			if((dp>0 && modellinks[i].side) || (dp<0 && !modellinks[i].side)) // flip
				{
				dp=dp*2*modellinks[i].Rref;
				p[0]+=dp*n[0]; p[1]+=dp*n[1]; p[2]+=dp*n[2];
				}
			}

		modellinks[i].f[0]=p[0]-Axis_model->GetPoint(i)[0];	modellinks[i].f[1]=p[1]-Axis_model->GetPoint(i)[1]; modellinks[i].f[2]=p[2]-Axis_model->GetPoint(i)[2]; 
		}
	if(GetInternalForce_AxialConstraint())
		for(i=0;i<nb_points;i++)
			{
			axiallinks[i].f[0]=0; axiallinks[i].f[1]=0; axiallinks[i].f[2]=0; 
			for(j=0;j<axiallinks[i].nb;j++) { axiallinks[i].f[0]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[0]; axiallinks[i].f[1]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[1]; axiallinks[i].f[2]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[2]; }
			if(axiallinks[i].nb!=0) {axiallinks[i].f[0]=-axiallinks[i].f[0]/axiallinks[i].nb; axiallinks[i].f[1]=-axiallinks[i].f[1]/axiallinks[i].nb; axiallinks[i].f[2]=-axiallinks[i].f[2]/axiallinks[i].nb;}
			}
	}
else if(AxialConstraint_mode==2) // Rref matching using weighted sum of forces
	{
	UpdateAxialLinks_Rmean();

	for(i=0;i<Axis_model->GetNumberOfPoints();i++)
		{
		GetAxialLinks_PaxisFromPmodel(i,p); u[0]=p[0]-Axis_model->GetPoint(i)[0]; u[1]=p[1]-Axis_model->GetPoint(i)[1]; u[2]=p[2]-Axis_model->GetPoint(i)[2];
		modellinks[i].dist=norm(u); if(modellinks[i].dist!=0) {u[0]=u[0]/modellinks[i].dist; u[1]=u[1]/modellinks[i].dist; u[2]=u[2]/modellinks[i].dist;}
		p[0]-=u[0]*modellinks[i].Rref; p[1]-=u[1]*modellinks[i].Rref; p[2]-=u[2]*modellinks[i].Rref;

		if(flipping)
			{
			GetAxialLinks_NaxisFromPmodel(i,n); 	dp=dotproduct(n,u); 
			if((dp>0 && modellinks[i].side) || (dp<0 && !modellinks[i].side)) // flip
				{
				dp=dp*2*modellinks[i].Rref;
				p[0]+=dp*n[0]; p[1]+=dp*n[1]; p[2]+=dp*n[2];
				}
			}

		modellinks[i].f[0]=p[0]-Axis_model->GetPoint(i)[0];	modellinks[i].f[1]=p[1]-Axis_model->GetPoint(i)[1]; modellinks[i].f[2]=p[2]-Axis_model->GetPoint(i)[2]; 
		}
	if(GetInternalForce_AxialConstraint())
		for(i=0;i<nb_points;i++)
			{
			axiallinks[i].f[0]=0; axiallinks[i].f[1]=0; axiallinks[i].f[2]=0; 
			for(j=0;j<axiallinks[i].nb;j++) { axiallinks[i].f[0]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[0]; axiallinks[i].f[1]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[1]; axiallinks[i].f[2]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[2]; }
			if(axiallinks[i].nb!=0) {axiallinks[i].f[0]=-axiallinks[i].f[0]/axiallinks[i].nb; axiallinks[i].f[1]=-axiallinks[i].f[1]/axiallinks[i].nb; axiallinks[i].f[2]=-axiallinks[i].f[2]/axiallinks[i].nb;}
			}
	}

else if(AxialConstraint_mode==3) // Update (R and links) and Rref matching using weighted sum of forces 
	{
	UpdateAxialLinks_links(); 
//UpdateAxialLinks_Ropt(100,.2);
	UpdateAxialLinks_Rmean(); 

/*
// Update and Rref matching using weighted sum of forces (radius and barycentric coord method)

	for(i=0;i<Axis_model->GetNumberOfPoints();i++)
		{
		GetAxialLinks_NaxisFromPmodel(i,u);
		u[0]=u[0]*modellinks[i].Rref; u[1]=u[1]*modellinks[i].Rref; u[2]=u[2]*modellinks[i].Rref;
		if(flipping) if(!modellinks[i].side) {u[0]=-u[0]; u[1]=-u[1]; u[2]=-u[2];}
		GetAxialLinks_PaxisFromPmodel(i,p); modellinks[i].f[0]=u[0]+p[0]-Axis_model->GetPoint(i)[0]; modellinks[i].f[1]=u[1]+p[1]-Axis_model->GetPoint(i)[1]; modellinks[i].f[2]=u[2]+p[2]-Axis_model->GetPoint(i)[2];
		}
*/

	for(i=0;i<Axis_model->GetNumberOfPoints();i++)
		{
		GetAxialLinks_PaxisFromPmodel(i,p); u[0]=p[0]-Axis_model->GetPoint(i)[0]; u[1]=p[1]-Axis_model->GetPoint(i)[1]; u[2]=p[2]-Axis_model->GetPoint(i)[2];
		modellinks[i].dist=norm(u); if(modellinks[i].dist!=0) {u[0]=u[0]/modellinks[i].dist; u[1]=u[1]/modellinks[i].dist; u[2]=u[2]/modellinks[i].dist;}
		p[0]-=u[0]*modellinks[i].Rref; p[1]-=u[1]*modellinks[i].Rref; p[2]-=u[2]*modellinks[i].Rref;

		if(flipping)
			{
			GetAxialLinks_NaxisFromPmodel(i,n); 	dp=dotproduct(n,u); 
			if((dp>0 && modellinks[i].side) || (dp<0 && !modellinks[i].side)) // flip
				{
				dp=dp*2*modellinks[i].Rref;
				p[0]+=dp*n[0]; p[1]+=dp*n[1]; p[2]+=dp*n[2];
				}
			}

		modellinks[i].f[0]=p[0]-Axis_model->GetPoint(i)[0];	modellinks[i].f[1]=p[1]-Axis_model->GetPoint(i)[1]; modellinks[i].f[2]=p[2]-Axis_model->GetPoint(i)[2]; 
		}
	if(GetInternalForce_AxialConstraint())
		for(i=0;i<nb_points;i++)
			{
			axiallinks[i].f[0]=0; axiallinks[i].f[1]=0; axiallinks[i].f[2]=0; 
			for(j=0;j<axiallinks[i].nb;j++) { axiallinks[i].f[0]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[0]; axiallinks[i].f[1]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[1]; axiallinks[i].f[2]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[2]; }
			if(axiallinks[i].nb!=0) {axiallinks[i].f[0]=-axiallinks[i].f[0]/axiallinks[i].nb; axiallinks[i].f[1]=-axiallinks[i].f[1]/axiallinks[i].nb; axiallinks[i].f[2]=-axiallinks[i].f[2]/axiallinks[i].nb;}
			}
	}

else if(AxialConstraint_mode==4) // Update (R) and Rref matching using weighted sum of forces 
	{
	UpdateAxialLinks_links(); 
	for(i=0;i<Axis_model->GetNumberOfPoints();i++)	{modellinks[i].dR-=modellinks[i].Rref; modellinks[i].Rref=0; 	for(j=0;j<modellinks[i].nb;j++) modellinks[i].Rref+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].Rref; modellinks[i].dR+=modellinks[i].Rref;}

	for(i=0;i<Axis_model->GetNumberOfPoints();i++)
		{
		GetAxialLinks_PaxisFromPmodel(i,p); u[0]=p[0]-Axis_model->GetPoint(i)[0]; u[1]=p[1]-Axis_model->GetPoint(i)[1]; u[2]=p[2]-Axis_model->GetPoint(i)[2];
		modellinks[i].dist=norm(u); if(modellinks[i].dist!=0) {u[0]=u[0]/modellinks[i].dist; u[1]=u[1]/modellinks[i].dist; u[2]=u[2]/modellinks[i].dist;}
		p[0]-=u[0]*modellinks[i].Rref; p[1]-=u[1]*modellinks[i].Rref; p[2]-=u[2]*modellinks[i].Rref;

		if(flipping)
			{
			GetAxialLinks_NaxisFromPmodel(i,n); 	dp=dotproduct(n,u); 
			if((dp>0 && modellinks[i].side) || (dp<0 && !modellinks[i].side)) // flip
				{
				dp=dp*2*modellinks[i].Rref;
				p[0]+=dp*n[0]; p[1]+=dp*n[1]; p[2]+=dp*n[2];
				}
			}

		modellinks[i].f[0]=p[0]-Axis_model->GetPoint(i)[0];	modellinks[i].f[1]=p[1]-Axis_model->GetPoint(i)[1]; modellinks[i].f[2]=p[2]-Axis_model->GetPoint(i)[2]; 
		}
	if(GetInternalForce_AxialConstraint())
		for(i=0;i<nb_points;i++)
			{
			axiallinks[i].f[0]=0; axiallinks[i].f[1]=0; axiallinks[i].f[2]=0; 
			for(j=0;j<axiallinks[i].nb;j++) { axiallinks[i].f[0]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[0]; axiallinks[i].f[1]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[1]; axiallinks[i].f[2]+=axiallinks[i].w[j]*modellinks[axiallinks[i].pts[j]].f[2]; }
			if(axiallinks[i].nb!=0) {axiallinks[i].f[0]=-axiallinks[i].f[0]/axiallinks[i].nb; axiallinks[i].f[1]=-axiallinks[i].f[1]/axiallinks[i].nb; axiallinks[i].f[2]=-axiallinks[i].f[2]/axiallinks[i].nb;}
			}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::GetNeighbors_c(int index,int nc[3]) {memcpy(nc,neighbors_c+3*index,3*sizeof(int));}
int* CSimplexSurf::GetNeighbors_c(int index) {return neighbors_c+3*index;}


void CSimplexSurf::SetNeighbors_c(int index, int i, int j, int k) {
	neighbors_c[3 * index + 0] = i;
	neighbors_c[3 * index + 1] = j;
	neighbors_c[3 * index + 2] = k;
}

void CSimplexSurf::Flip(bool reversepts) {
	int temp;
	P_float ftemp;
	for (int i = 0; i < nb_points; i++) {
		int *n = GetNeighbors(i);
		int new_n[3] = {n[0], n[1], n[2]};

		//temp = neighbors[3 * i + 1];		
		//neighbors[3 * i + 1] = neighbors[3 * i];		
		//neighbors[3 * i] = temp;
		
		temp = new_n[1];
		new_n[1] = new_n[0];
		new_n[0] = temp;
		SetNeighbors(i, new_n[0], new_n[1], new_n[2]);

	
		ftemp = params[3 * i + 1];		
		params[3 * i + 1] = params[3 * i];		
		params[3 * i] = ftemp;
		
		if(mass_inv[3 * i] != 0 && reversepts) 
			points[3 * i] = -points[3 * i];
	}

	UpdateFlipped();
	Equilibrium();
}

void CSimplexSurf::UpdateFlipped()
{
int n[3]; GetNeighbors(cells[0][2],n);
if(n[0]==cells[0][1] && n[2]==cells[0][3]) Flipped=true;
else if(n[1]==cells[0][1] && n[0]==cells[0][3]) Flipped=true;
else if(n[2]==cells[0][1] && n[1]==cells[0][3]) Flipped=true;
}

/**
 * Determines whether the cell is a cell of the shared boundary. 
 * Cells on the shared boundary will be made up of points whose
 * two material indices are non-zero. 
 */
bool CSimplexSurf::IsSharedBoundaryCell(vtkSmartPointer<vtkCell> cell) {
	bool sb = true;
	for (int i = 0; i < cell->GetNumberOfPoints(); i++) {
		int *mat = GetPointMaterialIndices(cell->GetPointId(i));
		
		// If any of the material indices of any point is zero
		// Then it means that the cell is not on the shared boundary.
		if (mat[0] == 0 || mat[1] == 0) { 
			sb = false;
			break;
		}
	}
	return sb;
}

/** 
 * This function returns the number of connected points for a given point. 
 * 
 */
vtkSmartPointer<vtkIdList> CSimplexSurf::GetNeighboringPoints(int index, vtkSmartPointer<vtkPolyData> mesh) {
	vtkSmartPointer<vtkIdList> neighboringPoints = vtkSmartPointer<vtkIdList>::New();

	//vtkSmartPointer<vtkPolyData> mesh = getAsVTKPolyData();
	//mesh->BuildCells();
	//mesh->BuildLinks();
	//mesh->Update();
	
	vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(index, cellList);

	for (int i = 0; i < cellList->GetNumberOfIds(); i++) {
		vtkSmartPointer<vtkCell> cell = mesh->GetCell(cellList->GetId(i));

		for (int j = 0; j < cell->GetNumberOfPoints(); j++) {
			int pt = cell->GetPointId(j);
			if (pt != index) { // Make sure the point is not the same as 'index'
				//mesh->BuildLinks();
				if (mesh->IsEdge(index, pt) != 0) { // If the point makes and edge with 'index'
					
					// Check in neighboringPoints to make sure that the point
					// is not already there. 
					bool found = false;
					for (int k = 0; k < neighboringPoints->GetNumberOfIds(); k++) {
						if (pt == neighboringPoints->GetId(k)) {
							found = true;
							//break;
						}
					}

					if (found == false) {
						neighboringPoints->InsertNextId(pt);
					}
				}
			}
		}
	}

	//cout << "Point " << index << " has neighboring points: ";
	//for (int i = 0; i < neighboringPoints->GetNumberOfIds(); i++) {
	//	cout << neighboringPoints->GetId(i) << ", ";
	//}
	//cout << endl;

	return neighboringPoints;
}

/**
 * Updates the neighbors of all points. This version is capable of producing oriented neighbors. 
 * The order of neighbors matterr a lot. 
 */
void CSimplexSurf::UpdateNeighbors_vtkPolyData() {	
	//ofstream f;
	//f.open("E:\\workspace\\CSimplexMesh\\CSimplexMesh\\data\\MultiMaterialNeighbors.txt");

	vtkSmartPointer<vtkPolyData> pdata = getAsVTKPolyData();
	pdata->BuildCells();
	pdata->BuildLinks();	

	// For every qth point
	for (int q = 0; q < pdata->GetNumberOfPoints(); q++) {
		vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
		pdata->GetPointCells(q, cellList);
		
		vtkSmartPointer<vtkIdList> qneighbors = GetNeighboringPoints(q, pdata);
		int no_neighboring_pts = qneighbors->GetNumberOfIds();

		// For a regular single material simplex mesh, each point will be shared by exactly 3 cells. 
		// In the multi-material case, if the point is shared by 3 or 4 or 5 cells, it means that the point
		// is on the non-manifold edge of the shared boundary, but this particular point can still
		// have only 3 neighboring points. 
		// We detect those neighbors in this part. 
		if (no_neighboring_pts == 3) {
			int qneigh[3] = {-1, -1, -1}; // To hold the neighbors of the qth point

			// Using one cell, obtain two neighbors for the qth point
			// In this way, the ordering of the neighbors is the same as the ordering of the cells. 
			vtkSmartPointer<vtkCell> cell0 = pdata->GetCell(cellList->GetId(0));
			for (int v = 0; v < cell0->GetNumberOfPoints(); v++) {
				int before_q = -1, after_q = -1;
				int qth_point = cell0->GetPointId(v);
				if (qth_point == q) {
					if (v == 0) {
						before_q = cell0->GetPointId(cell0->GetNumberOfPoints() - 1);
						after_q = cell0->GetPointId(1);
					}
					else if (v == cell0->GetNumberOfPoints() - 1) {
						before_q = cell0->GetPointId(cell0->GetNumberOfPoints() - 2);
						after_q = cell0->GetPointId(0);
					}
					else {
						before_q = cell0->GetPointId(v - 1);
						after_q = cell0->GetPointId(v + 1);
					}
				}

				if (before_q != -1 && after_q != -1) {
					qneigh[0] = before_q;
					qneigh[1] = after_q;

					break;
				}
			}

			// Using another cell, get the third neighbor. 
			// It is possible that the next cell would contain the same 2 neighbors 
			// as the first one. So use a while loop to go through the cell list until
			// the third neighbor is found. 
			bool found = false;
			int icell = 1; // Cell 0 in cellList was used to locate the first 2 neighbors. So start from 1.
			while (found == false && icell < cellList->GetNumberOfIds()) {
				vtkSmartPointer<vtkCell> cell1 = pdata->GetCell(cellList->GetId(icell));
				for (int b = 0; b < cell1->GetNumberOfPoints(); b++) {
					int before_q = -1, after_q = -1;
					if (cell1->GetPointId(b) == q) {
						if (b == 0) {
							before_q = cell1->GetPointId(cell1->GetNumberOfPoints() - 1);
							after_q = cell1->GetPointId(1);
						}
						else if (b == cell1->GetNumberOfPoints() - 1) {
							before_q = cell1->GetPointId(cell1->GetNumberOfPoints() - 2);
							after_q = cell1->GetPointId(0);
						}
						else {
							before_q = cell1->GetPointId(b - 1);
							after_q = cell1->GetPointId(b + 1);
						}
					}

					if (before_q != -1 && after_q != -1) {
						if ((before_q != qneigh[0]) && (before_q != qneigh[1])) {
							qneigh[2] = before_q;
							found = true;
							break;
						}
						else if ((after_q != qneigh[0]) && (after_q != qneigh[1])) {
							qneigh[2] = after_q;
							found = true;
							break;
						}

						// At this point, the two neighbors are the same as found before. So use the next cell in the list 
						else { 
							icell = icell + 1;
						}
					}
				}
			}
			
			//cout << "Point " << q << " has neighbors " << qneigh[0] << ", " << qneigh[1] << ", " << qneigh[2] << endl;

			if ((qneigh[0] == -1) || (qneigh[1] == -1) || (qneigh[2] == -1)) {
				cout << "\n\nERROR: UpdateNeighbors(): Neighbors not correctly found for point " << q << ".\n\n" << endl;
			}
			SetNeighbors(q, qneigh[0], qneigh[1], qneigh[2]);
		}

		// For multimaterial points having more than 3 neighbors,
		else {
			SetIsMultiMaterialPoint(q, true); // Set the qth point as a multimaterial point

			// Get the material indices of the qth point. 
			int *qMat = GetPointMaterialIndices(q);
			
			// Separate the neighbors wrt to qth point's material indices
			// One set of neighbors for qMat[0]
			int qMat0_neighbors[3] = {-1, -1, -1};
			int c0 = 0;
			for (int i = 0; i < no_neighboring_pts; i++) {
				int *n_mat = GetPointMaterialIndices(qneighbors->GetId(i));

				if (qMat[0] == n_mat[0] || qMat[0] == n_mat[1]) {
					qMat0_neighbors[c0] = qneighbors->GetId(i);
					c0 = c0 + 1;
				}
			}
			
			// Another set of neighbors for qMat[1]
			int qMat1_neighbors[3] = {-1, -1, -1};
			int c1 = 0;
			for (int i = 0; i < no_neighboring_pts; i++) {
				int *n_mat = GetPointMaterialIndices(qneighbors->GetId(i));

				if (qMat[1] == n_mat[0] || qMat[1] == n_mat[1]) {
					qMat1_neighbors[c1] = qneighbors->GetId(i);
					c1 = c1 + 1;
				}
			}
			

			// Reorder the neighbors for qMat[0]
			int i1 = 0;
			bool found1 = false;
			while (i1 < 3 && found1 == false) {
				// First find the neighbor that is not a point on the shared boundary. 
				int *nMat = GetPointMaterialIndices(qMat0_neighbors[i1]);
				if (nMat[0] == 0 || nMat[1] == 0) {
					int j = 0;
					while (j < cellList->GetNumberOfIds() && found1 == false) {
						// Then locate a cell that is used by both q and the neighbor that is not a point on the shared boundary
						// Using this cell, locate q in the cell, and the before and after points for q
						if ((pdata->IsPointUsedByCell(q, cellList->GetId(j)) == 1) && (pdata->IsPointUsedByCell(qMat0_neighbors[i1], cellList->GetId(j)) == 1)) {
							vtkSmartPointer<vtkCell> jCell = pdata->GetCell(cellList->GetId(j));
							for (int k = 0; k < jCell->GetNumberOfPoints(); k++) {
								int before_q = -1, after_q = -1;
								int qth_point = jCell->GetPointId(k);
								if (qth_point == q) {
									if (k == 0) {
										before_q = jCell->GetPointId(jCell->GetNumberOfPoints() - 1);
										after_q = jCell->GetPointId(1);
									}
									else if (k == jCell->GetNumberOfPoints() - 1) {
										before_q = jCell->GetPointId(jCell->GetNumberOfPoints() - 2);
										after_q = jCell->GetPointId(0);
									}
									else {
										before_q = jCell->GetPointId(k - 1);
										after_q = jCell->GetPointId(k + 1);
									}
								}

								// Then reorder the qMat0_neighbors array using before_q and after_q
								if (before_q != -1 && after_q != -1) {
									//cout << "Point " << q << endl;
									//cout << "\tNeighbors for material " << qMat[0] << " before: " << qMat0_neighbors[0] << ", " << qMat0_neighbors[1] << ", " << qMat0_neighbors[2] << endl;


									if (qMat0_neighbors[0] != before_q && qMat0_neighbors[0] != after_q) {
										qMat0_neighbors[2] = qMat0_neighbors[0];
										qMat0_neighbors[0] = before_q;
										qMat0_neighbors[1] = after_q;
									}
									else if (qMat0_neighbors[1] != before_q && qMat0_neighbors[1] != after_q) {
										qMat0_neighbors[2] = qMat0_neighbors[1];
										qMat0_neighbors[0] = before_q;
										qMat0_neighbors[1] = after_q;
									}
									else {
										qMat0_neighbors[0] = before_q;
										qMat0_neighbors[1] = after_q;
									}
									
									//cout << "\tNeighbors for material " << qMat[0] << " after: " << qMat0_neighbors[0] << ", " << qMat0_neighbors[1] << ", " << qMat0_neighbors[2] << endl;

									found1 = true;
									break;
								}
							}
						}

						j = j + 1;
					}
				}

				i1 = i1 + 1;
			}

			// Reorder the neighbors for qMat[1]
			int i2 = 0;
			bool found2 = false;
			while (i2 < 3 && found2 == false) {
				// First find the neighbor that is not a point on the shared boundary. 
				int *nMat = GetPointMaterialIndices(qMat1_neighbors[i2]);
				if (nMat[0] == 0 || nMat[1] == 0) {
					int j = 0;
					while (j < cellList->GetNumberOfIds() && found2 == false) {
						// Then locate a cell that is used by both q and the neighbor that is not a point on the shared boundary
						// Using this cell, locate q in the cell, and the before and after points for q
						if ((pdata->IsPointUsedByCell(q, cellList->GetId(j)) == 1) && (pdata->IsPointUsedByCell(qMat1_neighbors[i2], cellList->GetId(j)) == 1)) {
							vtkSmartPointer<vtkCell> jCell = pdata->GetCell(cellList->GetId(j));
							for (int k = 0; k < jCell->GetNumberOfPoints(); k++) {
								int before_q = -1, after_q = -1;
								int qth_point = jCell->GetPointId(k);
								if (qth_point == q) {
									if (k == 0) {
										before_q = jCell->GetPointId(jCell->GetNumberOfPoints() - 1);
										after_q = jCell->GetPointId(1);
									}
									else if (k == jCell->GetNumberOfPoints() - 1) {
										before_q = jCell->GetPointId(jCell->GetNumberOfPoints() - 2);
										after_q = jCell->GetPointId(0);
									}
									else {
										before_q = jCell->GetPointId(k - 1);
										after_q = jCell->GetPointId(k + 1);
									}
								}

								// Then reorder the qMat1_neighbors array using before_q and after_q
								if (before_q != -1 && after_q != -1) {
									//cout << "Point " << q << endl;
									//cout << "\tNeighbors for material " << qMat[1] << " before: " << qMat1_neighbors[0] << ", " << qMat1_neighbors[1] << ", " << qMat1_neighbors[2] << endl;


									if (qMat1_neighbors[0] != before_q && qMat1_neighbors[0] != after_q) {
										qMat1_neighbors[2] = qMat1_neighbors[0];
										qMat1_neighbors[0] = before_q;
										qMat1_neighbors[1] = after_q;
									}
									else if (qMat1_neighbors[1] != before_q && qMat1_neighbors[1] != after_q) {
										qMat1_neighbors[2] = qMat1_neighbors[1];
										qMat1_neighbors[0] = before_q;
										qMat1_neighbors[1] = after_q;
									}
									else {
										qMat1_neighbors[0] = before_q;
										qMat1_neighbors[1] = after_q;
									}
									
									//cout << "\tNeighbors for material " << qMat[1] << " after: " << qMat1_neighbors[0] << ", " << qMat1_neighbors[1] << ", " << qMat1_neighbors[2] << endl;

									found2 = true;
									break;
								}
							}
						}

						j = j + 1;
					}
				}

				i2 = i2 + 1;
			}


			if (qMat0_neighbors[0] == -1 || qMat0_neighbors[1] == -1  || qMat0_neighbors[2] == -1) {
				cout << "\n\nERROR: UpdateNeighbors(): Multimaterial neighbors fo qMat0 not correctly found for point " << q << ".\n\n" << endl;
			}
			if (qMat1_neighbors[0] == -1 || qMat1_neighbors[1] == -1  || qMat1_neighbors[2] == -1) {
				cout << "\n\nERROR: UpdateNeighbors(): Multimaterial neighbors fo qMat1 not correctly found for point " << q << ".\n\n" << endl;
			}
			
			//cout << "Point " << q << endl;
			//cout << "\tNeighbors for material " << qMat[0] << ": " << qMat0_neighbors[0] << ", " << qMat0_neighbors[1] << ", " << qMat0_neighbors[2] << endl;
			//cout << "\tNeighbors for material " << qMat[1] << ": " << qMat1_neighbors[0] << ", " << qMat1_neighbors[1] << ", " << qMat1_neighbors[2] << endl;
			
			
		
			SetMultiMaterialNeighbors(qMat[0], q, qMat0_neighbors[0], qMat0_neighbors[1], qMat0_neighbors[2]);
			SetMultiMaterialNeighbors(qMat[1], q, qMat1_neighbors[0], qMat1_neighbors[1], qMat1_neighbors[2]);
		}
	}
}

void CSimplexSurf::UpdateNeighbors_c_vtkPolyData() {
	vtkSmartPointer<vtkPolyData> pdata = getAsVTKPolyData();

	for (int i = 0; i < pdata->GetNumberOfPoints(); i++) {
		vtkSmartPointer<vtkIdList> cellList = vtkSmartPointer<vtkIdList>::New();
		pdata->GetPointCells(i, cellList);

		if (cellList->GetNumberOfIds() == 3) {
			SetNeighbors_c(i, cellList->GetId(0), cellList->GetId(1), cellList->GetId(2));
		}

		else if (cellList->GetNumberOfIds() > 3) {
			SetNeighbors_c(i, cellList->GetId(0), cellList->GetId(1), cellList->GetId(2));
		}
	}

}

void CSimplexSurf::UpdateNeighbors() {	
	int i;
	//for (i = 0; i < nb_cells; i++) 
	//	UpdateNeighbors(i);
	
	UpdateNeighbors_vtkPolyData();
	
	//for (i = 0; i < nb_points; i++) 
	//	UpdateNeighbors_c(i);

	UpdateNeighbors_c_vtkPolyData();

	for (i = 0; i < nb_points; i++) 
		UpdateNeighbors2(i);

	nb_points_border = 0; 
	for (i = 0; i < nb_points; i++) 
		if(neighbors_c[3 * i + 2] == -1) 
			nb_points_border++;

	if (nb_points_border != 0) 
		OrderNeighbors();
}

void CSimplexSurf::UpdateNeighbors_c(int pt_index) {
	int c[3]; 
	SearchCell(pt_index, c);
	//memcpy(neighbors_c + 3 * pt_index, c, 3 * sizeof(int));
	
	SetNeighbors_c(pt_index, c[0], c[1], c[2]);
	//cout << neighbors_c[3 * pt_index] << ", " << neighbors_c[3 * pt_index + 1] << ", " << neighbors_c[3 * pt_index + 2] << endl;
}

void CSimplexSurf::UpdateNeighbors(int cell_index) {
	if (cell_index==-1) 
		return;
	
	int j;
	int id1, id2, id3;
	int nb_p = cells[cell_index][0];
	
	for (j = 1; j < nb_p + 1; j++) {
		id1 = cells[cell_index][j];
		id2 = (j == 1) ? cells[cell_index][nb_p] : cells[cell_index][j - 1];
		id3 = (j == nb_p) ? cells[cell_index][1] : cells[cell_index][j + 1];

		//if(id1<0 || id2<0 || id3<0)  /// debug
		//this->SaveMesh("F:\\temp.defm",true,true);  ///

		SetNeighbors(id1, id2, id3);
	}
}

void CSimplexSurf::UpdateNeighbors2(int pt_index)
{
	int j, k, l, m, o, n[100], n2[100], nb_n = 3, nb_n2; 

	
	//n[0]=pt_index; 
	//n[1]=neighbors[3*pt_index]; 
	//n[2]=neighbors[3*pt_index+1]; 
	//n[3]=neighbors[3*pt_index+2];

	int *q = GetNeighbors(pt_index);
	n[0] = pt_index; 
	n[1] = q[0];
	n[2] = q[1];
	n[3] = q[2];

	//free(neighbors2[pt_index][0]); 
	//neighbors2[pt_index][0] = new int[4]; 
	//neighbors2[pt_index][0][0] = 3; 
	//memcpy(neighbors2[pt_index][0] + 1, neighbors + 3 * pt_index, 3 * sizeof(int));

	//cout << neighbors2[pt_index][0][0] << " " << neighbors2[pt_index][0][1] << " " << neighbors2[pt_index][0][2] << " " << neighbors2[pt_index][0][3] << endl;

	free(neighbors2[pt_index][0]); 
	neighbors2[pt_index][0] = new int[4]; 
	neighbors2[pt_index][0][0] = 3; 
	neighbors2[pt_index][0][1] = q[0]; 
	neighbors2[pt_index][0][2] = q[1]; 
	neighbors2[pt_index][0][3] = q[2]; 
	
	for (j = 1; j < NB_NEIGHBORHOOD; j++) {
		nb_n2 = 0;
		for (k = 0; k < neighbors2[pt_index][j - 1][0]; k++) {
			for (m = 0; m < 3; m++) {
				o = neighbors2[pt_index][j - 1][k + 1];
				for (l = 0; l < nb_n;l++) {
					//if (neighbors[3 * o + m] == n[l]) 
						//l = nb_n + 10; // check if already present

					int *qw = GetNeighbors(o);
					if (qw[m] == n[l]) 
						l = nb_n + 10; // check if already present
				}
				if (l == nb_n) {        // add neighbor in n and n2 
					//n2[nb_n2] = neighbors[3 * o + m]; 
					//nb_n2++;
					//n[nb_n] = neighbors[3 * o + m]; 
					//nb_n++;

					int *qwe = GetNeighbors(o);
					n2[nb_n2] = qwe[m]; 
					nb_n2++;
					n[nb_n] = qwe[m]; 
					nb_n++;
				}
			}
		}
		free(neighbors2[pt_index][j]); 
		neighbors2[pt_index][j] = new int[nb_n2 + 1];
		neighbors2[pt_index][j][0] = nb_n2;
		memcpy(neighbors2[pt_index][j] + 1, n2, nb_n2 * sizeof(int));
	}
}

void CSimplexSurf::OrderNeighbors()
{
int pt,n[3];
for(int i=0;i<nb_points;i++)
	if(neighbors_c[3*i+2]==-1)
		{
		GetNeighbors(i,n);
		pt=-1;
		if(neighbors_c[3*n[0]+2]!=-1) pt=0;
		else if(neighbors_c[3*n[1]+2]!=-1) pt=1;
		else if(neighbors_c[3*n[2]+2]!=-1) pt=2;

		if(pt==-1) // two sides
			{
			if((neighbors_c[3*n[0]]==neighbors_c[3*i+1] && neighbors_c[3*n[0]+1]==neighbors_c[3*i]) || (neighbors_c[3*n[0]]==neighbors_c[3*i] && neighbors_c[3*n[0]+1]==neighbors_c[3*i+1])) pt=0;
			if((neighbors_c[3*n[1]]==neighbors_c[3*i+1] && neighbors_c[3*n[1]+1]==neighbors_c[3*i]) || (neighbors_c[3*n[1]]==neighbors_c[3*i] && neighbors_c[3*n[1]+1]==neighbors_c[3*i+1])) pt=1;
			else pt=2;
			}

		if(pt==0) SetNeighbors(i,n[1],n[2],n[0]);
		else if(pt==1) SetNeighbors(i,n[2],n[0],n[1]);
		}
}

int* CSimplexSurf::GetBorder()
{
bool stop=false; 
int i,ni[3];

nb_points_border=0; for(i=0;i<nb_points;i++) if(neighbors_c[3*i+2]==-1) nb_points_border++;
int* ret=new int[nb_points_border+1];

ret[1]=0; 
while(!stop)
	{
	GetNeighbors(ret[1],ni);
	if(neighbors_c[3*(ret[1])+2]==-1 && (neighbors_c[3*(ni[0])+2]!=-1 || neighbors_c[3*(ni[1])+2]!=-1 || neighbors_c[3*(ni[2])+2]!=-1)) stop=true;
	else ret[1]++;
	}
if(neighbors_c[3*(ni[0])+2]!=-1) ret[2]=ni[2];
else if(neighbors_c[3*(ni[1])+2]!=-1) ret[2]=ni[0];
else ret[2]=ni[1];

ret[0]=2;
while(ret[0]!=nb_points_border)
	{
	GetNeighbors(ret[ret[0]],ni);
	if(ni[0]==ret[ret[0]-1]) ret[ret[0]+1]=ni[1];
	else if(ni[1]==ret[ret[0]-1]) ret[ret[0]+1]=ni[2];
	else if(ni[2]==ret[ret[0]-1]) ret[ret[0]+1]=ni[0];
	ret[0]++;
	}

return ret;
}

P_float CSimplexSurf::GetShortestDist(const int istart,const int iend)
{
P_float* SD=GetShortestDistGraph(istart,false,NULL);
P_float ret=SD[iend];
free(SD);
return ret;
}

P_float CSimplexSurf::GetShortestPath(const int istart,const int iend,int* path)
{
int i,*predecessors=new int[nb_points];
P_float* SD=GetShortestDistGraph(istart,true,predecessors);
P_float ret=SD[iend];
int count=2; 
i=iend;	while(predecessors[i]!=istart) {i=predecessors[i]; count++;}
path=new int[count+1]; path[0]=count; 
path[1]=istart; path[count]=iend; count--;
i=iend; while(predecessors[i]!=istart) {i=predecessors[i]; path[count]=i; count++;}
free(SD); free(predecessors);
return ret;
}

P_float* CSimplexSurf::GetShortestDistGraph(const int index,bool returnpredecessors,int* predecessors)
{
bool* flag=new bool[nb_points];
P_float *ND=GetNeighborsDistGraph(),*ret=new P_float[nb_points],dmin;
int i,iselect,ni[3];

for(i=0;i<nb_points;i++) {ret[i]=1E10; flag[i]=false;}
flag[index]=true;
GetNeighbors(index,ni);
ret[ni[0]]=ND[3*index]; ret[ni[1]]=ND[3*index+1]; ret[ni[2]]=ND[3*index+2];
if(returnpredecessors) {predecessors[ni[0]]=index; predecessors[ni[1]]=index; predecessors[ni[2]]=index;}

bool stop=false,stop2;
while(!stop)
	{
	stop2=true; dmin=1E10;
	for(i=0;i<nb_points;i++)
		if(!flag[i] && ret[i]<dmin)
			{
			dmin=ret[i];
			stop2=false;
			iselect=i;
			}
	if(!stop2)
		{
		flag[iselect]=true;
		GetNeighbors(iselect,ni);
		if(ret[ni[0]]>(ret[iselect]+ND[3*iselect])) {ret[ni[0]]=ret[iselect]+ND[3*iselect]; if(returnpredecessors) predecessors[ni[0]]=iselect;}
		if(ret[ni[1]]>(ret[iselect]+ND[3*iselect+1])) {ret[ni[1]]=ret[iselect]+ND[3*iselect+1]; if(returnpredecessors) predecessors[ni[1]]=iselect;}
		if(ret[ni[2]]>(ret[iselect]+ND[3*iselect+2])) {ret[ni[2]]=ret[iselect]+ND[3*iselect+2]; if(returnpredecessors) predecessors[ni[2]]=iselect;}
		}
	stop=stop2;
	}

free(ND); free(flag);
return ret;
}

P_float* CSimplexSurf::GetNeighborsDistGraph() {
	P_float* ret = new P_float[3 * nb_points];
	for (int i = 0; i < nb_points; i++) { 
		//ret[3 * i] = dist3D(points + 3 * i, points + 3 * neighbors[3 * i]);
		//ret[3 * i + 1] = dist3D(points + 3 * i, points + 3 * neighbors[3 * i + 1]);
		//ret[3 * i + 2] = dist3D(points + 3 * i, points + 3 * neighbors[3 * i + 2]);

		int *n = GetNeighbors(i);
		ret[3 * i] = dist3D(points + 3 * i, points + 3 * n[0]);
		ret[3 * i + 1] = dist3D(points + 3 * i, points + 3 * n[1]);
		ret[3 * i + 2] = dist3D(points + 3 * i, points + 3 * n[2]);
	}
	return ret;
}


void CSimplexSurf::SaveDistCurve(const char* filename,bool isdR)
{
P_float res=0.04;

P_float max,lm,lp,*coords=GetDistCoordinates();
int maxj,i,j,n=(int)floor(1./res); res=1./(P_float)n;
bool flag;

FILE* f=fopen(filename,"wt");
for(i=0;i<n;i++)
	{
    max=0; lp=((P_float)i+1.)*res; lm=(P_float)i*res; flag=false;
	for(j=0;j<nb_points;j++)
		if(coords[2*j+1]<lp &&  coords[2*j+1]>lm)
			{
			if(isdR && axiallinks[j].dR>max) {max=axiallinks[j].dR; maxj=j;}
			if(!isdR && axiallinks[j].Rref>max) {max=axiallinks[j].Rref; maxj=j;}
			flag=true;
			}
			 
	if(flag) fprintf(f,"%lf %lf %d\n",coords[2*maxj+1],max,maxj);
	}
fprintf(f,"\n");
for(i=0;i<n;i++)
	{
    max=0; lp=((P_float)i+1.)*res; lm=(P_float)i*res; flag=false;
	for(j=0;j<nb_points;j++)
		if(coords[2*j]<lp &&  coords[2*j]>lm)
			{
			if(isdR && axiallinks[j].dR>max) {max=axiallinks[j].dR; maxj=j;}
			if(!isdR && axiallinks[j].Rref>max) {max=axiallinks[j].Rref; maxj=j;}
			flag=true;
			}
			 
	if(flag) fprintf(f,"%lf %lf %d\n",coords[2*maxj],max,maxj);
	}


fclose(f);
free(coords);
}

void CSimplexSurf::LoadRadiiRef(const char* filename)
{
if(Axis_model==NULL || axiallinks==NULL) return;
if(Axis_model->GetAxialLinks()==NULL) return;
AXIALLINK *modellinks=Axis_model->GetAxialLinks();

FILE* f=fopen(filename,"rt");
P_float ref; for(int i=0;i<nb_points;i++) {fscanf(f,"%lf\n",&ref); axiallinks[i].dR=ref-axiallinks[i].Rref;}
for(int i=0;i<Axis_model->GetNumberOfPoints();i++)	{	modellinks[i].dR=0; for(int j=0;j<modellinks[i].nb;j++) modellinks[i].dR+=modellinks[i].w[j]*axiallinks[modellinks[i].pts[j]].dR; 	}
fclose(f);
}

void CSimplexSurf::SaveRadii(const char* filename)
{
FILE* f=fopen(filename,"wt");
for(int i=0;i<nb_points;i++) fprintf(f,"%lf\n",axiallinks[i].Rref);
fclose(f);
}

void CSimplexSurf::SaveDistCoordinates(const char* filename,int dim[2],bool normalizeVal,bool isdR)
{
P_float* coords=GetDistCoordinates();
P_float* vals=new P_float[nb_points];

int i;
for(i=0;i<nb_points;i++) if(isdR) vals[i]=axiallinks[i].dR; else vals[i]=axiallinks[i].Rref;

vtkStructuredPoints* im=RBF(nb_points,coords,vals,dim);
vtkStructuredPointsWriter* writer=vtkStructuredPointsWriter::New();
    writer->SetInput(im);
    writer->SetFileTypeToBinary();
    writer->SetFileName(filename);
    writer->Write();
    writer->Delete();

im->Delete();
free(coords);
free(vals);
}

P_float* CSimplexSurf::GetDistCoordinates()
{
int i,j,k,ksav,attachedval_start,attachedval_end;
P_float dmax,dmin,d,d1,d2,**SD=new P_float*[nb_points];
P_float* coords=new P_float[2*nb_points];

// select attachedval_start and attachedval_end
dmax=0;
for(i=0;i<nb_points;i++) 
	if(attached_point[i]!=0) 
		for(j=0;j<nb_points;j++) 
			if(attached_point[j]!=0 && abs(attached_point[j])!=abs(attached_point[i]))
				{
				d=dist3D(GetPoint(i),GetPoint(j));
				if(d>dmax) {dmax=d; attachedval_start=abs(attached_point[j]); attachedval_end=abs(attached_point[i]);}
				}
// compute shortest surface dists
for(i=0;i<nb_points;i++) SD[i]=GetShortestDistGraph(i,false,NULL);

// compute y
for(i=0;i<nb_points;i++) 
	{
	d1=1E10; for(j=0;j<nb_points;j++) if(abs(attached_point[j])==attachedval_start) if(SD[i][j]<d1) d1=SD[i][j];
	d2=1E10; for(j=0;j<nb_points;j++) if(abs(attached_point[j])==attachedval_end) if(SD[i][j]<d2) d2=SD[i][j];
	coords[2*i+1]=d1/(d1+d2);
	}

// compute resolutiony
P_float resolutiony=0;
for(i=0;i<nb_points;i++) 
	if(neighbors_c[3*i+2]==-1)
		{
		if(abs(coords[2*i+1]-coords[2*neighbors[3*i]+1])>resolutiony) resolutiony=abs(coords[2*i+1]-coords[2*neighbors[3*i]+1]);
		if(abs(coords[2*i+1]-coords[2*neighbors[3*i+1]+1])>resolutiony) resolutiony=abs(coords[2*i+1]-coords[2*neighbors[3*i+1]+1]);
		}

int nb_classes=ceil(1./resolutiony);
P_float ry=1./(P_float)nb_classes;
	
// compute x
for(i=0;i<nb_classes;i++)
	{
	dmax=0; dmin=1E10;
	// max distance in a class, for normalization
	for(j=0;j<nb_points;j++)
		if(coords[2*j+1]>i*ry && coords[2*j+1]<(i+1)*ry && neighbors_c[3*j+2]==-1)
			for(k=0;k<nb_points;k++)
				if(coords[2*k+1]>i*ry && coords[2*k+1]<(i+1)*ry && neighbors_c[3*k+2]==-1)
					if(SD[j][k]>dmax) dmax=SD[j][k];
	for(j=0;j<nb_points;j++)
		if(coords[2*j+1]>i*ry && coords[2*j+1]<(i+1)*ry)
			{
			// distance to border 
			for(k=0;k<nb_points;k++) if(coords[2*k+1]>i*ry && coords[2*k+1]<(i+1)*ry && neighbors_c[3*k+2]==-1)		if(SD[j][k]<dmin) {dmin=SD[j][k]; ksav=k;}
			// check side
			if(coords[2*neighbors[3*ksav]+1]>coords[2*neighbors[3*ksav+1]+1]) {dmin=dmax-dmin;} 
			coords[2*j]=dmin/dmax;
			}
	}

for(i=0;i<nb_points;i++) free(SD[i]); free(SD);
return coords;
}

void CSimplexSurf::SetElongationRef(const char* filename)
{
CSimplexSurf* maref;
if(filename!=NULL) {maref=new CSimplexSurf; maref->LoadMesh(filename);} else maref=this;

int i,j,attachedval_start,attachedval_end;
P_float dmax,d,d1ref,d2ref,**SDref=new P_float*[nb_points];

if(Elongations_ref!=NULL) free(Elongations_ref); Elongations_ref=new P_float[nb_points];
if(Elongations_refindex!=NULL) free(Elongations_refindex); Elongations_refindex=new int[2*nb_points];

// select attachedval_start and attachedval_end
dmax=0;
for(i=0;i<nb_points;i++) 
	if(attached_point[i]!=0) 
		for(j=0;j<nb_points;j++) 
			if(attached_point[j]!=0 && abs(attached_point[j])!=abs(attached_point[i]))
				{
				d=dist3D(GetPoint(i),GetPoint(j));
				if(d>dmax) {dmax=d; attachedval_start=abs(attached_point[j]); attachedval_end=abs(attached_point[i]);}
				}
// compute shortest surface dists
for(i=0;i<nb_points;i++) SDref[i]=maref->GetShortestDistGraph(i,false,NULL);

// compute y
for(i=0;i<nb_points;i++) 
	{
	d1ref=1E10; for(j=0;j<nb_points;j++) if(abs(attached_point[j])==attachedval_start) if(SDref[i][j]<d1ref) {d1ref=SDref[i][j]; Elongations_refindex[2*i]=j;}
	d2ref=1E10; for(j=0;j<nb_points;j++) if(abs(attached_point[j])==attachedval_end) if(SDref[i][j]<d2ref)  {d2ref=SDref[i][j]; Elongations_refindex[2*i+1]=j;}
	Elongations_ref[i]=d1ref+d2ref;
	}

for(i=0;i<nb_points;i++) free(SDref[i]); free(SDref);
if(filename!=NULL) delete(maref);
}

P_float CSimplexSurf::GetElongation(int index)
{
if(Elongations_refindex==NULL || Elongations_ref==NULL) return 0;
P_float elongation=GetShortestDist(index,Elongations_refindex[2*index])+GetShortestDist(index,Elongations_refindex[2*index+1])-Elongations_ref[index]; 
return elongation;
}

P_float CSimplexSurf::GetElongationPercentage(int index)
{
if(Elongations_refindex==NULL || Elongations_ref==NULL) return 0;
P_float elongation=(GetShortestDist(index,Elongations_refindex[2*index])+GetShortestDist(index,Elongations_refindex[2*index+1]))/Elongations_ref[index]; 
elongation=100.*elongation - 100.;
return elongation;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
int CSimplexSurf::TO_DeleteCell(int index)
{
int nc[3],count,j,*c=new int[GetCell(index)[0]+1];
c[0]=GetCell(index,c+1);
int *pi=new int[c[0]],*ci=new int[c[0]],*npi=new int[3*c[0]],*nci=new int[9*c[0]];

// delete border points
count=0;
for(j=0;j<c[0];j++)	
	{
	GetNeighbors_c(c[j+1],nc);
	if(nc[0]==-1) {pi[count]=c[j+1]; ci[count]=(nc[1]==index)?nc[2]:nc[1]; GetNeighbors(c[j+1],npi+3*count); count++;}
	if(nc[1]==-1) {pi[count]=c[j+1]; ci[count]=(nc[0]==index)?nc[2]:nc[0]; GetNeighbors(c[j+1],npi+3*count); count++;}
	if(nc[2]==-1) {pi[count]=c[j+1]; ci[count]=(nc[1]==index)?nc[0]:nc[1]; GetNeighbors(c[j+1],npi+3*count); count++;}
	}

//if(count!=2) return 0;
if(c[0]<4) return 0;
for(j=0;j<count;j++) if(GetCell(ci[j])[0]<4) return 0;


for(j=0;j<3*count;j++) SetNeighbors(npi[j],-1,-1,-1);
for(j=0;j<3*count;j++) if(npi[j]==-1) {nci[3*j]=-1; nci[3*j+1]=-1; nci[3*j+2]=-1;} else GetNeighbors_c(npi[j],nci+3*j);
for(j=0;j<count;j++) {DeletePointInCell(ci[j],pi[j]); DeletePointInCell(index,pi[j]);}

//for(j=0;j<9*count;j++) // debug
//	if(nci[j]>nb_cells)
//		int ion=0;


for(j=0;j<9*count;j++) if(nci[j]!=index && nci[j]!=-1) UpdateNeighbors(nci[j]);

// Delete Cell and update Neighbors_c
DeleteCell(index); for(j=0;j<count;j++) if(ci[j]>index) ci[j]--;
for(j=0;j<c[0];j++)	UpdateNeighbors_c(c[j+1]);

for(j=0;j<count;j++)
	{
	DeletePoint(pi[j]);
	for(int k=j+1;k<count;k++) if(pi[k]>pi[j]) pi[k]--;
	}

//for(j=0;j<count;j++) if(GetCell(ci[j])[0]<3) 
//	TO_DeleteCell(ci[j]);

//for(j=0;j<nb_points;j++) // debug
//	if(neighbors[3*j]==-1 || neighbors[3*j+1]==-1 || neighbors[3*j+2]==-1)
//		{
//		int ion=0;
//		int ijd=j;
//		}

free(c); free(ci);free(pi);free(npi);free(nci);
return 1;
}


void CSimplexSurf::DeleteCell(int index)
{
for(int i=0;i<nb_points;i++) 
	{
	if(neighbors_c[3*i]==index) {neighbors_c[3*i]=neighbors_c[3*i+1]; neighbors_c[3*i+1]=neighbors_c[3*i+2]; neighbors_c[3*i+2]=-1;}
	else if(neighbors_c[3*i+1]==index) {neighbors_c[3*i+1]=neighbors_c[3*i]; neighbors_c[3*i]=neighbors_c[3*i+2]; neighbors_c[3*i+2]=-1;}
	else if(neighbors_c[3*i+2]==index) {neighbors_c[3*i+2]=-1;}
	if(neighbors_c[3*i]>index) neighbors_c[3*i]-=1; 
	if(neighbors_c[3*i+1]>index) neighbors_c[3*i+1]-=1; 
	if(neighbors_c[3*i+2]>index) neighbors_c[3*i+2]-=1; 
	}

nb_cells--;

int** buff=(int**)malloc((nb_cells-index)*sizeof(int*));
memcpy(buff,cells+(index+1),(nb_cells-index)*sizeof(int*)); memcpy(cells+index,buff,(nb_cells-index)*sizeof(int*)); cells=(int**)realloc(cells,nb_cells*sizeof(int*));
free(buff);

int* buffi=(int*)malloc((nb_cells-index)*sizeof(int));
memcpy(buffi,attached_cell+(index+1),(nb_cells-index)*sizeof(int)); memcpy(attached_cell+index,buffi,(nb_cells-index)*sizeof(int)); attached_cell=(int*)realloc(attached_cell,nb_cells*sizeof(int));  
free(buffi);

P_float* bufff=(P_float*)malloc(3*(nb_cells-index)*sizeof(P_float));
memcpy(bufff,cellcenters+3*(index+1),3*(nb_cells-index)*sizeof(P_float)); memcpy(cellcenters+3*index,bufff,3*(nb_cells-index)*sizeof(P_float)); cellcenters=(P_float*)realloc(cellcenters,3*nb_cells*sizeof(P_float));  
memcpy(bufff,Val_c+3*(index+1),3*(nb_cells-index)*sizeof(P_float)); memcpy(Val_c+3*index,bufff,3*(nb_cells-index)*sizeof(P_float)); Val_c=(P_float*)realloc(Val_c,3*nb_cells*sizeof(P_float));  
free(bufff);

surfaces_cell=(P_float*)realloc(surfaces_cell,nb_cells*sizeof(P_float));  
}

void CSimplexSurf::InsertNextCell(int* f)
{
nb_cells++; cells=(int**)realloc(cells,nb_cells*sizeof(int*)); cells[nb_cells-1]=NULL;
SetCell(nb_cells-1,f);

cellcenters=(P_float*)realloc(cellcenters,3*nb_cells*sizeof(P_float));  
attached_cell=(int*)realloc(attached_cell,nb_cells*sizeof(int)); attached_cell[nb_cells-1]=0; 
Val_c=(P_float*)realloc(Val_c,3*nb_cells*sizeof(P_float));  P_float O[3]={0,0,0}; SetVal_c(nb_cells-1,O);
surfaces_cell=(P_float*)realloc(surfaces_cell,nb_cells*sizeof(P_float));  
}

void CSimplexSurf::SearchCell(int p,int c[3]) {
	int count=0; c[0] = -1; c[1] = -1; c[2] = -1;
	
	for (int i = 0; i < nb_cells; i++) { 
		for (int j = 0; j < cells[i][0]; j++) { 
			if (cells[i][j + 1] == p) {
				c[count] = i; 
				count++;
			}
		}
	}
	//if(count==2) 
	//c[2]=-1;
	//if(count==1) // debug
	//	{c[1]=-1; c[2]=-1;}
}

int CSimplexSurf::SearchCell(int p1,int p2,int c[2][2]) 
{
// c[0]: not adjacent cells
// c[1]: adjacent cells
int c1[3],c2[3];
SearchCell(p1,c1);
SearchCell(p2,c2);
int count =0;


//if(p1!=neighbors[3*p2] && p1!=neighbors[3*p2+1] && p1!=neighbors[3*p2+2]) // not adjacent points


int *n = GetNeighbors(p2);
if(p1!=n[0] && p1!=n[1] && p1!=n[2]) // not adjacent points
    {
    if(c1[0]==c2[0] || c1[0]==c2[1] || c1[0]==c2[2]) c[1][0]=c1[0];
    if(c1[1]==c2[0] || c1[1]==c2[1] || c1[1]==c2[2]) c[1][0]=c1[1];
    if(c1[2]==c2[0] || c1[2]==c2[1] || c1[2]==c2[2]) c[1][0]=c1[2];
    return 1;
    }

if(c1[0]!=c2[0] && c1[0]!=c2[1] && c1[0]!=c2[2]) c[0][0]=c1[0]; else if(count<2) {c[1][count]=c1[0];count++;}
if(c1[1]!=c2[0] && c1[1]!=c2[1] && c1[1]!=c2[2]) c[0][0]=c1[1]; else if(count<2) {c[1][count]=c1[1];count++;}
if(c1[2]!=c2[0] && c1[2]!=c2[1] && c1[2]!=c2[2]) c[0][0]=c1[2]; else if(count<2) {c[1][count]=c1[2];count++;}
if(c2[0]!=c1[0] && c2[0]!=c1[1] && c2[0]!=c1[2]) c[0][1]=c2[0]; 
if(c2[1]!=c1[0] && c2[1]!=c1[1] && c2[1]!=c1[2]) c[0][1]=c2[1];
if(c2[2]!=c1[0] && c2[2]!=c1[1] && c2[2]!=c1[2]) c[0][1]=c2[2]; 
return 0;
}

int CSimplexSurf::SearchCell(int p1, int p2,int p3) 
{
int c1[3];SearchCell(p1,c1);
int c2[3];SearchCell(p2,c2);
int c3[3];SearchCell(p3,c3);

if((c1[0]==c2[0] || c1[0]==c2[1] || c1[0]==c2[2])&&(c1[0]==c3[0] || c1[0]==c3[1] || c1[0]==c3[2])) return c1[0]; 
else if((c1[1]==c2[0] || c1[1]==c2[1] || c1[1]==c2[2])&&(c1[1]==c3[0] || c1[1]==c3[1] || c1[1]==c3[2])) return c1[1]; 
else return c1[2]; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::SaveMesh(const char* filename, bool writeneighbors, bool writemergedpoints) {
	int i, j, k;
	FILE* f = fopen(filename, "wt");
	fprintf(f, "%d %d %4.2f\n", nb_points, nb_cells, volume);
	
	for (i = 0; i < nb_points; i++) 
		fprintf(f, "%4.2f %4.2f %4.2f\n", points[3 * i], points[3 * i + 1], points[3 * i + 2]);
	
	for (i = 0; i < nb_cells; i++) {
		fprintf(f, "\n"); 
		
		for (j = 0; j < cells[i][0] + 1; j++) 
			fprintf(f, "%d ", cells[i][j]); 
	}
	 //fprintf(f,"\n"); for(i=0;i<nb_cells;i++) {fprintf(f,"\n"); fprintf(f,"%d ",flagmultires[i]);}

	if (writeneighbors) {
		fprintf(f, "\n\n%d\n", NB_NEIGHBORHOOD);
		fprintf(f, "\n");	
		
		for (i = 0; i < nb_points; i++) {
			//fprintf(f, "%d %d %d\n", neighbors[3 * i], neighbors[3 * i + 1], neighbors[3 * i + 2]);
			int *nn = GetNeighbors(i);
			fprintf(f, "%d %d %d\n", nn[0], nn[1], nn[2]);
		}

		ofstream file1, file2;
		file1.open("neighbors_c.txt");
		
		fprintf(f, "\n");	
		for (i = 0; i < nb_points; i++) {
			fprintf(f, "%d %d %d\n", neighbors_c[3 * i], neighbors_c[3 * i + 1], neighbors_c[3 * i + 2]);
			file1 << neighbors_c[3 * i] << ", " << neighbors_c[3 * i + 1] << ", " << neighbors_c[3 * i + 2] << endl;
		}
		file1.close();
		fprintf(f, "\n\n\n\n");

		file2.open("neighbors2.txt");
		for (i = 0; i < nb_points; i++) {
			for (j = 0; j < NB_NEIGHBORHOOD; j++) {	
				fprintf(f, "\n    ");
				
				for (k = 0; k <= neighbors2[i][j][0]; k++) {
					fprintf(f, "%d ", neighbors2[i][j][k]);

					file2 << "[" << i << ", " << j << ", " << k << "]: " << neighbors2[i][j][k] << ", ";
				}
				file2 << endl;
			}
		}
		file2.close();
	}

	if (writemergedpoints && MergedPoints != NULL) {
		fprintf(f, "\n\n");
		
		for (i = 0; i < 2 * MergedPoints[0] + 1; i++) 
			fprintf(f, "%d ", MergedPoints[i]);
	}

	fclose(f);
}

void CSimplexSurf::SaveParams(const char* filename)
{
int i;
FILE* f=fopen(filename,"wt");
fprintf(f,"%d\n",nb_points);
for(i=0;i<nb_points;i++) fprintf(f,"%4.2f %4.2f %4.2f\n",params[3*i],params[3*i+1],params[3*i+2]);
fclose(f);
}

void CSimplexSurf::LoadParams(const char* filename)
{
int i;
FILE* f=fopen(filename,"rt");
int tmp; fscanf(f,"%d\n",&tmp);
for(i=0;i<nb_points;i++) fscanf(f,"%lf %lf %lf\n",params+3*i,params+3*i+1,params+3*i+2);
fclose(f);
}


void CSimplexSurf::SaveOFFMesh(const char* filename,int lift)
{
vtkPolyData* polydata=GetPolyData(lift);
int i;
FILE* f=fopen(filename,"wt");
fprintf(f,"OFF\n%d %d %d\n",polydata->GetNumberOfPoints(),polydata->GetNumberOfCells(),polydata->GetNumberOfPoints()+polydata->GetNumberOfCells());
for(i=0;i<polydata->GetNumberOfPoints();i++) fprintf(f,"%4.2f %4.2f %4.2f\n",polydata->GetPoint(i)[0],polydata->GetPoint(i)[1],polydata->GetPoint(i)[2]);
for(i=0;i<polydata->GetNumberOfCells();i++)  fprintf(f,"3 %d %d %d\n",polydata->GetCell(i)->GetPointId(0),polydata->GetCell(i)->GetPointId(1),polydata->GetCell(i)->GetPointId(2));
fclose(f);
polydata->Delete();
}

void CSimplexSurf::LoadMesh(const char* filename)
{
int i,j,k;
FILE* f=fopen(filename,"rt");
Free();
fscanf(f,"%d %d %lf\n",&nb_points,&nb_cells,&volume);
Allocate(nb_points,nb_cells);
for(i=0;i<nb_points;i++) fscanf(f,"%lf %lf %lf\n",&points[3*i],&points[3*i+1],&points[3*i+2]);
int *cellp,n;
for(i=0;i<nb_cells;i++) {fscanf(f,"%d",&n); cellp=new int[n]; for(j=0;j<n;j++) fscanf(f,"%d ",&cellp[j]); SetCell(i,n,cellp); free(cellp);}
//int temp; for(i=0;i<nb_cells;i++) fscanf(f,"%d",&temp); // flagmultires

// Neighbors
if(!feof(f))
	{
    int nb,nb_neighbors; fscanf(f,"%d",&nb_neighbors);
	for(i=0;i<nb_points;i++) {
		//fscanf(f,"%d %d %d ",&neighbors[3*i],&neighbors[3*i+1],&neighbors[3*i+2]);
		int n1, n2, n3;
		fscanf(f,"%d %d %d ",&n1,&n2,&n3);
		SetNeighbors(i, n1, n2, n3);
	}
	for(i=0;i<nb_points;i++) fscanf(f,"%d %d %d ",&neighbors_c[3*i],&neighbors_c[3*i+1],&neighbors_c[3*i+2]); 
	if(nb_neighbors==NB_NEIGHBORHOOD)
		for(i=0;i<nb_points;i++) 
			{
			for(j=0;j<NB_NEIGHBORHOOD;j++)
				{	
				fscanf(f,"%d ",&nb);
				free(neighbors2[i][j]); neighbors2[i][j]=new int[nb+1]; neighbors2[i][j][0]=nb;
				for(k=0;k<neighbors2[i][j][0];k++) fscanf(f,"%d ",&neighbors2[i][j][k+1]);
				}
			}
	else for(i=0;i<nb_points;i++) UpdateNeighbors2(i);
	}
else UpdateNeighbors();

// axis?
bool flag=false; for(i=0;i<nb_points;i++) if(neighbors_c[3*i+2]==-1) flag=true; 
if(flag) SetIsAxis(true);

// MergedPoints
if(!feof(f))
	{
	fscanf(f,"\n\n");
	fscanf(f,"%d ",&i); MergedPoints=new int[2*i+1];  MergedPoints[0]=i;
	for(i=0;i<MergedPoints[0];i++) fscanf(f,"%d %d ",MergedPoints+2*i+1,MergedPoints+2*i+2);
	}

// nb_border_points
nb_points_border=0; for(i=0;i<nb_points;i++) if(neighbors_c[3*i+2]==-1) nb_points_border++;
if(nb_points_border!=0) OrderNeighbors();

UpdateFlipped();

Equilibrium();
UpdateParams();
UpdateMass();  
volume_ref=volume;
surface_ref=surface;
fclose(f);
}


void CSimplexSurf::RegisterToSimplexMesh(const char * filename)
{
FILE* f=fopen(filename,"rt");
P_float v; int i,nbp,nbc; fscanf(f,"%d %d %lf\n",&nbp,&nbc,&v);
P_float *p=new P_float[3*nbp];
for(i=0;i<nbp;i++) fscanf(f,"%lf %lf %lf\n",&p[3*i],&p[3*i+1],&p[3*i+2]);
fclose(f);
P_float dmin,d; int j,jmin;
for(i=0;i<nb_points;i++)
	{
	dmin=1E10;
	for(j=0;j<nbp;j++) {d=dist3D(points+3*i,p+3*j); if(d<dmin) {dmin=d; jmin=j;}}
	SetPoint(i,p+3*jmin);
	}
free(p);
Equilibrium();
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateUnstructuredGrid() {UGrid->Delete();UGrid=vtkUnstructuredGrid::New();  UGrid->Allocate(nb_points,nb_cells);UpdateUnstructuredGridPoints();UpdateUnstructuredGridCells();}
void CSimplexSurf::UpdateUnstructuredGridCells()
{
int i,j;
vtkPolygon* polygon=vtkPolygon::New();
for(i=0;i<nb_cells;i++) 
    { 
    polygon->GetPointIds()->SetNumberOfIds(cells[i][0]);
    for(j=0;j<cells[i][0];j++) polygon->GetPointIds()->SetId(j,cells[i][j+1]);
    UGrid->InsertNextCell(VTK_POLYGON,polygon->GetPointIds());
    }
polygon->Delete();
}

vtkPolyData* CSimplexSurf::GetPolyData(int lift)
{
if(lift!=-1) LiftCellcenters(lift);

int i,j;
P_float p[3];

j=0; for(i=0;i<nb_cells;i++) j+=cells[i][0]; 

vtkPolyData* PolyData=vtkPolyData::New();
vtkPoints* pts=vtkPoints::New();
    pts->SetNumberOfPoints(nb_points+nb_cells);
PolyData->SetPoints(pts);

PolyData->Allocate(j,j);

vtkTriangle* triangle=vtkTriangle::New();

for(i=0;i<nb_points;i++) 
    {
    GetPoint(i,p); 
    pts->SetPoint(i,p);
    } // copy points


for(i=0;i<nb_cells;i++) 
    { 
    triangle->GetPointIds()->SetId(0,nb_points+i);
    for(j=0;j<cells[i][0];j++) 
        {
        triangle->GetPointIds()->SetId(1,cells[i][j+1]);
        if(j==cells[i][0]-1) triangle->GetPointIds()->SetId(2,cells[i][1]); else triangle->GetPointIds()->SetId(2,cells[i][j+2]);
        PolyData->InsertNextCell(triangle->GetCellType(),triangle->GetPointIds()); // Insert triangle
        }
    pts->SetPoint(nb_points+i,cellcenters+3*i); // insert cell center
    }

PolyData->SetPoints(pts);
PolyData->Update();
triangle->Delete();

vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
    normals->SetInput(PolyData);
    normals->SplittingOff();
    normals->ConsistencyOn();
    normals->ComputePointNormalsOn();
    normals->Update();

PolyData->Delete();
PolyData=vtkPolyData::New();
    PolyData->DeepCopy(normals->GetOutput());
normals->Delete();

//if(lift!=-1) MapCurvatureOnPolyData(PolyData,5);
if(IsAxis) MapRadiusOnPolyData(PolyData);

return(PolyData);
}


void CSimplexSurf::MapCurvatureOnPolyData(vtkPolyData* model,int reject_percentage)
{
int i,nbp=model->GetNumberOfPoints(); 
P_float Kg_Kh_S[3],n[3],nm[3];
double* uca=new double[nbp]; for(i=0;i<nbp;i++) uca[i]=0;

for(i=0;i<nb_points;i++) {GetTriCurvature_Point(i,Kg_Kh_S,n); uca[i]=Kg_Kh_S[1];}
for(i=0;i<nb_cells;i++) {GetTriCurvature_Cell(i,Kg_Kh_S,n,nm); uca[i+nb_points]=Kg_Kh_S[1];}

// normalisation
P_float max=-1E10,min=1E10; for(i=0;i<nbp;i++) {if(uca[i]>max) max=uca[i]; if(uca[i]<min) min=uca[i];}
for(i=0;i<nbp;i++) uca[i]=(uca[i]-min)/(max-min);

Reject_Histo(nbp,uca,reject_percentage);

vtkDoubleArray* da=vtkDoubleArray::New();	da->SetArray(uca,nbp,0);
if(model->GetPointData()->GetScalars()!=NULL) model->GetPointData()->GetScalars()->Delete();	model->GetPointData()->SetScalars(da);
}


void CSimplexSurf::MapRadiusOnPolyData(vtkPolyData* model)
{
if(axiallinks==NULL) return;
int i,j;
double *uca=new double[nb_points+nb_cells]; 
for(i=0;i<nb_points;i++) uca[i]=axiallinks[i].Rref;
for(i=0;i<nb_cells;i++) 
	{
	uca[nb_points+i]=0;     
	for(j=0;j<cells[i][0];j++) uca[nb_points+i]+=axiallinks[cells[i][j+1]].Rref;
    uca[nb_points+i]=uca[nb_points+i]/(P_float)cells[i][0]; 
	}

vtkDoubleArray* da=vtkDoubleArray::New();	da->SetArray(uca,nb_points+nb_cells,0);
if(model->GetPointData()->GetScalars()!=NULL) model->GetPointData()->GetScalars()->Delete();	model->GetPointData()->SetScalars(da);
}

P_float CSimplexSurf::LiftCellcenters(int nb_it)
{
int i,j,k;
P_float Cm,C,Kg_Kh_S[3],n[3],nm[3],dev,mean,means;

P_float *steps=new P_float[nb_cells]; for(i=0;i<nb_cells;i++) steps[i]=1;

//FILE* f=fopen("C:\\Temp\\data2\\generic_male2\\tests\\testbed\\lifting.txt","rt+"); // DEBUG
//if(f==NULL) f=fopen("C:\\Temp\\data2\\generic_male2\\tests\\testbed\\lifting.txt","wt"); // DEBUG
//else fseek(f,0,SEEK_END); // DEBUG

for(k=0;k<nb_it;k++)
	{
	mean=0; dev=0; means=0; 
//	P_float scheck=0;	for(i=0;i<nb_cells;i++) {GetTriCurvature_Cell(i,Kg_Kh_S,n,nm); scheck+=abs(Kg_Kh_S[0])*Kg_Kh_S[2]; } 	for(i=0;i<nb_points;i++) {GetTriCurvature_Point(i,Kg_Kh_S,n); scheck+=abs(Kg_Kh_S[0])*Kg_Kh_S[2]; }

	for(i=0;i<nb_cells;i++)
		{
		// compute mean curvature
		Cm=0;	for(j=0;j<cells[i][0];j++) {GetTriCurvature_Point(cells[i][j+1],Kg_Kh_S,n); Cm+=Kg_Kh_S[1]; } Cm=Cm/(P_float)cells[i][0];
		// compute curvature
		GetTriCurvature_Cell(i,Kg_Kh_S,nm,n); C=Kg_Kh_S[1]; 
		// lift
		if(k==0) {if((Cm-C)*steps[i]<0) steps[i]=-steps[i];} else {if((Cm-C)*steps[i]<0) steps[i]=-steps[i]/2;}
		cellcenters[3*i]+=steps[i]*n[0];	cellcenters[3*i+1]+=steps[i]*n[1];	cellcenters[3*i+2]+=steps[i]*n[2];
		mean+=abs(Cm-C); dev+=(Cm-C)*(Cm-C); means+=abs(steps[i]);
		}
	mean=mean/(P_float)(nb_cells);	means=means/(P_float)(nb_cells);
	dev=sqrt(dev/(P_float)(nb_cells)-mean*mean);

//	fprintf(f,"%lf %lf %lf\n",mean,dev,means); // DEBUG
	}

//fclose(f); // DEBUG

free(steps);

// return curvature deviation
return dev;
}

P_float CSimplexSurf::GetTriCurvature_Point(int index,P_float Kg_Kh_S[3],P_float n[3])
{
int *ne=neighbors+3*index,*nc=neighbors_c+3*index;
int nb; if(nc[2]==-1) nb=5; else nb=6; P_float *pn=new P_float[3*nb]; 
int order[3]={-2,-2,-2};
if((nc[0]==neighbors_c[3*ne[0]] || nc[0]==neighbors_c[3*ne[0]+1] || nc[0]==neighbors_c[3*ne[0]+2])	&& (nc[0]==neighbors_c[3*ne[1]] || nc[0]==neighbors_c[3*ne[1]+1] || nc[0]==neighbors_c[3*ne[1]+2]))		order[0]=nc[0];	else if((nc[1]==neighbors_c[3*ne[0]] || nc[1]==neighbors_c[3*ne[0]+1] || nc[1]==neighbors_c[3*ne[0]+2]) && (nc[1]==neighbors_c[3*ne[1]] || nc[1]==neighbors_c[3*ne[1]+1] || nc[1]==neighbors_c[3*ne[1]+2]))		order[0]=nc[1];		else order[0]=nc[2];
if((nc[0]==neighbors_c[3*ne[1]] || nc[0]==neighbors_c[3*ne[1]+1] || nc[0]==neighbors_c[3*ne[1]+2]) && (nc[0]==neighbors_c[3*ne[2]] || nc[0]==neighbors_c[3*ne[2]+1] || nc[0]==neighbors_c[3*ne[2]+2])) order[1]=nc[0];	else if((nc[1]==neighbors_c[3*ne[1]] || nc[1]==neighbors_c[3*ne[1]+1] || nc[1]==neighbors_c[3*ne[1]+2])	&& (nc[1]==neighbors_c[3*ne[2]] || nc[1]==neighbors_c[3*ne[2]+1] || nc[1]==neighbors_c[3*ne[2]+2])) order[1]=nc[1];		else order[1]=nc[2];
if(nc[0]!=order[0] && nc[0]!=order[1]) order[2]=nc[0];	else if(nc[1]!=order[0] && nc[1]!=order[1]) order[2]=nc[1];		else order[2]=nc[2];
memcpy(pn,points+3*ne[0],3*sizeof(P_float));	memcpy(pn+3,cellcenters+3*order[0],3*sizeof(P_float));	memcpy(pn+6,points+3*ne[1],3*sizeof(P_float));	memcpy(pn+9,cellcenters+3*order[1],3*sizeof(P_float)); memcpy(pn+12,points+3*ne[2],3*sizeof(P_float)); if(order[2]!=-1) memcpy(pn+15,cellcenters+3*order[2],3*sizeof(P_float));
GetCurvatureTriangle(Kg_Kh_S,n,points+3*index,nb,pn,normals+3*index);
free(pn);
return Kg_Kh_S[1];
}

P_float CSimplexSurf::GetTriCurvature_Cell(int index,P_float Kg_Kh_S[3],P_float n[3],P_float nm[3])
{
int nb=cells[index][0]; 
P_float *pn=new P_float[3*nb]; 
// get mean normal (to check curv sign)
nm[0]=0; nm[1]=0; nm[2]=0; for(int j=0;j<cells[index][0];j++)	{memcpy(pn+3*j,points+3*cells[index][j+1],3*sizeof(P_float)); nm[0]+=normals[3*cells[index][j+1]]; nm[1]+=normals[3*cells[index][j+1]+1]; nm[2]+=normals[3*cells[index][j+1]+2];} P_float nrm=norm(nm); nm[0]=nm[0]/nrm; nm[1]=nm[1]/nrm; nm[2]=nm[2]/nrm;
GetCurvatureTriangle(Kg_Kh_S,n,cellcenters+3*index,nb,pn,nm);
free(pn);
return Kg_Kh_S[1];
}

P_float CSimplexSurf::GetCurvature(int index) {	
	//P_float C = 1. / GetCircumscribedRadius(points + 3 * index, points + 3 * neighbors[3 * index], points + 3 * neighbors[3 * index + 1], points + 3 * neighbors[3 * index + 2]);
	int *n = GetNeighbors(index);
	P_float C = 1. / GetCircumscribedRadius(points + 3 * index, points + 3 * n[0], points + 3 * n[1], points + 3 * n[2]);
	if (h[index] > 0) C = -C;
	return C;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateAll() {
	P_float n[3], v[3], *P, *P1, *P2, *P3, P1P3[3], P2P3[3], P1P2[3], PP1[3], s, dp;

	for (int i = 0; i < nb_points; i++) {
		P = GetPoint(i); 
		//P1 = GetPoint(neighbors[3 * i]); 
		//P2 = GetPoint(neighbors[3 * i + 1]); 
		//P3 = GetPoint(neighbors[3 * i + 2]); 
		
		int *pointNeighbors = GetNeighbors(i);
		P1 = GetPoint(pointNeighbors[0]); 
		P2 = GetPoint(pointNeighbors[1]); 
		P3 = GetPoint(pointNeighbors[2]); 

		P1P3[0] = P3[0] - P1[0];
		P1P3[1] = P3[1] - P1[1]; 
		P1P3[2] = P3[2] - P1[2];
		
		P2P3[0] = P3[0] - P2[0];
		P2P3[1] = P3[1] - P2[1];
		P2P3[2] = P3[2] - P2[2];
		
		PP1[0] = P1[0] - P[0];
		PP1[1] = P1[1] - P[1]; PP1[2] = P1[2] - P[2];


		// normals2
		if (neighbors_c[3 * i + 2] == -1) 	{
			P1P2[0] = P2[0] - P1[0]; 
			P1P2[1] = P2[1] - P1[1];
			P1P2[2] = P2[2] - P1[2]; 

			s = dotproduct(P1P2, P1P2); 
			dp = dotproduct(PP1, P1P2); 
			n[0] = -PP1[0] + dp * P1P2[0] / s; 
			n[1] = -PP1[1] + dp * P1P2[1] / s; 
			n[2] = -PP1[2] + dp * P1P2[2] / s;	
			
			s = norm(n);  
			normals2[3 * i] = n[0] / s; 
			normals2[3 * i + 1] = n[1] / s; 
			normals2[3 * i + 2] = n[2] / s;   
		}


		// normals
		crossproduct(n, P2P3, P1P3); 
		s = norm(n); 

		if (s == 0) { // colinear neighbors 
			
			//n=p1p3 n (p1p3 n pp1)
			crossproduct(v, P1P3, PP1); 
			crossproduct(n, P1P3, v);
			s = norm(n);
			if (s == 0) 
				SetNormal(i, 1, 0, 0); // colinear neighbors and point
			else 
				SetNormal(i, n[0] / s, n[1] / s, n[2] / s);
		}
		else 
			SetNormal(i, n[0] / s, n[1] / s, n[2] / s);	

	//debug
	//	if(neighbors_c[3*i+2]==-1) 	{if(dotproduct(normals2+3*i,normals+3*i)<0) {normals2[3*i]=-normals2[3*i]; normals2[3*i+1]=-normals2[3*i+1]; normals2[3*i+2]=-normals2[3*i+2]; } }


		// h
		if (neighbors_c[3 * i + 2] == -1) 
			h[i] = dotproduct(PP1, normals2 + 3 * i); 
		else 
			h[i] = dotproduct(PP1, normals + 3 * i); 

	//	volumes[i]=abs(h[i])*s/6.;
	}
	UpdateCellCenters();
	ComputeVolume();
}


/////////////////////////////////////////////////////////////////////////////////////////////////
P_float CSimplexSurf::GetHdiff(int index) {
	P_float P[3], P1[3], P2[3], P3[3], s;
	GetPoint(index, P); 

	//GetPoint(neighbors[3 * index], P1); 
	//GetPoint(neighbors[3 * index + 1], P2);  
	//GetPoint(neighbors[3 * index + 2], P3); 

	int *pointNeighbors = GetNeighbors(index);
	GetPoint(pointNeighbors[0], P1); 
	GetPoint(pointNeighbors[1], P2);  
	GetPoint(pointNeighbors[2], P3); 

	P_float PP1[3] = {P1[0] - P[0], P1[1] - P[1], P1[2] - P[2]};
	P_float *n; 
	if(neighbors_c[3 * index + 2] == -1) {
		P_float P1P2[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
		s = dotproduct(P1P2, P1P2); 
		n = normals2 + 3 * index;
	}
	else {
		P_float P1P3[3] = {P3[0] - P1[0], P3[1] - P1[1], P3[2] - P1[2]}; 
		P_float P2P3[3] = {P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2]}; 
		P_float cp[3]; 
		crossproduct(cp, P2P3, P1P3); 
		s = norm(cp); 
		n = GetNormal(index);
	}

	P_float h = dotproduct(PP1, n); 

	return h - sqrt(s) / params[3 * index + 2];	
}

void CSimplexSurf::UpdateParams() {
	for (int i = 0; i < nb_points; i++) 
		UpdateParams(i);
}

void CSimplexSurf::UpdateParams(int index) {
	P_float P[3], e[3], P1[3], P2[3], P3[3], s;
	GetPoint(index, P); 
	//GetPoint(neighbors[3 * index], P1); 
	//GetPoint(neighbors[3 * index + 1], P2);  
	//GetPoint(neighbors[3 * index + 2], P3); 
	
	int *pointNeighbors = GetNeighbors(index);
	GetPoint(pointNeighbors[0], P1); 
	GetPoint(pointNeighbors[1], P2);  
	GetPoint(pointNeighbors[2], P3); 

	P_float PP1[3] = {P1[0] - P[0], P1[1] - P[1], P1[2] - P[2]};
	// ei params: barycentric coordinates sum(ei.PiPproj)=0 sum(ei)=1 -> C=e1P1+e2P2+(1-e1-e2)P3
	// params(3)=sqrt(surf)/h
	P_float *n; 

	if (Forcehandleborder && neighbors_c[3 * index + 2] == -1) {
		P_float P1P2[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
		s = dotproduct(P1P2, P1P2); 
		n = normals2 + 3 * index;
	}
	else {
		P_float P1P3[3] = {P3[0] - P1[0], P3[1] - P1[1], P3[2] - P1[2]};
		P_float P2P3[3] = {P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2]}; 
		
		P_float cp[3]; 
		crossproduct(cp, P2P3, P1P3); 
		s = norm(cp); 
		n = GetNormal(index);
	}

	h[index] = dotproduct(PP1, n); 
	P_float Pproj[3]; 
	Pproj[0] = P[0] + h[index] * n[0]; 
	Pproj[1] = P[1] + h[index] * n[1];
	Pproj[2] = P[2] + h[index] * n[2];

	if (h[index] == 0) 
		h[index] = 1E-10; 

	Barycenter(e, Pproj, P1, P2, P3); 
	params[3 * index] = e[0];
	params[3 * index + 1] = e[1]; 
	params[3 * index + 2] = pow(s, SCALEINVARIANT_REFSHAPE) / h[index];	 
}

/*
void CSimplexSurf::UpdateParams_multires(CSimplexSurf* mesh)
{
int i,j,nb,c[200],pt_index=mesh->GetNumberOfPoints();
P_float P[3],P1[3],P2[3],P3[3],n[3],P1P3[3],P2P3[3],PP1[3],Pproj[3],e[3],s,h;
for(i=0;i<mesh->GetNumberOfPoints();i++) {params_multires[3*i]=0; params_multires[3*i+1]=1; params_multires[3*i+2]=0;}

for(i=0;i<mesh->GetNumberOfCells();i++)
    {
    nb=mesh->GetCell(i,c);
    for(j=0;j<nb;j++)
        {
        GetPoint(pt_index,P);
        mesh->GetPoint(c[(j==0)?(nb-1):(j-1)],P1);      mesh->GetPoint(c[j],P2);        mesh->GetPoint(c[(j==nb-1)?(0):(j+1)],P3); 

        P1P3[0]=P3[0]-P1[0];P1P3[1]=P3[1]-P1[1];P1P3[2]=P3[2]-P1[2];
        P2P3[0]=P3[0]-P2[0];P2P3[1]=P3[1]-P2[1];P2P3[2]=P3[2]-P2[2];
        PP1[0]=P1[0]-P[0];PP1[1]=P1[1]-P[1];PP1[2]=P1[2]-P[2];
        crossproduct(n,P2P3,P1P3); s=norm(n); if(s==0) {n[0]=0; n[1]=0;n[2]=0;} else {n[0]=n[0]/s; n[1]=n[1]/s; n[2]=n[2]/s;}
        h=dotproduct(PP1,n); 

        Pproj[0]=P[0]+h*n[0]; Pproj[1]=P[1]+h*n[1]; Pproj[2]=P[2]+h*n[2];
        Barycenter(e,Pproj,P1,P2,P3);
        params_multires[3*pt_index]=e[0]; params_multires[3*pt_index+1]=e[1]; params_multires[3*pt_index+2]=sqrt(s)/h;

        pt_index++;
        }
    }
}*/

P_float* CSimplexSurf::GetParams(int index) {return params+3*index; }
void CSimplexSurf::GetParams(int index,P_float p[3]) {memcpy(p,params+3*index,3*sizeof(P_float)); }
void CSimplexSurf::SetParams(int index,P_float p[3]) {memcpy(params+3*index,p,3*sizeof(P_float)); }
void CSimplexSurf::SetParams(int index,P_float p1,P_float p2,P_float p3) {params[3*index]=p1;params[3*index+1]=p2;params[3*index+2]=p3; }

P_float CSimplexSurf::GetHmean(int index)
{
P_float H=0;
// compute h=average(hi) i E neighborhood
int nb_val=0,n;

if(IsAxis) 
	{
	if(neighbors_c[3*index+2]==-1)  
		{
		//compute n (projection of n2 onto the neighbors plane)
		if(Forcehandleborder) {for(int j=0;j<NB_NEIGHBORHOOD;j++)	for(int k=0;k<neighbors2[index][j][0];k++) {n=neighbors2[index][j][k+1];	if(neighbors_c[3*n+2]==-1) {H+=h[n]; nb_val++;}}}
		else {H=0; nb_val=0;}
		}
	else 	{for(int j=0;j<NB_NEIGHBORHOOD;j++)	for(int k=0;k<neighbors2[index][j][0];k++)	{n=neighbors2[index][j][k+1]; if(neighbors_c[3*n+2]!=-1) {H+=h[n]; nb_val++;}}}
	}
else {for(int j=0;j<NB_NEIGHBORHOOD;j++)	for(int k=0;k<neighbors2[index][j][0];k++)		{n=neighbors2[index][j][k+1];	H+=h[n]; nb_val++; }}
if(nb_val!=0) H=H/(P_float)nb_val;

return H;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateNormals() {
	for (int i = 0; i < GetNumberOfPoints(); i++) {
		P_float n[3], v[3], *P, *P1, *P2, *P3, P1P3[3], P2P3[3], P1P2[3], PP1[3], s, dp;

		P = GetPoint(i); 
		//P1 = GetPoint(neighbors[3 * i]); 
		//P2 = GetPoint(neighbors[3 * i + 1]);  
		//P3 = GetPoint(neighbors[3 * i + 2]); 
				
		// If the point is a regular simplex point with only 3 neighbors
		if (IsMultiMaterialPoint(i) == false) {
			int *pointNeighbors = GetNeighbors(i);

			P1 = GetPoint(pointNeighbors[0]); 
			P2 = GetPoint(pointNeighbors[1]);  
			P3 = GetPoint(pointNeighbors[2]); 

			P1P3[0] = P3[0] - P1[0];
			P1P3[1] = P3[1] - P1[1];
			P1P3[2] = P3[2] - P1[2];

			P2P3[0] = P3[0] - P2[0];
			P2P3[1] = P3[1] - P2[1];
			P2P3[2] = P3[2] - P2[2];

			PP1[0] = P1[0] - P[0];
			PP1[1] = P1[1] - P[1];
			PP1[2] = P1[2] - P[2];

			// normals2
			if (neighbors_c[3 * i + 2] == -1) {
				P1P2[0] = P2[0] - P1[0];
				P1P2[1] = P2[1] - P1[1];
				P1P2[2] = P2[2] - P1[2]; 

				s = dotproduct(P1P2, P1P2); 
				dp = dotproduct(PP1, P1P2); 
			
				n[0] = -PP1[0] + dp * P1P2[0] / s; 
				n[1] = -PP1[1] + dp * P1P2[1] / s; 
				n[2] = -PP1[2] + dp * P1P2[2] / s;	
			
				s = norm(n);  
				normals2[3 * i] = n[0] / s; 
				normals2[3 * i + 1] = n[1] / s; 
				normals2[3 * i + 2] = n[2] / s;  
			}

			// normals
			crossproduct(n, P2P3, P1P3); 
			s = norm(n); 
			if (s == 0) {// colinear neighbors 		
				//n=p1p3 n (p1p3 n pp1)
				crossproduct(v, P1P3, PP1); 
				crossproduct(n, P1P3, v);
				s = norm(n);
				if (s == 0) 
					SetNormal(i, 1, 0, 0); // colinear neighbors and point
				else 
					SetNormal(i, n[0] / s, n[1] / s, n[2] / s);
			}
			else 
				SetNormal(i, n[0] / s, n[1] / s, n[2] / s);	
		}

		// Compute averaged normals for multimaterial points having 2 sets of neighbors.
		else if (IsMultiMaterialPoint(i) == true) {
			double normals_mm[3] = {0, 0, 0}, normals2_mm[3] = {0, 0, 0};

			int *mat = GetPointMaterialIndices(i);
			for (int iter = 0; iter < 2; iter++) {

				int *pointNeighbors = GetMultiMaterialNeighbors(mat[iter], i);
				P1 = GetPoint(pointNeighbors[0]); 
				P2 = GetPoint(pointNeighbors[1]);  
				P3 = GetPoint(pointNeighbors[2]); 

				P1P3[0] = P3[0] - P1[0];
				P1P3[1] = P3[1] - P1[1];
				P1P3[2] = P3[2] - P1[2];

				P2P3[0] = P3[0] - P2[0];
				P2P3[1] = P3[1] - P2[1];
				P2P3[2] = P3[2] - P2[2];

				PP1[0] = P1[0] - P[0];
				PP1[1] = P1[1] - P[1];
				PP1[2] = P1[2] - P[2];

				// normals2
				if (neighbors_c[3 * i + 2] == -1) {
					P1P2[0] = P2[0] - P1[0];
					P1P2[1] = P2[1] - P1[1];
					P1P2[2] = P2[2] - P1[2]; 

					s = dotproduct(P1P2, P1P2); 
					dp = dotproduct(PP1, P1P2); 
			
					n[0] = -PP1[0] + dp * P1P2[0] / s; 
					n[1] = -PP1[1] + dp * P1P2[1] / s; 
					n[2] = -PP1[2] + dp * P1P2[2] / s;	
			
					s = norm(n);  
					normals2_mm[0] = normals2_mm[0] + n[0] / s; 
					normals2_mm[1] = normals2_mm[1] + n[1] / s; 
					normals2_mm[2] = normals2_mm[2] + n[2] / s;  
				}

				// normals
				crossproduct(n, P2P3, P1P3); 
				s = norm(n); 
				if (s == 0) {// colinear neighbors 		
					//n=p1p3 n (p1p3 n pp1)
					crossproduct(v, P1P3, PP1); 
					crossproduct(n, P1P3, v);
					s = norm(n);
					if (s == 0) {
						//SetNormal(i, 1, 0, 0); // colinear neighbors and point
						normals_mm[0] = normals_mm[0] + 1;
						normals_mm[1] = normals_mm[1] + 0;
						normals_mm[2] = normals_mm[2] + 0;
					}
					else {
						//SetNormal(i, n[0] / s, n[1] / s, n[2] / s);
						normals_mm[0] = normals_mm[0] + n[0] / s; 
						normals_mm[1] = normals_mm[1] + n[1] / s; 
						normals_mm[2] = normals_mm[2] + n[2] / s;  
					}
				}
				else {
					//SetNormal(i, n[0] / s, n[1] / s, n[2] / s);
					normals_mm[0] = normals_mm[0] + n[0] / s; 
					normals_mm[1] = normals_mm[1] + n[1] / s; 
					normals_mm[2] = normals_mm[2] + n[2] / s;  

					
				}
			}


			// Average the normals
			normals2_mm[0] = normals2_mm[0] / 2;
			normals2_mm[1] = normals2_mm[1] / 2;
			normals2_mm[2] = normals2_mm[2] / 2;

			normals_mm[0] = normals_mm[0] / 2;
			normals_mm[1] = normals_mm[1] / 2;
			normals_mm[2] = normals_mm[2] / 2;

			/*normals2[3 * i + 0] = normals2_mm[0];
			normals2[3 * i + 1] = normals2_mm[1];
			normals2[3 * i + 2] = normals2_mm[2];*/

			SetNormal(i, normals_mm);
		}
	}
}
void CSimplexSurf::SetNormal(int index,P_float x,P_float y, P_float z) {normals[3*index]=x;normals[3*index+1]=y;normals[3*index+2]=z;}
void CSimplexSurf::SetNormal(int index,P_float p[3]) {memcpy(normals+3*index,p,3*sizeof(P_float));}
void CSimplexSurf::GetNormal(int index,P_float p[3]) {memcpy(p,normals+3*index,3*sizeof(P_float));}
P_float* CSimplexSurf::GetNormal(int index) {return normals+3*index;}

void CSimplexSurf::GetTangentVectors(int index,P_float t1[3],P_float t2[3]) {
	t1[0]=points[3*neighbors[3*index]]-points[3*neighbors[3*index+1]]; 
	t1[1]=points[3*neighbors[3*index]+1]-points[3*neighbors[3*index+1]+1]; 
	t1[2]=points[3*neighbors[3*index]+2]-points[3*neighbors[3*index+1]+2];
	
	P_float nt1=norm(t1); 
	t1[0]=t1[0]/nt1; 
	t1[1]=t1[1]/nt1; 
	t1[2]=t1[2]/nt1;
	crossproduct(t2,normals+3*index,t1);
}
void CSimplexSurf::GetTangentVectors(P_float *t1,P_float *t2) {for(int i=0;i<nb_points;i++) GetTangentVectors(i,t1+3*i,t2+3*i);}
/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::ComputeVolume() {
	int i, j, i2, i3;
	P_float cp[3], s, *p1, *p2, *p3, p1p3[3], p1p2[3];
	P_float volume2 = 0;
	surface = 0;

	for (i = 0; i < nb_points; i++) surfaces[i] = 0;

	for (i = 0; i < nb_cells; i++) { 
		surfaces_cell[i] = 0;
		p1 = cellcenters + 3 * i;
		// for each face calcul volume and surface
		for (j = 0; j < cells[i][0]; j++) {
			if (Flipped) {
				i3 = cells[i][j + 1]; 
				i2 = cells[i][(j == cells[i][0] - 1) ? 1 : (j + 2)];
			}
			else {
				i2 = cells[i][j + 1]; 
				i3 = cells[i][(j == cells[i][0] - 1) ? 1 : (j + 2)];
			}
			
			p2 = GetPoint(i2); 
			p3 = GetPoint(i3);

			p1p3[0] = p3[0] - p1[0]; 
			p1p3[1] = p3[1] - p1[1]; 
			p1p3[2] = p3[2] - p1[2];  
			
			p1p2[0] = p2[0] - p1[0]; 
			p1p2[1] = p2[1] - p1[1]; 
			p1p2[2] = p2[2] - p1[2];

			crossproduct(cp, p1p2, p1p3); 
			s = norm(cp) / 2.;
			
			volume2 += cp[2] * (p1[2] + p2[2] + p3[2]) / 6.;
			
			surface += s; 
			surfaces_cell[i] += s;
			
			surfaces[i2] += s / 2.; 
			surfaces[i3] += s / 2.; // 1 triangle <-> 2 points
		}
	}
	
	dv = volume2 - volume; 
	volume = volume2;	
}

void CSimplexSurf::SetVolume(P_float vol) {volume=vol;}
P_float CSimplexSurf::GetVolume() {return volume;}
P_float CSimplexSurf::GetVolume_ref() {return volume_ref;}
void CSimplexSurf::SetVolume_ref(P_float val,int mode) 
{
if(mode==0) volume_ref=val; 
else if(mode==1) volume_ref=RefModelVolume*val; 
else if(mode==2) volume_ref=volume*val;
}

void CSimplexSurf::SetSurface_ref(P_float val,int mode)
{
if(mode==0) surface_ref=val; 
else if(mode==1) surface_ref=RefModelSurface*val; 
else if(mode==2) surface_ref=surface*val;
}


P_float  CSimplexSurf::GetDv() {return dv;}

void CSimplexSurf::UpdateCellCenters()
{
int i,j;
P_float *p;

for(i=0;i<nb_cells;i++) 
    { 
    cellcenters[3*i]=0; cellcenters[3*i+1]=0; cellcenters[3*i+2]=0; 
    for(j=0;j<cells[i][0];j++) 
        {
        p=GetPoint(cells[i][j+1]);
        cellcenters[3*i]+=p[0]; cellcenters[3*i+1]+=p[1]; cellcenters[3*i+2]+=p[2];
        }
    cellcenters[3*i]=cellcenters[3*i]/(P_float)cells[i][0]; cellcenters[3*i+1]=cellcenters[3*i+1]/(P_float)cells[i][0]; cellcenters[3*i+2]=cellcenters[3*i+2]/(P_float)cells[i][0];
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::SetMRI_GradientMRI(vtkStructuredPoints* mri) {
	if(mri == NULL) return;

	MRI_GradientMRI = mri;
	MRI_ROI_GradientMRI[0] = mri->GetOrigin()[0]; 
	MRI_ROI_GradientMRI[1] = mri->GetOrigin()[0] + (P_float)(mri->GetDimensions()[0] - 1) * mri->GetSpacing()[0];
	MRI_ROI_GradientMRI[2] = mri->GetOrigin()[1]; 
	MRI_ROI_GradientMRI[3] = mri->GetOrigin()[1] + (P_float)(mri->GetDimensions()[1] - 1) * mri->GetSpacing()[1];
	MRI_ROI_GradientMRI[4] = mri->GetOrigin()[2]; 
	MRI_ROI_GradientMRI[5] = mri->GetOrigin()[2] + (P_float)(mri->GetDimensions()[2] - 1) * mri->GetSpacing()[2];
}

void CSimplexSurf::AddMRI_GradientMRI(vtkStructuredPoints* mri) 
{
if(mri==NULL) return;
MRI_GradientMRI2=mri;
MRI_ROI_GradientMRI2[0]=mri->GetOrigin()[0]; MRI_ROI_GradientMRI2[1]=mri->GetOrigin()[0]+(P_float)(mri->GetDimensions()[0]-1)*mri->GetSpacing()[0];
MRI_ROI_GradientMRI2[2]=mri->GetOrigin()[1]; MRI_ROI_GradientMRI2[3]=mri->GetOrigin()[1]+(P_float)(mri->GetDimensions()[1]-1)*mri->GetSpacing()[1];
MRI_ROI_GradientMRI2[4]=mri->GetOrigin()[2]; MRI_ROI_GradientMRI2[5]=mri->GetOrigin()[2]+(P_float)(mri->GetDimensions()[2]-1)*mri->GetSpacing()[2];
}


void CSimplexSurf::FreeMRI(int MRI_index) 
{
if(MRI_Type[MRI_index]==2) {free(MRI_Transform[MRI_index]); free(MRI_Spacing[MRI_index]);}
if(MRI_Normals[MRI_index]!=NULL) free(MRI_Normals[MRI_index]);
if(MRI_OffsetTransform[MRI_index]!=NULL) free(MRI_OffsetTransform[MRI_index]);
if(MRI_ROI[MRI_index]!=NULL) free(MRI_ROI[MRI_index]);
}

void CSimplexSurf::SetMRI(vtkStructuredPoints* mri,vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]) {
	if (MRI_nb == 0) 
		AddMRI();
	
	FreeMRI(MRI_index);
	MRI_Type[MRI_index] = 0; 
	MRI[MRI_index] = mri; 
	MRI_grad[MRI_index] = mri_grad;
	SetMRI_ROI(MRI_index, bounds);
}

void CSimplexSurf::SetMRI(vtkStructuredPoints* mri,P_float offset[16],vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]) 
{
if(MRI_nb==0) AddMRI();
FreeMRI(MRI_index);
MRI_Type[MRI_index]=0; MRI[MRI_index]=mri; MRI_grad[MRI_index]=mri_grad;
MRI_OffsetTransform[MRI_index]=new P_float[16]; memcpy(MRI_OffsetTransform[MRI_index],offset,16*sizeof(P_float));
SetMRI_ROI(MRI_index,bounds);
}

void CSimplexSurf::SetMRI(vtkStructuredPoints* mri,P_float **transform_RTMRI,P_float **spacing_RTMRI,vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]) 
{
if(MRI_nb==0) AddMRI();
FreeMRI(MRI_index);
MRI_Type[MRI_index]=1; MRI[MRI_index]=mri; MRI_grad[MRI_index]=mri_grad; 
MRI_Transform[MRI_index]=transform_RTMRI; MRI_Spacing[MRI_index]=spacing_RTMRI; 
SetMRI_ROI(MRI_index,bounds);
}

void CSimplexSurf::SetMRI(vtkStructuredPoints* mri,P_float *transform_RadialMRI,P_float *spacing_RadialMRI,vtkStructuredPoints* mri_grad,int MRI_index,P_float bounds[3][2]) 
{
if(MRI_nb==0) AddMRI();
FreeMRI(MRI_index);
MRI_Type[MRI_index]=2; MRI[MRI_index]=mri; MRI_grad[MRI_index]=mri_grad; 
MRI_Transform[MRI_index]=new P_float*[1]; MRI_Spacing[MRI_index]=new P_float*[1]; 
MRI_Transform[MRI_index][0]=transform_RadialMRI; MRI_Spacing[MRI_index][0]=spacing_RadialMRI; 
SetMRI_ROI(MRI_index,bounds);
}


void CSimplexSurf::AddMRI() {
	if (MRI_nb == 0) {
		MRI = new vtkStructuredPoints*[1]; 
		MRI_grad = new vtkStructuredPoints*[1]; 
		MRI_Transform = new P_float **[1]; 
		MRI_Spacing = new P_float **[1]; 
		MRI_Normals = new P_float*[1]; 
		MRI_OffsetTransform = new P_float*[1];	
		MRI_ROI = new P_float*[1]; 
		MRI_Type = new int[1]; 
	}
	else {
		MRI = (vtkStructuredPoints**)realloc(MRI, (MRI_nb + 1) * sizeof(vtkStructuredPoints*)); 
		MRI_grad = (vtkStructuredPoints**)realloc(MRI_grad, (MRI_nb + 1) * sizeof(vtkStructuredPoints*)); 
		MRI_Transform = (P_float ***)realloc(MRI_Transform, (MRI_nb + 1) * sizeof(P_float **)); 
		MRI_Spacing = (P_float ***)realloc(MRI_Spacing, (MRI_nb + 1) * sizeof(P_float **)); 
		MRI_Normals = (P_float**)realloc(MRI_Normals, (MRI_nb + 1) * sizeof(P_float*)); 
		MRI_OffsetTransform = (P_float**)realloc(MRI_OffsetTransform, (MRI_nb + 1) * sizeof(P_float*));	
		MRI_ROI = (P_float**)realloc(MRI_ROI, (MRI_nb + 1) * sizeof(P_float*)); 
		MRI_Type = (int*)realloc(MRI_Type, (MRI_nb + 1) * sizeof(int)); 
	}
	
	MRI[MRI_nb] = NULL; 
	MRI_grad[MRI_nb] = NULL; 
	MRI_Transform[MRI_nb] = NULL; 
	MRI_Spacing[MRI_nb] = NULL; 
	MRI_Normals[MRI_nb] = NULL; 
	MRI_OffsetTransform[MRI_nb] = NULL;	
	MRI_ROI[MRI_nb] = NULL; 
	MRI_Type[MRI_nb] = -1;
	MRI_nb++;
}

void CSimplexSurf::AddMRI(vtkStructuredPoints* mri,vtkStructuredPoints* mri_grad,P_float bounds[3][2]) 
{
AddMRI();
MRI_Type[MRI_nb-1]=0; MRI[MRI_nb-1]=mri; MRI_grad[MRI_nb-1]=mri_grad;
SetMRI_ROI(MRI_nb-1,bounds);
}

void CSimplexSurf::AddMRI(vtkStructuredPoints* mri,P_float offset[16],vtkStructuredPoints* mri_grad,P_float bounds[3][2]) 
{
AddMRI();
MRI_Type[MRI_nb-1]=0; MRI[MRI_nb-1]=mri; MRI_grad[MRI_nb-1]=mri_grad;
MRI_OffsetTransform[MRI_nb-1]=new P_float[16]; memcpy(MRI_OffsetTransform[MRI_nb-1],offset,16*sizeof(P_float));
SetMRI_ROI(MRI_nb-1,bounds);
}

void CSimplexSurf::AddMRI(vtkStructuredPoints* mri,P_float **transform_RTMRI,P_float **spacing_RTMRI,vtkStructuredPoints* mri_grad,P_float bounds[3][2]) 
{
AddMRI();
MRI_Type[MRI_nb-1]=1; MRI[MRI_nb-1]=mri; MRI_grad[MRI_nb-1]=mri_grad; 
MRI_Transform[MRI_nb-1]=transform_RTMRI; MRI_Spacing[MRI_nb-1]=spacing_RTMRI; 
SetMRI_ROI(MRI_nb-1,bounds);
}

void CSimplexSurf::AddMRI(vtkStructuredPoints* mri,P_float *transform_RadialMRI,P_float *spacing_RadialMRI,vtkStructuredPoints* mri_grad,P_float bounds[3][2]) 
{
AddMRI();
MRI_Type[MRI_nb-1]=2; MRI[MRI_nb-1]=mri; MRI_grad[MRI_nb-1]=mri_grad; 
MRI_Transform[MRI_nb-1]=new P_float*[1]; MRI_Spacing[MRI_nb-1]=new P_float*[1]; 
MRI_Transform[MRI_nb-1][0]=transform_RadialMRI; MRI_Spacing[MRI_nb-1][0]=spacing_RadialMRI; 
SetMRI_ROI(MRI_nb-1,bounds);
}

void CSimplexSurf::SetMRI_ROI(int MRI_index,P_float bounds[3][2])
{
MRI_ROI[MRI_index]=new P_float[6];  
MRI_ROI[MRI_index][0]=bounds[0][0]; MRI_ROI[MRI_index][1]=bounds[0][1];
MRI_ROI[MRI_index][2]=bounds[1][0]; MRI_ROI[MRI_index][3]=bounds[1][1];
MRI_ROI[MRI_index][4]=bounds[2][0]; MRI_ROI[MRI_index][5]=bounds[2][1];
}

//void CSimplexSurf::SetMRI2(vtkStructuredPoints* mri,P_float transform[16],bool* flagpt,vtkStructuredPoints* mri_grad) 
//{
//MRI2_grad=mri_grad; 
//MRI2=mri;
//memcpy(MRI2_transform,transform,16*sizeof(P_float));
//MRI_flagpt=flagpt;
//}

void CSimplexSurf::Analyse_IntensityProfile(const char* filename)
{
P_float tresh=1E10; //if (Metric_IntensityProfile==0 || Metric_IntensityProfile==2) tresh=1E10; else tresh=-.2;

int i;
FILE* f=fopen(filename,"rt+"); if(f==NULL) f=fopen(filename,"wt"); else fseek(f,0,SEEK_END);
P_float ACC_DO_NOM_RON_CR[5]; 

/*
if(ExternalForce_Demon) 
	{
	for(i=0;i<3*nb_points;i++) externalforces[i]=0; for(i=0;i<nb_points;i++) externalforces_ign[i]=false;
	P_float *nml=normals; if(Transform_RTMRI!=NULL) nml=Normals_RTMRI;
	GetDemonForces(externalforces,externalforces_ign,nb_points,nml,IntensityProfile,IntensityProfileRef[0],IntensityProfile_s,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref,1);
	ACC_DO_NOM_RON_CR[0]=0; for(i=0;i<nb_points;i++) ACC_DO_NOM_RON_CR[0]+=dotproduct(externalforces+3*i,externalforces+3*i); ACC_DO_NOM_RON_CR[0]=sqrt(ACC_DO_NOM_RON_CR[0]/(P_float)nb_points);
	}
else*/

bool grad=false; if(MRI_grad[0]!=NULL) grad=true;
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;

ComputeIntensityProfile(IntensityProfile_mn_ref+Depth_IntensityProfile,IntensityProfile_pn_ref+Depth_IntensityProfile,IntensityProfile_s,grad); 

P_float *E=new P_float[nb_points*(2*Depth_IntensityProfile+1)]; for(i=0;i<nb_points*(2*Depth_IntensityProfile+1);i++) E[i]=1E10;
if(ExternalForce_GradientMRI) {GetGradEnergies(E,nb_points,points,normals,MRI_GradientMRI,Size_GradientMRI,depth_GradientMRI,Interpolationmode,MRI_ROI_GradientMRI,oppositedirection_GradientMRI); GetGradEnergies(E,nb_points,points,normals,MRI_GradientMRI2,Size_GradientMRI,depth_GradientMRI,Interpolationmode,MRI_ROI_GradientMRI2,oppositedirection_GradientMRI);}
else if(vec) GetIPEnergies(E,Metric_IntensityProfile-2,nb_points,IntensityProfile_grad,IntensityProfileRef_grad,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
else if(grad) GetIPEnergies(E,Metric_IntensityProfile-2,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
else GetIPEnergies(E,Metric_IntensityProfile,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);

AnalyseEnergy(ACC_DO_NOM_RON_CR,E,nb_points,Depth_IntensityProfile,IntensityProfile_s,tresh);

if(ExternalForce_Demon && grad && !vec && IntensityProfile_Proc_Mask!=NULL) fprintf(f,"DemonGRADproc: \tACC %4.2f Size %4.2f Depth %d\n",ACC_DO_NOM_RON_CR[0],IntensityProfile_s,Depth_IntensityProfile);
else if(ExternalForce_Demon && grad && !vec) fprintf(f,"DemonGRAD: \tACC %4.2f Size %4.2f Depth %d\n",ACC_DO_NOM_RON_CR[0],IntensityProfile_s,Depth_IntensityProfile);
else if(ExternalForce_Demon && IntensityProfile_Proc_Mask!=NULL) fprintf(f,"DemonProc: \tACC %4.2f Size %4.2f Depth %d\n",ACC_DO_NOM_RON_CR[0],IntensityProfile_s,Depth_IntensityProfile);
else if(ExternalForce_Demon) fprintf(f,"Demon: \tACC %4.2f Size %4.2f Depth %d\n",ACC_DO_NOM_RON_CR[0],IntensityProfile_s,Depth_IntensityProfile);
else if(ExternalForce_GradientMRI && MRI_GradientMRI->GetNumberOfScalarComponents()==3) fprintf(f,"GRADV: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],Size_GradientMRI,depth_GradientMRI);
else if(ExternalForce_GradientMRI && MRI_GradientMRI->GetNumberOfScalarComponents()==1) fprintf(f,"GRAD: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],Size_GradientMRI,depth_GradientMRI);
else if(Metric_IntensityProfile==2 && vec) fprintf(f,"ADGRADV: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==3 && vec) fprintf(f,"NCCGRADV: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==2 && IntensityProfile_Proc_Mask!=NULL) fprintf(f,"ADGRADproc: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==3 && IntensityProfile_Proc_Mask!=NULL) fprintf(f,"NCCGRADproc: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==2) fprintf(f,"ADGRAD: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==3) fprintf(f,"NCCGRAD: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==0 && IntensityProfile_Proc_Mask!=NULL) fprintf(f,"ADproc: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==1 && IntensityProfile_Proc_Mask!=NULL) fprintf(f,"NCCproc: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==0) fprintf(f,"AD: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
else if(Metric_IntensityProfile==1) fprintf(f,"NCC: \tACC %4.2f DO %4.2f NOM %4.2f RON %4.2f CR %4.2f Size %4.2f Depth %d mn %d pn %d\n",ACC_DO_NOM_RON_CR[0],ACC_DO_NOM_RON_CR[1],ACC_DO_NOM_RON_CR[2],ACC_DO_NOM_RON_CR[3],ACC_DO_NOM_RON_CR[4],IntensityProfile_s,Depth_IntensityProfile,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
fclose(f);
}



void CSimplexSurf::Test_IntensityProfile(const char* filename)
{
bool grad=false; if(Metric_IntensityProfile==2 || Metric_IntensityProfile==3) grad=true;
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;

P_float tresh=1E10; //if (Metric_IntensityProfile==0 || Metric_IntensityProfile==2) tresh=1E10; else tresh=-.2;

int i,k;
FILE* f=fopen(filename,"wt"); 

// save IPref
int** prof_ref=NULL; P_float** vprof_ref=NULL;
if(vec) {vprof_ref=new P_float*[nb_points]; for(i=0;i<nb_points;i++) {vprof_ref[i]=new P_float[3*(IntensityProfile_mn_ref+IntensityProfile_pn_ref+1)]; memcpy(vprof_ref[i],IntensityProfileRef_grad[i],3*(IntensityProfile_mn_ref+IntensityProfile_pn_ref+1)*sizeof(P_float));}}
else {prof_ref=new int*[nb_points]; for(i=0;i<nb_points;i++) {prof_ref[i]=new int[IntensityProfile_mn_ref+IntensityProfile_pn_ref+1]; memcpy(prof_ref[i],IntensityProfileRef[i],(IntensityProfile_mn_ref+IntensityProfile_pn_ref+1)*sizeof(int));}}

P_float stddev,E,Emin=1E10,stddevmin=1E10;
int mn_min,pn_min,mn,pn,mn_ref=IntensityProfile_mn_ref,pn_ref=IntensityProfile_pn_ref;
P_float *nml=normals; if(MRI_Type[0]==2) nml=MRI_Normals[0];

for(mn=1;mn<=mn_ref;mn++)
    {   
    for(pn=1;pn<=pn_ref;pn++)
        {
		// restore IPref
		if(vec) {for(i=0;i<nb_points;i++) {free(IntensityProfileRef_grad[i]); IntensityProfileRef_grad[i]=new P_float[3*(mn+pn+1)]; memcpy(IntensityProfileRef_grad[i],vprof_ref[i]+3*(mn_ref-mn),3*(mn+pn+1)*sizeof(P_float));}      IntensityProfile_mn_ref=mn; IntensityProfile_pn_ref=pn;}
		else {for(i=0;i<nb_points;i++) {free(IntensityProfileRef[i]); IntensityProfileRef[i]=new int[mn+pn+1]; memcpy(IntensityProfileRef[i],prof_ref[i]+mn_ref-mn,(mn+pn+1)*sizeof(int));}      IntensityProfile_mn_ref=mn; IntensityProfile_pn_ref=pn;}

		//compute ip
		ComputeIntensityProfile(IntensityProfile_mn_ref+Depth_IntensityProfile,IntensityProfile_pn_ref+Depth_IntensityProfile,IntensityProfile_s,grad); 

		// compute forces
		for(k=0;k<3*nb_points;k++) externalforces[k]=0; for(k=0;k<nb_points;k++) externalforces_ign[k]=true;
		if(ExternalForce_Demon) GetDemonForces(externalforces,externalforces_ign,nb_points,nml,IntensityProfile,IntensityProfileRef,IntensityProfile_s,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref,1);
		else
			{
			P_float *pE=new P_float[nb_points*(2*Depth_IntensityProfile+1)]; for(i=0;i<nb_points*(2*Depth_IntensityProfile+1);i++) pE[i]=1E10;
			if(vec) GetIPEnergies(pE,Metric_IntensityProfile-2,nb_points,IntensityProfile_grad,IntensityProfileRef_grad,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
			else if(grad) GetIPEnergies(pE,Metric_IntensityProfile-2,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
			else GetIPEnergies(pE,Metric_IntensityProfile,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);
			GetIPForces(externalforces,externalforces_ign,pE,nb_points,IntensityProfile_mn,IntensityProfile_mn_ref,IntensityProfile_s,nml,1,tresh);
			}
		// compute norm/ stddev
		stddev=0; E=0; for(i=0;i<nb_points;i++) {E+=norm(externalforces+3*i); stddev+=dotproduct(externalforces+3*i,externalforces+3*i);} E=E/(P_float)nb_points; stddev=sqrt(stddev/(P_float)nb_points-E*E);
		fprintf(f,"%lf %d %d %lf %lf\n",IntensityProfile_s,IntensityProfile_mn_ref,IntensityProfile_pn_ref,E,stddev);      if((E+stddev)<(Emin+stddevmin)) {Emin=E; mn_min=mn; pn_min=pn; stddevmin=stddev;}
		}
	}

fprintf(f,"\n\n%d %d %lf %lf\n",mn_min,pn_min,Emin,stddevmin);
fclose(f);

// restore IPref
mn=mn_ref; pn=pn_ref;
IntensityProfile_mn_ref=mn; IntensityProfile_pn_ref=pn;

if(vec) {for(i=0;i<nb_points;i++) {free(IntensityProfileRef_grad[i]); IntensityProfileRef_grad[i]=new P_float[3*(mn+pn+1)]; memcpy(IntensityProfileRef_grad[i],vprof_ref[i]+3*(mn_ref-mn),3*(mn+pn+1)*sizeof(int));}}
else {for(i=0;i<nb_points;i++) {free(IntensityProfileRef[i]); IntensityProfileRef[i]=new int[mn+pn+1]; memcpy(IntensityProfileRef[i],prof_ref[i]+mn_ref-mn,(mn+pn+1)*sizeof(int));}}

if(vec) {for(i=0;i<nb_points;i++) free(vprof_ref[i]); free(vprof_ref);}
else {for(i=0;i<nb_points;i++) free(prof_ref[i]); free(prof_ref);}
}


void CSimplexSurf::ComputeSimilarityImage(int depth,const char* filename)
{
int i,j;
P_float tresh; if (Metric_IntensityProfile==0 || Metric_IntensityProfile==2) tresh=1E10; else tresh=-.2;
bool grad=false; if(Metric_IntensityProfile==2 || Metric_IntensityProfile==3) grad=true;
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;

ComputeIntensityProfile(IntensityProfile_mn_ref+depth,IntensityProfile_pn_ref+depth,IntensityProfile_s,grad); 

P_float *E=new P_float[nb_points*(2*Depth_IntensityProfile+1)]; for(i=0;i<nb_points*(2*Depth_IntensityProfile+1);i++) E[i]=1E10;
if(vec) GetIPEnergies(E,Metric_IntensityProfile-2,nb_points,IntensityProfile_grad,IntensityProfileRef_grad,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
else if(grad) GetIPEnergies(E,Metric_IntensityProfile-2,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
else GetIPEnergies(E,Metric_IntensityProfile,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);

vtkStructuredPoints* im=vtkStructuredPoints::New();
im->SetDimensions(2*depth+1,nb_points,1);
im->SetScalarTypeToUnsignedChar();
im->AllocateScalars();
unsigned char* pt=(unsigned char*)im->GetScalarPointer();
P_float min,max;

for(j=0;j<nb_points;j++)
    {
	// line normalization
	min=1E20;max=-1E10;
	for(i=0;i<2*depth+1;i++) if(E[i+j*(2*depth+1)]!=1E10) {if(E[i+j*(2*depth+1)]>max) max=E[i+j*(2*depth+1)]; if(E[i+j*(2*depth+1)]<min) min=E[i+j*(2*depth+1)]; }
	for(i=0;i<2*depth+1;i++) if(E[i+j*(2*depth+1)]!=1E10) pt[i+j*(2*depth+1)]=(unsigned char)floor( ((P_float)255*(E[i+j*(2*depth+1)]-min)) / (max-min) ); else pt[i+j*(2*depth+1)]=255;
	}
free(E);

std::string str=filename;
if(str.find(".vtk")!=-1)
    {
    vtkStructuredPointsWriter* writer=vtkStructuredPointsWriter::New();
    writer->SetInput(im);
    writer->SetFileTypeToBinary();
    writer->SetFileName(filename);
    writer->Write();
    writer->Delete();
    }
else
    {
    vtkJPEGWriter* writer2=vtkJPEGWriter::New();
    writer2->SetInput(im);
    writer2->SetFileDimensionality(2);
    writer2->SetFileName(filename);
    writer2->Write();
    writer2->Delete();
    }
}

void CSimplexSurf::ComputeIntensityProfile_GaussianWeights(P_float r,P_float stddev,P_float s)
{
int x,y;
int R=1+(int)ceil(r/s);
P_float l2,stddev2=stddev*stddev,sum=0;

if(IntensityProfile_GaussianWeights!=NULL) {for(y=0;y<IntensityProfile_GaussianR;y++) free(IntensityProfile_GaussianWeights[y]); free(IntensityProfile_GaussianWeights);}
IntensityProfile_GaussianWeights=new P_float*[R]; for(y=0;y<R;y++) IntensityProfile_GaussianWeights[y]=new P_float[R];
IntensityProfile_GaussianS=s; IntensityProfile_GaussianR=R;

for(y=0;y<R;y++)
	for(x=0;x<R;x++)
		{
		l2=s*s*(x*x+y*y);
		IntensityProfile_GaussianWeights[y][x]=exp(-l2/stddev2);
		}

sum=IntensityProfile_GaussianWeights[0][0];
for(y=1;y<R;y++) sum+=2*IntensityProfile_GaussianWeights[y][0];
for(x=1;x<R;x++) sum+=2*IntensityProfile_GaussianWeights[0][x];
for(y=1;y<R;y++) for(x=1;x<R;x++) sum+=4*IntensityProfile_GaussianWeights[y][x];
for(y=0;y<R;y++) for(x=0;x<R;x++) IntensityProfile_GaussianWeights[y][x]=IntensityProfile_GaussianWeights[y][x]/sum;
}


void CSimplexSurf::IniIntensityProfile(int mn,int pn,P_float s,bool grad)
{
int i,j;

bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;
if(!vec) {if(IntensityProfile==NULL) {IntensityProfile=(int**)malloc(nb_points*sizeof(int*));				IntensityProfile_s=s; IntensityProfile_pn=pn; IntensityProfile_mn=mn; for(j=0;j<nb_points;j++)  IntensityProfile[j]=new int[IntensityProfile_pn+IntensityProfile_mn+1]; }}
else {if(IntensityProfile_grad==NULL) {IntensityProfile_grad=(P_float**)malloc(nb_points*sizeof(P_float*)); IntensityProfile_s=s; IntensityProfile_pn=pn; IntensityProfile_mn=mn; for(j=0;j<nb_points;j++)  IntensityProfile_grad[j]=new P_float[3*(IntensityProfile_pn+IntensityProfile_mn+1)]; }}

if(IntensityProfile_mn!=mn || IntensityProfile_pn!=pn || IntensityProfile_s!=s)
    {
    IntensityProfile_s=s; IntensityProfile_pn=pn; IntensityProfile_mn=mn;
    if(!vec) for(j=0;j<nb_points;j++)   {free(IntensityProfile[j]);		 IntensityProfile[j]=new int[IntensityProfile_pn+IntensityProfile_mn+1];}
	else for(j=0;j<nb_points;j++)		{free(IntensityProfile_grad[j]); IntensityProfile_grad[j]=new P_float[3*(IntensityProfile_pn+IntensityProfile_mn+1)];}
    }

if(!vec) for(i=0;i<nb_points;i++) for(j=0;j<mn+pn+1;j++)  IntensityProfile[i][j]=-1E5;
else for(i=0;i<nb_points;i++) for(j=0;j<3*(mn+pn+1);j++)  IntensityProfile_grad[i][j]=-1E5;
}

void CSimplexSurf::ComputeIntensityProfile(int mn,int pn,P_float s,bool grad)
{
IniIntensityProfile(mn,pn,s,grad);
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;
for(int j=MRI_nb-1;j>=0;j--)
	{
	if(MRI_Normals[j]!=NULL) {free(MRI_Normals[j]); MRI_Normals[j]=NULL;}

	if(MRI_Type[j]==2) //radial
		{
		ComputeIP(nb_points,points,normals,MRI[j],MRI_Transform[j][0],MRI_Spacing[j][0],IntensityProfile,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s,Interpolationmode,MRI_ROI[j]);
		}
	else if(MRI_Type[j]==1) //rt
		{
		MRI_Normals[j]=ComputeIP(nb_points,points,normals,MRI[j],MRI_Transform[j],MRI_Spacing[j],10,IntensityProfile,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s,Interpolationmode,MRI_ROI[j]);
		}
	else if(IntensityProfile_GaussianWeights!=NULL) 
		{
		P_float *tangent1=new P_float[3*nb_points],*tangent2=new P_float[3*nb_points]; GetTangentVectors(tangent1,tangent2);
		ComputeIP(nb_points,points,normals,tangent1,tangent2,MRI[j],IntensityProfile,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s,IntensityProfile_GaussianWeights,IntensityProfile_GaussianR,Interpolationmode,MRI_ROI[j],MRI_OffsetTransform[j]);
		free(tangent1); free(tangent2);
		}
	else //cartesian
		{
		if(!grad) ComputeIP(nb_points,points,normals,MRI[j],IntensityProfile,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s,Interpolationmode,MRI_ROI[j],MRI_OffsetTransform[j]);
		else if(MRI_grad[j]->GetNumberOfScalarComponents()==1) 		ComputeIP(nb_points,points,normals,MRI_grad[j],IntensityProfile,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s,Interpolationmode,MRI_ROI[j],MRI_OffsetTransform[j]);
		else ComputeIP(nb_points,points,normals,MRI_grad[j],IntensityProfile_grad,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s,Interpolationmode,MRI_ROI[j],MRI_OffsetTransform[j]);
		}
	}
if(!vec && IntensityProfile_Proc_Mask!=NULL) ProcessIntensityProfile(); 
}



void CSimplexSurf::SaveIntensityProfile(const char* filename,bool grad)
{
int j; 
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;
std::string str=filename;

if(str.find(".vtk")!=-1)
    {
	vtkStructuredPoints* im=vtkStructuredPoints::New();
	im->SetDimensions(IntensityProfile_mn+IntensityProfile_pn+1,nb_points,1); 
	if(!vec) im->SetScalarTypeToUnsignedShort(); else im->SetScalarTypeToFloat();
	if(!vec) im->SetNumberOfScalarComponents(1); else im->SetNumberOfScalarComponents(3);
	im->AllocateScalars();
	unsigned short* pt; float* ptf;	if(!vec) pt=(unsigned short*)im->GetScalarPointer(); else ptf=(float*)im->GetScalarPointer();

	if(!vec) {for(j=0;j<nb_points;j++)   for(int i=0;i<IntensityProfile_mn+IntensityProfile_pn+1;i++) if(IntensityProfile[j][i]!=-1E5) *(pt+i+j*(IntensityProfile_mn+IntensityProfile_pn+1))=(unsigned short)IntensityProfile[j][i]; else *(pt+i+j*(IntensityProfile_mn+IntensityProfile_pn+1))=0;}
	else {for(j=0;j<nb_points;j++)   for(int i=0;i<IntensityProfile_mn+IntensityProfile_pn+1;i++) if(IntensityProfile_grad[j][3*i]!=-1E5) { *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1)))=(float)IntensityProfile_grad[j][3*i]; *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+1)=(float)IntensityProfile_grad[j][3*i+1]; *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+2)=(float)IntensityProfile_grad[j][3*i+2]; } else {*(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1)))=0;  *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+1)=0;  *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+2)=0; }}

	vtkStructuredPointsWriter* writer=vtkStructuredPointsWriter::New();
    writer->SetInput(im);
    writer->SetFileTypeToBinary();
    writer->SetFileName(filename);
    writer->Write();
    writer->Delete();
	im->Delete();
	}
else if(str.find(".jpg")!=-1)
    {
	vtkStructuredPoints* im=vtkStructuredPoints::New();
	im->SetDimensions(IntensityProfile_mn+IntensityProfile_pn+1,nb_points,1); 
	if(!vec) im->SetScalarTypeToUnsignedShort(); else im->SetScalarTypeToFloat();
	if(!vec) im->SetNumberOfScalarComponents(1); else im->SetNumberOfScalarComponents(3);
	im->AllocateScalars();
	unsigned short* pt; float* ptf;	if(!vec) pt=(unsigned short*)im->GetScalarPointer(); else ptf=(float*)im->GetScalarPointer();

	if(!vec) {for(j=0;j<nb_points;j++)   for(int i=0;i<IntensityProfile_mn+IntensityProfile_pn+1;i++) if(IntensityProfile[j][i]!=-1E5) *(pt+i+j*(IntensityProfile_mn+IntensityProfile_pn+1))=(unsigned short)IntensityProfile[j][i]; else *(pt+i+j*(IntensityProfile_mn+IntensityProfile_pn+1))=0;}
	else {for(j=0;j<nb_points;j++)   for(int i=0;i<IntensityProfile_mn+IntensityProfile_pn+1;i++) if(IntensityProfile_grad[j][3*i]!=-1E5) { *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1)))=(float)IntensityProfile_grad[j][3*i]; *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+1)=(float)IntensityProfile_grad[j][3*i+1]; *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+2)=(float)IntensityProfile_grad[j][3*i+2]; } else {*(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1)))=0;  *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+1)=0;  *(ptf+3*(i+j*(IntensityProfile_mn+IntensityProfile_pn+1))+2)=0; }}

	if(!vec) Normalize(im,255,0);
	vtkImageCast* cast=vtkImageCast::New();
		cast->SetInput(im);
		cast->SetOutputScalarTypeToUnsignedChar();
	vtkJPEGWriter* writer2=vtkJPEGWriter::New();
		writer2->SetInput(cast->GetOutput());
		writer2->SetFileDimensionality(2);
		writer2->SetFileName(filename);
		writer2->Write();
		writer2->Delete();
	cast->Delete();
	im->Delete();
	}
else
	{
	FILE* f=fopen(filename,"wt");
	fprintf(f,"%d %d %lf \n",IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_s);
	for(int i=0;i<nb_points;i++) 
		{
		if(!vec) for(int j=0;j<IntensityProfile_mn+IntensityProfile_pn+1;j++) fprintf(f,"%d ",IntensityProfile[i][j]);
		else for(int j=0;j<IntensityProfile_mn+IntensityProfile_pn+1;j++) {fprintf(f,"%4.2f %4.2f %4.2f ",IntensityProfile_grad[i][3*j],IntensityProfile_grad[i][3*j+1],IntensityProfile_grad[i][3*j+2]);}
		fprintf(f,"\n"); 
		}
	fclose(f);
	}

}


void CSimplexSurf::IniIntensityProfileRef(int mn,int pn,P_float s,bool grad)
{
int i,j;

bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;
if(!vec) {if(IntensityProfileRef==NULL) {IntensityProfileRef=(int**)malloc(nb_points*sizeof(int*));				IntensityProfile_s=s; IntensityProfile_pn_ref=pn; IntensityProfile_mn_ref=mn; for(j=0;j<nb_points;j++)  IntensityProfileRef[j]=new int[IntensityProfile_pn_ref+IntensityProfile_mn_ref+1]; }}
else {if(IntensityProfileRef_grad==NULL) {IntensityProfileRef_grad=(P_float**)malloc(nb_points*sizeof(P_float*)); IntensityProfile_s=s; IntensityProfile_pn_ref=pn; IntensityProfile_mn_ref=mn; for(j=0;j<nb_points;j++)  IntensityProfileRef_grad[j]=new P_float[3*(IntensityProfile_pn_ref+IntensityProfile_mn_ref+1)]; }}

if(IntensityProfile_pn_ref!=pn || IntensityProfile_mn_ref!=mn)
    {
    IntensityProfile_s=s; IntensityProfile_pn_ref=pn; IntensityProfile_mn_ref=mn;
	if(!vec)  for(j=0;j<nb_points;j++)  {free(IntensityProfileRef[j]); IntensityProfileRef[j]=new int[IntensityProfile_pn_ref+IntensityProfile_mn_ref+1]; }
	else	   for(j=0;j<nb_points;j++)  {free(IntensityProfileRef_grad[j]); IntensityProfileRef_grad[j]=new P_float[3*(IntensityProfile_pn_ref+IntensityProfile_mn_ref+1)];}
    }

if(!vec) for(i=0;i<nb_points;i++) for(j=0;j<mn+pn+1;j++)  IntensityProfileRef[i][j]=-1E5;
else for(i=0;i<nb_points;i++) for(j=0;j<3*(mn+pn+1);j++)  IntensityProfileRef_grad[i][j]=-1E5;
}

void CSimplexSurf::LoadIntensityProfileRef(const char* filename,bool grad)
{
int i,j;
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;
FILE* f=fopen(filename,"rt");
int pn,mn; P_float s; fscanf(f,"%d %d %lf \n",&mn,&pn,&s);
IniIntensityProfileRef(mn,pn,s,grad);
if(!vec)   for(i=0;i<nb_points;i++) for(j=0;j<IntensityProfile_mn_ref+IntensityProfile_pn_ref+1;j++)  fscanf(f,"%d ",IntensityProfileRef[i]+j); 
else	   for(i=0;i<nb_points;i++) for(j=0;j<IntensityProfile_mn_ref+IntensityProfile_pn_ref+1;j++)  fscanf(f,"%lf %lf %lf ",IntensityProfileRef_grad[i]+3*j,IntensityProfileRef_grad[i]+3*j+1,IntensityProfileRef_grad[i]+3*j+2); 
fclose(f);
}

void CSimplexSurf::LoadIntensityProfileRefFromHigherRes(CSimplexSurf* mesh,bool grad)
{
int i;
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;
P_float mnpns[3]; mesh->GetIntensityProfileRef(0,grad,mnpns); P_float s=mnpns[2]; int mn=(int)mnpns[0],pn=(int)mnpns[1];
IniIntensityProfileRef(mn,pn,s,grad);
if(!vec)   for(i=0;i<nb_points;i++) memcpy(IntensityProfileRef[i],(int*)mesh->GetIntensityProfileRef(i,grad,mnpns),(mn+pn+1)*sizeof(int));
else	   for(i=0;i<nb_points;i++) memcpy(IntensityProfileRef_grad[i],(P_float*)mesh->GetIntensityProfileRef(i,grad,mnpns),3*(mn+pn+1)*sizeof(P_float));
}

void CSimplexSurf::GetIntensityProfile(int indexpt,int* prof)
{
memcpy(prof,IntensityProfile[indexpt],(IntensityProfile_pn+IntensityProfile_mn+1)*sizeof(int));
}

int* CSimplexSurf::GetIntensityProfile(int indexpt) {return IntensityProfile[indexpt];}
void* CSimplexSurf::GetIntensityProfileRef(int indexpt,bool grad,P_float mnpns[2]) {mnpns[0]=IntensityProfile_mn_ref; mnpns[1]=IntensityProfile_pn_ref; mnpns[2]=IntensityProfile_s; if(grad) {if(MRI_grad[0]->GetNumberOfScalarComponents()) return IntensityProfileRef_grad[indexpt]; else return IntensityProfileRef[indexpt];} else return IntensityProfileRef[indexpt];}

void CSimplexSurf::SetIntensityProfileRef(int pn, int mn,P_float s,int* prof)
{
IniIntensityProfileRef(mn,pn,s,false);
for(int i=0;i<nb_points;i++) memcpy(IntensityProfileRef[i],prof,(IntensityProfile_mn_ref+IntensityProfile_pn_ref+1)*sizeof(int));
}

void CSimplexSurf::SetIntensityProfileProcessingMask(int* mask,int masksize) // convolution
{
if(IntensityProfile_Proc_Mask!=NULL) free(IntensityProfile_Proc_Mask);
IntensityProfile_Proc_Mask=mask;
IntensityProfile_Proc_MaskSize=masksize;
}

void CSimplexSurf::ProcessIntensityProfile() {ProcessIP(nb_points,IntensityProfile,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_Proc_Mask,IntensityProfile_Proc_MaskSize);}
void CSimplexSurf::ProcessIntensityProfileRef() {ProcessIP(nb_points,IntensityProfileRef,IntensityProfile_mn_ref,IntensityProfile_pn_ref,IntensityProfile_Proc_Mask,IntensityProfile_Proc_MaskSize); }



/////////////////////////////////////////////////////////////////////////////////////////////////
//desuet
void CSimplexSurf::ComputeIntensityProfileFromHigherRes(int mn,int pn,P_float s,CSimplexSurf* mesh)
{
int i,j,neighb[3];

IniIntensityProfile(mn,pn,s,false);

int *prof1,*prof2,*prof3,*prof4;

for(i=0;i<nb_points;i++)
    {
    for(j=0;j<IntensityProfile_mn+IntensityProfile_pn+1;j++) IntensityProfile[i][j]=0;
    mesh->GetNeighbors(i,neighb);
    prof1=mesh->GetIntensityProfile(i);  prof2=mesh->GetIntensityProfile(neighb[0]); prof3=mesh->GetIntensityProfile(neighb[1]);  prof4=mesh->GetIntensityProfile(neighb[2]); 
	for(j=0;j<IntensityProfile_mn+IntensityProfile_pn+1;j++) 
		if(prof1[j]==-1E5 || prof2[j]==-1E5 || prof3[j]==-1E5 || prof4[j]==-1E5) IntensityProfile[i][j]=-1E5;
		else IntensityProfile[i][j]=(prof1[j]+prof2[j]+prof3[j]+prof4[j])/4;
    }
}

void CSimplexSurf::ResizeIntensityProfile(int mn,int pn)
{
int* prof=new int[mn+pn+1];

int i,j;
for(j=0;j<nb_points;j++)  
    {
    for(i=-mn;i<=pn;i++)  prof[i+mn]=IntensityProfile[j][i+IntensityProfile_mn];
    free(IntensityProfile[j]); IntensityProfile[j]=new int[mn+pn+1];
    memcpy(IntensityProfile[j],prof,(mn+pn+1)*sizeof(int));
    }

free(prof); 
IntensityProfile_pn=pn;IntensityProfile_mn=mn;
}

void CSimplexSurf::ResizeIntensityProfileRef(int mn,int pn)
{
int* profi=new int[mn+pn+1];

int i,j;
for(j=0;j<nb_points;j++)  
    {
    for(i=-mn;i<=pn;i++)  profi[i+mn]=IntensityProfileRef[j][i+IntensityProfile_mn];
    free(IntensityProfileRef[j]); IntensityProfileRef[j]=new int[mn+pn+1];
    memcpy(IntensityProfileRef[j],profi,(mn+pn+1)*sizeof(int));

    for(i=-mn;i<=pn;i++)  profi[i+mn]=IntensityProfileRef[j][i+IntensityProfile_mn];
    free(IntensityProfileRef[j]); IntensityProfileRef[j]=new int[mn+pn+1];
    memcpy(IntensityProfileRef[j],profi,(mn+pn+1)*sizeof(int));
    }
free(profi);
IntensityProfile_pn_ref=pn;IntensityProfile_mn_ref=mn;
}

void CSimplexSurf::SimplifyIntensityProfile(unsigned short* levels,int nb_levels)
{
int j,nb_cortbone,nb_spongiousbone;
for(int i=0;i<nb_points;i++)
    {
    nb_cortbone=0;
    for(j=0;j<IntensityProfile_pn+IntensityProfile_mn+1;j++) if(IntensityProfile[i][j]<levels[1]) nb_cortbone++;
    if(nb_cortbone==0) nb_cortbone=1;
    if(nb_cortbone>IntensityProfile_mn+1) nb_cortbone=IntensityProfile_mn+1;
    nb_spongiousbone=IntensityProfile_mn+1-nb_cortbone; if(nb_spongiousbone<0) nb_spongiousbone=0;
    for(j=0;j<nb_spongiousbone;j++)  IntensityProfile[i][j]=levels[4];
    for(j=nb_spongiousbone;j<=IntensityProfile_mn;j++)  IntensityProfile[i][j]=levels[0];
    for(j=IntensityProfile_mn+1;j<IntensityProfile_pn+IntensityProfile_mn+1;j++)  IntensityProfile[i][j]=levels[2];
    }
}

unsigned short* CSimplexSurf::ComputeLevelsIntensityProfile(int nb_levels)
{
int i,j;
int histo[256]; for(i=0;i<256;i++) histo[i]=0;
for(i=0;i<nb_points;i++)  for(j=0;j<IntensityProfile_pn+IntensityProfile_mn+1;j++)  histo[IntensityProfile[i][j]]++;

// compute 3 maxima (cortical bone, muscles, fat/bone) and 2 minima
unsigned short* ret=new unsigned short[nb_levels*2-1];  
int max1[2]={0,0},max2[2]={0,0},max3[2]={0,0},min1[2]={1E10,0},min2[2]={1E10,0};
for(i=0;i<25;i++) if(histo[i]>max1[0]) {max1[0]=histo[i]; max1[1]=i;}
for(i=25;i<80;i++) if(histo[i]>max2[0]) {max2[0]=histo[i]; max2[1]=i;}
for(i=80;i<256;i++) if(histo[i]>max3[0]) {max3[0]=histo[i]; max3[1]=i;}
for(i=max1[1];i<max2[1];i++) if(histo[i]<min1[0]) {min1[0]=histo[i]; min1[1]=i;}
for(i=max2[1];i<max3[1];i++) if(histo[i]<min2[0]) {min2[0]=histo[i]; min2[1]=i;}

ret[0]=max1[1]; ret[1]=min1[1]; ret[2]=max2[1]; ret[3]=min2[1]; ret[4]=max3[1];
return ret;
}


void CSimplexSurf::RegisterIntensityProfile(int mn,int pn,P_float s,P_float searchdepth, P_float searchstep)
{
int i,j,norm;

P_float* ip0=new P_float[mn+pn+1];P_float* ip1=new P_float[mn+pn+1];P_float* ip2=new P_float[mn+pn+1];P_float* ip3=new P_float[mn+pn+1];
P_float S1,S2,S3,S1min,S2min,S3min,i1,i2,i3,w,val1,val2,val3,ind,ptemp[3],p0[3],p1[3],p2[3],p3[3],n0[3],n1[3],n2[3],n3[3];
int neighb[3];

vtkStructuredPoints* vol;
for(int k=0;k<nb_points;k++)  
    {
    GetNeighbors(k,neighb);
    GetPoint(k,p0); GetPoint(neighb[0],p1);GetPoint(neighb[1],p2);GetPoint(neighb[2],p3);
    GetNormal(k,n0);    GetNormal(neighb[0],n1);GetNormal(neighb[1],n2);GetNormal(neighb[2],n3);

    // handle possible MRI transformation
	vol=MRI[0];   

    // compute profiles
    for(j=-mn;j<=pn;j++)  ip0[j+mn]=Interpolation(vol,p0[0]+j*s*n0[0],p0[1]+j*s*n0[1],p0[2]+j*s*n0[2],Interpolationmode);
    for(j=-mn;j<=pn;j++)   ip1[j+mn]=Interpolation(vol,p1[0]+j*s*n1[0],p1[1]+j*s*n1[1],p1[2]+j*s*n1[2],Interpolationmode);
    for(j=-mn;j<=pn;j++)   ip2[j+mn]=Interpolation(vol,p2[0]+j*s*n2[0],p2[1]+j*s*n2[1],p2[2]+j*s*n2[2],Interpolationmode);
    for(j=-mn;j<=pn;j++)   ip3[j+mn]=Interpolation(vol,p3[0]+j*s*n3[0],p3[1]+j*s*n3[1],p3[2]+j*s*n3[2],Interpolationmode);

    // compute initial similarity
    i1=0;i2=0;i3=0;

    S1min=0; S2min=0; S3min=0;
    for(j=0;j<pn+mn+1;j++)
        {
        S1min+=fabs(ip1[j]-ip0[j]); S2min+=fabs(ip2[j]-ip0[j]); S3min+=fabs(ip3[j]-ip0[j]);
        }
    S1min=S1min/((P_float)(pn+mn+1)+0.0001); S2min=S2min/((P_float)(pn+mn+1)+0.0001); S3min=S3min/((P_float)(pn+mn+1)+0.0001);

    // register neighbors profiles to the ref profile
    for(i=-searchdepth;i<searchdepth;i++)
        {
        S1=0; S2=0; S3=0;
        norm=0;
        for(j=-mn;j<=pn;j++)
            {
            ind=(i*searchstep)/s+j;
            if(ind>=-mn && ind<=pn)
                {
                w=ind-floor(ind);
                val1=(1-w)*ip1[(int)floor(ind)+mn]+w*ip1[(int)ceil(ind)+mn]; val2=(1-w)*ip2[(int)floor(ind)+mn]+w*ip2[(int)ceil(ind)+mn]; val3=(1-w)*ip3[(int)floor(ind)+mn]+w*ip3[(int)ceil(ind)+mn];
                S1+=fabs(val1-ip0[j+mn]); S2+=fabs(val2-ip0[j+mn]); S3+=fabs(val3-ip0[j+mn]);
                norm++;
                }
            }
        S1=S1/((P_float)norm+0.0001); S2=S2/((P_float)norm+0.0001); S3=S3/((P_float)norm+0.0001);
        if(S1<S1min) {S1min=S1; i1=i;} if(S2<S2min) {S2min=S2; i2=i;} if(S3<S3min) {S3min=S3; i3=i;}
        }
    // apply transform
    ptemp[0]=p1[0]+i1*searchstep*n1[0];ptemp[1]=p1[1]+i1*searchstep*n1[1];ptemp[2]=p1[2]+i1*searchstep*n1[2]; SetPoint(neighb[0],ptemp);
    ptemp[0]=p2[0]+i2*searchstep*n2[0];ptemp[1]=p2[1]+i2*searchstep*n2[1];ptemp[2]=p2[2]+i2*searchstep*n2[2]; SetPoint(neighb[1],ptemp);
    ptemp[0]=p3[0]+i3*searchstep*n3[0];ptemp[1]=p3[1]+i3*searchstep*n3[1];ptemp[2]=p3[2]+i3*searchstep*n3[2]; SetPoint(neighb[2],ptemp);
    }

free(ip0); free(ip1); free(ip2); free(ip3); 
UpdateAll();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::InsertInternalConstraintPoint(P_float p[3]) {InsertInternalConstraintPoint(p[0],p[1],p[2]);}
void CSimplexSurf::InsertExternalConstraintPoint(P_float p[3]) {InsertExternalConstraintPoint(p[0],p[1],p[2]);}
void CSimplexSurf::InsertFrontierConstraintPoint(P_float p[3]) {InsertFrontierConstraintPoint(p[0],p[1],p[2]);}
void CSimplexSurf::InsertInternalConstraintPoint(P_float x, P_float y, P_float z) { InternalConstraintPoints->InsertNextPoint((double)x,(double)y,(double)z);}
void CSimplexSurf::InsertExternalConstraintPoint(P_float x, P_float y, P_float z) { ExternalConstraintPoints->InsertNextPoint((double)x,(double)y,(double)z);}
void CSimplexSurf::InsertFrontierConstraintPoint(P_float x, P_float y, P_float z) { FrontierConstraintPoints->InsertNextPoint((double)x,(double)y,(double)z);}
void CSimplexSurf::SetInternalConstraintPointsFromRefModel() {InternalConstraintPoints->DeepCopy(RefModel->GetPoints());}
void CSimplexSurf::InsertConstraintPointForce(int index,P_float f[3]) {constraintpointsforces[4*index]++; constraintpointsforces[4*index+1]+=f[0]; constraintpointsforces[4*index+2]+=f[1]; constraintpointsforces[4*index+3]+=f[2];}
void CSimplexSurf::InsertConstraintPointForce(int index,P_float x,P_float y,P_float z) {constraintpointsforces[4*index]++; constraintpointsforces[4*index+1]+=x; constraintpointsforces[4*index+2]+=y; constraintpointsforces[4*index+3]+=z;}

void CSimplexSurf::DeleteConstraintPoint(P_float x, P_float y, P_float z,P_float radius) 
{ 
int i;
double p[3],l;
vtkPoints* intpts=vtkPoints::New();
vtkPoints* extpts=vtkPoints::New();
vtkPoints* frpts=vtkPoints::New();

for(i=0;i<InternalConstraintPoints->GetNumberOfPoints();i++)
	{
	InternalConstraintPoints->GetPoint(i,p);
	l=sqrt((p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y)+(p[2]-z)*(p[2]-z));
	if(l>radius) intpts->InsertNextPoint(p);
	}
for(i=0;i<ExternalConstraintPoints->GetNumberOfPoints();i++)
	{
	ExternalConstraintPoints->GetPoint(i,p);
	l=sqrt((p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y)+(p[2]-z)*(p[2]-z));
	if(l>radius) extpts->InsertNextPoint(p);
	}
for(i=0;i<FrontierConstraintPoints->GetNumberOfPoints();i++)
	{
	FrontierConstraintPoints->GetPoint(i,p);
	l=sqrt((p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y)+(p[2]-z)*(p[2]-z));
	if(l>radius) frpts->InsertNextPoint(p);
	}
InternalConstraintPoints->DeepCopy(intpts);ExternalConstraintPoints->DeepCopy(extpts);FrontierConstraintPoints->DeepCopy(frpts);
intpts->Delete();extpts->Delete();frpts->Delete();
}

void CSimplexSurf::SetInternalConstraintPoints(vtkPoints* pts) {if(InternalConstraintPoints!=NULL) InternalConstraintPoints->Delete(); InternalConstraintPoints=pts;}
void CSimplexSurf::SetExternalConstraintPoints(vtkPoints* pts) {if(ExternalConstraintPoints!=NULL) ExternalConstraintPoints->Delete(); ExternalConstraintPoints=pts;}
void CSimplexSurf::SetFrontierConstraintPoints(vtkPoints* pts) {if(FrontierConstraintPoints!=NULL) FrontierConstraintPoints->Delete(); FrontierConstraintPoints=pts;}

vtkPoints* CSimplexSurf::GetInternalConstraintPoints() { return InternalConstraintPoints;}
vtkPoints* CSimplexSurf::GetExternalConstraintPoints() { return ExternalConstraintPoints;}
vtkPoints* CSimplexSurf::GetFrontierConstraintPoints() { return FrontierConstraintPoints;}
void CSimplexSurf::RemoveAllInternalConstraintPoints() {InternalConstraintPoints->Delete(); InternalConstraintPoints=vtkPoints::New(); InternalConstraintPoints->SetDataTypeToFloat(); }
void CSimplexSurf::RemoveAllExternalConstraintPoints() {ExternalConstraintPoints->Delete(); ExternalConstraintPoints=vtkPoints::New(); ExternalConstraintPoints->SetDataTypeToFloat(); }
void CSimplexSurf::RemoveAllFrontierConstraintPoints() {FrontierConstraintPoints->Delete(); FrontierConstraintPoints=vtkPoints::New(); FrontierConstraintPoints->SetDataTypeToFloat(); }

void CSimplexSurf::SaveInternalConstraintPoints(const char* filename)
{
int i; double dp[3];
FILE* f=fopen(filename,"wt");
fprintf(f,"%d \n",InternalConstraintPoints->GetNumberOfPoints());
for(i=0;i<InternalConstraintPoints->GetNumberOfPoints();i++) 
    {
    InternalConstraintPoints->GetPoint(i,dp);
    fprintf(f,"%4.2f %4.2f %4.2f\n",dp[0],dp[1],dp[2]);
    }
fclose(f);
}
void CSimplexSurf::LoadInternalConstraintPoints(const char* filename)
{
std::string name=filename;
if(name.find(".vtk")!=-1)
    {
    vtkPolyDataReader* reader=vtkPolyDataReader::New();
        reader->SetFileName(filename);
        reader->Update();
    InternalConstraintPoints->DeepCopy(reader->GetOutput()->GetPoints());
    reader->Delete();
    }
else
    {
    int i,nb; P_float p[3];
    FILE* f=fopen(filename,"rt");
    InternalConstraintPoints->Reset();
    fscanf(f,"%d \n",&nb);
    for(i=0;i<nb;i++) 
        {
        fscanf(f,"%lf %lf %lf\n",&p[0],&p[1],&p[2]);
        InsertInternalConstraintPoint(p);
        }
    fclose(f);
    }
}

void CSimplexSurf::SaveExternalConstraintPoints(const char* filename)
{
int i; double dp[3];
FILE* f=fopen(filename,"wt");
fprintf(f,"%d \n",ExternalConstraintPoints->GetNumberOfPoints());
for(i=0;i<ExternalConstraintPoints->GetNumberOfPoints();i++) 
    {
    ExternalConstraintPoints->GetPoint(i,dp);
    fprintf(f,"%4.2f %4.2f %4.2f\n",dp[0],dp[1],dp[2]);
    }
fclose(f);
}
void CSimplexSurf::LoadExternalConstraintPoints(const char* filename)
{
std::string name=filename;
if(name.find(".vtk")!=-1)
    {
    vtkPolyDataReader* reader=vtkPolyDataReader::New();
        reader->SetFileName(filename);
        reader->Update();
    ExternalConstraintPoints->DeepCopy(reader->GetOutput()->GetPoints());
    reader->Delete();
    }
else
    {
    int i,nb; P_float p[3];
    FILE* f=fopen(filename,"rt");
    ExternalConstraintPoints->Reset();
    fscanf(f,"%d \n",&nb);
    for(i=0;i<nb;i++) 
        {
        fscanf(f,"%lf %lf %lf\n",&p[0],&p[1],&p[2]);
        InsertExternalConstraintPoint(p);
        }
    fclose(f);
    }
}
void CSimplexSurf::SaveFrontierConstraintPoints(const char* filename)
{
int i; double dp[3];
FILE* f=fopen(filename,"wt");
fprintf(f,"%d \n",FrontierConstraintPoints->GetNumberOfPoints());
for(i=0;i<FrontierConstraintPoints->GetNumberOfPoints();i++) 
    {
    FrontierConstraintPoints->GetPoint(i,dp);
    fprintf(f,"%4.2f %4.2f %4.2f\n",dp[0],dp[1],dp[2]);
    }
fclose(f);
}
void CSimplexSurf::LoadFrontierConstraintPoints(const char* filename)
{
std::string name=filename;
if(name.find(".vtk")!=-1)
    {
    vtkPolyDataReader* reader=vtkPolyDataReader::New();
        reader->SetFileName(filename);
        reader->Update();
    FrontierConstraintPoints->DeepCopy(reader->GetOutput()->GetPoints());
    reader->Delete();
    }
else
    {
    int i,nb; P_float p[3];
    FILE* f=fopen(filename,"rt");
    FrontierConstraintPoints->Reset();
    fscanf(f,"%d \n",&nb);
    for(i=0;i<nb;i++) 
        {
        fscanf(f,"%lf %lf %lf\n",&p[0],&p[1],&p[2]);
        InsertFrontierConstraintPoint(p);
        }
    fclose(f);
    }
}


void CSimplexSurf::SetRefModel(vtkPolyData* refmodel,int inimode)
{
RefPointsLocator=vtkPointLocator::New();
RefModel=vtkPolyData::New();
RefModel->DeepCopy(refmodel);
RefPointsLocator->SetDataSet(RefModel);
RefPointsLocator->BuildLocator();
IniRefModel(inimode);
SetInternalConstraintPointsFromRefModel();
}

void CSimplexSurf::IniRefModel(int mode)
{
int i; double dp[3],cref[3]={0,0,0},c[3]={0,0,0}; P_float p[3]; P_float nm;
// calcul normals
vtkPolyDataNormals *normals=vtkPolyDataNormals::New();
    normals->SetInput(RefModel);
    normals->SplittingOff();
    normals->ConsistencyOn();
    normals->ComputePointNormalsOn();
    normals->SetFlipNormals(0);
    normals->Update();

RefModel->Delete();
RefModel=vtkPolyData::New();
    RefModel->DeepCopy(normals->GetOutput());

// calcul centroides
for(i=0;i<RefModel->GetNumberOfPoints();i++) {RefModel->GetPoint(i,dp); cref[0]+=dp[0]; cref[1]+=dp[1]; cref[2]+=dp[2];}
cref[0]=cref[0]/RefModel->GetNumberOfPoints();cref[1]=cref[1]/RefModel->GetNumberOfPoints();cref[2]=cref[2]/RefModel->GetNumberOfPoints();

// calcul Radius/direction
P_float r=0;
for(i=0;i<RefModel->GetNumberOfPoints();i++) {RefModel->GetPoint(i,dp); p[0]=cref[0]-dp[0]; p[1]=cref[1]-dp[1]; p[2]=cref[2]-dp[2]; nm=norm(p); if(nm>r) r=nm;}
// scale model
if(mode==1)  /// SCALE SPHERE TO BOUND THE MODEL
{ for(i=0;i<GetNumberOfPoints();i++) {GetPoint(i,p); nm=norm(p); p[0]=p[0]*r/nm; p[1]=p[1]*r/nm; p[2]=p[2]*r/nm; SetPoint(i,p);} }

//ComputeVolume(); SetVolume_ref(volume,0);

// translate model
for(i=0;i<GetNumberOfPoints();i++) {GetPoint(i,p); c[0]+=p[0]; c[1]+=p[1]; c[2]+=p[2];}
c[0]=cref[0]-c[0]/GetNumberOfPoints();c[1]=cref[1]-c[1]/GetNumberOfPoints();c[2]=cref[2]-c[2]/GetNumberOfPoints();

if(mode==1) for(i=0;i<GetNumberOfPoints();i++) {GetPoint(i,p); SetPoint(i,p[0]+c[0],p[1]+c[1],p[2]+c[2]);}
Equilibrium();

// pointcells 
RefModel->BuildCells(); 
RefModel->Modified();
RefModelPointCells=new vtkIdList*[RefModel->GetNumberOfPoints()];
for(i=0;i<RefModel->GetNumberOfPoints();i++) 
    {
    RefModelPointCells[i]=vtkIdList::New();
    RefModel->GetPointCells(i,RefModelPointCells[i]);
    }   
// cell normals & volume
RefModelCellNormals=new P_float[3*RefModel->GetNumberOfCells()];
P_float p1[3],p2[3],p3[3],p1p2[3],p1p3[3],n[3];
int ids[3];
RefModelVolume=0; RefModelSurface=0;
for(i=0;i<RefModel->GetNumberOfCells();i++) 
    {
    ids[0]=RefModel->GetCell(i)->GetPointIds()->GetId(0);       ids[1]=RefModel->GetCell(i)->GetPointIds()->GetId(1);           ids[2]=RefModel->GetCell(i)->GetPointIds()->GetId(2);
    p1[0]=RefModel->GetPoint(ids[0])[0]; p1[1]=RefModel->GetPoint(ids[0])[1]; p1[2]=RefModel->GetPoint(ids[0])[2];
    p2[0]=RefModel->GetPoint(ids[1])[0]; p2[1]=RefModel->GetPoint(ids[1])[1]; p2[2]=RefModel->GetPoint(ids[1])[2];
    p3[0]=RefModel->GetPoint(ids[2])[0]; p3[1]=RefModel->GetPoint(ids[2])[1]; p3[2]=RefModel->GetPoint(ids[2])[2];
    p1p2[0]=p2[0]-p1[0];    p1p2[1]=p2[1]-p1[1];    p1p2[2]=p2[2]-p1[2];
    p1p3[0]=p3[0]-p1[0];    p1p3[1]=p3[1]-p1[1];    p1p3[2]=p3[2]-p1[2];
    crossproduct(n,p1p2,p1p3);  nm=norm(n);	RefModelSurface+=nm/2;
    RefModelVolume+=n[2] * (p1[2]+p2[2]+p3[2]) /6;
    RefModelCellNormals[3*i]=n[0]/nm;   RefModelCellNormals[3*i+1]=n[1]/nm;     RefModelCellNormals[3*i+2]=n[2]/nm;
    }   
SetVolume_ref(1,1);	SetSurface_ref(1,1);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::SwicthOffAllForces() {InternalForce_AxialConstraint=false;InternalForce_SurfEqu=false;InternalForce_SurfEquFlexion=false;InternalForce_Laplacian=false;InternalForce_LaplacianFlexion=false;InternalForce_RefShape=false;InternalForce_VolPres=false;InternalForce_SurfPres=false;ExternalForce_ICP=false;ExternalForce_GradientMRI=false;ExternalForce_Demon=false;ExternalForce_IntensityProfile=false;ExternalForce_ConstraintPoints=false;ExternalForce_ConstraintMesh=false;}
void CSimplexSurf::SetInternalForceAxialConstraint(bool i,int mode,P_float alpha,P_float alpha_lim,P_float alpha_incr) {AxialConstraint_mode=mode; InternalForce_AxialConstraint=i;Alpha_AxialConstraint=alpha;Alpha_AxialConstraint_lim=alpha_lim;Alpha_AxialConstraint_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceLaplacian(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_Laplacian=i;Alpha_Laplacian=alpha;Alpha_Laplacian_lim=alpha_lim;Alpha_Laplacian_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceSurfEqu(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_SurfEqu=i;Alpha_SurfEqu=alpha;Alpha_SurfEqu_lim=alpha_lim;Alpha_SurfEqu_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceLaplacianFlexion(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_LaplacianFlexion=i;Alpha_LaplacianFlexion=alpha;Alpha_LaplacianFlexion_lim=alpha_lim;Alpha_LaplacianFlexion_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceSurfEquFlexion(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_SurfEquFlexion=i;Alpha_SurfEquFlexion=alpha;Alpha_SurfEquFlexion_lim=alpha_lim;Alpha_SurfEquFlexion_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceRefShape(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_RefShape=i;Alpha_RefShape=alpha;Alpha_RefShape_lim=alpha_lim;Alpha_RefShape_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceVolPres(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_VolPres=i;Alpha_VolPres=alpha;Alpha_VolPres_lim=alpha_lim;Alpha_VolPres_incr=alpha_incr;}
void CSimplexSurf::SetInternalForceSurfPres(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {InternalForce_SurfPres=i;Alpha_SurfPres=alpha;Alpha_SurfPres_lim=alpha_lim;Alpha_SurfPres_incr=alpha_incr;}
//void CSimplexSurf::SetInternalForceGaussianFilterDisplacement(bool i,P_float rigidity) {InternalForce_GaussianFilterDisplacement=i; if(rigidity!=-1) SetRigidity(rigidity);}
void CSimplexSurf::SetExternalForceICP(bool i,double dmax,P_float alpha,P_float alpha_lim,P_float alpha_incr) {  ExternalForce_ICP=i;Alpha_ICP=alpha;DMax_ICP=dmax;Alpha_ICP_lim=alpha_lim;Alpha_ICP_incr=alpha_incr; }

void CSimplexSurf::SetExternalForceGradientMRI(bool i, double size, int depth, P_float alpha, P_float alpha_lim, P_float alpha_incr, vtkStructuredPoints* gradvol, bool oppositedirection) {
	if (gradvol != NULL) 
		SetMRI_GradientMRI(gradvol); 

	if (MRI_grad != NULL && i && MRI_GradientMRI == NULL) 
		SetMRI_GradientMRI(MRI_grad[0]); 
	
	ExternalForce_GradientMRI = i; 
	Alpha_GradientMRI = alpha; 
	Alpha_GradientMRI_lim = alpha_lim; 
	Alpha_GradientMRI_incr = alpha_incr; 
	Size_GradientMRI = size; 
	depth_GradientMRI = depth;
	oppositedirection_GradientMRI = oppositedirection;
}

void CSimplexSurf::SetExternalForceDemon(bool i,int depth,P_float alpha,P_float alpha_lim,P_float alpha_incr) {ExternalForce_Demon=i;Alpha_Demon=alpha;Alpha_Demon_lim=alpha_lim;Alpha_Demon_incr=alpha_incr;depth_Demon=depth;}
void CSimplexSurf::SetExternalForceIntensityProfile(bool i,int metric,P_float depth,P_float alpha,P_float alpha_lim,P_float alpha_incr,P_float depth_lim,P_float depth_incr) {ExternalForce_IntensityProfile=i;Alpha_IntensityProfile=alpha;Alpha_IntensityProfile_lim=alpha_lim;Alpha_IntensityProfile_incr=alpha_incr;Metric_IntensityProfile=metric;Depth_IntensityProfile=depth;Depth_IntensityProfile_float=depth; Depth_IntensityProfile_incr=depth_incr;Depth_IntensityProfile_lim=depth_lim;}
void CSimplexSurf::SetExternalForceConstraintPoints(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr) {ExternalForce_ConstraintPoints=i;Alpha_ConstraintPoints=alpha;Alpha_ConstraintPoints_lim=alpha_lim;Alpha_ConstraintPoints_incr=alpha_incr;}
void CSimplexSurf::SetExternalForceConstraintMesh(bool i,P_float alpha,P_float alpha_lim,P_float alpha_incr,P_float tresh,CSimplexSurf* mesh) {ExternalForce_ConstraintMesh_mesh=mesh; ExternalForce_ConstraintMesh=i;Alpha_ConstraintMesh=alpha;Alpha_ConstraintMesh_lim=alpha_lim;Alpha_ConstraintMesh_incr=alpha_incr;ExternalForce_ConstraintMesh_tresh=tresh; Init_ConstraintMeshForces();}

bool CSimplexSurf::GetInternalForce_AxialConstraint() {return InternalForce_AxialConstraint;}
void CSimplexSurf::SetAxialConstraint_mode(int mode) {AxialConstraint_mode=mode;}

void CSimplexSurf::SetHandleborder(bool hb) {Forcehandleborder=hb; UpdateParams();}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::Init_ConstraintMeshForces()
{
// ad-hoc technique for cartilage
if(ExternalForce_ConstraintMesh_pts!=NULL) free(ExternalForce_ConstraintMesh_pts);
if(!ExternalForce_ConstraintMesh) return;
int i,j,ID; P_float d,dmin;
ExternalForce_ConstraintMesh_pts=new int[nb_points];
for(i=0;i<nb_points;i++)
	{
	dmin=ExternalForce_ConstraintMesh_tresh; ID=-1;
	// get closest facing point
	for(j=0;j<ExternalForce_ConstraintMesh_mesh->GetNumberOfPoints();j++)
		{
		d=dist3D(GetPoint(i),ExternalForce_ConstraintMesh_mesh->GetPoint(j));
		if(d<dmin && dotproduct(GetNormal(i),ExternalForce_ConstraintMesh_mesh->GetNormal(j))<0) {dmin=d; ID=j;}
		}
	ExternalForce_ConstraintMesh_pts[i]=ID;
	}
}

void CSimplexSurf::UpdateForcesConstraintMesh()
{
// ad-hoc technique for cartilage
CSimplexSurf *mesh=ExternalForce_ConstraintMesh_mesh;
int jmax; P_float wmax;
CLOSESTSTRUCT *closest=new CLOSESTSTRUCT[1]; 
P_float PPref[3],alpha=Alpha_ConstraintMesh,extrapolation=ISOTROPICEXTRAPOLATION;

for(int i=0;i<nb_points;i++)
	if(ExternalForce_ConstraintMesh_pts[i]!=-1)
		{
		closest->weights=new P_float[1]; closest->pts=new int[1];     closest->dist2=1E10; closest->nb=0; 
		mesh->GetClosest_P(closest,GetPoint(i),ExternalForce_ConstraintMesh_pts[i]); 
		wmax=0; for(int j=0;j<closest->nb;j++) if(closest->weights[j]>wmax) {wmax=closest->weights[j]; jmax=closest->pts[j];} ExternalForce_ConstraintMesh_pts[i]=jmax;

		PPref[0]=closest->c[0]-GetPoint(i)[0]; PPref[1]=closest->c[1]-GetPoint(i)[1]; PPref[2]=closest->c[2]-GetPoint(i)[2]; 
		
		Add_indep(i,PPref,alpha,extrapolation);
		
		free(closest->pts); free(closest->weights);
		}
free(closest);
}



void CSimplexSurf::GetClosestRefModel(P_float p[3],P_float n[3],P_float u[3])
{
int i,index_cell,ids[3];
u[0]=0; u[1]=0; u[2]=0; 

P_float ntest; 
double dp[3]={(double)p[0],(double)p[1],(double)p[2]};
double d; int ID=RefPointsLocator->FindClosestPointWithinRadius(DMax_ICP,dp,d);
P_float p1[3],p2[3],p3[3],cp[3];
P_float cpp[3],*norml;
if(ID==-1) return;

P_float fa,fb,fc,fd,fe,ff,fout[2],fproj[3],fdist,fdist_sav=1E10;
P_float p1p2[3],p1p3[3],pp1[3];

// get the closest point on the surface
cp[0]=RefModel->GetPoint(ID)[0]; cp[1]=RefModel->GetPoint(ID)[1]; cp[2]=RefModel->GetPoint(ID)[2];
cpp[0]=p[0]-cp[0];cpp[1]=p[1]-cp[1];cpp[2]=p[2]-cp[2];  
for(i=0;i<RefModelPointCells[ID]->GetNumberOfIds();i++)
    {
    index_cell=RefModelPointCells[ID]->GetId(i);
    norml=RefModelCellNormals+3*index_cell;
    
    ids[0]=RefModel->GetCell(index_cell)->GetPointIds()->GetId(0);     ids[1]=RefModel->GetCell(index_cell)->GetPointIds()->GetId(1);         ids[2]=RefModel->GetCell(index_cell)->GetPointIds()->GetId(2);
    p1[0]=RefModel->GetPoint(ids[0])[0]; p1[1]=RefModel->GetPoint(ids[0])[1]; p1[2]=RefModel->GetPoint(ids[0])[2]; p2[0]=RefModel->GetPoint(ids[1])[0]; p2[1]=RefModel->GetPoint(ids[1])[1]; p2[2]=RefModel->GetPoint(ids[1])[2]; p3[0]=RefModel->GetPoint(ids[2])[0]; p3[1]=RefModel->GetPoint(ids[2])[1]; p3[2]=RefModel->GetPoint(ids[2])[2];
    p1p2[0]=p2[0]-p1[0];p1p2[1]=p2[1]-p1[1];p1p2[2]=p2[2]-p1[2]; p1p3[0]=p3[0]-p1[0];p1p3[1]=p3[1]-p1[1];p1p3[2]=p3[2]-p1[2]; pp1[0]=p1[0]-p[0];pp1[1]=p1[1]-p[1];pp1[2]=p1[2]-p[2];

//  dot=dotproduct(cpp,norml); 
//  f[0]=-dot*ntest*n[0]; f[1]=-dot*ntest*n[0]; f[2]=-dot*ntest*n[0]; 
//  f[0]=-dot*norml[0]; f[1]=-dot*norml[1]; f[2]=-dot*norml[2]; 

    ntest=dotproduct(n,norml); 
    fa=dotproduct(p1p2,p1p2); fb=dotproduct(p1p2,p1p3); fc=dotproduct(p1p3,p1p3); fd=dotproduct(p1p2,pp1); fe=dotproduct(p1p3,pp1); ff=dotproduct(pp1,pp1);
    fdist=Closest(fa,fb,fc,fd,fe,ff,fout);

    if(ntest>0 && fdist<fdist_sav)  
        {
        fproj[0]=p1[0]+fout[0]*p1p2[0]+fout[1]*p1p3[0];fproj[1]=p1[1]+fout[0]*p1p2[1]+fout[1]*p1p3[1];fproj[2]=p1[2]+fout[0]*p1p2[2]+fout[1]*p1p3[2];
        u[0]=fproj[0]-p[0]; u[1]=fproj[1]-p[1]; u[2]=fproj[2]-p[2]; 
        fdist_sav=fdist;
        } 
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateForcesLaplacianFlexion()
{
int *neighb;
P_float H,PPref[3],pref[3]={1./3.,1./3.,1./3.},C[3],*P,*P1,*P2,*P3,*n,alpha=Alpha_LaplacianFlexion,extrapolation=INTERNALFORCEEXTRAPOLATION;
for(int i=0;i<nb_points;i++)
    {
	neighb=GetNeighbors(i); 
	P=GetPoint(i); P1=GetPoint(neighb[0]); P2=GetPoint(neighb[1]); P3=GetPoint(neighb[2]);
	
	if(neighbors_c[3*i+2]==-1) {n=normals2+3*i; C[0]=(P1[0]+P2[0])/2; C[1]=(P1[1]+P2[1])/2; C[2]=(P1[2]+P2[2])/2; }
	else {n=normals+3*i; C[0]=(P1[0]+P2[0]+P3[0])/3; C[1]=(P1[1]+P2[1]+P3[1])/3; C[2]=(P1[2]+P2[2]+P3[2])/3;}

	H=GetHmean(i);

	PPref[0]=C[0]-P[0]-H*n[0];PPref[1]=C[1]-P[1]-H*n[1];PPref[2]=C[2]-P[2]-H*n[2];
	
	if(USEMULTIPLEPARTICLES) Add_multi(i,pref,H,alpha,extrapolation);
	else Add_indep(i,PPref,alpha,extrapolation);
	}
}

void CSimplexSurf::UpdateForcesLaplacian()
{
int *neighb;
P_float PPref[3],pref[3]={1./3.,1./3.,1./3.},C[3],*P,*P1,*P2,*P3,alpha=Alpha_Laplacian,extrapolation=INTERNALFORCEEXTRAPOLATION;
for(int i=0;i<nb_points;i++)
    {
	neighb=GetNeighbors(i); 
	P=GetPoint(i); 
	P1=GetPoint(neighb[0]); 
	P2=GetPoint(neighb[1]); 
	P3=GetPoint(neighb[2]);
	
	if(neighbors_c[3*i+2]==-1) {
		C[0]=(P1[0]+P2[0])/2; 
		C[1]=(P1[1]+P2[1])/2; 
		C[2]=(P1[2]+P2[2])/2; 
	}
	else {
		C[0]=(P1[0]+P2[0]+P3[0])/3; 
		C[1]=(P1[1]+P2[1]+P3[1])/3;
		C[2]=(P1[2]+P2[2]+P3[2])/3;
	}
	
	PPref[0]=C[0]-P[0];
	PPref[1]=C[1]-P[1];
	PPref[2]=C[2]-P[2];
	
	if(USEMULTIPLEPARTICLES) Add_multi(i,pref,0,alpha,extrapolation);
	else Add_indep(i,PPref,alpha,extrapolation);
	}
}

void CSimplexSurf::UpdateForcesSurfEquFlexion()
{
int *neighb;
P_float H,PPref[3],pref[3],prefm,*P,*P1,*P2,*P3,*n,alpha=Alpha_SurfEquFlexion,extrapolation=INTERNALFORCEEXTRAPOLATION;
for(int i=0;i<nb_points;i++)
    {
	neighb=GetNeighbors(i); 
	P=GetPoint(i); P1=GetPoint(neighb[0]); P2=GetPoint(neighb[1]); P3=GetPoint(neighb[2]);

	pref[0]=surfaces[neighb[0]]; pref[1]=surfaces[neighb[1]]; pref[2]=surfaces[neighb[2]];

	if(neighbors_c[3*i+2]==-1) {n=normals2+3*i; pref[2]=0;}
	else {n=normals+3*i;}

	prefm=(pref[0]+pref[1]+pref[2]); pref[0]=pref[0]/prefm;pref[1]=pref[1]/prefm;pref[2]=pref[2]/prefm;

	H=GetHmean(i);

	PPref[0]=P1[0]*pref[0]+P2[0]*pref[1]+P3[0]*pref[2]-P[0]-H*n[0]; PPref[1]=P1[1]*pref[0]+P2[1]*pref[1]+P3[1]*pref[2]-P[1]-H*n[1]; PPref[2]=P1[2]*pref[0]+P2[2]*pref[1]+P3[2]*pref[2]-P[2]-H*n[2];
	
	if(USEMULTIPLEPARTICLES) Add_multi(i,pref,H,alpha,extrapolation);
	else Add_indep(i,PPref,alpha,extrapolation);
	}
}


void CSimplexSurf::UpdateForcesSurfEqu() {
	int *neighb;
	P_float PPref[3], pref[3], prefm, *P, *P1, *P2, *P3, alpha = Alpha_SurfEqu, extrapolation = INTERNALFORCEEXTRAPOLATION;
	
	for (int i = 0; i < nb_points; i++) { 
		neighb = GetNeighbors(i); 
		P = GetPoint(i); 
		P1 = GetPoint(neighb[0]); 
		P2 = GetPoint(neighb[1]); 
		P3 = GetPoint(neighb[2]);

		pref[0] = surfaces[neighb[0]]; 
		pref[1] = surfaces[neighb[1]]; 
		pref[2] = surfaces[neighb[2]];

		if (neighbors_c[3 * i + 2] == -1) pref[2] = 0;

		prefm = (pref[0] + pref[1] + pref[2]); 
		pref[0] = pref[0] / prefm;
		pref[1] = pref[1] / prefm;
		pref[2] = pref[2] / prefm;

		PPref[0] = P1[0] * pref[0] + P2[0] * pref[1] + P3[0] * pref[2] - P[0]; 
		PPref[1] = P1[1] * pref[0] + P2[1] * pref[1] + P3[1] * pref[2] - P[1]; 
		PPref[2] = P1[2] * pref[0] + P2[2] * pref[1] + P3[2] * pref[2] - P[2];
	
		if (extrapolation != 0) {
			cout << " ";
		}
		
		if (USEMULTIPLEPARTICLES) Add_multi(i, pref, 0, alpha, extrapolation);
		else Add_indep(i, PPref, alpha, extrapolation);
	}
}

void CSimplexSurf::UpdateForcesRefShape() {
	int *neighb;
	P_float ht, s, P1P2[3], P1P3[3], P2P3[3], S[3], PPref[3], pref[3], *n, *P, *P1, *P2, *P3, alpha = Alpha_RefShape, extrapolation = INTERNALFORCEEXTRAPOLATION;
	
	for (int i = 0; i < nb_points; i++) {
		if (mass_constraint[i] == NULL) {// hack to prevent skin border from instbilities		
			neighb = GetNeighbors(i); 
			GetParams(i, pref); 
			
			P = GetPoint(i); 
			P1 = GetPoint(neighb[0]);
			P2 = GetPoint(neighb[1]); 
			P3 = GetPoint(neighb[2]);

			// get surface
			if (Forcehandleborder && neighbors_c[3 * i + 2] == -1) {
				n = normals2 + 3 * i;
				P1P2[0] = P2[0] - P1[0]; P1P2[1] = P2[1] - P1[1]; P1P2[2] = P2[2] - P1[2];
				s = dotproduct(P1P2, P1P2); 
			}
			else { 
				n = normals + 3 * i; 
				P1P3[0] = P3[0] - P1[0]; 
				P1P3[1] = P3[1] - P1[1]; 
				P1P3[2] = P3[2] - P1[2]; 

				P2P3[0] = P3[0] - P2[0]; 
				P2P3[1] = P3[1] - P2[1];
				P2P3[2] = P3[2] - P2[2]; 

				crossproduct(S, P2P3, P1P3);  s = norm(S);
			}

			if (pref[2] == 0) {
				ht = 0;
			} 
			else {
				ht = pow(s, SCALEINVARIANT_REFSHAPE) / pref[2]; 
			}

			pref[2] = 1 - pref[0] - pref[1];
			PPref[0] = P1[0] * pref[0] + P2[0] * pref[1] + P3[0] * pref[2] - P[0] - ht * n[0]; 
			PPref[1] = P1[1] * pref[0] + P2[1] * pref[1] + P3[1] * pref[2] - P[1] - ht * n[1]; 
			PPref[2] = P1[2] * pref[0] + P2[2] * pref[1] + P3[2] * pref[2] - P[2] - ht * n[2];

			if (extrapolation != 0) {
				cout << " ";
			}
			
			if (USEMULTIPLEPARTICLES) {
				Add_multi(i, pref, ht, alpha, extrapolation);
			}
			else {
				Add_indep(i, PPref, alpha, extrapolation);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateForcesVolPres() {
	P_float dV = volume_ref - volume;
	P_float PPref[3], *n;
	P_float dH = dV / surface; // dH constant

	if (abs(dH) > Threshold_volsurfpres) {
		dH = (dH < 0)?(-Threshold_volsurfpres):(Threshold_volsurfpres);
	}

	P_float alpha = Alpha_VolPres, extrapolation = ANISOTROPICEXTRAPOLATION;
	// dV constant
	//if(surfaces[i]!=0) dH=dV/(surfaces[i]*nb_points); else dH=0;
	for (int i = 0; i < nb_points; i++) {
		n = GetNormal(i);
		PPref[0] = dH * n[0]; PPref[1] = dH * n[1]; PPref[2] = dH * n[2];
		Add_indep(i, PPref, alpha, extrapolation);
	}
}

void CSimplexSurf::UpdateForcesSurfPres()
{
P_float L,dS=surface_ref-surface;
L=0;  for(int i=0;i<nb_points;i++)	if(neighbors_c[3*i+2]==-1) 		{if(GetNeighbors(i)[0]<i) L+=dist3D(GetPoint(i),GetPoint(GetNeighbors(i)[0])); if(GetNeighbors(i)[1]<i) L+=dist3D(GetPoint(i),GetPoint(GetNeighbors(i)[1]));}
P_float PPref[3],*n,dH=dS/L; // dH constant
if(abs(dH)>Threshold_volsurfpres) dH=(dH<0)?(-Threshold_volsurfpres):(Threshold_volsurfpres);
P_float alpha=Alpha_SurfPres,extrapolation=ANISOTROPICEXTRAPOLATION;
//P_float nrmv,*P1=GetPoint(neighbors[3*i]),*P2=GetPoint(neighbors[3*i+1]),*P3=GetPoint(neighbors[3*i+2]); 
//P_float P1P2[3]={P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]}; nrmv=norm(P1P2);
//P_float dH=2*dS/(nrmv*nb_points_border);
for(int i=0;i<nb_points;i++)
	if(neighbors_c[3*i+2]==-1)
		{
		n=normals2+3*i;
		PPref[0]=dH*n[0]; PPref[1]=dH*n[1]; PPref[2]=dH*n[2];
		Add_indep(i,PPref,alpha,extrapolation);
		}
}

void CSimplexSurf::UpdateForcesAxialConstraint()
{
if(axiallinks==NULL) return;
P_float alpha=Alpha_AxialConstraint,extrapolation=ANISOTROPICEXTRAPOLATION;
for(int i=0;i<nb_points;i++)
	{
	Add_indep(i,axiallinks[i].f,alpha,extrapolation);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateForcesConstraintPoints()
{
int i,j; CLOSESTSTRUCT* closest=new CLOSESTSTRUCT[1]; double dblp[3]; P_float p[3],resp[3],W;

for(i=0;i<4*nb_points;i++) constraintpointsforces[i]=0;

P_float inf=10000000;
UpdateDOPCellTree(0);
for(i=0;i<ExternalConstraintPoints->GetNumberOfPoints();i++)
        {
        GetExternalConstraintPoints()->GetPoint(i,dblp); p[0]=dblp[0]; p[1]=dblp[1]; p[2]=dblp[2];
        closest->pts=new int[1];closest->weights=new P_float[1]; closest->dist2=1E10; closest->t=1E10; closest->nb=0; closest->crossing=false; closest->side=1;
		if(p[0]>=GetDOPCellTree()->ext[1] && p[1]<=GetDOPCellTree()->ext[2] && p[1]>=GetDOPCellTree()->ext[3] &&   p[2]<=GetDOPCellTree()->ext[4] && p[2]>=GetDOPCellTree()->ext[5] &&   (p[0]+p[1])>=GetDOPCellTree()->ext[7] && (p[0]-p[1])>=GetDOPCellTree()->ext[9] &&  (p[0]+p[2])>=GetDOPCellTree()->ext[11] &&   (p[0]-p[2])>=GetDOPCellTree()->ext[13] &&  (p[1]+p[2])<=GetDOPCellTree()->ext[14] && (p[1]+p[2])>=GetDOPCellTree()->ext[15] &&  (p[1]-p[2])<=GetDOPCellTree()->ext[16] && (p[1]-p[2])>=GetDOPCellTree()->ext[17] )  
			DOPInclusion(p,closest,DOPCellTree);
		free(closest->pts); free(closest->weights);
		if(closest->side==-1)       
			{
			closest->pts=new int[1];closest->weights=new P_float[1]; closest->dist2=1E10; closest->nb=0;
			GetClosest(closest,p); 
			if(closest->side==-1)       
				{
				resp[0]=p[0]-closest->c[0]; resp[1]=p[1]-closest->c[1]; resp[2]=p[2]-closest->c[2]; 
				W=0; for(j=0;j<closest->nb;j++) W+=closest->weights[j]*closest->weights[j]; if(W!=0) W=1/W;
				for(j=0;j<closest->nb;j++) InsertConstraintPointForce(closest->pts[j],closest->weights[j]*W*resp[0],closest->weights[j]*W*resp[1],closest->weights[j]*W*resp[2]);
				}
			free(closest->pts); free(closest->weights); 
			}
        }

for(i=0;i<InternalConstraintPoints->GetNumberOfPoints();i++)
        {
        GetInternalConstraintPoints()->GetPoint(i,dblp); p[0]=dblp[0]; p[1]=dblp[1]; p[2]=dblp[2];
        closest->pts=new int[1];closest->weights=new P_float[1]; closest->dist2=1E10; closest->t=1E10; closest->nb=0; closest->crossing=false; closest->side=1;
//		if(p[0]>=GetDOPCellTree()->ext[1] && p[1]<=GetDOPCellTree()->ext[2] && p[1]>=GetDOPCellTree()->ext[3] &&   p[2]<=GetDOPCellTree()->ext[4] && p[2]>=GetDOPCellTree()->ext[5] &&   (p[0]+p[1])>=GetDOPCellTree()->ext[7] && (p[0]-p[1])>=GetDOPCellTree()->ext[9] &&  (p[0]+p[2])>=GetDOPCellTree()->ext[11] &&   (p[0]-p[2])>=GetDOPCellTree()->ext[13] &&  (p[1]+p[2])<=GetDOPCellTree()->ext[14] && (p[1]+p[2])>=GetDOPCellTree()->ext[15] &&  (p[1]-p[2])<=GetDOPCellTree()->ext[16] && (p[1]-p[2])>=GetDOPCellTree()->ext[17] )  
            DOPInclusion(p,closest,DOPCellTree);
        free(closest->pts); free(closest->weights);
		if(closest->side==1)        
            {
			closest->pts=new int[1];closest->weights=new P_float[1]; closest->dist2=1E10; closest->nb=0;
            GetClosest(closest,p); 
	        if(closest->side==1)       
				{
				resp[0]=p[0]-closest->c[0]; resp[1]=p[1]-closest->c[1]; resp[2]=p[2]-closest->c[2]; 
				W=0; for(j=0;j<closest->nb;j++) W+=closest->weights[j]*closest->weights[j]; if(W!=0) W=1/W;
				for(j=0;j<closest->nb;j++) InsertConstraintPointForce(closest->pts[j],closest->weights[j]*W*resp[0],closest->weights[j]*W*resp[1],closest->weights[j]*W*resp[2]);
				}
            free(closest->pts); free(closest->weights); 
			}
        }

for(i=0;i<FrontierConstraintPoints->GetNumberOfPoints();i++)
        {
        GetFrontierConstraintPoints()->GetPoint(i,dblp); p[0]=dblp[0]; p[1]=dblp[1]; p[2]=dblp[2];
        closest->pts=new int[1];closest->weights=new P_float[1]; closest->dist2=1E10; closest->t=1E10; closest->nb=0; closest->crossing=false; closest->side=1;
        GetClosest(closest,p); resp[0]=p[0]-closest->c[0]; resp[1]=p[1]-closest->c[1]; resp[2]=p[2]-closest->c[2];  
        W=0; for(j=0;j<closest->nb;j++) W+=closest->weights[j]*closest->weights[j]; if(W!=0) W=1/W;
        for(j=0;j<closest->nb;j++) InsertConstraintPointForce(closest->pts[j],closest->weights[j]*W*resp[0],closest->weights[j]*W*resp[1],closest->weights[j]*W*resp[2]);
        free(closest->pts); free(closest->weights);
        }

for(i=0;i<nb_points;i++) if(constraintpointsforces[4*i]!=0) {constraintpointsforces[4*i+1]=constraintpointsforces[4*i+1]/constraintpointsforces[4*i]; constraintpointsforces[4*i+2]=constraintpointsforces[4*i+2]/constraintpointsforces[4*i];      constraintpointsforces[4*i+3]=constraintpointsforces[4*i+3]/constraintpointsforces[4*i];}
free(closest);

P_float alpha=Alpha_ConstraintPoints,extrapolation=ANISOTROPICEXTRAPOLATION;
for(int i=0;i<nb_points;i++)
	{
	Add_indep(i,constraintpointsforces+4*i+1,alpha,extrapolation);
//	memcpy(externalforces+3*i,constraintpointsforces+4*i+1,3*sizeof(P_float));
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::SetExternalForceRegularization(bool en, bool smooth, int transform,P_float lambda,P_float lambda_lim,P_float lambda_incr)
{
Regularization_Lambda=lambda;
Regularization_Lambda_lim=lambda_lim;
Regularization_Lambda_incr=lambda_incr;
Regularization_En=en;
Regularization_Smooth=smooth;
if(transform==0) {Regularization_Rigid=false; Regularization_Simi=false; Regularization_Affine=false;}
else if(transform==1) {Regularization_Rigid=true; Regularization_Simi=false; Regularization_Affine=false;}
else if(transform==2) {Regularization_Rigid=false; Regularization_Simi=true; Regularization_Affine=false;}
else if(transform==3) {Regularization_Rigid=false; Regularization_Simi=false; Regularization_Affine=true;}
}

void CSimplexSurf::UpdateForcesExternal() {
	
	int j;
	for (j = 0; j < 3 * nb_points; j++) externalforces[j] = 0;
	for (j = 0; j < nb_points; j++) externalforces_ign[j] = true;

	if (ExternalForce_IntensityProfile) UpdateIntensityProfileForces(Metric_IntensityProfile);	
	if (ExternalForce_GradientMRI) UpdateGradientMRIForces();
	if (ExternalForce_Demon) UpdateDemonForces();
	if (ExternalForce_ICP) UpdateForcesICP();
	//if(ExternalForce_ConstraintPoints)	UpdateForcesConstraintPoints();
	RegularizeExternalForces();	

	P_float alpha = 0,extrapolation = ANISOTROPICEXTRAPOLATION;
	if (ExternalForce_IntensityProfile) alpha += Alpha_IntensityProfile;
	if (ExternalForce_GradientMRI) alpha += Alpha_GradientMRI; 
	if (ExternalForce_Demon) alpha += Alpha_Demon; 
	if (ExternalForce_ICP) alpha += Alpha_ICP; 
	//if(ExternalForce_ConstraintPoints) alpha+=Alpha_ConstraintPoints; 

	//cout << "External forces" << endl;
	for (int i = 0; i < nb_points; i++){
		if (attached_point[i] == 0 && mass_constraint[i] == NULL) {// hack to prevent attached points / skin border from instbilities
			externalforces[3 * i] = externalforces[3 * i] / alpha; 
			externalforces[3 * i + 1] = externalforces[3 * i + 1] / alpha; 
			externalforces[3 * i + 2] = externalforces[3 * i + 2] / alpha;
			

			//cout << "  External force for point " << i << ": (" << externalforces[3 * i] << ", " << externalforces[3 * i + 1] << ", " << externalforces[3 * i + 2] << ")" << endl;
			Add_indep(i, externalforces + 3 * i, alpha, extrapolation);
		}
	}
	
}

void CSimplexSurf::UpdateForcesICP()
{
P_float u[3],alpha=Alpha_ICP,extrapolation=ANISOTROPICEXTRAPOLATION;
for(int i=0;i<nb_points;i++)
	if(attached_point[i]==0 && mass_constraint[i]==NULL)
		{
		GetClosestRefModel(GetPoint(i),GetNormal(i),u);
		externalforces[3*i]+=u[0];	externalforces[3*i+1]+=u[1];	externalforces[3*i+2]+=u[2];
		externalforces_ign[i]=false;
		}
}

void CSimplexSurf::UpdateIntensityProfileForces(int metric)
{
if(MRI_nb==0) return;

int i;
P_float *E=new P_float[nb_points*(2*Depth_IntensityProfile+1)]; for(i=0;i<nb_points*(2*Depth_IntensityProfile+1);i++) E[i]=1E10;
P_float tresh; if (metric==0 || metric==2) tresh=1E10; else tresh=-.2;
bool grad=false; if(metric==2 || metric==3) grad=true;
bool vec=false; if(grad) if(MRI_grad[0]->GetNumberOfScalarComponents()==3) vec=true;

ComputeIntensityProfile(IntensityProfile_mn_ref+Depth_IntensityProfile,IntensityProfile_pn_ref+Depth_IntensityProfile,IntensityProfile_s,grad); 
if(vec) GetIPEnergies(E,metric-2,nb_points,IntensityProfile_grad,IntensityProfileRef_grad,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
else if(grad) GetIPEnergies(E,metric-2,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);	
else GetIPEnergies(E,metric,nb_points,IntensityProfile,IntensityProfileRef,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref);

if(Regularization_En) RegularizeEnergies(E);
P_float *nml=normals; if(MRI_Type[0]==1) nml=MRI_Normals[0]; //rt
GetIPForces(externalforces,externalforces_ign,E,nb_points,IntensityProfile_mn,IntensityProfile_mn_ref,IntensityProfile_s,nml,Alpha_IntensityProfile,tresh);
}

void CSimplexSurf::UpdateDemonForces()
{
bool grad=false; if(MRI_grad[0]!=NULL) if(MRI_grad[0]->GetNumberOfScalarComponents()==1) grad=true;
if(!ExternalForce_IntensityProfile) ComputeIntensityProfile(IntensityProfile_mn_ref+(int)depth_Demon,IntensityProfile_pn_ref+(int)depth_Demon,IntensityProfile_s,grad); 
P_float *nml=normals; if(MRI_Type[0]==1) nml=MRI_Normals[0]; //rt
GetDemonForces(externalforces,externalforces_ign,nb_points,nml,IntensityProfile,IntensityProfileRef,IntensityProfile_s,IntensityProfile_mn,IntensityProfile_pn,IntensityProfile_mn_ref,IntensityProfile_pn_ref,Alpha_Demon);
}


// The depth in this function is set from the scripts. 
void CSimplexSurf::UpdateGradientMRIForces() {
	if (MRI_GradientMRI == NULL) return;

	P_float *E = new P_float[nb_points * (2 * depth_GradientMRI + 1)]; 
	
	for (int i = 0; i < nb_points * (2 * depth_GradientMRI + 1); i++) {
		E[i] = 1E10;
	}

	P_float tresh = 1E10;
	
	GetGradEnergies(E, nb_points, points, normals, MRI_GradientMRI, Size_GradientMRI, depth_GradientMRI, Interpolationmode, MRI_ROI_GradientMRI, oppositedirection_GradientMRI);
	
	if (MRI_GradientMRI2 != NULL) {
		GetGradEnergies(E, nb_points, points, normals, MRI_GradientMRI2, Size_GradientMRI, depth_GradientMRI, Interpolationmode, MRI_ROI_GradientMRI2, oppositedirection_GradientMRI);
	}
	if (Regularization_En) {
		RegularizeEnergies(E);
	}

	GetIPForces(externalforces, externalforces_ign, E, nb_points, depth_GradientMRI, Size_GradientMRI, normals, Alpha_GradientMRI, tresh);
}


void CSimplexSurf::RegularizeEnergies(P_float *E) {
	int i1, j1, j2, i, j, d = 2 * Depth_IntensityProfile + 1;
	P_float *E2 = new P_float[nb_points * d], w[3], dp1, dp2, dp3;
	
	for (i = 0; i < nb_points; i++) {
		dp1 = dotproduct(normals + 3 * i, normals + 3 * neighbors[3 * i]); 
		dp2 = dotproduct(normals + 3 * i, normals + 3 * neighbors[3 * i + 1]); 
		dp3 = dotproduct(normals + 3 * i, normals + 3 * neighbors[3 * i + 2]); 

		int *neighb = GetNeighbors(i);
		dp1 = dotproduct(normals + 3 * i, normals + 3 * neighb[0]); 
		dp2 = dotproduct(normals + 3 * i, normals + 3 * neighb[1]); 
		dp3 = dotproduct(normals + 3 * i, normals + 3 * neighb[2]); 
		
		for (j = -Depth_IntensityProfile; j <= Depth_IntensityProfile; j++) {
			i1 = j + Depth_IntensityProfile;
			E2[i * d + i1] = E[i * d + i1];
			
			int *Eneighb = GetNeighbors(i);

			w[0] = dp1 * j; 
			w[1] = w[0] - floor(w[0]); 
			w[2] = 1 - w[1]; 
			j1 = (int)floor(w[0]) + Depth_IntensityProfile; 
			j2 = (int)ceil(w[0]) + Depth_IntensityProfile;
			//E2[i * d + i1] += E[neighbors[3 * i + 0] * d + j1] * w[2] + E[neighbors[3 * i + 0] * d + j2] * w[1];
			E2[i * d + i1] += E[Eneighb[0] * d + j1] * w[2] + E[Eneighb[0] * d + j2] * w[1];

			w[0] = dp2 * j; 
			w[1] = w[0] - floor(w[0]); 
			w[2] = 1 - w[1]; 
			j1 = (int)floor(w[0]) + Depth_IntensityProfile; 
			j2 = (int)ceil(w[0]) + Depth_IntensityProfile; 		
			//E2[i * d + i1] += E[neighbors[3 * i + 1] * d + j1] * w[2] + E[neighbors[3 * i + 1] * d + j2] * w[1];
			E2[i * d + i1] += E[Eneighb[1] * d + j1] * w[2] + E[Eneighb[1] * d + j2] * w[1];

			w[0] = dp3 * j; 
			w[1] = w[0] - floor(w[0]); 
			w[2] = 1 - w[1]; 
			j1 = (int)floor(w[0]) + Depth_IntensityProfile; 
			j2 = (int)ceil(w[0]) + Depth_IntensityProfile; 		
			//E2[i * d + i1] += E[neighbors[3 * i + 2] * d + j1] * w[2] + E[neighbors[3 * i + 2] * d + j2] * w[1];
			E2[i * d + i1] += E[Eneighb[2] * d + j1] * w[2] + E[Eneighb[2] * d + j2] * w[1];
			
			E2[i * d + i1] = E2[i * d + i1] / 4.;
		}
	}
	memcpy(E, E2, nb_points * d * sizeof(P_float));
	free(E2);
}


void CSimplexSurf::RegularizeExternalForces() {
	int i, neighb[3];
	bool ign = true; 
	for (i = 0; i < nb_points; i++) 
		if (!externalforces_ign[i]) 
			ign = false; 
	if (ign) return;

	P_float n[3];
	if (Regularization_Smooth) {
		for (i = 0; i < nb_points; i++) {
			if (!externalforces_ign[i]) {
				n[0] = externalforces[3 * i]; 
				n[1] = externalforces[3 * i + 1]; 
				n[2] = externalforces[3 * i + 2];
				
				GetNeighbors(i, neighb);
				n[0] += externalforces[3 * neighb[0]];
				n[1] += externalforces[3 * neighb[0] + 1];
				n[2] += externalforces[3 * neighb[0] + 2];

				n[0] += externalforces[3 * neighb[1]];
				n[1] += externalforces[3 * neighb[1] + 1];
				n[2] += externalforces[3 * neighb[1] + 2];

				n[0] += externalforces[3 * neighb[2]]; 
				n[1] += externalforces[3 * neighb[2] + 1];
				n[2] += externalforces[3 * neighb[2] + 2];

				n[0] = n[0] / 4; 
				n[1] = n[1] / 4; 
				n[2] = n[2] / 4;
				
				memcpy(externalforces + 3 * i, n, 3 * sizeof(P_float));
			}
		}
	}

	P_float T[16];
	if (Regularization_Rigid) 
		Regularize_Rigid(T, nb_points, points, externalforces, externalforces_ign);
	else if (Regularization_Simi) 
		Regularize_Simi(T, nb_points, points, externalforces, externalforces_ign);
	else if (Regularization_Affine) 
		Regularize_Affine(T, nb_points, points, externalforces, externalforces_ign);

	P_float fg[3];
	if (Regularization_Rigid || Regularization_Affine || Regularization_Simi) {
		// fi= lfi + (1-l)(T-I)pi
		T[0] = T[0] - 1; 
		T[5] = T[5] - 1; 
		T[10] = T[10] - 1;
		for (i = 0; i < nb_points; i++) {
			Transform(points + 3 * i, fg, T);
			externalforces[3 * i] = Regularization_Lambda * externalforces[3 * i] + (1 - Regularization_Lambda) * fg[0]; 	
			externalforces[3 * i + 1] = Regularization_Lambda * externalforces[3 * i + 1] + (1 - Regularization_Lambda) * fg[1];		
			externalforces[3 * i + 2] = Regularization_Lambda * externalforces[3 * i + 2] + (1 - Regularization_Lambda) * fg[2];
		}
	}
}




/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::UpdateAlphas()
{
if(ExternalForce_IntensityProfile || ExternalForce_GradientMRI)
	{
	if((Regularization_Lambda<Regularization_Lambda_lim && Regularization_Lambda_incr>0) ||(Regularization_Lambda>Regularization_Lambda_lim && Regularization_Lambda_incr<0))  Regularization_Lambda+=Regularization_Lambda_incr;
	if(Depth_IntensityProfile_lim!=Depth_IntensityProfile_float)
		{
		if((Depth_IntensityProfile_float<Depth_IntensityProfile_lim && Depth_IntensityProfile_incr>0) ||(Depth_IntensityProfile_float>Depth_IntensityProfile_lim && Depth_IntensityProfile_incr<0))  Depth_IntensityProfile_float+=Depth_IntensityProfile_incr;
		Depth_IntensityProfile=(int)Depth_IntensityProfile_float;
		}
	}
if(ExternalForce_ConstraintPoints) if((Alpha_ConstraintPoints<Alpha_ConstraintPoints_lim && Alpha_ConstraintPoints_incr>0) ||(Alpha_ConstraintPoints>Alpha_ConstraintPoints_lim && Alpha_ConstraintPoints_incr<0))  Alpha_ConstraintPoints+=Alpha_ConstraintPoints_incr; 
if(ExternalForce_ConstraintMesh) if((Alpha_ConstraintMesh<Alpha_ConstraintMesh_lim && Alpha_ConstraintMesh_incr>0) ||(Alpha_ConstraintMesh>Alpha_ConstraintMesh_lim && Alpha_ConstraintMesh_incr<0))  Alpha_ConstraintMesh+=Alpha_ConstraintMesh_incr; 
if(ExternalForce_ICP) if((Alpha_ICP<Alpha_ICP_lim && Alpha_ICP_incr>0) ||(Alpha_ICP>Alpha_ICP_lim && Alpha_ICP_incr<0))  Alpha_ICP+=Alpha_ICP_incr; 
if(InternalForce_AxialConstraint) if((Alpha_AxialConstraint<Alpha_AxialConstraint_lim && Alpha_AxialConstraint_incr>0) ||(Alpha_AxialConstraint>Alpha_AxialConstraint_lim && Alpha_AxialConstraint_incr<0))  Alpha_AxialConstraint+=Alpha_AxialConstraint_incr; 
if(InternalForce_VolPres) if((Alpha_VolPres<Alpha_VolPres_lim && Alpha_VolPres_incr>0) ||(Alpha_VolPres>Alpha_VolPres_lim && Alpha_VolPres_incr<0))  Alpha_VolPres+=Alpha_VolPres_incr; 
if(InternalForce_SurfPres) if((Alpha_SurfPres<Alpha_SurfPres_lim && Alpha_SurfPres_incr>0) ||(Alpha_SurfPres>Alpha_SurfPres_lim && Alpha_SurfPres_incr<0))  Alpha_SurfPres+=Alpha_SurfPres_incr; 
if(InternalForce_Laplacian) if((Alpha_Laplacian<Alpha_Laplacian_lim && Alpha_Laplacian_incr>0) ||(Alpha_Laplacian>Alpha_Laplacian_lim && Alpha_Laplacian_incr<0))  Alpha_Laplacian+=Alpha_Laplacian_incr; 
if(InternalForce_LaplacianFlexion) if((Alpha_LaplacianFlexion<Alpha_LaplacianFlexion_lim && Alpha_LaplacianFlexion_incr>0) ||(Alpha_LaplacianFlexion>Alpha_LaplacianFlexion_lim && Alpha_LaplacianFlexion_incr<0))  Alpha_LaplacianFlexion+=Alpha_LaplacianFlexion_incr; 
if(InternalForce_SurfEqu) if((Alpha_SurfEqu<Alpha_SurfEqu_lim && Alpha_SurfEqu_incr>0) ||(Alpha_SurfEqu>Alpha_SurfEqu_lim && Alpha_SurfEqu_incr<0))  Alpha_SurfEqu+=Alpha_SurfEqu_incr; 
if(InternalForce_SurfEquFlexion) if((Alpha_SurfEquFlexion<Alpha_SurfEquFlexion_lim && Alpha_SurfEquFlexion_incr>0) ||(Alpha_SurfEquFlexion>Alpha_SurfEquFlexion_lim && Alpha_SurfEquFlexion_incr<0))  Alpha_SurfEquFlexion+=Alpha_SurfEquFlexion_incr; 
if(InternalForce_RefShape) if((Alpha_RefShape<Alpha_RefShape_lim && Alpha_RefShape_incr>0) ||(Alpha_RefShape>Alpha_RefShape_lim && Alpha_RefShape_incr<0))  Alpha_RefShape+=Alpha_RefShape_incr; 
if(ExternalForce_GradientMRI) if((Alpha_GradientMRI<Alpha_GradientMRI_lim && Alpha_GradientMRI_incr>0) ||(Alpha_GradientMRI>Alpha_GradientMRI_lim && Alpha_GradientMRI_incr<0))  Alpha_GradientMRI+=Alpha_GradientMRI_incr; 
if(ExternalForce_Demon) if((Alpha_Demon<Alpha_Demon_lim && Alpha_Demon_incr>0) ||(Alpha_Demon>Alpha_Demon_lim && Alpha_Demon_incr<0))  Alpha_Demon+=Alpha_Demon_incr; 
if(ExternalForce_IntensityProfile) if((Alpha_IntensityProfile<Alpha_IntensityProfile_lim && Alpha_IntensityProfile_incr>0) ||(Alpha_IntensityProfile>Alpha_IntensityProfile_lim && Alpha_IntensityProfile_incr<0))  Alpha_IntensityProfile+=Alpha_IntensityProfile_incr; 
}


void CSimplexSurf::UpdateExternalForces() {
	
	int j;
	UpdateAll();
	UpdateAlphas();
	for (j = 0; j < 6 * nb_points; j++) df_p[j] = 0; 
	for (j = 0; j < 3 * nb_points; j++) forces[j] = 0; 

	//if(ExternalForce_IntensityProfile || ExternalForce_GradientMRI  || ExternalForce_Demon || ExternalForce_ICP || ExternalForce_ConstraintPoints) UpdateForcesExternal();
	if (ExternalForce_IntensityProfile || ExternalForce_GradientMRI  || ExternalForce_Demon || ExternalForce_ICP) UpdateForcesExternal();
	if (ExternalForce_ConstraintPoints)	UpdateForcesConstraintPoints();
	if (ExternalForce_ConstraintMesh)	UpdateForcesConstraintMesh();
	if (Axis_model != NULL && IsAxis == true) if (InternalForce_AxialConstraint || Axis_model->GetInternalForce_AxialConstraint()) UpdateAxialLinks_forces(); 
	
}


void CSimplexSurf::UpdateInternalForces() {
	
	if (InternalForce_AxialConstraint)	{
		UpdateForcesAxialConstraint();
	}
	if (InternalForce_VolPres) {
		UpdateForcesVolPres();
	}
	if (InternalForce_SurfPres) {
		UpdateForcesSurfPres();
	}
	if (InternalForce_Laplacian)	{
		UpdateForcesLaplacian();
	}
	if (InternalForce_LaplacianFlexion) {
		UpdateForcesLaplacianFlexion();
	}
	if (InternalForce_SurfEqu) {
		UpdateForcesSurfEqu();
	}
	if (InternalForce_SurfEquFlexion) {
		UpdateForcesSurfEquFlexion();
	}
	if (InternalForce_RefShape) {
		UpdateForcesRefShape();
	}
	
	//verifyPointsAndParams();
	

	//mergedpoints
	int j;
	if (MergedPoints != NULL) {
		if (USEFORCES) { 
			for (j = 0; j < MergedPoints[0]; j++) {
				AddForce(MergedPoints[2 * j + 1], GetForce(MergedPoints[2 * j + 2]));		
				AddDfp(MergedPoints[2 * j + 1], GetDfp(MergedPoints[2 * j + 2]));
			}
		}
		else {
			for (j = 0; j < MergedPoints[0]; j++) {
				P_float p[3], p2[3]; 
				GetPoint(MergedPoints[2 * j + 1], p); 
				GetPoint(MergedPoints[2 * j + 2], p2);	
				
				p[0] = (p[0] + p2[0]) / 2.; 
				p[1] = (p[1] + p2[1]) / 2.; 
				p[2] = (p[2] + p2[2]) / 2.;

				SetPoint(MergedPoints[2 * j + 1], p); 
				SetPoint(MergedPoints[2 * j + 2], p);	
			}
		}
	}
	
}

void CSimplexSurf::Add_multi(int index,P_float *pref,P_float h,P_float alpha,P_float extrapolation) {if(USEFORCES) AddForce_multi(index,pref,h,alpha,extrapolation); else AddDisp_multi(index,pref,h,alpha); }

void CSimplexSurf::AddForce_multi(int index,P_float *pref,P_float h,P_float alpha,P_float extrapolation)
{
int k,*ne=neighbors+3*index;
P_float F[3],F1[3],F2[3],F3[3],DFP[6],DFP1[6],DFP2[6],DFP3[6],lambda[4]; 
GetForces_Shape(GetPoint(index),GetPoint(ne[0]),GetPoint(ne[1]),GetPoint(ne[2]),GetNormal(index),pref,h,F,F1,F2,F3,lambda);

for(k=0;k<3;k++) {F[k]=F[k]*alpha; F1[k]=F1[k]*alpha; F2[k]=F2[k]*alpha; F3[k]=F3[k]*alpha;}
GetDerivatives(extrapolation,F,DFP); for(k=0;k<6;k++) {DFP[k]=DFP[k]*alpha; DFP1[k]=pref[0]*lambda[1]*DFP[k]; DFP2[k]=pref[1]*lambda[2]*DFP[k];  DFP3[k]=pref[2]*lambda[3]*DFP[k]; }

if(GetMass_constraint(index)!=NULL)
	{
	// to complete
	}
else if(mass_modified)
	{
	// to complete
	}
else
	{
	P_float M=pref[0]*lambda[1]*mass_inv[3*ne[0]] + pref[1]*lambda[2]*mass_inv[3*ne[1]] + pref[2]*lambda[3]*mass_inv[3*ne[2]] + mass_inv[3*index];
	if(M!=0) M=1/M;

	for(k=0;k<3;k++) {F[k]=F[k]*M; F1[k]=F1[k]*M; F2[k]=F2[k]*M; F3[k]=F3[k]*M;}
	for(k=0;k<6;k++) {DFP[k]=DFP[k]*M; DFP1[k]=DFP1[k]*M; DFP2[k]=DFP2[k]*M; DFP3[k]=DFP3[k]*M;}
	}

for(k=0;k<3;k++) {F[k]=F[k]*timestep_inv2; F1[k]=F1[k]*timestep_inv2; F2[k]=F2[k]*timestep_inv2; F3[k]=F3[k]*timestep_inv2;}
for(k=0;k<6;k++) {DFP[k]=DFP[k]*timestep_inv2; DFP1[k]=DFP1[k]*timestep_inv2; DFP2[k]=DFP2[k]*timestep_inv2; DFP3[k]=DFP3[k]*timestep_inv2;}

AddForce(index,F); AddForce(ne[0],F1); AddForce(ne[1],F2); AddForce(ne[2],F3);
AddDfp(index,DFP); AddDfp(ne[0],DFP1); AddDfp(ne[1],DFP2); AddDfp(ne[2],DFP3);
}

void CSimplexSurf::AddDisp_multi(int index,P_float *pref,P_float h,P_float alpha)
{
int k,*ne=neighbors+3*index;
P_float F[3],F1[3],F2[3],F3[3],lambda[4]; GetForces_Shape(GetPoint(index),GetPoint(ne[0]),GetPoint(ne[1]),GetPoint(ne[2]),GetNormal(index),pref,h,F,F1,F2,F3,lambda);

for(k=0;k<3;k++) {F[k]=F[k]*alpha; F1[k]=F1[k]*alpha; F2[k]=F2[k]*alpha; F3[k]=F3[k]*alpha;}

if(GetMass_constraint(index)!=NULL)
	{
	// to complete
	}
else if(mass_modified)
	{
	// to complete
	}
else
	{
	P_float M=pref[0]*lambda[1]*mass_inv[3*ne[0]] + pref[1]*lambda[2]*mass_inv[3*ne[1]] + pref[2]*lambda[3]*mass_inv[3*ne[2]] + mass_inv[3*index];
	if(M!=0) M=1/M;
	for(k=0;k<3;k++) {F[k]=F[k]*mass_inv[3*index]*M; F1[k]=F1[k]*mass_inv[3*ne[0]]*M; F2[k]=F2[k]*mass_inv[3*ne[1]]*M; F3[k]=F3[k]*mass_inv[3*ne[2]]*M;}
	}

if(GetAttachedPoint(index)==0) AddDisp(index,F); 
if(GetAttachedPoint(ne[0])==0) AddDisp(ne[0],F1); 
if(GetAttachedPoint(ne[1])==0) AddDisp(ne[1],F2); 
if(GetAttachedPoint(ne[2])==0) AddDisp(ne[2],F3);
}

void CSimplexSurf::AddDisp(int index,P_float *F) {
	P_float n[3], cp[3], *P, *P1, *P2, *P3, P1P3[3], P2P3[3], P1P2[3], PP1[3], s, dp;
	int i, *ne = GetNeighbors(index), *nc = GetNeighbors_c(index);
	/*
	P_float v;
	int j,i2,i3,*c;

	// before update, substract surface/ volume
	for(i=0;i<3;i++) 
		if(nc[i]!=-1)
			{ 
			surfaces_cell[nc[i]]=0; P1=cellcenters+3*nc[i];	c=cells[nc[i]];
			// for each face calcul volume and surface
			for(j=0;j<c[0];j++) 
				{
				if(Flipped) {i3=c[j+1]; i2=c[(j==c[0]-1)?1:(j+2)];} else {i2=c[j+1]; i3=c[(j==c[0]-1)?1:(j+2)];}
				P2=GetPoint(i2); P3=GetPoint(i3); P1P3[0]=P3[0]-P1[0]; P1P3[1]=P3[1]-P1[1]; P1P3[2]=P3[2]-P1[2];  P1P2[0]=P2[0]-P1[0]; P1P2[1]=P2[1]-P1[1]; P1P2[2]=P2[2]-P1[2];  
				crossproduct(cp,P1P2,P1P3); s=norm(cp)/2.;
				v= cp[2] * (P1[2]+P2[2]+P3[2]) /6.; 
				volume-=v; dv-=v; surface-=s; surfaces_cell[nc[i]]-=s; surfaces[i2]-=s/2.; surfaces[i3]-=s/2.; 
				}
			}*/
	
	// update point and ext forces
	GetPoint(index)[0] += F[0]; GetPoint(index)[1] += F[1]; GetPoint(index)[2] += F[2];
	if(axiallinks!=NULL) {axiallinks[index].f[0]-=F[0]; axiallinks[index].f[1]-=F[1]; axiallinks[index].f[2]-=F[2];}
	externalforces[3*index]-=F[0]; externalforces[3*index+1]-=F[1]; externalforces[3*index+2]-=F[2];
	constraintpointsforces[4*index+1]-=F[0]; constraintpointsforces[4*index+2]-=F[1]; constraintpointsforces[4*index+3]-=F[2];
	
	// update h
	P=GetPoint(index); P1=GetPoint(ne[0]); PP1[0]=P1[0]-P[0];PP1[1]=P1[1]-P[1];PP1[2]=P1[2]-P[2];
	if(nc[2]==-1) h[index]=dotproduct(PP1,normals2+3*index); 	else h[index]=dotproduct(PP1,normals+3*index); 
	// Update normals
	for(i=0;i<3;i++) {
		P=GetPoint(ne[i]); 
		//P1=GetPoint(neighbors[3*ne[i]]); 
		//P2=GetPoint(neighbors[3*ne[i]+1]);  
		//P3=GetPoint(neighbors[3*ne[i]+2]);    

		int *neig = GetNeighbors(ne[i]);
		P1=GetPoint(neig[0]); 
		P2=GetPoint(neig[1]);  
		P3=GetPoint(neig[2]);    
		
		P1P3[0]=P3[0]-P1[0];P1P3[1]=P3[1]-P1[1];P1P3[2]=P3[2]-P1[2];    P2P3[0]=P3[0]-P2[0];P2P3[1]=P3[1]-P2[1];P2P3[2]=P3[2]-P2[2];    PP1[0]=P1[0]-P[0];PP1[1]=P1[1]-P[1];PP1[2]=P1[2]-P[2];
		// normals2
		if(neighbors_c[3*ne[i]+2]==-1) 	{P1P2[0]=P2[0]-P1[0];P1P2[1]=P2[1]-P1[1];P1P2[2]=P2[2]-P1[2]; s=dotproduct(P1P2,P1P2); dp=dotproduct(PP1,P1P2); n[0]=-PP1[0]+dp*P1P2[0]/s; n[1]=-PP1[1]+dp*P1P2[1]/s; n[2]=-PP1[2]+dp*P1P2[2]/s;	s=norm(n);  normals2[3*ne[i]]=n[0]/s; normals2[3*ne[i]+1]=n[1]/s; normals2[3*ne[i]+2]=n[2]/s;  }
		// normals
		crossproduct(n,P2P3,P1P3); 	s=norm(n); 	if(s==0) {crossproduct(cp,P1P3,PP1); crossproduct(n,P1P3,cp);	s=norm(n); if(s==0) SetNormal(ne[i],1,0,0); else SetNormal(ne[i],n[0]/s,n[1]/s,n[2]/s);}	else SetNormal(ne[i],n[0]/s,n[1]/s,n[2]/s);	
		// h
		if(neighbors_c[3*ne[i]+2]==-1) h[ne[i]]=dotproduct(PP1,normals2+3*ne[i]); 	else h[ne[i]]=dotproduct(PP1,normals+3*ne[i]); 
		}
	// update cell centers
	cellcenters[3*nc[0]]+=F[0]/(P_float)cells[nc[0]][0]; cellcenters[3*nc[0]+1]+=F[1]/(P_float)cells[nc[0]][0]; cellcenters[3*nc[0]+2]+=F[2]/(P_float)cells[nc[0]][0];
	cellcenters[3*nc[1]]+=F[0]/(P_float)cells[nc[1]][0]; cellcenters[3*nc[1]+1]+=F[1]/(P_float)cells[nc[1]][0]; cellcenters[3*nc[1]+2]+=F[2]/(P_float)cells[nc[1]][0];
	if(nc[2]!=-1) {cellcenters[3*nc[2]]+=F[0]/(P_float)cells[nc[2]][0]; cellcenters[3*nc[2]+1]+=F[1]/(P_float)cells[nc[2]][0]; cellcenters[3*nc[2]+2]+=F[2]/(P_float)cells[nc[2]][0];}
	/*// update surfaces/ volume
	for(i=0;i<3;i++) 
		if(nc[i]!=-1)
			{ 
			surfaces_cell[nc[i]]=0; P1=cellcenters+3*nc[i];	c=cells[nc[i]];
			// for each face calcul volume and surface
			for(j=0;j<c[0];j++) 
				{
				if(Flipped) {i3=c[j+1]; i2=c[(j==c[0]-1)?1:(j+2)];} else {i2=c[j+1]; i3=c[(j==c[0]-1)?1:(j+2)];}
				P2=GetPoint(i2); P3=GetPoint(i3); P1P3[0]=P3[0]-P1[0]; P1P3[1]=P3[1]-P1[1]; P1P3[2]=P3[2]-P1[2];  P1P2[0]=P2[0]-P1[0]; P1P2[1]=P2[1]-P1[1]; P1P2[2]=P2[2]-P1[2];  
				crossproduct(cp,P1P2,P1P3); s=norm(cp)/2.;
				v= cp[2] * (P1[2]+P2[2]+P3[2]) /6.; 
				volume+=v; dv+=v; surface+=s; surfaces_cell[nc[i]]+=s; surfaces[i2]+=s/2.; surfaces[i3]+=s/2.; 
				}
			}*/
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::TimeStep()
{
int i;

TimeCounter++;

// mergedpoints
if(MergedPoints!=NULL) if(USEFORCES) for(i=0;i<MergedPoints[0];i++) SetPoint(MergedPoints[2*i+2],GetPoint(MergedPoints[2*i+1]));		

// Topo
if(Topo!=0)
    {
	//CPerfTimer start; start.Start();

	bool flag=false;
	if(Topo_Elongation[0]!=0) if(fmod((double)TimeCounter+Topo_Increase_L[0]/4,Topo_Elongation[0])==0) {Subdivide_Elongation(); flag=true;}
	if(Topo_Increase_L[0]!=0) if(fmod((double)TimeCounter+Topo_Increase_L[0]/2,Topo_Increase_L[0])==0) {Subdivide_Increase_L(); flag=true;}
	if(Topo_Decrease_L[0]!=0) if(fmod((double)TimeCounter,Topo_Decrease_L[0])==0) {Subdivide_Decrease_L(); flag=true;}
	if(Topo_Exchange[0]!=0) if(fmod((double)TimeCounter,Topo_Exchange[0])==0) {Subdivide_Exchange(); flag=true;}
//	if(Topo_Increase_R[0]!=0) if(fmod((double)TimeCounter+3,Topo_Increase_R[0])==0) {Subdivide_Increase_R(); flag=true;}
//	if(Topo_Decrease_R[0]!=0) if(fmod((double)TimeCounter+3+Topo_Increase_R[0]/2,Topo_Decrease_R[0])==0) {Subdivide_Decrease_R(); flag=true;}
//	if(Topo_Increase_Err[0]!=0) if(fmod((double)TimeCounter+6,Topo_Increase_Err[0])==0) {Subdivide_Increase_Err(); flag=true;}
//	if(Topo_Decrease_Err[0]!=0) if(fmod((double)TimeCounter+6+Topo_Increase_Err[0]/2,Topo_Decrease_Err[0])==0) {Subdivide_Decrease_Err(); flag=true;}

	//start.Stop(); Monitor_TOduration+= start.Elapsedms(); 
	if(flag) 
		{
		BuildDOPCellTree(-1,0);
	    if(Monitor_surfaces) Update_Monitor_surfaces();
		if(Monitor_edgelenght) Update_Monitor_edgelenght();
		}
	}

// Monitoring
if(Monitor)
    {
    if(Monitor_axiallinks_error) Update_Monitor_axiallinks_error();
    if(Monitor_density) Update_Monitor_density();
    if(Monitor_displacement) Update_Monitor_displacement();
    }

}

void CSimplexSurf::Equilibrium()
{
for(int i=0; i<3*nb_points;i++) {speeds[i]=0; speeds_tm1[i]=0; forces[i]=0; points_tm1[i]=points[i];}
for(int i=0; i<6*nb_points;i++) {df_p[i]=0; df_s[i]=0;}
dv=0;
UpdateAll();
FillVal0();

BuildDOPCellTree(-1,0);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::SetMonitor_axiallinks_error(bool val,const char* filename) {Monitor_axiallinks_error=val; Monitor_axiallinks_error_filename=filename; Update_Monitor();}
void CSimplexSurf::SetMonitor_surfaces(bool val,const char* filename) {Monitor_surfaces=val; Monitor_surfaces_filename=filename; Update_Monitor();}
void CSimplexSurf::Update_Monitor() {Monitor=Monitor_axiallinks_error || Monitor_surfaces || Monitor_edgelenght || Monitor_density || Monitor_displacement;}

void CSimplexSurf::Update_Monitor_axiallinks_error()
{
FILE* f=fopen(Monitor_axiallinks_error_filename.c_str(),"rt+");
if(f==NULL) f=fopen(Monitor_axiallinks_error_filename.c_str(),"wt"); 
else fseek(f,0,SEEK_END);
P_float *err=ComputeAxialLinks_error();
fprintf(f,"%lf ( %lf )\n",err[0],err[1]);
free(err);
fclose(f);
}

void CSimplexSurf::Update_Monitor_surfaces()
{
FILE* f=fopen(Monitor_surfaces_filename.c_str(),"rt+");
if(f==NULL) f=fopen(Monitor_surfaces_filename.c_str(),"wt"); 
else fseek(f,0,SEEK_END);
/*  // HISTO
int tolerance=1,max=200,i;
if(Monitor_surfaces==1) // initialisation
    {
    for(i=0;i<max/tolerance;i++) fprintf(f,"%d ",i*tolerance);
    Monitor_surfaces=2;
    }
fprintf(f,"\n");
int* histo=new int[max/tolerance]; for(i=0;i<max/tolerance;i++) histo[i]=0;
for(i=0;i<nb_points;i++) histo[(int)floor(surfaces[i]/(P_float)tolerance)]++;
for(i=0;i<max/tolerance;i++) fprintf(f,"%d ",histo[i]);
free(histo);
*/
// MEAN / std dev
P_float mean=0,stddev=0;
//int i; for(i=0;i<nb_points;i++) {mean+=Topo_sref-surfaces[i]; stddev+=(Topo_sref-surfaces[i])*(Topo_sref-surfaces[i]);} // difference with target
int i; for(i=0;i<nb_points;i++) {mean+=surfaces[i]; stddev+=(surfaces[i])*(surfaces[i]);} // absolute
mean=mean/(P_float)(nb_points);
stddev=sqrt(stddev/(P_float)nb_points-mean*mean);
fprintf(f,"%lf %lf %lf\n",mean,stddev,Monitor_TOduration);
fclose(f);
}


void CSimplexSurf::Update_Monitor_edgelenght()
{
FILE* f=fopen(Monitor_edgelenght_filename.c_str(),"rt+");
if(f==NULL) f=fopen(Monitor_edgelenght_filename.c_str(),"wt"); 
else fseek(f,0,SEEK_END);
// MEAN / std dev
P_float mean=0,stddev=0,d;
int i,j,count=0;
for(i=0;i<nb_points;i++) 
	for(j=0;j<3;j++) 
		if(neighbors[3*i+j]>i)
			{
			//d=Topo_lref-dist3D(GetPoint(i),GetPoint(neighbors[3*i+j]));  // difference with target
			d=dist3D(GetPoint(i),GetPoint(neighbors[3*i+j])); // absolute
			mean+=d; stddev+=d*d;
			count++;
			}
mean=mean/((P_float)(count));
stddev=sqrt(stddev/((P_float)(count))-mean*mean);
fprintf(f,"%lf %lf %lf\n",mean,stddev,Monitor_TOduration);
fclose(f);
}

void CSimplexSurf::Update_Monitor_density()
{
FILE* f=fopen(Monitor_density_filename.c_str(),"rt+");
if(f==NULL) f=fopen(Monitor_density_filename.c_str(),"wt"); 
else fseek(f,0,SEEK_END);
P_float surf=0; for(int i=0;i<nb_points;i++) surf+=surfaces[i];
fprintf(f,"%lf\n",surf/(P_float)nb_cells);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::SetTopo_Increase_L(int frequ, int nb, bool ignoreborder) {
	Topo_ignoreborder = ignoreborder; 
	Topo_Increase_L[0] = frequ; 
	Topo_Increase_L[1] = nb;  
	Update_Topo();
}

void CSimplexSurf::SetTopo_Decrease_L(int frequ, int nb, bool ignoreborder) {
	Topo_ignoreborder = ignoreborder; 
	Topo_Decrease_L[0] = frequ; 
	Topo_Decrease_L[1] = nb; 
	Update_Topo();
}

void CSimplexSurf::SetTopo_Elongation(int frequ, int nb, bool ignoreborder) {
	Topo_ignoreborder = ignoreborder; 
	Topo_Elongation[0] = frequ; 
	Topo_Elongation[1] = nb; 
	Update_Topo();
}

void CSimplexSurf::SetTopo_Exchange(int frequ, int nb, bool ignoreborder) {
	Topo_ignoreborder = ignoreborder; 
	Topo_Exchange[0] = frequ; 
	Topo_Exchange[1] = nb; 
	Update_Topo();
}

void CSimplexSurf::Update_Topo() {
	Topo = Topo_Increase_L[0] + Topo_Decrease_L[0] + Topo_Elongation[0] + Topo_Exchange[0];
}

bool CSimplexSurf::GetTopo() {
	if (Topo == 0) return false; 
	else return true;
}

void CSimplexSurf::SetTopo_Tol(P_float tol) {
	Topo_Tol = tol;
}

void CSimplexSurf::SetTopo_Lref(P_float lref) {
	Topo_lref = lref; 
	Topo_sref = 3. * sqrt((double)3) * lref * lref / 4.;
	Topo_scellref = 2. * Topo_sref;
	Topo_nref = (int)((double)surface / (double)Topo_sref);
}

void CSimplexSurf::Subdivide_Update() {
	Topo_nref = (int)((double)surface / (double)Topo_sref); 
	for (int i = 0; i < nb_points; i++) 
		UpdateNeighbors2(i); 
	
	OrderNeighbors(); 
	UpdateAll();
}

void CSimplexSurf::Subdivide_Exchange() {
	int i, j, count = 0, countmax = (Topo_Exchange[1] == -1) ? (2.5 * nb_points):(Topo_Exchange[1]);
	while (count < countmax) {
		i = rand() * nb_points / RAND_MAX; 
		if (i == nb_points) 
			i = nb_points - 1; 
		
		if (TOech2_test(i)) {
			TOech2(i); 
			Subdivide_Update();
		}    
		
		count++;
		i = rand() * nb_points / RAND_MAX; 
		if (i == nb_points) 
			i = nb_points - 1; 
		
		j = rand() * 3 / RAND_MAX; 
		if (j == 3) 
			j = 2;	
		
		//if (TOech_test(i, neighbors[3 * i + j])) {
		//	TOech(i, neighbors[3 * i + j]);
		//	Subdivide_Update();
		//}    
		int *neigh = GetNeighbors(i);
		if (TOech_test(i, neigh[j])) {
			TOech(i, neigh[j]);
			Subdivide_Update();
		}
		
		count++;
	}
}


void CSimplexSurf::Subdivide_Increase_L() { // increase res with TO2 at maximum surface_cell 
	int i, j, count = 0, countmax, pt[4]; 
	P_float l[2];

	if (Topo_Increase_L[1] == -1) 
		countmax = nb_points; 
	else if (Topo_Increase_L[1] == -2) {
		countmax = (Topo_nref - nb_points) / 2; 
		if (countmax <= 0) return;
	} 
	else countmax = Topo_Increase_L[1];
	
	while (count < countmax) {
		j = rand() * nb_cells / RAND_MAX; 

		if (j == nb_cells) 
			j = nb_cells - 1; 
		
		if (TO2_test(j, 0, pt, l)) {
			TO2(j, pt[0], pt[1], pt[2], pt[3]);
			Subdivide_Update();	
		} 
		
		count++;
		
		i = rand() * nb_points / RAND_MAX; 
		if (i == nb_points) 
			i = nb_points - 1; 
		
		if (TOvertex_i_test(i)) {
			TOvertex_i(i); 	
			Subdivide_Update();	
		}
		
		count++;
		
		j = rand() * nb_cells / RAND_MAX; 
		if (j == nb_cells) 
			j = nb_cells - 1; 
		
		if (TOcell_i_test(j)) {
			TOcell_i(j); 	
			Subdivide_Update();	
		}
		
		count++;
		count++;
	}
}


void CSimplexSurf::Subdivide_Elongation() { // increase res with TO2 if 1st_principal_length>2.2nd_principal_length
	int j, count = 0, countmax, pt[4]; 
	P_float l[2];
	
	if (Topo_Elongation[1] == -1) 
		countmax = nb_points; 
	else 
		countmax = Topo_Elongation[1];

	while (count < countmax) {
		j = rand() * nb_cells / RAND_MAX; 
		
		if (j == nb_cells) 
			j = nb_cells - 1;
		if (TO2_test(j, sqrt(2.) * Topo_Tol, pt, l)) {
			TO2(j, pt[0], pt[1], pt[2], pt[3]); 
			Subdivide_Update();
		}
		count++;   
	}
}


void CSimplexSurf::Subdivide_Decrease_L() { // decrease res with TO1 at minimum joined surface_cells 
	bool *first = new bool[1];
	int i, j, k, count = 0, countmax;
	
	if (Topo_Decrease_L[1] == -1) 
		countmax = nb_points; 
	else if (Topo_Decrease_L[1] == -2) {
		countmax = (nb_points - Topo_nref) / 2; 
		if (countmax <= 0) return;
	} 
	else countmax = Topo_Decrease_L[1];
	
	while (count < countmax) {
		i = rand() * nb_points / RAND_MAX; 
		if (i == nb_points) 
			i = nb_points - 1; 
		
		j = rand() * 3 / RAND_MAX; 
		if (j == 3) 
			j = 2; 
		
		j = neighbors[3 * i + j]; 	
		if (TO1_test(i, j)) {
			TO1(i, j); 
			Subdivide_Update(); 
		} 
		
		count++;
		
		k = rand() * nb_cells / RAND_MAX;	
		if (k == nb_cells) 
			k = nb_cells - 1; 
		
		if (TOcell_d_test(k)) {
			TOcell_d(k); 
			Subdivide_Update(); 
		} 
		
		count++;
		
		k = rand() * nb_cells / RAND_MAX; 
		if (k == nb_cells) 
			k = nb_cells - 1; 
		
		if (TOvertex_d_test(k, first)) {
			TOvertex_d(k, *first); 
			Subdivide_Update(); 
		} 
		count++;
		count++;
	}
	free(first);
}

//
//void CSimplexSurf::Subdivide_Elongation()
//{
//// split cells/ elongation
//int i1,i2,i3,i4,ign1=0,ign2,test; P_float *l,e[3],p1[3],p2[3],p3[3],p4[3],lmax,lmax2;
//int i,count=0,countmax=(Topo_Elongation[1]==-1)?(nb_cells):(Topo_Elongation[1]);
//
//while(count!=countmax)
//    {
//    i=rand()*nb_cells/RAND_MAX;
//    l=new P_float[cells[i][0]];
//    lmax=0; lmax2=0;
//    for(int k=1;k<=cells[i][0];k++) 
//        {
//        i1=cells[i][k]; i2=cells[i][(k==cells[i][0])?1:(k+1)];e[0]=points[3*i1]-points[3*i2]; e[1]=points[3*i1+1]-points[3*i2+1];e[2]=points[3*i1+2]-points[3*i2+2]; l[k-1]=norm(e); 
//        if(l[k-1]>lmax) {lmax2=lmax; ign2=ign1; lmax=l[k-1]; ign1=k-1;}
//        else if(l[k-1]>lmax2) {lmax2=l[k-1]; ign2=k-1;}
//        }
//    i1=cells[i][ign1+1];
//    i2=cells[i][(ign1==(cells[i][0]-1))?1:(ign1+2)];
//    i3=cells[i][ign2+1];
//    i4=cells[i][(ign2==(cells[i][0]-1))?1:(ign2+2)];
//
//	// abs(l11'/2-lref)+abs(l12'-lref)+abs(l21'/2-lref)+abs(l22'-lref)+abs(ln-lref) <? abs(l1-lref) + abs(l2-lref)
//	GetPoint(i1,p1);GetPoint(i2,p2);GetPoint(i3,p3);GetPoint(i4,p4); 
//	P_float A[3]={(0.5*p3[0]+0.5*p4[0]+p1[0]+p2[0])/3,(0.5*p3[1]+0.5*p4[1]+p1[1]+p2[1])/3,(0.5*p3[2]+0.5*p4[2]+p1[2]+p2[2])/3};
//	P_float B[3]={(0.5*p1[0]+0.5*p2[0]+p3[0]+p4[0])/3,(0.5*p1[1]+0.5*p2[1]+p3[1]+p4[1])/3,(0.5*p1[2]+0.5*p2[2]+p3[2]+p4[2])/3};
//	e[0]=(p1[0]-A[0]); e[1]=(p1[1]-A[1]); e[2]=(p1[2]-A[2]); P_float l11=norm(e);
//	e[0]=(p2[0]-A[0]); e[1]=(p2[1]-A[1]); e[2]=(p2[2]-A[2]); P_float l12=norm(e);
//	e[0]=(p3[0]-B[0]); e[1]=(p3[1]-B[1]); e[2]=(p3[2]-B[2]); P_float l21=norm(e);
//	e[0]=(p4[0]-B[0]); e[1]=(p4[1]-B[1]); e[2]=(p4[2]-B[2]); P_float l22=norm(e);
//	e[0]=(A[0]-B[0]); e[1]=(A[1]-B[1]); e[2]=(A[2]-B[2]); P_float ln=norm(e);
//	test= (abs(l[ign1]-Topo_lref) + abs(l[ign2]-Topo_lref))/2 - (abs(l11-Topo_lref) + abs(l12-Topo_lref) + abs(l21-Topo_lref) + abs(l22-Topo_lref) + abs(ln-Topo_lref))/5;
//	if(test>0 && abs(ign2-ign1)>1 && abs(ign2-ign1)!=(cells[i][0]-1))   TO2(i,i1,i2,i3,i4);
//    free(l);
//    count++;
//	}
//}

/*
void CSimplexSurf::Subdivide_Elongation()
{
// split cells/ elongation
int i,i1,i2,i3,i4,imax,imax1,imax2,imax3,imax4,ign1=0,ign2,count=0,countmax=(Topo_Elongation[1]==-1)?(nb_cells):(Topo_Elongation[1]);
P_float *l,e[3],p1[3],p2[3],p3[3],p4[3],lmax,lmax2,test,testmax;

while(count!=countmax)
    {
	testmax=-1E10;

	i=rand()*nb_cells/RAND_MAX; //for(i=0;i<nb_cells;i++)
		{
		l=new P_float[cells[i][0]];
		lmax=0; lmax2=0;
		for(int k=1;k<=cells[i][0];k++) 
			{
			i1=cells[i][k]; i2=cells[i][(k==cells[i][0])?1:(k+1)];e[0]=points[3*i1]-points[3*i2]; e[1]=points[3*i1+1]-points[3*i2+1];e[2]=points[3*i1+2]-points[3*i2+2]; l[k-1]=norm(e); 
			if(l[k-1]>lmax) {lmax2=lmax; ign2=ign1; lmax=l[k-1]; ign1=k-1;}
			else if(l[k-1]>lmax2) {lmax2=l[k-1]; ign2=k-1;}
			}
		if(abs(ign2-ign1)>1 && abs(ign2-ign1)!=(cells[i][0]-1))
			{
			i1=cells[i][ign1+1];
			i2=cells[i][(ign1==(cells[i][0]-1))?1:(ign1+2)];
			i3=cells[i][ign2+1];
			i4=cells[i][(ign2==(cells[i][0]-1))?1:(ign2+2)];

			// abs(l11'/2-lref)+abs(l12'-lref)+abs(l21'/2-lref)+abs(l22'-lref)+abs(ln-lref) <? abs(l1-lref) + abs(l2-lref)
			GetPoint(i1,p1);GetPoint(i2,p2);GetPoint(i3,p3);GetPoint(i4,p4); 
			P_float A[3]={(0.5*p3[0]+0.5*p4[0]+p1[0]+p2[0])/3,(0.5*p3[1]+0.5*p4[1]+p1[1]+p2[1])/3,(0.5*p3[2]+0.5*p4[2]+p1[2]+p2[2])/3};
			P_float B[3]={(0.5*p1[0]+0.5*p2[0]+p3[0]+p4[0])/3,(0.5*p1[1]+0.5*p2[1]+p3[1]+p4[1])/3,(0.5*p1[2]+0.5*p2[2]+p3[2]+p4[2])/3};
			e[0]=(p1[0]-A[0]); e[1]=(p1[1]-A[1]); e[2]=(p1[2]-A[2]); P_float l11=norm(e);
			e[0]=(p2[0]-A[0]); e[1]=(p2[1]-A[1]); e[2]=(p2[2]-A[2]); P_float l12=norm(e);
			e[0]=(p3[0]-B[0]); e[1]=(p3[1]-B[1]); e[2]=(p3[2]-B[2]); P_float l21=norm(e);
			e[0]=(p4[0]-B[0]); e[1]=(p4[1]-B[1]); e[2]=(p4[2]-B[2]); P_float l22=norm(e);
			e[0]=(A[0]-B[0]); e[1]=(A[1]-B[1]); e[2]=(A[2]-B[2]); P_float ln=norm(e);
			test= (abs(l[ign1]-Topo_lref) + abs(l[ign2]-Topo_lref))/2 - Topo_Tol*(abs(l11-Topo_lref) + abs(l12-Topo_lref) + abs(l21-Topo_lref) + abs(l22-Topo_lref) + abs(ln-Topo_lref))/5;
			if(test>0 && test>testmax) if(TO2_test(i,i1,i2,i3,i4)) {testmax=test; imax=i; imax1=i1; imax2=i2; imax3=i3; imax4=i4;}
			}
	    free(l);
		}
	if(testmax!=-1E10) 	{TO2(imax,imax1,imax2,imax3,imax4); Subdivide_Update();}
    count++;
	}
}
*/

void CSimplexSurf::Subdivide_Merge(P_float p1[3],P_float p2[3])
{
CLOSESTSTRUCT *closest=new CLOSESTSTRUCT[1]; 
closest->weights=new P_float[1]; closest->pts=new int[1];     closest->dist2=1E10; closest->nb=0; 
GetClosest(closest,p1); int cc1=closest->cell;
free(closest->pts); free(closest->weights); closest->weights=new P_float[1]; closest->pts=new int[1];     closest->dist2=1E10; closest->nb=0; 
GetClosest(closest,p2); int cc2=closest->cell;

int offset;
P_float q1[3],q2[3],v[3],l,min=1E10; GetPoint(cells[cc1][1],q1); for(int i=0;i<cells[cc2][0];i++) {GetPoint(cells[cc2][i+1],q2); v[0]=q1[0]-q2[0];v[1]=q1[1]-q2[1]; v[2]=q1[2]-q2[2]; l=norm(v); if(l<min) {min=l;offset=i;}}
TOmerge(cc1,cc2,offset);
Equilibrium();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CSimplexSurf::Subdivide_test(P_float p[3]) //subdivide from an input point 
{
//Monito_TO("C:\\Temp\\data2\\generic_male2\\tests\\testbed\\topo_femur_analyse.txt");

int pt=GetClosestPoint(p); int pt2=neighbors[3*pt];
CLOSESTSTRUCT *closest=new CLOSESTSTRUCT[1]; closest->pts=new int[1]; closest->weights=new P_float[1]; closest->nb=0; closest->dist2=1E10; GetClosest(closest,p); int c=closest->cell; free(closest->pts); free(closest->weights); free(closest);
TOcell_d(c);
TO_test(); 
}


void CSimplexSurf::Monito_TO(const char* filename)
{
FILE* f=fopen(filename,"rt+");
if(f==NULL) f=fopen(filename,"wt"); 
else fseek(f,0,SEEK_END);

fprintf(f,"\n");
// write number of operations 
fprintf(f,"TO2: %d (%d tests)\n",Monitor_TO2,Monitor_TO2_test);
fprintf(f,"TO1: %d (%d tests)\n",Monitor_TO1,Monitor_TO1_test);
fprintf(f,"TOech: %d (%d tests)\n",Monitor_TOech,Monitor_TOech_test);
fprintf(f,"TOech2: %d (%d tests)\n",Monitor_TOech2,Monitor_TOech2_test);
fprintf(f,"TOcell_i: %d (%d tests)\n",Monitor_TOcell_i,Monitor_TOcell_i_test);
fprintf(f,"TOcell_d: %d (%d tests)\n",Monitor_TOcell_d,Monitor_TOcell_d_test);
fprintf(f,"TOvertex_i: %d (%d tests)\n",Monitor_TOvertex_i,Monitor_TOvertex_i_test);
fprintf(f,"TOvertex_d: %d (%d tests)\n",Monitor_TOvertex_d,Monitor_TOvertex_d_test);

fprintf(f,"\n");
// write triangle mesh quality (radius ratio hitogram)
int i,j,nb_BINS=100; P_float val,*p1,*p2,*p3;
int *histo=new int[nb_BINS]; for(i=0;i<nb_BINS;i++) histo[i]=0;
for(i=0;i<nb_cells;i++) 
	{
	p3=cellcenters+3*i;
	for(j=0;j<cells[i][0];j++) 
		{
		p1=points+3*cells[i][j+1];
		p2=(j==cells[i][0]-1)?points+3*cells[i][1]:points+3*cells[i][j+2];
		val=(P_float)(nb_BINS-1)*GetTriangleQuality(p1,p2,p3);
		histo[(int)val]++;
		}
	}
for(i=0;i<nb_BINS;i++) fprintf(f,"%d\n",histo[i]);

fprintf(f,"\n");
fclose(f);
}

void CSimplexSurf::TO_test()
{/*
if(!IsAxis)
	{
	for(int i=0;i<3*nb_points;i++)
		if(neighbors_c[i]<0)
			int oihin=0;
	int c[3]; 
	for(int i=0;i<nb_points;i++)
		{	
		SearchCell(i,c);
		if(c[0]<0 || c[1]<0 || c[2]<0)
			int oihsdfo=0;
		}
	}
for(int i=0;i<3*nb_points;i++)
	if(neighbors[i]<0)
		int uiu=0;
*/
	
// test neighbors_c

int nc[3],*nnc;
for(int i=0;i<nb_points;i++)
	{
	nnc=neighbors_c+3*i;
	SearchCell(i,nc);
	if(nnc[0]==-1 || nnc[1]==-1)
		int oisgno=0;

	if(nnc[0]!=nc[0] && nnc[0]!=nc[1] && nnc[0]!=nc[2])
		int oisgno=0;
	if(nnc[1]!=nc[0] && nnc[1]!=nc[1] && nnc[1]!=nc[2])
		int oisgno=0;
	if(nnc[2]!=nc[0] && nnc[2]!=nc[1] && nnc[2]!=nc[2])
		int oisgno=0;
	}

// test if points are in exactly 3 cells;
int c[10],count=0; 
for(int i=0;i<nb_points;i++)
	{
	count=0;
	for(int j=0;j<nb_cells;j++) 
	    for(int k=0;k<cells[j][0];k++) 
		    if(cells[j][k+1]==i) {c[count]=j;count++;}
	if(count!=3) 
			int oisgno=0;
	}

}

/////////////////////////////////////////////////////////////////////
bool CSimplexSurf::TO1_test(int p1,int p2)
{
Monitor_TO1_test++;
int c[2][2]={-1,-1,-1,-1},*c1=neighbors_c+3*p1,*c2=neighbors_c+3*p2,count=0;;
if(c1[0]!=c2[0] && c1[0]!=c2[1] && c1[0]!=c2[2]) c[0][0]=c1[0]; else if(count<2) {c[1][count]=c1[0];count++;}
if(c1[1]!=c2[0] && c1[1]!=c2[1] && c1[1]!=c2[2]) c[0][0]=c1[1]; else if(count<2) {c[1][count]=c1[1];count++;}
if(c1[2]!=c2[0] && c1[2]!=c2[1] && c1[2]!=c2[2]) c[0][0]=c1[2]; else if(count<2) {c[1][count]=c1[2];count++;}
if(c2[0]!=c1[0] && c2[0]!=c1[1] && c2[0]!=c1[2]) c[0][1]=c2[0]; 
if(c2[1]!=c1[0] && c2[1]!=c1[1] && c2[1]!=c1[2]) c[0][1]=c2[1];
if(c2[2]!=c1[0] && c2[2]!=c1[1] && c2[2]!=c1[2]) c[0][1]=c2[2]; 

if(c[1][1]==-1 || c[1][0]==-1)  return false; // don't delete border edge

if(surfaces_cell[c[1][0]]*surfaces_cell[c[1][1]]*Topo_Tol*Topo_Tol>Topo_scellref*Topo_scellref/2.) return false; //criterion

if(attached_point[p1]!=0 || attached_point[p2]!=0) return false;
if(Topo_ignoreborder && (c[1][0]==-1 || c[0][0]==-1 || c[1][1]==-1 || c[0][1]==-1)) return false;
if(c[0][0]!=-1) if(cells[c[0][0]][0]<=4)  return false;
if(c[0][1]!=-1) if(cells[c[0][1]][0]<=4)  return false;
if(c[1][1]!=-1 && c[1][0]!=-1) if(cells[c[1][1]][0]==3 && cells[c[1][0]][0]==3) return false;
return true;
}

void CSimplexSurf::TO1(int p1,int p2)
{
Monitor_TO1++;
// search for neighbors and adjacent cells
int c[2][2]={-1,-1,-1,-1},*c1=neighbors_c+3*p1,*c2=neighbors_c+3*p2,count=0;;
if(c1[0]!=c2[0] && c1[0]!=c2[1] && c1[0]!=c2[2]) c[0][0]=c1[0]; else if(count<2) {c[1][count]=c1[0];count++;}
if(c1[1]!=c2[0] && c1[1]!=c2[1] && c1[1]!=c2[2]) c[0][0]=c1[1]; else if(count<2) {c[1][count]=c1[1];count++;}
if(c1[2]!=c2[0] && c1[2]!=c2[1] && c1[2]!=c2[2]) c[0][0]=c1[2]; else if(count<2) {c[1][count]=c1[2];count++;}
if(c2[0]!=c1[0] && c2[0]!=c1[1] && c2[0]!=c1[2]) c[0][1]=c2[0]; 
if(c2[1]!=c1[0] && c2[1]!=c1[1] && c2[1]!=c1[2]) c[0][1]=c2[1];
if(c2[2]!=c1[0] && c2[2]!=c1[1] && c2[2]!=c1[2]) c[0][1]=c2[2]; 


// modify cells
if(c[0][0]!=-1) DeletePointInCell(c[0][0],p1);
if(c[0][1]!=-1) DeletePointInCell(c[0][1],p2);

if(c[1][1]==-1) {DeleteCell(c[1][0]); }
else if(c[1][0]==-1) {DeleteCell(c[1][1]); }
else
	{
	int* cell=Concatenate(cells[c[1][0]],cells[c[1][1]],p1,p2); 
	SetCell(c[1][0],cell); 
	free(cell);
	for(int i=0;i<3*nb_points;i++) if(neighbors_c[i]==c[1][1]) neighbors_c[i]=c[1][0];
	DeleteCell(c[1][1]);
	}

// Update Neighbors of order 1
int *np1=GetNeighbors(p1),*np2=GetNeighbors(p2),*n;
if(np1[0]==p2) {n=GetNeighbors(np1[1]); if(n[0]==p1) n[0]=np1[2]; if(n[1]==p1) n[1]=np1[2]; if(n[2]==p1) n[2]=np1[2]; n=GetNeighbors(np1[2]); if(n[0]==p1) n[0]=np1[1]; if(n[1]==p1) n[1]=np1[1]; if(n[2]==p1) n[2]=np1[1];}
if(np1[1]==p2) {n=GetNeighbors(np1[0]); if(n[0]==p1) n[0]=np1[2]; if(n[1]==p1) n[1]=np1[2]; if(n[2]==p1) n[2]=np1[2]; n=GetNeighbors(np1[2]); if(n[0]==p1) n[0]=np1[0]; if(n[1]==p1) n[1]=np1[0]; if(n[2]==p1) n[2]=np1[0];}
if(np1[2]==p2) {n=GetNeighbors(np1[1]); if(n[0]==p1) n[0]=np1[0]; if(n[1]==p1) n[1]=np1[0]; if(n[2]==p1) n[2]=np1[0]; n=GetNeighbors(np1[0]); if(n[0]==p1) n[0]=np1[1]; if(n[1]==p1) n[1]=np1[1]; if(n[2]==p1) n[2]=np1[1];}
if(np2[0]==p1) {n=GetNeighbors(np2[1]); if(n[0]==p2) n[0]=np2[2]; if(n[1]==p2) n[1]=np2[2]; if(n[2]==p2) n[2]=np2[2]; n=GetNeighbors(np2[2]); if(n[0]==p2) n[0]=np2[1]; if(n[1]==p2) n[1]=np2[1]; if(n[2]==p2) n[2]=np2[1];}
if(np2[1]==p1) {n=GetNeighbors(np2[0]); if(n[0]==p2) n[0]=np2[2]; if(n[1]==p2) n[1]=np2[2]; if(n[2]==p2) n[2]=np2[2]; n=GetNeighbors(np2[2]); if(n[0]==p2) n[0]=np2[0]; if(n[1]==p2) n[1]=np2[0]; if(n[2]==p2) n[2]=np2[0];}
if(np2[2]==p1) {n=GetNeighbors(np2[1]); if(n[0]==p2) n[0]=np2[0]; if(n[1]==p2) n[1]=np2[0]; if(n[2]==p2) n[2]=np2[0]; n=GetNeighbors(np2[0]); if(n[0]==p2) n[0]=np2[1]; if(n[1]==p2) n[1]=np2[1]; if(n[2]==p2) n[2]=np2[1];}

// Delete points 
if(p2>p1) {DeletePoint(p2); DeletePoint(p1);} else {DeletePoint(p1); DeletePoint(p2);}
}



bool CSimplexSurf::TO2_test(const int c,const P_float tresh,int pt[4],P_float l[2])
{
Monitor_TO2_test++;
if(c==-1) return false;
if(tresh<=0) 
	if(surfaces_cell[c]<sqrt(2.)*Topo_scellref*Topo_Tol) return false;

if(!GetCutEdges_PrincipalDirection(c,tresh,pt,l)) return false;

if(attached_point[pt[0]]!=0 || attached_point[pt[1]]!=0 || attached_point[pt[2]]!=0 || attached_point[pt[3]]!=0) return false;
if(Topo_ignoreborder)
	{
	int ca,cb;
	if(neighbors_c[3*pt[0]]!=c && (neighbors_c[3*pt[0]]==neighbors_c[3*pt[1]] || neighbors_c[3*pt[0]]==neighbors_c[3*pt[1]+1] || neighbors_c[3*pt[0]]==neighbors_c[3*pt[1]+2])) ca=neighbors_c[3*pt[0]]; else if(neighbors_c[3*pt[0]+1]!=c && (neighbors_c[3*pt[0]+1]==neighbors_c[3*pt[1]] || neighbors_c[3*pt[0]+1]==neighbors_c[3*pt[1]+1] || neighbors_c[3*pt[0]+1]==neighbors_c[3*pt[1]+2])) ca=neighbors_c[3*pt[0]+1]; else ca=neighbors_c[3*pt[0]+2];
	if(neighbors_c[3*pt[2]]!=c && (neighbors_c[3*pt[2]]==neighbors_c[3*pt[3]] || neighbors_c[3*pt[2]]==neighbors_c[3*pt[3]+1] || neighbors_c[3*pt[2]]==neighbors_c[3*pt[3]+2])) cb=neighbors_c[3*pt[2]]; else if(neighbors_c[3*pt[2]+1]!=c && (neighbors_c[3*pt[2]+1]==neighbors_c[3*pt[3]] || neighbors_c[3*pt[2]+1]==neighbors_c[3*pt[3]+1] || neighbors_c[3*pt[2]+1]==neighbors_c[3*pt[3]+2])) cb=neighbors_c[3*pt[2]+1]; else cb=neighbors_c[3*pt[2]+2];
	if(ca==-1 || cb==-1) return false;
	}

return true;
}

// p1, p2, p3 and p4 must be entered in counterclockwise or clockwise orientation
// Probably the same cell orientation as in the mesh. 
void CSimplexSurf::TO2(int c,int p1,int p2,int p3,int p4)
{
Monitor_TO2++;
int i;

int ca,cb;
if(neighbors_c[3*p1]!=c && (neighbors_c[3*p1]==neighbors_c[3*p2] || neighbors_c[3*p1]==neighbors_c[3*p2+1] || neighbors_c[3*p1]==neighbors_c[3*p2+2])) ca=neighbors_c[3*p1]; else if(neighbors_c[3*p1+1]!=c && (neighbors_c[3*p1+1]==neighbors_c[3*p2] || neighbors_c[3*p1+1]==neighbors_c[3*p2+1] || neighbors_c[3*p1+1]==neighbors_c[3*p2+2])) ca=neighbors_c[3*p1+1]; else ca=neighbors_c[3*p1+2];
if(neighbors_c[3*p3]!=c && (neighbors_c[3*p3]==neighbors_c[3*p4] || neighbors_c[3*p3]==neighbors_c[3*p4+1] || neighbors_c[3*p3]==neighbors_c[3*p4+2])) cb=neighbors_c[3*p3]; else if(neighbors_c[3*p3+1]!=c && (neighbors_c[3*p3+1]==neighbors_c[3*p4] || neighbors_c[3*p3+1]==neighbors_c[3*p4+1] || neighbors_c[3*p3+1]==neighbors_c[3*p4+2])) cb=neighbors_c[3*p3+1]; else cb=neighbors_c[3*p3+2];

int a=nb_points,b=nb_points+1;
P_float A[3]={(0.5*GetPoint(p3)[0]+0.5*GetPoint(p4)[0]+GetPoint(p1)[0]+GetPoint(p2)[0])/3,(0.5*GetPoint(p3)[1]+0.5*GetPoint(p4)[1]+GetPoint(p1)[1]+GetPoint(p2)[1])/3,(0.5*GetPoint(p3)[2]+0.5*GetPoint(p4)[2]+GetPoint(p1)[2]+GetPoint(p2)[2])/3};
P_float B[3]={(0.5*GetPoint(p1)[0]+0.5*GetPoint(p2)[0]+GetPoint(p3)[0]+GetPoint(p4)[0])/3,(0.5*GetPoint(p1)[1]+0.5*GetPoint(p2)[1]+GetPoint(p3)[1]+GetPoint(p4)[1])/3,(0.5*GetPoint(p1)[2]+0.5*GetPoint(p2)[2]+GetPoint(p3)[2]+GetPoint(p4)[2])/3};
P_float SA[3]={(0.5*GetSpeed(p3)[0]+0.5*GetSpeed(p4)[0]+GetSpeed(p1)[0]+GetSpeed(p2)[0])/3,(0.5*GetSpeed(p3)[1]+0.5*GetSpeed(p4)[1]+GetSpeed(p1)[1]+GetSpeed(p2)[1])/3,(0.5*GetSpeed(p3)[2]+0.5*GetSpeed(p4)[2]+GetSpeed(p1)[2]+GetSpeed(p2)[2])/3};
P_float SB[3]={(0.5*GetSpeed(p1)[0]+0.5*GetSpeed(p2)[0]+GetSpeed(p3)[0]+GetSpeed(p4)[0])/3,(0.5*GetSpeed(p1)[1]+0.5*GetSpeed(p2)[1]+GetSpeed(p3)[1]+GetSpeed(p4)[1])/3,(0.5*GetSpeed(p1)[2]+0.5*GetSpeed(p2)[2]+GetSpeed(p3)[2]+GetSpeed(p4)[2])/3};

// modify adjacent cells
InsertPointInCell(ca,a,p1); InsertPointInCell(cb,b,p3);

// add cells
int* cell=cells[c];
int* cell1=new int[cell[0]+2];
int* cell2=new int[cell[0]+2];
Separate(cell,p1,a,p2,p3,b,p4,cell1,cell2);
InsertNextCell(cell1); InsertNextCell(cell2);
InsertNextPoint(A);InsertNextPoint(B);
SetSpeed(a,SA); SetSpeed(b,SB);
free(cell1); free(cell2);
DeleteCell(c);

// Update Neighbors of order 1
neighbors[3*a]=p1; neighbors[3*a+1]=b; neighbors[3*a+2]=p2;
neighbors[3*b]=p3; neighbors[3*b+1]=a; neighbors[3*b+2]=p4;
for(i=0;i<3;i++)
	{
	if(neighbors[3*p1+i]==p2) neighbors[3*p1+i]=a;
	if(neighbors[3*p2+i]==p1) neighbors[3*p2+i]=a;
	if(neighbors[3*p4+i]==p3) neighbors[3*p4+i]=b;
	if(neighbors[3*p3+i]==p4) neighbors[3*p3+i]=b;
	}

cell=cells[nb_cells-2];	for(i=0;i<cell[0];i++) UpdateNeighbors_c(cell[i+1]);
cell=cells[nb_cells-1];	for(i=0;i<cell[0];i++) UpdateNeighbors_c(cell[i+1]);
}

/////////////////////////////////////////////
bool CSimplexSurf::TOcell_i_test(int c)
{
Monitor_TOcell_i_test++;
if(c==-1) return false;
if(surfaces_cell[c]<sqrt((P_float)cells[c][0]+1.)*Topo_scellref*Topo_Tol) return false; // criterion
for(int i=0;i<cells[c][0];i++) if(attached_point[cells[c][i+1]]!=0) return false;
if(Topo_ignoreborder) {for(int i=0;i<cells[c][0];i++) if(neighbors_c[3*cells[c][i+1]+2]==-1) return false;}
return true;
}

void CSimplexSurf::TOcell_i(int c) // insert a cell inside c
{
Monitor_TOcell_i++;
int i,j,*cell=cells[c],pt1,pt2,ce,cb[6],A,B,Ap,Bp,Am,Bm,*cc;
P_float p1[3],p2[3];
cb[0]=5;

// insert points
for(i=0;i<cell[0];i++)
	{
	pt1=cell[i+1]; pt2=(i==cell[0]-1)?cell[1]:cell[i+2];
	// compute the two points to add per edge
	p1[0]=(points[3*pt1]+points[3*pt2])/2.; p1[1]=(points[3*pt1+1]+points[3*pt2+1])/2.;	p1[2]=(points[3*pt1+2]+points[3*pt2+2])/2.;	InsertNextPoint(p1); 
	p2[0]=(p1[0]+cellcenters[3*c])/2.; p2[1]=(p1[1]+cellcenters[3*c+1])/2.;	p2[2]=(p1[2]+cellcenters[3*c+2])/2.;	InsertNextPoint(p2); 
	}

// insert cells and compute neighbors
cc=new int[cell[0]+1]; cc[0]=cell[0];
for(i=0;i<cell[0];i++)
	{
	pt1=cell[i+1]; pt2=(i==cell[0]-1)?cell[1]:cell[i+2];
	A=nb_points-2*cell[0]+2*i;	B=nb_points-2*cell[0]+2*i+1;
	Ap=(i==cell[0]-1)?(nb_points-2*cell[0]):(nb_points-2*cell[0]+2*i+2);	Bp=(i==cell[0]-1)?(nb_points-2*cell[0]+1):(nb_points-2*cell[0]+2*i+3);
	Am=(i==0)?(nb_points-2):(nb_points-2*cell[0]+2*i-2);	Bm=(i==0)?(nb_points-1):(nb_points-2*cell[0]+2*i-1);

	// update neighbors
	neighbors[3*A]=pt1; neighbors[3*A+1]=B; neighbors[3*A+2]=pt2;
	neighbors[3*B]=Bp; neighbors[3*B+1]=A; neighbors[3*B+2]=Bm;
	for(j=0;j<3;j++) {if(neighbors[3*pt1+j]==pt2) neighbors[3*pt1+j]=A; if(neighbors[3*pt2+j]==pt1) neighbors[3*pt2+j]=A;}

	// update neighbor cell (adjacent cell != c)
	for(j=0;j<3;j++) {if(neighbors_c[3*pt1+j]!=c) if(neighbors_c[3*pt1+j]==neighbors_c[3*pt2] || neighbors_c[3*pt1+j]==neighbors_c[3*pt2+1] ||neighbors_c[3*pt1+j]==neighbors_c[3*pt2+2]) ce=neighbors_c[3*pt1+j];}
	if(ce!=-1) InsertPointInCell(ce,A,pt1);
	
	// insert border cell 
	cb[1]=Am; cb[2]=pt1; cb[3]=A; cb[4]=B; cb[5]=Bm; InsertNextCell(cb);
	
	// update center cell
	cc[i+1]=B;

	// update pt1 neighbors_c
	if(neighbors_c[3*pt1]==c) neighbors_c[3*pt1]=nb_cells-1;	else if(neighbors_c[3*pt1+1]==c) neighbors_c[3*pt1+1]=nb_cells-1;	else neighbors_c[3*pt1+2]=nb_cells-1;
	}
SetCell(c,cc); free(cc);

// update A/B neighbors_c
cell=cells[c]; for(i=0;i<cell[0];i++)	{UpdateNeighbors_c(nb_points-2*cell[0]+2*i); UpdateNeighbors_c(nb_points-2*cell[0]+2*i+1);}
}

bool CSimplexSurf::TOcell_d_test(int c) // delete cell c
{
Monitor_TOcell_d_test++;
int i;
if(c==-1) return false;
int *cell=cells[c],ptp,pt,ptm,np,nc;
P_float sum=surfaces_cell[c];
for(i=0;i<cell[0];i++) if(attached_point[cell[i+1]]!=0 || neighbors_c[3*cell[i+1]+2]==-1) return false;
for(i=0;i<cell[0];i++)
	{
	ptm=(i==0)?cell[cell[0]]:cell[i]; pt=cell[i+1]; ptp=(i==cell[0]-1)?cell[1]:cell[i+2];
	if(neighbors[3*pt]!=ptp && neighbors[3*pt]!=ptm) np=neighbors[3*pt];	else if(neighbors[3*pt+1]!=ptp && neighbors[3*pt+1]!=ptm) np=neighbors[3*pt+1];	else if(neighbors[3*pt+2]!=ptp && neighbors[3*pt+2]!=ptm) np=neighbors[3*pt+2];
	if(attached_point[np]!=0 || neighbors_c[3*np+2]==-1) return false;
	if(neighbors_c[3*np]!=neighbors_c[3*pt] && neighbors_c[3*np]!=neighbors_c[3*pt+1] && neighbors_c[3*np]!=neighbors_c[3*pt+2]) nc=neighbors_c[3*np];	else if(neighbors_c[3*np+1]!=neighbors_c[3*pt] && neighbors_c[3*np+1]!=neighbors_c[3*pt+1] && neighbors_c[3*np+1]!=neighbors_c[3*pt+2]) nc=neighbors_c[3*np+1];	else if(neighbors_c[3*np+2]!=neighbors_c[3*pt] && neighbors_c[3*np+2]!=neighbors_c[3*pt+1] && neighbors_c[3*np+2]!=neighbors_c[3*pt+2]) nc=neighbors_c[3*np+2];
	if(cells[nc][0]==3) return false;
	sum+=surfaces_cell[nc];
	}
if(sum*Topo_Tol>sqrt((P_float)cells[c][0]+1.)*Topo_scellref) return false; // criterion
return true;
}

void CSimplexSurf::TOcell_d(int c)
{
Monitor_TOcell_d++;
int i,j,*cell=cells[c],*nc,*nc2,*np,ptm,pt,ptp,pt1,pt2,*ci,*cl,np2,c0;
ci=new int[1]; nc=new int[cell[0]+1]; nc2=new int[cell[0]]; np=new int[2*cell[0]];
// cfind adjacent point/cells
ci[0]=0;
for(i=0;i<cell[0];i++)
	{
	ptm=(i==0)?cell[cell[0]]:cell[i]; pt=cell[i+1]; ptp=(i==cell[0]-1)?cell[1]:cell[i+2];
	if(neighbors[3*pt]!=ptp && neighbors[3*pt]!=ptm) np[i]=neighbors[3*pt];
	else if(neighbors[3*pt+1]!=ptp && neighbors[3*pt+1]!=ptm) np[i]=neighbors[3*pt+1];
	else if(neighbors[3*pt+2]!=ptp && neighbors[3*pt+2]!=ptm) np[i]=neighbors[3*pt+2];
	for(j=0;j<3;j++) {if(neighbors_c[3*pt+j]!=c) if(neighbors_c[3*pt+j]==neighbors_c[3*ptp] || neighbors_c[3*pt+j]==neighbors_c[3*ptp+1] ||neighbors_c[3*pt+j]==neighbors_c[3*ptp+2]) nc[i]=neighbors_c[3*pt+j];}
	if(neighbors_c[3*np[i]]!=neighbors_c[3*pt] && neighbors_c[3*np[i]]!=neighbors_c[3*pt+1] && neighbors_c[3*np[i]]!=neighbors_c[3*pt+2]) nc2[i]=neighbors_c[3*np[i]];
	else if(neighbors_c[3*np[i]+1]!=neighbors_c[3*pt] && neighbors_c[3*np[i]+1]!=neighbors_c[3*pt+1] && neighbors_c[3*np[i]+1]!=neighbors_c[3*pt+2]) nc2[i]=neighbors_c[3*np[i]+1];
	else if(neighbors_c[3*np[i]+2]!=neighbors_c[3*pt] && neighbors_c[3*np[i]+2]!=neighbors_c[3*pt+1] && neighbors_c[3*np[i]+2]!=neighbors_c[3*pt+2]) nc2[i]=neighbors_c[3*np[i]+2];
	np[i+cell[0]]=pt; // for ordering
	}
nc[cell[0]]=c;
// update neighbors add delete point in cells
for(i=0;i<cell[0];i++)
	{
	pt=cell[i+1];
	// find the neighbors of np[i] != pt
	if(neighbors[3*np[i]]==pt) {pt1=neighbors[3*np[i]+1]; pt2=neighbors[3*np[i]+2];}	else if(neighbors[3*np[i]+1]==pt) {pt1=neighbors[3*np[i]+2]; pt2=neighbors[3*np[i]];}	else if(neighbors[3*np[i]+2]==pt) {pt1=neighbors[3*np[i]]; pt2=neighbors[3*np[i]+1];}
	// update pt1 and pt2 's neighbors
	if(neighbors[3*pt1]==np[i]) neighbors[3*pt1]=pt2;	else if(neighbors[3*pt1+1]==np[i]) neighbors[3*pt1+1]=pt2;	else if(neighbors[3*pt1+2]==np[i]) neighbors[3*pt1+2]=pt2;
	if(neighbors[3*pt2]==np[i]) neighbors[3*pt2]=pt1;	else if(neighbors[3*pt2+1]==np[i]) neighbors[3*pt2+1]=pt1;	else if(neighbors[3*pt2+2]==np[i]) neighbors[3*pt2+2]=pt1;
	DeletePointInCell(nc2[i],np[i]);	
	
	cl=cells[nc[i]]; np2=(i==cell[0]-1)?np[0]:np[i+1];
	j=1; 
	while(cl[j]!=np[i]) 
		j++;	
	
	j++; 
	
	if(j==cl[0]+1) 
		j=1; // put j after np[i]

	while(cl[j]!=np2) {ci[0]++; ci=(int*)realloc(ci,(ci[0]+1)*sizeof(int)); ci[ci[0]]=cl[j]; if(j==cl[0]) j=1; else j++;} 
	}
// Delete cells and points
c0=cell[0];
qsort(nc, cell[0]+1, sizeof(int),qsortintmaxtomin);
qsort(np, 2*cell[0], sizeof(int),qsortintmaxtomin);
for(i=0;i<c0+1;i++) DeleteCell(nc[i]);
InsertNextCell(ci); 
for(i=0;i<2*c0;i++) DeletePoint(np[i]);
free(nc); free(nc2); free(np); free(ci);

//Update neighbors_c
for(i=0;i<cells[nb_cells-1][0];i++)	UpdateNeighbors_c(cells[nb_cells-1][i+1]);
}

/////////////////////////////////////////////
bool CSimplexSurf::TOvertex_i_test(int p)
{
Monitor_TOvertex_i_test++;
int* np=neighbors+3*p,*nc=neighbors_c+3*p;
int *nc1=neighbors_c+3*np[0],*nc2=neighbors_c+3*np[1],*nc3=neighbors_c+3*np[2],ncc1,ncc2,ncc3;

if(nc[2]==-1 || nc1[2]==-1 || nc2[2]==-1 || nc3[2]==-1) return false;

if(attached_point[p]!=0 || attached_point[np[0]]!=0 || attached_point[np[1]]!=0 || attached_point[np[2]]!=0) return false;
if(cells[nc[0]][0]==3 || cells[nc[1]][0]==3 || cells[nc[2]][0]==3) return false;

if(nc1[0]!=nc[0] && nc1[0]!=nc[1] && nc1[0]!=nc[2]) ncc1=nc1[0]; else if(nc1[1]!=nc[0] && nc1[1]!=nc[1] && nc1[1]!=nc[2]) ncc1=nc1[1]; else ncc1=nc1[2];
if(nc2[0]!=nc[0] && nc2[0]!=nc[1] && nc2[0]!=nc[2]) ncc2=nc2[0]; else if(nc2[1]!=nc[0] && nc2[1]!=nc[1] && nc2[1]!=nc[2]) ncc2=nc2[1]; else ncc2=nc2[2];
if(nc3[0]!=nc[0] && nc3[0]!=nc[1] && nc3[0]!=nc[2]) ncc3=nc3[0]; else if(nc3[1]!=nc[0] && nc3[1]!=nc[1] && nc3[1]!=nc[2]) ncc3=nc3[1]; else ncc3=nc3[2];
P_float sum=surfaces_cell[nc[0]]+surfaces_cell[nc[1]]+surfaces_cell[nc[2]]+surfaces_cell[ncc1]+surfaces_cell[ncc2]+surfaces_cell[ncc3];
if(sum<sqrt(42.)*Topo_Tol*Topo_scellref) return false; // criterion

return true;
}

void CSimplexSurf::TOvertex_i(int p) // replace p and its neighbors by a 6-cell
{
Monitor_TOvertex_i++;
int* np=neighbors+3*p,*nc=neighbors_c+3*p,n1,n2,n3,p1,p2,i;
for(i=0;i<cells[nc[0]][0];i++)
	if(cells[nc[0]][i+1]==p)
		{
		n1=(i==0)?cells[nc[0]][cells[nc[0]][0]]:cells[nc[0]][i];		if(i==0) p1=cells[nc[0]][cells[nc[0]][0]-1];	else if(i==1) p1=cells[nc[0]][cells[nc[0]][0]];	else p1=cells[nc[0]][i-1];
		n2=(i==cells[nc[0]][0]-1)?cells[nc[0]][1]:cells[nc[0]][i+2];	if(i==cells[nc[0]][0]-1) p2=cells[nc[0]][2];	else if(i==cells[nc[0]][0]-2) p2=cells[nc[0]][1];	else p2=cells[nc[0]][i+3];
		}
if(np[0]!=n1 && np[0]!=n2) n3=np[0]; else if(np[1]!=n1 && np[1]!=n2) n3=np[1]; else n3=np[2];
TO2(nc[0],p1,n1,n2,p2);
TOech(p,n3);
}

bool CSimplexSurf::TOvertex_d_test(int c,bool *first) 
{
Monitor_TOvertex_d_test++;
int i;
if(c==-1) return false;
if(cells[c][0]!=6) return false;
for(i=0;i<cells[c][0];i++) if(attached_point[cells[c][i+1]]!=0) return false;
P_float sum=surfaces_cell[c];
int *cell=cells[c],ptp,pt,nc,n1[3],n2[3];
for(i=0;i<cell[0];i++)
	{
	pt=cell[i+1]; ptp=(i==cell[0]-1)?cell[1]:cell[i+2];
	if(neighbors_c[3*pt]!=c && (neighbors_c[3*pt]==neighbors_c[3*ptp] || neighbors_c[3*pt]==neighbors_c[3*ptp+1] || neighbors_c[3*pt]==neighbors_c[3*ptp+2])) nc=neighbors_c[3*pt];	else if(neighbors_c[3*pt+1]!=c && (neighbors_c[3*pt+1]==neighbors_c[3*ptp] || neighbors_c[3*pt+1]==neighbors_c[3*ptp+1] || neighbors_c[3*pt+1]==neighbors_c[3*ptp+2])) nc=neighbors_c[3*pt+1];	else nc=neighbors_c[3*pt+2];
	if(nc==-1) return false; //Topo_ignoreborder
	if(i==0 || i==2 || i==4) n1[i/2]=nc;
	else n2[i/2]=nc;
	sum+=surfaces_cell[nc];
	}

if(cells[n1[0]][0]+cells[n1[1]][0]+cells[n1[2]][0]>cells[n2[0]][0]+cells[n2[1]][0]+cells[n2[2]][0]) *first=true; else *first=false;
if(first) {if(cells[n1[0]][0]==3 || cells[n1[1]][0]==3 || cells[n1[2]][0]==3) return false; } else {if(cells[n2[0]][0]==3 || cells[n2[1]][0]==3 || cells[n2[2]][0]==3) return false; }
	
if(sum*Topo_Tol>sqrt(42.)*Topo_scellref) return false; // criterion
return true;
}

void CSimplexSurf::TOvertex_d(int c,bool first) // replace a 6-cell by a vertex
{
Monitor_TOvertex_d++;
int *cell=cells[c],p1,p2,p3,p4;
if(first) {p1=cell[1]; p2=cell[2]; p3=cell[4]; p4=cell[5];}
else {p1=cell[2]; p2=cell[3]; p3=cell[5]; p4=cell[6];}
TOech(p1,p2);
TO1(p3,p4);
}


/////////////////////////////////////////////////
bool CSimplexSurf::TOech2_test(int p)
{
Monitor_TOech2_test++;
int i,*np=neighbors+3*p,*nc=neighbors_c+3*p,nc2[3];
for(i=0;i<3;i++) 
	{
	if(attached_point[np[i]]!=0) return false; 
	if(neighbors_c[3*np[i]]!=nc[0] && neighbors_c[3*np[i]]!=nc[1] && neighbors_c[3*np[i]]!=nc[2]) nc2[i]=neighbors_c[3*np[i]]; else if(neighbors_c[3*np[i]+1]!=nc[0] && neighbors_c[3*np[i]+1]!=nc[1] && neighbors_c[3*np[i]+1]!=np[2]) nc2[i]=neighbors_c[3*np[i]+1]; else nc2[i]=neighbors_c[3*np[i]+2];
	}

//if(Topo_ignoreborder) //debug
	if(nc[0]==-1 || nc[1]==-1 || nc[2]==-1 || nc2[0]==-1 || nc2[1]==-1 || nc2[2]==-1) return false;

if(cells[nc[0]][0]<=4 || cells[nc[1]][0]<=4 || cells[nc[2]][0]<=4) return false;
if(cells[nc[0]][0] + cells[nc[1]][0] + cells[nc[2]][0] - cells[nc2[0]][0] - cells[nc2[1]][0] - cells[nc2[2]][0] <=6) return false;
return true;
}

void CSimplexSurf::TOech2(int p) // 2nd-order exchange
{
Monitor_TOech2++;
int np[3]; GetNeighbors(p,np);

TOech(np[0],p);
TOech(np[0],np[2]);
TOech(np[1],p);
TOech(p,np[2]);
}


bool CSimplexSurf::TOech_test(int p1,int p2)
{
Monitor_TOech_test++;
if(attached_point[p1]!=0 || attached_point[p2]!=0) return false; 

int c[2][2]={-1,-1,-1,-1},*c1=neighbors_c+3*p1,*c2=neighbors_c+3*p2,count=0;;
if(c1[0]!=c2[0] && c1[0]!=c2[1] && c1[0]!=c2[2]) c[0][0]=c1[0]; else if(count<2) {c[1][count]=c1[0];count++;}
if(c1[1]!=c2[0] && c1[1]!=c2[1] && c1[1]!=c2[2]) c[0][0]=c1[1]; else if(count<2) {c[1][count]=c1[1];count++;}
if(c1[2]!=c2[0] && c1[2]!=c2[1] && c1[2]!=c2[2]) c[0][0]=c1[2]; else if(count<2) {c[1][count]=c1[2];count++;}
if(c2[0]!=c1[0] && c2[0]!=c1[1] && c2[0]!=c1[2]) c[0][1]=c2[0]; 
if(c2[1]!=c1[0] && c2[1]!=c1[1] && c2[1]!=c1[2]) c[0][1]=c2[1];
if(c2[2]!=c1[0] && c2[2]!=c1[1] && c2[2]!=c1[2]) c[0][1]=c2[2]; 

if(c[0][0]==-1 && c[0][1]==-1) return false;
if(Topo_ignoreborder) if(c[0][0]==-1 || c[0][1]==-1  || c[1][0]==-1  || c[1][1]==-1) return false;
if(c[1][0]!=-1) if(cells[c[1][0]][0]==3) return false; // 4?
if(c[1][1]!=-1) if(cells[c[1][1]][0]==3) return false; // 4?
	
int n1,n2,n3,n4;
if(c[1][0]==-1) n1=cells[c[1][1]][0]; else n1=cells[c[1][0]][0];
if(c[1][1]==-1) n2=cells[c[1][0]][0]; else n2=cells[c[1][1]][0];
if(c[0][0]==-1) n3=cells[c[0][1]][0]; else n3=cells[c[0][0]][0];
if(c[0][1]==-1) n4=cells[c[0][0]][0]; else n4=cells[c[0][1]][0];

if(n1+n2-n3-n4<=2) return false;

// prevent from deleting small tubes
int i,j;
if(c[0][0]!=-1 && c[0][1]!=-1)
	for(i=0;i<cells[c[0][0]][0];i++)
		for(j=0;j<cells[c[0][1]][0];j++)
			if(cells[c[0][0]][i+1]==cells[c[0][1]][j+1])
				return false;
return true;

}

void CSimplexSurf::TOech(int p1,int p2) // 1st-order exchange
{
Monitor_TOech++;
int i,j;

// search for neighbors and adjacent cells
int c[2][2]={-1,-1,-1,-1},*c1=neighbors_c+3*p1,*c2=neighbors_c+3*p2,count=0;
if(c1[0]!=c2[0] && c1[0]!=c2[1] && c1[0]!=c2[2]) c[0][0]=c1[0]; else if(count<2) {c[1][count]=c1[0];count++;}
if(c1[1]!=c2[0] && c1[1]!=c2[1] && c1[1]!=c2[2]) c[0][0]=c1[1]; else if(count<2) {c[1][count]=c1[1];count++;}
if(c1[2]!=c2[0] && c1[2]!=c2[1] && c1[2]!=c2[2]) c[0][0]=c1[2]; else if(count<2) {c[1][count]=c1[2];count++;}
if(c2[0]!=c1[0] && c2[0]!=c1[1] && c2[0]!=c1[2]) c[0][1]=c2[0]; 
if(c2[1]!=c1[0] && c2[1]!=c1[1] && c2[1]!=c1[2]) c[0][1]=c2[1];
if(c2[2]!=c1[0] && c2[2]!=c1[1] && c2[2]!=c1[2]) c[0][1]=c2[2]; 

// calcul new position of p1 and p2
int p3,p4,p5,p6;
if(neighbors[3*p1]==p2) {p3=neighbors[3*p1+1]; p4=neighbors[3*p1+2];}
if(neighbors[3*p1+1]==p2) {p3=neighbors[3*p1+2]; p4=neighbors[3*p1];}
if(neighbors[3*p1+2]==p2) {p3=neighbors[3*p1]; p4=neighbors[3*p1+1];}
if(neighbors[3*p2]==p1) {p5=neighbors[3*p2+1]; p6=neighbors[3*p2+2];}
if(neighbors[3*p2+1]==p1) {p5=neighbors[3*p2+2]; p6=neighbors[3*p2];}
if(neighbors[3*p2+2]==p1) {p5=neighbors[3*p2]; p6=neighbors[3*p2+1];}

// c[1][0] should contain p3 , otherwise switch c[1][0] and c[1][1]
bool swtch=false;
if(c[1][1]!=-1) {for(i=0; i<cells[c[1][1]][0]; i++)	if(cells[c[1][1]][i+1]==p3) swtch=true;}
else {for(i=0; i<cells[c[1][0]][0]; i++)	if(cells[c[1][0]][i+1]==p4) swtch=true;}
if(swtch) {j=c[1][0]; c[1][0]=c[1][1]; c[1][1]=j;}

P_float P1[3]={(points[3*p1]+points[3*p2]+points[3*p3]+points[3*p6])/4,(points[3*p1+1]+points[3*p2+1]+points[3*p3+1]+points[3*p6+1])/4,(points[3*p1+2]+points[3*p2+2]+points[3*p3+2]+points[3*p6+2])/4};
P_float P2[3]={(points[3*p1]+points[3*p2]+points[3*p4]+points[3*p5])/4,(points[3*p1+1]+points[3*p2+1]+points[3*p4+1]+points[3*p5+1])/4,(points[3*p1+2]+points[3*p2+2]+points[3*p4+2]+points[3*p5+2])/4};
P_float SP1[3]={(speeds[3*p1]+speeds[3*p2]+speeds[3*p3]+speeds[3*p6])/4,(speeds[3*p1+1]+speeds[3*p2+1]+speeds[3*p3+1]+speeds[3*p6+1])/4,(speeds[3*p1+2]+speeds[3*p2+2]+speeds[3*p3+2]+speeds[3*p6+2])/4};
P_float SP2[3]={(speeds[3*p1]+speeds[3*p2]+speeds[3*p4]+speeds[3*p5])/4,(speeds[3*p1+1]+speeds[3*p2+1]+speeds[3*p4+1]+speeds[3*p5+1])/4,(speeds[3*p1+2]+speeds[3*p2+2]+speeds[3*p4+2]+speeds[3*p5+2])/4};

SetPoint(p1,P1);SetPoint(p2,P2);
SetSpeed(p1,SP1);SetSpeed(p2,SP2);
 
// modify cells
InsertPointInCell(c[0][0],p2,p4); InsertPointInCell(c[0][1],p1,p6);
DeletePointInCell(c[1][0],p2); DeletePointInCell(c[1][1],p1);

neighbors[3*p1]=p3; neighbors[3*p1+1]=p2; neighbors[3*p1+2]=p6;
neighbors[3*p2]=p5; neighbors[3*p2+1]=p1; neighbors[3*p2+2]=p4;
if(neighbors[3*p4]==p1)  neighbors[3*p4]=p2; else if(neighbors[3*p4+1]==p1)  neighbors[3*p4+1]=p2; else  neighbors[3*p4+2]=p2;
if(neighbors[3*p6]==p2)  neighbors[3*p6]=p1; else if(neighbors[3*p6+1]==p2)  neighbors[3*p6+1]=p1; else  neighbors[3*p6+2]=p1;

if(c[0][0]==-1) {neighbors_c[3*p1]=c[0][1]; neighbors_c[3*p1+1]=c[1][0]; neighbors_c[3*p1+2]=c[0][0];} else if(c[0][1]==-1) {neighbors_c[3*p1]=c[1][0]; neighbors_c[3*p1+1]=c[0][0]; neighbors_c[3*p1+2]=c[0][1];} else {neighbors_c[3*p1]=c[0][0]; neighbors_c[3*p1+1]=c[0][1]; neighbors_c[3*p1+2]=c[1][0];}
if(c[0][0]==-1) {neighbors_c[3*p2]=c[1][1]; neighbors_c[3*p2+1]=c[0][1]; neighbors_c[3*p2+2]=c[0][0];} else if(c[0][1]==-1) {neighbors_c[3*p2]=c[0][0]; neighbors_c[3*p2+1]=c[1][1]; neighbors_c[3*p2+2]=c[0][1];} else {neighbors_c[3*p2]=c[0][1]; neighbors_c[3*p2+1]=c[0][0]; neighbors_c[3*p2+2]=c[1][1];}
}

///////////////////////////////////////////////////////////////////
void CSimplexSurf::TOmerge(int c1,int c2,int offset)
{
if(cells[c1][0]!=cells[c2][0]) return;

int i1,i2,p1,i,j,c[2][2],nb=cells[c1][0]; int* corres= new int[2*nb]; int* corres_c= new int[2*nb];
for(i=0;i<nb;i++) { corres[2*i]=cells[c1][i+1]; corres[(i>offset)?(2*(offset-i+nb)+1):(2*(offset-i)+1)]=cells[c2][(i+1+offset>nb)?(i+1+offset-nb):(i+1+offset)];}
for(i=0;i<nb;i++) { SearchCell(corres[2*i],corres[(i==nb-1)?0:(2*i+2)],c); if(c[1][0]!=c1) corres_c[2*i]=c[1][0]; else corres_c[2*i]=c[1][1];  SearchCell(corres[2*i+1],corres[(i==nb-1)?1:(2*i+3)],c); if(c[1][0]!=c2) corres_c[2*i+1]=c[1][0]; else corres_c[2*i+1]=c[1][1];}

for(i=0;i<nb;i++)
    {
    i1=1; while(cells[corres_c[2*i]][i1]!=corres[2*i]) i1++; i1=(i1==1)?(cells[corres_c[2*i]][0]):(i1-1); p1=cells[corres_c[2*i]][i1];
    i2=1; while(cells[corres_c[2*i+1]][i2]!=corres[2*i+1]) i2++; i2=(i2>cells[corres_c[2*i+1]][0]-2)?(i2-cells[corres_c[2*i+1]][0]+2):(i2+2);
    for(j=0;j<cells[corres_c[2*i+1]][0]-2;j++)
        {
        InsertPointInCell(corres_c[2*i],cells[corres_c[2*i+1]][i2],p1);
        i2=i2+1; if(i2==cells[corres_c[2*i+1]][0]+1) i2=1;
        }
    DeletePointInCell(corres_c[2*i],corres[2*i]);
    DeletePointInCell(corres_c[2*i],p1);
    }

// delete nb+2 cells
int* values=new int[nb+2]; for(i=0;i<nb;i++) values[i+2]=corres_c[2*i+1]; values[0]=c1; values[1]=c2;
qsort(values, nb+2, sizeof(int),qsortintmaxtomin);
for(i=0;i<nb+2;i++) DeleteCell(values[i]);
// delete 2nb points
qsort(corres, 2*nb, sizeof(int),qsortintmaxtomin);
for(i=0;i<nb*2;i++) DeletePoint(corres[i]);

free(values); free(corres);free(corres_c);
UpdateNeighbors();
}



/////////////////////////////////////////////////////////////////////////////////////////////////
int qsortintmaxtomin (const void * a, const void * b) {  return ( *(int*)b - *(int*)a ); }
int qsortP_floatmaxtomin (const void * a, const void * b) {  P_float f1 = *( (P_float *) a);    P_float f2 = *( (P_float *) b);    if (f1 < f2) return 1;  if (f1 == f2) return 0;  return -1;}

void CSimplexSurf::decompose(int index,P_float f[3],P_float ftg[3],P_float fn[3]) 
{
P_float n[3]; GetNormal(index,n);
P_float cp=dotproduct(f,n);
fn[0]=cp*n[0];fn[1]=cp*n[1];fn[2]=cp*n[2];
ftg[0]=f[0]-fn[0]; ftg[1]=f[1]-fn[1]; ftg[2]=f[2]-fn[2];
}

/*
void CSimplexSurf::GetClosest(CLOSESTSTRUCT* closest,P_float p[3],P_float p1[3],P_float p2[3],P_float p3[3]) 
{

closest->cell=-1;
P_float p1p2[3]={p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
P_float p1p3[3]={p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2]};
P_float pp1[3]={p1[0]-p[0],p1[1]-p[1],p1[2]-p[2]};
P_float norml[3]; crossproduct(norml,p1p2,p1p3);

if(dotproduct(pp1,norml)<=0) closest->side=1; else closest->side=-1;

P_float fout[2];
closest->dist2=GetClosest(dotproduct(p1p2,p1p2),dotproduct(p1p2,p1p3),dotproduct(p1p3,p1p3),dotproduct(p1p2,pp1),dotproduct(p1p3,pp1),dotproduct(pp1,pp1),fout);
P_float fu=1-fout[0]-fout[1];

closest->c[0]=p1[0]+fout[0]*p1p2[0]+fout[1]*p1p3[0];closest->c[1]=p1[1]+fout[0]*p1p2[1]+fout[1]*p1p3[1];closest->c[2]=p1[2]+fout[0]*p1p2[2]+fout[1]*p1p3[2];
closest->nb=(fu==0)?(0):(1);closest->nb+=(fout[0]==0)?(0):(1);closest->nb+=(fout[1]==0)?(0):(1);
closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb]; 
int count=0;
if(fu!=0) {closest->weights[0]=fu; closest->pts[0]=0; count++;}
if(fout[0]!=0) {closest->weights[count]=fout[0]; closest->pts[count]=1; count++;}
if(fout[1]!=0) {closest->weights[count]=fout[1]; closest->pts[count]=2; count++;}
}

void CSimplexSurf::GetClosest(CLOSESTSTRUCT* closest,P_float p[3],P_float p1[3],P_float p1p2[3],P_float p1p3[3],P_float norml[3],P_float P1P2, P_float P1P3,P_float P1P2P1P3) 
{

closest->cell=-1;
P_float pp1[3]={p1[0]-p[0],p1[1]-p[1],p1[2]-p[2]};
if(dotproduct(pp1,norml)<=0) closest->side=1; else closest->side=-1;

P_float fout[2];
closest->dist2=GetClosest(dotproduct(p1p2,p1p2),dotproduct(p1p2,p1p3),dotproduct(p1p3,p1p3),dotproduct(p1p2,pp1),dotproduct(p1p3,pp1),dotproduct(pp1,pp1),fout);
P_float fu=1-fout[0]-fout[1];

closest->c[0]=p1[0]+fout[0]*p1p2[0]+fout[1]*p1p3[0];closest->c[1]=p1[1]+fout[0]*p1p2[1]+fout[1]*p1p3[1];closest->c[2]=p1[2]+fout[0]*p1p2[2]+fout[1]*p1p3[2];
closest->nb=(fu==0)?(0):(1);closest->nb+=(fout[0]==0)?(0):(1);closest->nb+=(fout[1]==0)?(0):(1);
closest->weights=new P_float[closest->nb]; closest->pts=new int[closest->nb]; 
int count=0;
if(fu!=0) {closest->weights[0]=fu; closest->pts[0]=0; count++;}
if(fout[0]!=0) {closest->weights[count]=fout[0]; closest->pts[count]=1; count++;}
if(fout[1]!=0) {closest->weights[count]=fout[1]; closest->pts[count]=2; count++;}
}
*/

bool CSimplexSurf::GetCutEdges_PrincipalDirection(const int index_cell,P_float tresh,int pt[4],P_float l[2]) 
{
pt[0]=-1; pt[1]=-1; pt[2]=-1; pt[3]=-1;
int i,cross=0; 
int *c=GetCell(index_cell); 
P_float *pts=new P_float[3*c[0]]; for(i=0;i<c[0];i++) GetPoint(c[i+1],pts+3*i);
P_float pm[3],axis[3][3]; GetPrincipalDirections(pm,axis,c[0],pts,false); l[0]=norm(axis[0]); l[1]=norm(axis[1]); 
if(l[0]/l[1]<tresh) {free(pts); return false;}

P_float v,d=dotproduct(pm,axis[0]);
int pt2=c[0]; bool sgn=(dotproduct(pts+3*pt2,axis[0])>d)?true:false;
for(i=0;i<c[0];i++)
	{
	v=dotproduct(pts+3*i,axis[0]); 
	if((sgn && v<d)|| (!sgn && v>d)) 
		{
		if(cross<2)
			{
			pt[2*cross]=c[pt2]; 
			pt[2*cross+1]=c[i+1];
			}
		cross++;
		}
	sgn=(v>d)?true:false; pt2=i+1;
	}
free(pts);
//if(pt[0]==pt[2] || pt[1]==pt[3] || neighbors[3*pt[0]]==pt[2] || neighbors[3*pt[0]+1]==pt[2] || neighbors[3*pt[0]+2]==pt[2] || neighbors[3*pt[1]]==pt[3] || neighbors[3*pt[1]+1]==pt[3] || neighbors[3*pt[1]+2]==pt[3])	return false; // return true if edges are not neighboring edges
//else 

if(cross<2) return false;

return true;
}




/////////////////////////////////////////////////////////////////////////////
/*
bool CSimplexSurf::TOstar_test(int p)
{
if(attached_point[p]!=0) return false; 
if(Topo_ignoreborder && neighbors_c[3*p+2]==-1)	{return false;}
return true;
}
void CSimplexSurf::TOstar(int p)
{
// increase resolution locally
// star->triangle transformation
int n1=neighbors[3*p],n2=neighbors[3*p+1],n3=neighbors[3*p+2];
int c[3]={SearchCell(n1,p,n2),SearchCell(n3,p,n1),SearchCell(n2,p,n3)};

P_float A[3]={(points[3*p]+points[3*n1])/2,(points[3*p+1]+points[3*n1+1])/2,(points[3*p+2]+points[3*n1+2])/2};    
P_float B[3]={(points[3*p]+points[3*n2])/2,(points[3*p+1]+points[3*n2+1])/2,(points[3*p+2]+points[3*n2+2])/2};    
P_float C[3]={(points[3*p]+points[3*n3])/2,(points[3*p+1]+points[3*n3+1])/2,(points[3*p+2]+points[3*n3+2])/2};    
P_float SA[3]={(speeds[3*p]+speeds[3*n1])/2,(speeds[3*p+1]+speeds[3*n1+1])/2,(speeds[3*p+2]+speeds[3*n1+2])/2};    
P_float SB[3]={(speeds[3*p]+speeds[3*n2])/2,(speeds[3*p+1]+speeds[3*n2+1])/2,(speeds[3*p+2]+speeds[3*n2+2])/2};    
P_float SC[3]={(speeds[3*p]+speeds[3*n3])/2,(speeds[3*p+1]+speeds[3*n3+1])/2,(speeds[3*p+2]+speeds[3*n3+2])/2};    

int a=nb_points,b=nb_points+1;
InsertNextPoint(A);InsertNextPoint(B); SetPoint(p,C);
SetSpeed(a,SA); SetSpeed(b,SB);	SetSpeed(p,SC);

DeletePointInCell(c[0],p);InsertPointInCell(c[0],a,n2);InsertPointInCell(c[0],b,n2);InsertPointInCell(c[2],b,p);InsertPointInCell(c[1],a,n1);
int cell[4]={3,p,b,a};
InsertNextCell(cell);

// Update Neighbors of order 1
SetNeighbors(p,-1,-1,-1);SetNeighbors(n1,-1,-1,-1);SetNeighbors(n2,-1,-1,-1);SetNeighbors(n3,-1,-1,-1);
UpdateNeighbors(nb_cells-1);UpdateNeighbors(c[0]);UpdateNeighbors(c[1]);UpdateNeighbors(c[2]);

if(c[0]==-1 || c[1]==-1 || c[2]==-1) 
	{
	int cn[3];
	SearchCell(n1,cn); if(cn[0]!=c[0] && cn[0]!=c[1] && cn[0]!=c[2]) UpdateNeighbors(cn[0]); if(cn[1]!=c[0] && cn[1]!=c[1] && cn[1]!=c[2]) UpdateNeighbors(cn[1]); if(cn[2]!=c[0] && cn[2]!=c[1] && cn[2]!=c[2]) UpdateNeighbors(cn[2]);
	SearchCell(n2,cn); if(cn[0]!=c[0] && cn[0]!=c[1] && cn[0]!=c[2]) UpdateNeighbors(cn[0]); if(cn[1]!=c[0] && cn[1]!=c[1] && cn[1]!=c[2]) UpdateNeighbors(cn[1]); if(cn[2]!=c[0] && cn[2]!=c[1] && cn[2]!=c[2]) UpdateNeighbors(cn[2]);
	SearchCell(n3,cn); if(cn[0]!=c[0] && cn[0]!=c[1] && cn[0]!=c[2]) UpdateNeighbors(cn[0]); if(cn[1]!=c[0] && cn[1]!=c[1] && cn[1]!=c[2]) UpdateNeighbors(cn[1]); if(cn[2]!=c[0] && cn[2]!=c[1] && cn[2]!=c[2]) UpdateNeighbors(cn[2]);
	}
if(TOech_test(a,n1,-1)) TOech(a,n1);
if(TOech_test(b,n2,-1)) TOech(b,n2);
if(TOech_test(p,n3,-1)) TOech(p,n3);
}


bool CSimplexSurf::TOtri_test(int c)
{
int i;
for(i=0;i<cells[c][0];i++) if(attached_point[cells[c][i+1]]!=0) return false;
// condition on adjacent cells: nb_point(adj_c)>4
int adjc[3];
for(i=0;i<cells[c][0]-1;i++)
    {
    SearchCell(cells[c][i+1],adjc);
	if(adjc[0]!=-1) {if(cells[adjc[0]][0]<5) return false;} else if(Topo_ignoreborder) return false; 
	if(adjc[1]!=-1) {if(cells[adjc[1]][0]<5) return false;} else if(Topo_ignoreborder) return false; 
	if(adjc[2]!=-1) {if(cells[adjc[2]][0]<5) return false;} else if(Topo_ignoreborder) return false; 
    }
return true;
}

void CSimplexSurf::TOtri(int c)
{
// decrease resolution locally
// cell->triangle->star transformation
int i,j=1; 

// cell->triangle
int count=0;
while(cells[c][0]!=3) 
    {
    if(TOech_test(cells[c][j],cells[c][j+1],c)) TOech(cells[c][j],cells[c][j+1]); 
    j++; count++;
    if(j>=cells[c][0]) j=1;
	if(count==200) 
		return;
	}

int p1=cells[c][1]; int p2=cells[c][2]; int p3=cells[c][3];

int cs[3]; 
int ac[2][2]; SearchCell(p1,p2,ac); cs[0]=(ac[1][0]==c)?ac[1][1]:ac[1][0]; cs[1]=ac[0][0]; cs[2]=ac[0][1];

if(cs[0]==cs[1] || cs[0]==cs[2] || cs[1]==cs[2]) return; //to prevent from object topology change (tubular structure)

int cn2,cn3;
if(cs[0]==-1 || cs[1]==-1 || cs[2]==-1)
	{
	int np2; if(neighbors[3*p2]!=p3 && neighbors[3*p2]!=p1) np2=neighbors[3*p2];	else if(neighbors[3*p2+1]!=p3 && neighbors[3*p2+1]!=p1) np2=neighbors[3*p2+1]; else	np2=neighbors[3*p2+2];
	int np3; if(neighbors[3*p3]!=p2 && neighbors[3*p3]!=p1) np3=neighbors[3*p3];	else if(neighbors[3*p3+1]!=p2 && neighbors[3*p3+1]!=p1) np3=neighbors[3*p3+1]; else	np3=neighbors[3*p3+2];
	SearchCell(np2,p2,ac); cn2=ac[0][0];
	SearchCell(np3,p3,ac); cn3=ac[0][0];
	}

DeletePointInCell(cs[0],p2); DeletePointInCell(cs[1],p3); DeletePointInCell(cs[2],p2); if(cs[2]!=-1) for(i=0;i<cells[cs[2]][0];i++) if(cells[cs[2]][i+1]==p3) cells[cs[2]][i+1]=p1;
P_float P[3]={(points[3*p1]+points[3*p2]+points[3*p3])/3,(points[3*p1+1]+points[3*p2+1]+points[3*p3+1])/3,(points[3*p1+2]+points[3*p2+2]+points[3*p3+2])/3};
P_float SP[3]={(speeds[3*p1]+speeds[3*p2]+speeds[3*p3])/3,(speeds[3*p1+1]+speeds[3*p2+1]+speeds[3*p3+1])/3,(speeds[3*p1+2]+speeds[3*p2+2]+speeds[3*p3+2])/3};
SetPoint(p1,P); SetSpeed(p1,SP);

// Update Neighbors of order 1

if(neighbors[3*p2]==p3) {SetNeighbors(neighbors[3*p2+1],-1,-1,-1);SetNeighbors(neighbors[3*p2+2],-1,-1,-1);}
if(neighbors[3*p2+1]==p3) {SetNeighbors(neighbors[3*p2],-1,-1,-1);SetNeighbors(neighbors[3*p2+2],-1,-1,-1);}
if(neighbors[3*p2+2]==p3) {SetNeighbors(neighbors[3*p2],-1,-1,-1);SetNeighbors(neighbors[3*p2+1],-1,-1,-1);}
SetNeighbors(neighbors[3*p3],-1,-1,-1);SetNeighbors(neighbors[3*p3+1],-1,-1,-1);SetNeighbors(neighbors[3*p3+2],-1,-1,-1);
UpdateNeighbors(cs[0]);UpdateNeighbors(cs[1]);UpdateNeighbors(cs[2]);

if(cs[0]==-1 || cs[1]==-1 || cs[2]==-1) {UpdateNeighbors(cn2); UpdateNeighbors(cn3);}

DeleteCell(c);
if(p3>p2) {DeletePoint(p3); DeletePoint(p2);} else {DeletePoint(p2); DeletePoint(p3);}
}

void CSimplexSurf::Subdivide_Increase_star() // increase res locally if surface>surface_ref
{
int i,j,count=0,countmax=(Topo_Increase_L[1]==-1)?(nb_points):(Topo_Increase_L[1]);
P_float s,a,min;
while(count!=countmax)
    {
	min=1E10;
	for(j=0;j<nb_points;j++) //j=rand()*nb_points/RAND_MAX; //
		{	
		s=surfaces[j]+surfaces[GetNeighbors(j)[0]]+surfaces[GetNeighbors(j)[1]]+surfaces[GetNeighbors(j)[2]];
		a=5*Topo_sref-s;
		if(a<min && a<0) if(TOstar_test(j)) {min=a; i=j;}
		}
	if(min!=1E10) {TOstar(i); Subdivide_Update();}
    count++;
    }
}


void CSimplexSurf::Subdivide_Increase_R()
{
// increase res locally
int i,imax,count=0,countmax=(Topo_Increase_R[1]==-1)?(nb_points):(Topo_Increase_R[1]);
P_float max,Rm,stdev;
while(count!=countmax)
    {
	max=-1E10;
	for(i=0;i<nb_points;i++)
		{
		Rm=(axiallinks[i].Rref+axiallinks[neighbors[3*i]].Rref+axiallinks[neighbors[3*i+1]].Rref+axiallinks[neighbors[3*i+2]].Rref)/4; 
		stdev=(axiallinks[i].Rref-Rm)*(axiallinks[i].Rref-Rm) + (axiallinks[neighbors[3*i]].Rref-Rm)*(axiallinks[neighbors[3*i]].Rref-Rm)  + (axiallinks[neighbors[3*i+1]].Rref-Rm)*(axiallinks[neighbors[3*i+1]].Rref-Rm) + (axiallinks[neighbors[3*i+2]].Rref-Rm)*(axiallinks[neighbors[3*i+2]].Rref-Rm);
		stdev=sqrt(stdev/(P_float)(nb_points));
		if(stdev>max) if(TOstar_test(i)) {max=stdev; imax=i;}
		}
	if(max!=-1E10) {TOstar(imax); Subdivide_Update();}
    count++;
    }
}


void CSimplexSurf::Subdivide_Increase_Err()
{
// increase res locally
int i,imax,count=0,countmax=(Topo_Increase_Err[1]==-1)?(nb_points):(Topo_Increase_Err[1]);
P_float max;
while(count!=countmax)
    {
	max=-1E10;
	for(i=0;i<nb_points;i++)	if(axiallinks[i].err>max) if(TOstar_test(i)) {max=axiallinks[i].err; imax=i;}
	if(max!=-1E10) {TOstar(imax); Subdivide_Update();}
    count++;
    }
}

void CSimplexSurf::Subdivide_Decrease_tri()
{
// decrease res locally
int i,j,k,count=0,countmax=(Topo_Decrease_L[1]==-1)?(nb_cells):(Topo_Decrease_L[1]);
P_float min,a,s;
while(count!=countmax)
    {
	min=1E10;
	for(j=0;j<nb_cells;j++) //j=rand()*nb_cells/RAND_MAX; //
		{	
		s=0; for(k=0;k<cells[j][0];k++) s+=surfaces[cells[j][k+1]]; s=6*s/cells[j][0]; 
		a=Topo_Tol*abs(s/(P_float)4-Topo_sref)-abs(s/(P_float)6-Topo_sref);
		if(a<min && a<0) if(TOtri_test(j)) {min=a; i=j;}
		}
    if(s<Topo_sref-Topo_Tol) {TOtri(i); Subdivide_Exchange(); } //debug
    if(s<Topo_sref) {TOtri(i); Subdivide_Exchange(); } 
	if(min!=1E10) {TOtri(i); Subdivide_Update();} 
	count++;
    }
}



void CSimplexSurf::Subdivide_Decrease_R()
{
// decrease res locally
int i,imin,count=0,countmax=(Topo_Decrease_R[1]==-1)?(nb_points):(Topo_Decrease_R[1]);
P_float min,Rm,stdev;
while(count!=countmax)
    {
	min=1E10;
	for(i=0;i<nb_points;i++)
		{
		Rm=(axiallinks[i].Rref+axiallinks[neighbors[3*i]].Rref+axiallinks[neighbors[3*i+1]].Rref+axiallinks[neighbors[3*i+2]].Rref)/4; 
		stdev=(axiallinks[i].Rref-Rm)*(axiallinks[i].Rref-Rm) + (axiallinks[neighbors[3*i]].Rref-Rm)*(axiallinks[neighbors[3*i]].Rref-Rm)  + (axiallinks[neighbors[3*i+1]].Rref-Rm)*(axiallinks[neighbors[3*i+1]].Rref-Rm) + (axiallinks[neighbors[3*i+2]].Rref-Rm)*(axiallinks[neighbors[3*i+2]].Rref-Rm);
		if(stdev<min) if(TOtri_test(i)) {min=stdev; imin=i;}
		}
	if(min!=1E10) {TOtri(imin); Subdivide_Update();}
    count++;
    }
}

void CSimplexSurf::Subdivide_Decrease_Err()
{
// increase res locally
int i,imin,count=0,countmax=(Topo_Decrease_Err[1]==-1)?(nb_points):(Topo_Decrease_Err[1]);
P_float min;
while(count!=countmax)
    {
	min=1E10;
	for(i=0;i<nb_points;i++)	if(axiallinks[i].err<min) if(TOtri_test(i)) {min=axiallinks[i].err; imin=i;}
	if(min!=1E10) {TOtri(imin); Subdivide_Update();}
    count++;
    }
}

void CSimplexSurf::SetTopo_Increase_R(int frequ,int nb,bool ignoreborder) {Topo_ignoreborder=ignoreborder; Topo_Increase_R[0]=frequ; Topo_Increase_R[1]=nb;  Update_Topo();}
void CSimplexSurf::SetTopo_Decrease_R(int frequ,int nb,bool ignoreborder) {Topo_ignoreborder=ignoreborder; Topo_Decrease_R[0]=frequ; Topo_Decrease_R[1]=nb; Update_Topo();}
void CSimplexSurf::SetTopo_Increase_Err(int frequ,int nb,bool ignoreborder) {Topo_ignoreborder=ignoreborder; Topo_Increase_Err[0]=frequ; Topo_Increase_Err[1]=nb;  Update_Topo();}
void CSimplexSurf::SetTopo_Decrease_Err(int frequ,int nb,bool ignoreborder) {Topo_ignoreborder=ignoreborder; Topo_Decrease_Err[0]=frequ; Topo_Decrease_Err[1]=nb; Update_Topo();}


*/

P_float CSimplexSurf::getSurface() {
	return surface;
}

P_float CSimplexSurf::getSurfaces_cell(int index) {
	return surfaces_cell[index];
}

vtkSmartPointer<vtkPolyData> CSimplexSurf::getAsVTKPolyData() {
	vtkSmartPointer<vtkPolyData> pdata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

	for (int i = 0; i < nb_points; i++) {
		double *d = GetPoint(i);
		points->InsertNextPoint(d[0], d[1], d[2]);
	}

	for (int j = 0; j < nb_cells; j++) {		
		vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();

		int *c = GetCell(j);
		for (int k = 1; k <= c[0]; k++) {
			list->InsertNextId(c[k]);
		}

		cells->InsertNextCell(list);
	}

	pdata->SetPoints(points);
	pdata->SetPolys(cells);
	
	// Material indices and scalars are stored as vtkPointData, in three arrays. 
	// The first two arrays contain the two material indices, 
	// and the third array contain the scalar values. 
	vtkSmartPointer<vtkIntArray> matvals1 = vtkSmartPointer<vtkIntArray>::New();
	vtkSmartPointer<vtkIntArray> matvals2 = vtkSmartPointer<vtkIntArray>::New();
	vtkSmartPointer<vtkIntArray> scalars = vtkSmartPointer<vtkIntArray>::New();

	matvals1->SetNumberOfValues(GetNumberOfPoints());
	matvals1->SetNumberOfComponents(1);
	matvals1->SetName("PointMaterialIndices1");

	matvals2->SetNumberOfValues(GetNumberOfPoints());
	matvals2->SetNumberOfComponents(1);
	matvals2->SetName("PointMaterialIndices2");

	scalars->SetNumberOfValues(GetNumberOfPoints());
	scalars->SetNumberOfComponents(1);
	scalars->SetName("scalars");
	
	for (int mi = 0; mi < GetNumberOfPoints(); mi++) {
		int *mat = GetPointMaterialIndices(mi);
		matvals1->SetValue(mi, mat[0]);
		matvals2->SetValue(mi, mat[1]);
		scalars->SetValue(mi, mat[2]);
	}
 
	pdata->GetPointData()->AddArray(matvals1);
	pdata->GetPointData()->AddArray(matvals2);
	pdata->GetPointData()->AddArray(scalars);


	// Extra info. If the mesh is a split of a larger, 
	// multi-material mesh, then there is extra info
	// i.e. original mesh point ids. 
	if (GetSplitMesh() == true) {
		vtkSmartPointer<vtkIntArray> originalIds = vtkSmartPointer<vtkIntArray>::New();
		originalIds->SetNumberOfValues(GetNumberOfPoints());
		originalIds->SetNumberOfComponents(1);
		originalIds->SetName("OriginalIds");

		for (int mi = 0; mi < GetNumberOfPoints(); mi++) {
			originalIds->SetValue(mi, GetOriginalPointIndex(mi));
		}

		pdata->GetPointData()->AddArray(originalIds);
	}


	

	pdata->BuildCells();
	pdata->BuildLinks();
	pdata->Update();

	return pdata;
}

void CSimplexSurf::writeCSimplexMeshAsVTKPolyData(std::string filename) {
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInput(getAsVTKPolyData());
	writer->Write();
}

/**
 * Checks the connectivity, i.e. neighbors of each point 
 * to verify that the 2-Simplex criterion is true. 
*/
void CSimplexSurf::verify3PointConnectivity() {
	vtkSmartPointer<vtkPolyData> pdata = getAsVTKPolyData();

	for (int i = 0; i < pdata->GetNumberOfPoints(); i++) {
		
	}
}

/**
 * A triangular cell in the simplex mesh can be replaced with a 
 * single point using TO1(int p1, int p2) operation where p1 and
 * p2 are two points of the triangular cell. This function locates 
 * all triangular cells in the simplex mesh and performs a TO1() 
 * operation on two applicable points of the triangular cell. 
*/
void CSimplexSurf::decimateByRemovingTriangularCells() {
	cout << "Decimating by replacing triangular cells with a single point..." << endl;

	int nb_pts = GetNumberOfPoints();
	int nb_cells = GetNumberOfCells();

	cout << "\tNumber of points in mesh before decimation: " << GetNumberOfPoints() << endl;
	cout << "\tNumber of cells in mesh before decimation: " << GetNumberOfCells() << endl;

	int i = 0;
	while (i < GetNumberOfCells()) {
		int *cell = GetCell(i);
		
		if (cell[0] == 3) {
			if (TO1_test(cell[1], cell[2]) == true) {
				TO1(cell[1], cell[2]);
				i = i - 1;
			}
			else if (TO1_test(cell[2], cell[3]) == true) {
				TO1(cell[2], cell[3]);
				i = i - 1;
			}
			else if (TO1_test(cell[3], cell[1]) == true) {
				TO1(cell[3], cell[1]);
				i = i - 1;
			} 
			else {
				cout << "WARNING: Cannot safely apply TO1() on cell " << i << endl;
			}
		}

		i = i + 1;
	}

	UpdateNeighbors();
	UpdateAll();

	cout << endl;
	cout << "\tNumber of points in mesh after decimation: " << GetNumberOfPoints() << endl;
	cout << "\tNumber of cells in mesh after decimation: " << GetNumberOfCells() << endl;
}

void CSimplexSurf::decimateByMergingSmallSurfaceAreaCells(double factor) {
	double averageCellSurface = getSurface() / (double)GetNumberOfCells();
	
	vtkSmartPointer<vtkPolyData> pdata = getAsVTKPolyData();

	// For all cells
	for (int i = 0; i < GetNumberOfCells(); i++) {
		double d = getSurfaces_cell(i);

		if (i == 406) {
			writeCSimplexMeshAsVTKPolyData("C:\\Users\\Tanweer Rashid\\Desktop\\TEMP.vtk");
		}
		

		// If the surface area of the ith cell is less than (factor * averageSurfaceArea)
		if (d < averageCellSurface * factor) {
			cout << "Cell " << i << " has small surface area" << endl;

			vtkSmartPointer<vtkCell> icell = pdata->GetCell(i);
			vtkSmartPointer<vtkIdList> pointList = icell->GetPointIds(); // Get the points that make up the cell. 
			
			int *neighboringCellsIndex = new int[icell->GetNumberOfPoints()]; // To keep neighboring cells indices
			double *neighboringCellsSurface = new double[icell->GetNumberOfPoints()]; // To keep the surface areas of the neighboring cells
			int **neighboringCellsEdge_pts = new int*[icell->GetNumberOfPoints()]; // To keep the points making up the neighboring cell's edge. 

			int c = 0;			
			for (int j = 0; j < pointList->GetNumberOfIds(); j++) {
				int pt0 = pointList->GetId(j % pointList->GetNumberOfIds());
				int pt1 = pointList->GetId((j + 1) % pointList->GetNumberOfIds());

				vtkSmartPointer<vtkIdList> neighborList = vtkSmartPointer<vtkIdList>::New();
				pdata->GetCellEdgeNeighbors(i, pt0, pt1, neighborList);
				
				if (neighborList->GetNumberOfIds() > 1) {
					cout << "\tERROR: Too many neighbors detected in decimateByMergingSmallSurfaceAreaCells()" << endl; 
				}
				else if (neighborList->GetNumberOfIds() < 1) {
					cout << "\tERROR: 0 neighbors detected in decimateByMergingSmallSurfaceAreaCells()" << endl; 
				}
				else {
					neighboringCellsIndex[c] = neighborList->GetId(0);
					neighboringCellsSurface[c] = getSurfaces_cell(neighborList->GetId(0));
					
					neighboringCellsEdge_pts[c] = new int[2];
					neighboringCellsEdge_pts[c][0] = pt0;
					neighboringCellsEdge_pts[c][1] = pt1;

					c = c  + 1;
				}
			}
			
			cout << "\n\tCell " << i << "'s neighbors unsorted:" << endl;
			for (int e = 0; e < icell->GetNumberOfPoints(); e++) {
				cout << "\t\t" << neighboringCellsIndex[e] << "   " << neighboringCellsSurface[e] << "   " << neighboringCellsEdge_pts[e][0] << "   " << neighboringCellsEdge_pts[e][1] << endl;
			}

			// Sorting (basic sort) in ascending order wrt surface area. 
			for (int a = 0; a < icell->GetNumberOfPoints() - 1; a++) {
				for (int b = a + 1; b < icell->GetNumberOfPoints(); b++) {
					if (neighboringCellsSurface[b] < neighboringCellsSurface[a]) { // Swap
						int tempIndex = neighboringCellsIndex[a];
						double tempSurface = neighboringCellsSurface[a];
						int temp_pt0 = neighboringCellsEdge_pts[a][0];
						int temp_pt1 = neighboringCellsEdge_pts[a][1];

						neighboringCellsIndex[a] = neighboringCellsIndex[b];
						neighboringCellsSurface[a] = neighboringCellsSurface[b];
						neighboringCellsEdge_pts[a][0] = neighboringCellsEdge_pts[b][0];
						neighboringCellsEdge_pts[a][1] = neighboringCellsEdge_pts[b][1];

						neighboringCellsIndex[b] = tempIndex;
						neighboringCellsSurface[b] = tempSurface;
						neighboringCellsEdge_pts[b][0] = temp_pt0;
						neighboringCellsEdge_pts[b][1] = temp_pt1;
					}
				}
			}

			cout << "\n\tCell " << i << "'s neighbors sorted:" << endl;
			for (int e = 0; e < icell->GetNumberOfPoints(); e++) {
				cout << "\t\t" << neighboringCellsIndex[e] << "   " << neighboringCellsSurface[e] << "   " << neighboringCellsEdge_pts[e][0] << "   " << neighboringCellsEdge_pts[e][1] << endl;
			}

			// After sorting, perform safe T0() with the cell having smallest surface area
			for (int d = 0; d < icell->GetNumberOfPoints(); d++) {
				if (TO1_test(neighboringCellsEdge_pts[d][0], neighboringCellsEdge_pts[d][1]) == true) {
					TO1(neighboringCellsEdge_pts[d][0], neighboringCellsEdge_pts[d][1]);

					cout << "\tTO1() applied safely on edge made of points " << neighboringCellsEdge_pts[d][0] << " and " << neighboringCellsEdge_pts[d][1] << endl;

					UpdateNeighbors();
					UpdateAll();
					ComputeVolume();

					break;
				}
				else {
					cout << "\tTO1() cannot be safely applied to edge made of points " << neighboringCellsEdge_pts[d][0] << " and " << neighboringCellsEdge_pts[d][1] << ". Trying neighbor with next smallest surface. " << endl;

					if (d == (icell->GetNumberOfPoints() - 1)) {
						cout << "\tWARNING: TO1() operator could not be safely applied with any neighbors of cell " << i << endl;
					}
				}
			}

			//if (TO1_test(smallest_pt0, smallest_pt1) == true) {
			//	cout << "\n\tSmallest neighboring cell shares points " << smallest_pt0 << " and " << smallest_pt1 << endl;
			//	TO1(smallest_pt0, smallest_pt1);
			//}
			//else {
			//	cout << "\n\tWARNING: TO1() cannot be safely performed for points " << smallest_pt0 << " and " << smallest_pt1 << endl;
			//}
			


			cout << endl;
		}
	}

	UpdateNeighbors();
	UpdateAll();
}

void CSimplexSurf::verifyPointsAndParams() {
	//UpdateNormals();
	/*
	ofstream fpoints("C:\\Users\\Tanweer Rashid\\Desktop\\MRI_SEG_Points.txt");
	ofstream fnormals("C:\\Users\\Tanweer Rashid\\Desktop\\MRI_SEG_Normals.txt");

	for (int i = 0; i < GetNumberOfPoints(); i++) {
		double *p = GetPoint(i);
		
		fpoints << p[0] << ", " << p[1] << ", " << p[2] << "\n";
		fnormals << normals[3 * i + 0] << ", " << normals[3 * i + 1] << ", " << normals[3 * i + 2] << "\n";
	}

	fpoints.close();
	fnormals.close();
	*/
	

	ofstream fparams("C:\\Users\\Tanweer Rashid\\Desktop\\params.txt", ios::app);
	for (int k = 0; k < GetNumberOfPoints(); k++) {
		double *params = GetParams(k);
		double e1 = params[0]; 
		double e2 = params[1];
		double e3 = 1 - e1 - e2;
		double hh = params[2];
		//double h = this->h[k];


		double *P, *P1, *P2, *P3;
		P = GetPoint(k); 
		P1 = GetPoint(neighbors[3*k]); 
		P2 = GetPoint(neighbors[3*k+1]); 
		P3 = GetPoint(neighbors[3*k+2]); 

		double np[3];
	
		np[0] = (e1 * P1[0]) + (e2 * P2[0]) + (e3 * P3[0]) + (normals[3 * k + 0] * hh);
		np[1] = (e1 * P1[1]) + (e2 * P2[1]) + (e3 * P3[1]) + (normals[3 * k + 1] * hh);
		np[2] = (e1 * P1[2]) + (e2 * P2[2]) + (e3 * P3[2]) + (normals[3 * k + 2] * hh);



		fparams << (P[0] - np[0]) << "    " << (P[1] - np[1]) << "    " << (P[2] - np[2]) << "\n";
	}
	fparams << "================================================================\n";
	fparams.close();

} 