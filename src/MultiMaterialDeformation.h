
#include "StdFunctions.h"
#include "SimplexSurf.h"

class MultiMaterialDeformation {

	public:
		MultiMaterialDeformation(CSimplexSurf **newMesh, double timestep, double alpha, int num_meshes);
		~MultiMaterialDeformation();

		void ImpEuler();
				
		void setMeshes(CSimplexSurf **newMesh);
		CSimplexSurf** getMeshes();

		P_float getTimestep();
		void setTimestep(double ts);

		P_float getAlpha();
		void setAlpha(double a);

		int getNumberOfDeformableMeshes();
		void setNumberOfDeformableMeshes(int n);

	private:
		CSimplexSurf **meshes;		
		P_float Timestep;	
		P_float Timestep2;
		P_float Alpha; 
		P_float Alpha2;
		int numberOfDeformableMeshes;
		
		void MMDeform();
		void ImpEuler_ApplyResults(P_float* x);
		void ImpEuler_SetParams(P_float* r, P_float* h);
};
