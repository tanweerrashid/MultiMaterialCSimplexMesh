
#include "StdFunctions.h"
#include "SimplexSurf.h"

class Deformation {

	public:
		Deformation(CSimplexSurf *newMesh, double timestep, double alpha);
		~Deformation();

		void ImpEuler();
				
		void setMesh(CSimplexSurf *newMesh);
		CSimplexSurf* getMesh();

		P_float getTimestep();
		void setTimestep(double ts);

		P_float getAlpha();
		void setAlpha(double a);

	private:
		CSimplexSurf *mesh;		
		P_float Timestep;	
		P_float Timestep2;
		P_float Alpha; 
		P_float Alpha2;
		
		void ImpEuler_ApplyResults(P_float* x);
		void ImpEuler_SetParams(P_float* r, P_float* h);
};