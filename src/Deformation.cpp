
#include "Deformation.h"

Deformation::Deformation(CSimplexSurf *newMesh, double timestep, double alpha) {
	mesh = newMesh;

	Timestep = timestep;
	Timestep2 = Timestep * Timestep;

	Alpha = alpha;
	Alpha2 = Alpha * Alpha;
}

Deformation::~Deformation() {
}

void Deformation::setMesh(CSimplexSurf *newMesh) {
	mesh = newMesh;
}

CSimplexSurf* Deformation::getMesh() {
	return mesh;
}

P_float Deformation::getTimestep() {
	return Timestep;
}

void Deformation::setTimestep(double ts) {
	Timestep = ts;
	Timestep2 = Timestep * Timestep;
}

P_float Deformation::getAlpha() {
	return Alpha;
}

void Deformation::setAlpha(double a) {
	Alpha = a;
	Alpha2 = Alpha * Alpha;
}

void Deformation::ImpEuler() {
	int i, count = 0; P_float epsilon = 0.05;

	//this->UpdateForces();
	mesh->UpdateExternalForces();
	mesh->UpdateInternalForces();
	mesh->UpdateForcesDamping();

	int nb_Points = mesh->GetNumberOfPoints();

	//initialisation
	P_float bet = 0, alph, asb;
	P_float *x = new P_float[nb_Points*3], *h = new P_float[nb_Points * 6];
	P_float *t = new P_float[nb_Points * 3], *ht = new P_float[nb_Points * 3], *r = new P_float[nb_Points * 3];

	for (i = 0; i < 3 * nb_Points; i++) {
		x[i] = 0;
	}

	ImpEuler_SetParams(r, h); 
	//AddCollision_ImplIniParams(r);

	// cg
	while (-1) {
		alph = 0; 
		for (i = 0; i < 3 * nb_Points; i++) {
			alph += r[i] * r[i];
		}

		if (bet != 0) {
			asb = alph / bet; 
			for (i = 0; i < 3 * nb_Points; i++) {
				t[i] = r[i] + asb * t[i];
			}
		}
		else {
			memcpy(t, r, 3 * nb_Points * sizeof(P_float));
		}

		bet = 0;
		for (i = 0; i < nb_Points; i++) {
			Multi_MsT(ht + 3 * i, h + 6 * i, t + 3 * i); 
			bet += ht[3 * i] * t[3 * i] + ht[3 * i + 1] * t[3 * i + 1] + ht[3 * i + 2] * t[3 * i + 2];
		}

		//	bet+=AddCollision_ImplParams(ht,t);
		if (bet == 0) break;
		
		asb = alph / bet;
		for (i = 0; i < nb_Points; i++)	{
			r[3 * i] = r[3 * i] - asb * ht[3 * i]; 
			r[3 * i + 1] = r[3 * i + 1] - asb * ht[3 * i + 1];    
			r[3 * i + 2] = r[3 * i + 2] - asb * ht[3 * i + 2];

			x[3 * i] = x[3 * i] + asb * t[3 * i]; 
			x[3 * i + 1] = x[3 * i + 1] + asb * t[3 * i + 1]; 
			x[3 * i + 2] = x[3 * i + 2] + asb * t[3 * i + 2];
		}

		bet = alph; count++;
		if (bet < epsilon || count > 20) break;
	}

	// apply results
	ImpEuler_ApplyResults(x);

	// end
	free(x); free(h); free(t); free(r); free(ht); 

}

void Deformation::ImpEuler_ApplyResults(P_float* x) {
	P_float dep = 0, ds[3], *s, *p, disp[3]; 
	int i, count = 0;

	for (i = 0; i < mesh->GetNumberOfPoints(); i++) {
		P_float M = mesh->GetMass_inv(i)[0]; 
		s = mesh->GetSpeed(i); 
		p = mesh->GetPoint(i);

		// ds=L.x
		ds[0] = M * x[3 * count];
		ds[1] = M * x[3 * count + 1];
		ds[2] = M * x[3 * count + 2];

		// P+= (So + alpha*ds)*dt
		mesh->SetPoint_tm1(i, p);
		disp[0] = (s[0] + Alpha * ds[0]) * Timestep; 
		disp[1] = (s[1] + Alpha * ds[1]) * Timestep; 
		disp[2] = (s[2] + Alpha * ds[2]) * Timestep; 

		p[0] += disp[0];
		p[1] += disp[1];
		p[2] += disp[2];
		dep += norm(disp);

		// S+= ds
		mesh->SetSpeed_tm1(i, s);
		s[0] += ds[0];	
		s[1] += ds[1];	
		s[2] += ds[2];

		count++;
	}
		
}

void Deformation::ImpEuler_SetParams(P_float* r, P_float* h) {
	P_float *m, *dfp, *dfs, *s, *f;
	int i, count = 0;

	for (i = 0; i < mesh->GetNumberOfPoints(); i++) {
		m = mesh->GetMass_inv(i); 
		
		f = mesh->GetForce(i);
		s = mesh->GetSpeed(i); 
		dfs = mesh->GetDfs(i);  
		dfp = mesh->GetDfp(i);  

		P_float M = m[0], M2 = m[0] * m[0];
		//H= m.I -Alpha2*m^2.df/dp.dt^2 -Alpha*m^2.df/ds.dt
		h[6 * count] =- M2 * (Alpha2 * dfp[0] * Timestep2 + Alpha * dfs[0] * Timestep) + M; 			
		h[6 * count + 1] =- M2 * (Alpha2 * dfp[1] * Timestep2 + Alpha * dfs[1] * Timestep);   			
		h[6 * count + 2] =- M2 * (Alpha2 * dfp[2] * Timestep2 + Alpha * dfs[2] * Timestep) + M;   			
		h[6 * count + 3] =- M2 * (Alpha2 * dfp[3] * Timestep2 + Alpha * dfs[3] * Timestep);   			
		h[6 * count + 4] =- M2 * (Alpha2 * dfp[4] * Timestep2 + Alpha * dfs[4] * Timestep);   			
		h[6 * count + 5] =- M2 * (Alpha2 * dfp[5] * Timestep2 + Alpha * dfs[5] * Timestep) + M;

		// R = m.F.dt + Alpha*m.df/dp.S.dt^2 
		r[3 * count] = M * (f[0] * Timestep + Alpha * (dfp[0] * s[0] + dfp[1] * s[1] + dfp[3] * s[2]) * Timestep2);
		r[3 * count + 1] = M * (f[1] * Timestep + Alpha * (dfp[1] * s[0] + dfp[2] * s[1] + dfp[4] * s[2]) * Timestep2);
		r[3 * count + 2] = M * (f[2] * Timestep + Alpha * (dfp[3] * s[0] + dfp[4] * s[1] + dfp[5] * s[2]) * Timestep2);

		count++;
	}
}