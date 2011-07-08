#include "sys_param.h"
void D_matrix(complex<Doub> D[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha],complex<Doub> Sig_fe[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha],
			  Doub r[N_LAT][3][3])
{
	complex<Doub> (*F)[N_Alpha][3][N_LAT][N_LAT] = new complex<Doub>[N_Alpha][N_Alpha][3][N_LAT][N_LAT];
	complex<Doub> (*D_k)[N_LAT*N_Alpha] = new complex<Doub>[N_LAT*N_Alpha][N_LAT*N_Alpha];
	Doub a = 3.0;
	Doub koef = 547.5133;

	for (int i = 0; i < N_Alpha; i++)
	{
		for (int j = 0; j < N_Alpha; j++)
		{
			Doub is_diag = (i==j)?1:0;
	  		for (int l = 0; l < 3; l++)
			{
	   			for (int i2 = 0; i2 < N_LAT; i2++)
				{
	    			for (int j2 = 0; j2 < N_LAT; j2++)
					{
						if ((i2 != j2) || (l != 1))
						{
		Doub dist = (r[i2][0][1] - r[j2][0][l])*(r[i2][0][1] - r[j2][0][l]) +(r[i2][1][1]-r[j2][1][l])*(r[i2][1][1]-r[j2][1][l]) +
			(r[i2][2][1]-r[j2][2][l])*(r[i2][2][1]-r[j2][2][l]);
		F[i][j][l][i2][j2] = -2*(3*(r[i2][i][1]-r[j2][i][l])*(r[i2][j][1]-r[j2][j][l])- dist*is_diag)/(19.342*sqrt(dist*dist*dist*dist*dist))
			+ real(Sig_fe[1][l][i2][j2][i][j]);
		F[i][j][1][i2][i2] -= F[i][j][l][i2][j2];
						}
					}
				}
			}
		}
	}
	for (int k11 = -Nkp; k11 < Nkp; k11++)
	{
		Doub k1 = 2*Pi*k11/(3.0*Nkp);
		//for (int i = 0; i < N_LAT; i++)
		//{
		//	for (int j = 0; j < N_LAT; j++)
		//	{
		//		for (int i1 = 0; i1 < N_Alpha; i1++)
		//		{
		//			for (int j1 = 0; j1 < N_Alpha; j1++)
		//			{
		//				if (j*N_Alpha+j1 >= i*N_Alpha+i1)
		//				{
		//					for (int l = 0; l < 3; l++)
		//					{
		//						D_k[i*N_Alpha+i1][j*N_Alpha+j1] += (F[i1][j1][l][i][j])*cdexp(CI*k1*(r[j][0][l]-r[i][0][1]));
		//					}
		//				}
		//			}
		//		}
		//	}
		//}
		//for (int i = 0; i < N_LAT; i++)
		//{
		//	for (int j = 0; j < N_LAT; j++)
		//	{
		//		for (int i1 = 0; i1 < N_Alpha; i1++)
		//		{
		//			for (int j1 = 0; j1 < N_Alpha; j1++)
		//			{
		//				if (j*N_Alpha+j1 <= i*N_Alpha+i1)
		//				{
		//					D_k[i*N_Alpha+i1][j*N_Alpha+j1] = conj(D_k[j*N_Alpha+j1][i*N_Alpha+i1]);
		//				}
		//			}
		//		}
		//	}
		//}
		for (int i = 0; i < N_LAT; i++)
		{
			for (int j = 0; j < N_LAT; j++)
			{
				for (int i1 = 0; i1 < N_Alpha; i1++)
				{
					for (int j1 = 0; j1 < N_Alpha; j1++)
					{
						for (int l = 0; l < 3; l++)
						{
							D[k11+Nkp][i][j][i1][j1] += (F[i1][j1][l][i][j])*exp(CI*k1*(r[j][0][l]-r[i][0][1]));
						}
					}
				}
			}
		}
	}
	delete [] D_k;
	delete [] F;
}
void G_phonon(complex<Doub> G_p[N_LAT][N_LAT][N_Alpha][N_Alpha],complex<Doub> D_k[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha], Doub E_r,
			  complex<Doub> Coh_p[N_LAT][N_Alpha][N_Alpha], Doub r[N_LAT][3][3])
{
	complex<Doub> (*Im_d)[N_Alpha] = new complex<Doub>[N_Alpha][N_Alpha];
	complex<Doub> (*C_P)[N_LAT][N_Alpha][N_Alpha] = new complex<Doub>[N_LAT][N_LAT][N_Alpha][N_Alpha];
	complex<Doub> (*I1)[N_LAT][N_Alpha][N_Alpha] = new complex<Doub>[N_LAT][N_LAT][N_Alpha][N_Alpha];

	for (int i = 0; i < N_LAT; i++)
	{
		for (int i_1 = 0; i_1 < N_Alpha; i_1++)
		{
			I1[i][i][i_1][i_1] = 1.0;
		}
	}
	for (int i = 0; i < N_LAT; i++)
	{
		for (int p1 = 0; p1 < N_Alpha; p1++)
		{
			for (int p2 = 0; p2 < N_Alpha; p2++)
			{
				C_P[i][i][p1][p2]=Coh_p[i][p1][p2];
			}
		}
	}
	Doub koef = 547.5133;
	Doub E = koef*E_r*E_r*19.93826;
	for (int k1 = -Nkp; k1 < Nkp; k1++)
	{
		Doub k11 = Pi*k1/(3.0*Nkp);
		Doub Ro_k = 1.0;
		complex<Doub> (*G_k_tem)[N_LAT*N_Alpha] = new complex<Doub>[N_LAT*N_Alpha][N_LAT*N_Alpha];
		complex<Doub> (*D_k_tem)[N_LAT*N_Alpha] = new complex<Doub>[N_LAT*N_Alpha][N_LAT*N_Alpha];
		for (int i = 0; i < N_LAT; i++)
		{
			for (int is = 0; is < N_LAT; is++)
			{
				for (int a = 0; a < N_Alpha; a++)
				{
					for (int as = 0; as < N_Alpha; as++)
					{
						D_k_tem[i*N_Alpha+a][is*N_Alpha+as] = E*I1[i][is][a][as] - 
							D_k[k1+Nkp][i][is][a][as] - C_P[i][is][a][as];
					}
				}
			}
		}
		matinv36(D_k_tem, G_k_tem);
		if ((k1 == -Nkp) || (k1 == Nkp-1)) {Ro_k=0.5;}
		for (int i = 0; i < N_LAT; i++)
		{
			for (int i11 = 0; i11 < N_Alpha; i11++)
			{
				for (int j1 = 0; j1 < N_Alpha; j1++)
				{
					G_p[i][i][i11][j1] += Ro_k*G_k_tem[i*N_Alpha+i11][i*N_Alpha+j1]/(2.0*Nkp);
				}
			}
		}
		delete [] D_k_tem;
		delete [] G_k_tem;
	}
	delete [] Im_d;
	delete [] C_P;
	delete [] I1;
}
void G_phonon2(complex<Doub> G_p[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha], complex<Doub> D_k[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha],
			   Doub E_r, complex<Doub> Coh_p[N_LAT][N_Alpha][N_Alpha], Doub r[N_LAT][3][3])
{
	complex<Doub> (*Im_d)[N_Alpha] = new complex<Doub>[N_Alpha][N_Alpha];
	complex<Doub> (*C_P)[N_LAT][N_Alpha][N_Alpha] = new complex<Doub>[N_LAT][N_LAT][N_Alpha][N_Alpha];
	complex<Doub> (*I1)[N_LAT][N_Alpha][N_Alpha] = new complex<Doub>[N_LAT][N_LAT][N_Alpha][N_Alpha];

	for (int i = 0; i < N_LAT; i++)
	{
		for (int i_1 = 0; i_1 < N_Alpha; i_1++)
		{
			I1[i][i][i_1][i_1] = 1.0;
		}
	}
	for (int i = 0; i < N_LAT; i++)
	{
		for (int p1 = 0; p1 < N_Alpha; p1++)
		{
			for (int p2 = 0; p2 < N_Alpha; p2++)
			{
				C_P[i][i][p1][p2]=Coh_p[i][p1][p2];
			}
		}
	}

	Doub koef = 547.5133;
	Doub E = koef*E_r*E_r*19.93826;
	for (int k1 = -Nkp; k1 < Nkp; k1++)
	{
		Doub k11 = Pi*k1/(3.0*Nkp);
		Doub Ro_k = 1.0;
		complex<Doub> (*G_k_tem)[N_LAT*N_Alpha] = new complex<Doub>[N_LAT*N_Alpha][N_LAT*N_Alpha];
		complex<Doub> (*D_k_tem)[N_LAT*N_Alpha] = new complex<Doub>[N_LAT*N_Alpha][N_LAT*N_Alpha];

		for (int i = 0; i < N_LAT; i++)
		{
			for (int j = 0; j < N_LAT; j++)
			{
				for (int a = 0; a < N_Alpha; a++)
				{
					for (int as = 0; as < N_Alpha; as++)
					{
						D_k_tem[i*N_Alpha+a][j*N_Alpha+as] = E*I1[i][j][a][as] - D_k[k1+Nkp][i][j][a][as] - C_P[i][j][a][as];
					}
				}
			}
		}

		matinv36(D_k_tem, G_k_tem);
		
		if ((k1 == -Nkp) || (k1 == Nkp-1)) {Ro_k=0.5;}
		for (int n = 0; n < 3; n++)
		{
			for (int ns = 0; ns < 3; ns++)
			{
				for (int i = 0; i < N_LAT; i++)
				{
					for (int j = 0; j < N_LAT; j++)
					{
						for (int i11 = 0; i11 < N_Alpha; i11++)
						{
							for (int j1 = 0; j1 < N_Alpha; j1++)
							{
								G_p[n][ns][i][j][i11][j1] += Ro_k*G_k_tem[i*N_Alpha+i11][j*N_Alpha+j1]*exp(CI*k11*(r[i][0][n]-r[j][0][ns]))/(2.0*Nkp);
							}
						}
					}
				}
			}
		}
		delete [] D_k_tem;
		delete [] G_k_tem;
	}
}
