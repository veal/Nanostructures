#include "sys_param.h"
void Conductivity(complex<Doub> sigma_coh[2], complex<Doub> sigma_comp[18], complex<Doub> sigmam[N_Com][N_LAT][2][N_Band], complex<Doub> sigma_m[N_Com][N_LAT][2][N_Band],
				  complex<Doub> W[N_Com][N_LAT][N_Band][N_Band], complex<Doub> vm[2][N_Com][N_LAT][2][N_Band], Doub Pm[2][N_LAT], Doub Cl[N_Com][N_LAT],
				  complex<Doub> Hr[3][3][N_LAT][N_LAT][N_Band][N_Band],int Nkp, complex<Doub> H_Difr[3][3][N_LAT][N_LAT][N_Band][N_Band], Doub Eps,
				  Doub Eps_sec, Doub ee[N_En],int niE,Doub E_N_Start,Doub E_N_End,Doub T_K,int spin,Doub dE,
				  Doub Ef,int imm, string Path, complex<Doub> Coh_p[N_En][2][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3])
{     
	Doub E3[N_point] = {0};
	Doub E3_S, E3_E, dE3;
	complex<Doub> c1(0.0, 1.0);
	complex<Doub> (*C_P)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*G_k)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GVdGVG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GV)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*dGVG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*Go_k)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*dG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*VdG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GoV)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*VdGVdG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*VG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GVG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GoVG)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GVGo)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GoVGo)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*dGVGo)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*GoVdGVGo)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*Ho_k)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*H_k)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*H_Dk)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*VGo)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
	complex<Doub> (*I1)[N_LAT][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][N_Band][N_Band];

	complex<Doub> (*Gs)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*G_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*Gos)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*Go_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GVdGVGs)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GoVGs)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GVGs)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*VdGVdGs)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GoVGos)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GVGos)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GVG_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GVGo_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GoVG_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GoVGo_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];
	complex<Doub> (*GoVdGVGos)[3][N_LAT][N_LAT][N_Band][N_Band] = new complex<Doub>[3][3][N_LAT][N_LAT][N_Band][N_Band];

	complex<Doub> (*t)[N_Com][N_LAT][N_Band][N_Band] = new complex<Doub>[3][N_Com][N_LAT][N_Band][N_Band];
	complex<Doub> (*t0)[N_Com][N_LAT][N_Band][N_Band] = new complex<Doub>[3][N_Com][N_LAT][N_Band][N_Band];
	complex<Doub> (*Wm_D)[N_Com][N_LAT][N_Band][N_Band] = new complex<Doub>[2][N_Com][N_LAT][N_Band][N_Band];

	Doub Eps_a, Eps_a_sec, Eps_a1, Eps1;
	complex<Doub> (*I2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x3)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x5)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x6)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x7)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x8)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x9)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x11)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo3)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo7)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo5)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo6)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo8)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo9)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*xo11)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*Td)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*Tdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma1_1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma1_2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma1_3)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma1_4)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma1_5)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_6)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_7)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_8)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_9)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_10)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_11)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o7)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o6)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o8)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o9)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o10)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*x_o11)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T_d)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T_do)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T_n)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T_no)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma1)[N_Band] = new complex<Doub>[N_Band][N_Band];

	complex<Doub> (*summa1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*summa1_1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*Kt)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*tKt)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*KTd)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*tKT2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T2Kt)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T2KT2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*KTn)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*TnKTn)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*summa)[N_Band] = new complex<Doub>[N_Band][N_Band];
	
	complex<Doub> (*summa2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*summa2_1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*oKt)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*t0oKt)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*oKTd)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*t0oKT2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*To2oKt)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*To2oKT2)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*oKTn)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*TnoKoTn)[N_Band] = new complex<Doub>[N_Band][N_Band];

	complex<Doub> (*summa3)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*summa3_1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*Kot0)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*tKot0)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*KoTdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*tKoTdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T2Kot0)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*T2KoTdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*KoTno)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*TnKoTno)[N_Band] = new complex<Doub>[N_Band][N_Band];

	complex<Doub> (*summa4)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*summa4_1)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*oKot0)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*t0oKot0)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*oKoTdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*t0oKoTdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*TdooKot0)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*TdooKoTdo)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*oKoTno)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*TnooKoTno)[N_Band] = new complex<Doub>[N_Band][N_Band];
	complex<Doub> (*sigma)[N_Band] = new complex<Doub>[N_Band][N_Band];

	complex<Doub> (*sigma_cpa)[N_Band] = new complex<Doub>[N_LAT][N_Band];
	complex<Doub> (*sigmaEsm)[N_Com][N_LAT][N_Band] = new complex<Doub>[2][N_Com][N_LAT][N_Band];

	complex<Doub> (*Sig_ef)[2][N_LAT][N_Band][N_Band] = new complex<Doub>[N_En][2][N_LAT][N_Band][N_Band];
	complex<Doub> (*Sig_ef_e_l)[N_LAT][2][N_Band][N_Band] = new complex<Doub>[N_LAT][N_LAT][2][N_Band][N_Band];
	
	complex<Doub> sigma_coh_old[2], sigma_coh_new[2];
	complex<Doub> (*sigmam_old)[N_LAT][N_Band] = new complex<Doub>[N_Com][N_LAT][N_Band];
	complex<Doub> (*sigma_m_old)[N_LAT][N_Band] = new complex<Doub>[N_Com][N_LAT][N_Band];
	complex<Doub> (*sigmam_new)[N_LAT][N_Band] = new complex<Doub>[N_Com][N_LAT][N_Band];
	complex<Doub> (*sigma_m_new)[N_LAT][N_Band] = new complex<Doub>[N_Com][N_LAT][N_Band];

	complex<Doub> (*sigmaEsm_old)[N_Com][N_LAT][N_Band] = new complex<Doub>[2][N_Com][N_LAT][N_Band];
	complex<Doub> (*sigmaEsm_new)[N_Com][N_LAT][N_Band] = new complex<Doub>[2][N_Com][N_LAT][N_Band];

	int le_n, N_S;
	Doub N_start, N_end;
//!********************
//	common/Parameters/ eps_a, eps_a_sec, key_sec_sph  
//	common/Omega/      a
//	Common/Nkp_tot/ Ntot
//	Common/Out_it/  it, f2
//	Common/Sig_ef/		Sig_ef_s
//	Common/Coh_p/		Coh_p_e
//	Common/Temp/		iii
//!********************

//	Sig_ef(:,:,:,:,:) = Sig_ef_s(:,spin,:,:,:,:)
//	Coh_p(:,:,:,:) = Coh_p_e(:,spin,:,:,:)

	if (T_K == 0.0)
	{
		N_start = niE;
		N_end = niE;
	}

	Doub Koef;// = -(0.193e-4)/a	// (e**2.d0)/2*h = 0.193d-4; 1/(a**2) sokratilos s razmernostu prishedshei iz skorosti
	Doub K_int_old = 0;
	Doub K_int_new = 0;

	for (int i = 0; i < N_LAT; i++)
	{
		for (int j = 0; j < N_Band; j++)
		{
			I1[i][i][j][j] = 1.0;
		}
	}
	E3_S = E_N_Start; 
	E3_E = E_N_End;
	dE3 = (E3_E - E3_S)/(N_point-1);
	for (int nE3 = 0; nE3 < N_point; nE3++)
	{
		E3[nE3] = E3_S + dE3*nE3;
	}
	if (E_N_Start == E_N_End)
	{
		cout << "Warning in conduct: E_N_Start=E_N_End" << '\n';
	}
	if (E_N_Start > E_N_End)
	{
		cout << "Warning in conduct: E_N_Start>E_N_End" << '\n';
	}
	for (int N_integr = 0; N_integr < N_point; N_integr++)
	{
		cout << "Number of point in conduct=" << N_integr << '\n';
		for (int i = 0; i < N_En; i++)
		{
			if ((ee[i] - E3[N_integr])*(ee[i+1] - E3[N_integr]) < 0)
			{
				le_n = (fabs(ee[i] - E3[N_integr]) < fabs(ee[i+1] - E3[N_integr]))?i:i+1;
			}
		}
		for (int j = 0; j < N_LAT; j++)
		{
			for (int b1 = 0; b1 < N_Band; b1++)
			{
				for (int b2 = 0; b2 < N_Band; b2++)
				{
					C_P[j][j][b1][b2] = Coh_p[le_n][spin][j][b1][b2];
					Sig_ef_e_l[j][j][0][b1][b2]=Sig_ef[le_n][0][j][b1][b2];
					Sig_ef_e_l[j][j][1][b1][b2]=Sig_ef[le_n][1][j][b1][b2];
				}
			}
		}
		for (int k1 = -Nkp; k1 < Nkp; k1++)
		{
			Doub kx = 2.0*Pi*static_cast<Doub>(k1)/(3*Nkp);
			Hij_k(kx, Hr, H_k, r);
			Hij_k (kx, H_Difr, H_Dk, r);
			complex<Doub> (*Ho_k1)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
			complex<Doub> (*H_k1)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
			complex<Doub> (*G_k1)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
			complex<Doub> (*Go_k1)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
			complex<Doub> (*H_Dk1)[N_LAT*N_Band] = new complex<Doub>[N_LAT*N_Band][N_LAT*N_Band];
			for (int i1 = 0; i1 < N_LAT; i1++)
			{
				for (int i2 = 0; i2 < N_LAT; i2++)
				{
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							Ho_k1[i1*N_Band+b1][i2*N_Band+b2] = ee[le_n]*I1[i1][i2][b1][b2] - (conj(C_P[i1][i2][b1][b2]) + H_k[i1][i2][b1][b2] + 
								Sig_ef_e_l[i1][i2][spin][b1][b2]);
							H_k1[i1*N_Band+b1][i2*N_Band+b2] = ee[le_n]*I1[i1][i2][b1][b2] - (C_P[i1][i2][b1][b2] + H_k[i1][i2][b1][b2] + 
								Sig_ef_e_l[i1][i2][spin][b1][b2]);
							H_Dk1[i1*N_Band+b1][i2*N_Band+b2] = H_Dk[i1][i2][b1][b2];
						}
					}
				}
			}
			matinv(H_k1, G_k1);
			matinv(Ho_k1, Go_k1);
			delete [] Ho_k1;
			delete [] H_k1;
//          G - ������������� ����. �����, Go - ����������� 
//********************VdGVdG****************************************************
			for (int i1 = 0; i1 < N_LAT; i1++)
			{
				for (int i2 = 0; i2 < N_LAT; i2++)
				{
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							dG[i1*N_Band+b1][i2*N_Band+b2] = G_k1[i1*N_Band+b1][i2*N_Band+b2] - Go_k1[i1*N_Band+b1][i2*N_Band+b2];
						}
					}
				}
			}
			matmul(H_Dk1, dG, VdG);
			matmul(VdG, VdG, VdGVdG); // sigm_zz - diagonal component of tensor

			matmul(H_Dk1, G_k1, VG);
			matmul(G_k1, VG, GVG);
			matmul(Go_k1, VG, GoVG);
			matmul(H_Dk1, Go_k1, VGo);
			matmul(G_k1, VGo, GVGo);
			matmul(Go_k1, VGo, GoVGo);

			matmul(dG, VG, dGVG);
			matmul(G_k1, H_Dk1, GV);
			matmul(GV, dGVG, GVdGVG);

			matmul(dG, VGo, dGVGo);
			matmul(Go_k1, H_Dk1, GoV);
			matmul(GoV, dGVGo, GoVdGVGo);
			delete [] H_Dk1;
			delete [] Go_k1;
			delete [] G_k1;
//********************End VdGVdG************************************************
			Doub rok = ((k1 == -Nkp) || (k1 == Nkp-1))?0.5:1;
			for (int n = 0; n < 3; n++)
			{
				for (int n1 = 0; n1 < 3; n1++)
				{
					for (int i = 0; i < N_LAT; i++)
					{
						for (int j = 0; j < N_LAT; j++)
						{
							complex<Doub> p = c1*kx*(r[i][0][n]-r[j][0][n1]);
//********************Gij0m and Gijm0*******************************************
							for (int b1 = 0; b1 < N_Band; b1++)
							{
								for (int b2 = 0; b2 < N_Band; b2++)
								{
									Gs[n][n1][i][j][b1][b2] += G_k[i][j][b1][b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
									Gos[n][n1][i][j][b1][b2] += Go_k[i][j][b1][b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (VdGVdG)ij0m************************************************
									VdGVdGs[n][n1][i][j][b1][b2] += VdGVdG[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (GVG)ij0m***************************************************
									GVGs[n][n1][i][j][b1][b2] += GVG[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (GoVG)ij0m**************************************************
									GoVGs[n][n1][i][j][b1][b2] += GoVG[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (GVGo)ij0m***************************************************
									GVGos[n][n1][i][j][b1][b2] += GVGo[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (GoVGo)ij0m************************************************
									GoVGos[n][n1][i][j][b1][b2] += GoVGo[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (GVdGVG)ij0m************************************************
									GVdGVGs[n][n1][i][j][b1][b2] += GVdGVG[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
//********* ������� (GoVdGVGo)ij0m************************************************
									GoVdGVGos[n][n1][i][j][b1][b2] += GoVdGVGo[i*N_Band+b1][j*N_Band+b2]*rok*exp(p)/static_cast<Doub>(2*Nkp);
								}
							}
						}
					}
				}
			}
		}
//*********************End Green function***************************************
//*********************T-matrix*************************************************
		for (int i = 0; i < N_Band; i++)
		{
			I2[i][i] = 1.0;
			for (int nn = 0; nn < 2; nn++)
			{
				for (int nc = 0; nc < N_LAT; nc++)
				{
					for (int lat = 0; lat < N_LAT; lat++)
					{
						Wm_D[nn][nc][lat][i][i] = vm[nn][nc][lat][spin][i];
					}
				}
			}
		}
//**** set up single-site scattering matrix ***********	
		for (int i = 0; i < N_LAT; i++)
		{
			for (int j = 0; j < N_Com; j++)
			{
				for (int m = 0; m < 2; m++)
				{
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							x1[b1][b2] = W[j][i][b1][b2] + Wm_D[m][j][i][b1][b2] - Coh_p[le_n][spin][i][b1][b2];
						}
					}
					matmul4(x1, Gs[1][1][i][i], x2);
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							x3[b1][b2] = I2[b1][b2] - x2[b1][b2];
						}
					}
					matinv4(x3, x2);
					matmul4(x2, x1, t[m][j][i]);
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							xo1[b1][b2] = W[j][i][b1][b2] + Wm_D[m][j][i][b1][b2] - conj(Coh_p[le_n][spin][i][b1][b2]);
						}
					}
					matmul4(xo1, Gos[1][1][i][i], xo2);
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							xo3[b1][b2] = I2[b1][b2] - xo2[b1][b2];
						}
					}
					matinv4(xo3, xo2);
					matmul4(xo2, xo1, t0[m][j][i]);
				}
			}
		}
//*********************conductivity*********************************************
//	sigma_coh = 0.0
//	sigmaEsm = 0.d0
//	sigmaTemp = 0.d0
		for (int i = 0; i < N_LAT; i++)
		{
			for (int ib = 0; ib < N_Band; ib++)
			{
				sigma_cpa[i][ib] = VdGVdGs[1][1][i][i][ib][ib];
				sigma_coh_new[i] += sigma_cpa[i][ib];
			}
		}
		ofstream file_out("E:\\MyBComp\\Study\\Clean Nanotube\\Progy\\mine\\Program\\sigma_minus.dat");
		for (int i = 0; i < N_LAT; i++)
		{
			for (int ib = 0; ib < N_Band; ib++)
			{
				file_out << imm << "   " << ee[le_n] << spin << i << ib << real(sigma_cpa[i][ib]) << '\n';
			}
		}
		file_out.close();
//******************************************************************************
		for (int m = 0; m < 2; m++)
		{
			for (int il = 0; il < N_LAT; il++)
			{
				for (int i = 0; i < N_Com; i++)
				{
					matmul4(GVdGVGs[1][1][il][il],t[m][i][il], x5);
					matmul4(GoVdGVGos[1][1][il][il],t0[m][i][il], xo5);
					for (int b1 = 0; b1 < N_Band; b1++)
					{
						for (int b2 = 0; b2 < N_Band; b2++)
						{
							x5[b1][b2] *= 2.0;
							xo5[b1][b2] *= -2.0;
							sigma1_1[b1][b2] = x5[b1][b2] + xo5[b1][b2];
						}
					}
					for (int jl = 0; jl < N_LAT; jl++)
					{
					//if (il==jl) then
					//	Eps1 = Eps_sec
					//	Eps_a1 = Eps_a_sec
					//else
					//	Eps1 = Eps
					//	Eps_a1 = Eps_a
					//end if
						for (int m1 = 0; m1 < 2; m1++)
						{
							for (int j = 0; j < N_Com; j++)
							{
								for (int n = 0; n < 3; n++)
								{
									if ((n != 1) || (il != jl))
									{
	//************GVdGVG Tdiag(E+)***************************************************
	//         i <-> il,m; j <-> jl,m1
								matmul4(t[m][i][il], Gs[1][n][il][jl], x1);  
								matmul4(t[m1][j][jl], Gs[n][1][jl][il], x7);			   
								matmul4(x1, x7, x6);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										x2[b1][b2] = I2[b1][b2] - x6[b1][b2];
									}
								}
								matinv4(x2, x8);
								matmul4(x8, x6, x9);
								matmul4(x9, t[m][i][il], Td);
								matmul4(GVdGVGs[1][1][il][il], Td, x11);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										sigma1_2[b1][b2] = 2.0*x11[b1][b2];
									}
								}
	//***********GVdGVG Tdiag(E-)***************************************************
								matmul4(t0[m][i][il], Gos[1][n][il][jl], xo1);  
								matmul4(t0[m1][j][jl], Go_s[n][1][jl][il], xo7);			   
								matmul4(xo1, xo7, xo6);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										xo2[b1][b2] = I2[b1][b2] - xo6[b1][b2];
									}
								}
								matinv4(xo2, xo8);							   
								matmul4(xo8, xo6, xo9);
								matmul4(xo9, t0[m][i][il], Tdo);
								matmul4(GoVdGVGos[1][1][il][il], Tdo, xo11);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										sigma1_3[b1][b2] = -2.0 * xo11[b1][b2];
									}
								}
	//***********GVdGVG Tnediag(E+)*�������� ������� �������************************

								matmul4(t[m1][j][jl], Gs[n][1][jl][il], x_1);  
								matmul4(t[m][i][il], Gs[1][n][il][jl], x_7);			   
								matmul4(x_1, x_7, x_6);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										x_2[b1][b2] = I2[b1][b2] - x_6[b1][b2];
									}
								}
								matinv4(x_2, x_8);							   
								matmul4(x_8, x_1, x_9);
								matmul4(x_9, t[m][i][il], T_n);
								matmul4(GVdGVGs[1][n][il][jl], T_n, x_11);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										sigma1_4[b1][b2] = 2.0*x_11[b1][b2];
									}
								}

	//***********GVdGVG Tnediag(E-)*�������� ������� �������************************

								matmul4(t0[m1][j][jl], Gos[n][1][jl][il], x_o1);  
								matmul4(t0[m][i][il], Gos[1][n][il][jl], x_o7);			   
								matmul4(x_o1, x_o7, x_o6);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										x_o2[b1][b2] = I2[b1][b2] - x_o6[b1][b2];
									}
								}
								matinv4(x_o2, x_o8);							   
								matmul4(x_o8, x_o1, x_o9);
								matmul4(x_o9, t0[m][i][il], T_no);
								matmul4(GoVdGVGos[1][n][il][jl], T_no, x_o11);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										sigma1_5[b1][b2] = -2.0 * x_o11[b1][b2];
									}
								}
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										sigma1[b1][b2] = sigma1_2[b1][b2] + sigma1_3[b1][b2]
										+ sigma1_4[b1][b2] + sigma1_5[b1][b2];
									}
								}

	//******************************************************************************
	//***********GVdGVG Tdiag(E+)*�������� ������� �������**************************

								matmul4(T_n, Gs[1][n][il][jl], x_10);
								matmul4(x_10, t[m1][j][jl], T_d);

	//************GVdGVG Tdiag(E-)*�������� ������� �������**************************

								matmul4(T_no, Gos[1][n][il][jl], x_o10);
								matmul4(x_o10, t0[m1][j][jl], T_do);
	//*******************************************************************************

	//*******************E+ E+*******************************************************

								matmul4(GVG_s[n][1][jl][il], t[m][i][il], Kt);  
								matmul4(t[m1][j][jl], Kt, tKt); 
								matmul4(GVG_s[n][1][jl][il], Td, KTd);  
								matmul4(t[m1][j][jl], KTd, tKT2);
								matmul4(T_d, Kt, T2Kt);  
								matmul4(T_d, KTd, T2KT2);  
								matmul4(GVGs[1][n][il][jl], T_n, KTn);  
								matmul4(T_n, KTn, TnKTn);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										summa1[b1][b2] = tKt[b1][b2] + tKT2[b1][b2] + 
											T2Kt[b1][b2] + T2KT2[b1][b2] + TnKTn[b1][b2];
									}
								}
								matmul4(GVGs[1][n][il][jl], summa1, summa1_1);  

	//*******************E+ E-*******************************************************
								matmul4(GoVG_s[n][1][jl][il], t[m][i][il], oKt);  
								matmul4(t0[m1][j][jl], oKt, t0oKt);  
								matmul4(GoVG_s[n][1][jl][il], Td, oKTd);  
								matmul4(t0[m1][j][jl], oKTd, t0oKT2);
								matmul4(T_do, oKt, To2oKt);  
								matmul4(T_do, oKTd, To2oKT2);  
								matmul4(GoVGs[1][n][il][jl], T_n, oKTn);  
								matmul4(T_no, oKTn, TnoKoTn);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										summa2[b1][b2] = t0oKt[b1][b2] + t0oKT2[b1][b2]
										+ To2oKt[b1][b2] + To2oKT2[b1][b2] + TnoKoTn[b1][b2];
									}
								}
								matmul4(GVGos[1][n][il][jl], summa2, summa2_1);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										summa2_1[b1][b2] *= -1.0;
									}
								}

	//*******************E- E+*******************************************************

								matmul4(GVGos[n][1][jl][il], t0[m][i][il],Kot0);  
								matmul4(t[m1][j][jl], Kot0, tKot0);  
								matmul4(GVGos[n][1][jl][il], Tdo, KoTdo);  
								matmul4(t[m1][j][jl], KoTdo, tKoTdo);
								matmul4(T_d, Kot0, T2Kot0);
								matmul4(T_d, KoTdo, T2KoTdo);
								matmul4(GVGos[1][n][il][jl], T_no, KoTno);
								matmul4(T_n, KoTno, TnKoTno);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										summa3[b1][b2] = tKot0[b1][b2] + tKoTdo[b1][b2]
										+ T2Kot0[b1][b2] + T2KoTdo[b1][b2] + TnKoTno[b1][b2];
									}
								}
								matmul4(GoVGs[1][n][il][jl], summa3, summa3_1);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										summa3_1[b1][b2] *= -1.0;
									}
								}

	//*******************E- E-*******************************************************
								matmul4(GoVGos[n][1][jl][il],t0[m][i][il], oKot0);  
								matmul4(t0[m1][j][jl], oKot0, t0oKot0);
								matmul4(GoVGos[n][1][jl][il], Tdo, oKoTdo);  
								matmul4(t0[m1][j][jl], oKoTdo, t0oKoTdo);
								matmul4(T_do,  oKot0, TdooKot0);
								matmul4(T_do, oKoTdo, TdooKoTdo);
								matmul4(GoVGos[1][n][il][jl], T_no, oKoTno);
								matmul4(T_no, oKoTno, TnooKoTno);
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										summa4[b1][b2] = t0oKot0[b1][b2] + t0oKoTdo[b1][b2]
										+ TdooKot0[b1][b2] + TdooKoTdo[b1][b2] + TnooKoTno[b1][b2];
									}
								}
								matmul4(GoVGos[1][n][il][jl], summa4, summa4_1);
	//*******************************************************************************
								for (int b1 = 0; b1 < N_Band; b1++)
								{
									for (int b2 = 0; b2 < N_Band; b2++)
									{
										sigma2[b1][b2] = summa1_1[b1][b2] + summa2_1[b1][b2]
										+ summa3_1[b1][b2] + summa4_1[b1][b2];
										sigma[b1][b2] = sigma1[b1][b2] + sigma2[b1][b2];
									}
								}
								for (int ib = 0; ib < N_Band; ib++)
								{
									complex<Doub> sigma_ib = sigma[ib][ib]*Pmm(m1,m,jl,il,Pm,Eps1)*Cll(j,i,jl,il,Cl,Eps_a1);
									sigmaEsm_new[m][i][il][ib] += sigma_ib;
								}
									}
								}
							}
						}
					}
					for (int ib = 0; ib < N_Band; ib++)
					{
						complex<Doub> sigma_ib = sigma1_1[ib][ib] +  sigma_cpa[il][ib];
						sigmaEsm_new[m][i][il][ib] += sigma_ib;
					}
				}
			}
			if(m == 0)
			{
				for (int nc = 0; nc < N_Com; nc++)
				{
					for (int nl = 0; nl < N_LAT; nl++)
					{
						for (int nb = 0; nb < N_Band; nb++)
						{
							sigmam_new[nc][nl][nb]  = sigmaEsm_new[m][nc][nl][nb];
						}
					}
				}
			}
			else
			{
				for (int nc = 0; nc < N_Com; nc++)
				{
					for (int nl = 0; nl < N_LAT; nl++)
					{
						for (int nb = 0; nb < N_Band; nb++)
						{
							sigma_m_new[nc][nl][nb] = sigmaEsm_new[m][nc][nl][nb];
						}
					}
				}
			}
		}

		if (T_K == 0.0)
		{
			K_int_new = 1.0;
		}
		else
		{
			K_int_new = 0.5*(-DFun_Fermi(E3[N_integr],Ef,T_K))*dE3;
		}
		for (int nc = 0; nc < N_Com; nc++)
		{
			sigma_coh[nc] += (K_int_old*sigma_coh_old[nc] + K_int_new*sigma_coh_new[nc]);
			sigma_coh_old[nc] = sigma_coh_new[nc];
			sigma_coh_new[nc] = 0.0;
			for (int nl = 0; nl < N_LAT; nl++)
			{
				for (int nb = 0; nb < N_Band; nb++)
				{
					sigmam[nc][nl][spin][nb] += (K_int_old*sigmam_old[nc][nl][nb] + K_int_new*sigmam_new[nc][nl][nb]);
					sigma_m[nc][nl][spin][nb] += (K_int_old*sigma_m_old[nc][nl][nb] + K_int_new*sigma_m_new[nc][nl][nb]);
					sigmaEsm[0][nc][nl][nb] += (K_int_old*sigmaEsm_old[0][nc][nl][nb] + K_int_new*sigmaEsm_new[0][nc][nl][nb]);
					sigmaEsm[1][nc][nl][nb] += (K_int_old*sigmaEsm_old[1][nc][nl][nb] + K_int_new*sigmaEsm_new[1][nc][nl][nb]);
					sigmam_old[nc][nl][nb] = sigmam_new[nc][nl][nb];
					sigma_m_old[nc][nl][nb] = sigma_m_new[nc][nl][nb];
					sigmaEsm_old[0][nc][nl][nb] = sigmaEsm_new[0][nc][nl][nb];
					sigmaEsm_new[0][nc][nl][nb] = 0.0;
					sigmaEsm_old[1][nc][nl][nb] = sigmaEsm_new[1][nc][nl][nb];
					sigmaEsm_new[1][nc][nl][nb] = 0.0;
					sigmam_new[nc][nl][nb] = 0.0;
					sigma_m_new[nc][nl][nb]	= 0.0;
				}
			}
		}
		K_int_old = K_int_new;
	}
	for (int nc = 0; nc < N_Com; nc++)
	{
		for (int nl = 0; nl < N_LAT; nl++)
		{
			for (int nb = 0; nb < N_Band; nb++)
			{
				sigmam[nc][nl][spin][nb]  = Koef*sigmam[nc][nl][spin][nb];
				sigma_m[nc][nl][spin][nb] = Koef*sigma_m[nc][nl][2][nb];
				for (int sp = 0; sp < 2; sp++)
				{
					sigmaEsm[sp][nc][nl][nb] = Koef*sigmaEsm[sp][nc][nl][nb];
				}
			}
		}
	}
	sigma_coh[0] = Koef*sigma_coh[0];
	sigma_coh[1] = Koef*sigma_coh[1];

	for (int il = 0; il < N_LAT; il++)
	{
		for (int ib = 0; ib < N_Band; ib++)
		{
			for (int m = 0; m < 2; m++)
			{
				for (int i = 0; i < N_Com; i++)
				{
					complex<Doub> dd = (1.0/12.0) * Cl[i][il] * Pm[m][il] * sigmaEsm[m][i][il][ib];
					if(m == 0) {sigma_comp[0] += dd;}
					if(m == 1) {sigma_comp[1] += dd;}
					if(i == 0)  {sigma_comp[2] += dd;}
					if(i == 1)  {sigma_comp[3] += dd;}
					if(i == 2)  {sigma_comp[4] += dd;}
					if(i == 3)  {sigma_comp[5] += dd;}
					if(i == 4)  {sigma_comp[6] += dd;}
					if(i == 5)  {sigma_comp[7] += dd;}
					if(i == 6)  {sigma_comp[8] += dd;}
					if(i == 7)  {sigma_comp[9] += dd;}
					if(i == 8)  {sigma_comp[10] += dd;}
					if(i == 9)  {sigma_comp[11] += dd;}
					if(i == 10)  {sigma_comp[12] += dd;}
					if(i == 11)  {sigma_comp[13] += dd;}
					if(ib == 0) {sigma_comp[14] += dd;}
					if(ib <= 3 && ib > 0) {sigma_comp[15] += dd;}
					if(ib <= 8 && ib > 3) {sigma_comp[16] += dd;}
					if(ib == 9)	{sigma_comp[17] += dd;}
				}
			}
		}
	}
}