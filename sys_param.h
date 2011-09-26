#ifndef SYS_PARAMS_H
#define	SYS_PARAMS_H

#include "nr3.h"
#include "PWF.h"

#include <stdio.h>
#include <pthread.h>
#include <fstream>
#include <string>
#include <iostream>
#include <complex>

#include "W_matrix.h"

class Hop_integrals;
typedef complex<Doub> Comp;
using namespace std;
const int N_LAT=12, N_Band=10, N_Com = 2, N_En_Pr = 74, Nkp = 200, Nkp2 = 400, MAX_IT = 400, N_En = 200, N_point = 3, N_Alpha = 3, n_Cores = 6;
const Doub EPS_CP = 1e-3;
extern Comp CI;
const Doub Pi = 3.14;
extern Doub l_CPA;

const bool VERBOSE = false;

extern Doub (*r)[3][3];
extern Doub (*nm_1)[N_Com][N_LAT][N_Band];
extern Doub (*mm_1)[N_Com][N_LAT][N_Band];
extern Doub (*ni0)[4];
extern Doub (*mm0)[4];
extern Doub U_matr[N_Band];
extern Comp (*Coh_p)[2][N_LAT][N_LAT][N_Band][N_Band];
extern Comp (*W)[N_LAT][N_LAT][N_LAT][N_Band][N_Band];
extern Comp (*sigmam)[N_LAT][2][N_Band];
extern Comp (*sigma_m)[N_LAT][2][N_Band];
extern Comp (*sigma_coh)[N_LAT];
extern Comp (*sigma_comp)[18];
extern Comp (*H_k)[N_LAT][N_LAT][N_Band][N_Band];
extern Comp (*S_k)[N_LAT][N_LAT][N_Band][N_Band];
extern Comp (*D_k)[N_LAT][N_LAT][N_Alpha][N_Alpha];
extern Comp (*H_Difr)[3][N_LAT][N_LAT][N_Band][N_Band];
extern Comp (*Gs)[3][3][N_LAT][N_LAT][N_Band][N_Band];
extern Comp (*Sig_ef_e)[N_LAT][N_Band][N_Band];
extern Comp (*Sig_fe)[3][N_LAT][N_LAT][N_Alpha][N_Alpha];

extern Doub (*vm)[N_Com][N_LAT][2][N_Band];

extern Doub (*Pm)[N_LAT];
extern Doub (*Cl)[N_LAT];
extern Doub ee[N_En];

extern bool calculateHamiltonian;
extern bool calculateWamiltonian;
extern bool calculateGreen;
extern bool calculateSig_fe;
extern bool calculate_nm_mm;

extern Doub En_S;
extern Doub En_E;
extern Doub Et;
extern Doub Et_a;
extern Doub gmi;
extern Doub Esmr;
extern Doub N;
extern Doub Ef, T, E_N_Start, E_N_End;
extern int N_f_int, N_start, N_end;
extern int itmS, it_m;
extern int mm_stop, l_data;
extern Doub del_m;
extern pthread_mutex_t job_queue;

void CP_Calculation(Comp Coh_p[N_LAT][N_LAT][N_Band][N_Band], Doub E, int niE, Comp W[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band],
        Doub vm[2][N_Com][N_LAT][2][N_Band], Doub Pm[2][N_LAT], Doub Cl[N_Com][N_LAT],
        Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3], int spin, double *f2);
void Eigen_values(Comp S[N_LAT][N_LAT][N_Band][N_Band], double*,  bool isOutputNeeded);
void matmul(Comp A[N_LAT*N_Band][N_LAT*N_Band], Comp B[N_LAT*N_Band][N_LAT*N_Band],
			Comp res[N_LAT*N_Band][N_LAT*N_Band]);
void matmul3(Comp A[N_Alpha][N_Alpha], Comp B[N_Alpha][N_Alpha], Comp res[N_Alpha][N_Alpha]);
void matmul36(Comp A[N_LAT*N_Alpha][N_LAT*N_Alpha], Comp B[N_LAT*N_Alpha][N_LAT*N_Alpha],
			  Comp res[N_LAT*N_Alpha][N_LAT*N_Alpha]);
void matmul4(Comp A[N_Band][N_Band], Comp B[N_Band][N_Band], Comp res[N_Band][N_Band]);
void matinv(Comp B[N_LAT*N_Band][N_LAT*N_Band],Comp res[N_LAT*N_Band][N_LAT*N_Band]);
void matinv36(Comp B[N_LAT*N_Alpha][N_LAT*N_Alpha],Comp res[N_LAT*N_Alpha][N_LAT*N_Alpha]);
void matinv4(Comp B[N_Band][N_Band],Comp res[N_Band][N_Band]);
void matinv3(Comp B[N_Alpha][N_Alpha],Comp res[N_Alpha][N_Alpha]);
void Trans_Hk(Comp S[N_LAT][N_LAT][N_Band][N_Band], Comp H[N_LAT][N_LAT][N_Band][N_Band],
			  Comp H1[N_LAT][N_LAT][N_Band][N_Band], bool);
void checkIfMatrixIsHermitian(Comp[N_LAT][N_LAT][N_Band][N_Band], Doub ACCURACY);
void G_phonon(Comp G_p[N_LAT][N_LAT][N_Alpha][N_Alpha],Comp D_k[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha], Doub E_r,
			  Comp Coh_p[N_LAT][N_Alpha][N_Alpha], Doub r[N_LAT][3][3]);
Doub dI(Doub y_2, Doub y_1, Doub y, Doub h);
void Calculate_Wamiltonian(PWF*, PWF*);
struct parameters
{
	int num;
	Comp (*Gs)[N_LAT][N_LAT][N_Band][N_Band];
	Comp (*Hk)[N_LAT][N_LAT][N_Band][N_Band];
	Doub Energy;
	Comp (*CP)[N_LAT][N_Band][N_Band];
};

struct parameters2
{
	int num;
	Comp (*Gs)[3][3][N_LAT][N_LAT][N_Band][N_Band];
	Comp (*Hk)[N_LAT][N_LAT][N_Band][N_Band];
	Doub Energy;
	Comp (*CP)[N_LAT][N_Band][N_Band];
};

struct Sig_fe_param
{
	int num;
	int niE_S;
	int niE_E;
	Doub (*V_f)[N_Alpha][N_Band][N_Band];
	Comp (*Sigma)[3][N_LAT][N_LAT][N_Alpha][N_Alpha];
	Comp (*G_s)[3][3][N_LAT][N_LAT][N_Band][N_Band];
	Doub T;
	Doub Ef;
	Doub (*Cl)[N_LAT];
	Doub* ee;
};
void* CalculateGreen( void* );
void* CalculateGreenFull( void* );
void* Sigma_fe_thread( void* );
void Speed(Doub kx, Doub V_p[N_Com][N_En_Pr], Comp H[N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]);
void Set_coord(Doub[N_LAT][3][3]);

void Read_Energy_Integrals(Doub W[N_Com][N_En_Pr], Doub V[N_Com][N_En_Pr], Doub S[N_Com][N_En_Pr]);
void Read_Coulomb_integral(Doub U_matr[N_Band]);
void Hij_k (Doub, Comp[3][3][N_LAT][N_LAT][N_Band][N_Band],
			Comp[N_LAT][N_LAT][N_Band][N_Band], Doub[N_LAT][3][3]);
void cd_Invert(const Comp[N_LAT*N_Band][N_LAT*N_Band], Comp[N_LAT*N_Band][N_LAT*N_Band]);
void cd_Invert(const Comp**, Comp**, int dimension);
void G_site(Comp[N_LAT][N_LAT][N_Band][N_Band], Doub, Comp[N_LAT][N_LAT][N_Band][N_Band],
			 Comp[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub[N_LAT][3][3]);
void G_site_all(Comp Gs[3][3][N_LAT][N_LAT][N_Band][N_Band], Doub E, int spin, Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band],
				Comp[2][N_LAT][N_LAT][N_Band][N_Band], Comp[N_LAT][N_Band][N_Band], Doub[N_LAT][3][3]);
void Scat_Pr(Doub[2][4], Doub[2][4], Doub[2][N_LAT], Doub[N_Com][N_LAT]);
void Read_input_file(Doub[2][4], Doub[2][4], Doub[2][N_Com][N_LAT][N_Band], Doub[2][N_Com][N_LAT][N_Band]);
void Density(Comp G_s[N_En][3][3][N_LAT][N_LAT][N_Band][N_Band], Doub g_coh[N_LAT], Doub g_cpa[N_LAT][N_Band],
        Doub g_comp[8], Doub gm[N_Com][N_LAT][2][N_Band], Doub g_m[N_Com][N_LAT][2][N_Band], int niE,
        int spin, double E, Comp W[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band], Doub vm[2][N_Com][N_LAT][2][N_Band],
        double Pm[2][N_LAT], double Cl[N_Com][N_LAT], double Eps, Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band],
        Comp Coh_p_e[2][N_LAT][N_LAT][N_Band][N_Band], double r[N_LAT][3][3]);
void Conductivity(Comp sigma_coh[2], Comp sigma_comp[18], Comp sigmam[N_Com][N_LAT][2][N_Band], Comp sigma_m[N_Com][N_LAT][2][N_Band],
				  Comp W[N_Com][N_LAT][N_Band][N_Band], Comp vm[2][N_Com][N_LAT][2][N_Band], Doub Pm[2][N_LAT], Doub Cl[N_Com][N_LAT],
				  Comp Hr[3][3][N_LAT][N_LAT][N_Band][N_Band],int Nkp, Comp H_Difr[3][3][N_LAT][N_LAT][N_Band][N_Band], Doub Eps,
				  Doub Eps_sec, Doub ee[N_En],int niE,Doub E_N_start,Doub E_N_end,Doub T_K,int spin,Doub dE,
				  Doub Ef,int imm, Comp Coh_p[N_En][2][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]);
Doub DFun_Fermi(Doub E, Doub Ef, Doub kT);
Doub Cll(int, int, int, int, Doub[N_Com][N_LAT], Doub);
Doub Pmm(int, int, int, int, Doub[2][N_LAT], Doub);
void Int_Simpson(Doub&, int, Doub, Doub, Doub, Doub);
void Sigma_fe(Doub Ef, Doub T, Comp Gs[N_En][3][3][N_LAT][N_LAT][N_Band][N_Band],
			  Comp Sig_fe[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha], Doub ee[N_En], Doub Cl[N_Com][N_LAT]);
void D_matrix(Comp D[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha],Comp Sig_fe[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha],
			  Doub r[N_LAT][3][3]);
Doub Fun_Fermi(Doub E, Doub Ef, Doub T);
void Mm_nm(Doub gsm[N_Com][N_LAT][2][N_Band], Doub gs_m[N_Com][N_LAT][2][N_Band], Doub gs_m_1[N_Com][N_LAT][2][N_Band],
		   Doub gsm_1[N_Com][N_LAT][2][N_Band], Doub gsm_2[N_Com][N_LAT][2][N_Band], Doub gs_m_2[N_Com][N_LAT][2][N_Band],
		   Doub Pm[2][N_LAT], Doub nm[2][N_Com][N_LAT][N_Band], Doub mm[2][N_Com][N_LAT][N_Band],
		   Doub ee[N_En], int niE, Doub dE, Doub Ef, Doub T);
extern void allocate();

#endif
