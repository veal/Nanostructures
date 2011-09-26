#ifndef GLOBALCONSTVAR_H
#define GLOBALCONSTVAR_H

#endif // GLOBALCONSTVAR_H

#include "nr3.h"
#include <complex>

typedef complex<Doub> Comp;
const int N_LAT=12, N_Band=10, N_Com = 2, N_En_Pr = 74, Nkp = 200, Nkp2 = 400, MAX_IT = 400, N_En = 200, N_point = 3, N_Alpha = 3, n_Cores = 6;
const Doub EPS_CP = 1e-3;
const Doub Pi = 3.14;
Doub l_CPA;
Comp CI(0.0, 1.0);

Doub(*r)[3][3];
Doub(*nm_1)[N_Com][N_LAT][N_Band];
Doub(*mm_1)[N_Com][N_LAT][N_Band];
Doub(*ni0)[4];
Doub(*mm0)[4];
Doub U_matr[N_Band];
Comp (*Coh_p)[2][N_LAT][N_LAT][N_Band][N_Band];
Comp (*W)[N_LAT][N_LAT][N_LAT][N_Band][N_Band];
Comp (*sigmam)[N_LAT][2][N_Band];
Comp (*sigma_m)[N_LAT][2][N_Band];
Comp (*sigma_coh)[N_LAT];
Comp (*sigma_comp)[18];
Comp (*Hr)[3][N_LAT][N_LAT][N_Band][N_Band];
Comp (*H_k)[N_LAT][N_LAT][N_Band][N_Band];
Comp (*S_k)[N_LAT][N_LAT][N_Band][N_Band];
Comp (*D_k)[N_LAT][N_LAT][N_Alpha][N_Alpha];
Comp (*H_Difr)[3][N_LAT][N_LAT][N_Band][N_Band];
Comp (*Gs)[3][3][N_LAT][N_LAT][N_Band][N_Band];
Comp (*Sig_ef_e)[N_LAT][N_Band][N_Band];
Comp (*Sig_fe)[3][N_LAT][N_LAT][N_Alpha][N_Alpha];
Doub(*Pm)[N_LAT];
Doub(*Cl)[N_LAT];
Doub(*vm)[N_Com][N_LAT][2][N_Band];
Doub ee[N_En];
bool calculateHamiltonian;
bool calculateWamiltonian;
bool calculateGreen;
bool calculateSig_fe;
bool calculate_nm_mm;

Doub En_S;
Doub En_E;
Doub Et;
Doub Et_a;
Doub gmi;
Doub Esmr;
Doub N;
Doub Ef, T, E_N_Start, E_N_End;
int N_f_int, N_start, N_end;
int itmS, it_m;
int mm_stop, l_data;
Doub del_m;
