#include "sys_param.h"
Doub l_CPA;
Comp CI(0.0, 1.0);

Doub(*r)[3][3];
Doub* R;
Doub* translationVectors;
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

void allocate() {
    l_CPA = 0;
    r = new Doub[N_LAT][3][3];
    R = new Doub[N_LAT * 3];
    translationVectors = new Doub[DIMENSION * 3];
    nm_1 = new Doub[2][N_Com][N_LAT][N_Band];
    mm_1 = new Doub[2][N_Com][N_LAT][N_Band];
    ni0 = new Doub[2][4];
    mm0 = new Doub[2][4];
    Coh_p = new Comp[N_En][2][N_LAT][N_LAT][N_Band][N_Band];
    W = new Comp[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band];
    sigmam = new Comp[N_Com][N_LAT][2][N_Band];
    sigma_m = new Comp[N_Com][N_LAT][2][N_Band];
    sigma_coh = new Comp[2][N_LAT];
    sigma_comp = new Comp[2][18];
    D_k = new Comp[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha];
    H_Difr = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    Sig_ef_e = new Comp[N_En][N_LAT][N_Band][N_Band];
    Sig_fe = new Comp[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha];
    H_Difr = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    Gs = new Comp[N_En][3][3][N_LAT][N_LAT][N_Band][N_Band];
    Sig_ef_e = new Comp[N_En][N_LAT][N_Band][N_Band];
    Sig_fe = new Comp[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha];
    vm = new Doub[2][N_Com][N_LAT][2][N_Band];
    r = new Doub[N_LAT][3][3];
    Pm = new Doub[2][N_LAT];
    Cl = new Doub[N_Com][N_LAT];
    calculateHamiltonian = true; //Caution! if false S_k won't be calculated!!!
    calculateWamiltonian = true;
    calculateGreen = true;
    calculateSig_fe = false;
    calculate_nm_mm = true;
    En_S = -20.0;
    En_E = 15.0;
    Et = 0;
    Et_a = 0;
    gmi = 0;
    Esmr = 0;
    N = 4.1;
    T = 300;
    itmS = 0;
    it_m = 1;
    mm_stop = 0;
    l_data = 1;
    del_m = 0.1;
}
