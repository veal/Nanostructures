#include "sys_param.h"

Doub Pmm(int m1, int m2, int i, int j, Doub Pm[2][N_LAT], Doub Eps) {
    //**************      m1->i ; m2->j ; I->(-m); II->(+m)
    Doub Sign_P;
    //****
    if (Pm[m2][j] == 0.0) {
        return 0.0;
    } else {
        Sign_P = (m1 == m2) ? 1.0 : -1.0;
        return Pm[m1][i] + Eps * Sign_P / Pm[m2][j];
    }
}

Doub Cll(int l1, int l2, int i, int j, Doub Cl[N_Com][N_LAT], Doub Eps_a) {
    //**************      l1->i ; l2->j ; I->(B); II->(A)
    Doub Sign_C;
    //****
    if (Cl[l2][j] == 0.0) {
        return 0.0;
    } else {
        Sign_C = (l1 == l2) ? 1.0 : -1.0;
        return Cl[l1][i] + Eps_a * Sign_C / Cl[l2][j];
    }
}

void Mm_nm(Doub gsm[N_Com][N_LAT][2][N_Band], Doub gs_m[N_Com][N_LAT][2][N_Band], Doub gs_m_1[N_Com][N_LAT][2][N_Band],
        Doub gsm_1[N_Com][N_LAT][2][N_Band], Doub gsm_2[N_Com][N_LAT][2][N_Band], Doub gs_m_2[N_Com][N_LAT][2][N_Band],
        Doub Pm[2][N_LAT], Doub nm[2][N_Com][N_LAT][N_Band], Doub mm[2][N_Com][N_LAT][N_Band],
        Doub ee[N_En], int niE, Doub dE, Doub Ef, Doub T) {
    for (int i1 = 0; i1 < N_LAT; i1++) {
        for (int ib = 0; ib < N_Band; ib++) {
            for (int j = 0; j < N_Com; j++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 2; n++) {
                        Doub gmi = 0.0;
                        Doub gmi_1 = 0.0;
                        Doub gmi_2 = 0.0;
                        //						spin = 1 <=> spin_down(-);  spin = 2 <=> spin_up(+)
                        for (int spin = 0; spin < 2; spin++) {
                            int sign = 1;
                            if (spin == 0 && n == 1) {
                                sign = -1;
                            }
                            if (m == 0) {
                                Doub k = (gs_m[j][i1][spin][ib] > 0) ? gs_m[j][i1][spin][ib] : 0;
                                gmi += sign*k;
                                k = (gs_m_1[j][i1][spin][ib] > 0) ? gs_m_1[j][i1][spin][ib] : 0;
                                gmi_1 += sign*k;
                                k = (gs_m_2[j][i1][spin][ib] > 0) ? gs_m_2[j][i1][spin][ib] : 0;
                                gmi_2 += sign*k;
                            } else {
                                Doub k = (gsm[j][i1][spin][ib] > 0) ? gsm[j][i1][spin][ib] : 0;
                                gmi += sign*k;
                                k = (gsm_1[j][i1][spin][ib] > 0) ? gsm_1[j][i1][spin][ib] : 0;
                                gmi_1 += sign*k;
                                k = (gsm_2[j][i1][spin][ib] > 0) ? gsm_2[j][i1][spin][ib] : 0;
                                gmi_2 += sign*k;
                            }
                        }
                        if (n == 0) {
                            //if (niE >= 2)
                            //{
                            //	Int_Simpson( nm[m][j][i1][ib], niE, gmi_2.r, gmi_1.r, gmi.r, dE);
                            //}
                            if (niE >= 1) {
                                nm[m][j][i1][ib] += (gmi_1 + gmi) * dE / 2;
                                //Int_Simpson( nm[m][j][i1][ib], niE, 0, gmi_1.r, gmi.r, dE);
                            }
                            if (niE == 0) {
                                nm[m][j][i1][ib] += gmi * dE / 2;
                                //Int_Simpson( nm[m][j][i1][ib], niE, 0, 0, gmi.r, dE);
                            }
                        } else {
                            //if (niE >= 2)
                            //{
                            //	Int_Simpson( mm[m][j][i1][ib], niE, gmi_2.r,gmi_1.r, gmi.r, dE);
                            //}
                            if (niE >= 1) {
                                mm[m][j][i1][ib] += (gmi_1 + gmi) * dE / 2;
                                //Int_Simpson( mm[m][j][i1][ib], niE, 0, gmi_1.r, gmi.r, dE);
                            }
                            if (niE == 0) {
                                mm[m][j][i1][ib] += gmi * dE / 2;
                                //Int_Simpson( mm[m][j][i1][ib], niE, 0, 0, gmi.r, dE);
                            }
                        }
                    }
                }
            }
        }
    }
}
