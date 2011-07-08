#include "sys_param.h"

void CP_Calculation(Comp Coh_p[N_LAT][N_LAT][N_Band][N_Band], Doub E, int niE, Comp W[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band],
        Doub vm[2][N_Com][N_LAT][2][N_Band], Doub Pm[2][N_LAT], Doub Cl[N_Com][N_LAT],
        Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3], int spin, double *f2) {

    Comp (*dcp)[N_LAT*N_Band][N_LAT*N_Band] = new Comp[N_LAT][N_LAT*N_Band][N_LAT*N_Band];
    Comp (*Wab)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*G_s)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*cp1)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*cp2)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*I1)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*Im_d)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*x1)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*x2)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*x3)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*Wm_D)[2][N_LAT][N_Band][N_Band] = new Comp[2][2][N_LAT][N_Band][N_Band];
    bool isCalculated = false;

    for (int i = 0; i < N_LAT; i++) {
        for (int i1 = 0; i1 < N_Band; i1++) {
            I1[i][i][i1][i1] = 1.0;
        }
    }

    for (int i = 0; i < N_Band; i++) {
        for (int n = 0; n < 2; n++) {
            for (int n1 = 0; n1 < 2; n1++) {
                for (int n2 = 0; n2 < N_LAT; n2++) {
                    Wm_D[n][n1][n2][i][i] = vm[n][n1][n2][spin][i];
                }
            }
        }
    }
    for (int i = 0; i < N_Band; i++) {
        Im_d[i][i] = Comp(0.0, -1.0e-3);
    }
    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_Com; j++) {
            for (int m = 0; m < 2; m++) {
                for (int i1 = 0; i1 < N_LAT; i1++) {
                    for (int j1 = 0; j1 < N_Band; j1++) {
                        for (int j2 = 0; j2 < N_Band; j2++) {
                            Wab[i][i1][j1][j2] += (real(W[j][i][i1][i1][j1][j2]) + Wm_D[m][j][i][j1][j2]) * Pm[m][i] * Cl[j][i];
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < N_LAT; i++) {
        for (int i1 = 0; i1 < N_LAT; i1++) {
            for (int j1 = 0; j1 < N_Band; j1++) {
                for (int j2 = 0; j2 < N_Band; j2++) {
                    Coh_p[i][i1][j1][j2] = Wab[i][i1][j1][j2] + Im_d[j1][j2];
                }
            }
        }
    }

    for (int it = 0; it < MAX_IT; it++) {
        l_CPA = 1;
        *f2 = 0.0;
        for (int i = 0; i < N_LAT; i++) {
            for (int i1 = 0; i1 < N_LAT; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    for (int j2 = 0; j2 < N_Band; j2++) {
                        Coh_p[i][i1][j1][j2] = Comp (real(Coh_p[i][i1][j1][j2]), -abs(imag(Coh_p[i][i1][j1][j2])));
                    }
                }
            }
        }
//        if (!isCalculated) {
            G_site(G_s, E, Coh_p, Hr, r);
//            isCalculated = true;
//        }
        
        Comp (*G_temp)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
        for (int i1 = 0; i1 < N_LAT; i1++) {
            for (int i2 = 0; i2 < N_LAT; i2++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    for (int j2 = 0; j2 < N_Band; j2++) {
                        G_temp[i1*N_Band+j1][i2*N_Band+j2] = G_s[i1][i2][j1][j2];
                    }
                }
            }
        }
        for (int i = 0; i < N_LAT; i++) {
            for (int j1 = 0; j1 < N_Band; j1++) {
                for (int j2 = 0; j2 < N_Band; j2++) {
                    cp1[j1][j2] = 0.0;
                    cp2[j1][j2] = 0.0;
                }
            }
            for (int j = 0; j < N_Com; j++) {
                for (int m = 0; m < 2; m++) {
                    for (int i1 = 0; i1 < N_LAT; i1++) {
                        for (int j1 = 0; j1 < N_Band; j1++) {
                            for (int j2 = 0; j2 < N_Band; j2++) {
                                x1[i1*N_Band+j1][i1*N_Band+j2] = real(W[j][i][i1][i1][j1][j2]) + Wm_D[m][j][i][j1][j2] - real(Coh_p[i][i1][j1][j2]);
                            }
                        }
                    }
                    matmul(x1, G_temp, x3);
                    for (int i1 = 0; i1 < N_LAT; i1++) {
                        for (int j1 = 0; j1 < N_Band; j1++) {
                            for (int j2 = 0; j2 < N_Band; j2++) {
                                x2[i1*N_Band+j1][i1*N_Band+j2] = I1[i1][i1][j1][j2] - x3[i1*N_Band+j1][i1*N_Band+j2];
                            }
                        }
                    }
                    cd_Invert((const Comp**) x2, (Comp**) x3, N_LAT*N_Band);
//                    matinv(x2, x3);
                    matmul(x3, x1, x2);
                    for (int i1 = 0; i1 < N_LAT; i1++) {
                        for (int i2 = 0; i2 < N_LAT; i2++) {
                            for (int j1 = 0; j1 < N_Band; j1++) {
                                for (int j2 = 0; j2 < N_Band; j2++) {
                                    cp1[i1*N_Band+j1][i2*N_Band+j2] += Pm[m][i] * Cl[j][i] * x2[i1*N_Band+j1][i2*N_Band+j2];
                                    cp2[i1*N_Band+j1][i2*N_Band+j2] += Pm[m][i] * Cl[j][i] * x3[i1*N_Band+j1][i2*N_Band+j2];
                                }
                            }
                        }
                    }
                }
            }
//            for (int i1 = 0; i1 < N_LAT; i1++) {
//                for (int j1 = 0; j1 < N_Band; j1++) {
//                    double rratio = abs(Coh_p[i][i1][j1][j1] - cp1[i1*N_Band+j1][i1*N_Band+j1])/
//                            abs(Coh_p[i][i1][j1][j1]);
//                    if (rratio > 0.0005)
//                        cout << i1 << "   " << j1 << "   " << cp1[i1*N_Band+j1][i1*N_Band+j1] << "   " << Coh_p[i][i1][j1][j1]
//                                << "   " << rratio << '\n';
//                }
//            }
            cd_Invert((const Comp**) cp2, (Comp**) x3, N_LAT*N_Band);
//            matinv(cp2, x3);
            matmul(x3, cp1, dcp[i]);
            
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j = 0; j < N_Band; j++) {
                    Doub f1 = abs(dcp[i][i1*N_Band+j][i1*N_Band+j])/abs(Coh_p[i][i1][j][j]); //.r*dcp[i][j][j].r+dcp[i][j][j].i*dcp[i][j][j].i);
                    if (f1 > EPS_CP)//.r*Coh_p[i][j][j].r+Coh_p[i][j][j].i*Coh_p[i][j][j].i))
                    {
                        l_CPA = 0;
                    }
                    if (f1 > (*f2)) {
                        (*f2) = f1;
                    }
                }
            }
        }

        for (int i = 0; i < N_LAT; i++) {
            for (int i1 = 0; i1 < N_LAT; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    for (int j2 = 0; j2 < N_Band; j2++) {
                        Coh_p[i][i1][j1][j2] += dcp[i][i1*N_Band+j1][i1*N_Band+j2];
                        Coh_p[i][i1][j1][j2] = Comp (real(Coh_p[i][i1][j1][j2]),
                                -abs(imag(Coh_p[i][i1][j1][j2])));
                    }
                }
            }
        }
        delete[] G_temp;
        if (l_CPA == 1 && it > 0) {
            delete [] dcp;
            delete [] Wab;
            delete [] G_s;
            delete [] cp1;
            delete [] cp2;
            delete [] I1;
            delete [] Im_d;
            delete [] x1;
            delete [] x2;
            delete [] x3;
            delete [] Wm_D;
            return;
        }
        cout << "CP differs from average t in " << *f2 << " times" << '\n';
    }
    if ((*f2) <= 1.0e-1) {
        l_CPA = 1;
    }
    delete [] dcp;
    delete [] Wab;
    delete [] G_s;
    delete [] cp1;
    delete [] cp2;
    delete [] I1;
    delete [] Im_d;
    delete [] x1;
    delete [] x2;
    delete [] x3;
    delete [] Wm_D;
}

void C_Ph_calculation(Comp Coh_p[N_LAT][N_Alpha][N_Alpha], Doub E, int niE, Doub Cl[N_Com][N_LAT],
        Comp Dk[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha], Doub r[N_LAT][3][3],
        Comp Sig_fe[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha]) {
    int N_E = 360;
    Doub E_s = -0.02, E_e = 0.02;
    int Nkph = 200;
    Doub EPS_CPh = 1e-5;
    complex<Doub > (*dcp)[N_Alpha][N_Alpha] = new Comp[N_LAT][N_Alpha][N_Alpha];
    complex<Doub > (*W)[N_LAT][N_Alpha][N_Alpha] = new Comp[N_Com][N_LAT][N_Alpha][N_Alpha];
    complex<Doub > (*Wab)[N_Alpha][N_Alpha] = new Comp[N_LAT][N_Alpha][N_Alpha];
    complex<Doub > (*G_s)[N_LAT][N_Alpha][N_Alpha] = new Comp[N_LAT][N_LAT][N_Alpha][N_Alpha];
    complex<Doub > (*cp1)[N_Alpha] = new Comp[N_Band][N_Alpha];
    complex<Doub > (*cp2)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
    complex<Doub > (*I1)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
    complex<Doub > (*Im_d)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
    complex<Doub > (*x1)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
    complex<Doub > (*x2)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
    complex<Doub > (*x3)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
    Doub f2 = 0.0;
    l_CPA = 0;

    for (int i = 0; i < N_Band; i++) {
        I1[i][i] = 1.0;
        Im_d[i][i] = (E < 0) ? complex<Doub > (0.0, 1.0e-7) : complex<Doub > (0.0, -1.0e-7);
    }

    //**** Define start value for coherent potential ******
    Doub koef = 547.5133;
    for (int i = 0; i < N_Alpha; i++) {
        for (int il = 0; il < N_LAT; il++) {
            W[0][il][i][i] = 0;
            W[1][il][i][i] = -koef * E * E * 3.312862;
        }
    }
    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_Com; j++) {
            for (int t = 0; t < N_Alpha; t++) {
                for (int t1 = 0; t1 < N_Alpha; t1++) {
                    Wab[i][t][t1] += Cl[j][i] * W[j][i][t][t1];
                }
            }
        }
        for (int t = 0; t < N_Alpha; t++) {
            for (int t1 = 0; t1 < N_Alpha; t1++) {
                Coh_p[i][t][t1] = Wab[i][t][t1] + Im_d[t][t1];
            }
        }
    }
    for (int it = 0; it < MAX_IT; it++) {
        l_CPA = 1;
        for (int i = 0; i < N_LAT; i++) {
            for (int j1 = 0; j1 < N_Alpha; j1++) {
                for (int j2 = 0; j2 < N_Alpha; j2++) {
                    Coh_p[i][j1][j2] = (E < 0) ? complex<Doub > (real(Coh_p[i][j1][j2]), fabs(imag(Coh_p[i][j1][j2]))) :
                            complex<Doub > (real(Coh_p[i][j1][j2]), -fabs(imag(Coh_p[i][j1][j2])));
                }
            }
        }
        G_phonon(G_s, Dk, E, Coh_p, r);
        //*******  third  variant  of  Soven  equation ********
        for (int i = 0; i < N_LAT; i++) {
            for (int j1 = 0; j1 < N_Alpha; j1++) {
                for (int j2 = 0; j2 < N_Alpha; j2++) {
                    cp1[j1][j2] = 0.0;
                    cp2[j1][j2] = 0.0;
                }
            }
            for (int j = 0; j < N_Com; j++) {
                for (int b1 = 0; b1 < N_Alpha; b1++) {
                    for (int b2 = 0; b2 < N_Alpha; b2++) {
                        x1[b1][b2] = W[j][i][b1][b2] - Coh_p[i][b1][b2];
                    }
                }
                matmul3(x1, G_s[i][i], x3);
                for (int j1 = 0; j1 < N_Alpha; j1++) {
                    for (int j2 = 0; j2 < N_Alpha; j2++) {
                        x2[j1][j2] = I1[j1][j2] - x3[j1][j2];
                    }
                }
                matinv3(x2, x3);
                matmul3(x3, x1, x2);
                for (int j1 = 0; j1 < N_Alpha; j1++) {
                    for (int j2 = 0; j2 < N_Alpha; j2++) {
                        cp1[j1][j2] += Cl[j][i] * x2[j1][j2];
                        cp2[j1][j2] += Cl[j][i] * x3[j1][j2];
                    }
                }
            }
            matinv3(cp2, x3);
            matmul3(x3, cp1, dcp[i]);
            for (int j = 0; j < N_Alpha; j++) {
                Doub f1 = abs(dcp[i][j][j]); //sqrt(dcp[i][j][j].r*dcp[i][j][j].r+dcp[i][j][j].i*dcp[i][j][j].i);
                //				cout << sqrt(dcp[i][j][j].r*dcp[i][j][j].r+dcp[i][j][j].i*dcp[i][j][j].i)/sqrt(Coh_p[i][j][j].r*Coh_p[i][j][j].r+Coh_p[i][j][j].i*Coh_p[i][j][j].i) << '\n';
                if (f1 > EPS_CP * abs(Coh_p[i][j][j])) //.r*Coh_p[i][j][j].r+Coh_p[i][j][j].i*Coh_p[i][j][j].i))
                {
                    l_CPA = 0;
                }
                if (f1 > f2) {
                    f2 = f1;
                }
            }
        }
        for (int i = 0; i < N_LAT; i++) {
            for (int j1 = 0; j1 < N_Alpha; j1++) {
                for (int j2 = 0; j2 < N_Alpha; j2++) {
                    Coh_p[i][j1][j2] += dcp[i][j1][j2];
                }
            }
        }
        for (int i = 0; i < N_LAT; i++) {
            for (int j1 = 0; j1 < N_Alpha; j1++) {
                for (int j2 = 0; j2 < N_Alpha; j2++) {
                    Coh_p[i][j1][j2] = (E < 0) ? complex<Doub > (real(Coh_p[i][j1][j2]), fabs(imag(Coh_p[i][j1][j2]))) :
                            complex<Doub > (real(Coh_p[i][j1][j2]), -fabs(imag(Coh_p[i][j1][j2])));
                }
            }
        }
        if (l_CPA == 1 && it > 0) {
            delete [] dcp;
            delete [] W;
            delete [] Wab;
            delete [] G_s;
            delete [] cp1;
            delete [] cp2;
            delete [] I1;
            delete [] Im_d;
            delete [] x1;
            delete [] x2;
            delete [] x3;
            return;
        }
    }
    if (f2 <= 1.0e-1) {
        l_CPA = 1;
    }
    delete [] dcp;
    delete [] W;
    delete [] Wab;
    delete [] G_s;
    delete [] cp1;
    delete [] cp2;
    delete [] I1;
    delete [] Im_d;
    delete [] x1;
    delete [] x2;
    delete [] x3;
}