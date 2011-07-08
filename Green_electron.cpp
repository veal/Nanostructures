#include "sys_param.h"

void G_site(Comp G_s[N_LAT][N_LAT][N_Band][N_Band], Doub E, Comp Coh_p[N_LAT][N_LAT][N_Band][N_Band],
        Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]) {
    pthread_t aThread[n_Cores];
    Comp (*Gs)[N_LAT][N_LAT][N_Band][N_Band] = new Comp[n_Cores][N_LAT][N_LAT][N_Band][N_Band];
    parameters *param = new parameters[n_Cores];
    for (int iter = 0; iter < n_Cores; iter++) {
        param[iter].num = iter;
        param[iter].CP = Coh_p;
        param[iter].Energy = E;
        param[iter].Hk = Hr;
        param[iter].Gs = Gs;
    }

    for (int iter = 0; iter < n_Cores; iter++) {
        pthread_create(&aThread[iter], NULL, &CalculateGreen, &param[iter]);
    }

    for (int iter = 0; iter < n_Cores; iter++) {
        pthread_join(aThread[iter], NULL);
    }

    for (int j = 0; j < N_LAT; j++) {
        for (int b1 = 0; b1 < N_Band; b1++) {
            for (int b2 = 0; b2 < N_Band; b2++) {
                G_s[j][j][b1][b2] = 0;
                for (int n_core = 0; n_core < n_Cores; n_core++) {
                    G_s[j][j][b1][b2] += Gs[n_core][j][j][b1][b2];
                }
            }
        }
    }
    delete [] Gs;
    delete [] param;
}

void G_site_all(Comp Gs[3][3][N_LAT][N_LAT][N_Band][N_Band], Doub E, int spin, Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band],
        Comp Coh_p_e[2][N_LAT][N_LAT][N_Band][N_Band], Comp Sig_ef_e[N_LAT][N_Band][N_Band],
        Doub r[N_LAT][3][3]) {
    Comp (*G_s)[3][3][N_LAT][N_LAT][N_Band][N_Band] = new Comp[n_Cores][3][3][N_LAT][N_LAT][N_Band][N_Band];
    //	Comp (*G2s)[3][N_LAT][N_LAT][N_Band][N_Band] = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    //	Comp (*Gos)[3][N_LAT][N_LAT][N_Band][N_Band] = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    //	Comp (*Go_s)[3][N_LAT][N_LAT][N_Band][N_Band] = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    Comp (*H_k)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*H_k_tem)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*G_k_tem)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];

    Comp (*C_P)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*G_k)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    //	Comp (*G2_k)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*I1)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*Sig_ef_e_l)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*Coh_p)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_Band; j++) {
            I1[i][i][j][j] = 1.0;
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int i1 = 0; i1 < N_LAT; i1++) {
            for (int j1 = 0; j1 < N_Band; j1++) {
                for (int j2 = 0; j2 < N_Band; j2++) {
                    Coh_p[i][i1][j1][j2] = Coh_p_e[spin][i][i1][j1][j2];
                }
            }
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int i1 = 0; i1 < N_LAT; i1++) {
            for (int j1 = 0; j1 < N_Band; j1++) {
                for (int j2 = 0; j2 < N_Band; j2++) {
                    C_P[i1][i1][j1][j2] += Coh_p[i][i1][j1][j2];
                    Sig_ef_e_l[i][i][j1][j2] = Sig_ef_e[i][j1][j2];
                }
            }
        }
    }

    pthread_t aThread[n_Cores];
    parameters2 *param = new parameters2[n_Cores];
    for (int iter = 0; iter < n_Cores; iter++) {
        param[iter].num = iter;
        param[iter].CP = Coh_p;
        param[iter].Energy = E;
        param[iter].Hk = Hr;
        param[iter].Gs = G_s;
    }

    for (int iter = 0; iter < n_Cores; iter++) {
        pthread_create(&aThread[iter], NULL, &CalculateGreenFull, &param[iter]);
    }

    for (int iter = 0; iter < n_Cores; iter++) {
        pthread_join(aThread[iter], NULL);
    }

    for (int n = 0; n < 3; n++) {
        for (int ns = 0; ns < 3; ns++) {
            for (int j1 = 0; j1 < N_LAT; j1++) {
                for (int j2 = 0; j2 < N_LAT; j2++) {
                    for (int b1 = 0; b1 < N_Band; b1++) {
                        for (int b2 = 0; b2 < N_Band; b2++) {
                            Gs[n][ns][j1][j2][b1][b2] = 0;
                            for (int n_core = 0; n_core < n_Cores; n_core++) {
//                                if (n == ns && j1 == j2 && b1 == b2 && imag(G_s[n_core][n][ns][j1][j2][b1][b2]) > 0.0)
//                                    cout << imag(G_s[n_core][n][ns][j1][j2][b1][b2]) << '\t' << n_core << '\n';
                                Gs[n][ns][j1][j2][b1][b2] += G_s[n_core][n][ns][j1][j2][b1][b2];
                            }
                        }
                    }
                }
            }
        }
    }
    delete [] G_s;
    delete [] H_k;
    delete [] H_k_tem;
    delete [] G_k_tem;

    delete [] C_P;
    delete [] G_k;
    //	delete [] G2_k;
    delete [] I1;
    delete [] Sig_ef_e_l;
    delete [] Coh_p;
}

void* CalculateGreen(void *param) {
    parameters* t = (parameters*) (param);
    int num = t->num;
    Comp (*G_s)[N_LAT][N_Band][N_Band];
    Comp (*Hr)[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*C_P)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*C_P1)[N_Band][N_Band] = new Comp[N_LAT][N_Band][N_Band];
    Comp (*I)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*H_k_tem)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*G_k_tem)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*Coh_p)[N_LAT][N_Band][N_Band];
    Coh_p = t->CP;
    G_s = t->Gs[num];
    Hr = t->Hk;
    Doub E = t->Energy;
    int NKP1 = -Nkp + 2 * num * Nkp / n_Cores;
    int NKP2 = -Nkp + 2 * (num + 1) * Nkp / n_Cores;
    for (int j = 0; j < N_LAT; j++) {
        for (int j1 = 0; j1 < N_Band; j1++) {
            I[j][j][j1][j1] = 1.0;
        }
    }
    pthread_mutex_lock(&job_queue);
    //	WaitForSingleObject( ghMutex, INFINITE );
    
//    for (int i1 = 0; i1 < N_LAT; i1++) {
//        for (int j1 = 0; j1 < N_Band; j1++) {
//            for (int j2 = 0; j2 < N_Band; j2++) {
//                C_P1[i1][j1][j2] = 0.0;
//                for (int j = 0; j < N_LAT; j++) {
//                    C_P1[i1][j1][j2] += Coh_p[j][i1][j1][j2];
//                }
//            }
//        }
//    }
    
    for (int j = 0; j < N_LAT; j++) {
        for (int j1 = 0; j1 < N_Band; j1++) {
            for (int j2 = 0; j2 < N_Band; j2++) {
                C_P[j][j][j1][j2] = Coh_p[j][j][j1][j2];
            }
        }
    }
    
    Comp gg(0.0, -0.001); //!!!!!!!!!!!!!!!!!!
    pthread_mutex_unlock(&job_queue);
    //	ReleaseMutex(ghMutex);
    for (int k1 = NKP1; k1 < NKP2; k1++) {
        Doub ro_loc = 1;
        Doub kx = 2 * Pi * static_cast<Doub> (k1) / static_cast<Doub> (3 * Nkp);
        for (int i = 0; i < N_LAT; i++) {
            for (int i1 = 0; i1 < N_LAT; i1++) {
                for (int j = 0; j < N_Band; j++) {
                    for (int j1 = 0; j1 < N_Band; j1++) {
//                        Comp p = CI * kx * (r[i1][0][1] - r[i][0][1]);
                        H_k_tem[i * N_Band + j][i1 * N_Band + j1] = E * I[i][i1][j][j1] -
                                (C_P[i][i1][j][j1]/**exp(p)*/ + Hr[k1 + Nkp][i][i1][j][j1]);
//						H_k_tem[i*N_Band+j][i1*N_Band+j1] = E*I[i][i1][j][j1] -
//							(gg*I[i][i1][j][j1] + Hr[k1+Nkp][i][i1][j][j1]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    }
                }
            }
        }
        cd_Invert(H_k_tem, G_k_tem);
        if ((k1 == Nkp - 1) || (k1 == -Nkp)) {
            ro_loc = 0.5;
        }
        for (int j = 0; j < N_LAT; j++) {
            for (int b1 = 0; b1 < N_Band; b1++) {
                for (int b2 = 0; b2 < N_Band; b2++) {
                    G_s[j][j][b1][b2] += G_k_tem[j * N_Band + b1][j * N_Band + b2] * ro_loc / static_cast<Doub> (Nkp2);
                }
            }
        }
    }
    delete [] I;
    delete [] C_P;
    delete [] C_P1;
    delete [] G_k_tem;
    delete [] H_k_tem;
    return NULL;
}

void* CalculateGreenFull(void *param) {
    parameters2* t = (parameters2*) (param);
    int num = t->num;
    Comp (*G_s)[3][N_LAT][N_LAT][N_Band][N_Band];
    Comp (*Hr)[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*C_P)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*C_P1)[N_Band][N_Band] = new Comp[N_LAT][N_Band][N_Band];
    Comp (*I)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
    Comp (*H_k_tem)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*G_k_tem)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*Coh_p)[N_LAT][N_Band][N_Band];
    Coh_p = t->CP;
    G_s = t->Gs[num];
    Hr = t->Hk;
    Doub E = t->Energy;
    int NKP1 = -Nkp + 2 * num * Nkp / n_Cores;
    int NKP2 = -Nkp + 2 * (num + 1) * Nkp / n_Cores;
    for (int j = 0; j < N_LAT; j++) {
        for (int j1 = 0; j1 < N_Band; j1++) {
            I[j][j][j1][j1] = 1.0;
        }
    }
    pthread_mutex_lock(&job_queue);

    for (int i1 = 0; i1 < N_LAT; i1++) {
        for (int j1 = 0; j1 < N_Band; j1++) {
            for (int j2 = 0; j2 < N_Band; j2++) {
                C_P1[i1][j1][j2] = 0.0;
                for (int j = 0; j < N_LAT; j++) {
                    C_P1[i1][j1][j2] += Coh_p[j][i1][j1][j2];
                }
            }
        }
    }

    for (int j = 0; j < N_LAT; j++) {
        for (int j1 = 0; j1 < N_Band; j1++) {
            for (int j2 = 0; j2 < N_Band; j2++) {
                C_P[j][j][j1][j2] = C_P1[j][j1][j2];
//                if (imag(Coh_p[j][j1][j2]) > 0)
//                    cout << j << "  " << j1 << "  " << j2 << "  " << real(Coh_p[j][j1][j2])
//                            << "  " << imag(Coh_p[j][j1][j2]) << '\n';
            }
        }
    }
    
    Comp gg(0.0, -0.01); //!!!!!!!!!!!!!!!!!!
    pthread_mutex_unlock(&job_queue);

    for (int k1 = NKP1; k1 < NKP2; k1++) {
        Doub kx = 2 * Pi * static_cast<Doub> (k1) / static_cast<Doub> (3 * Nkp);
        Doub ro_loc = 1;
//        double* diag_val = new double[N_LAT*N_Band];
//        Eigen_values(Hr[k1+Nkp], diag_val, true);
        for (int i = 0; i < N_LAT; i++) {
            for (int i1 = 0; i1 < N_LAT; i1++) {
                for (int j = 0; j < N_Band; j++) {
                    for (int j1 = 0; j1 < N_Band; j1++) {
//                        H_k_tem[i * N_Band + j][i1 * N_Band + j1] = E * I[i][i1][j][j1] -
//                                (C_P[i][i1][j][j1] + Hr[k1 + Nkp][i][i1][j][j1]);
                        H_k_tem[i * N_Band + j][i1 * N_Band + j1] = E * I[i][i1][j][j1] -
                                (gg*I[i][i1][j][j1] + Hr[k1 + Nkp][i][i1][j][j1]);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    }
                }
            }
        }
        cd_Invert(H_k_tem, G_k_tem);

        if ((k1 == Nkp - 1) || (k1 == -Nkp)) {
            ro_loc = 0.5;
        }
        for (int n = 0; n < 3; n++) {
            for (int ns = 0; ns < 3; ns++) {
                for (int i1 = 0; i1 < N_LAT; i1++) {
                    for (int i2 = 0; i2 < N_LAT; i2++) {
                        for (int j1 = 0; j1 < N_Band; j1++) {
                            for (int j2 = 0; j2 < N_Band; j2++) {
                                Comp p = CI * kx * (r[i1][0][n] - r[i2][0][ns]);
                                G_s[n][ns][i1][i2][j1][j2] += G_k_tem[i1 * N_Band + j1][i2 * N_Band + j2] * ro_loc * exp(p) / static_cast<Doub> (2 * Nkp);
                            }
                        }
                    }
                }
            }
        }
    }
    delete [] I;
    delete [] C_P;
    delete [] C_P1;
    delete [] G_k_tem;
    delete [] H_k_tem;
    return NULL;
}