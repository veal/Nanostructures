#include "sys_param.h"

void Density(Comp G_s[N_En][3][3][N_LAT][N_LAT][N_Band][N_Band], Doub g_coh[N_LAT], Doub g_cpa[N_LAT][N_Band],
        Doub g_comp[8], Doub gm[N_Com][N_LAT][2][N_Band], Doub g_m[N_Com][N_LAT][2][N_Band], int niE,
        int spin, double E, Comp W[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band], Doub vm[2][N_Com][N_LAT][2][N_Band],
        double Pm[2][N_LAT], double Cl[N_Com][N_LAT], double Eps, Comp Hr[Nkp2][N_LAT][N_LAT][N_Band][N_Band],
        Comp Coh_p_e[2][N_LAT][N_LAT][N_Band][N_Band], double r[N_LAT][3][3]) {

    Doub(*gEsm)[N_Com][N_LAT][N_Band] = new Doub[2][N_Com][N_LAT][N_Band];
    double Eps1 = 0, Eps_a1 = 0;
    Comp (*Coh_p)[N_Band][N_Band] = new Comp[N_LAT][N_Band][N_Band];
    Comp (*t)[N_Com][N_LAT][N_Band][N_Band] = new Comp[2][N_Com][N_LAT][N_Band][N_Band];
    Comp (*Wm_D)[N_Com][N_LAT][N_Band][N_Band] = new Comp[2][N_Com][N_LAT][N_Band][N_Band];
    Comp (*Gs)[3][N_LAT][N_LAT][N_Band][N_Band] = G_s[niE];
    Comp (*x1)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*x2)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*x3)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
    Comp (*x4)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*x5)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*x6)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*x7)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*x8)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*I1)[N_Band] = new Comp[N_Band][N_Band];
    Comp (*Sig_ef_e)[N_Band][N_Band] = new Comp[N_LAT][N_Band][N_Band];
    Comp (*sigmam1)[N_LAT][N_Band] = new Comp[N_Com][N_LAT][N_Band];
    Comp (*sigma_m1)[N_LAT][N_Band] = new Comp[N_Com][N_LAT][N_Band];
    Doub(*vm1)[N_Com][N_LAT][N_Band] = new Doub[2][N_Com][N_LAT][N_Band];
    Doub(*gm1)[N_LAT][N_Band] = new Doub[N_Com][N_LAT][N_Band];
    Doub(*g_m1)[N_LAT][N_Band] = new Doub[N_Com][N_LAT][N_Band];
    for (int nc = 0; nc < N_Com; nc++) {
        for (int lat = 0; lat < N_LAT; lat++) {
            for (int band = 0; band < N_Band; band++) {
                vm1[0][nc][lat][band] = vm[0][nc][lat][spin][band];
                vm1[1][nc][lat][band] = vm[1][nc][lat][spin][band];
            }
        }
    }
    
    G_site_all(Gs, E, spin, Hr, Coh_p_e, Sig_ef_e, r);
    
    for (int i = 0; i < N_LAT; i++) {
        g_coh[i] = 0.0;
        for (int ib = 0; ib < N_Band; ib++) {
            g_cpa[i][ib] = imag(Gs[1][1][i][i][ib][ib]);
            g_coh[i] -= g_cpa[i][ib] / Pi;
        }
    }

//    //************************
//
//    for (int i = 0; i < N_Band; i++) {
//        I1[i][i] = 1.0;
//        for (int m = 0; m < 2; m++) {
//            for (int j1 = 0; j1 < N_Com; j1++) {
//                for (int j2 = 0; j2 < N_LAT; j2++) {
//                    Wm_D[m][j1][j2][i][i] = vm1[m][j1][j2][i];
//                }
//            }
//        }
//    }
//
//    //***** set up single-site scattering matrix ***********
//    for (int i = 0; i < N_LAT; i++) {
//        for (int j = 0; j < N_Com; j++) {
//            for (int m = 0; m < 2; m++) {
//                for (int i1 = 0; i1 < N_LAT; i1++) {
//                    for (int j1 = 0; j1 < N_Band; j1++) {
//                        for (int j2 = 0; j2 < N_Band; j2++) {
//                            x1[i1*N_Band+j1][i1*N_Band+j2] = W[j][i][i1][i1][j1][j2] + Wm_D[m][j][i][j1][j2] - Coh_p_e[spin][i][i1][j1][j2];
//                        }
//                    }
//                }
//                matmul(x1, Gs[1][1], x2);
//                for (int i1 = 0; i1 < N_LAT; i1++) {
//                    for (int i2 = 0; i2 < N_LAT; i2++) {
//                        for (int j1 = 0; j1 < N_Band; j1++) {
//                            for (int j2 = 0; j2 < N_Band; j2++) {
//                                x3[i1*N_Band+j1][i1*N_Band+j2] = I1[i1][i2][j1][j2] - x2[i1*N_Band+j1][i1*N_Band+j2];
//                            }
//                        }
//                    }
//                }
//                matinv(x3, x2);
//                matmul(x2, x1, t[m][j][i]);
//            }
//        }
//    }
    //for (int j1 = 0; j1 < N_Band; j1++)
    //{
    //	for(int j2 = 0; j2 < N_Band; j2++)
    //	{
    //		Comp t_average;
    //		for (int i = 0; i < N_LAT; i++)
    //		{
    //			for (int j = 0; j < N_Com; j++)
    //			{
    //				for (int m = 0; m < 2; m++)
    //				{
    //					t_average += Cl[j][i]*Pm[m][i]*t[m][j][i][j1][j2];
    //				}
    //			}
    //		}
    //		if (sqrt(t_average.r*t_average.r + t_average.i*t_average.i) > EPS_CP)
    //		{
    //			cout << j1 << "   " << j2 << "   " << sqrt(t_average.r*t_average.r + t_average.i*t_average.i) << '\n';
    //		}
    //	}
    //}
    //********  Pairs  **********
    
//    for (int m = 0; m < 2; m++) {
//        for (int il = 0; il < N_LAT; il++) {
//            for (int i = 0; i < N_Com; i++) {
//                for (int ib = 0; ib < N_Band; ib++) {
//                        gEsm[m][i][il][ib] = 0.0;
//                }
//                matmul4(t[m][i][il][0][0], Gs[1][1][il][il], x5);
//                matmul4(Gs[1][1][il][il], x5, x4);
//                for (int n = 0; n < 3; n++) {
//                    for (int jl = 0; jl < N_LAT; jl++) {
//                        for (int m1 = 0; m1 < 2; m1++) {
//                            for (int j = 0; j < N_Com; j++) {
//                                if ((n != 1) || (il != jl)) {
//                                    matmul4(t[m][i][il][0][0], Gs[1][n][il][jl], x1);
//                                    matmul4(t[m1][j][jl][0][0], Gs[n][1][jl][il], x7);
//                                    matmul4(x1, x7, x6);
//                                    for (int b1 = 0; b1 < N_Band; b1++) {
//                                        for (int b2 = 0; b2 < N_Band; b2++) {
//                                            x2[b1][b2] = I1[b1][b2] - x6[b1][b2];
//                                        }
//                                    }
//                                    matinv4(x2, x8);
//                                    matmul4(Gs[1][1][il][il], x8, x2);
//                                    matmul4(x2, x6, x8);
//                                    for (int b1 = 0; b1 < N_Band; b1++) {
//                                        for (int b2 = 0; b2 < N_Band; b2++) {
//                                            x3[b1][b2] = I1[b1][b2] + x5[b1][b2];
//                                        }
//                                    }
//                                    matmul4(x8, x3, x2);
//                                    matmul4(x7, x1, x6);
//                                    for (int b1 = 0; b1 < N_Band; b1++) {
//                                        for (int b2 = 0; b2 < N_Band; b2++) {
//                                            x8[b1][b2] = I1[b1][b2] - x6[b1][b2];
//                                        }
//                                    }
//                                    matinv4(x8, x1);
//                                    matmul4(Gs[1][n][il][jl], x1, x8);
//                                    matmul4(x7, t[m][i][il][0][0], x1);
//                                    matmul4(x8, x1, x6);
//                                    matmul4(Gs[1][n][il][jl], x7, x1);
//                                    for (int b1 = 0; b1 < N_Band; b1++) {
//                                        for (int b2 = 0; b2 < N_Band; b2++) {
//                                            x3[b1][b2] = x1[b1][b2] + Gs[1][1][il][il][b1][b2];
//                                        }
//                                    }
//                                    matmul4(x6, x3, x1);
//                                    for (int ib = 0; ib < N_Band; ib++) {
//                                        double g_ib = (imag(x1[ib][ib]) + imag(x2[ib][ib])) * Pmm(m1, m, jl, il, Pm, Eps) *
//                                                Cll(j, i, jl, il, Cl, Eps_a1);
//                                        gEsm[m][i][il][ib] -= g_ib / Pi;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                for (int ib = 0; ib < N_Band; ib++) {
//                    double g_ib = imag(x4[ib][ib]) + g_cpa[il][ib];
//                    gEsm[m][i][il][ib] -= g_ib / Pi;
//                }
//            }
//        }
//        if (m == 1) {
//            for (int jn = 0; jn < N_Com; jn++) {
//                for (int ji = 0; ji < N_LAT; ji++) {
//                    for (int jb = 0; jb < N_Band; jb++) {
//                        gm1[jn][ji][jb] = gEsm[m][jn][ji][jb];
//                    }
//                }
//            }
//        } else {
//            for (int jn = 0; jn < N_Com; jn++) {
//                for (int ji = 0; ji < N_LAT; ji++) {
//                    for (int jb = 0; jb < N_Band; jb++) {
//                        g_m1[jn][ji][jb] = gEsm[m][jn][ji][jb];
//                    }
//                }
//            }
//        }
//    }
//    for (int jl = 0; jl < N_LAT; jl++) {
//        for (int jb = 0; jb < N_Band; jb++) {
//            g_cpa[jl][jb] /= -Pi;
//        }
//    }
//    for (int il = 0; il < N_LAT; il++) {
//        for (int ib = 0; ib < N_Band; ib++) {
//            for (int m = 0; m < 2; m++) {
//                for (int i = 0; i < N_Com; i++) {
//                    Doub dd = gEsm[m][i][il][ib] * Cl[i][il] * Pm[m][il] / N_LAT;
//                    if (m == 1) {
//                        g_comp[0] += dd;
//                    }
//                    if (m == 0) {
//                        g_comp[1] += dd;
//                    }
//                    if (i == 0) {
//                        g_comp[2] += dd;
//                    }
//                    if (i == 1) {
//                        g_comp[3] += dd;
//                    }
//                    if (ib == 1) {
//                        g_comp[4] += dd;
//                    }
//                    if (ib > 1 && ib <= 3) {
//                        g_comp[5] += dd;
//                    }
//                    if (ib > 3 && ib <= 8) {
//                        g_comp[6] += dd;
//                    }
//                    if (ib == 9) {
//                        g_comp[7] += dd;
//                    }
//                }
//            }
//        }
//    }
//    for (int nc = 0; nc < N_Com; nc++) {
//        for (int lat = 0; lat < N_LAT; lat++) {
//            for (int band = 0; band < N_Band; band++) {
//                gm[nc][lat][spin][band] = gm1[nc][lat][band];
//                g_m[nc][lat][spin][band] = g_m1[nc][lat][band];
//            }
//        }
//    }
    delete [] gEsm;
    delete [] Coh_p;
    delete [] t;
    delete [] Wm_D;
    delete [] x1;
    delete [] x2;
    delete [] x3;
    delete [] x4;
    delete [] x5;
    delete [] x6;
    delete [] x7;
    delete [] x8;
    delete [] I1;
    delete [] Sig_ef_e;
    delete [] sigmam1;
    delete [] sigma_m1;
    delete [] vm1;
    delete [] gm1;
    delete [] g_m1;
}