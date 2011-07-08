#include "lab.h"

Lab::Lab(){
}
void Lab::run() {
    runCalculation();
}

void Lab::runCalculation() {
    //�������� ������ �� ���� ��� �������������� ��������
    allocate();
    //������������� ���������� ����� � ������
    Set_coord(r);

    Scat_Pr(ni0, mm0, Pm, Cl);
    //��� ������� ���� �� �������������
    Read_input_file(ni0, mm0, nm_1, mm_1);
    //�������� ����������� ��������� (39,39)
    Read_Coulomb_integral(U_matr);
    Hop_integrals* hopIntegral = new Hop_integrals(matrixInputFile.toStdString(), impurityInputFile.toStdString());
    Calculate_Hamiltonian(hopIntegral);
    return;
    Calculate_Wamiltonian(hopIntegral->getMatrixPWF(), hopIntegral->getImpurityPWF());

    //    return 0;
    ofstream f_tot_start[2], f_tot;
    if (calculateGreen == true) {
        Doub ni_E[100] = {0};
        Doub ni_E_1[100] = {0};
        if (calculate_nm_mm == true) {
            // ���� �� ������� ��������� (37,2)
            for (int nEps = 0; nEps < 1; nEps++) {
                ni_E[nEps] = 0.0;
                ni_E_1[nEps] = 0.0;
                Doub Eps = 0.20 + 0.005 * nEps;
                Doub free_energy = 0;
                Doub(*nm)[N_Com][N_LAT][N_Band] = new Doub[2][N_Com][N_LAT][N_Band];
                Doub(*mm)[N_Com][N_LAT][N_Band] = new Doub[2][N_Com][N_LAT][N_Band];
                // �������� ��������� ���������� ��������� ������� �� ����� ������� (36,2)
                for (int imm = itmS; imm < it_m; imm++) {
                    int Stop_intE = 0;
                    Doub Int_E = 0;
                    Doub Int_E_1 = 0;
                    if (/*(mm_stop == 0) || */(l_data == 1)) {
                        Doub Int_E = 0.0, Ee2 = 0.0, Ee1 = 0.0;
                        for (int i1 = 0; i1 < N_LAT; i1++) {
                            for (int j = 0; j < N_Com; j++) {
                                for (int m = 0; m < 2; m++) {
                                    for (int ib = 0; ib < N_Band; ib++) {
                                        nm[m][j][i1][ib] = 0;
                                        mm[m][j][i1][ib] = 0;
                                    }
                                }
                            }
                        }
                        Doub Ee10 = 0.0;
                        Doub ni2[N_LAT] = {0.0};
                        for (int ib = 0; ib < N_Band; ib++) {
                            for (int i1 = 0; i1 < N_LAT; i1++) {
                                for (int j = 0; j < N_Com; j++) {
                                    for (int m = 0; m < 2; m++) {
                                        for (int spin = 0; spin < 2; spin++) {
                                            if (spin == 1) // spin_up(+)
                                            {
                                                vm[m][j][i1][spin][ib] = U_matr[ib]*(nm_1[m][j][i1][ib] - mm_1[m][j][i1][ib]) / 2.0;
                                            } else // spin_down(-)
                                            {
                                                vm[m][j][i1][spin][ib] = U_matr[ib]*(nm_1[m][j][i1][ib] + mm_1[m][j][i1][ib]) / 2.0;
                                            }
                                            //                                              vm[m][j][i1][spin][ib] -= (8.525e-6)*100*0.5*2*(spin-0.5);
                                        }
                                        ni2[i1] += Cl[j][i1] * Pm[m][i1]*(mm_1[m][j][i1][ib] * mm_1[m][j][i1][ib] - nm_1[m][j][i1][ib] * nm_1[m][j][i1][ib]);
                                    }
                                    Ee10 += U_matr[ib] * ni2[i1] / (4.0 * N_LAT);
                                }
                            }
                        }

                        Doub (*g_pkp)[2] = new Doub[N_En][2];
                        Doub (*g)[2] = new Doub[N_En][2];
                        Doub (*gsm)[N_LAT][2][N_Band] = new Doub[N_Com][N_LAT][2][N_Band];
                        Doub (*gs_m)[N_LAT][2][N_Band] = new Doub[N_Com][N_LAT][2][N_Band];
                        Doub (*gsm_1)[N_LAT][2][N_Band] = new Doub[N_Com][N_LAT][2][N_Band];
                        Doub (*gs_m_1)[N_LAT][2][N_Band] = new Doub[N_Com][N_LAT][2][N_Band];
                        Doub (*gsm_2)[N_LAT][2][N_Band] = new Doub[N_Com][N_LAT][2][N_Band];
                        Doub (*gs_m_2)[N_LAT][2][N_Band] = new Doub[N_Com][N_LAT][2][N_Band];
                        //���� �� �������
                        for (int niE = 0; niE < N_En; niE++) {
                            Doub *g_coh = new Doub[N_LAT];
                            Doub (*g_cpa)[N_LAT][N_Band] = new Doub[2][N_LAT][N_Band];
                            Doub (*g_comp)[8] = new Doub[2][8];
                            Doub dE = (En_E - En_S) / N_En;
                            Doub E = En_S + dE*niE;
                            for (int spin = 0; spin < 2; spin++) {
                                ee[niE] = E;
                                // ������� ����������� ��������� (32,51)
                                double CP_ratio = 0;
                                //                                CP_Calculation(Coh_p[niE][spin], E, niE, W, vm, Pm, Cl, H_k, r, spin, &CP_ratio);
                                //                                cout << "CP differs from average t in " << CP_ratio << " times" << '\n';
                                //                                return 0;
                                // ��������� ��������� ����������� ��������� (36,4) �� � �����
                                // ���������� ����� � ������������ ����������� ����������� (42 ���. Gaa)
                                for (int num = 0; num < 8; num++)
                                    g_comp[spin][num] = 0.0;

                                Density(Gs, g_coh, g_cpa[spin], g_comp[spin], gsm, gs_m, niE, spin, E, W, vm,
                                        Pm, Cl, Eps, H_k, Coh_p[niE], r);
                                g_pkp[niE][spin] = 0.0;
                                g[niE][spin] = 0.0;
                                for (int i1 = 0; i1 < N_LAT; i1++) {
                                    for (int ib = 0; ib < N_Band; ib++) {
                                        g_cpa[spin][i1][ib] /= N_LAT;
                                    }
                                    if (g_coh[i1] < 0.0)
                                        cout << g_coh[i1] << "  " << '\n';
                                    g_coh[i1] /= N_LAT;
                                    g_pkp[niE][spin] += g_coh[i1];
                                    for (int j = 0; j < N_Com; j++) {
                                        for (int m = 0; m < 2; m++) {
                                            for (int ib = 0; ib < N_Band; ib++) {
                                                if (m == -1) {
                                                    gmi = gs_m[j][i1][spin][ib];
                                                } else {
                                                    gmi = gsm[j][i1][spin][ib];
                                                }
                                                g[niE][spin] += Cl[j][i1] * Pm[m][i1] * gmi / N_LAT;
                                            }
                                        }
                                    }
                                }
                                //                                if (g_pkp[niE][spin] < 1.0e-25) {
                                //                                    g_pkp[niE][spin] = 0.0;
                                //                                }
                                if (g[niE][spin] < 1.0e-25) {
                                    g[niE][spin] = 0.0;
                                }
                            }
                            f_tot.open("Rezults/density.dat", ios::app);
                            f_tot << E - Esmr << '\t' << g_pkp[niE][0] + g_pkp[niE][1] << '\t' << g[niE][0] + g[niE][1] <<
                                    '\t' << l_CPA << '\n';
                            if (niE == N_En - 1) {
                                f_tot << "END OF ITERATION...   " << imm << "   " << nEps << '\n';
                            }
                            f_tot.close();
                            f_tot_start[0].open("Rezults/g-.dat", ios::app);
                            f_tot_start[1].open("Rezults/g+.dat", ios::app);
                            for (int spin = 0; spin < 2; spin++) {
                                f_tot_start[spin] << E - Esmr << '\t' <<
                                        g_comp[spin][0] + g_comp[spin][1] << '\t' <<
                                        g_comp[spin][0] << '\t' <<
                                        g_comp[spin][1] << '\t' <<
                                        g_comp[spin][2] << '\t' <<
                                        g_comp[spin][3] << '\t' <<
                                        g_comp[spin][4] << '\t' <<
                                        g_comp[spin][5] << '\t' <<
                                        g_comp[spin][6] << '\t' <<
                                        g_comp[spin][7] << '\t' << '\n';
                                if (niE == N_En - 1) {
                                    f_tot_start[spin] << "END OF ITERATION...   " << imm << '\t' << nEps << '\n';
                                }
                            }
                            f_tot_start[0].close();
                            f_tot_start[1].close();
                            cout << "E = " << E << '\t' << "l_CPA = " << l_CPA << '\n';
                            cout << g_pkp[niE][0] + g_pkp[niE][1] << '\t' << g[niE][0] + g[niE][1] << '\n';
                            //������ ���-�� ������� �� (33,17) ��� ���������� ������ ����� � ��������� �������� ��. (36,2)
                            if (Stop_intE == 0) {
                                Int_E_1 = Int_E;
                                //if (niE >= 2)
                                //{
                                //	Int_Simpson(Int_E, niE, g[niE-2][0].r + g[niE-2][1].r, g[niE-1][0].r + g[niE-1][1].r, g[niE][0].r + g[niE][1].r, dE);
                                //}
                                if (niE == 0) {
                                    Int_E += (g[niE][0] + g[niE][1]) * dE / 2;
                                    if (imm == it_m - 1) {
                                        free_energy += (ee[niE] * g[niE][0] + ee[niE] * g[niE][1]) * dE / 2;
                                    }
                                    //Int_Simpson(Int_E, niE, 0.0, 0.0, g[niE][0].r + g[niE][1].r, dE);
                                }
                                if (niE >= 1) {
                                    if (imm == it_m - 1) {
                                        free_energy += ee[niE] * (g[niE - 1][0] + g[niE - 1][1] + g[niE][0] + g[niE][1]) * dE / 2;
                                    }
                                    Int_E += (g[niE - 1][0] + g[niE - 1][1] + g[niE][0] + g[niE][1]) * dE / 2;
                                    //Int_Simpson(Int_E,niE,0.0,g[niE-1][0].r + g[niE-1][1].r, g[niE][0].r + g[niE][1].r, dE);
                                }
                                cout << "Int_E = " << Int_E << '\n';
                                if (Int_E >= N) {
                                    Stop_intE = 1;
                                    if (fabs(N - Int_E_1) > fabs(Int_E - N)) {
                                        Ef = E;
                                        N_f_int = niE;
                                    } else {
                                        Ef = E - dE;
                                        N_f_int = niE - 1;
                                    }
                                }
                                Mm_nm(gsm, gs_m, gs_m_1, gsm_1, gsm_2, gs_m_2, Pm, nm, mm, ee, niE, dE, Ef, T);
                                for (int i1 = 0; i1 < N_LAT; i1++) {
                                    for (int j = 0; j < N_Com; j++) {
                                        for (int spin = 0; spin < 2; spin++) {
                                            for (int ib = 0; ib < N_Band; ib++) {
                                                gsm_2[j][i1][spin][ib] = gsm_1[j][i1][spin][ib];
                                                gs_m_2[j][i1][spin][ib] = gs_m_1[j][i1][spin][ib];
                                                gsm_1[j][i1][spin][ib] = gsm[j][i1][spin][ib];
                                                gs_m_1[j][i1][spin][ib] = gs_m[j][i1][spin][ib];
                                                gsm[j][i1][spin][ib] = 0.0;
                                                gs_m[j][i1][spin][ib] = 0.0;
                                            }
                                        }
                                    }
                                }
                            }
                            delete [] g_coh;
                            delete [] g_cpa;
                            delete [] g_comp;
                        }
                        for (int i1 = 0; i1 < N_LAT; i1++) {
                            for (int j = 0; j < N_Com; j++) {
                                for (int m = 0; m < 2; m++) {
                                    for (int ib = 0; ib < N_Band; ib++) {
                                        if (abs(mm[m][j][i1][ib] - mm_1[m][j][i1][ib]) > abs(mm[m][j][i1][ib]) * del_m) {
                                            mm_stop = 0;
                                        } else {
                                            mm_stop = 1;
                                        }
                                        //cout << "m=" << m << "   " << "j=" << j << "   " << "i1=" << i1 << "   " <<
                                        //	"ib=" << ib << "   " << mm[m][j][i1][ib] << '\n';
                                        //cout << "m=" << m << "   " << "j=" << j << "   " << "i1=" << i1 << "   " <<
                                        //	"ib=" << ib << "   " << nm[m][j][i1][ib] << '\n';
                                        mm_1[m][j][i1][ib] = mm[m][j][i1][ib];
                                        nm_1[m][j][i1][ib] = nm[m][j][i1][ib];
                                        nm[m][j][i1][ib] = 0;
                                        mm[m][j][i1][ib] = 0;
                                    }
                                }
                            }
                        }
                        delete [] gsm;
                        delete [] gs_m;
                        delete [] g;
                        delete [] g_pkp;
                        delete [] gsm_1;
                        delete [] gs_m_1;
                        delete [] gsm_2;
                        delete [] gs_m_2;
                        ofstream file_out_mm("Rezults/mm.dat", ios::app);
                        ofstream file_out_nm("Rezults/nm.dat", ios::app);
                        Doub nm_tot1 = 0;
                        for (int i1 = 0; i1 < N_LAT; i1++) {
                            for (int j = 0; j < N_Com; j++) {
                                for (int m = 0; m < 2; m++) {
                                    Doub nm_tot = 0;
                                    Doub mm_tot = 0;
                                    for (int ib = 0; ib < N_Band; ib++) {
                                        nm_tot += nm_1[m][j][i1][ib];
                                        mm_tot += mm_1[m][j][i1][ib];
                                        nm_tot1 += Cl[j][i1] * Pm[m][i1] * nm_1[m][j][i1][ib];
                                        file_out_nm << m << "   " << j << "   " << i1 << "   " << ib << "   " << nm_1[m][j][i1][ib] << '\n';
                                        file_out_mm << m << "   " << j << "   " << i1 << "   " << ib << "   " << mm_1[m][j][i1][ib] << '\n';
                                    }
                                    file_out_nm << m << "   " << j << "   " << i1 << "         " << nm_tot << "    TOTAL!!!" << '\n';
                                    file_out_mm << m << "   " << j << "   " << i1 << "         " << mm_tot << "    TOTAL!!!" << '\n';
                                }
                            }
                        }
                        file_out_nm << "END OF ITERATION...   " << imm << "   " << nEps << '\n';
                        file_out_nm << "Total number of electrons   " << nm_tot1 << '\n';
                        file_out_mm.close();
                        file_out_nm.close();
                    }
                }
                delete [] nm;
                delete [] mm;
                ofstream file_out_En("Rezults/free_energy.dat", ios::app);
                file_out_En << Eps << "    " << free_energy << '\n';
                file_out_En.close();
                ofstream file_out;
                file_out.open("Rezults/Green.dat", ios::binary);
                file_out.write((char*) Gs, N_En * 9 * N_LAT * N_LAT * N_Band * N_Band * sizeof (complex<Doub>));
                file_out.close();
                file_out.open("Rezults/Coh_pot+.dat"/*, ios::binary*/);
                for (int eee = 0; eee < N_En; eee++) {
                    for (int nl = 0; nl < N_LAT; nl++) {
                        for (int nlsh = 0; nlsh < N_LAT; nlsh++) {
                            for (int nb1 = 0; nb1 < N_Band; nb1++) {
                                for (int nb2 = 0; nb2 < N_Band; nb2++) {
                                    file_out << real(Coh_p[eee][1][nl][nlsh][nb1][nb2]) << "   " << imag(Coh_p[eee][1][nl][nlsh][nb1][nb2]) << '\n';
                                }
                            }
                        }
                    }
                }
                //				file_out.write((char*)Coh_p, N_En*2*N_LAT*N_Band*N_Band*sizeof(complex<Doub>));
                file_out.close();
                file_out.open("Rezults/Coh_pot-.dat"/*ios::binary*/);
                for (int eee = 0; eee < N_En; eee++) {
                    for (int nl = 0; nl < N_LAT; nl++) {
                        for (int nlsh = 0; nlsh < N_LAT; nlsh++) {
                            for (int nb1 = 0; nb1 < N_Band; nb1++) {
                                for (int nb2 = 0; nb2 < N_Band; nb2++) {
                                    file_out << real(Coh_p[eee][0][nl][nlsh][nb1][nb2]) << "   " << imag(Coh_p[eee][0][nl][nlsh][nb1][nb2]) << '\n';
                                }
                            }
                        }
                    }
                }
                //				file_out.write((char*)Coh_p, N_En*2*N_LAT*N_Band*N_Band*sizeof(complex<Doub>));
                file_out.close();
                file_out.open("Rezults/nm.bin", ios::binary);
                file_out.write((char*) nm_1, 2 * N_Com * N_LAT * N_Band * sizeof (complex<Doub>));
                file_out.close();
                file_out.open("Rezults/mm.bin", ios::binary);
                file_out.write((char*) mm_1, 2 * N_Com * N_LAT * N_Band * sizeof (complex<Doub>));
                file_out.close();
            }
            return;
        } else {
            ifstream file_out;
            file_out.open("Rezults/Green.dat", ios::binary);
            file_out.read((char*) Gs, N_En * 9 * N_LAT * N_LAT * N_Band * N_Band * sizeof (complex<Doub>));
            file_out.close();
            file_out.open("Rezults/Coh_pot.dat", ios::binary);
            file_out.read((char*) Coh_p, N_En * 2 * N_LAT * N_Band * N_Band * sizeof (complex<Doub>));
            file_out.close();
            file_out.open("Rezults/nm.bin", ios::binary);
            file_out.read((char*) nm_1, 2 * N_Com * N_LAT * N_Band * sizeof (complex<Doub>));
            file_out.close();
            file_out.open("Rezults/mm.bin", ios::binary);
            file_out.read((char*) mm_1, 2 * N_Com * N_LAT * N_Band * sizeof (complex<Doub>));
            file_out.close();
        }
    } else {
        cout << "Reading Green function.." << '\n';
        for (int niE = 0; niE < N_En; niE++) {
            Doub dE = (En_E - En_S) / N_En;
            Doub E = En_S + dE*niE;
            ee[niE] = E;
        }
        Ef = -0.28;
        ifstream file_out;
        file_out.open("Rezults/Green.dat", ios::binary);
        file_out.read((char*) Gs, N_En * 9 * N_LAT * N_LAT * N_Band * N_Band * sizeof (complex<Doub>));
        file_out.close();
    }
    if (calculateSig_fe == true) {
        Sigma_fe(Ef, T, Gs, Sig_fe, ee, Cl);
        ofstream file_out;
        file_out.open("Rezults/Sig_fe.dat", ios::binary);
        file_out.write((char*) Sig_fe, 9 * N_LAT * N_LAT * N_Alpha * N_Alpha * sizeof (complex<Doub>));
        file_out.close();
    } else {
        cout << "Reading Sigma_fe.." << '\n';
        ifstream file_out;
        file_out.open("Rezults/Sig_fe.dat", ios::binary);
        file_out.read((char*) Sig_fe, 9 * N_LAT * N_LAT * N_Alpha * N_Alpha * sizeof (complex<Doub>));
        file_out.close();
    }
    cout << "Calculating D_k.." << '\n';
    D_matrix(D_k, Sig_fe, r);
    ofstream file_out;
    file_out.open("Rezults/D_k.dat");
    for (int k11 = -Nkp; k11 < Nkp; k11++) {
        for (int i = 0; i < N_LAT; i++) {
            for (int j = 0; j < N_LAT; j++) {
                for (int i1 = 0; i1 < N_Alpha; i1++) {
                    for (int j1 = 0; j1 < N_Alpha; j1++) {
                        file_out << real(D_k[k11 + Nkp][i][j][i1][j1]) << "   " << imag(D_k[k11 + Nkp][i][j][i1][j1]) << '\n';
                        //if ((D_k[k11+Nkp][i][j][i1][j1].r != D_k[k11+Nkp][j][i][j1][i1].r) || (D_k[k11+Nkp][i][j][i1][j1].i != -D_k[k11+Nkp][j][i][j1][i1].i))
                        //{
                        //	cout << i << "   " << j << "   " << i1 << "   " << j1 << '\n';
                        //}
                    }
                }
            }
        }
    }

    //file_out.close();
    //*********************Interval for Conductivity*********************************
    //*******************************************************************************

    //		l_conduct = 1;
    //		if (T == 0.0)
    //		{
    //			l_T = 1;
    //			N_start = N_f_int;
    //			N_end = N_f_int;
    //		}
    //		else
    //		{
    //			E_N_Start = Ef - 3.01008e-4*T;
    //			E_N_End = Ef + 3.01008e-4*T;
    //		}
    ////*******************************************************************************
    ////*******************************************************************************
    //
    //		for (int spin = 0; spin < 2; spin++)
    //		{
    //			Conductivity(sigma_coh[spin],sigma_comp[spin],
    //				sigmam, sigma_m, W, vm,	Pm, Cl, Hr, Nkp, H_Difr, Eps, Eps_sec, Et,
    //				Et_a, ee, N_f_int, E_N_Start,E_N_End,T,spin,dE,Ef,imm, Coh_p, r);
    //			Doub sigma_pkp = 0;
    //			Doub sigma[2] = {0};
    //			for (int i1 = 0; i1 < N_LAT; i1)
    //			{
    //				sigma_coh[spin][i1] /= N_LAT;
    //				sigma_pkp += sigma_coh[spin][i1].r;
    //				for (int j = 0; j < N_Com; j++)
    //				{
    //					for (int m = 0; m < 2; m++)
    //					{
    //						for (int ib = 0; ib < N_Band; ib++)
    //						{
    //							Doub sigmami;
    //							if (m == -1)
    //							{
    //								sigmami = sigma_m[j][i1][spin][ib].r;
    //							}
    //							else
    //							{
    //								sigmami = sigmam[j][i1][spin][ib].r;
    //							}
    //							sigma[spin] += Cl[j][i1]*Pm[m][i1]*sigmami/N_LAT;
    //						}
    //					}
    //				}
    //			}
    //			if (spin == 1)
    //			{
    //				ofstream f_sig_tot_start("H:\\MyBComp\\Study\\Clean Nanotube\\Progy\\mine\\Program\\Conduct.dat");
    //				f_sig_tot_start << E - Esmr << "    " << sigma_pkp << "    " << sigma_coh[0][0].r+sigma_coh[0][1].r << "    " <<
    //					sigma_coh[1][0].r+sigma_coh[1][1].r << "    " << sigma[0]+sigma[1] << "    " << sigma[0] << "    " <<
    //					sigma[1] << "    " << sigma_comp[0][0].r << "    " << sigma_comp[1][0].r << "    " << sigma_comp[0][1].r << "    " <<
    //					sigma_comp[1][1].r << "    " << sigma_comp[0][2].r << "    " << sigma_comp[1][2].r << "    " <<
    //					sigma_comp[0][3].r << "    " << sigma_comp[1][3].r << "    " << sigma_comp[0][4].r << "    " <<
    //					sigma_comp[1][4].r << "    " << sigma_comp[0][5].r << "    " << sigma_comp[1][5].r << "    " <<
    //					sigma_comp[0][6].r << "    " << sigma_comp[1][6].r << "    " << sigma_comp[0][7].r << "    " <<
    //					sigma_comp[1][7].r << "    " << sigma_comp[0][8].r << "    " << sigma_comp[1][8].r << '\n';
    //				f_sig_tot_start.close();
    //				cout << "Conductivity calculated";
    //			}
    //		}

    delete [] H_k;
    delete [] H_Difr;
    delete [] Gs;
    delete [] Sig_ef_e;
    delete [] Sig_fe;

    delete [] r;
    delete [] nm_1;
    delete [] mm_1;
    delete [] ni0;
    delete [] mm0;
    delete [] vm;
    delete [] Coh_p;
    delete [] W;
    delete [] Pm;
    delete [] Cl;
    delete [] sigmam;
    delete [] sigma_m;
    delete [] sigma_coh;
    delete [] sigma_comp;
}
void Lab::set_matrixInputFile(QString str) {
    matrixInputFile = str;
}

void Lab::set_impurityInputFile(QString str) {
    impurityInputFile = str;
}
