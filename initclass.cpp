#include "initclass.h"

InitClass::InitClass()
{
}

void InitClass::Set_coord()
{
        double a, a1, a2, r2, a1_2, a2_2, r2_2;
        a=3.0;				//����� ����� ������������ ������
        a1=3.0/4.0;
        a2 = sqrt(3.0)/4.0;
        r2 = sqrt(3.0)/2.0;
        a1_2=a1*5.0*3.0;
        a2_2 = a2*5.0*3.0;
        r2_2 = r2*5.0*3.0;
        r[0][0][1]=0.5;
        r[0][1][1]=0.0;
        r[0][2][1]=-r2;
        r[1][0][1]=0.0;
        r[1][1][1]=a1;
        r[1][2][1]=-a2;
        r[2][0][1]=0.5;
        r[2][1][1]=a1;
        r[2][2][1]=a2;
        r[3][0][1]=0.0;
        r[3][1][1]=0.0;
        r[3][2][1]=r2;
        r[4][0][1]=0.5;
        r[4][1][1]=-a1;
        r[4][2][1]=a2;
        r[5][0][1]=0.0;
        r[5][1][1]=-a1;
        r[5][2][1]=-a2;
        r[6][0][1]=3.0/2.0;
        r[6][1][1]=0.0;
        r[6][2][1]=-r2;
        r[7][0][1]=2.0;
        r[7][1][1]=a1;
        r[7][2][1]=-a2;
        r[8][0][1]=3.0/2.0;
        r[8][1][1]=a1;
        r[8][2][1]=a2;
        r[9][0][1]=2.0;
        r[9][1][1]=0.0;
        r[9][2][1]=r2;
        r[10][0][1]=3.0/2.0;
        r[10][1][1]=-a1;
        r[10][2][1]=a2;
        r[11][0][1]=2.0;
        r[11][1][1]=-a1;
        r[11][2][1]=-a2;
//        r[12][0][1]=0.5;
//	r[12][1][1]=0.0;
//	r[12][2][1]=-r2_2;
//	r[13][0][1]=0.0;
//	r[13][1][1]=a1_2;
//	r[13][2][1]=-a2_2;
//	r[14][0][1]=0.5;
//	r[14][1][1]=a1_2;
//	r[14][2][1]=a2_2;
//	r[15][0][1]=0.0;
//	r[15][1][1]=0.0;
//	r[15][2][1]=r2_2;
//	r[16][0][1]=0.5;
//	r[16][1][1]=-a1_2;
//	r[16][2][1]=a2_2;
//	r[17][0][1]=0.0;
//	r[17][1][1]=-a1_2;
//	r[17][2][1]=-a2_2;
//	r[18][0][1]=3.0/2.0;
//	r[18][1][1]=0.0;
//	r[18][2][1]=-r2_2;
//	r[19][0][1]=2.0;
//	r[19][1][1]=a1_2;
//	r[19][2][1]=-a2_2;
//	r[20][0][1]=3.0/2.0;
//	r[20][1][1]=a1_2;
//	r[20][2][1]=a2_2;
//	r[21][0][1]=2.0;
//	r[21][1][1]=0.0;
//	r[21][2][1]=r2_2;
//	r[22][0][1]=3.0/2.0;
//	r[22][1][1]=-a1_2;
//	r[22][2][1]=a2_2;
//	r[23][0][1]=2.0;
//	r[23][1][1]=-a1_2;
//	r[23][2][1]=-a2_2;
//**********************************************************

    for(int l = 0; l < 3; l++)
    {
        for(int i = 0; i < N_LAT; i++)
        {
            r[i][0][l]=r[i][0][1]+a*(l-1);
            r[i][1][l]=r[i][1][1];
            r[i][2][l]=r[i][2][1];
        }
    }
}

void InitClass::Read_Coulomb_integral()
{
    QFile file_desc("Coulomb/C_ee.dat");
    if (!file_desc.open(QFile::ReadOnly)) {
        cout << "Can't find file for Hop_integrals" << '\n';
        return;
    }

    QTextStream Coulomb_file(&file_desc);

    int is;
    Coulomb_file.readLine();
    Coulomb_file.readLine();
    Coulomb_file >> is >> is >> U_matr[0];
    Coulomb_file >> is >> is >> U_matr[1];
    U_matr[2] = U_matr[1];
    U_matr[3] = U_matr[1];
    Coulomb_file >> is >> is >> U_matr[4];
    U_matr[5] = U_matr[4];
    U_matr[6] = U_matr[4];
    U_matr[7] = U_matr[4];
    U_matr[8] = U_matr[4];
    Coulomb_file >> is >> is >> U_matr[9];

    file_desc.close();
}
void InitClass::Read_input_file() {
    for(int i1 = 0; i1 < N_LAT; i1++) {
        for(int j = 0; j < N_Com; j++) {
            for(int m = 0; m < 2; m++) {
                for(int ib = 0; ib < N_Band; ib++) {
                    if (ib == 0) {
                        mm_1[m][j][i1][ib] = (2*static_cast<double>(m)-1)*mm0[j][0];
                        nm_1[m][j][i1][ib] = ni0[j][0];
                    }
                    else {
                        if(ib < 4 && ib >= 1) {
                            mm_1[m][j][i1][ib] = (2*static_cast<double>(m)-1)*mm0[j][1]/3.0;
                            nm_1[m][j][i1][ib] = ni0[j][1]/3.0;
                        }
                        else {
                            if(ib >=4 && ib < 9) {
                                mm_1[m][j][i1][ib] = (2*static_cast<double>(m)-1)*mm0[j][2]/5.0;
                                nm_1[m][j][i1][ib] = ni0[j][2]/5.0;
                            }
                            else {
                                mm_1[m][j][i1][ib] = (2*static_cast<double>(m)-1)*mm0[j][3];
                                nm_1[m][j][i1][ib] = ni0[j][3];
                            }
                        }
                    }
                }
            }
        }
    }
}
void InitClass::Scat_Pr() {
//    Doub _eps = 0.0;
    ni0[0][0] = 2;
    ni0[0][1] = 2;
    ni0[0][2] = 0;
    ni0[0][3] = 0;
    ni0[1][0] = 2;
    ni0[1][1] = 3;
    ni0[1][2] = 0;
    ni0[1][3] = 0;
    mm0[0][0] = 0;
    mm0[0][1] = 0;
    mm0[0][2] = 0;
    mm0[0][3] = 0;
    mm0[1][0] = 0;
    mm0[1][1] = 1;
    mm0[1][2] = 0;
    mm0[1][3] = 0;
    double x2 = 0.05;
    for(int i = 0; i < N_LAT; i++) {
//        if (i < 12) {
//            Pm[0][i] = 1.0;
//            Pm[1][i] = 0.0;
//            Pm[2][i] = 0.0;
//            Cl[0][i] = 1.0;
//            Cl[1][i] = 0.0;
//            Cl[2][i] = 0.0;
//        } else {

    Pm[0][i] = 1.0-0.5;
    Pm[1][i] = 0.5;
    Cl[0][i] = 1.0-x2;
    Cl[1][i] = x2;

//            Pm[0][i] = 0.0;
//            Pm[1][i] = x2;
//            Pm[2][i] = 1.0-x2;
//            Cl[0][i] = 0.0;
//            Cl[1][i] = x2;
//            Cl[2][i] = 1.0-x2;
//        }

//        for (int j = 0; j < N_LAT; j++) {
//            for (int l = 0; l < N_LAT; l++) {
//                for (int l1 = 0; l1 < N_LAT; l1++) {
//                    Doub d1 = l1 == 1?1.0:0.0;
//                    Doub d2 = l1 == 2?1.0:0.0;
//                    Doub d3 = l == 1?1.0:0.0;
//                    Doub d4 = l == 2?1.0:0.0;
//                    if (j < 12) {
//                        Pmm[0][l1][j][i] = 1.0;
//                        Pmm[1][l1][j][i] = 0.0;
//                        Pmm[2][l1][j][i] = 0.0;
//                    } else {
//                        Pmm[0][l1][j][i] = 0.0;
//                    }
//                    if (i >= 12 && j >= 12 && l > 0 && l1 > 0) {
//                        Pmm[l][l1][j][i] = Pm[l1][j] + _eps/Pm[l][i]*(d1-d2)*(d3-d4);
//                    }
//                }
//            }
//        }

        }
}
void InitClass::allocate() {
    l_CPA = 0;
    r = new Doub[N_LAT][3][3];
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
    Hr = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    H_k = new Comp[Nkp2][N_LAT][N_LAT][N_Band][N_Band];
    S_k = new Comp[Nkp2][N_LAT][N_LAT][N_Band][N_Band];
    D_k = new Comp[Nkp2][N_LAT][N_LAT][N_Alpha][N_Alpha];
    H_Difr = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
    Sig_ef_e = new Comp[N_En][N_LAT][N_Band][N_Band];
    Sig_fe = new Comp[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha];
    Hr = new Comp[3][3][N_LAT][N_LAT][N_Band][N_Band];
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
    En_S = -4.7;
    En_E = 0.0;
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
