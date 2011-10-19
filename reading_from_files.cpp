#include "sys_param.h"
void Set_coord(Doub* r) {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;

    r[3] = 1.0;
    r[4] = 0.0;
    r[5] = 0.0;

//    r[6] = 0.5;
//    r[7] = -0.866;
//    r[8] = 0.3;

    translationVectors[0] = 1.5;
    translationVectors[1] = 0.866;
    translationVectors[2] = 0.0;

    translationVectors[3] = 0.0;
    translationVectors[4] = -1.732;
    translationVectors[5] = 0;

}


void Read_Coulomb_integral(Doub U_matr[N_Band])
{
    ifstream Coulomb_file("/home/veal/Sandbox/GUINano/Rezults/Coulomb/C_ee.dat");
    if (!Coulomb_file.is_open()) {
        cout << "Can't find file for Hop_integrals" << '\n';
        return;
    }

//    QTextStream Coulomb_file(&file_desc);

    string line;
    int is;
    getline(Coulomb_file, line);
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

    Coulomb_file.close();
}
void Read_input_file(Doub ni0[2][4], Doub mm0[2][4], Doub nm_1[2][N_Com][N_LAT][N_Band], Doub mm_1[2][N_Com][N_LAT][N_Band])
{
        for(int i1 = 0; i1 < N_LAT; i1++)
        {
                for(int j = 0; j < N_Com; j++)
                {
                        for(int m = 0; m < 2; m++)
                        {
                                for(int ib = 0; ib < N_Band; ib++)
                                {
                                        if (ib == 0)
                                        {
                                                mm_1[m][j][i1][ib] = (2*static_cast<Doub>(m)-1)*mm0[j][0];
                                                nm_1[m][j][i1][ib] = ni0[j][0];
                                        }
                                        else
                                        {
                                                if(ib < 4 && ib >= 1)
                                                {
                                                        mm_1[m][j][i1][ib] = (2*static_cast<Doub>(m)-1)*mm0[j][1]/3.0;
                                                        nm_1[m][j][i1][ib] = ni0[j][1]/3.0;
                                                }
                                                else
                                                {
                                                        if(ib >=4 && ib < 9)
                                                        {
                                                                mm_1[m][j][i1][ib] = (2*static_cast<Doub>(m)-1)*mm0[j][2]/5.0;
                                                                nm_1[m][j][i1][ib] = ni0[j][2]/5.0;
                                                        }
                                                        else
                                                        {
                                                                mm_1[m][j][i1][ib] = (2*static_cast<Doub>(m)-1)*mm0[j][3];
                                                                nm_1[m][j][i1][ib] = ni0[j][3];
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}
void Scat_Pr(Doub ni0[2][4], Doub mm0[2][4], Doub Pm[2][N_LAT], Doub Cl[N_Com][N_LAT]) {
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
    Doub x2 = 0.05;
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
