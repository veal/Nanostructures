#include "W_matrix.h"

W_matrix::W_matrix(PWF* matrixPWF, PWF* impurityPWF) {

    _dBasisWF = new double[4][N_tot];
    _dteta = new double[N_Max];
    _dteta2 = new double[N_tot][N_Max];
    _dR2 = new double[N_tot][N_Max];
    _mk_atom2 = new int[N_tot][N_Max];
    d_DeltaR = new double[N_tot];
    _dR = new double[N_tot];
    PARAM_FILE = "Rezults/W_params.txt";
    _dstep_teta = Pi/(N_Max-1);

    _matrixPWF = matrixPWF;
    setBasis(_matrixPWF);
    _impurityPWF = impurityPWF;
    if (_matrixPWF == NULL || _impurityPWF == NULL) {
        cout << "Failed to read files with data for W_matrix\n";
    }
    setGrid(_matrixPWF->nuclearCharge);
}

int W_matrix::get_integrals(double** W_container, double** Wsh_container,
        const double distance) {
    for (int iter = 0; iter < integralList.size(); iter++) {
        if (abs(integralList[iter].distance - distance) <= 0.01*distance) {
            *W_container = integralList[iter].W_container;
            *Wsh_container = integralList[iter].Wsh_container;
            return 0;
        }
    }

    integralContainer newIntegral(distance);
    calculateIntegrals(&newIntegral, distance);

    integralList.push_back(newIntegral);
    for (int iter = 0; iter < integralList.size(); iter++) {
        if (abs(integralList[iter].distance - distance) <= 0.01*distance) {
            *W_container = integralList[iter].W_container;
            *Wsh_container = integralList[iter].Wsh_container;
            return 0;
        }
    }
    return 0;
}

int W_matrix::calculateIntegrals(integralContainer *intCon, double distance) {

    double lclDistance = distance / 5.29e-11;    

    anotherCoordSys(lclDistance);

    ofstream VSFile;
    VSFile.open("Rezults/W_container.dat", ios::app);

    VSFile << "Computation result :" << '\n';
    VSFile <<  "\tW\t" << "Distance = " << distance << '\n';
    VSFile.precision(5);
    double* cleanPotential = _matrixPWF->getScaledPotential(_matrixPWF->nuclearCharge);
    double* impurityPotential = _impurityPWF->getScaledPotential(_matrixPWF->nuclearCharge);

    for (int m_orb = 0; m_orb < 14; m_orb++) {
        int n, m, lm;
        if (m_orb == 0) {
            n = 0;
            m = 0;
            lm = 0;
        }
        if (m_orb == 1) {
            n = 1;
            m = 1;
            lm = 0;
        }
        if (m_orb == 2) {
            n = 1;
            m = 1;
            lm = 1;
        }
        if (m_orb == 3) {
            n = 2;
            m = 2;
            lm = 0;
        }
        if (m_orb == 4) {
            n = 2;
            m = 2;
            lm = 1;
        }
        if (m_orb == 5) {
            n = 2;
            m = 2;
            lm = 2;
        }
        if (m_orb == 6) {
            n = 0;
            m = 1;
            lm = 0;
        }
        if (m_orb == 7) {
            n = 0;
            m = 2;
            lm = 0;
        }
        if (m_orb == 8) {
            n = 1;
            m = 2;
            lm = 0;
        }
        if (m_orb == 9) {
            n = 1;
            m = 2;
            lm = 1;
        }
        if (m_orb == 10) {
            n = 3;
            m = 0;
            lm = 0;
        }
        if (m_orb == 11) {
            n = 3;
            m = 1;
            lm = 0;
        }
        if (m_orb == 12) {
            n = 3;
            m = 2;
            lm = 0;
        }
        if (m_orb == 13) {
            n = 3;
            m = 3;
            lm = 0;
        }
        double d_Int2 = 0, d_Int3 = 0;
//        double Max = -1e+6, Min = 1e+6;
//        int imin, jmin;
//        int imax, jmax;
        for (int i = 1; i < N_tot; i++) {
            for (int j = 1; j < N_Max; j++) {
//                if (FunV(n,m,i,j) < Min) {
//                    imin = i; jmin = j;
//                    Min = FunV(n,m,i,j);
//                }
//                if (FunV(n,m,i,j) > Max) {
//                    imax = i; jmax = j;
//                    Max = FunV(n,m,i,j);
//                }    
                d_Int2 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunW(cleanPotential,n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i-1][j-1])+
                FunW(cleanPotential,n,m,i-1,j)*PolLagr(m_orb,_dteta[j],_dteta2[i-1][j]) +
                FunW(cleanPotential,n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i][j-1]) +
                FunW(cleanPotential,n,m,i,j)*PolLagr(m_orb,_dteta[j],_dteta2[i][j]));

                d_Int3 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunW(impurityPotential,n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i-1][j-1])+
                FunW(impurityPotential,n,m,i-1,j)*PolLagr(m_orb,_dteta[j],_dteta2[i-1][j]) +
                FunW(impurityPotential,n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i][j-1]) +
                FunW(impurityPotential,n,m,i,j)*PolLagr(m_orb,_dteta[j],_dteta2[i][j]));
            }
//            cout << FunW(impurityPotential, 2, 2, i, 220) - FunW(cleanPotential, 2, 2, i, 220) << '\n'; //!!!!!!
        }

        intCon->W_container[m_orb] = d_Int3 - d_Int2;

        d_Int3 = 0;
        d_Int2 = 0;

        for (int i = 1; i < N_tot; i++) {
            for (int j = 1; j < N_Max; j++) {
                d_Int2 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunWsh(cleanPotential,n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i-1][j-1])+
                FunWsh(cleanPotential,n,m,i-1,j)*PolLagr(m_orb,_dteta[j],_dteta2[i-1][j]) +
                FunWsh(cleanPotential,n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i][j-1]) +
                FunWsh(cleanPotential,n,m,i,j)*PolLagr(m_orb,_dteta[j],_dteta2[i][j]));

                d_Int3 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunWsh(impurityPotential,n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i-1][j-1])+
                FunWsh(impurityPotential,n,m,i-1,j)*PolLagr(m_orb,_dteta[j],_dteta2[i-1][j]) +
                FunWsh(impurityPotential,n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i][j-1]) +
                FunWsh(impurityPotential,n,m,i,j)*PolLagr(m_orb,_dteta[j],_dteta2[i][j]));
            }
        }

        intCon->Wsh_container[m_orb] = d_Int3 - d_Int2;

        VSFile << n << '\t' << m << '\t' << lm << '\t' << scientific << intCon->W_container[m_orb] <<
                '\t' << intCon->Wsh_container[m_orb] << '\n';
//        VSFile << m_orb << "-------------------------------------------------------" << '\n';
//        VSFile << imin << '\t' << jmin << '\t' << scientific << Min << '\n';
//        VSFile << imax << '\t' << jmax << '\t' << scientific << Max << '\n';  
    }
    VSFile.close();

    return 0;
}

double W_matrix::FunW(double* pot, int n_, int m_, int i_, int j_) {
    return _dBasisWF[n_][i_] * _dBasisWF[m_][_mk_atom2[i_][j_]]
            * _dR[i_] * _dR[i_] * sin(_dteta[j_])*(pot[i_]);
}

double W_matrix::FunWsh(double* pot, int n_, int m_, int i_, int j_) {
    return _dBasisWF[n_][i_] * _dBasisWF[m_][_mk_atom2[i_][j_]]
            * _dR[i_] * _dR[i_] * sin(_dteta[j_])*(pot[_mk_atom2[i_][j_]]);
}

double W_matrix::PolLagr(int orb, double teta, double teta2) {
        switch (orb) {
            case 0:
                return 1.0 / 2.0;
            case 1:
                return (3.0 / 2.0) * cos(teta) * cos(teta2);
            case 2: return (3.0 / 4.0) * sin(teta) * sin(teta2);
            case 3: return (5.0 / 8.0)*(3.0 * cos(teta) * cos(teta) - 1)*
                        (3.0 * cos(teta2) * cos(teta2) - 1);
            case 4: return (15.0 / 16.0) * sin(2.0 * teta) * sin(2.0 * teta2);
            case 5: return (15.0 / 16.0) * sin(teta) * sin(teta) * sin(teta2) * sin(teta2);
            case 6: return (sqrt(3.0) / 2.0) * cos(teta2);
            case 7: return (sqrt(5.0) / 4.0)*
                        (3.0 * cos(teta2) * cos(teta2) - 1);
            case 8: return (sqrt(15.0) / 4.0) *
                        cos(teta)*(3.0 * cos(teta2) * cos(teta2) - 1);
            case 9: return (sqrt(45.0) / 8.0) * sin(teta) * sin(2.0 * teta2);
            case 10: return 1.0 / 2.0;
            case 11: return (sqrt(3.0) / 2.0) * cos(teta2);
            case 12: return (sqrt(5.0) / 4.0)*(3.0 * cos(teta2) * cos(teta2) - 1);
            case 13: return 1.0 / 2.0;
        }
    }

W_matrix::W_matrix(const W_matrix& orig) {
}

W_matrix::~W_matrix() {
    for (int iter = 0; iter < integralList.size(); iter++) {
        delete[] integralList[iter].W_container;
        delete[] integralList[iter].Wsh_container;
    }
    delete[] d_DeltaR;
    delete[] _dR;
    delete[] _dBasisWF;
    delete[] _dteta;
    delete[] _dteta2;
    delete[] _dR2;
    delete[] _mk_atom2;
}

void W_matrix::setGrid(int m_Z) {
    int m_Nblock = N_tot/40;
    int iter = 0;
    double d_C = 0.88534138/pow(m_Z, 1.0/3.0);
    double d_DeltaX = 0.0025;
    
    ofstream file(PARAM_FILE.c_str(), ios::app);

    file << "Nuclear charge = " << m_Z << '\n'
            << "Number of blocks = " << m_Nblock << '\n';

    file << "\td_X\t" << "\td_R\t" << "\td_DeltaR" << '\n';
    file.precision(5);
    double* d_X = new double[N_tot];

    for (int j = 0; j < m_Nblock; j++) {
        for (int k = 0; k < 40; k++) {
            iter++;
            d_X[iter] = d_X[iter-1] + d_DeltaX;
            _dR[iter-1] = d_X[iter-1]*d_C;
            d_DeltaR[iter] = d_C*d_DeltaX;
            file << iter-1 << '\t' << scientific << d_X[iter-1] << '\t' << _dR[iter-1] << '\t' <<
                    d_DeltaR[iter] << '\n';
        }
        d_DeltaX += d_DeltaX;
    }

    _dR[iter] = d_X[iter]*d_C;
    file << iter << '\t' << scientific << d_X[iter] << '\t' << _dR[iter] << '\t' <<
                    d_DeltaR[iter] << '\n';
    file.close();

    delete[] d_X;
}

void W_matrix::anotherCoordSys(double distance) {
//    ofstream VSFile("Rezults/V_S_container.dat", ios::app);
//    VSFile << "Another coord sys indices" << '\n';
    for (int i = 1; i < N_Max; i++) {
        _dteta[i] = _dteta[i-1] + _dstep_teta;
    }

    for (int i = 0; i < N_tot; i++) {
        for (int j = 0; j < N_Max; j++) {
            if (j == 0) {
                _dR2[i][j] = abs(_dR[i] - distance);
            } else {
                _dR2[i][j] = pow((_dR[i]*_dR[i] + distance*distance -
                    2.0*_dR[i]*distance*cos(_dteta[j])), 0.5);
                if (distance == 0) {
                    _dteta2[i][j] = _dteta[j];
                } else {
                    double d_tmp = (_dR2[i][j]*_dR2[i][j] + distance*distance - _dR[i]*_dR[i])
                        /(2.0*_dR2[i][j]*distance);
                    if (abs(d_tmp) - 1.0 > 0.0)
                        d_tmp /= abs(d_tmp);
                    _dteta2[i][j] = acos(d_tmp);
                    _dteta2[i][j] = Pi - _dteta2[i][j];
                }
            }
            double d_max_number = 10000;
            for(int tin = 0; tin < N_tot; tin++) {
                if (abs(_dR[tin] - _dR2[i][j]) < d_max_number) {
                    d_max_number = abs(_dR[tin] - _dR2[i][j]);
                    _mk_atom2[i][j] = tin;
                }
            }
//            VSFile << i << '\t' << j << '\t' << _dteta2[i][j] << '\n';
        }
    }
//    VSFile.close();
}

void W_matrix::setBasis(PWF* pwf) {
    double *tmp0, *tmp1, *tmp2, *tmp3;
    tmp0 = pwf->getScaledWF0(_matrixPWF->nuclearCharge);
    tmp1 = pwf->getScaledWF1(_matrixPWF->nuclearCharge);
    tmp2 = pwf->getScaledWF2(_matrixPWF->nuclearCharge);
    tmp3 = pwf->getScaledWF3(_matrixPWF->nuclearCharge);
    for (int iter = 0; iter < N_tot; iter++) {
        _dBasisWF[0][iter] = tmp0[iter];
        _dBasisWF[1][iter] = tmp1[iter];
        _dBasisWF[2][iter] = tmp2[iter];
        _dBasisWF[3][iter] = tmp3[iter];
    }
}
void W_matrix::F_Wk(Doub kx, Comp H[N_LAT][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]) {
    Doub *V1, *V2, *V, *Vsh;
    Doub a, d, Pi, x, y, z;

    a = 3.0; //����� ������������ ������
    Pi = 3.1415926536;

    for (int it = 0; it < N_LAT; it++) {
        for (int i = 0; i < N_LAT; i++) {
            for (int j = 0; j < N_LAT; j++) {
//                for (int n = 0; n < 3; n++) {

//                    x = (r[j][0][n] - r[i][0][1]);
//                    y = (r[j][1][n] - r[i][1][1]);
//                    z = (r[j][2][n] - r[i][2][1]);

                    x = (r[j][0][1] - r[i][0][1]);
                    y = (r[j][1][1] - r[i][1][1]);
                    z = (r[j][2][1] - r[i][2][1]);

                    d = sqrt(x * x + y * y + z * z);

                    if (d < 2.1 && d > 0.0) {

                        get_integrals(&V, &Vsh, 1.42e-10 * d);

                        if (it == i) {
                            V1 = V;
                            V2 = Vsh;
                        } else if (it == j) {
                            V2 = V;
                            V1 = Vsh;
                        } else {
                            break;
                        }

                        H[it][i][j][0][0] += V1[0] * exp(CI * kx * x);
                        H[it][i][j][0][1] += (x / d) * V1[6] * exp(CI * kx * x);
                        H[it][i][j][0][2] += (y / d) * V1[6] * exp(CI * kx * x);
                        H[it][i][j][0][3] += (z / d) * V1[6] * exp(CI * kx * x);
                        H[it][i][j][0][4] += 1.732 * (x / d)*(y / d) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][0][5] += 1.732 * (y / d)*(z / d) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][0][6] += 1.732 * (z / d)*(x / d) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][0][7] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][0][8] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][0][9] += V1[10] * exp(CI * kx * x);
                        H[it][i][j][1][0] -= (x / d) * V2[6] * exp(CI * kx * x);   ///!!!!!!
                        H[it][i][j][1][1] += ((x / d)*(x / d) * V1[1] + (1 - (x / d)*(x / d)) * V1[2]) * exp(CI * kx * x);
                        H[it][i][j][1][2] += (x / d)*(y / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                        H[it][i][j][1][3] += (x / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                        H[it][i][j][1][4] += (1.732 * (x / d)*(x / d)*(y / d) * V1[8]+(y / d)*(1 - 2 * (x / d)*(x / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][1][5] += (x / d)*(y / d)*(z / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][1][6] += (1.732 * (x / d)*(x / d)*(z / d) * V1[8]+(z / d)*(1 - 2 * (x / d)*(x / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][1][7] += (0.5 * 1.732 * (x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]+
                                (x / d)*(1 - (x / d)*(x / d)+ (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][1][8] += ((x / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8]
                                - 1.732 * (x / d)*(z / d)*(z / d) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][1][9] -= (x / d) * V2[11] * exp(CI * kx * x);
                        H[it][i][j][2][0] -= (y / d) * V2[6] * exp(CI * kx * x);
                        H[it][i][j][2][1] += (x / d)*(y / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                        H[it][i][j][2][2] += ((y / d)*(y / d) * V1[1] + (1 - (y / d)*(y / d)) * V1[2]) * exp(CI * kx * x);
                        H[it][i][j][2][3] += (y / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                        H[it][i][j][2][4] += (1.732 * (y / d)*(y / d)*(x / d) * V1[8]+(x / d)*(1 - 2 * (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][2][5] += (1.732 * (y / d)*(y / d)*(z / d) * V1[8]+(z / d)*(1 - 2 * (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][2][6] += (y / d)*(z / d)*(x / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][2][7] += (0.5 * 1.732 * (y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]-
                                (y / d)*(1 + (x / d)*(x / d)-(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][2][8] += ((y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8] -
                                1.732 * (y / d)*(z / d)*(z / d) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][2][9] -= (y / d) * V2[11] * exp(CI * kx * x);
                        H[it][i][j][3][0] -= (z / d) * V2[6] * exp(CI * kx * x);
                        H[it][i][j][3][1] += (x / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                        H[it][i][j][3][2] += (y / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                        H[it][i][j][3][3] += ((z / d)*(z / d) * V1[1] + (1 - (z / d)*(z / d)) * V1[2]) * exp(CI * kx * x);
                        H[it][i][j][3][4] += (x / d)*(y / d)*(z / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][3][5] += (1.732 * (z / d)*(z / d)*(y / d) * V1[8]+(y / d)*(1 - 2 * (z / d)*(z / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][3][6] += (1.732 * (z / d)*(z / d)*(x / d) * V1[8]+(x / d)*(1 - 2 * (z / d)*(z / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][3][7] += (0.5 * 1.732 * (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]-
                                (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][3][8] += ((z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8] +
                                1.732 * (z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                        H[it][i][j][3][9] -= (z / d) * V2[11] * exp(CI * kx * x);
                        H[it][i][j][4][0] += 1.732 * (x / d)*(y / d) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][4][1] -= (1.732 * (x / d)*(x / d)*(y / d) * V2[8]+(y / d)*(1 - 2.0 * (x / d)*(x / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][4][2] -= (1.732 * (y / d)*(y / d)*(x / d) * V2[8]+(x / d)*(1 - 2.0 * (y / d)*(y / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][4][3] -= (x / d)*(y / d)*(z / d)*(1.732 * V2[8] - 2.0 * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][4][4] += (3.0 * (x / d)*(x / d)*(y / d)*(y / d) * V1[3]+((x / d)*(x / d)+(y / d)*(y / d) - 4.0 * (x / d)*(x / d)*(y / d)*(y / d)) * V1[4]+
                                ((z / d)*(z / d)+(x / d)*(x / d)*(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][4][5] += (3.0 * (x / d)*(y / d)*(y / d)*(z / d) * V1[3]+(x / d)*(z / d)*(1 - 4 * (y / d)*(y / d)) * V1[4]+
                                (x / d)*(z / d)*((y / d)*(y / d) - 1) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][4][6] += (3 * (x / d)*(x / d)*(y / d)*(z / d) * V1[3]+(y / d)*(z / d)*(1 - 4.0 * (x / d)*(x / d)) * V1[4]+
                                (y / d)*(z / d)*((x / d)*(x / d) - 1) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][4][7] += (1.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3] + 2.0 * (x / d)*(y / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                                0.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][4][8] += (1.732 * (x / d)*(y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                2.0 * 1.732 * (x / d)*(y / d)*(z / d)*(z / d) * V1[4] + 0.5 * 1.732 * (x / d)*(y / d)*(1 + (z / d)*(z / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][4][9] += 1.732 * (x / d)*(y / d) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][5][0] += 1.732 * (y / d)*(z / d) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][5][1] -= (x / d)*(y / d)*(z / d)*(1.732 * V2[8] - 2.0 * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][5][2] -= (1.732 * (y / d)*(y / d)*(z / d) * V2[8]+(z / d)*(1 - 2 * (y / d)*(y / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][5][3] -= (1.732 * (z / d)*(z / d)*(y / d) * V2[8]+(y / d)*(1 - 2 * (z / d)*(z / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][5][4] += (3.0 * (x / d)*(y / d)*(y / d)*(z / d) * V1[3]+(x / d)*(z / d)*(1 - 4 * (y / d)*(y / d)) * V1[4]+
                                (x / d)*(z / d)*((y / d)*(y / d) - 1) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][5][5] += (3.0 * (y / d)*(y / d)*(z / d)*(z / d) * V1[3]+((y / d)*(y / d)+(z / d)*(z / d) - 4.0 * (y / d)*(y / d)*(z / d)*(z / d)) * V1[4]+
                                ((x / d)*(x / d)+(y / d)*(y / d)*(z / d)*(z / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][5][6] += (3.0 * (y / d)*(z / d)*(z / d)*(x / d) * V1[3]+(y / d)*(x / d)*(1 - 4 * (z / d)*(z / d)) * V1[4]+
                                (y / d)*(x / d)*((z / d)*(z / d) - 1) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][5][7] += (1.5 * (y / d)*(z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]-
                                (y / d)*(z / d)*(1 + 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]+
                                (y / d)*(z / d)*(1.0 + 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][5][8] += (1.732 * (y / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                                0.5 * 1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][5][9] += (1.732 * (y / d)*(z / d) * V1[7]) * exp(CI * kx * x);
                        H[it][i][j][6][0] += 1.732 * (z / d)*(x / d) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][6][1] -= (1.732 * (x / d)*(x / d)*(z / d) * V2[8]+(z / d)*(1 - 2 * (x / d)*(x / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][6][2] -= (y / d)*(z / d)*(x / d)*(1.732 * V2[8] - 2.0 * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][6][3] -= (1.732 * (z / d)*(z / d)*(x / d) * V2[8]+(x / d)*(1 - 2 * (z / d)*(z / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][6][4] += (3 * (x / d)*(x / d)*(y / d)*(z / d) * V1[3]+(y / d)*(z / d)*(1 - 4.0 * (x / d)*(x / d)) * V1[4]+
                                (y / d)*(z / d)*((x / d)*(x / d) - 1) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][6][5] += (3.0 * (y / d)*(z / d)*(z / d)*(x / d) * V1[3]+(y / d)*(x / d)*(1 - 4 * (z / d)*(z / d)) * V1[4]+
                                (y / d)*(x / d)*((z / d)*(z / d) - 1) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][6][6] += (3.0 * (z / d)*(z / d)*(x / d)*(x / d) * V1[3]+((z / d)*(z / d)+(x / d)*(x / d) - 4.0 * (z / d)*(z / d)*(x / d)*(x / d)) * V1[4]+
                                ((y / d)*(y / d)+(z / d)*(z / d)*(x / d)*(x / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][6][7] += (1.5 * (z / d)*(x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]+
                                (z / d)*(x / d)*(1.0 - 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]-
                                (z / d)*(x / d)*(1.0 - 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][6][8] += (1.732 * (x / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                                0.5 * 1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][6][9] += 1.732 * (z / d)*(x / d) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][7][0] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][7][1] -= (0.5 * 1.732 * (x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V2[8]+
                                (x / d)*(1 - (x / d)*(x / d) + (y / d)*(y / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][7][2] -= (0.5 * 1.732 * (y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V2[8]-
                                (y / d)*(1 + (x / d)*(x / d)-(y / d)*(y / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][7][3] -= (0.5 * 1.732 * (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V2[8]-
                                (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][7][4] += (1.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3] + 2.0 * (x / d)*(y / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                                0.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][7][5] += (1.5 * (y / d)*(z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]-
                                (y / d)*(z / d)*(1 + 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]+
                                (y / d)*(z / d)*(1.0 + 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][7][6] += (1.5 * (z / d)*(x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]+
                                (z / d)*(x / d)*(1.0 - 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]-
                                (z / d)*(x / d)*(1.0 - 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][7][7] += (0.75 * ((x / d)*(x / d)-(y / d)*(y / d))*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]+
                                ((x / d)*(x / d)+(y / d)*(y / d)-((x / d)*(x / d)-(y / d)*(y / d))*((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]+
                                ((z / d)*(z / d) + 0.25 * ((x / d)*(x / d)-(y / d)*(y / d))*((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][7][8] += (0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d))*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                1.732 * (z / d)*(z / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                                0.25 * 1.732 * (1 + (z / d)*(z / d))*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][7][9] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][8][0] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[7] * exp(CI * kx * x);
                        H[it][i][j][8][1] -= ((x / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V2[8]
                                - 1.732 * (x / d)*(z / d)*(z / d) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][8][2] -= ((y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V2[8] -
                                1.732 * (y / d)*(z / d)*(z / d) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][8][3] -= ((z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V2[8] +
                                1.732 * (z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V2[9]) * exp(CI * kx * x);
                        H[it][i][j][8][4] += (1.732 * (x / d)*(y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                2.0 * 1.732 * (x / d)*(y / d)*(z / d)*(z / d) * V1[4] + 0.5 * 1.732 * (x / d)*(y / d)*(1 + (z / d)*(z / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][8][5] += (1.732 * (y / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                                0.5 * 1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][8][6] += (1.732 * (x / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                                0.5 * 1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][8][7] += (0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d))*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                1.732 * (z / d)*(z / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                                0.25 * 1.732 * (1 + (z / d)*(z / d))*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][8][8] += (((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d)))*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                                3.0 * (z / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[4] +
                                0.75 * ((x / d)*(x / d)+(y / d)*(y / d))*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                        H[it][i][j][8][9] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][9][0] += V1[10] * exp(CI * kx * x);
                        H[it][i][j][9][1] += (x / d) * V1[11] * exp(CI * kx * x);
                        H[it][i][j][9][2] += (y / d) * V1[11] * exp(CI * kx * x);
                        H[it][i][j][9][3] += (z / d) * V1[11] * exp(CI * kx * x);
                        H[it][i][j][9][4] += 1.732 * (x / d)*(y / d) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][9][5] += (1.732 * (y / d)*(z / d) * V1[7]) * exp(CI * kx * x);
                        H[it][i][j][9][6] += 1.732 * (z / d)*(x / d) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][9][7] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][9][8] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[12] * exp(CI * kx * x);
                        H[it][i][j][9][9] += V1[13] * exp(CI * kx * x);
                    }
                    if (d == 0 && i == it && j == it) {

                        get_integrals(&V, &Vsh, 1.42e-10 * d);

                        V1 = V;

                        H[it][i][j][0][0] = V1[0];
                        H[it][i][j][0][1] = 0.0;
                        H[it][i][j][0][2] = 0.0;
                        H[it][i][j][0][3] = 0.0;
                        H[it][i][j][0][4] = 0.0;
                        H[it][i][j][0][5] = 0.0;
                        H[it][i][j][0][6] = 0.0;
                        H[it][i][j][0][7] = 0.0;
                        H[it][i][j][0][8] = 0.0;
                        H[it][i][j][0][9] = 0.0;
                        H[it][i][j][1][0] = 0.0;
                        H[it][i][j][1][1] = V1[1];
                        H[it][i][j][1][2] = 0.0;
                        H[it][i][j][1][3] = 0.0;
                        H[it][i][j][1][4] = 0.0;
                        H[it][i][j][1][5] = 0.0;
                        H[it][i][j][1][6] = 0.0;
                        H[it][i][j][1][7] = 0.0;
                        H[it][i][j][1][8] = 0.0;
                        H[it][i][j][1][9] = 0.0;
                        H[it][i][j][2][0] = 0.0;
                        H[it][i][j][2][1] = 0.0;
                        H[it][i][j][2][2] = V1[1];
                        H[it][i][j][2][3] = 0.0;
                        H[it][i][j][2][4] = 0.0;
                        H[it][i][j][2][5] = 0.0;
                        H[it][i][j][2][6] = 0.0;
                        H[it][i][j][2][7] = 0.0;
                        H[it][i][j][2][8] = 0.0;
                        H[it][i][j][2][9] = 0.0;
                        H[it][i][j][3][0] = 0.0;
                        H[it][i][j][3][1] = 0.0;
                        H[it][i][j][3][2] = 0.0;
                        H[it][i][j][3][3] = V1[1];
                        H[it][i][j][3][4] = 0.0;
                        H[it][i][j][3][5] = 0.0;
                        H[it][i][j][3][6] = 0.0;
                        H[it][i][j][3][7] = 0.0;
                        H[it][i][j][3][8] = 0.0;
                        H[it][i][j][3][9] = 0.0;
                        H[it][i][j][4][0] = 0.0;
                        H[it][i][j][4][1] = 0.0;
                        H[it][i][j][4][2] = 0.0;
                        H[it][i][j][4][3] = 0.0;
                        H[it][i][j][4][4] = V1[4];
                        H[it][i][j][4][5] = 0.0;
                        H[it][i][j][4][6] = 0.0;
                        H[it][i][j][4][7] = 0.0;
                        H[it][i][j][4][8] = 0.0;
                        H[it][i][j][4][9] = 0.0;
                        H[it][i][j][5][0] = 0.0;
                        H[it][i][j][5][1] = 0.0;
                        H[it][i][j][5][2] = 0.0;
                        H[it][i][j][5][3] = 0.0;
                        H[it][i][j][5][4] = 0.0;
                        H[it][i][j][5][5] = V1[4];
                        H[it][i][j][5][6] = 0.0;
                        H[it][i][j][5][7] = 0.0;
                        H[it][i][j][5][8] = 0.0;
                        H[it][i][j][5][9] = 0.0;
                        H[it][i][j][6][0] = 0.0;
                        H[it][i][j][6][1] = 0.0;
                        H[it][i][j][6][2] = 0.0;
                        H[it][i][j][6][3] = 0.0;
                        H[it][i][j][6][4] = 0.0;
                        H[it][i][j][6][5] = 0.0;
                        H[it][i][j][6][6] = V1[4];
                        H[it][i][j][6][7] = 0.0;
                        H[it][i][j][6][8] = 0.0;
                        H[it][i][j][6][9] = 0.0;
                        H[it][i][j][7][0] = 0.0;
                        H[it][i][j][7][1] = 0.0;
                        H[it][i][j][7][2] = 0.0;
                        H[it][i][j][7][3] = 0.0;
                        H[it][i][j][7][4] = 0.0;
                        H[it][i][j][7][5] = 0.0;
                        H[it][i][j][7][6] = 0.0;
                        H[it][i][j][7][7] = V1[4];
                        H[it][i][j][7][8] = 0.0;
                        H[it][i][j][7][9] = 0.0;
                        H[it][i][j][8][0] = 0.0;
                        H[it][i][j][8][1] = 0.0;
                        H[it][i][j][8][2] = 0.0;
                        H[it][i][j][8][3] = 0.0;
                        H[it][i][j][8][4] = 0.0;
                        H[it][i][j][8][5] = 0.0;
                        H[it][i][j][8][6] = 0.0;
                        H[it][i][j][8][7] = 0.0;
                        H[it][i][j][8][8] = V1[4];
                        H[it][i][j][8][9] = 0.0;
                        H[it][i][j][9][0] = 0.0;
                        H[it][i][j][9][1] = 0.0;
                        H[it][i][j][9][2] = 0.0;
                        H[it][i][j][9][3] = 0.0;
                        H[it][i][j][9][4] = 0.0;
                        H[it][i][j][9][5] = 0.0;
                        H[it][i][j][9][6] = 0.0;
                        H[it][i][j][9][7] = 0.0;
                        H[it][i][j][9][8] = 0.0;
                        H[it][i][j][9][9] = V1[13];
                    }
//                }
            }
        }
    }
}

void W_matrix::Wk(int Nkp, Comp Wr[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band],
        Comp S_k[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3], PWF* _matrixPWF, PWF* _impurityPWF) {

    W_matrix WIntegral(_matrixPWF, _impurityPWF);

    for (int k = -Nkp; k < Nkp; k++) {
        Comp (*W_temp)[N_LAT][N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_LAT][N_Band][N_Band];
        cout << k << '\n';

        Doub kx = 2 * Pi * static_cast<Doub> (k) / (3 * Nkp);

        F_Wk(kx, W_temp, r);

        for (int i1 = 0; i1 < N_Band; i1++) {
            for (int j1 = 0; j1 < N_Band; j1++) {
                for (int i = 0; i < N_LAT; i++) {
                    for (int ish = 0; ish < N_LAT; ish++) {
                        for (int iter = 1; iter < N_LAT; iter++) {
                            W_temp[0][i][ish][i1][j1] += W_temp[iter][i][ish][i1][j1];
//                            cout << iter << '\t' << i << '\t' << ish << '\t' << i1 <<
//                                    '\t' << j1 << '\t' << W[1][iter][i][i][i1][j1] << '\n';
                        }
//                            cout << i << '\t' << ish << '\t' << i1 <<
//                                    '\t' << j1 << '\t' << W_temp[0][i][ish][i1][j1] << '\n';
                    }
                }
            }
        }

//        for (int iter = 0; iter < 1; iter++) {
//            Trans_Hk(S_k[k + Nkp], W_temp[iter], W_temp[iter], false);
//        }


        for (int i1 = 0; i1 < N_Band; i1++) {
            for (int j1 = 0; j1 < N_Band; j1++) {
                for (int i = 0; i < N_LAT; i++) {
                    for (int ish = 0; ish < N_LAT; ish++) {
                        cout << i << '\t' << ish << '\t' << i1 <<
                                '\t' << j1 << '\t' << W_temp[0][i][ish][i1][j1] << '\n';

                    }
                }
            }
        }
        break;
        Doub ro_loc = ((k == -Nkp) || (k == Nkp - 1)) ? 0.5 : 1;

        for (int iter = 0; iter < N_LAT; iter++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    for (int i = 0; i < N_LAT; i++) {
                        for (int j = 0; j < N_LAT; j++) {
                            Comp p = exp(CI * kx * (r[i][0][1] - r[j][0][1]));
                            Wr[1][iter][i][j][i1][j1] += W_temp[iter][i][j][i1][j1] * p * ro_loc / (static_cast<Doub> (Nkp2));
                            Wr[0][iter][i][j][i1][j1] = 0;
                        }
                    }
                }
            }
        }
        delete [] W_temp;
    }
}

void W_matrix::Calculate_Wamiltonian(PWF* _matrixPWF, PWF* _impurityPWF) {
    if (calculateWamiltonian == true) {
        ofstream Ham_r_file;
        Ham_r_file.open("Rezults/W_k.dat", ios::binary);

//        Wk(Nkp, W, S_k, r, _matrixPWF, _impurityPWF);

        Ham_r_file.write((char*) W, N_Com * N_LAT * N_LAT * N_LAT * N_Band * N_Band * sizeof (Comp));
        Ham_r_file.close();
    } else {
        cout << "Reading Wamiltonian.." << '\n';
        ifstream Ham_r_file;
        Ham_r_file.open("Rezults/W_k.dat", ios::binary);
        Ham_r_file.read((char*) W, N_Com * N_LAT * N_LAT * N_LAT * N_Band * N_Band * sizeof (Comp));
        Ham_r_file.close();

//        for (int iter = 0; iter < N_LAT; iter++) {
//            for (int i1 = 0; i1 < N_Band; i1++) {
//                for (int j1 = 0; j1 < N_Band; j1++) {
//                    for (int i = 0; i < N_LAT; i++) {
//                        for (int ish = 0; ish < N_LAT; ish++) {
//                            cout << iter << '\t' << i << '\t' << ish << '\t' << i1 <<
//                                    '\t' << j1 << '\t' << W[1][iter][i][i][i1][j1] << '\n';
//                        }
//                    }
//                }
//            }
//        }

    }
}
