#include "W_matrix.h"

W_matrix::W_matrix(PWF* matrixPWF, PWF* impurityPWF) {

    _dBasisWF = new double[4][N_tot];
    _dteta = new double[N_Max];
    _dteta2 = new double[N_tot][N_Max];
    _dR2 = new double[N_tot][N_Max];
    _mk_atom2 = new int[N_tot][N_Max];
    d_DeltaR = new double[N_tot];
    _dR = new double[N_tot];
    PARAM_FILE = "/home/veal/Sandbox/GUINano/Rezults/W_params.txt";
    WF_FILE = "/home/veal/Sandbox/GUINano/Rezults/";
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

    QString str;
    QTextStream sstr(&str);
    sstr << WF_FILE << "W_container.dat";
    QFile VSFile(str);

    if (!VSFile.open(QFile::WriteOnly | QFile::Append)) {
        QString temp = "Can't find pwf file :";
        temp += "Rezults/W_container.dat";
        cout << temp.toStdString() << '\n';
        return NULL;
    }

    QTextStream in(&VSFile);
    in.setRealNumberNotation(QTextStream::ScientificNotation);

    in << "Computation result :" << '\n';
    in <<  "\tW\t" << "Distance = " << distance << '\n';

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
        for (int i = 1; i < N_tot; i++) {
            for (int j = 1; j < N_Max; j++) {   
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

        in << n << '\t' << m << '\t' << lm << '\t' << scientific << intCon->W_container[m_orb] <<
                '\t' << intCon->Wsh_container[m_orb] << '\n'; 
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

    QFile VSFile(PARAM_FILE);

    if (!VSFile.open(QFile::WriteOnly | QFile::Append)) {
        QString temp = "Can't find pwf file :";
        temp += "Rezults/W_param.dat";
        cout << temp.toStdString() << '\n';
        return NULL;
    }

    QTextStream in(&VSFile);
    in.setRealNumberNotation(QTextStream::ScientificNotation);

    in << "Nuclear charge = " << m_Z << '\n'
            << "Number of blocks = " << m_Nblock << '\n';

    in << "\td_X\t" << "\td_R\t" << "\td_DeltaR" << '\n';

    double* d_X = new double[N_tot];

    for (int j = 0; j < m_Nblock; j++) {
        for (int k = 0; k < 40; k++) {
            iter++;
            d_X[iter] = d_X[iter-1] + d_DeltaX;
            _dR[iter-1] = d_X[iter-1]*d_C;
            d_DeltaR[iter] = d_C*d_DeltaX;
            in << iter-1 << '\t' << scientific << d_X[iter-1] << '\t' << _dR[iter-1] << '\t' <<
                    d_DeltaR[iter] << '\n';
        }
        d_DeltaX += d_DeltaX;
    }

    _dR[iter] = d_X[iter]*d_C;
    in << iter << '\t' << scientific << d_X[iter] << '\t' << _dR[iter] << '\t' <<
                    d_DeltaR[iter] << '\n';
    VSFile.close();

    delete[] d_X;
}

void W_matrix::anotherCoordSys(double distance) {
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
        }
    }
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

