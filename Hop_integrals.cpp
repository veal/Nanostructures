#include "Hop_integrals.h"
Hop_integrals::Hop_integrals() {
    _dBasisWF = new double[4][N_tot];
    _dteta = new double[N_Max];
    _dteta2 = new double[N_tot][N_Max];
    _dR2 = new double[N_tot][N_Max];
    _mk_atom2 = new int[N_tot][N_Max];
    d_DeltaR = new double[N_tot];
    _dR = new double[N_tot];
    PARAM_FILE = "/home/veal/Sandbox/GUINano/Rezults/Hop_intergral_params.txt";
    WF_FILE = "/home/veal/Sandbox/GUINano/Rezults/";
    _dstep_teta = Pi/(N_Max-1);
}

Hop_integrals::Hop_integrals(QString matrixFilename, QString impurityFilename) {

    _matrixFilename = matrixFilename;
    _impurityFilename = impurityFilename;

    _dBasisWF = new double[4][N_tot];
    _dteta = new double[N_Max];
    _dteta2 = new double[N_tot][N_Max];
    _dR2 = new double[N_tot][N_Max];
    _mk_atom2 = new int[N_tot][N_Max];
    d_DeltaR = new double[N_tot];
    _dR = new double[N_tot];
    PARAM_FILE = "/home/veal/Sandbox/GUINano/Rezults/Hop_intergral_params.txt";
    WF_FILE = "/home/veal/Sandbox/GUINano/Rezults/";
    _dstep_teta = Pi/(N_Max-1);

    _matrixParam = readParamFile(_matrixFilename);
//    assert(_matrixParam);
    _impurityParam = readParamFile(_impurityFilename);
//    assert(_impurityParam);
    _matrixPWF = readPWFFile(_matrixFilename);
    setBasis(_matrixPWF);
    _impurityPWF = readPWFFile(_impurityFilename);
    if (_matrixPWF == NULL || _impurityPWF == NULL) {
        cout << "Failed to read files with data for Hop_integrals\n";
//        assert(_matrixPWF == NULL);
//        assert(_impurityPWF == NULL);
    }
    setGrid(_matrixPWF->nuclearCharge);
}

int Hop_integrals::get_integrals(double** V_container, double** S_container, double** W_container,
        const double distance) {
    for (int iter = 0; iter < integralList.size(); iter++) {
        if (abs(integralList[iter].distance - distance) <= 0.01*distance) {
            *V_container = integralList[iter].Vcontainer;
            *S_container = integralList[iter].Scontainer;
            *W_container = integralList[iter].Wcontainer;
            return 0;
        }
    }

    integralContainer newIntegral(distance); ///!!!!!!!!!!!!
    calculateIntegrals(
            &newIntegral, distance);

//    *V_container = newIntegral.Vcontainer;
//    *S_container = newIntegral.Scontainer;
//    *W_container = newIntegral.Wcontainer;
    integralList.push_back(newIntegral);
    for (int iter = 0; iter < integralList.size(); iter++) {
        if (abs(integralList[iter].distance - distance) <= 0.01*distance) {
            *V_container = integralList[iter].Vcontainer;
            *S_container = integralList[iter].Scontainer;
            *W_container = integralList[iter].Wcontainer;
            return 0;
        }
    }

    return 0;
}
int Hop_integrals::calculateIntegrals(integralContainer *intCon, double distance) {

    double lclDistance = distance / 5.29e-11;    

    anotherCoordSys(lclDistance);

    QString str;
    QTextStream sstr(&str);
    sstr << WF_FILE << "V_S_container.dat";
    QFile VSFile(str);

    if (!VSFile.open(QFile::WriteOnly | QFile::Append)) {
        QString temp = "Can't find pwf file :";
        temp += "Rezults/V_S_container.dat";
        cout << temp.toStdString() << '\n';
        return NULL;
    }

    QTextStream in(&VSFile);
    in.setRealNumberNotation(QTextStream::ScientificNotation);

    in << "Computation result for " <<_matrixFilename << '\n';
    in << "Cut radius = " << _matrixParam->dr_max << "\tV\t" << "\tS\t"
            << "\tW\t" << "Distance = " << distance << '\n';
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
        double d_Int1 = 0, d_Int2 = 0, d_Int3 = 0;
        for (int i = 1; i < N_tot; i++) {
            for (int j = 1; j < N_Max; j++) {   
                d_Int1 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunS(n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i-1][j-1]) +
                FunS(n,m,i-1,j)*PolLagr(m_orb, _dteta[j], _dteta2[i-1][j]) +
                FunS(n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i][j-1]) +
                FunS(n,m,i,j)*PolLagr(m_orb,_dteta[j], _dteta2[i][j]));
                d_Int2 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunV(cleanPotential,n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i-1][j-1])+
                FunV(cleanPotential,n,m,i-1,j)*PolLagr(m_orb,_dteta[j],_dteta2[i-1][j]) +
                FunV(cleanPotential,n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i][j-1]) +
                FunV(cleanPotential,n,m,i,j)*PolLagr(m_orb,_dteta[j],_dteta2[i][j]));
                d_Int3 += (1.0/4.0)*_dstep_teta*d_DeltaR[i] *
                (FunW(impurityPotential,n,m,i-1,j-1)*PolLagr(m_orb,_dteta[j-1], _dteta2[i-1][j-1])+
                FunW(impurityPotential,n,m,i-1,j)*PolLagr(m_orb,_dteta[j],_dteta2[i-1][j]) +
                FunW(impurityPotential,n,m,i,j-1)*PolLagr(m_orb,_dteta[j-1],_dteta2[i][j-1]) +
                FunW(impurityPotential,n,m,i,j)*PolLagr(m_orb,_dteta[j],_dteta2[i][j]));
            }
        }
        intCon->Scontainer[m_orb] = d_Int1;
        intCon->Vcontainer[m_orb] = (1.0/2.0)*(_matrixPWF->getEnergyLevel()[n]+_matrixPWF->getEnergyLevel()[m])*
        d_Int1 + d_Int2;

        intCon->Wcontainer[m_orb] = d_Int3 - d_Int2;
        in << n << '\t' << m << '\t' << lm << '\t' << scientific << intCon->Vcontainer[m_orb] <<
                '\t' << intCon->Scontainer[m_orb] << '\t' << intCon->Wcontainer[m_orb] << '\n';   
    }

    VSFile.close();

    return 0;
}

double Hop_integrals::FunS(int n_, int m_, int i_, int j_) {
    double result = _dBasisWF[n_][i_] * _dBasisWF[m_][_mk_atom2[i_][j_]]
                    * _dR[i_] * _dR[i_] * sin(_dteta[j_]);
    if (result < 0.0)
        result = 0.0;
    return result;
}

double Hop_integrals::FunV(double* pot, int n_, int m_, int i_, int j_) {
    return (1.0 / 2.0)*_dBasisWF[n_][i_] * _dBasisWF[m_][_mk_atom2[i_][j_]]
            * _dR[i_] * _dR[i_] * sin(_dteta[j_])*(pot[i_] + pot[_mk_atom2[i_][j_]]);
}

double Hop_integrals::FunW(double* pot, int n_, int m_, int i_, int j_) {
    return _dBasisWF[n_][i_] * _dBasisWF[m_][_mk_atom2[i_][j_]]
            * _dR[i_] * _dR[i_] * sin(_dteta[j_])*(pot[i_] + pot[_mk_atom2[i_][j_]]);
}

double Hop_integrals::PolLagr(int orb, double teta, double teta2) {
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

double Hop_integrals::FunV1(int n_, int i_) {
    return (1.0 / 2.0)*_dBasisWF[n_][i_] * _dBasisWF[n_][i_]
            * _dR[i_] * _dR[i_] * _impurityPWF->getScaledPotential(_matrixPWF->nuclearCharge)[i_];
}

Hop_integrals::Hop_integrals(const Hop_integrals& orig) {
}

Hop_integrals::~Hop_integrals() {
    for (int iter = 0; iter < integralList.size(); iter++) {
        delete[] integralList[iter].Scontainer;
        delete[] integralList[iter].Vcontainer;
        delete[] integralList[iter].Wcontainer;
    }
    delete[] d_DeltaR;
    delete[] _dR;
    delete[] _dBasisWF;
    delete[] _dteta;
    delete[] _dteta2;
    delete[] _dR2;
    delete[] _mk_atom2;
}

void Hop_integrals::setGrid(int m_Z) {
    int m_Nblock = N_tot/40;
    int iter = 0;
    double d_C = 0.88534138/pow(m_Z, 1.0/3.0);
    double d_DeltaX = 0.0025;
    
    QFile file(PARAM_FILE);
    if (!file.open(QFile::WriteOnly | QFile::Append)) {
        cout << "couldn't open paramFile for writing.. Now exiting..\n";
        return;
    }
    QTextStream fileStream(&file);

    fileStream << "Nuclear charge = " << m_Z << '\n'
            << "Number of blocks = " << m_Nblock << '\n';

    fileStream << "\td_X\t" << "\td_R\t" << "\td_DeltaR" << '\n';
    double* d_X = new double[N_tot];

    for (int j = 0; j < m_Nblock; j++) {
        for (int k = 0; k < 40; k++) {
            iter++;
            d_X[iter] = d_X[iter-1] + d_DeltaX;
            _dR[iter-1] = d_X[iter-1]*d_C;
            d_DeltaR[iter] = d_C*d_DeltaX;
            fileStream << iter-1 << '\t' << scientific << d_X[iter-1] << '\t' << _dR[iter-1] << '\t' <<
                    d_DeltaR[iter] << '\n';
        }
        d_DeltaX += d_DeltaX;
    }

    _dR[iter] = d_X[iter]*d_C;
    fileStream << iter << '\t' << scientific << d_X[iter] << '\t' << _dR[iter] << '\t' <<
                    d_DeltaR[iter] << '\n';
    file.close();

    delete[] d_X;
}

Hop_integrals::paramFile* Hop_integrals::readParamFile(QString filename) {
    paramFile* newParam = new paramFile();
    newParam->initFile = filename;
    QFile file_desc(filename);
//    ifstream file_desc();
    if (!file_desc.open(QFile::ReadOnly)) {
        cout << "Can't find file for Hop_integrals" << '\n';
        return NULL;
    }

    QTextStream in(&file_desc);

    QString sz_temp;

    sz_temp = in.readLine();
    newParam->elementLabel = in.readLine();
    sz_temp = in.readLine();

    in >> newParam->m_Norbital;
    in >> newParam->m_BandNumber[0];
    in >> newParam->m_BandNumber[1];
    in >> newParam->m_BandNumber[2];
    in >> newParam->m_BandNumber[3];

    sz_temp = in.readLine();
    sz_temp = in.readLine();
    sz_temp = in.readLine();
    newParam->latticeConst = sz_temp.trimmed().toDouble();
    sz_temp = in.readLine();

    in >> newParam->dr_max;
    sz_temp = in.readLine();
    sz_temp = in.readLine();
    sz_temp = in.readLine();
    newParam->sz_InputFile = sz_temp;

//    newParam->sz_InputFile.resize(newParam->sz_InputFile.size() - 1);

    file_desc.close();

    QFile file(PARAM_FILE);
    if (!file.open(QFile::WriteOnly | QFile::Append)) {
        cout << "Couldn't open paramFile for writing.. Now exiting..\n";
        return NULL;
    }

    QTextStream fileStream(&file);
    fileStream << "Number of orbitals : " << newParam->m_Norbital << '\n' <<
            "First Band : " << newParam->m_BandNumber[0] << '\n' <<
            "Second Band : " << newParam->m_BandNumber[1] << '\n' <<
            "Third Band : " << newParam->m_BandNumber[2] << '\n' <<
            "Fourth Band : " << newParam->m_BandNumber[3] << '\n' <<
            "Cutting radius : " << newParam->dr_max << '\n' <<
            "Data file : " << newParam->sz_InputFile << '\n';

    file.close();
    paramFileList.push_back(*newParam);
    return newParam;
}

PWF* Hop_integrals::readPWFFile(QString initFile) {
    paramFile currParamFile;
    double* dPotential = new double[N_tot];
    double* dEnergyLevels = new double[4];
    int m_Z;
    double (*_dRadialPartWF)[N_tot] = new double[4][N_tot];
    for (int iter = 0; iter < paramFileList.size(); iter++) {
        if (!paramFileList[iter].initFile.compare(initFile)) {
            currParamFile = paramFileList[iter];
            break;
        }
    }
    QFile inputFile(currParamFile.sz_InputFile);

    if (!inputFile.open(QFile::ReadOnly)) {
        QString temp = "Can't find pwf file :";
        temp += currParamFile.sz_InputFile;
        cout << temp.toStdString() << '\n';
        return NULL;
    }

    QTextStream in(&inputFile);

    in >> m_Z >> m_temp >> m_temp >> m_temp >> m_temp;

    sz_temp = in.readLine();
    sz_temp = in.readLine();

    setGrid(m_Z);

    for (int iter = 0; iter < N_tot; iter++) {
        in >> dPotential[iter];
    }

    for (int j = 1; j < N_tot; j++) {
        dPotential[j] /= _dR[j];
    }
    dPotential[0] = dPotential[1];

    for (int i = 0; i < N_tot; i++) {
        _dRadialPartWF[0][i] = 0.0;
        _dRadialPartWF[1][i] = 0.0;
        _dRadialPartWF[2][i] = 0.0;
        _dRadialPartWF[3][i] = 0.0;
    }

    for (int iter1 = 0; iter1 < currParamFile.m_Norbital; iter1++) {
        int m_Orbital_number;
        double d_EnergyofState;
        double d_temp;
        int m_NumberofPoints;
        in >> m_Orbital_number >> d_temp >> d_EnergyofState >> d_temp >> m_NumberofPoints >> d_temp;
        sz_temp = in.readLine();
        if (m_Orbital_number == currParamFile.m_BandNumber[0]) {
            for (int iter = 0; iter < m_NumberofPoints; iter++) {
                in >> _dRadialPartWF[0][iter];
                if (iter != 0)
                    _dRadialPartWF[0][iter] /= _dR[iter];
            }
            _dRadialPartWF[0][0] = _dRadialPartWF[0][1];
            dEnergyLevels[0] = d_EnergyofState;
        } else
        if (m_Orbital_number == currParamFile.m_BandNumber[1]) {
                for (int iter = 0; iter < m_NumberofPoints; iter++) {
                    in >> _dRadialPartWF[1][iter];
                    if (iter != 0)
                        _dRadialPartWF[1][iter] /= _dR[iter];
                }
                _dRadialPartWF[1][0] = _dRadialPartWF[1][1];
                dEnergyLevels[1] = d_EnergyofState;
        } else
        if (m_Orbital_number == currParamFile.m_BandNumber[2]) {
                for (int iter = 0; iter < m_NumberofPoints; iter++) {
                    in >> _dRadialPartWF[2][iter];
                    if (iter != 0)
                        _dRadialPartWF[2][iter] /= _dR[iter];
                }
                _dRadialPartWF[2][0] = _dRadialPartWF[2][1];
                dEnergyLevels[2] = d_EnergyofState;
        } else
        if (m_Orbital_number == currParamFile.m_BandNumber[3]) {
                for (int iter = 0; iter < m_NumberofPoints; iter++) {
                    in >> _dRadialPartWF[3][iter];
                    if (iter != 0)
                        _dRadialPartWF[3][iter] /= _dR[iter];
                }
                _dRadialPartWF[3][0] = _dRadialPartWF[3][1];
                dEnergyLevels[3] = d_EnergyofState;
        }
        else {
            for (int iter = 0; iter < m_NumberofPoints; iter++) {
                in >> d_temp;
            }
            sz_temp = in.readLine();
        }

    }

    inputFile.close();

    for (int i = 0; i < N_tot; i++) {
        if ((currParamFile.dr_max*LATTICE_CONSTANT)<_dR[i]) {
            _dRadialPartWF[0][i] = 0.0;
            _dRadialPartWF[1][i] = 0.0;
            _dRadialPartWF[2][i] = 0.0;
            _dRadialPartWF[3][i] = 0.0;
        }
    }

    publishRezults(0, _dRadialPartWF[0]);
    publishRezults(1, _dRadialPartWF[1]);
    publishRezults(2, _dRadialPartWF[2]);
    publishRezults(3, _dRadialPartWF[3]);
    publishRezults(-1, dPotential);

    return new PWF(dPotential, _dRadialPartWF[0], _dRadialPartWF[1], _dRadialPartWF[2],
            _dRadialPartWF[3], N_tot, m_Z, dEnergyLevels);

}

void Hop_integrals::anotherCoordSys(double distance) {
//    ofstream VSFile("Rezults/V_S_container.dat", ios::app);
//    VSFile << "Another coord sys indices" << '\n';
    _dteta[0] = 0;
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

void Hop_integrals::publishRezults(int m_band, double* source) {
    QString temp_path;
    QTextStream sstr(&temp_path);
    if (m_band >= 0) {
        sstr << WF_FILE <<  _matrixParam->m_BandNumber[m_band] << "_WF_values_.dat";
    } else {
        sstr << WF_FILE << "potential_values_.dat";
    }
    QFile file(temp_path);

    if (!file.open(QFile::WriteOnly)) {
        QString temp = "Can't write to ";
        temp += temp_path;
        cout << temp.toStdString() <<  "file.\n";
        return;
    }

    QTextStream filestream(&file);

    filestream << "\n_______________________________\n";

//    file.precision(5);
    for (int i = 0; i < N_tot; i++) {
        filestream << i << '\t' << _dR[i] << '\t' << scientific << source[i] << '\n';
    }

    file.close();
}

void Hop_integrals::setBasis(PWF* pwf) {
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
