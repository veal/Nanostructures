#ifndef HOP_INTEGRALS_H
#define	HOP_INTEGRALS_H

#include <QFile>
#include <QTextStream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <complex>
#include "PWF.h"

class Hop_integrals {

    class integralContainer {
        public:
            double* Vcontainer;
            double* Scontainer;
            double* Wcontainer;
            double distance;
            integralContainer() {
                Vcontainer = new double[14];
                Scontainer = new double[14];
                Wcontainer = new double[14];
                distance = 0.0;
            }
            integralContainer(double dist) {
                Vcontainer = new double[14];
                Scontainer = new double[14];
                Wcontainer = new double[14];
                distance = dist;
            }
            ~integralContainer() {
//                delete[] Vcontainer;
//                delete[] Scontainer;
//                delete[] Wcontainer;
            }
    };

    static const int N_tot = 441;
    static const int N_Max = 441;
    static const double LATTICE_CONSTANT = 2.684310019; //honeycomb graphene structures Bohr radii

public:

    struct paramFile {
        QString initFile;
        double dr_max;
        int m_BandNumber[4], m_Norbital;
        QString sz_InputFile;
        QString elementLabel;
        double latticeConst;
    };

    Hop_integrals();
    Hop_integrals(QString, QString);
    Hop_integrals(const Hop_integrals& orig);
    int get_integrals(double** V_container, double** S_container, double** W_container, const double distance);
    virtual ~Hop_integrals();

    PWF* getImpurityPWF() const {
        return _impurityPWF;
    }

    PWF* getMatrixPWF() const {
        return _matrixPWF;
    }

    paramFile* readParamFile(QString fileName);

private:

    QString _matrixFilename;
    QString _impurityFilename;
    QString sz_matrixInputFile;
    QString sz_ImpurityInputFile;
    QString sz_temp;
    QString PARAM_FILE;
    QString WF_FILE;

    PWF* _matrixPWF;
    PWF* _impurityPWF;
    paramFile* _matrixParam;
    paramFile* _impurityParam;
    double (*_dBasisWF)[N_tot];
    double* _dteta;
    double (*_dteta2)[N_Max];
    double (*_dR2)[N_tot];
    double* d_DeltaR;
    double* _dR;
    double _dstep_teta;
    int (*_mk_atom2)[N_Max];
    int m_temp;
    
    vector<integralContainer> integralList;
    vector<paramFile> paramFileList;

    double FunS(int n_, int m_, int i_, int j_);

    double FunV(double* pot,int n_, int m_, int i_, int j_);

    double FunW(double* pot,int n_, int m_, int i_, int j_);

    double FunV1(int n_, int i_);

    double PolLagr(int orb, double teta, double teta2);

    void setGrid(const int);

    PWF* readPWFFile(QString initFile);

    int readImpurityPotentialFromPWFFile();

    void anotherCoordSys(double distance);

    void publishRezults(int, double*);

    int calculateIntegrals(integralContainer*, double distance);

    void setBasis(PWF* pwf);
};

#endif	/* HOP_INTEGRALS_H */
