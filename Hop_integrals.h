#include "PWF.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <complex>

#ifndef HOP_INTEGRALS_H
#define	HOP_INTEGRALS_H

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

    struct paramFile {
        string initFile;
        double dr_max;
        int m_BandNumber[4], m_Norbital;
        string sz_InputFile;
    };

    vector<integralContainer> integralList;
    vector<paramFile> paramFileList;

    static const int N_tot = 441;
    static const int N_Max = 441;
    static const double LATTICE_CONSTANT = 2.684310019; //honeycomb graphene structures Bohr radii

public:
    Hop_integrals(string, string);
    Hop_integrals(const Hop_integrals& orig);
    int get_integrals(double** V_container, double** S_container, double** W_container, const double distance);
    virtual ~Hop_integrals();

    PWF* getImpurityPWF() const {
        return _impurityPWF;
    }

    PWF* getMatrixPWF() const {
        return _matrixPWF;
    }


private:
    
    string _matrixFilename;
    string _impurityFilename;
    string sz_matrixInputFile;
    string sz_ImpurityInputFile;
    string sz_temp;
    string PARAM_FILE;
    string WF_FILE;

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
    
    
    double FunS(int n_, int m_, int i_, int j_);

    double FunV(double* pot,int n_, int m_, int i_, int j_);

    double FunW(double* pot,int n_, int m_, int i_, int j_);

    double FunV1(int n_, int i_);

    double PolLagr(int orb, double teta, double teta2);

    void setGrid(const int);

    paramFile* readParamFile(string fileName);

    PWF* readPWFFile(string initFile);

    int readImpurityPotentialFromPWFFile();

    void anotherCoordSys(double distance);

    void publishRezults(int, double*);

    int calculateIntegrals(integralContainer*, double distance);

    void setBasis(PWF* pwf);
};

#endif	/* HOP_INTEGRALS_H */