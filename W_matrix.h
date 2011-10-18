#ifndef W_MATRIX_H
#define	W_MATRIX_H

#include "PWF.h"

class W_matrix {

    class integralContainer {
        public:
            double* W_container;
            double* Wsh_container;
            double distance;
            integralContainer() {
                W_container = new double[14];
                Wsh_container = new double[14];
                distance = 0.0;
            }

            integralContainer(double dist) {
                W_container = new double[14];
                Wsh_container = new double[14];
                distance = dist;
            }

            ~integralContainer() {
            }
    };

    vector<integralContainer> integralList;

    static const int N_tot = 441;
    static const int N_Max = 441;
    static const double LATTICE_CONSTANT = 2.684310019; //honeycomb graphene structures Bohr radii

public:
    W_matrix(PWF*, PWF*);
    W_matrix(const W_matrix& orig);
    int get_integrals(double** W_container, double** Wsh_container, const double distance);
    virtual ~W_matrix();

    PWF* getImpurityPWF() const {
        return _impurityPWF;
    }

    PWF* getMatrixPWF() const {
        return _matrixPWF;
    }
    void F_Wk(Doub kx, Comp H[N_LAT][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]);
    void Wk(int Nkp, Comp Wr[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band],
            Comp S_k[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3], PWF* _matrixPWF, PWF* _impurityPWF);
    void Calculate_Wamiltonian(PWF* _matrixPWF, PWF* _impurityPWF);


private:
    
    string PARAM_FILE;

    PWF* _matrixPWF;
    PWF* _impurityPWF;

    double (*_dBasisWF)[N_tot];
    double* _dteta;
    double (*_dteta2)[N_Max];
    double (*_dR2)[N_tot];
    double* d_DeltaR;
    double* _dR;
    double _dstep_teta;
    int (*_mk_atom2)[N_Max];
    int m_temp;

    double FunWsh(double* pot,int n_, int m_, int i_, int j_);

    double FunW(double* pot, int n_, int m_, int i_, int j_);

    double PolLagr(int orb, double teta, double teta2);

    void setGrid(const int);

    void anotherCoordSys(double distance);

    void publishRezults(int, double*);

    int calculateIntegrals(integralContainer*, double distance);

    void setBasis(PWF* pwf);
};

#endif	/* W_MATRIX_H */
