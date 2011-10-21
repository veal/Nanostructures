#ifndef WAMILTONIAN_H
#define WAMILTONIAN_H

#include "W_matrix.h"
#include "unitcell.h"

class Wamiltonian
{
public:
    Wamiltonian();
    void Calculate_Wamiltonian(W_matrix* wMatrix, UnitCell* cell);

private:
    void F_Wk(Doub kx, Doub ky, Comp H[N_LAT][N_LAT][N_LAT][N_Band][N_Band], UnitCell* cell, W_matrix* wMatrix);
    void Wk(int Nkp, UnitCell* cell, W_matrix* wMatrix);

    Comp (*Wr)[N_Com][N_LAT][N_LAT][N_LAT][N_Band][N_Band];
    Comp (*S_k)[Nkp2][N_LAT][N_LAT][N_Band][N_Band];
};

#endif // WAMILTONIAN_H
