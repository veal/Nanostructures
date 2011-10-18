#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Hop_integrals.h"
#include "unitcell.h"

class Hamiltonian {
public:
    Hamiltonian();
    ~Hamiltonian();
    void F_Hk(Doub kx, Doub ky, Comp H[N_LAT][N_LAT][N_Band][N_Band], UnitCell* cell,
            Hop_integrals* hopIntegral, uint integralType);
    void Hop_int(int Nkp, UnitCell* cell, Hop_integrals* hopIntegral);
    void Hij_k(Doub kx, Comp Hr[3][3][N_LAT][N_LAT][N_Band][N_Band],
            Comp Hij[N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]);
    void Calculate_Hamiltonian(Hop_integrals* hopIntegral, UnitCell* cell);

    Comp (*H_k)[Nkp2][N_LAT][N_LAT][N_Band][N_Band];
    Comp (*S_k)[Nkp2][N_LAT][N_LAT][N_Band][N_Band];

};

#endif // HAMILTONIAN_H
