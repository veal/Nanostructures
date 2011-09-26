#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

//#include "sys_param.h"
//#include "lab.h"
#include "Hop_integrals.h"

class Hamiltonian {
public:
    Hamiltonian();
    void F_Hk(Doub kx, Comp H[N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3],
            Hop_integrals* hopIntegral, uint integralType);
    void Hop_int(int Nkp, Comp H_k[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3], Hop_integrals* hopIntegral);
    void Hij_k(Doub kx, Comp Hr[3][3][N_LAT][N_LAT][N_Band][N_Band],
            Comp Hij[N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]);
    void Calculate_Hamiltonian(Hop_integrals* hopIntegral);

};

#endif // HAMILTONIAN_H
