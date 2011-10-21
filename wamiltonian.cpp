#include "wamiltonian.h"

Wamiltonian::Wamiltonian()
{

}
void Wamiltonian::F_Wk(Doub kx, Doub ky, Comp H[N_LAT][N_LAT][N_LAT][N_Band][N_Band], UnitCell* cell, W_matrix* wMatrix) {
    Doub *V1, *V2, *V, *Vsh;
    Doub d, Pi, x, y, z;
    Doub *V01, *V02;
    wMatrix->get_integrals(&V01, &V02, 0);
    V = wMatrix->getMatrixPWF()->getEnergyLevel();

    Pi = 3.1415926536;

    for (int it = 0; it < N_LAT; it++) {
        for (int i = 0; i < N_LAT; i++) {
            for (int j = 0; j < N_LAT; j++) {

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

void Wamiltonian::Wk(int Nkp, UnitCell* cell, W_matrix* wMatrix) {
    for (int k1 = -Nkp; k1 < Nkp + 1; k1++) {
        for (int k2 = -Nkp; k2 < Nkp + 1; k2++) {
            Comp (*W_temp)[N_LAT][N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_LAT][N_Band][N_Band];
            cout << (k1 * 2 * Nkp + k2) * 100 / (4 * Nkp * Nkp) << '\n';

            Doub kx = static_cast<Doub> (k1) / (2.0 * Nkp);
            Doub ky = static_cast<Doub> (k2) / (2.0 * Nkp);
            F_Wk(kx, ky, W_temp, cell, wMatrix);


            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    for (int i = 0; i < N_LAT; i++) {
                        for (int ish = 0; ish < N_LAT; ish++) {
                            for (int iter = 1; iter < N_LAT; iter++) {
                                W_temp[0][i][ish][i1][j1] += W_temp[iter][i][ish][i1][j1];
                            }
                        }
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

void Wamiltonian::Calculate_Wamiltonian(W_matrix* wMatrix, UnitCell* cell) {
    if (calculateWamiltonian == true) {
        ofstream Ham_r_file;
        Ham_r_file.open("Rezults/W_k.dat", ios::binary);

        Wk(Nkp, cell, wMatrix);

        Ham_r_file.write((char*) W, N_Com * N_LAT * N_LAT * N_LAT * N_Band * N_Band * sizeof (Comp));
        Ham_r_file.close();
    } else {
        cout << "Reading Wamiltonian.." << '\n';
        ifstream Ham_r_file;
        Ham_r_file.open("Rezults/W_k.dat", ios::binary);
        Ham_r_file.read((char*) W, N_Com * N_LAT * N_LAT * N_LAT * N_Band * N_Band * sizeof (Comp));
        Ham_r_file.close();

    }
}
