#include "sys_param.h"
#include "Hop_integrals.h"

void F_Hk(Doub kx, Comp H[N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3],
        Hop_integrals* hopIntegral, uint integralType) {
    if (VERBOSE)
        cout << "F_Hk started.." << '\n';
    Doub *V1, *V2, *V3;
    Doub *V;
    Doub V_S[] = {1.0, 1.0, 0.231147, 0.591972}; //!!! Change 3-rd and 4-th element on cut radius change
    V = hopIntegral->getMatrixPWF()->getEnergyLevel();
    Doub a, d, Pi, x, y, z;

    a = 3.0; //����� ������������ ������
    Pi = 3.1415926536;

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int n = 0; n < 3; n++) {
                x = (r[j][0][n] - r[i][0][1]);
                y = (r[j][1][n] - r[i][1][1]);
                z = (r[j][2][n] - r[i][2][1]);
                d = sqrt(x * x + y * y + z * z);
                if (d < 2.1) {
                    hopIntegral->get_integrals(&V1, &V2, &V3, 1.42e-10 * d);
                    if (VERBOSE)
                        cout << "Still working.." << '\n';
                    if (integralType == 1) {
                        V1 = V2;
                        V = V_S;
                    }
                }
                if (d < 2.1 && d > 0.0) {
                    H[i][j][0][0] += V1[0] * exp(CI * kx * x);
                    H[i][j][0][1] += (x / d) * V1[6] * exp(CI * kx * x);
                    H[i][j][0][2] += (y / d) * V1[6] * exp(CI * kx * x);
                    H[i][j][0][3] += (z / d) * V1[6] * exp(CI * kx * x);
                    H[i][j][0][4] += 1.732 * (x / d)*(y / d) * V1[7] * exp(CI * kx * x);
                    H[i][j][0][5] += 1.732 * (y / d)*(z / d) * V1[7] * exp(CI * kx * x);
                    H[i][j][0][6] += 1.732 * (z / d)*(x / d) * V1[7] * exp(CI * kx * x);
                    H[i][j][0][7] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[7] * exp(CI * kx * x);
                    H[i][j][0][8] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[7] * exp(CI * kx * x);
                    H[i][j][0][9] += V1[10] * exp(CI * kx * x);
                    H[i][j][1][0] -= (x / d) * V1[6] * exp(CI * kx * x);
                    H[i][j][1][1] += ((x / d)*(x / d) * V1[1] + (1 - (x / d)*(x / d)) * V1[2]) * exp(CI * kx * x);
                    H[i][j][1][2] += (x / d)*(y / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                    H[i][j][1][3] += (x / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                    H[i][j][1][4] += (1.732 * (x / d)*(x / d)*(y / d) * V1[8]+(y / d)*(1 - 2 * (x / d)*(x / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][1][5] += (x / d)*(y / d)*(z / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                    H[i][j][1][6] += (1.732 * (x / d)*(x / d)*(z / d) * V1[8]+(z / d)*(1 - 2 * (x / d)*(x / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][1][7] += (0.5 * 1.732 * (x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]+
                            (x / d)*(1 - (x / d)*(x / d)+ (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][1][8] += ((x / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8]
                            - 1.732 * (x / d)*(z / d)*(z / d) * V1[9]) * exp(CI * kx * x);
                    H[i][j][1][9] -= (x / d) * V1[11] * exp(CI * kx * x);
                    H[i][j][2][0] -= (y / d) * V1[6] * exp(CI * kx * x);
                    H[i][j][2][1] += (x / d)*(y / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                    H[i][j][2][2] += ((y / d)*(y / d) * V1[1] + (1 - (y / d)*(y / d)) * V1[2]) * exp(CI * kx * x);
                    H[i][j][2][3] += (y / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                    H[i][j][2][4] += (1.732 * (y / d)*(y / d)*(x / d) * V1[8]+(x / d)*(1 - 2 * (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][2][5] += (1.732 * (y / d)*(y / d)*(z / d) * V1[8]+(z / d)*(1 - 2 * (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][2][6] += (y / d)*(z / d)*(x / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                    H[i][j][2][7] += (0.5 * 1.732 * (y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]-
                            (y / d)*(1 + (x / d)*(x / d)-(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][2][8] += ((y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8] -
                            1.732 * (y / d)*(z / d)*(z / d) * V1[9]) * exp(CI * kx * x);
                    H[i][j][2][9] -= (y / d) * V1[11] * exp(CI * kx * x);
                    H[i][j][3][0] -= (z / d) * V1[6] * exp(CI * kx * x);
                    H[i][j][3][1] += (x / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                    H[i][j][3][2] += (y / d)*(z / d)*(V1[1] - V1[2]) * exp(CI * kx * x);
                    H[i][j][3][3] += ((z / d)*(z / d) * V1[1] + (1 - (z / d)*(z / d)) * V1[2]) * exp(CI * kx * x);
                    H[i][j][3][4] += (x / d)*(y / d)*(z / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                    H[i][j][3][5] += (1.732 * (z / d)*(z / d)*(y / d) * V1[8]+(y / d)*(1 - 2 * (z / d)*(z / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][3][6] += (1.732 * (z / d)*(z / d)*(x / d) * V1[8]+(x / d)*(1 - 2 * (z / d)*(z / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][3][7] += (0.5 * 1.732 * (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]-
                            (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][3][8] += ((z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8] +
                            1.732 * (z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][3][9] -= (z / d) * V1[11] * exp(CI * kx * x);
                    H[i][j][4][0] += 1.732 * (x / d)*(y / d) * V1[7] * exp(CI * kx * x);
                    H[i][j][4][1] -= (1.732 * (x / d)*(x / d)*(y / d) * V1[8]+(y / d)*(1 - 2.0 * (x / d)*(x / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][4][2] -= (1.732 * (y / d)*(y / d)*(x / d) * V1[8]+(x / d)*(1 - 2.0 * (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][4][3] -= (x / d)*(y / d)*(z / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                    H[i][j][4][4] += (3.0 * (x / d)*(x / d)*(y / d)*(y / d) * V1[3]+((x / d)*(x / d)+(y / d)*(y / d) - 4.0 * (x / d)*(x / d)*(y / d)*(y / d)) * V1[4]+
                            ((z / d)*(z / d)+(x / d)*(x / d)*(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][4][5] += (3.0 * (x / d)*(y / d)*(y / d)*(z / d) * V1[3]+(x / d)*(z / d)*(1 - 4 * (y / d)*(y / d)) * V1[4]+
                            (x / d)*(z / d)*((y / d)*(y / d) - 1) * V1[5]) * exp(CI * kx * x);
                    H[i][j][4][6] += (3 * (x / d)*(x / d)*(y / d)*(z / d) * V1[3]+(y / d)*(z / d)*(1 - 4.0 * (x / d)*(x / d)) * V1[4]+
                            (y / d)*(z / d)*((x / d)*(x / d) - 1) * V1[5]) * exp(CI * kx * x);
                    H[i][j][4][7] += (1.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3] + 2.0 * (x / d)*(y / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                            0.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][4][8] += (1.732 * (x / d)*(y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            2.0 * 1.732 * (x / d)*(y / d)*(z / d)*(z / d) * V1[4] + 0.5 * 1.732 * (x / d)*(y / d)*(1 + (z / d)*(z / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][4][9] += 1.732 * (x / d)*(y / d) * V1[12] * exp(CI * kx * x);
                    H[i][j][5][0] += 1.732 * (y / d)*(z / d) * V1[7] * exp(CI * kx * x);
                    H[i][j][5][1] -= (x / d)*(y / d)*(z / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                    H[i][j][5][2] -= (1.732 * (y / d)*(y / d)*(z / d) * V1[8]+(z / d)*(1 - 2 * (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][5][3] -= (1.732 * (z / d)*(z / d)*(y / d) * V1[8]+(y / d)*(1 - 2 * (z / d)*(z / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][5][4] += (3.0 * (x / d)*(y / d)*(y / d)*(z / d) * V1[3]+(x / d)*(z / d)*(1 - 4 * (y / d)*(y / d)) * V1[4]+
                            (x / d)*(z / d)*((y / d)*(y / d) - 1) * V1[5]) * exp(CI * kx * x);
                    H[i][j][5][5] += (3.0 * (y / d)*(y / d)*(z / d)*(z / d) * V1[3]+((y / d)*(y / d)+(z / d)*(z / d) - 4.0 * (y / d)*(y / d)*(z / d)*(z / d)) * V1[4]+
                            ((x / d)*(x / d)+(y / d)*(y / d)*(z / d)*(z / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][5][6] += (3.0 * (y / d)*(z / d)*(z / d)*(x / d) * V1[3]+(y / d)*(x / d)*(1 - 4 * (z / d)*(z / d)) * V1[4]+
                            (y / d)*(x / d)*((z / d)*(z / d) - 1) * V1[5]) * exp(CI * kx * x);
                    H[i][j][5][7] += (1.5 * (y / d)*(z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]-
                            (y / d)*(z / d)*(1 + 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]+
                            (y / d)*(z / d)*(1.0 + 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                    H[i][j][5][8] += (1.732 * (y / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                            0.5 * 1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][5][9] += (1.732 * (y / d)*(z / d) * V1[7]) * exp(CI * kx * x);
                    H[i][j][6][0] += 1.732 * (z / d)*(x / d) * V1[7] * exp(CI * kx * x);
                    H[i][j][6][1] -= (1.732 * (x / d)*(x / d)*(z / d) * V1[8]+(z / d)*(1 - 2 * (x / d)*(x / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][6][2] -= (y / d)*(z / d)*(x / d)*(1.732 * V1[8] - 2.0 * V1[9]) * exp(CI * kx * x);
                    H[i][j][6][3] -= (1.732 * (z / d)*(z / d)*(x / d) * V1[8]+(x / d)*(1 - 2 * (z / d)*(z / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][6][4] += (3 * (x / d)*(x / d)*(y / d)*(z / d) * V1[3]+(y / d)*(z / d)*(1 - 4.0 * (x / d)*(x / d)) * V1[4]+
                            (y / d)*(z / d)*((x / d)*(x / d) - 1) * V1[5]) * exp(CI * kx * x);
                    H[i][j][6][5] += (3.0 * (y / d)*(z / d)*(z / d)*(x / d) * V1[3]+(y / d)*(x / d)*(1 - 4 * (z / d)*(z / d)) * V1[4]+
                            (y / d)*(x / d)*((z / d)*(z / d) - 1) * V1[5]) * exp(CI * kx * x);
                    H[i][j][6][6] += (3.0 * (z / d)*(z / d)*(x / d)*(x / d) * V1[3]+((z / d)*(z / d)+(x / d)*(x / d) - 4.0 * (z / d)*(z / d)*(x / d)*(x / d)) * V1[4]+
                            ((y / d)*(y / d)+(z / d)*(z / d)*(x / d)*(x / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][6][7] += (1.5 * (z / d)*(x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]+
                            (z / d)*(x / d)*(1.0 - 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]-
                            (z / d)*(x / d)*(1.0 - 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                    H[i][j][6][8] += (1.732 * (x / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                            0.5 * 1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][6][9] += 1.732 * (z / d)*(x / d) * V1[12] * exp(CI * kx * x);
                    H[i][j][7][0] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[7] * exp(CI * kx * x);
                    H[i][j][7][1] -= (0.5 * 1.732 * (x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]+
                            (x / d)*(1 - (x / d)*(x / d) + (y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][7][2] -= (0.5 * 1.732 * (y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]-
                            (y / d)*(1 + (x / d)*(x / d)-(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][7][3] -= (0.5 * 1.732 * (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[8]-
                            (z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][7][4] += (1.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3] + 2.0 * (x / d)*(y / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                            0.5 * (x / d)*(y / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][7][5] += (1.5 * (y / d)*(z / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]-
                            (y / d)*(z / d)*(1 + 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]+
                            (y / d)*(z / d)*(1.0 + 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                    H[i][j][7][6] += (1.5 * (z / d)*(x / d)*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]+
                            (z / d)*(x / d)*(1.0 - 2.0 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]-
                            (z / d)*(x / d)*(1.0 - 0.5 * ((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                    H[i][j][7][7] += (0.75 * ((x / d)*(x / d)-(y / d)*(y / d))*((x / d)*(x / d)-(y / d)*(y / d)) * V1[3]+
                            ((x / d)*(x / d)+(y / d)*(y / d)-((x / d)*(x / d)-(y / d)*(y / d))*((x / d)*(x / d)-(y / d)*(y / d))) * V1[4]+
                            ((z / d)*(z / d) + 0.25 * ((x / d)*(x / d)-(y / d)*(y / d))*((x / d)*(x / d)-(y / d)*(y / d))) * V1[5]) * exp(CI * kx * x);
                    H[i][j][7][8] += (0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d))*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            1.732 * (z / d)*(z / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                            0.25 * 1.732 * (1 + (z / d)*(z / d))*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][7][9] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[12] * exp(CI * kx * x);
                    H[i][j][8][0] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[7] * exp(CI * kx * x);
                    H[i][j][8][1] -= ((x / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8]
                            - 1.732 * (x / d)*(z / d)*(z / d) * V1[9]) * exp(CI * kx * x);
                    H[i][j][8][2] -= ((y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8] -
                            1.732 * (y / d)*(z / d)*(z / d) * V1[9]) * exp(CI * kx * x);
                    H[i][j][8][3] -= ((z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[8] +
                            1.732 * (z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[9]) * exp(CI * kx * x);
                    H[i][j][8][4] += (1.732 * (x / d)*(y / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            2.0 * 1.732 * (x / d)*(y / d)*(z / d)*(z / d) * V1[4] + 0.5 * 1.732 * (x / d)*(y / d)*(1 + (z / d)*(z / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][8][5] += (1.732 * (y / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                            0.5 * 1.732 * (y / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][8][6] += (1.732 * (x / d)*(z / d)*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)-(z / d)*(z / d)) * V1[4] -
                            0.5 * 1.732 * (x / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][8][7] += (0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d))*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            1.732 * (z / d)*(z / d)*((y / d)*(y / d)-(x / d)*(x / d)) * V1[4] +
                            0.25 * 1.732 * (1 + (z / d)*(z / d))*((x / d)*(x / d)-(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][8][8] += (((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d)))*((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[3] +
                            3.0 * (z / d)*(z / d)*((x / d)*(x / d)+(y / d)*(y / d)) * V1[4] +
                            0.75 * ((x / d)*(x / d)+(y / d)*(y / d))*((x / d)*(x / d)+(y / d)*(y / d)) * V1[5]) * exp(CI * kx * x);
                    H[i][j][8][9] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[12] * exp(CI * kx * x);
                    H[i][j][9][0] += V1[10] * exp(CI * kx * x);
                    H[i][j][9][1] += (x / d) * V1[11] * exp(CI * kx * x);
                    H[i][j][9][2] += (y / d) * V1[11] * exp(CI * kx * x);
                    H[i][j][9][3] += (z / d) * V1[11] * exp(CI * kx * x);
                    H[i][j][9][4] += 1.732 * (x / d)*(y / d) * V1[12] * exp(CI * kx * x);
                    H[i][j][9][5] += (1.732 * (y / d)*(z / d) * V1[7]) * exp(CI * kx * x);
                    H[i][j][9][6] += 1.732 * (z / d)*(x / d) * V1[12] * exp(CI * kx * x);
                    H[i][j][9][7] += 0.5 * 1.732 * ((x / d)*(x / d)-(y / d)*(y / d)) * V1[12] * exp(CI * kx * x);
                    H[i][j][9][8] += ((z / d)*(z / d) - 0.5 * ((x / d)*(x / d)+(y / d)*(y / d))) * V1[12] * exp(CI * kx * x);
                    H[i][j][9][9] += V1[13] * exp(CI * kx * x);
                }
                if (d == 0) {
                    H[i][j][0][0] = V[0];
                    H[i][j][0][1] = 0.0;
                    H[i][j][0][2] = 0.0;
                    H[i][j][0][3] = 0.0;
                    H[i][j][0][4] = 0.0;
                    H[i][j][0][5] = 0.0;
                    H[i][j][0][6] = 0.0;
                    H[i][j][0][7] = 0.0;
                    H[i][j][0][8] = 0.0;
                    H[i][j][0][9] = 0.0;
                    H[i][j][1][0] = 0.0;
                    H[i][j][1][1] = V[1];
                    H[i][j][1][2] = 0.0;
                    H[i][j][1][3] = 0.0;
                    H[i][j][1][4] = 0.0;
                    H[i][j][1][5] = 0.0;
                    H[i][j][1][6] = 0.0;
                    H[i][j][1][7] = 0.0;
                    H[i][j][1][8] = 0.0;
                    H[i][j][1][9] = 0.0;
                    H[i][j][2][0] = 0.0;
                    H[i][j][2][1] = 0.0;
                    H[i][j][2][2] = V[1];
                    H[i][j][2][3] = 0.0;
                    H[i][j][2][4] = 0.0;
                    H[i][j][2][5] = 0.0;
                    H[i][j][2][6] = 0.0;
                    H[i][j][2][7] = 0.0;
                    H[i][j][2][8] = 0.0;
                    H[i][j][2][9] = 0.0;
                    H[i][j][3][0] = 0.0;
                    H[i][j][3][1] = 0.0;
                    H[i][j][3][2] = 0.0;
                    H[i][j][3][3] = V[1];
                    H[i][j][3][4] = 0.0;
                    H[i][j][3][5] = 0.0;
                    H[i][j][3][6] = 0.0;
                    H[i][j][3][7] = 0.0;
                    H[i][j][3][8] = 0.0;
                    H[i][j][3][9] = 0.0;
                    H[i][j][4][0] = 0.0;
                    H[i][j][4][1] = 0.0;
                    H[i][j][4][2] = 0.0;
                    H[i][j][4][3] = 0.0;
                    H[i][j][4][4] = V[2];
                    H[i][j][4][5] = 0.0;
                    H[i][j][4][6] = 0.0;
                    H[i][j][4][7] = 0.0;
                    H[i][j][4][8] = 0.0;
                    H[i][j][4][9] = 0.0;
                    H[i][j][5][0] = 0.0;
                    H[i][j][5][1] = 0.0;
                    H[i][j][5][2] = 0.0;
                    H[i][j][5][3] = 0.0;
                    H[i][j][5][4] = 0.0;
                    H[i][j][5][5] = V[2];
                    H[i][j][5][6] = 0.0;
                    H[i][j][5][7] = 0.0;
                    H[i][j][5][8] = 0.0;
                    H[i][j][5][9] = 0.0;
                    H[i][j][6][0] = 0.0;
                    H[i][j][6][1] = 0.0;
                    H[i][j][6][2] = 0.0;
                    H[i][j][6][3] = 0.0;
                    H[i][j][6][4] = 0.0;
                    H[i][j][6][5] = 0.0;
                    H[i][j][6][6] = V[2];
                    H[i][j][6][7] = 0.0;
                    H[i][j][6][8] = 0.0;
                    H[i][j][6][9] = 0.0;
                    H[i][j][7][0] = 0.0;
                    H[i][j][7][1] = 0.0;
                    H[i][j][7][2] = 0.0;
                    H[i][j][7][3] = 0.0;
                    H[i][j][7][4] = 0.0;
                    H[i][j][7][5] = 0.0;
                    H[i][j][7][6] = 0.0;
                    H[i][j][7][7] = V[2];
                    H[i][j][7][8] = 0.0;
                    H[i][j][7][9] = 0.0;
                    H[i][j][8][0] = 0.0;
                    H[i][j][8][1] = 0.0;
                    H[i][j][8][2] = 0.0;
                    H[i][j][8][3] = 0.0;
                    H[i][j][8][4] = 0.0;
                    H[i][j][8][5] = 0.0;
                    H[i][j][8][6] = 0.0;
                    H[i][j][8][7] = 0.0;
                    H[i][j][8][8] = V[2];
                    H[i][j][8][9] = 0.0;
                    H[i][j][9][0] = 0.0;
                    H[i][j][9][1] = 0.0;
                    H[i][j][9][2] = 0.0;
                    H[i][j][9][3] = 0.0;
                    H[i][j][9][4] = 0.0;
                    H[i][j][9][5] = 0.0;
                    H[i][j][9][6] = 0.0;
                    H[i][j][9][7] = 0.0;
                    H[i][j][9][8] = 0.0;
                    H[i][j][9][9] = V[3];
                }
            }
        }
    }
    if (VERBOSE)
        cout << "Hop_int finished.." << '\n';
}

void Hop_int(int Nkp, Comp H_k[Nkp2][N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3], Hop_integrals* hopIntegral) {
    if (VERBOSE)
        cout << "Hop_int started.." << '\n';
    ofstream Ham_k_file("Rezults/Ham.dat");
    ofstream HDif_k_file("Rezults/HDif.dat");

    for (int k = -Nkp; k < Nkp; k++) { //!!!!
        complex<Doub > (*H)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
        complex<Doub > (*W)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
        complex<Doub > (*S)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
        complex<Doub > (*H_Dif)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
        cout << k << '\n';

        Doub kx = 2 * Pi * static_cast<Doub> (k) / (3 * Nkp);
        F_Hk(kx, H, r, hopIntegral, 0);
        F_Hk(kx, S, r, hopIntegral, 1);

//		Speed(kx, V_p, H_Dif, r);
//		Trans_Hk(S, H_Dif, H_Dif);
        Trans_Hk(S, H, H, false);

//        Eigen_values(H, NULL, true);

        Doub ro_loc = ((k == -Nkp) || (k == Nkp - 1)) ? 0.5 : 1;

        for (int i = 0; i < N_LAT; i++) {
            for (int j = 0; j < N_LAT; j++) {
                for (int i1 = 0; i1 < N_Band; i1++) {
                    for (int j1 = 0; j1 < N_Band; j1++) {
                        H_k[k + Nkp][i][j][i1][j1] = H[i][j][i1][j1];
                        S_k[k + Nkp][i][j][i1][j1] = S[i][j][i1][j1];
//                        if (abs(imag(H[i][j][i1][j1]) + imag(H[j][i][j1][i1])) > 0.01*abs(imag(H[i][j][i1][j1]))
//                                && abs(real(H[i][j][i1][j1]) - real(H[j][i][j1][i1])) > 0.01*abs(real(H[i][j][i1][j1])))
//                        {
//                            cout << i << "  " << j << "  " << i1 << "  " << j1  << '\n';
//                            cout << real(H[i][j][i1][j1]) << "  " << imag(H[i][j][i1][j1]) << '\n';
//                            cout << real(H[j][i][j1][i1]) << "  " << imag(H[j][i][j1][i1]) << '\n';
//                        }
                        Ham_k_file << real(H[i][j][i1][j1]) << "       " << imag(H[i][j][i1][j1]) << '\n';
                    }
                }
            }
        }
        delete [] W;
        delete [] H;
        delete [] S;
        delete [] H_Dif;
    }

    Ham_k_file.close();
    HDif_k_file.close();
    if (VERBOSE)
        cout << "Hop_int finished.." << '\n';
}

void Hij_k(Doub kx, Comp Hr[3][3][N_LAT][N_LAT][N_Band][N_Band],
        Comp Hij[N_LAT][N_LAT][N_Band][N_Band], Doub r[N_LAT][3][3]) {
    Comp c1(0.0, 1.0), z;
    for (int i = 0; i < N_LAT; i++) {
        for (int i1 = 0; i1 < N_LAT; i1++) {
            for (int j = 0; j < N_Band; j++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    Hij[i][i1][j][j1] = 0;
                    for (int n = 0; n < 3; n++) {
                        for (int n1 = 0; n1 < 3; n1++) {
                            z = c1 * kx * (r[i][0][n] - r[i1][0][n1]);
                            Hij[i][i1][j][j1] = Hij[i][i1][j][j1] + Hr[n][n1][i][i1][j][j1] * exp(z);
                        }
                    }
                }
            }
        }
    }
}

void Calculate_Hamiltonian(Hop_integrals* hopIntegral) {
    if (VERBOSE)
        cout << "Calculate_hamiltonian started.." << '\n';
    if (calculateHamiltonian == true) {
        ofstream Ham_r_file;
        Ham_r_file.open("Rezults/H_k.dat", ios::binary);
        Hop_int(Nkp, H_k, r, hopIntegral);
        Ham_r_file.write((char*) H_k, Nkp2 * N_LAT * N_LAT * N_Band * N_Band * sizeof (Comp));
        Ham_r_file.close();
    } else {
        cout << "Reading Hamiltonian.." << '\n';
        ifstream Ham_r_file;
        Ham_r_file.open("Rezults/H_k.dat", ios::binary);
        Ham_r_file.read((char*) H_k, Nkp2 * N_LAT * N_LAT * N_Band * N_Band * sizeof (Comp));
        Ham_r_file.close();
    }
    if (VERBOSE)
        cout << "Calculate_hamiltonian finished.." << '\n';
}