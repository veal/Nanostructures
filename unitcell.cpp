#include "unitcell.h"

UnitCell::UnitCell(int nlat, int dimension)
{
    N_LAT = nlat;
    DIMENSION = dimension;
    numberOfCells = pow(3, DIMENSION);
    r = new double[N_LAT * numberOfCells * 3];
    translationVectors = new double[DIMENSION * 3];
    inverseVectors = new double[DIMENSION * 3];
    hamVectors = new double[N_LAT * 3 * numberOfCells];
    switch (DIMENSION) {
    case 1:
        zeroCell = 1;
        break;
    case 2:
        zeroCell = 4;
        break;
    case 3:
        zeroCell = 13;
    }
}
UnitCell::~UnitCell() {
    delete r;
    delete translationVectors;
    delete inverseVectors;
    delete hamVectors;
}

void UnitCell::setAtomPositionsInCell(double* r) {
    if (DIMENSION == 1) {
        int num = 0;
        for (int x = -1; x < 2; x++) {
            for (int i = 0; i < N_LAT; i++) {
                for (int j = 0; j < 3; j++) {
                    this->r[(i * 3 + j) * numberOfCells + num] = r[i * 3 + j] + x * translationVectors[j];
                }
            }
            num++;
        }
    }
    if (DIMENSION == 2) {
        int num = 0;
        for (int x = -1; x < 2; x++) {
            for (int y = -1; y < 2; y++) {
                for (int i = 0; i < N_LAT; i++) {
                    for (int j = 0; j < 3; j++) {
                        this->r[(i * 3 + j) * numberOfCells + num] = r[i * 3 + j] + x * translationVectors[j] + y * translationVectors[3 + j];
                    }
                }
                num++;
            }
        }
    }
    if (DIMENSION == 3) {
        int num = 0;
        for (int x = -1; x < 2; x++) {
            for (int y = -1; y < 2; y++) {
                for (int z = -1; z < 2; z++) {
                    for (int i = 0; i < N_LAT; i++) {
                        for (int j = 0; j < 3; j++) {
                            this->r[(i * 3 + j) * numberOfCells + num] = r[i * 3 + j] + x * translationVectors[j] +
                                    y * translationVectors[3 + j] + z * translationVectors[6 + j];
                        }
                    }
                    num++;
                }
            }
        }
    }
    getAtomPositionsForHam();
}

void UnitCell::setTranslationVectors(double* vectors) {
    for (int i = 0; i < DIMENSION; i++) {
        for (int j = 0; j < 3; j++) {
            translationVectors[i * 3 + j] = vectors[i * 3 + j];
        }
    }
}

void UnitCell::getAtomPositionsForHam() {
    obtainInverseVectors();
    QString str;
    QTextStream sstr(&str);
    sstr <<  "/home/veal/Sandbox/GUINano/Rezults/Lattice_container.dat";
    QFile VSFile(str);

    if (!VSFile.open(QFile::WriteOnly | QFile::Append)) {
        QString temp = "Can't open file :";
        temp += "Rezults/Lattice_container.dat. Exiting..";
        std::cout << temp.toStdString() << '\n';
        return;
    }

    QTextStream in(&VSFile);
    in.setRealNumberNotation(QTextStream::ScientificNotation);

    in << "Computation result for lattice sites:\n";
    for (int w = 0; w < numberOfCells; w++) {                   //loop for iterating neighbouring unit cells
        std::cout << "cell = " << w << "; \n";
        for (int k = 0; k < N_LAT; k++) {                       //loop for iterating all hamvectors
            std::cout << "nlat = " << k << "; \n";
            for (int j = 0; j < DIMENSION; j++) {               //loop for iterating all components of a single hamvector
                std::cout << "j = " << j << "; \n";
                hamVectors[(k * 3 + j) * numberOfCells + w] = 0.0;
                for (int i = 0; i < DIMENSION; i++) {
                    hamVectors[(k * 3 + j) * numberOfCells + w] += inverseVectors[j * 3 + i] * r[(k * 3 + i) * numberOfCells + w];
                }
                double a = hamVectors[(k * 3 + j) * numberOfCells + w];
                double absVal = a > 0.0 ? a : -a;
                if (absVal < 0.1)
                    hamVectors[(k * 3 + j) * numberOfCells + w] = 0.0;
            }
            in << "origVector = (" << r[k * 3 * numberOfCells + w] << ", " <<
                  r[(k * 3 + 1) * numberOfCells + w] << ", " <<  r[(k * 3 + 2) * numberOfCells + w] << "); \n";
            in << "hamVector = (" << hamVectors[k * 3 * numberOfCells + w] << ", " <<
                  hamVectors[(k * 3 + 1) * numberOfCells + w] << ", " << hamVectors[(k * 3 + 2) * numberOfCells + w] << "); \n";
            in << "\n";
        }
    }
    VSFile.close();
}

void UnitCell::obtainInverseVectors() {
    if (DIMENSION == 2) {
        inverseVectors[ind(0, 0)] = translationVectors[ind(1, 1)];
        inverseVectors[ind(0, 1)] = -translationVectors[ind(1, 0)];

        inverseVectors[ind(1, 0)] = translationVectors[ind(0, 1)];
        inverseVectors[ind(1, 1)] = -translationVectors[ind(0, 0)];

        double sign = inverseVectors[ind(0, 0)] * translationVectors[ind(0, 0)] + inverseVectors[ind(0, 1)] * translationVectors[ind(0, 1)];
        double norm = 2 * 3.14 / sign;
        inverseVectors[ind(0, 0)] *= norm;
        inverseVectors[ind(0, 1)] *= norm;

//        std::cout << "j = " << 0 << "; kx = " << inverseVectors[ind(0, 0)] << ";\n";
//        std::cout << "j = " << 0 << "; ky = " << inverseVectors[ind(0, 1)] << ";\n";

        sign = inverseVectors[ind(1, 0)] * translationVectors[ind(1, 0)] + inverseVectors[ind(1, 1)] * translationVectors[ind(1, 1)];
        norm = 2 * 3.14 / sign;
        inverseVectors[ind(1, 0)] *= norm;
        inverseVectors[ind(1, 1)] *= norm;

//        std::cout << "j = " << 1 << "; kx = " << inverseVectors[ind(1, 0)] << ";\n";
//        std::cout << "j = " << 1 << "; ky = " << inverseVectors[ind(1, 1)] << ";\n";
    }
    if (DIMENSION == 3) {
        inverseVectors[ind(0, 0)] = translationVectors[ind(1, 1)] * translationVectors[ind(2, 2)] - translationVectors[ind(1, 2)] * translationVectors[ind(2, 1)];
        inverseVectors[ind(0, 1)] = translationVectors[ind(1, 2)] * translationVectors[ind(2, 0)] - translationVectors[ind(1, 0)] * translationVectors[ind(2, 2)];
        inverseVectors[ind(0, 2)] = translationVectors[ind(1, 0)] * translationVectors[ind(2, 1)] - translationVectors[ind(1, 1)] * translationVectors[ind(2, 0)];

        inverseVectors[ind(1, 0)] = translationVectors[ind(0, 1)] * translationVectors[ind(2, 2)] - translationVectors[ind(0, 2)] * translationVectors[ind(2, 1)];
        inverseVectors[ind(1, 1)] = translationVectors[ind(0, 2)] * translationVectors[ind(2, 0)] - translationVectors[ind(0, 0)] * translationVectors[ind(2, 2)];
        inverseVectors[ind(1, 2)] = translationVectors[ind(0, 0)] * translationVectors[ind(2, 1)] - translationVectors[ind(0, 1)] * translationVectors[ind(2, 0)];

        inverseVectors[ind(2, 0)] = translationVectors[ind(0, 1)] * translationVectors[ind(1, 2)] - translationVectors[ind(0, 2)] * translationVectors[ind(1, 1)];
        inverseVectors[ind(2, 1)] = translationVectors[ind(0, 2)] * translationVectors[ind(1, 0)] - translationVectors[ind(0, 0)] * translationVectors[ind(1, 2)];
        inverseVectors[ind(2, 2)] = translationVectors[ind(0, 0)] * translationVectors[ind(1, 1)] - translationVectors[ind(0, 1)] * translationVectors[ind(1, 0)];

        double sign = inverseVectors[ind(0, 0)] * translationVectors[ind(0, 0)] + inverseVectors[ind(0, 1)] * translationVectors[ind(0, 1)] + inverseVectors[ind(0, 2)] * inverseVectors[ind(0, 2)];
        double norm = 2 * 3.14 / sqrt(inverseVectors[ind(0, 0)] * inverseVectors[ind(0, 0)] + inverseVectors[ind(0, 1)] * inverseVectors[ind(0, 1)] + inverseVectors[ind(0, 2)] * inverseVectors[ind(0, 2)]);
        inverseVectors[ind(0, 0)] *= norm;
        inverseVectors[ind(0, 1)] *= norm;
        inverseVectors[ind(0, 2)] *= norm;

        sign = inverseVectors[ind(1, 0)] * translationVectors[ind(1, 0)] + inverseVectors[ind(1, 1)] * translationVectors[ind(1, 1)] + inverseVectors[ind(1, 2)] * inverseVectors[ind(1, 2)];
        norm = 2 * 3.14 / sqrt(inverseVectors[ind(1, 0)] * inverseVectors[ind(1, 0)] + inverseVectors[ind(1, 1)] * inverseVectors[ind(1, 1)] + inverseVectors[ind(1, 2)] * inverseVectors[ind(1, 2)]);
        inverseVectors[ind(1, 0)] *= norm;
        inverseVectors[ind(1, 1)] *= norm;
        inverseVectors[ind(1, 2)] *= norm;

        sign = inverseVectors[ind(2, 0)] * translationVectors[ind(2, 0)] + inverseVectors[ind(2, 1)] * translationVectors[ind(2, 1)] + inverseVectors[ind(2, 2)] * inverseVectors[ind(2, 2)];
        norm = 2 * 3.14 / sqrt(inverseVectors[ind(2, 0)] * inverseVectors[ind(2, 0)] + inverseVectors[ind(2, 1)] * inverseVectors[ind(2, 1)] + inverseVectors[ind(2, 2)] * inverseVectors[ind(2, 2)]);
        inverseVectors[ind(2, 0)] *= norm;
        inverseVectors[ind(2, 1)] *= norm;
        inverseVectors[ind(2, 2)] *= norm;
    }
}
double UnitCell::getExpPositions(int cellNumber, int lattice, int component) {
    return hamVectors[(lattice * 3 + component) * numberOfCells + cellNumber];
}
double UnitCell::getCosPositions(int cellNumber, int lattice, int component) {
    return r[(lattice * 3 + component) * numberOfCells + cellNumber];
}
