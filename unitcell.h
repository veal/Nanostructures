#ifndef UNITCELL_H
#define UNITCELL_H

#include <math.h>
#include <iostream>
#include <QFile>
#include <QTextStream>

#define ind(x, y)   (x * 3 + y)

class UnitCell
{
public:
    UnitCell(int, int);
    void setAtomPositionsInCell(double* r);
    void setTranslationVectors(double* vectors);
    double getExpPositions(int cellNumber, int lattice, int component);
    double getCosPositions(int cellNumber, int lattice, int component);
    int numberOfCells;
    int zeroCell;
    ~UnitCell();
private:
    int N_LAT;
    int DIMENSION;
    double* r;
    double* hamVectors;
    double* translationVectors;
    double* inverseVectors;
    void obtainInverseVectors();
    void getAtomPositionsForHam();
};


#endif // UNITCELL_H
