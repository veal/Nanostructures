#ifndef LAB_H
#define LAB_H

#include <QString>
#include <QRunnable>
#include "unitcell.h"
#include "Hop_integrals.h"
#include "Hamiltonian.h"
#include "wamiltonian.h"
#include "W_matrix.h"
#include "controller.h"

class Lab : public QRunnable
{
public:
    Lab();
    void runCalculation();
    void set_matrixInputFile(QString);
    void set_impurityInputFile(QString);
    Hop_integrals* getHopIntegrals();
    void run();

    void setProgressUpdateCallback(Controller *c);
    void updateProgress(int);

private:
    Controller *controller;
    QString matrixInputFile;
    QString impurityInputFile;
    Hamiltonian *hamiltonian;
    Wamiltonian *wamiltonian;
    Hop_integrals hopIntegrals;
};

#endif // LAB_H
