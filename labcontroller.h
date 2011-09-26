#ifndef LABCONTROLLER_H
#define LABCONTROLLER_H

#include <QThreadPool>
#include "lab.h"
#include "controller.h"
#include "mainwindow.h"

class LabController : public Controller
{
public:
    LabController();
    void set_matrixInputFileStr(QString);
    void set_impurityInputFileStr(QString);
    void notifyStartcalculating();
    void matrixFileChosen(QString);
    void impurityFileChosen(QString);
    void fillMatrixFields(Hop_integrals::paramFile*);
    void fillImpurFields(Hop_integrals::paramFile *param);
    void setMainWindow(void*);
    void HamCalcProgressUpdate(int);

    Lab lab;
    QThreadPool threadPool;
    QString matrixInputFileStr;
    QString impurityInputFileStr;
    MainWindow *window;
    static const int SIZE = 5;
};

#endif // LABCONTROLLER_H
