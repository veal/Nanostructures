#ifndef LABCONTROLLER_H
#define LABCONTROLLER_H

#include "controller.h"
#include <QThreadPool>
#include "lab.h"

class LabController : public Controller
{
public:
    LabController();
    void set_matrixInputFileStr(QString);
    void set_impurityInputFileStr(QString);
    void notifyStartcalculating();

    Lab lab;
    QThreadPool threadPool;
    static const int SIZE = 5;
};

#endif // LABCONTROLLER_H
