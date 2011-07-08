#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <QString>

class Controller
{
public:
    Controller();

    QString matrixInputFileStr;
    QString impurityInputFileStr;

    virtual void set_matrixInputFileStr(QString) = 0;
    virtual void set_impurityInputFileStr(QString) = 0;
    virtual void notifyStartcalculating() = 0;

};

#endif // CONTROLLER_H
