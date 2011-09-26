#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <QString>

class Controller
{
public:
    Controller();

    virtual void set_matrixInputFileStr(QString) = 0;
    virtual void set_impurityInputFileStr(QString) = 0;
    virtual void notifyStartcalculating() = 0;
    virtual void matrixFileChosen(QString) = 0;
    virtual void impurityFileChosen(QString) = 0;
    virtual void setMainWindow(void*) = 0;
    virtual void HamCalcProgressUpdate(int) = 0;

};

#endif // CONTROLLER_H
