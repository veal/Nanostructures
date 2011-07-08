#ifndef LAB_H
#define LAB_H

#include "sys_param.h"
#include <QString>
#include <QRunnable>

class Lab : public QRunnable
{
public:
    Lab();
    void runCalculation();
    void set_matrixInputFile(QString);
    void set_impurityInputFile(QString);

    void run();

private:
    QString matrixInputFile;
    QString impurityInputFile;
};

#endif // LAB_H
