#ifndef INITCLASS_H
#define INITCLASS_H

#include <complex>
#include <QFile>
#include <QTextStream>
#include "GlobalConstVar.h"

class InitClass
{
public:
    InitClass();
    void Set_coord();
    void Read_Coulomb_integral();
    void Read_input_file();
    void Scat_Pr();
    void allocate();

};

#endif // INITCLASS_H
