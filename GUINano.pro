#-------------------------------------------------
#
# Project created by QtCreator 2011-07-07T12:41:51
#
#-------------------------------------------------

QT       += core gui

TARGET = GUINano
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    W_matrix.cpp \
    Wamiltonian.cpp \
    Speed.cpp \
    sigma_fe.cpp \
    reading_from_files.cpp \
    PWF.cpp \
    Math.cpp \
    Magnetic_moments.cpp \
    Log.cpp \
    lab.cpp \
    Hop_integrals.cpp \
    Hamiltonian.cpp \
    Green_phonon.cpp \
    Green_electron.cpp \
    Density.cpp \
    Conductivity.cpp \
    Coherent_pot.cpp \
    allocation.cpp \
    controller.cpp \
    labcontroller.cpp \
    task.cpp

HEADERS  += mainwindow.h \
    W_matrix.h \
    sys_param.h \
    PWF.h \
    nr3.h \
    ludcmp.h \
    LUcomplex.h \
    Log.h \
    lab.h \
    Hop_integrals.h \
    eigen_sym.h \
    ceigen_sym.h \
    controller.h \
    labcontroller.h \
    task.h

FORMS    += mainwindow.ui
