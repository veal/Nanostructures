#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTextStream>
#include <QtGui/QFileDialog>
#include "controller.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void setCallBackListener(Controller*);
    ~MainWindow();
    QString get_matrixStr();
    QString get_impurityStr();

    void setMatrixLabel(QString);
    void setMatrixNumOrb(int);
    void setMatrixTruncRad(double);
    void setMatrixDist(double);
    void setMatrixPWFPath(QString);

    void setImpurLabel(QString);
    void setImpurNumOrb(int);
    void setImpurTruncRad(double);
    void setImpurDist(double);
    void setImpurPWFPath(QString);

    void setHamProgress(int);

private:
    Ui::MainWindow *ui;
    Controller *_controller;
    QString matrixStr;
    QString impurityStr;

private slots:
    void on_pushButton_3_clicked();
    void on_pushButton_2_clicked();
    void on_pushButton_clicked();
};

#endif // MAINWINDOW_H
