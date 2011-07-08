#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
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
