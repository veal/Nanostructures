#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    QFileDialog::Options options;
    QString selectedFilter;
    QString fileName = QFileDialog::getOpenFileName(this,
                                tr("QFileDialog::getOpenFileName()"),
                                "/home/veal/",
                                tr("All Files (*);;Data Files (*.dat)"),
                                &selectedFilter,
                                options);
    if (!fileName.isEmpty()) {
        ui->lineEdit->setText(fileName);
        matrixStr = fileName;
    }
}
void MainWindow::setCallBackListener(Controller* controller) {
    _controller = controller;
}

void MainWindow::on_pushButton_2_clicked() {
    QFileDialog::Options options;
    QString selectedFilter;
    QString fileName = QFileDialog::getOpenFileName(this,
                                tr("QFileDialog::getOpenFileName()"),
                                "/home/veal/",
                                tr("All Files (*);;Data Files (*.dat)"),
                                &selectedFilter,
                                options);
    if (!fileName.isEmpty()) {
        ui->lineEdit_2->setText(fileName);
        impurityStr = fileName;
    }
}

void MainWindow::on_pushButton_3_clicked() {
    _controller->set_matrixInputFileStr(ui->lineEdit->text());
    _controller->set_impurityInputFileStr(ui->lineEdit_2->text());
    _controller->notifyStartcalculating();
}
QString MainWindow::get_matrixStr(){
    return matrixStr;
}
QString MainWindow::get_impurityStr(){
    return impurityStr;
}
