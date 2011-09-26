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
                                "/home/veal/Sandbox/GUINano",
                                tr("All Files (*);;Data Files (*.dat)"),
                                &selectedFilter,
                                options);
    if (!fileName.isEmpty()) {
        _controller->matrixFileChosen(fileName);
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
                                "/home/veal/Sandbox/GUINano",
                                tr("All Files (*);;Data Files (*.dat)"),
                                &selectedFilter,
                                options);
    if (!fileName.isEmpty()) {
        _controller->impurityFileChosen(fileName);
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

void MainWindow::setMatrixLabel(QString str) {
   ui->matrixElementLabel->setText(str);
}
void MainWindow::setMatrixNumOrb(int num) {
    QString str;
    QTextStream sstr(&str);
    sstr << num;
    ui->matrixNumOrbLabel->setText(str);
}
void MainWindow::setMatrixDist(double num) {
    QString str;
    QTextStream sstr(&str);
    sstr << num;
    ui->matrixDistLabel->setText(str);
}
void MainWindow::setMatrixTruncRad(double num) {
    QString str;
    QTextStream sstr(&str);
    sstr << num;
    ui->matrixRadLabel->setText(str);
}
void MainWindow::setMatrixPWFPath(QString str) {
   ui->matrixPwfPathLabel->setText(str);
}
void MainWindow::setImpurLabel(QString str) {
   ui->impurElementLabel->setText(str);
}
void MainWindow::setImpurNumOrb(int num) {
    QString str;
    QTextStream sstr(&str);
    sstr << num;
    ui->impurNumOrbLabel->setText(str);
}
void MainWindow::setImpurDist(double num) {
    QString str;
    QTextStream sstr(&str);
    sstr << num;
    ui->impurDistLabel->setText(str);
}
void MainWindow::setImpurTruncRad(double num) {
    QString str;
    QTextStream sstr(&str);
    sstr << num;
    ui->impurRadLabel->setText(str);
}
void MainWindow::setImpurPWFPath(QString str) {
    ui->impurPwfPathLabel->setText(str);
}
void MainWindow::setHamProgress(int percent) {
    ui->progressBar->setValue(percent);
}
