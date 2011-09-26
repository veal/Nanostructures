#include <QtGui/QApplication>
#include "mainwindow.h"
#include "labcontroller.h"

pthread_mutex_t job_queue = PTHREAD_MUTEX_INITIALIZER;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;

    Controller* controller = new LabController();
    w.setCallBackListener(controller);
    controller->setMainWindow(&w);
    w.show();

    return a.exec();
}
