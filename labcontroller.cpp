#include "labcontroller.h"

LabController::LabController() {
    threadPool.setMaxThreadCount(SIZE);
}
void LabController::set_matrixInputFileStr(QString str) {
    matrixInputFileStr = str;
}
void LabController::set_impurityInputFileStr(QString str) {
    impurityInputFileStr = str;
}
void LabController::notifyStartcalculating() {
    threadPool.start(&lab);
}
