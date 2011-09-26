#include "labcontroller.h"

LabController::LabController() {
    threadPool.setMaxThreadCount(SIZE);
    lab.setProgressUpdateCallback(this);
}
void LabController::set_matrixInputFileStr(QString str) {
    matrixInputFileStr = str;
}
void LabController::set_impurityInputFileStr(QString str) {
    impurityInputFileStr = str;
}
void LabController::notifyStartcalculating() {
    lab.set_matrixInputFile(matrixInputFileStr);
    lab.set_impurityInputFile(impurityInputFileStr);
    threadPool.start(&lab);
}
void LabController::matrixFileChosen(QString str) {
    fillMatrixFields(lab.getHopIntegrals()->readParamFile(str));
}
void LabController::impurityFileChosen(QString str) {
    fillImpurFields(lab.getHopIntegrals()->readParamFile(str));
}
void LabController::fillMatrixFields(Hop_integrals::paramFile *param) {
    window->setMatrixLabel(param->elementLabel);
    window->setMatrixNumOrb(param->m_Norbital);
    window->setMatrixDist(param->latticeConst);
    window->setMatrixTruncRad(param->dr_max);
    window->setMatrixPWFPath(param->sz_InputFile);
}
void LabController::fillImpurFields(Hop_integrals::paramFile *param) {
    window->setImpurLabel(param->elementLabel);
    window->setImpurNumOrb(param->m_Norbital);
    window->setImpurDist(param->latticeConst);
    window->setImpurTruncRad(param->dr_max);
    window->setImpurPWFPath(param->sz_InputFile);
}
void LabController::setMainWindow(void *w) {
    window = static_cast<MainWindow*> (w);
}
void LabController::HamCalcProgressUpdate(int percent) {
    window->setHamProgress(percent);
}
