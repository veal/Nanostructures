#include "PWF.h"

PWF::PWF(double* potential, double* wf0, double* wf1, double* wf2, double* wf3,
            int numPoints, int nuclearCharge, double* energyLevel) {
    nativePWF.potentialPtr = potential;
    nativePWF.wFBand0Ptr = wf0;
    nativePWF.wFBand1Ptr = wf1;
    nativePWF.wFBand2Ptr = wf2;
    nativePWF.wFBand3Ptr = wf3;
    numberofPoints = numPoints;
    nativePWF.radialPoints = new double[numPoints];
    nativePWF.scale = getScale(nuclearCharge);
    setGrid(nativePWF.scale, &nativePWF);

    scaleContainer.push_back(nativePWF);
    _dEnergyLevel = energyLevel;
    this->nuclearCharge = nuclearCharge;
}

PWF::PWF(const PWF& orig) {

}

PWF::~PWF() {
    for (int iter = 0; iter < scaleContainer.size(); iter++) {
        delete scaleContainer[iter].potentialPtr;
        delete scaleContainer[iter].radialPoints;
        delete scaleContainer[iter].wFBand0Ptr;
        delete scaleContainer[iter].wFBand1Ptr;
        delete scaleContainer[iter].wFBand2Ptr;
        delete scaleContainer[iter].wFBand3Ptr;
    }
    scaleContainer.clear();
}

double PWF::getScale(const int numberOfElectrons){
    return 0.88534138/pow(numberOfElectrons, 1.0/3.0);
}

PWF::scaledPWF* PWF::getScaledPWF(const int nuclearCharge) {
    double scale = getScale(nuclearCharge);
    for (int iter = 0; iter < scaleContainer.size(); iter++) {
        if (abs(scaleContainer[iter].scale - scale) < scale*0.01) {
            return &(scaleContainer[iter]);
        }
    }

    return makeRescalePWF(scale);
}

PWF::scaledPWF* PWF::makeRescalePWF(const double scale) {
    scaledPWF futureScaledPWF;
    futureScaledPWF.scale = scale;
    futureScaledPWF.potentialPtr = new double[numberofPoints];
    futureScaledPWF.radialPoints = new double[numberofPoints];
    futureScaledPWF.wFBand0Ptr = new double[numberofPoints];
    futureScaledPWF.wFBand1Ptr = new double[numberofPoints];
    futureScaledPWF.wFBand2Ptr = new double[numberofPoints];
    futureScaledPWF.wFBand3Ptr = new double[numberofPoints];
    setGrid(scale, &futureScaledPWF);
    
    for (int iter = 0; iter < numberofPoints; iter++) {
        for (int iter2 = 0; iter2 < numberofPoints-1; iter2++) {
            if (futureScaledPWF.radialPoints[iter] > nativePWF.radialPoints[iter2] &&
                    futureScaledPWF.radialPoints[iter] < nativePWF.radialPoints[iter2+1]) {
                double alpha = (futureScaledPWF.radialPoints[iter] - nativePWF.radialPoints[iter2]) /
                    (nativePWF.radialPoints[iter2+1] - nativePWF.radialPoints[iter2]);
                if (alpha < 0) cout << "alpha < 0!!!\n";
                futureScaledPWF.potentialPtr[iter] = alpha*nativePWF.potentialPtr[iter2] +
                        (1 - alpha)*nativePWF.potentialPtr[iter2+1];
                futureScaledPWF.wFBand0Ptr[iter] = alpha*nativePWF.wFBand0Ptr[iter2] +
                        (1 - alpha)*nativePWF.wFBand0Ptr[iter2+1];
                futureScaledPWF.wFBand1Ptr[iter] = alpha*nativePWF.wFBand1Ptr[iter2] +
                        (1 - alpha)*nativePWF.wFBand1Ptr[iter2+1];
                futureScaledPWF.wFBand2Ptr[iter] = alpha*nativePWF.wFBand2Ptr[iter2] +
                        (1 - alpha)*nativePWF.wFBand2Ptr[iter2+1];
                futureScaledPWF.wFBand3Ptr[iter] = alpha*nativePWF.wFBand3Ptr[iter2] +
                        (1 - alpha)*nativePWF.wFBand3Ptr[iter2+1];
            }
        }
    }

    futureScaledPWF.potentialPtr[0] = nativePWF.potentialPtr[1];
    futureScaledPWF.wFBand0Ptr[0] = nativePWF.wFBand0Ptr[1];
    futureScaledPWF.wFBand1Ptr[0] = nativePWF.wFBand1Ptr[1];
    futureScaledPWF.wFBand2Ptr[0] = nativePWF.wFBand2Ptr[1];
    futureScaledPWF.wFBand3Ptr[0] = nativePWF.wFBand3Ptr[1];

    scaleContainer.push_back(futureScaledPWF);
    return &(scaleContainer.back());
}
void* PWF::setGrid(double scale, scaledPWF* futureScaledPWF) {
    int m_Nblock = numberofPoints/40;
    int d_R = 0;
    double d_DeltaX = 0.0025;
    double* d_X = new double[numberofPoints];

    d_X[0] = 0.0;

    for (int j = 0; j < m_Nblock; j++) {
        for (int k = 0; k < 40; k++) {
            d_R++;
            d_X[d_R] = d_X[d_R-1] + d_DeltaX;
            futureScaledPWF->radialPoints[d_R-1] = d_X[d_R-1]*scale;
        }
        d_DeltaX += d_DeltaX;
    }

    futureScaledPWF->radialPoints[d_R] = d_X[d_R]*scale;

    delete[] d_X;
}
double* PWF::getScaledPotential(const int nuclearCharge) {
    return getScaledPWF(nuclearCharge)->potentialPtr;
}
double* PWF::getScaledWF0(const int nuclearCharge) {
    return getScaledPWF(nuclearCharge)->wFBand0Ptr;
}
double* PWF::getScaledWF1(const int nuclearCharge) {
    return getScaledPWF(nuclearCharge)->wFBand1Ptr;
}
double* PWF::getScaledWF2(const int nuclearCharge) {
    return getScaledPWF(nuclearCharge)->wFBand2Ptr;
}
double* PWF::getScaledWF3(const int nuclearCharge) {
    return getScaledPWF(nuclearCharge)->wFBand3Ptr;
}
double* PWF::getEnergyLevel() {
    return _dEnergyLevel;
}
