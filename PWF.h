#ifndef PWF_H
#define	PWF_H

#include "sys_param.h"
class PWF {
public:
    PWF(double* potential, double* wf1, double* wf2, double* wf3, double* wf4,
            int numPoints, int nuclearCharge, double* energyLevel);
    PWF(const PWF& orig);
    virtual ~PWF();
    double* getScaledPotential(const int nuclearCharge);
    double* getScaledWF0(const int nuclearCharge);
    double* getScaledWF1(const int nuclearCharge);
    double* getScaledWF2(const int nuclearCharge);
    double* getScaledWF3(const int nuclearCharge);
    double* getEnergyLevel();
    int nuclearCharge;

private:
    class scaledPWF {
    public:
        scaledPWF() {
            scale = -1.0;
        }
        scaledPWF(const scaledPWF& orig) {
            scale = orig.scale;
            potentialPtr = orig.potentialPtr;
            wFBand0Ptr = orig.wFBand0Ptr;
            wFBand1Ptr = orig.wFBand1Ptr;
            wFBand2Ptr = orig.wFBand2Ptr;
            wFBand3Ptr = orig.wFBand3Ptr;
            radialPoints = orig.radialPoints;
        }
        double scale;
        double* potentialPtr;
        double* wFBand0Ptr;
        double* wFBand1Ptr;
        double* wFBand2Ptr;
        double* wFBand3Ptr;
        double* radialPoints;

        ~scaledPWF() {
        }
    };
    scaledPWF nativePWF;
    vector<scaledPWF> scaleContainer;
    scaledPWF* getScaledPWF(const int nuclearCharge);
    scaledPWF* makeRescalePWF(const double scale);
    double getScale(const int);
    void* setGrid(double scale, scaledPWF* futureScaledPWF);
    
    int numberofPoints;
    double* _dEnergyLevel;

};

#endif	/* PWF_H */

