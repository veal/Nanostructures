/* 
 * File:   Log.h
 * Author: veal
 *
 * Created on May 22, 2011, 10:23 AM
 */

#include <fstream>
#include <string>
#include <iostream>
#include <complex>

#ifndef LOG_H
#define	LOG_H

class Log {
public:
    Log();
    Log(const Log& orig);
    void WLog();
    virtual ~Log();
private:

};

#endif	/* LOG_H */

