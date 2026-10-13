//=============================================================================================
// MadgwickAHRS.h
//=============================================================================================
//
// Implementation of Madgwick's IMU and AHRS algorithms.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// From https://web.archive.org/web/20190712041101/http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/:
// "Open-source resources available on this website are provided under the
// [GNU General Public Licence](https://web.archive.org/web/20190712041101/http://www.gnu.org/licenses/gpl.html)
// unless an alternative licence is provided in source."
//
// Date			Author          Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
//
//=============================================================================================
#ifndef MadgwickAHRS_h
#define MadgwickAHRS_h
#include <math.h>
#include "stdbool.h"
//--------------------------------------------------------------------------------------------
// Variable declaration
static float invSqrt(float x);
typedef struct {
    float q0;
    float q1;
    float q2;
    float q3;
}Quaternion;

typedef struct
{
    float beta;             // algorithm gain
    float invSampleFreq;
    Quaternion q0;
    Quaternion qCurr;
    Quaternion twist;
    Quaternion bending;
    bool anglesComputed;
    float bending_angle, twist_angle;
}MadgwickFilter;


void InitFilter(MadgwickFilter *Filter);
//void computeAngles();
void updateIMU(MadgwickFilter *Filter, float gx, float gy, float gz, float ax, float ay, float az);
void computeTwistBendAngles(MadgwickFilter *Filter);
Quaternion mulQuaternion(Quaternion q, Quaternion r);
Quaternion invQuaternion(Quaternion x);
float invSqrt(float x);
#endif

