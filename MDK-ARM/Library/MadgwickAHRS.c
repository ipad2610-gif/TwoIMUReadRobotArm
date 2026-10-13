//=============================================================================================
// MadgwickAHRS.c
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
// 19/02/2012	SOH Madgwick	Magnetometer measurement is normalised
//
//=============================================================================================

//-------------------------------------------------------------------------------------------
// Header files

#include "MadgwickAHRS.h"
#include <math.h>
#include "utilzfunc.h"
//-------------------------------------------------------------------------------------------
// Definitions

#define sampleFreqDef   100.0f          // sample frequency in Hz
#define betaDef         0.1f            // 2 * proportional gain


//============================================================================================
// Functions

//-------------------------------------------------------------------------------------------
// AHRS algorithm update

void InitFilter(MadgwickFilter *Filter){
    Filter->beta = betaDef;
    Filter->q0.q0 = 1.0f; Filter->q0.q1 = 0.0f; Filter->q0.q2 = 0.0f; Filter->q0.q3 = 0.0f;
    Filter->qCurr.q0 = 1.0f; Filter->qCurr.q1 = 0.0f; Filter->qCurr.q2 = 0.0f; Filter->qCurr.q3 = 0.0f;
    Filter->invSampleFreq = 0.01;
    Filter->anglesComputed = 0;
}


// IMU algorithm update
void updateIMU(MadgwickFilter *Filter, float gx, float gy, float gz, float ax, float ay, float az) {
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

	// Convert gyroscope degrees/sec to radians/sec
	gx *= 0.0174533f;
	gy *= 0.0174533f;
	gz *= 0.0174533f;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-Filter->qCurr.q1 * gx - Filter->qCurr.q2 * gy - Filter->qCurr.q3 * gz);
	qDot2 = 0.5f * (Filter->qCurr.q0 * gx + Filter->qCurr.q2 * gz - Filter->qCurr.q3 * gy);
	qDot3 = 0.5f * (Filter->qCurr.q0 * gy - Filter->qCurr.q1 * gz + Filter->qCurr.q3 * gx);
	qDot4 = 0.5f * (Filter->qCurr.q0 * gz + Filter->qCurr.q1 * gy - Filter->qCurr.q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0f * Filter->qCurr.q0;
		_2q1 = 2.0f * Filter->qCurr.q1;
		_2q2 = 2.0f * Filter->qCurr.q2;
		_2q3 = 2.0f * Filter->qCurr.q3;
		_4q0 = 4.0f * Filter->qCurr.q0;
		_4q1 = 4.0f * Filter->qCurr.q1;
		_4q2 = 4.0f * Filter->qCurr.q2;
		_8q1 = 8.0f * Filter->qCurr.q1;
		_8q2 = 8.0f * Filter->qCurr.q2;
		q0q0 = Filter->qCurr.q0 * Filter->qCurr.q0;
		q1q1 = Filter->qCurr.q1 * Filter->qCurr.q1;
		q2q2 = Filter->qCurr.q2 * Filter->qCurr.q2;
		q3q3 = Filter->qCurr.q3 * Filter->qCurr.q3;

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * Filter->qCurr.q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0f * q0q0 * Filter->qCurr.q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0f * q1q1 * Filter->qCurr.q3 - _2q1 * ax + 4.0f * q2q2 * Filter->qCurr.q3 - _2q2 * ay;
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= Filter->beta * s0;
		qDot2 -= Filter->beta * s1;
		qDot3 -= Filter->beta * s2;
		qDot4 -= Filter->beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	Filter->qCurr.q0 += qDot1 * Filter->invSampleFreq;
	Filter->qCurr.q1 += qDot2 * Filter->invSampleFreq;
	Filter->qCurr.q2 += qDot3 * Filter->invSampleFreq;
	Filter->qCurr.q3 += qDot4 * Filter->invSampleFreq;

	// Normalise quaternion
	recipNorm = invSqrt(Filter->qCurr.q0 * Filter->qCurr.q0 + Filter->qCurr.q1 * Filter->qCurr.q1 + Filter->qCurr.q2 * Filter->qCurr.q2 + Filter->qCurr.q3 * Filter->qCurr.q3);
	Filter->qCurr.q0 *= recipNorm;
	Filter->qCurr.q1 *= recipNorm;
	Filter->qCurr.q2 *= recipNorm;
	Filter->qCurr.q3 *= recipNorm;
	Filter->anglesComputed = 0;
}

float invSqrt(float x) {
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
	y = y * (1.5f - (halfx * y * y));
	return y;
}

//-------------------------------------------------------------------------------------------

//void computeAngles()
//{
//	roll = atan2f(qCurr.q0*qCurr.q1 + qCurr.q2*qCurr.q3, 0.5f - qCurr.q1*qCurr.q1 - qCurr.q2*qCurr.q2);
//	pitch = asinf(-2.0f * (qCurr.q1*qCurr.q3 - qCurr.q0*qCurr.q2));
//	yaw = atan2f(qCurr.q1*qCurr.q2 + qCurr.q0*qCurr.q3, 0.5f - qCurr.q2*qCurr.q2 - qCurr.q3*qCurr.q3);
//	anglesComputed = 1;
//}

void computeTwistBendAngles(MadgwickFilter *Filter)
{
	Quaternion q0_inv = invQuaternion(Filter->q0);
	Quaternion rel = mulQuaternion(q0_inv, Filter->qCurr);
	float x = 1; float y = 0; float z = 0; // Twist around x
    float u = x*rel.q1 + y*rel.q2 + z*rel.q3;
	float n = x*x + y*y + z*z;
	float m = rel.q0*n;
	float l = sqrt(m*m + u*u*n);
	if (l > 1e-6f) {
	    Filter->twist.q0 = m/l;
	    Filter->twist.q1 = x*u/l;
	    Filter->twist.q2 = y*u/l;
	    Filter->twist.q3 = z*u/l;
    } else {
        Filter->twist.q0 = 1.0f;
        Filter->twist.q1 = 0.0f;
        Filter->twist.q2 = 0.0f;
        Filter->twist.q3 = 0.0f;
    }
	Filter->bending = mulQuaternion(rel,invQuaternion(Filter->twist));
	float twist_rad = 2.0f * atan2(Filter->twist.q1, Filter->twist.q0);
	Filter->twist_angle = twist_rad * (180.0f / M_PI);

    float flexion_rad= 2.0f * atan2(Filter->bending.q2, Filter->bending.q0);
    Filter->bending_angle = flexion_rad * (180.0f / M_PI);
}

Quaternion mulQuaternion(Quaternion q, Quaternion r)
{
    Quaternion result;
    result.q0 = r.q0*q.q0 - r.q1*q.q1 - r.q2*q.q2 - r.q3*q.q3;
    result.q1 = r.q0*q.q1 + r.q1*q.q0 - r.q2*q.q3 + r.q3*q.q2;
    result.q2 = r.q0*q.q2 + r.q1*q.q3 + r.q2*q.q0 - r.q3*q.q1;
    result.q3 = r.q0*q.q3 - r.q1*q.q2 + r.q2*q.q1 + r.q3*q.q0;
    return result;
}
Quaternion invQuaternion(Quaternion x)
{
    Quaternion result;
    float norm = x.q0*x.q0 + x.q1*x.q1 + x.q2*x.q2 + x.q3*x.q3;
    if (norm < 1e-6f) return x;
    result.q0 = x.q0/norm; result.q1 = -x.q1/norm;
    result.q2 = -x.q2/norm; result.q3 = -x.q3/norm;
    return result;
}
