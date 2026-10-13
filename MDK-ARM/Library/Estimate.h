//=============================================================================================
#ifndef ESTIMATE_h
#define ESTIMATE_h
#include <math.h>
#include "stdbool.h"
#include "stdint.h"
#include "vqf6d.h"
//-------------------------------------------------------------------------------------------
// Definitions
#define sampleFreqDef   100.0f          // sample frequency in Hz
#define RAD2DEG 				57.2957795f
#define DEG2RAD					0.0174533f
#define ALPHA2_DEG 	87.92f //88.19f
#define PHI_DEG     1.81f

#define SA 0.99934112f  	 // sin(alpha2)
#define CA 0.03629488f   // cos(alpha2)

#define ALPHA2_FC_HZ   1.0f
#define Q2KF_SIGMA_RATE       0.0035f
#define Q2KF_SIGMA_ACC        0.0175f

#define U1_MIN_GYRO_RAD_S       0.1f
#define U1_MIN_SAMPLES          100U
#define U1_CONFIDENCE_MIN       0.95f
#define U1_POWER_ITERATIONS     8U
#define U1_EPS                  1e-9f

#define U2_Q1_STILL_RAD_S          0.08f
#define U2_MIN_GYRO_RAD_S          0.10f
#define U2_MIN_SAMPLES             100U
#define U2_CONFIDENCE_MIN          0.95f
#define U2_MIN_MOTION_RAD          0.5f

typedef struct
{
	float twist_angle, bending_angle;

	float q1_acc, q1_pred,q1dot_pred;
	float q2_acc, q2_bias;
	float	q2_pred, q2dot_pred;
	float dq2_gyro,dq2_corr;
	float ePred_g1, ePred_g2;
	float a1_norm, a2_norm, h2_norm, gyro1_norm, gyro2_norm;
	
	float P;	// Covariance  of Prediction
	float Q;	// Covariance  of measurement
	float q2_R;
	float q2_Kq;
	// Compensate linear acceleration IMU2
	float r2[3];
	float gyro2_prev[3];
	float alpha2_f[3];
	float a_dyn2[3];
	float a_tan2[3];
	float a_cent2[3];
	float acc_comp[3];
	float dyn_comp_gain;
	uint8_t dyn_comp_initialized;
	float c12;
	float j1[3],j2[3];
	float g10[3],g20[3],j12[3];
	float g2_ref[3];
	float invSampleFreq;
	float lambda_q1, lambda_q2, Kq1, Kq2,w1,w2;
	
}Estimator;

typedef struct
{
	float dt;
	// Calib j10
	uint32_t u1_sample_count;
	float u1_mean[3];
	float u1_M2[3][3];
	float u1[3];
	float u1_confidence;
	uint8_t u1_valid;
	uint8_t u1_locked;
	uint8_t u1_compute_counter;
	// Calib j20
	float u2[3];
	float u2_mean[3];
	float u2_M2[3][3];
	float u2_confidence;
	float u2_motion;
	uint8_t u2_valid;
	uint8_t u2_locked;
	uint8_t u2_compute_counter;
	uint32_t u2_sample_count;
		// Calib j12
	float j12_num[3];
	float j12_den;
	uint32_t j12_count;
	float j12[3];
	float j12_gyro2_energy;
	float j12_confidence;
	float j12_scale;
	uint8_t j12_valid;
	float j12_motion;
}CalibProcess;

void InitEstimator(Estimator *est);
void CompensateLinearAcc(Estimator *est, float *gyro2, float *gValue2);
void ComputeAngles(Estimator *est, float *gyro1, float *gyro2, float *gValue1, float *gValue2);

void InitCalibProcess(CalibProcess *calib);

void ResetJointAxis1(CalibProcess *calib);
void UpdateJointAxis1(CalibProcess *calib, const float *gyro);
uint8_t ComputeJointAxis1(CalibProcess *calib); 

void ResetJointAxis2(CalibProcess *calib);
void UpdateJointAxis2(CalibProcess *calib, const float *gyro1, const float *gyro2);
uint8_t ComputeJointAxis2(CalibProcess *calib);

void ResetJ12Axis(CalibProcess *calib);
void ComputeJ12Axis(CalibProcess *calib, float *gyro1, float *gyro2);

#endif /*ESTIMATE_H*/

