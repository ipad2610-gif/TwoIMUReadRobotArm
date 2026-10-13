//=============================================================================================
#ifndef ESTIMATE_h
#define ESTIMATE_h
#include <math.h>
#include "stdbool.h"
#include "stdint.h"
#include "vqf6d.h"
//--------------------------------------------------------------------------------------------
#define M_PI 3.1415926
#define U1_MIN_GYRO_RAD_S       0.05f
#define U1_MIN_SAMPLES          100U
#define U1_CONFIDENCE_MIN       0.95f
#define U1_POWER_ITERATIONS     8U
#define U1_EPS                  1e-9f

#define U2_Q1_STILL_RAD_S          0.08f
#define U2_MIN_GYRO_RAD_S          0.10f
#define U2_MIN_SAMPLES             100U
#define U2_CONFIDENCE_MIN          0.95f
#define U2_MIN_MOTION_RAD          0.5f

#define HEADING_BUFFER_SIZE        200U
#define HEADING_SAMPLE_DIV         5U

#define HEADING_INIT_COARSE_RAD    0.0872664626f
#define HEADING_INIT_FINE_RAD      0.00872664626f

#define HEADING_TRACK_RANGE_RAD    0.174532925f
#define HEADING_TRACK_STEP_RAD     0.00436332313f
#define HEADING_MAX_STEP_RAD 0.008726646f   // 0.5 deg/update

typedef struct
{
	float gyro1_norm;
	float gyro2_norm;

	float u1_axis_error_deg;
	float u1_perp_norm;
	float u1_perp_ratio;

	float j1_3D[3];
	float j1_6D[3];

	float j1_3D_h;
	float j1_3D_z;
	float j1_6D_h;
	float j1_6D_z;

	float omega1_E1[3];
	float omega2_E2[3];

	float omega1_h;
	float omega2_h;
	float norm_v,norm_v_ratio;
	float heading_direct;
	float R1_ref[3][3], R2_ref[3][3];
	float drift1, drift2, relative_drift;
} EstimatorDebug;

typedef struct
{
	float bending_angle, twist_angle;
	float q1_pred,q1dot_pred, q2_pred, q2dot_pred;
	float invSampleFreq;
	uint8_t gravity_valid;
	uint32_t u1_sample_count;
	float u1_mean[3];
	float u1_M2[3][3];
	float u1[3];
	float u1_confidence;
	uint8_t u1_valid;
	uint8_t u1_locked;
	uint8_t u1_compute_counter;
	
	float u2[3];
	float u2_mean[3];
	float u2_M2[3][3];
	float u2_confidence;
	float u2_motion;
	uint8_t u2_valid;
	uint8_t u2_locked;
	uint8_t u2_compute_counter;
	uint32_t u2_sample_count;
	float heading_delta;
	float heading_cost;
	float heading_excitation;
	uint8_t compute_heading_counter;
	uint8_t heading_valid;

	uint8_t heading_sample_div_counter;
	uint16_t heading_buffer_count;
	uint16_t heading_buffer_head;
	
	float heading_rel_norm;
	float heading_perp_norm;
	float heading_parallel;
	
	float heading_A[HEADING_BUFFER_SIZE][3];
	float heading_B[HEADING_BUFFER_SIZE][3];
	float heading_D[HEADING_BUFFER_SIZE][3];
	
	EstimatorDebug debug;
}Estimator;



void InitEstimator(Estimator *est);
void UpdateJointAxis1(Estimator *est, const float *gyro);
uint8_t ComputeJointAxis1(Estimator *est);
void ResetJointAxis1(Estimator *est);
void ResetJointAxis2(Estimator *est);
void UpdateJointAxis2(Estimator *est, const float *gyro1, const float *gyro2);
uint8_t ComputeJointAxis2(Estimator *est);
void UpdateHeadingData(Estimator *est, const float *gyro1, const float *gyro2, const float R1[3][3], const float R2[3][3]);
uint8_t EstimateInitialHeading(Estimator *est);
uint8_t TrackHeadingDelta(Estimator *est);
void ComputeHeadingDiagnostic(Estimator *est, const float *gyro1, const float *gyro2, const float R1[3][3], const float R2[3][3]);
void EstimatorDebugUpdate(Estimator *est, VQF6D *vqf1, const float *gyro1, const float *gyro2, const float R1[3][3], const float R2[3][3]);

#endif /*ESTIMATE_H*/

