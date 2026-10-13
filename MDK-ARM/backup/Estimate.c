//=============================================================================================
// ESTIMATE.c

#include "Estimate.h"
#include <math.h>

//-------------------------------------------------------------------------------------------
// Definitions

#define sampleFreqDef   100.0f          // sample frequency in Hz
#define ALPHA2_DEG 88.19f
#define PHI_DEG     1.81f

#define SA 0.99950106f   // sin(alpha2)
#define CA 0.03158521f   // cos(alpha2)

#define SP 0.03158521f   // sin(phi)
#define CP 0.99950106f

//============================================================================================
// Functions

//-------------------------------------------------------------------------------------------

static float Norm3(const float *v)
{
	return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static float Dot3(const float *a, const float *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static uint8_t Normalize3(float *v)
{
	float norm = Norm3(v);

	if (norm < U1_EPS) return 0U;

	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;

	return 1U;
}

static void Mat3Vec(const float A[3][3], const float *x, float *y)
{
	// y = A*x
	y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
	y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
	y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}
static void Mat3TVec(const float A[3][3], const float *x, float *y)
{
  // y = A^T*x
	y[0] = A[0][0]*x[0] + A[1][0]*x[1] + A[2][0]*x[2];
	y[1] = A[0][1]*x[0] + A[1][1]*x[1] + A[2][1]*x[2];
	y[2] = A[0][2]*x[0] + A[1][2]*x[1] + A[2][2]*x[2];
}
static float WrapPi(float angle)
{
	while (angle > M_PI) angle -= 2.0f*M_PI;
	while (angle < -M_PI) angle += 2.0f*M_PI;

	return angle;
}
static void Mat3MulTransB(const float A[3][3], const float B[3][3], float C[3][3])
{
    for (uint8_t i = 0U; i < 3U; i++)
    {
        for (uint8_t j = 0U; j < 3U; j++)
        {
            C[i][j] = A[i][0]*B[j][0] + A[i][1]*B[j][1] + A[i][2]*B[j][2];
        }
    }
}
static uint8_t PrincipalAxis3(const float C[3][3], float *axis, float *lambda)
{
	// Find axis based on eigenvalue of covariance matrix C
	float trace = C[0][0] + C[1][1] + C[2][2];

	if (trace < U1_EPS) return 0U;

	float v[3];

	if ((C[0][0] >= C[1][1]) && (C[0][0] >= C[2][2]))
	{
		v[0] = 1.0f;
		v[1] = 0.0f;
		v[2] = 0.0f;
	}
	else if (C[1][1] >= C[2][2])
	{
		v[0] = 0.0f;
		v[1] = 1.0f;
		v[2] = 0.0f;
	}
	else
	{
		v[0] = 0.0f;
		v[1] = 0.0f;
		v[2] = 1.0f;
	}

	for (uint8_t iteration = 0U; iteration < U1_POWER_ITERATIONS; iteration++)
	{
		float temp[3];

		Mat3Vec(C, v, temp);

		if (!Normalize3(temp)) return 0U;

		v[0] = temp[0];
		v[1] = temp[1];
		v[2] = temp[2];
	}

	float Cv[3];

	Mat3Vec(C, v, Cv);

	*lambda = Dot3(v, Cv);

	axis[0] = v[0];
	axis[1] = v[1];
	axis[2] = v[2];

	return 1U;
}

void InitEstimator(Estimator *est){
	est->invSampleFreq = 0.01f;
	ResetJointAxis1(est);
	ResetJointAxis2(est);
	est->compute_heading_counter = 0;

	est->heading_delta = 0.0f;
    est->heading_cost = 0.0f;
    est->heading_excitation = 0.0f;
    est->heading_valid = 0U;

    est->compute_heading_counter = 0U;
    est->heading_sample_div_counter = 0U;
    est->heading_buffer_count = 0U;
    est->heading_buffer_head = 0U;
}


void ResetJointAxis1(Estimator *est)
{
	est->u1_sample_count = 0U;

	est->u1_mean[0] = 0.0f;
	est->u1_mean[1] = 0.0f;
	est->u1_mean[2] = 0.0f;

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			est->u1_M2[i][j] = 0.0f;
		}
	}

	est->u1[0] = 0.0549092f;
	est->u1[1] = 0.991894947f;
	est->u1[2] = 0.114592113f;

	est->u1_confidence = 0.0f;
	est->u1_valid = 1U;
	est->u1_locked = 1U;
	est->u1_compute_counter = 0U;
}

void UpdateJointAxis1(Estimator *est, const float *gyro)
{
	// Builf Welford covariance
	if (est->u1_locked) return;

	float gyro_mag = Norm3(gyro);

	if (gyro_mag < U1_MIN_GYRO_RAD_S) return;
	
	// remove gyro bias
	est->u1_sample_count++;

	float n = (float)est->u1_sample_count;

	float delta[3];
	float delta2[3];

	delta[0] = gyro[0] - est->u1_mean[0];
	delta[1] = gyro[1] - est->u1_mean[1];
	delta[2] = gyro[2] - est->u1_mean[2];

	est->u1_mean[0] += delta[0] / n;
	est->u1_mean[1] += delta[1] / n;
	est->u1_mean[2] += delta[2] / n;

	delta2[0] = gyro[0] - est->u1_mean[0];
	delta2[1] = gyro[1] - est->u1_mean[1];
	delta2[2] = gyro[2] - est->u1_mean[2];
	// Compute Covarian Matrix of gyro to see which direction gyro distribute the most
	// Welford update
	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			est->u1_M2[i][j] += delta[i] * delta2[j];
		}
	}
}

uint8_t ComputeJointAxis1(Estimator *est)
{
	if (est->u1_sample_count < U1_MIN_SAMPLES) return 0U;

	float inv = 1.0f / (float)(est->u1_sample_count - 1U);
	float C[3][3];

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			C[i][j] = est->u1_M2[i][j] * inv;
		}
	}

	float trace = C[0][0] + C[1][1] + C[2][2];

	if (trace < U1_EPS)
	{
		est->u1_valid = 0U;
		return 0U;
	}

	float v[3];

	if (est->u1_valid)
	{
		v[0] = est->u1[0];
		v[1] = est->u1[1];
		v[2] = est->u1[2];
	}
	else
	{
		if ((C[0][0] >= C[1][1]) && (C[0][0] >= C[2][2]))
		{
			v[0] = 1.0f;
			v[1] = 0.0f;
			v[2] = 0.0f;
		}
		else if (C[1][1] >= C[2][2])
		{
			v[0] = 0.0f;
			v[1] = 1.0f;
			v[2] = 0.0f;
		}
		else
		{
			v[0] = 0.0f;
			v[1] = 0.0f;
			v[2] = 1.0f;
		}
	}

	for (uint8_t iteration = 0U; iteration < U1_POWER_ITERATIONS; iteration++)
	{
		float temp[3];

		Mat3Vec(C, v, temp);

		if (!Normalize3(temp))
		{
			est->u1_valid = 0U;
			return 0U;
		}

		v[0] = temp[0];
		v[1] = temp[1];
		v[2] = temp[2];
	}

	if (est->u1_valid)
	{
		if (Dot3(v, est->u1) < 0.0f)
		{
			v[0] = -v[0];
			v[1] = -v[1];
			v[2] = -v[2];
		}
	}

	float Cv[3];

	Mat3Vec(C, v, Cv);

	float lambda1 = Dot3(v, Cv);

	est->u1_confidence = lambda1 / trace;

	est->u1[0] = v[0];
	est->u1[1] = v[1];
	est->u1[2] = v[2];

	est->u1_valid = (est->u1_confidence >= U1_CONFIDENCE_MIN);

	return est->u1_valid;
}

void ResetJointAxis2(Estimator *est)
{
	est->u2_sample_count = 0U;

	est->u2_mean[0] = 0.0f;
	est->u2_mean[1] = 0.0f;
	est->u2_mean[2] = 0.0f;

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			est->u2_M2[i][j] = 0.0f;
		}
	}

	est->u2[0] = 0.99938f;
	est->u2[1] = -0.03507f;
	est->u2[2] = -0.00127f;

	est->u2_confidence = 0.0f;
	est->u2_motion = 0.0f;

	est->u2_valid = 1U;
	est->u2_locked = 1U;
	est->u2_compute_counter = 0U;
}
void UpdateJointAxis2(Estimator *est, const float *gyro1, const float *gyro2)
{
	if (est->u2_locked) return;
	
	float gyro1_mag = Norm3(gyro1);
	float gyro2_mag = Norm3(gyro2);

	if (gyro1_mag > U2_Q1_STILL_RAD_S) return;
	if (gyro2_mag < U2_MIN_GYRO_RAD_S) return;

	est->u2_sample_count++;

	float n = (float)est->u2_sample_count;

	float delta[3];
	float delta2[3];

	delta[0] = gyro2[0] - est->u2_mean[0];
	delta[1] = gyro2[1] - est->u2_mean[1];
	delta[2] = gyro2[2] - est->u2_mean[2];

	est->u2_mean[0] += delta[0] / n;
	est->u2_mean[1] += delta[1] / n;
	est->u2_mean[2] += delta[2] / n;

	delta2[0] = gyro2[0] - est->u2_mean[0];
	delta2[1] = gyro2[1] - est->u2_mean[1];
	delta2[2] = gyro2[2] - est->u2_mean[2];

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			est->u2_M2[i][j] += delta[i] * delta2[j];
		}
	}

	est->u2_motion += gyro2_mag * est->invSampleFreq;
}
uint8_t ComputeJointAxis2(Estimator *est)
{
	if (est->u2_sample_count < U2_MIN_SAMPLES) return 0U;

	float inv = 1.0f / (float)(est->u2_sample_count - 1U);
	float C[3][3];

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			C[i][j] = est->u2_M2[i][j] * inv;
		}
	}

	float trace = C[0][0] + C[1][1] + C[2][2];

	if (trace < U1_EPS)
	{
		est->u2_valid = 0U;
		return 0U;
	}

	float axis_new[3];
	float lambda = 0.0f;

	if (!PrincipalAxis3(C, axis_new, &lambda))
	{
		est->u2_valid = 0U;
		return 0U;
	}

	if (est->u2_valid)
	{
		if (Dot3(axis_new, est->u2) < 0.0f)
		{
			axis_new[0] = -axis_new[0];
			axis_new[1] = -axis_new[1];
			axis_new[2] = -axis_new[2];
		}
	}

	est->u2[0] = axis_new[0];
	est->u2[1] = axis_new[1];
	est->u2[2] = axis_new[2];

	est->u2_confidence = lambda / trace;

	est->u2_valid = (est->u2_confidence >= U2_CONFIDENCE_MIN) && (est->u2_motion >= U2_MIN_MOTION_RAD);

	return est->u2_valid;	
}
void UpdateHeadingData(Estimator *est, const float *gyro1, const float *gyro2, const float R1[3][3], const float R2[3][3])
{
	est->heading_sample_div_counter++;

	if (est->heading_sample_div_counter < HEADING_SAMPLE_DIV) return;

	est->heading_sample_div_counter = 0U;

	float omega1_E1[3];

	Mat3Vec(R1, gyro1, omega1_E1);

	float vecA[3];
	float vecB[3];
	float vecC[3];

	vecA[0] = omega1_E1[0];
	vecA[1] = omega1_E1[1];
	vecA[2] = 0.0f;

	vecB[0] = omega1_E1[1];
	vecB[1] = -omega1_E1[0];
	vecB[2] = 0.0f;

	vecC[0] = 0.0f;
	vecC[1] = 0.0f;
	vecC[2] = omega1_E1[2];

	float A[3];
	float B[3];
	float C[3];

	Mat3TVec(R2, vecA, A);
	Mat3TVec(R2, vecB, B);
	Mat3TVec(R2, vecC, C);

	uint16_t k = est->heading_buffer_head;

	est->heading_A[k][0] = A[0];
	est->heading_A[k][1] = A[1];
	est->heading_A[k][2] = A[2];

	est->heading_B[k][0] = B[0];
	est->heading_B[k][1] = B[1];
	est->heading_B[k][2] = B[2];

	est->heading_D[k][0] = gyro2[0] - C[0];
	est->heading_D[k][1] = gyro2[1] - C[1];
	est->heading_D[k][2] = gyro2[2] - C[2];

	est->heading_buffer_head++;

	if (est->heading_buffer_head >= HEADING_BUFFER_SIZE) est->heading_buffer_head = 0U;
	if (est->heading_buffer_count < HEADING_BUFFER_SIZE) est->heading_buffer_count++;
}
static uint8_t EvaluateHeadingFixedAxis(const Estimator *est, float delta, float *cost, float *excitation)
{
	if (!est->u2_locked) return 0U;
	if (est->heading_buffer_count < 20U) return 0U;

	float c = cosf(delta);
	float s = sinf(delta);

	float cost_sum = 0.0f;
	float excitation_sum = 0.0f;

	uint16_t start = (uint16_t)((est->heading_buffer_head + HEADING_BUFFER_SIZE - est->heading_buffer_count) % HEADING_BUFFER_SIZE);

	for (uint16_t n = 0U; n < est->heading_buffer_count; n++)
	{
		uint16_t k = (uint16_t)((start + n) % HEADING_BUFFER_SIZE);

		float v[3];

		v[0] = est->heading_D[k][0] - c*est->heading_A[k][0] - s*est->heading_B[k][0];
		v[1] = est->heading_D[k][1] - c*est->heading_A[k][1] - s*est->heading_B[k][1];
		v[2] = est->heading_D[k][2] - c*est->heading_A[k][2] - s*est->heading_B[k][2];

		float parallel = Dot3(v, est->u2);

		float rx = v[0] - parallel*est->u2[0];
		float ry = v[1] - parallel*est->u2[1];
		float rz = v[2] - parallel*est->u2[2];

		cost_sum += rx*rx + ry*ry + rz*rz;

		excitation_sum += est->heading_A[k][0]*est->heading_A[k][0] + est->heading_A[k][1]*est->heading_A[k][1] + est->heading_A[k][2]*est->heading_A[k][2];
	}

	*cost = cost_sum / (float)est->heading_buffer_count;
	*excitation = excitation_sum / (float)est->heading_buffer_count;

	return 1U;
}
uint8_t EstimateInitialHeading(Estimator *est)
{
	if (!est->u2_locked) return 0U;

	float best_delta = 0.0f;
	float best_cost = 1e30f;
	float best_excitation = 0.0f;

	for (float delta = -M_PI; delta <= M_PI; delta += HEADING_INIT_COARSE_RAD)
	{
		float cost;
		float excitation;

		if (!EvaluateHeadingFixedAxis(est, delta, &cost, &excitation)) continue;

		if (cost < best_cost)
		{
			best_cost = cost;
			best_delta = delta;
			best_excitation = excitation;
		}
	}

	float start = best_delta - HEADING_INIT_COARSE_RAD;
	float end = best_delta + HEADING_INIT_COARSE_RAD;

	for (float delta = start; delta <= end; delta += HEADING_INIT_FINE_RAD)
	{
		float test_delta = WrapPi(delta);

		float cost;
		float excitation;

		if (!EvaluateHeadingFixedAxis(est, test_delta, &cost, &excitation)) continue;

		if (cost < best_cost)
		{
			best_cost = cost;
			best_delta = test_delta;
			best_excitation = excitation;
		}
	}
	if (best_excitation < 0.1f)
	{
			return 0U;
	}
	est->heading_delta = WrapPi(best_delta);
	est->heading_cost = best_cost;
	est->heading_excitation = best_excitation;
	est->heading_valid = 1U;

	return 1U;
}
uint8_t TrackHeadingDelta(Estimator *est)
{
    if (!est->heading_valid) return 0U;
    if (!est->u2_locked) return 0U;

    float best_delta = est->heading_delta;
    float best_cost = 1e30f;
    float best_excitation = 0.0f;

    for (float offset = -HEADING_TRACK_RANGE_RAD; offset <= HEADING_TRACK_RANGE_RAD; offset += HEADING_TRACK_STEP_RAD)
    {
        float test_delta = WrapPi(est->heading_delta + offset);

        float cost;
        float excitation;

        if (!EvaluateHeadingFixedAxis(est, test_delta, &cost, &excitation)) continue;

        if (cost < best_cost)
        {
            best_cost = cost;
            best_delta = test_delta;
            best_excitation = excitation;
        }
    }

    est->heading_cost = best_cost;
    est->heading_excitation = best_excitation;

    if (best_excitation < 0.05f) return 0U;

    float delta_step = WrapPi(best_delta - est->heading_delta);
    float max_step = 0.008726646f;

    if (delta_step > max_step) delta_step = max_step;
    if (delta_step < -max_step) delta_step = -max_step;

    est->heading_delta = WrapPi(est->heading_delta + delta_step);

    return 1U;
}
void ComputeHeadingDiagnostic(Estimator *est, const float *gyro1, const float *gyro2, const float R1[3][3], const float R2[3][3])
{
	float omega1_E1[3];
	float omega1_E2[3];
	float omega1_S2[3];

	Mat3Vec(R1, gyro1, omega1_E1);

	float c = cosf(est->heading_delta);
	float s = sinf(est->heading_delta);

	omega1_E2[0] =  c*omega1_E1[0] + s*omega1_E1[1];
	omega1_E2[1] = -s*omega1_E1[0] + c*omega1_E1[1];
	omega1_E2[2] =  omega1_E1[2];

	Mat3TVec(R2, omega1_E2, omega1_S2);

	float v[3];

	v[0] = gyro2[0] - omega1_S2[0];
	v[1] = gyro2[1] - omega1_S2[1];
	v[2] = gyro2[2] - omega1_S2[2];
	float w1 = Norm3(gyro1);
	float w2 = Norm3(gyro2);
	float wavg = 0.5f*(w1 + w2);

	est->debug.norm_v = Norm3(v);

	if (wavg > 0.3f)
	{
			est->debug.norm_v_ratio = est->debug.norm_v / wavg;
	}
	est->heading_parallel = Dot3(v, est->u2);

	float perp[3];

	perp[0] = v[0] - est->heading_parallel*est->u2[0];
	perp[1] = v[1] - est->heading_parallel*est->u2[1];
	perp[2] = v[2] - est->heading_parallel*est->u2[2];

	est->heading_rel_norm = Norm3(v);
	est->heading_perp_norm = Norm3(perp);
}
void EstimatorDebugUpdate(Estimator *est, VQF6D *vqf1, const float *gyro1, const float *gyro2, const float R1[3][3], const float R2[3][3])
{
	float gyro1_mag = Norm3(gyro1);
	float gyro2_mag = Norm3(gyro2);

	est->debug.gyro1_norm = gyro1_mag;
	est->debug.gyro2_norm = gyro2_mag;

	if (gyro1_mag > 0.3f)
	{
		float dot_axis = fabsf(Dot3(gyro1, est->u1)) / gyro1_mag;

		if (dot_axis > 1.0f) dot_axis = 1.0f;

		est->debug.u1_axis_error_deg = acosf(dot_axis) * 57.2957795f;

		float parallel = Dot3(gyro1, est->u1);

		float perp[3];

		perp[0] = gyro1[0] - parallel*est->u1[0];
		perp[1] = gyro1[1] - parallel*est->u1[1];
		perp[2] = gyro1[2] - parallel*est->u1[2];

		est->debug.u1_perp_norm = Norm3(perp);
		est->debug.u1_perp_ratio = est->debug.u1_perp_norm / gyro1_mag;
	}

	Mat3Vec(R1, gyro1, est->debug.omega1_E1);
	Mat3Vec(R2, gyro2, est->debug.omega2_E2);

	est->debug.omega1_h = sqrtf(est->debug.omega1_E1[0]*est->debug.omega1_E1[0] + est->debug.omega1_E1[1]*est->debug.omega1_E1[1]);
	est->debug.omega2_h = sqrtf(est->debug.omega2_E2[0]*est->debug.omega2_E2[0] + est->debug.omega2_E2[1]*est->debug.omega2_E2[1]);

	if ((est->debug.omega1_h > 0.3f) && (est->debug.omega2_h > 0.3f))
	{
		est->debug.heading_direct = atan2f(est->debug.omega2_E2[0]*est->debug.omega1_E1[1] - est->debug.omega2_E2[1]*est->debug.omega1_E1[0], est->debug.omega2_E2[0]*est->debug.omega1_E1[0] + est->debug.omega2_E2[1]*est->debug.omega1_E1[1]);
	}

	float q3D[4];
	float R3D[3][3];

	VQF6D_GetQuat3D(vqf1, q3D);
	VQF6D_QuatToRotMat(q3D, R3D);

	Mat3Vec(R3D, est->u1, est->debug.j1_3D);
	Mat3Vec(R1, est->u1, est->debug.j1_6D);

	est->debug.j1_3D_h = sqrtf(est->debug.j1_3D[0]*est->debug.j1_3D[0] + est->debug.j1_3D[1]*est->debug.j1_3D[1]);
	est->debug.j1_3D_z = est->debug.j1_3D[2];

	est->debug.j1_6D_h = sqrtf(est->debug.j1_6D[0]*est->debug.j1_6D[0] + est->debug.j1_6D[1]*est->debug.j1_6D[1]);
	est->debug.j1_6D_z = est->debug.j1_6D[2];
	
	float D1[3][3];
	float D2[3][3];

	Mat3MulTransB(R1, est->debug.R1_ref, D1);
	Mat3MulTransB(R2, est->debug.R2_ref, D2);

	est->debug.drift1 = atan2f(D1[1][0], D1[0][0]);
	est->debug.drift2 = atan2f(D2[1][0], D2[0][0]);

	est->debug.relative_drift = WrapPi(est->debug.drift2 - est->debug.drift1);
}