//=============================================================================================
// ESTIMATE.c

#include "Estimate.h"
#include <math.h>
#include "utilzfunc.h"

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

static float AccMagWeight(float acc_norm)
{
    const float G = 1.0f;
    float rel_err = fabsf(acc_norm - G);
    if (rel_err <= 0.02f) return 1.0f;
    if (rel_err >= 0.15f)return 0.0f;
    return 1.0f - (rel_err - 0.02f)/ (0.15f - 0.02f);
}
static float MotionWeight(float omega, float omega_ref)
{
    float x = omega / omega_ref;
    float x2 = x*x;
    return 1.0f / (1.0f + x2*x2);
}
static float ObsWeightQ2(float h)
{
		float w_obs = 0.0f;
		if(h<=0.45f) w_obs = 0.0f;
		else if(h >= 0.65f) w_obs = 1.0f;
		else
		{
			float x = (h-0.45f)/0.2f;
			w_obs = x*x*(3.0f - 2.0f*x);
		}
		return w_obs;
}

void InitEstimator(Estimator *est)
{
	est->invSampleFreq = 0.01f;
	est->q1_pred = 0.0f, est->q2_pred = 0.0f;
	est->q1dot_pred = 0.0f, est->q2dot_pred = 0.0f;
	
	est->j1[0] = 0.0558833f;
	est->j1[1] = 0.9918949f;
	est->j1[2] = 0.1143946f;
	
	est->j2[0] = 0.999426f;
	est->j2[1] = -0.032874f;
	est->j2[2] = -0.0081791f;
	
	est->j12[0] = 0.0769084f;
	est->j12[1] = 0.9948321f;
	est->j12[2] = 0.06628749f;
	
	est->q2_bias = 0.0003f;
	est->P = 0.0027f;
	est->Q = Q2KF_SIGMA_RATE * Q2KF_SIGMA_RATE * est->invSampleFreq *est->invSampleFreq;
	
	est->q2_R = 0.0f;
	est->q2_Kq = 0.0f;
	est->h2_norm = 0.0f;

	est->r2[0] = 0.20746f;
	est->r2[1] = 0.0045f;
	est->r2[2] = 0.11862f;
	est->c12 = Dot3(est->j2, est->j12);
	est->lambda_q1 = 0.005f; est->lambda_q2 = 0.005f;
	est->Kq1 = 0.05f;	est->w1 = 0.0f;
	est->Kq2 = 0.01f; est->w2 = 0.0f;
	for (uint8_t i = 0U; i < 3U; i++)
	{
			est->gyro2_prev[i] = 0.0f;
			est->alpha2_f[i] = 0.0f;
			est->a_dyn2[i] = 0.0f;
			est->a_tan2[i] = 0.0f;
			est->a_cent2[i] = 0.0f;
	}

	est->dyn_comp_gain = 1.0f;
	est->dyn_comp_initialized = 0U;
}

void CompensateLinearAcc(Estimator *est, float *gyro2, float *gValue2)
{
	float dt = 0.01f;
	if (!est->dyn_comp_initialized)
	{
			for (uint8_t i = 0U; i < 3U; i++)
			{
					est->gyro2_prev[i] = gyro2[i];
					est->alpha2_f[i] = 0.0f;
			}

			est->dyn_comp_initialized = 1U;
	}
	float tau = 1.0f / (2.0f * M_PI * ALPHA2_FC_HZ);
	float k = dt / (tau + dt);
	for (uint8_t i = 0U; i < 3U; i++)
	{
			float alpha_raw = (gyro2[i] - est->gyro2_prev[i]) / dt;
			est->alpha2_f[i] += k * (alpha_raw - est->alpha2_f[i]);
			est->gyro2_prev[i] = gyro2[i];
	}
	/* ---------------------------------------------
	 * Tangential:
	 *
	 * a_tan = alpha x r
	 * --------------------------------------------- */
	Cross3(est->alpha2_f,est->r2, est->a_tan2);
	/* ---------------------------------------------
	 * Centripetal:
	 *
	 * a_cent = omega x (omega x r)
	 * --------------------------------------------- */
	float omega_cross_r[3];
	Cross3(gyro2, est->r2, omega_cross_r);
	Cross3(gyro2,omega_cross_r,est->a_cent2);
	for (uint8_t i = 0U; i < 3U; i++)
	{
			est->a_dyn2[i] = est->a_tan2[i] + est->a_cent2[i];
			est->acc_comp[i] = gValue2[i] - est->dyn_comp_gain * est->a_dyn2[i]/9.80665f;
	}
}

void ComputeAngles(Estimator *est, float *gyro1, float *gyro2, float *gValue1, float *gValue2)
{
	CompensateLinearAcc(est,gyro2, gValue2);
	float g1m[3] = {gValue1[0], gValue1[1], gValue1[2]};
	float g2m[3] = {est->acc_comp[0], est->acc_comp[1], est->acc_comp[2]};
	est->a1_norm = Norm3(gValue1);
	est->a2_norm = Norm3(est->acc_comp);
	est->gyro1_norm = Norm3(gyro1);
	est->gyro2_norm = Norm3(gyro2);
	uint8_t stationary =
    (est->gyro1_norm < 0.1f) &&
    (est->gyro2_norm < 0.1f) &&
    (fabsf(est->a1_norm - 1.0f) < 0.02f) &&
    (fabsf(est->a2_norm - 1.0f) < 0.02f);
	// Prediction phase
	est->q1dot_pred = Dot3(est->j1,gyro1);
	est->q1_pred += est->q1dot_pred*est->invSampleFreq;
	
	est->q2dot_pred = (stationary) ? 0.0f : Dot3(est->j2, gyro2)- est->c12 * est->q1dot_pred;
	est->q2_pred += 	est->q2dot_pred * est->invSampleFreq;
	
	Normalize3(g1m);
	Normalize3(g2m);
	
	// q1 correction
	est->q1_acc = -AngleAroundAxis(est->g10,g1m,est->j1);
	est->ePred_g1 = WrapPi(est->q1_acc - est->q1_pred);
	float wmag1 = AccMagWeight(est->a1_norm);
	float omega1 = Norm3(gyro1);
	float wmotion1 = MotionWeight(omega1, 0.60f);
	est->w1 = wmag1 * wmotion1;
	est->ePred_g1 = ClampF(est->ePred_g1,-10.0f*DEG2RAD,10.0f*DEG2RAD);
	est->q1_pred += est->Kq1 * est->w1 * est->ePred_g1;
	
	//q2 correction
	RotateAxis(est->j12,-est->q1_pred,est->g20, est->g2_ref);
	Normalize3(est->g2_ref);
	float h2[3];
	Cross3(est->j2, est->g2_ref, h2);
	h2[0] = -h2[0];
	h2[1] = -h2[1];
	h2[2] = -h2[2];
	est->h2_norm = Norm3(h2);	// Check observability
//	float w_unobs = 1.0f - est->h2_norm*est->h2_norm;
//	float wobs = ObsWeightQ2(est->h2_norm);		
//	if (est->h2_norm < 0.20f)
//	{
//			est->P += est->Q + 3e-7f;
//	}
//	else if (est->h2_norm < 0.65f)
//	{
//			est->P += est->Q+ 3e-7f * (1.0f - wobs);
//	}
//	else
//	{
//			est->P += est->Q;
//	}
	est->q2_acc = -AngleAroundAxis(est->g2_ref, g2m, est->j2);
	est->ePred_g2 = WrapPi(est->q2_acc - est->q2_pred);

	float wmag2 = AccMagWeight(est->a2_norm);
	float omega2 = Norm3(gyro2);
	float wmotion2 = MotionWeight(omega2, 0.5f);
	float wobs2 = ObsWeightQ2(est->h2_norm);
	float quality = wobs2*wmag2 * wmotion2;
	
	est->q2_pred +=	0.01f*quality * est->ePred_g2;
	
//	if(wobs > 0.0f & quality > 0.05f)
//	{
//		float sigma_meas =  Q2KF_SIGMA_ACC / (wobs*quality);
//		est->q2_R = sigma_meas * sigma_meas;;
//		
//		float S = est->P + est->q2_R;
//		if (S > 1e-12f)
//		{
//			est->q2_Kq = est->P/S;
//			est->q2_pred +=	est->q2_Kq * est->ePred_g2;
//			est->dq2_corr =
//			est->q2_Kq *
//			est->ePred_g2;
//			est->P = (1.0f - est->q2_Kq)*est->P;
//			if (est->P < 1e-12f) est->P = 1e-12f;
//		}
//	}
	est->bending_angle = -est->q1_pred*RAD2DEG;
	est->twist_angle = -est->q2_pred*RAD2DEG;

}

void InitCalibProcess(CalibProcess *calib)
{
	calib->dt = 0.01f;
	ResetJointAxis1(calib);
	ResetJointAxis2(calib);
	ResetJ12Axis(calib);
}

void ResetJ12Axis(CalibProcess *calib)
{
	calib->j12_num[0] = 0.0f; calib->j12_num[1] = 0.0f; calib->j12_num[2] = 0.0f;
	calib->j12_den = 0.0f;
	calib->j12_count = 0;
	calib->j12_gyro2_energy = 0.0f;
	calib->j12_confidence = 0.0f;
	calib->j12_scale = 0.0f;
	calib->j12_motion = 0.0f;
	calib->j12_valid = 0U;
}
void ComputeJ12Axis(CalibProcess *calib, float *gyro1, float *gyro2)
{
	float q1dot_pred = Dot3(calib->u1,gyro1);
	if(fabsf(q1dot_pred) > 0.3f)
	{
		calib->j12_num[0] += q1dot_pred * gyro2[0];
    calib->j12_num[1] += q1dot_pred* gyro2[1];
    calib->j12_num[2] += q1dot_pred * gyro2[2];

    calib->j12_den += q1dot_pred * q1dot_pred;
		calib->j12_gyro2_energy += Dot3(gyro2, gyro2);
		calib->j12_motion +=  fabsf(q1dot_pred) * calib->dt;
    calib->j12_count++;
	}
	if(calib->j12_den > 1e-6f && calib->j12_gyro2_energy > 1e-9f)
	{
			calib->j12[0] = calib->j12_num[0] / calib->j12_den;
			calib->j12[1] = calib->j12_num[1] / calib->j12_den;
			calib->j12[2] = calib->j12_num[2] / calib->j12_den;
			Normalize3(calib->j12);
			float num_energy = Dot3(calib->j12_num,calib->j12_num); 
			float explained = num_energy / calib->j12_den;
			calib->j12_confidence = explained / calib->j12_gyro2_energy;
			calib->j12_confidence = ClampF(calib->j12_confidence,  0.0f, 1.0f);
	}
}

void ResetJointAxis1(CalibProcess *calib)
{
	calib->u1_sample_count = 0U;

	calib->u1_mean[0] = 0.0f;
	calib->u1_mean[1] = 0.0f;
	calib->u1_mean[2] = 0.0f;

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			calib->u1_M2[i][j] = 0.0f;
		}
	}

	calib->u1[0] = 0.0558833f;
	calib->u1[1] = 0.9918949f;
	calib->u1[2] = 0.1143946f;

	calib->u1_confidence = 0.0f;
	calib->u1_valid = 0U;
	calib->u1_locked = 1U;
	calib->u1_compute_counter = 0U;
}
void UpdateJointAxis1(CalibProcess *calib, const float *gyro)
{
	// Builf Welford covariance
	if (calib->u1_locked) return;

	float gyro_mag = Norm3(gyro);

	if (gyro_mag < U1_MIN_GYRO_RAD_S) return;
	
	// remove gyro bias
	calib->u1_sample_count++;

	float n = (float)calib->u1_sample_count;

	float delta[3];
	float delta2[3];

	delta[0] = gyro[0] - calib->u1_mean[0];
	delta[1] = gyro[1] - calib->u1_mean[1];
	delta[2] = gyro[2] - calib->u1_mean[2];

	calib->u1_mean[0] += delta[0] / n;
	calib->u1_mean[1] += delta[1] / n;
	calib->u1_mean[2] += delta[2] / n;

	delta2[0] = gyro[0] - calib->u1_mean[0];
	delta2[1] = gyro[1] - calib->u1_mean[1];
	delta2[2] = gyro[2] - calib->u1_mean[2];
	// Compute Covarian Matrix of gyro to see which direction gyro distribute the most
	// Welford update
	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			calib->u1_M2[i][j] += delta[i] * delta2[j];
		}
	}
}

uint8_t ComputeJointAxis1(CalibProcess *calib)
{
	if (calib->u1_sample_count < U1_MIN_SAMPLES) return 0U;

	float inv = 1.0f / (float)(calib->u1_sample_count - 1U);
	float C[3][3];

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			C[i][j] = calib->u1_M2[i][j] * inv;
		}
	}

	float trace = C[0][0] + C[1][1] + C[2][2];

	if (trace < U1_EPS)
	{
		calib->u1_valid = 0U;
		return 0U;
	}

	float v[3];

	if (calib->u1_valid)
	{
		v[0] = calib->u1[0];
		v[1] = calib->u1[1];
		v[2] = calib->u1[2];
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
			calib->u1_valid = 0U;
			return 0U;
		}

		v[0] = temp[0];
		v[1] = temp[1];
		v[2] = temp[2];
	}

	if (calib->u1_valid)
	{
		if (Dot3(v, calib->u1) < 0.0f)
		{
			v[0] = -v[0];
			v[1] = -v[1];
			v[2] = -v[2];
		}
	}

	float Cv[3];

	Mat3Vec(C, v, Cv);

	float lambda1 = Dot3(v, Cv);

	calib->u1_confidence = lambda1 / trace;

	calib->u1[0] = v[0];
	calib->u1[1] = v[1];
	calib->u1[2] = v[2];

	calib->u1_valid = (calib->u1_confidence >= U1_CONFIDENCE_MIN);

	return calib->u1_valid;
}

void ResetJointAxis2(CalibProcess *calib)
{
	calib->u2_sample_count = 0U;

	calib->u2_mean[0] = 0.0f;
	calib->u2_mean[1] = 0.0f;
	calib->u2_mean[2] = 0.0f;

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			calib->u2_M2[i][j] = 0.0f;
		}
	}

	calib->u2[0] = 0.999426f;
	calib->u2[1] = -0.032874f;
	calib->u2[2] = -0.0081791f;

	calib->u2_confidence = 0.0f;
	calib->u2_motion = 0.0f;

	calib->u2_valid = 0U;
	calib->u2_locked = 1U;
	calib->u2_compute_counter = 0U;
}

void UpdateJointAxis2(CalibProcess *calib, const float *gyro1, const float *gyro2)
{
	if (calib->u2_locked) return;
	
	float gyro1_mag = Norm3(gyro1);
	float gyro2_mag = Norm3(gyro2);

	if (gyro1_mag > U2_Q1_STILL_RAD_S) return;
	if (gyro2_mag < U2_MIN_GYRO_RAD_S) return;

	calib->u2_sample_count++;

	float n = (float)calib->u2_sample_count;

	float delta[3];
	float delta2[3];

	delta[0] = gyro2[0] - calib->u2_mean[0];
	delta[1] = gyro2[1] - calib->u2_mean[1];
	delta[2] = gyro2[2] - calib->u2_mean[2];

	calib->u2_mean[0] += delta[0] / n;
	calib->u2_mean[1] += delta[1] / n;
	calib->u2_mean[2] += delta[2] / n;

	delta2[0] = gyro2[0] - calib->u2_mean[0];
	delta2[1] = gyro2[1] - calib->u2_mean[1];
	delta2[2] = gyro2[2] - calib->u2_mean[2];

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			calib->u2_M2[i][j] += delta[i] * delta2[j];
		}
	}

	calib->u2_motion += gyro2_mag * calib->dt;
}
uint8_t ComputeJointAxis2(CalibProcess *calib)
{
	if (calib->u2_sample_count < U2_MIN_SAMPLES) return 0U;

	float inv = 1.0f / (float)(calib->u2_sample_count - 1U);
	float C[3][3];

	for (uint8_t i = 0U; i < 3U; i++)
	{
		for (uint8_t j = 0U; j < 3U; j++)
		{
			C[i][j] = calib->u2_M2[i][j] * inv;
		}
	}

	float trace = C[0][0] + C[1][1] + C[2][2];

	if (trace < U1_EPS)
	{
		calib->u2_valid = 0U;
		return 0U;
	}

	float axis_new[3];
	float lambda = 0.0f;

	if (!PrincipalAxis3(C, axis_new, &lambda))
	{
		calib->u2_valid = 0U;
		return 0U;
	}

	if (calib->u2_valid)
	{
		if (Dot3(axis_new, calib->u2) < 0.0f)
		{
			axis_new[0] = -axis_new[0];
			axis_new[1] = -axis_new[1];
			axis_new[2] = -axis_new[2];
		}
	}

	calib->u2[0] = axis_new[0];
	calib->u2[1] = axis_new[1];
	calib->u2[2] = axis_new[2];

	calib->u2_confidence = lambda / trace;

	calib->u2_valid = (calib->u2_confidence >= U2_CONFIDENCE_MIN) && (calib->u2_motion >= U2_MIN_MOTION_RAD);

	return calib->u2_valid;	
}