/*
 * SPDX-FileCopyrightText: 2021 Daniel Laidig <laidig@control.tu-berlin.de>
 * SPDX-License-Identifier: MIT
 *
 * Embedded-C adaptation of the 6D path of VQF (gyro + accelerometer only).
 * The algorithm and default parameters follow the current VQF MATLAB implementation.
 */

#include "vqf6d.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#define VQF6D_PI 3.14159265358979323846f
#define VQF6D_SQRT2 1.41421356237309504880f
#define VQF6D_EPS FLT_EPSILON

static float vqf6d_clip_scalar(float x, float lo, float hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static float vqf6d_norm(const float *v, size_t n)
{
    float s = 0.0f;
    size_t i;
    for (i = 0; i < n; ++i) s += v[i] * v[i];
    return sqrtf(s);
}

static void vqf6d_normalize(float *v, size_t n)
{
    float nrm = vqf6d_norm(v, n);
    size_t i;
    if (nrm <= VQF6D_EPS) return;
    for (i = 0; i < n; ++i) v[i] /= nrm;
}

static void vqf6d_quat_multiply(const float q1[4], const float q2[4], float out[4])
{
    float w = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
    float x = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
    float y = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1];
    float z = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0];
    out[0] = w;
    out[1] = x;
    out[2] = y;
    out[3] = z;
}

void VQF6D_RotateVector(const float q[4], const float v[3], float out[3])
{
    out[0] = (1.0f - 2.0f*q[2]*q[2] - 2.0f*q[3]*q[3])*v[0] + 2.0f*v[1]*(q[2]*q[1] - q[0]*q[3]) + 2.0f*v[2]*(q[0]*q[2] + q[3]*q[1]);
    out[1] = 2.0f*v[0]*(q[0]*q[3] + q[2]*q[1]) + (1.0f - 2.0f*q[1]*q[1] - 2.0f*q[3]*q[3])*v[1] + 2.0f*v[2]*(q[2]*q[3] - q[1]*q[0]);
    out[2] = 2.0f*v[0]*(q[3]*q[1] - q[0]*q[2]) + 2.0f*v[1]*(q[0]*q[1] + q[3]*q[2]) + (1.0f - 2.0f*q[1]*q[1] - 2.0f*q[2]*q[2])*v[2];
}

void VQF6D_QuatToRotMat(const float q[4], float R[3][3])
{
    R[0][0] = 1.0f - 2.0f*q[2]*q[2] - 2.0f*q[3]*q[3];
    R[0][1] = 2.0f*(q[2]*q[1] - q[0]*q[3]);
    R[0][2] = 2.0f*(q[0]*q[2] + q[3]*q[1]);
    R[1][0] = 2.0f*(q[0]*q[3] + q[2]*q[1]);
    R[1][1] = 1.0f - 2.0f*q[1]*q[1] - 2.0f*q[3]*q[3];
    R[1][2] = 2.0f*(q[2]*q[3] - q[1]*q[0]);
    R[2][0] = 2.0f*(q[3]*q[1] - q[0]*q[2]);
    R[2][1] = 2.0f*(q[0]*q[1] + q[3]*q[2]);
    R[2][2] = 1.0f - 2.0f*q[1]*q[1] - 2.0f*q[2]*q[2];
}

static void vqf6d_filter_coeffs(float tau, float Ts, float b[3], float a[2])
{
    if (tau < Ts * 0.5f) {
        b[0] = 1.0f;
        b[1] = 0.0f;
        b[2] = 0.0f;
        a[0] = 0.0f;
        a[1] = 0.0f;
        return;
    }

    {
        float fc = VQF6D_SQRT2 / (2.0f * VQF6D_PI * tau);
        float C = tanf(VQF6D_PI * fc * Ts);
        float D = C*C + VQF6D_SQRT2*C + 1.0f;
        float b0 = C*C / D;
        b[0] = b0;
        b[1] = 2.0f*b0;
        b[2] = b0;
        a[0] = 2.0f*(C*C - 1.0f) / D;
        a[1] = (1.0f - VQF6D_SQRT2*C + C*C) / D;
    }
}

static void vqf6d_filter_reset(VQF6D_FilterState *state)
{
    memset(state, 0, sizeof(*state));
    state->initialized = false;
    state->count = 0U;
}

static void vqf6d_filter_initial_state(float x0, const float b[3], const float a[2], float *z1, float *z2)
{
    *z1 = x0 * (1.0f - b[0]);
    *z2 = x0 * (b[2] - a[1]);
}

static void vqf6d_filter_vec(const float *x, size_t n, float tau, float Ts, const float b[3], const float a[2], VQF6D_FilterState *state, float *out)
{
    size_t i;

    if (!state->initialized) {
        state->count++;
        for (i = 0; i < n; ++i) {
            state->sum[i] += x[i];
            out[i] = state->sum[i] / (float)state->count;
        }

        if ((float)state->count * Ts >= tau) {
            for (i = 0; i < n; ++i) vqf6d_filter_initial_state(out[i], b, a, &state->z1[i], &state->z2[i]);
            state->initialized = true;
        }
        return;
    }

    for (i = 0; i < n; ++i) {
        float y = b[0]*x[i] + state->z1[i];
        state->z1[i] = b[1]*x[i] - a[0]*y + state->z2[i];
        state->z2[i] = b[2]*x[i] - a[1]*y;
        out[i] = y;
    }
}

static bool vqf6d_inv3(const float A[3][3], float invA[3][3])
{
    float C00 = A[1][1]*A[2][2] - A[1][2]*A[2][1];
    float C01 = A[1][2]*A[2][0] - A[1][0]*A[2][2];
    float C02 = A[1][0]*A[2][1] - A[1][1]*A[2][0];
    float C10 = A[0][2]*A[2][1] - A[0][1]*A[2][2];
    float C11 = A[0][0]*A[2][2] - A[0][2]*A[2][0];
    float C12 = A[0][1]*A[2][0] - A[0][0]*A[2][1];
    float C20 = A[0][1]*A[1][2] - A[0][2]*A[1][1];
    float C21 = A[0][2]*A[1][0] - A[0][0]*A[1][2];
    float C22 = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    float det = A[0][0]*C00 + A[0][1]*C01 + A[0][2]*C02;

    if (fabsf(det) < 1e-12f) return false;

    {
        float id = 1.0f / det;
        invA[0][0] = C00*id; invA[0][1] = C10*id; invA[0][2] = C20*id;
        invA[1][0] = C01*id; invA[1][1] = C11*id; invA[1][2] = C21*id;
        invA[2][0] = C02*id; invA[2][1] = C12*id; invA[2][2] = C22*id;
    }
    return true;
}

static void vqf6d_mat3_mul(const float A[3][3], const float B[3][3], float C[3][3])
{
    size_t r, c, k;
    for (r = 0; r < 3; ++r) {
        for (c = 0; c < 3; ++c) {
            float sum = 0.0f;
            for (k = 0; k < 3; ++k) sum += A[r][k] * B[k][c];
            C[r][c] = sum;
        }
    }
}

static void vqf6d_mat3_transpose(const float A[3][3], float AT[3][3])
{
    AT[0][0] = A[0][0]; AT[0][1] = A[1][0]; AT[0][2] = A[2][0];
    AT[1][0] = A[0][1]; AT[1][1] = A[1][1]; AT[1][2] = A[2][1];
    AT[2][0] = A[0][2]; AT[2][1] = A[1][2]; AT[2][2] = A[2][2];
}

static void vqf6d_bias_update(VQF6D *vqf, const float accEarth[3], const float RlpIn[3][3], const float biasLpIn[2])
{
    float R[3][3];
    float e[3] = {0.0f, 0.0f, 0.0f};
    float w[3] = {-1.0f, -1.0f, -1.0f};
    float biasClip = vqf->params.biasClip * VQF6D_PI / 180.0f;
    float RT[3][3], PRt[3][3], RP[3][3], S[3][3], invS[3][3], K[3][3], KR[3][3], KRP[3][3];
    float newBias[3];
    size_t i;

    memcpy(R, RlpIn, sizeof(R));

    if (vqf->restDetected && vqf->params.restBiasEstEnabled) {
        e[0] = vqf->restLastGyrLp[0] - vqf->bias[0];
        e[1] = vqf->restLastGyrLp[1] - vqf->bias[1];
        e[2] = vqf->restLastGyrLp[2] - vqf->bias[2];
        memset(R, 0, sizeof(R));
        R[0][0] = R[1][1] = R[2][2] = 1.0f;
        w[0] = w[1] = w[2] = vqf->coeffs.biasRestW;
    } else if (vqf->params.motionBiasEstEnabled) {
        e[0] = -accEarth[1]/vqf->coeffs.accTs + biasLpIn[0] - (R[0][0]*vqf->bias[0] + R[0][1]*vqf->bias[1] + R[0][2]*vqf->bias[2]);
        e[1] =  accEarth[0]/vqf->coeffs.accTs + biasLpIn[1] - (R[1][0]*vqf->bias[0] + R[1][1]*vqf->bias[1] + R[1][2]*vqf->bias[2]);
        e[2] = -(R[2][0]*vqf->bias[0] + R[2][1]*vqf->bias[1] + R[2][2]*vqf->bias[2]);
        w[0] = vqf->coeffs.biasMotionW;
        w[1] = vqf->coeffs.biasMotionW;
        w[2] = vqf->coeffs.biasVerticalW;
    }

    if (vqf->biasP[0][0] < vqf->coeffs.biasP0) vqf->biasP[0][0] += vqf->coeffs.biasV;
    if (vqf->biasP[1][1] < vqf->coeffs.biasP0) vqf->biasP[1][1] += vqf->coeffs.biasV;
    if (vqf->biasP[2][2] < vqf->coeffs.biasP0) vqf->biasP[2][2] += vqf->coeffs.biasV;

    if (w[0] < 0.0f) return;

    for (i = 0; i < 3; ++i) e[i] = vqf6d_clip_scalar(e[i], -biasClip, biasClip);

    vqf6d_mat3_transpose(R, RT);
    vqf6d_mat3_mul(vqf->biasP, RT, PRt);
    vqf6d_mat3_mul(R, vqf->biasP, RP);
    vqf6d_mat3_mul(RP, RT, S);
    S[0][0] += w[0];
    S[1][1] += w[1];
    S[2][2] += w[2];

    if (!vqf6d_inv3(S, invS)) return;
    vqf6d_mat3_mul(PRt, invS, K);

    newBias[0] = vqf->bias[0] + K[0][0]*e[0] + K[0][1]*e[1] + K[0][2]*e[2];
    newBias[1] = vqf->bias[1] + K[1][0]*e[0] + K[1][1]*e[1] + K[1][2]*e[2];
    newBias[2] = vqf->bias[2] + K[2][0]*e[0] + K[2][1]*e[1] + K[2][2]*e[2];

    vqf6d_mat3_mul(K, R, KR);
    vqf6d_mat3_mul(KR, vqf->biasP, KRP);
    for (i = 0; i < 3; ++i) {
        size_t j;
        for (j = 0; j < 3; ++j) vqf->biasP[i][j] -= KRP[i][j];
    }

    /* Preserve symmetry against float round-off. */
    vqf->biasP[0][1] = vqf->biasP[1][0] = 0.5f*(vqf->biasP[0][1] + vqf->biasP[1][0]);
    vqf->biasP[0][2] = vqf->biasP[2][0] = 0.5f*(vqf->biasP[0][2] + vqf->biasP[2][0]);
    vqf->biasP[1][2] = vqf->biasP[2][1] = 0.5f*(vqf->biasP[1][2] + vqf->biasP[2][1]);

    vqf->bias[0] = vqf6d_clip_scalar(newBias[0], -biasClip, biasClip);
    vqf->bias[1] = vqf6d_clip_scalar(newBias[1], -biasClip, biasClip);
    vqf->bias[2] = vqf6d_clip_scalar(newBias[2], -biasClip, biasClip);
}

void VQF6D_GetDefaultParams(VQF6D_Params *params)
{
    if (params == NULL) return;
    params->tauAcc = 1.0f;
    params->motionBiasEstEnabled = true;
    params->restBiasEstEnabled = true;
    params->biasSigmaInit = 0.5f;
    params->biasForgettingTime = 100.0f;
    params->biasClip = 2.0f;
    params->biasSigmaMotion = 0.1f;
    params->biasVerticalForgettingFactor = 0.0001f;
    params->biasSigmaRest = 0.03f;
    params->restMinT = 1.5f;
    params->restFilterTau = 0.5f;
    params->restThGyr = 2.0f;
    params->restThAcc = 0.5f;
}

void VQF6D_Reset(VQF6D *vqf)
{
    if (vqf == NULL) return;

    vqf->gyrQuat[0]=1.0f; vqf->gyrQuat[1]=0.0f; vqf->gyrQuat[2]=0.0f; vqf->gyrQuat[3]=0.0f;
    vqf->accQuat[0]=1.0f; vqf->accQuat[1]=0.0f; vqf->accQuat[2]=0.0f; vqf->accQuat[3]=0.0f;

    memset(vqf->lastAccLp, 0, sizeof(vqf->lastAccLp));
    vqf6d_filter_reset(&vqf->accLpState);
    vqf->lastAccCorrAngularRate = 0.0f;

    memset(vqf->bias, 0, sizeof(vqf->bias));
    memset(vqf->biasP, 0, sizeof(vqf->biasP));
    vqf->biasP[0][0] = vqf->coeffs.biasP0;
    vqf->biasP[1][1] = vqf->coeffs.biasP0;
    vqf->biasP[2][2] = vqf->coeffs.biasP0;
    vqf6d_filter_reset(&vqf->motionBiasEstRLpState);
    vqf6d_filter_reset(&vqf->motionBiasEstBiasLpState);

    vqf->restDetected = false;
    vqf->restLastSquaredDeviations[0] = 0.0f;
    vqf->restLastSquaredDeviations[1] = 0.0f;
    vqf->restT = 0.0f;
    memset(vqf->restLastGyrLp, 0, sizeof(vqf->restLastGyrLp));
    memset(vqf->restLastAccLp, 0, sizeof(vqf->restLastAccLp));
    vqf6d_filter_reset(&vqf->restGyrLpState);
    vqf6d_filter_reset(&vqf->restAccLpState);
}

bool VQF6D_Init(VQF6D *vqf, float gyrTs, float accTs, const VQF6D_Params *params)
{
    float pMotion, pRest;

    if (vqf == NULL || gyrTs <= 0.0f) return false;
    if (accTs <= 0.0f) accTs = gyrTs;

    memset(vqf, 0, sizeof(*vqf));
    if (params != NULL) vqf->params = *params;
    else VQF6D_GetDefaultParams(&vqf->params);

    if (vqf->params.tauAcc <= 0.0f || vqf->params.restFilterTau <= 0.0f || vqf->params.biasForgettingTime <= 0.0f) return false;

    vqf->coeffs.gyrTs = gyrTs;
    vqf->coeffs.accTs = accTs;
    vqf6d_filter_coeffs(vqf->params.tauAcc, accTs, vqf->coeffs.accLpB, vqf->coeffs.accLpA);

    vqf->coeffs.biasP0 = (vqf->params.biasSigmaInit*100.0f)*(vqf->params.biasSigmaInit*100.0f);
    vqf->coeffs.biasV = (0.1f*100.0f)*(0.1f*100.0f)*accTs/vqf->params.biasForgettingTime;
    pMotion = (vqf->params.biasSigmaMotion*100.0f)*(vqf->params.biasSigmaMotion*100.0f);
    vqf->coeffs.biasMotionW = pMotion*pMotion/vqf->coeffs.biasV + pMotion;
    vqf->coeffs.biasVerticalW = vqf->coeffs.biasMotionW / fmaxf(vqf->params.biasVerticalForgettingFactor, 1e-10f);
    pRest = (vqf->params.biasSigmaRest*100.0f)*(vqf->params.biasSigmaRest*100.0f);
    vqf->coeffs.biasRestW = pRest*pRest/vqf->coeffs.biasV + pRest;

    vqf6d_filter_coeffs(vqf->params.restFilterTau, gyrTs, vqf->coeffs.restGyrLpB, vqf->coeffs.restGyrLpA);
    vqf6d_filter_coeffs(vqf->params.restFilterTau, accTs, vqf->coeffs.restAccLpB, vqf->coeffs.restAccLpA);

    VQF6D_Reset(vqf);
    return true;
}

void VQF6D_UpdateGyr(VQF6D *vqf, const float gyr[3])
{
    float gyrLp[3];
    float gyrNoBias[3];
    float gyrNorm;
    float angle;

    if (vqf == NULL || gyr == NULL) return;

    if (vqf->params.restBiasEstEnabled) {
        float dev[3];
        float squaredDeviation;
        float biasClip = vqf->params.biasClip * VQF6D_PI / 180.0f;
        float maxAbs;

        vqf6d_filter_vec(gyr, 3U, vqf->params.restFilterTau, vqf->coeffs.gyrTs, vqf->coeffs.restGyrLpB, vqf->coeffs.restGyrLpA, &vqf->restGyrLpState, gyrLp);
        dev[0]=gyr[0]-gyrLp[0]; dev[1]=gyr[1]-gyrLp[1]; dev[2]=gyr[2]-gyrLp[2];
        squaredDeviation = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
        maxAbs = fmaxf(fabsf(gyrLp[0]), fmaxf(fabsf(gyrLp[1]), fabsf(gyrLp[2])));
        if (squaredDeviation >= (vqf->params.restThGyr*VQF6D_PI/180.0f)*(vqf->params.restThGyr*VQF6D_PI/180.0f) || maxAbs > biasClip) {
            vqf->restT = 0.0f;
            vqf->restDetected = false;
        }
        memcpy(vqf->restLastGyrLp, gyrLp, sizeof(gyrLp));
        vqf->restLastSquaredDeviations[0] = squaredDeviation;
    }

    gyrNoBias[0] = gyr[0] - vqf->bias[0];
    gyrNoBias[1] = gyr[1] - vqf->bias[1];
    gyrNoBias[2] = gyr[2] - vqf->bias[2];
    gyrNorm = vqf6d_norm(gyrNoBias, 3U);
    angle = gyrNorm * vqf->coeffs.gyrTs;

    if (gyrNorm > VQF6D_EPS) {
        float step[4];
        float qNew[4];
        float c = cosf(0.5f*angle);
        float s = sinf(0.5f*angle)/gyrNorm;
        step[0]=c; step[1]=s*gyrNoBias[0]; step[2]=s*gyrNoBias[1]; step[3]=s*gyrNoBias[2];
        vqf6d_quat_multiply(vqf->gyrQuat, step, qNew);
        memcpy(vqf->gyrQuat, qNew, sizeof(qNew));
        vqf6d_normalize(vqf->gyrQuat, 4U);
    }
}

void VQF6D_UpdateAcc(VQF6D *vqf, const float acc[3])
{
    float accNormSq;
    float accEarth[3];

    if (vqf == NULL || acc == NULL) return;

    accNormSq = acc[0]*acc[0] + acc[1]*acc[1] + acc[2]*acc[2];
    if (accNormSq <= VQF6D_EPS) return;

    if (vqf->params.restBiasEstEnabled) {
        float accLp[3];
        float dev[3];
        float squaredDeviation;
        vqf6d_filter_vec(acc, 3U, vqf->params.restFilterTau, vqf->coeffs.accTs, vqf->coeffs.restAccLpB, vqf->coeffs.restAccLpA, &vqf->restAccLpState, accLp);
        dev[0]=acc[0]-accLp[0]; dev[1]=acc[1]-accLp[1]; dev[2]=acc[2]-accLp[2];
        squaredDeviation = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
        if (squaredDeviation >= vqf->params.restThAcc*vqf->params.restThAcc) {
            vqf->restT = 0.0f;
            vqf->restDetected = false;
        } else {
            vqf->restT += vqf->coeffs.accTs;
            if (vqf->restT >= vqf->params.restMinT) vqf->restDetected = true;
        }
        memcpy(vqf->restLastAccLp, accLp, sizeof(accLp));
        vqf->restLastSquaredDeviations[1] = squaredDeviation;
    }

    VQF6D_RotateVector(vqf->gyrQuat, acc, accEarth);
    vqf6d_filter_vec(accEarth, 3U, vqf->params.tauAcc, vqf->coeffs.accTs, vqf->coeffs.accLpB, vqf->coeffs.accLpA, &vqf->accLpState, vqf->lastAccLp);
    VQF6D_RotateVector(vqf->accQuat, vqf->lastAccLp, accEarth);
    vqf6d_normalize(accEarth, 3U);

    {
        float qwArg = 0.5f*(accEarth[2] + 1.0f);
        float qw;
        float accCorrQuat[4];
        float qNew[4];
        if (qwArg < 0.0f) qwArg = 0.0f;
        qw = sqrtf(qwArg);
        if (qw > 1e-6f) {
            accCorrQuat[0] = qw;
            accCorrQuat[1] = 0.5f*accEarth[1]/qw;
            accCorrQuat[2] = -0.5f*accEarth[0]/qw;
            accCorrQuat[3] = 0.0f;
        } else {
            accCorrQuat[0]=0.0f; accCorrQuat[1]=1.0f; accCorrQuat[2]=0.0f; accCorrQuat[3]=0.0f;
        }
        vqf6d_quat_multiply(accCorrQuat, vqf->accQuat, qNew);
        memcpy(vqf->accQuat, qNew, sizeof(qNew));
        vqf6d_normalize(vqf->accQuat, 4U);
        vqf->lastAccCorrAngularRate = acosf(vqf6d_clip_scalar(accEarth[2], -1.0f, 1.0f))/vqf->coeffs.accTs;
    }

    if (vqf->params.motionBiasEstEnabled || vqf->params.restBiasEstEnabled) {
        float q6D[4];
        float R[3][3];
        float Rlp[3][3];
        float biasLpRaw[2];
        float biasLp[2];

        VQF6D_GetQuat6D(vqf, q6D);
        VQF6D_QuatToRotMat(q6D, R);
        biasLpRaw[0] = R[0][0]*vqf->bias[0] + R[0][1]*vqf->bias[1] + R[0][2]*vqf->bias[2];
        biasLpRaw[1] = R[1][0]*vqf->bias[0] + R[1][1]*vqf->bias[1] + R[1][2]*vqf->bias[2];
        vqf6d_filter_vec(&R[0][0], 9U, vqf->params.tauAcc, vqf->coeffs.accTs, vqf->coeffs.accLpB, vqf->coeffs.accLpA, &vqf->motionBiasEstRLpState, &Rlp[0][0]);
        vqf6d_filter_vec(biasLpRaw, 2U, vqf->params.tauAcc, vqf->coeffs.accTs, vqf->coeffs.accLpB, vqf->coeffs.accLpA, &vqf->motionBiasEstBiasLpState, biasLp);
        vqf6d_bias_update(vqf, accEarth, Rlp, biasLp);
    }
}

void VQF6D_Update(VQF6D *vqf, const float gyr[3], const float acc[3])
{
    VQF6D_UpdateGyr(vqf, gyr);
    VQF6D_UpdateAcc(vqf, acc);
}

void VQF6D_GetQuat3D(const VQF6D *vqf, float q[4])
{
    if (vqf == NULL || q == NULL) return;
    memcpy(q, vqf->gyrQuat, 4U*sizeof(float));
}

void VQF6D_GetQuat6D(const VQF6D *vqf, float q[4])
{
    if (vqf == NULL || q == NULL) return;
    vqf6d_quat_multiply(vqf->accQuat, vqf->gyrQuat, q);
}

void VQF6D_GetBiasEstimate(const VQF6D *vqf, float bias[3], float *sigma)
{
    if (vqf == NULL) return;
    if (bias != NULL) memcpy(bias, vqf->bias, 3U*sizeof(float));
    if (sigma != NULL) {
        float row0 = fabsf(vqf->biasP[0][0]) + fabsf(vqf->biasP[0][1]) + fabsf(vqf->biasP[0][2]);
        float row1 = fabsf(vqf->biasP[1][0]) + fabsf(vqf->biasP[1][1]) + fabsf(vqf->biasP[1][2]);
        float row2 = fabsf(vqf->biasP[2][0]) + fabsf(vqf->biasP[2][1]) + fabsf(vqf->biasP[2][2]);
        float Pmax = fmaxf(row0, fmaxf(row1, row2));
        if (Pmax > vqf->coeffs.biasP0) Pmax = vqf->coeffs.biasP0;
        *sigma = sqrtf(Pmax) * VQF6D_PI / (100.0f*180.0f);
    }
}

bool VQF6D_GetRestDetected(const VQF6D *vqf)
{
    if (vqf == NULL) return false;
    return vqf->restDetected;
}
