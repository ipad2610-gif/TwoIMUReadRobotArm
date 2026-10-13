#ifndef VQF6D_H
#define VQF6D_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * 6D (gyro + accelerometer) subset of VQF, adapted for embedded C.
 * Quaternion convention: q = [w, x, y, z].
 * Gyroscope input: rad/s.
 * Accelerometer input: m/s^2.
 * Output q rotates vectors from sensor frame S into the VQF 6D earth frame E_i.
 *
 * Based on the MIT-licensed VQF implementation by Daniel Laidig and Thomas Seel:
 * https://github.com/dlaidig/vqf
 */

#define VQF6D_MAX_FILTER_DIM 9U

typedef struct {
    float tauAcc;
    bool motionBiasEstEnabled;
    bool restBiasEstEnabled;
    float biasSigmaInit;
    float biasForgettingTime;
    float biasClip;
    float biasSigmaMotion;
    float biasVerticalForgettingFactor;
    float biasSigmaRest;
    float restMinT;
    float restFilterTau;
    float restThGyr;
    float restThAcc;
} VQF6D_Params;

typedef struct {
    bool initialized;
    uint32_t count;
    float sum[VQF6D_MAX_FILTER_DIM];
    float z1[VQF6D_MAX_FILTER_DIM];
    float z2[VQF6D_MAX_FILTER_DIM];
} VQF6D_FilterState;

typedef struct {
    float gyrTs;
    float accTs;
    float accLpB[3];
    float accLpA[2];
    float biasP0;
    float biasV;
    float biasMotionW;
    float biasVerticalW;
    float biasRestW;
    float restGyrLpB[3];
    float restGyrLpA[2];
    float restAccLpB[3];
    float restAccLpA[2];
} VQF6D_Coeffs;

typedef struct {
    VQF6D_Params params;
    VQF6D_Coeffs coeffs;

    float gyrQuat[4];
    float accQuat[4];

    float lastAccLp[3];
    VQF6D_FilterState accLpState;
    float lastAccCorrAngularRate;

    float bias[3];
    float biasP[3][3];
    VQF6D_FilterState motionBiasEstRLpState;
    VQF6D_FilterState motionBiasEstBiasLpState;

    bool restDetected;
    float restLastSquaredDeviations[2];
    float restT;
    float restLastGyrLp[3];
    VQF6D_FilterState restGyrLpState;
    float restLastAccLp[3];
    VQF6D_FilterState restAccLpState;
} VQF6D;

void VQF6D_GetDefaultParams(VQF6D_Params *params);
bool VQF6D_Init(VQF6D *vqf, float gyrTs, float accTs, const VQF6D_Params *params);
void VQF6D_Reset(VQF6D *vqf);

void VQF6D_UpdateGyr(VQF6D *vqf, const float gyr[3]);
void VQF6D_UpdateAcc(VQF6D *vqf, const float acc[3]);
void VQF6D_Update(VQF6D *vqf, const float gyr[3], const float acc[3]);

void VQF6D_GetQuat3D(const VQF6D *vqf, float q[4]);
void VQF6D_GetQuat6D(const VQF6D *vqf, float q[4]);
void VQF6D_GetBiasEstimate(const VQF6D *vqf, float bias[3], float *sigma);
bool VQF6D_GetRestDetected(const VQF6D *vqf);

void VQF6D_QuatToRotMat(const float q[4], float R[3][3]);
void VQF6D_RotateVector(const float q[4], const float v[3], float out[3]);

#ifdef __cplusplus
}
#endif

#endif
