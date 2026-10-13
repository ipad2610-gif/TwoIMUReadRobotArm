#ifndef _UTILZFUNC_H
#define _UTILZFUNC_H

#include "stdint.h"
#define M_PI 3.1415926

float ClampF(float x, float xmin, float xmax);
float Norm3(const float *v);
float Dot3(const float *a, const float *b);
uint8_t Normalize3(float *v);
void Mat3Vec(const float A[3][3], const float *x, float *y);
void Mat3TVec(const float A[3][3], const float *x, float *y);
float WrapPi(float angle);
void Mat3MulTransB(const float A[3][3], const float B[3][3], float C[3][3]);
void Cross3(const float *a, const float *b, float *c);
float AngleAroundAxis(const float *ref,
                             const float *cur,
                             const float *axis);
void RotateAxis(const float *axis, float angle, const float *v,float *out);
														 
#endif /*UTILZFUNC_H*/