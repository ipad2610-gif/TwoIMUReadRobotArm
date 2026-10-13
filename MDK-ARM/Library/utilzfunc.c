#include "utilzfunc.h"
#include "math.h"

float Norm3(const float *v)
{
	return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

float Dot3(const float *a, const float *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

uint8_t Normalize3(float *v)
{
	float norm = Norm3(v);

	if (norm < 1e-3f) return 0U;

	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;

	return 1U;
}

void Mat3Vec(const float A[3][3], const float *x, float *y)
{
	// y = A*x
	y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
	y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
	y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}

void Mat3TVec(const float A[3][3], const float *x, float *y)
{
  // y = A^T*x
	y[0] = A[0][0]*x[0] + A[1][0]*x[1] + A[2][0]*x[2];
	y[1] = A[0][1]*x[0] + A[1][1]*x[1] + A[2][1]*x[2];
	y[2] = A[0][2]*x[0] + A[1][2]*x[1] + A[2][2]*x[2];
}
float WrapPi(float angle)
{
	while (angle > M_PI) angle -= 2.0f*M_PI;
	while (angle < -M_PI) angle += 2.0f*M_PI;

	return angle;
}
void Mat3MulTransB(const float A[3][3], const float B[3][3], float C[3][3])
{
    for (uint8_t i = 0U; i < 3U; i++)
    {
        for (uint8_t j = 0U; j < 3U; j++)
        {
            C[i][j] = A[i][0]*B[j][0] + A[i][1]*B[j][1] + A[i][2]*B[j][2];
        }
    }
}
void Cross3(const float *a, const float *b, float *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

float AngleAroundAxis(const float *ref,
                             const float *cur,
                             const float *axis)
{
    float p0[3], p[3], cr[3];

    float d0 = Dot3(ref, axis);
    float d  = Dot3(cur, axis);

    for (uint8_t i = 0; i < 3; i++)
    {
        p0[i] = ref[i] - d0*axis[i];
        p[i]  = cur[i] - d*axis[i];
    }

    if (!Normalize3(p0) || !Normalize3(p))
        return 0.0f;

    Cross3(p0, p, cr);

    return atan2f(Dot3(axis, cr), Dot3(p0, p));
}
void RotateAxis(const float *axis, float angle, const float *v,float *out)
{
    float c = cosf(angle);
    float s = sinf(angle);

    float axv[3];
    Cross3(axis, v, axv);

    float d = Dot3(axis, v);
    float k = (1.0f - c)*d;

    out[0] = c*v[0] + s*axv[0] + k*axis[0];
    out[1] = c*v[1] + s*axv[1] + k*axis[1];
    out[2] = c*v[2] + s*axv[2] + k*axis[2];
}
float ClampF(float x, float xmin, float xmax)
{
    if (x < xmin) return xmin;
    if (x > xmax) return xmax;
    return x;
}