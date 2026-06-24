/* This file contains modified functions originally from the cglm library.
 * Original Copyright (c) 2015-2024 Recep Aslantas
 * MIT License (see LICENSE-cglm or text below)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef QBBPSO_H
#define QBBPSO_H
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define min(a, b) ((a) <= (b) ? (a) : (b))
#define max(a, b) ((a) >= (b) ? (a) : (b))
#define sat(x, lb, ub) max(min((x), (ub)), (lb))

#define MAXITERWOODSAMPLER 10
#define MAXITERBETASAMPLER 10

typedef float vec3[3];
typedef float quat[4];
typedef vec3 mat3[3];
typedef quat mat4[4];

/* ordering labels for different inertia ratios */
typedef enum {
    WRONG,
    SPHERE,   // I1 == I2 == I3
    SYMZMAJ,  // I1 == I2 < I3
    SYMZMIN,  // I1 == I2 > I3
    SYMXMAJ,  // I1 > I2 == I3
    SYMXMIN,  // I1 < I2 == I3
    SYMYMAJ,  // I1 == I3 < I2
    SYMYMIN,  // I1 == I3 > I2
    IXIYIZ,   // I1 < I2 < I3
    IXIZIY,   // I1 < I3 < I2
    IYIXIZ,   // I2 < I1 < I3
    IYIZIX,   // i2 < I3 < I1
    IZIXIY,   // I3 < I1 < I2
    IZIYIX,   // I3 < I2 < I1
} c_order;

struct particle {
    float v;   /* particle optimfunc value */
    quat q0;   /* attitude at t = 0 */
    quat l;    /* angular momentum vector: first 3 components direction, 4th component magnitude */
    vec3 c;    /* inertia constants */
    c_order o; /* ordering of principal moments of inertia */
};

struct measurement {
    vec3 r;  /* inertial vector */
    vec3 b;  /* body vector */
    float t; /* timestamp in seconds */
};

/* vector and quaternion operations adapted from the cglm library */

static inline float vec3_dot(vec3 a, vec3 b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline float vec3_norm2(vec3 v) {
    return vec3_dot(v, v);
}

static inline float vec3_norm(vec3 v) {
    return sqrtf(vec3_norm2(v));
}

static inline void vec3_zero(vec3 v) {
    v[0] = v[1] = v[2] = 0.0f;
}

static inline void vec3_scale(vec3 v, float s, vec3 dest) {
    dest[0] = v[0] * s;
    dest[1] = v[1] * s;
    dest[2] = v[2] * s;
}

static inline void vec3_normalize_to(vec3 v, vec3 dest) {
    float norm;

    norm = vec3_norm(v);

    if (norm < FLT_EPSILON) {
        vec3_zero(dest);
        return;
    }

    vec3_scale(v, 1.0f / norm, dest);
}

static inline void vec3_normalize(vec3 v) {
    vec3_normalize_to(v, v);
}

static inline void vec3_add(vec3 a, vec3 b, vec3 dest) {
    dest[0] = a[0] + b[0];
    dest[1] = a[1] + b[1];
    dest[2] = a[2] + b[2];
}
static inline void vec3_sub(vec3 a, vec3 b, vec3 dest) {
    dest[0] = a[0] - b[0];
    dest[1] = a[1] - b[1];
    dest[2] = a[2] - b[2];
}

static inline void vec3_cross(vec3 a, vec3 b, vec3 dest) {
    /* (u2.v3 - u3.v2, u3.v1 - u1.v3, u1.v2 - u2.v1) */
    dest[0] = a[1] * b[2] - a[2] * b[1];
    dest[1] = a[2] * b[0] - a[0] * b[2];
    dest[2] = a[0] * b[1] - a[1] * b[0];
}

static inline float quat_dot(quat a, quat b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

static inline float quat_norm2(quat q) {
    return quat_dot(q, q);
}

static inline float quat_norm(quat q) {
    return sqrtf(quat_norm2(q));
}

static inline void quat_identity(quat q) {
    q[0] = q[1] = q[2] = 0.0f;
    q[3] = 1.0f;
}
static inline void quat_zero(quat q) {
    q[0] = q[1] = q[2] = q[3] = 0.0f;
}

static inline void quat_scale(quat q, float s, quat dest) {
    dest[0] = q[0] * s;
    dest[1] = q[1] * s;
    dest[2] = q[2] * s;
    dest[3] = q[3] * s;
}
static inline void quat_muladds(quat q, float s, quat dest) {
    dest[0] += q[0] * s;
    dest[1] += q[1] * s;
    dest[2] += q[2] * s;
    dest[3] += q[3] * s;
}

static inline void quat_normalize_to(quat q, quat dest) {
    float dot;

    dot = quat_norm2(q);

    if (dot <= 0.0f) {
        quat_identity(dest);
        return;
    }

    quat_scale(q, 1.0f / sqrtf(dot), dest);
}

static inline void quat_normalize(quat q) {
    quat_normalize_to(q, q);
}

static inline void quatv(quat q, float angle, vec3 axis) {
    vec3 k;
    float a, c, s;

    a = angle * 0.5f;
    c = cosf(a);
    s = sinf(a);

    vec3_normalize_to(axis, k);

    q[0] = s * k[0];
    q[1] = s * k[1];
    q[2] = s * k[2];
    q[3] = c;
}

static inline void quat_mul(quat p, quat q, quat dest) {
    /*
    + (a1 b2 + b1 a2 + c1 d2 − d1 c2)i
    + (a1 c2 − b1 d2 + c1 a2 + d1 b2)j
    + (a1 d2 + b1 c2 − c1 b2 + d1 a2)k
       a1 a2 − b1 b2 − c1 c2 − d1 d2
   */

    dest[0] = p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1];
    dest[1] = p[3] * q[1] - p[0] * q[2] + p[1] * q[3] + p[2] * q[0];
    dest[2] = p[3] * q[2] + p[0] * q[1] - p[1] * q[0] + p[2] * q[3];
    dest[3] = p[3] * q[3] - p[0] * q[0] - p[1] * q[1] - p[2] * q[2];
}

static inline void mat3_mulv(mat3 m, vec3 v, vec3 dest) {
    dest[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2];
    dest[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2];
    dest[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2];
}

static inline void mat4_mulv(mat4 m, quat v, quat dest) {
    dest[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0] * v[3];
    dest[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1] * v[3];
    dest[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2] * v[3];
    dest[3] = m[0][3] * v[0] + m[1][3] * v[1] + m[2][3] * v[2] + m[3][3] * v[3];
}

static inline void quat_imag(quat q, vec3 dest) {
    dest[0] = q[0];
    dest[1] = q[1];
    dest[2] = q[2];
}

static inline float quat_real(quat q) {
    return q[3];
}

static inline void quat_rotatev(quat q, vec3 v, vec3 dest) {
    quat p;
    vec3 u, v1, v2;
    float s;

    quat_normalize_to(q, p);
    quat_imag(p, u);
    s = quat_real(p);

    vec3_scale(u, 2.0f * vec3_dot(u, v), v1);
    vec3_scale(v, s * s - vec3_dot(u, u), v2);
    vec3_add(v1, v2, v1);

    vec3_cross(u, v, v2);
    vec3_scale(v2, 2.0f * s, v2);

    vec3_add(v1, v2, dest);
}

static inline void quat_negate_to(quat q, quat dest) {
    dest[0] = -q[0];
    dest[1] = -q[1];
    dest[2] = -q[2];
    dest[3] = -q[3];
}

static inline void quat_conjugate(quat q, quat dest) {
    quat_negate_to(q, dest);
    dest[3] = -dest[3];
}

static inline void quat_midpoint(quat a, quat b, quat c) {
    float dot;

    dot = quat_dot(a, b);

    float sign = (dot >= 0.0f) ? 1.0f : -1.0f;

    c[0] = a[0] + sign * b[0];
    c[1] = a[1] + sign * b[1];
    c[2] = a[2] + sign * b[2];
    c[3] = a[3] + sign * b[3];

    quat_normalize(c);
}

/* determine order of the principal moments of inertia based on the 
 * Euler equations inertia ratios c1,c2, and c3
 */
static inline c_order get_order(vec3 c) {
    int I2gtI3 = (c[0] > 0.0);
    int I3gtI1 = (c[1] > 0.0);
    int I1gtI2 = (c[2] > 0.0);

    int I2eqI3 = (c[0] == 0.0);
    int I3eqI1 = (c[1] == 0.0);
    int I1eqI2 = (c[2] == 0.0);

    int I2ltI3 = (c[0] < 0.0);
    int I3ltI1 = (c[1] < 0.0);
    int I1ltI2 = (c[2] < 0.0);

    if (I2eqI3 && I3eqI1 && I1eqI2)
        return SPHERE;
    else if (I1eqI2 && I2gtI3 && I3ltI1)
        return SYMZMIN;
    else if (I1eqI2 && I2ltI3 && I3gtI1)
        return SYMZMAJ;
    else if (I3eqI1 && I2ltI3 && I1gtI2)
        return SYMYMIN;
    else if (I3eqI1 && I2gtI3 && I1ltI2)
        return SYMYMAJ;
    else if (I2eqI3 && I3gtI1 && I1ltI2)
        return SYMXMIN;
    else if (I2eqI3 && I3ltI1 && I1gtI2)
        return SYMXMAJ;
    else if (I2ltI3 && I3gtI1 && I1ltI2)
        return IXIYIZ;  // I1 < I2 < I3
    else if (I2gtI3 && I3gtI1 && I1ltI2)
        return IXIZIY;  // I1 < I3 < I2
    else if (I2ltI3 && I3gtI1 && I1gtI2)
        return IYIXIZ;  // I2 < I1 < I3
    else if (I2ltI3 && I3ltI1 && I1gtI2)
        return IYIZIX;  // I2 < I3 < I1
    else if (I2gtI3 && I3ltI1 && I1ltI2)
        return IZIXIY;  // I3 < I1 < I2
    else if (I2gtI3 && I3ltI1 && I1gtI2)
        return IZIYIX;  // I3 < I2 < I1

    return WRONG;
}

/* uniform random number generator */
static inline float randf(void) {
    return (float)rand() / ((float)RAND_MAX);
}

/* normal distribution sampler */
static inline float randnrm() {
    float x, y, r2;
    do {
        x = -1.0f + 2.0f * randf();
        y = -1.0f + 2.0f * randf();
        r2 = x * x + y * y;
    } while (r2 > 1.0f || r2 == 0.0f);
    /* Box-Muller transform */
    return y * sqrtf(-2.0f * logf(r2) / r2);
}

/*
 * The following Thompson Sampling / BetaSampler implementation is adapted
 * and renamed from the MSDN Magazine article:
 * "Test Run - Thompson Sampling Using C#" by James McCaffrey (February 2018)
 * URL: https://learn.microsoft.com/en-us/archive/msdn-magazine/2018/february/test-run-thompson-sampling-using-csharp
 * Copyright (c) Microsoft Corporation. All rights
 */
static inline float samplebeta(float a, float b) {
    int iter = 0;
    float x;
    float alpha, beta, gamma;
    float u1 = 0.f, u2 = 0.f, w = 0.f, v = 0.f;
    float tmp;
    alpha = a + b;
    beta = 0.0;

    if (min(a, b) <= 1.0) {
        beta = max(1.f / a, 1.f / b);
    } else {
        beta = sqrtf((alpha - 2.f) / (2 * a * b - alpha));
    }
    gamma = a + 1.f / beta;

    while (iter++ < MAXITERBETASAMPLER) {
        u1 = randf();
        u2 = randf();
        v = beta * logf(u1 / (1.f - u1));
        w = a * expf(v);
        tmp = logf(alpha / (b + w));
        if (alpha * tmp + (gamma * v) - 1.3862944f >= logf(u1 * u1 * u2)) break;
    }

    x = w / (b + w);
    return x;
}

/* used to transform samples from the vMF distribution in 4D
 * from a distribution around the standard direction to the
 * central direction given by mu
 */
static inline void hypercrossm(vec3 v0, vec3 v1, vec3 v2, vec3 v3) {
    float xy01, xy02, xy03, xy12, xy13, xy23;

    xy01 = v0[0] * v1[1] - v0[1] * v1[0];
    xy02 = v0[0] * v1[2] - v0[2] * v1[0];
    xy03 = v0[0] * v1[3] - v0[3] * v1[0];
    xy12 = v0[1] * v1[2] - v0[2] * v1[1];
    xy13 = v0[1] * v1[3] - v0[3] * v1[1];
    xy23 = v0[2] * v1[3] - v0[3] * v1[2];

    v3[0] = -(+xy23 * v2[1] - xy13 * v2[2] + xy12 * v2[3]);
    v3[1] = -(-xy23 * v2[0] + xy03 * v2[2] - xy02 * v2[3]);
    v3[2] = -(+xy13 * v2[0] - xy03 * v2[1] + xy01 * v2[3]);
    v3[3] = -(-xy12 * v2[0] + xy02 * v2[1] - xy01 * v2[2]);

    return;
}

// uniform quaternion sample
// G. Marsaglia, "Choosing a Point from the Surface of a Sphere,"
// The Annals of Mathematical Statistics, Vol. 43, No. 2, 1972, pp. 645-646.
static inline void randquat(quat q) {
    float x, y, z, u, v, w, s;
    do {
        x = -1.f + randf() * 2.f;
        y = -1.f + randf() * 2.f;
        z = x * x + y * y;
    } while (z > 1);
    do {
        u = -1.f + randf() * 2.f;
        v = -1.f + randf() * 2.f;
        w = u * u + v * v;
    } while (w > 1 || w == 0.0f);
    s = sqrtf((1 - z) / w);

    q[0] = x;
    q[1] = y;
    q[2] = s * u;
    q[3] = s * v;

    return;
}

// uniform unit vector sample
// G. Marsaglia, "Choosing a Point from the Surface of a Sphere,"
// The Annals of Mathematical Statistics, Vol. 43, No. 2, 1972, pp. 645-646.
static inline void randaxis(vec3 n) {
    float x, y, z, s;
    do {
        x = -1.f + randf() * 2.f;
        y = -1.f + randf() * 2.f;
        z = x * x + y * y;
    } while (z > 1);
    s = sqrtf(1 - z);
    n[0] = 2 * x * s;
    n[1] = 2 * y * s;
    n[2] = 1.0 - 2 * z;
    return;
}

/* used to transform samples from the vMF distribution in 3D
 * from a distribution around the standard direction to the
 * central direction given by mu, orthonormal basis with the 
 * direction mu as z-axis
 */
static inline void onb3(mat3 v, vec3 u2) {
    vec3_normalize_to(u2, v[2]);

    float sign = copysignf(1.0f, v[2][2]);
    const float a = -1.0f / (sign + v[2][2]);
    const float b = v[2][0] * v[2][1] * a;
    v[0][0] = 1.0f + sign * v[2][0] * v[2][0] * a;
    v[0][1] = sign * b;
    v[0][2] = -sign * v[2][0];
    v[1][0] = b;
    v[1][1] = sign + v[2][1] * v[2][1] * a;
    v[1][2] = -v[2][1];
    return;
}

/* vmf sample in 3D
 *Jakob, W. (2012). Numerically stable sampling of the von Mises-Fisher distribution on Sˆ2 (and other tricks). 
 *Interactive Geometry Lab, ETH Zürich, Tech. Rep, 6.
 */
static inline void samplevmf3(vec3 n, vec3 mu, float kappa) {
    float u1, u2, phi, r, cos_theta, sin_theta;
    float thr = FLT_EPSILON / 4.f;
    vec3 tmp;
    mat3 onb;

    if (kappa > 1e7) {
        memcpy(n, mu, sizeof(vec3));
    } else {
        u1 = randf();
        u2 = randf();
        phi = 2.f * 3.14159265359f * u1;
        r = kappa > thr ? log1pf(u2 * expm1f(-2.f * kappa)) / kappa : -2.f * u2;

        cos_theta = 1.f + r;
        sin_theta = sqrtf(-fmaf(r, r, 2 * r));

        tmp[0] = cosf(phi) * sin_theta;
        tmp[1] = sinf(phi) * sin_theta;
        tmp[2] = cos_theta;

        onb3(onb, mu);
        mat3_mulv(onb, tmp, n);
        vec3_normalize(n);
    }
}
/* used to transform samples from the vMF distribution in 4D
 * from a distribution around the standard direction to the
 * central direction given by mu, orthonormal basis with the 
 * direction mu as last axis
 */
static inline void onb4(mat4 v, quat u3) {
    int i, maxidx;
    float absvalue, maxabsvalue;

    quat_normalize_to(u3, v[3]);
    maxidx = 0;
    maxabsvalue = fabsf(v[3][0]);
    for (i = 1; i < 4; i++) {
        absvalue = fabsf(v[3][i]);
        if (absvalue > maxabsvalue) {
            maxidx = i;
            maxabsvalue = absvalue;
        }
    }
    quat_zero(v[2]);
    if (maxidx < 2) {
        v[0][0] = -v[3][1];
        v[0][1] = +v[3][0];
        v[0][2] = 0.0f;
        v[0][3] = 0.0f;
        quat_normalize(v[0]);
        if (v[3][2] < v[3][3]) {
            v[2][2] = 1.0f;
            quat_muladds(v[3], -v[3][2], v[2]);  // v[2] = ez -  <ez,v[0]> * v[0]
        } else {
            v[2][3] = 1.0f;
            quat_muladds(v[3], -v[3][3], v[2]);  // v[2] = ew -  <ew,v[0]> * v[0]
        }
    } else {
        v[0][0] = 0.0f;
        v[0][1] = 0.0f;
        v[0][2] = +v[3][3];
        v[0][3] = -v[3][2];
        quat_normalize(v[0]);
        if (v[3][0] < v[3][1]) {
            v[2][0] = 1.0f;
            quat_muladds(v[3], -v[3][0], v[2]);  // v[2] = ex -  <ex,v[0]> * v[0]
        } else {
            v[2][1] = 1.0f;
            quat_muladds(v[3], -v[3][1], v[2]);  // v[2] = ey -  <ey,v[0]> * v[0]
        }
    }
    quat_normalize(v[2]);
    hypercrossm(v[2], v[3], v[0], v[1]);
    quat_normalize(v[1]);
    return;
}

/* vmf sample in 4D
 * Wood, A. T. (1994). Simulation of the von Mises Fisher distribution.
 * Communications in statistics-simulation and computation, 23(1), 157-164.
 */
static inline void samplevmf4(quat q, quat mu, float lambda) {
    float b, x0, c, z, u, w, nrm;
    quat tmp;
    mat4 onb;

    float thr = 1e-4;

    if (lambda < thr) {
        randquat(q);
    } else {
        // for small lambda b -> 1
        b = (-2 * lambda + sqrtf(4 * lambda * lambda + 9.f)) / 3.f;
        // for beta -> 1 x0 -> 0
        x0 = (1 - b) / (1 + b);
        // for x0 -> 0 x0^2 -> 0
        c = lambda * x0 + 3.f * log1pf(-x0 * x0);

        int iter = 0;
        while (iter++ < MAXITERWOODSAMPLER) {
            z = samplebeta(1.5f, 1.5f);
            u = (float)rand() / ((float)RAND_MAX);
            w = (1.f - (1.f + b) * z) / (1.f - (1.f - b) * z);

            if (lambda * w + 3.f * log1p(-x0 * w) - c >= log(u)) break;
        }

        tmp[0] = randnrm();
        tmp[1] = randnrm();
        tmp[2] = randnrm();
        tmp[3] = 0.f;
        nrm = sqrtf((1.f - w * w) / (tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]));
        tmp[0] = tmp[0] * nrm;
        tmp[1] = tmp[1] * nrm;
        tmp[2] = tmp[2] * nrm;
        tmp[3] = w;
        onb4(onb, mu);
        mat4_mulv(onb, tmp, q);
        quat_normalize(q);
    }

    return;
}

/* implementation of the quaternion integration using the Crouch-Grossman method
 * Andrle, M. S., & Crassidis, J. L. (2013). Geometric integration of quaternions. 
 * Journal of Guidance, Control, and Dynamics, 36(6), 1762-1767.*/
static inline void propagate_cg5(double qatt[4], double wb[3], double c[3], double dt) {
    static const double a[] = {
        0.8177227988124852,   //a21
        0.3199876375476427,   //a31
        0.0659864263556022,   //a32
        0.9214417194464946,   //a41
        0.4997857776773573,   //a42
        -1.0969984448371582,  //a43
        0.3552358559023322,   //a51
        0.2390958372307326,   //a52
        1.3918565724203246,   //a53
        -1.1092979392113565,  //a54
    };
    static const double b[] = {
        0.1370831520630755,
        -0.0183698531564020,
        0.7397813985370780,
        -0.1907142565505889,
        0.3322195591068374,
    };

    double qstep[4] = {0.0, 0.0, 0.0, 1.0};
    double qi[4], qtmp[4];
    double k[15] = {0};
    double dw[3] = {0};
    double wi[3], dwi[3];
    double wi_nrm, scale, ha_sq, cos_ha, sin_ha;

    int ij = 0;

    for (int i = 0; i < 5; i++) {
        dwi[0] = 0.0;
        dwi[1] = 0.0;
        dwi[2] = 0.0;

        for (int j = 0; j < i; j++, ij++) {
            dwi[0] += a[ij] * k[j * 3 + 0];
            dwi[1] += a[ij] * k[j * 3 + 1];
            dwi[2] += a[ij] * k[j * 3 + 2];
        }

        wi[0] = wb[0] + dwi[0] * dt;
        wi[1] = wb[1] + dwi[1] * dt;
        wi[2] = wb[2] + dwi[2] * dt;

        k[i * 3 + 0] = c[0] * wi[1] * wi[2];
        k[i * 3 + 1] = c[1] * wi[2] * wi[0];
        k[i * 3 + 2] = c[2] * wi[0] * wi[1];

        dwi[0] += k[i * 3 + 0] * b[i];
        dwi[1] += k[i * 3 + 1] * b[i];
        dwi[2] += k[i * 3 + 2] * b[i];

        wi_nrm = sqrt(wi[0] * wi[0] +
                      wi[1] * wi[1] +
                      wi[2] * wi[2]);

        scale = 0.5 * dt * b[i];

        ha_sq = (scale * wi_nrm) * (scale * wi_nrm);

        /* 4-th order Taylor expansion */
        sin_ha = scale * (1.0 - ha_sq / 6.0 + (ha_sq * ha_sq) / 120);
        cos_ha = (1.0 - ha_sq / 2.0 + (ha_sq * ha_sq) / 24);

        qi[0] = wi[0] * sin_ha;
        qi[1] = wi[1] * sin_ha;
        qi[2] = wi[2] * sin_ha;
        qi[3] = cos_ha;

        qtmp[0] = qstep[0];
        qtmp[1] = qstep[1];
        qtmp[2] = qstep[2];
        qtmp[3] = qstep[3];

        qstep[0] = qtmp[3] * qi[0] + qtmp[0] * qi[3] + qtmp[1] * qi[2] - qtmp[2] * qi[1];
        qstep[1] = qtmp[3] * qi[1] - qtmp[0] * qi[2] + qtmp[1] * qi[3] + qtmp[2] * qi[0];
        qstep[2] = qtmp[3] * qi[2] + qtmp[0] * qi[1] - qtmp[1] * qi[0] + qtmp[2] * qi[3];
        qstep[3] = qtmp[3] * qi[3] - qtmp[0] * qi[0] - qtmp[1] * qi[1] - qtmp[2] * qi[2];
    }

    qtmp[0] = qatt[0];
    qtmp[1] = qatt[1];
    qtmp[2] = qatt[2];
    qtmp[3] = qatt[3];

    qatt[0] = qtmp[3] * qstep[0] + qtmp[0] * qstep[3] + qtmp[1] * qstep[2] - qtmp[2] * qstep[1];
    qatt[1] = qtmp[3] * qstep[1] - qtmp[0] * qstep[2] + qtmp[1] * qstep[3] + qtmp[2] * qstep[0];
    qatt[2] = qtmp[3] * qstep[2] + qtmp[0] * qstep[1] - qtmp[1] * qstep[0] + qtmp[2] * qstep[3];
    qatt[3] = qtmp[3] * qstep[3] - qtmp[0] * qstep[0] - qtmp[1] * qstep[1] - qtmp[2] * qstep[2];

    wb[0] += dt * dw[0];
    wb[1] += dt * dw[1];
    wb[2] += dt * dw[2];
}

#endif