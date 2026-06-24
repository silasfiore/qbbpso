#include "qbbpso.h"

static inline void init_parts(struct particle *parts, int n_parts, float mag_l_prior, float std_l_prior, vec3 c) {
    /* initialize particles */
    for (int i = 0; i < n_parts; i++) {
        /* random initial attitude */
        randquat(parts[i].q0);
        /* random angular momentum direction */
        randaxis(parts[i].l);
        /* magnitude of angular momentum from prior distribution */
        parts[i].l[3] = mag_l_prior + randnrm() * std_l_prior;
        /* inertia constants are known */
        memcpy(parts[i].c, c, sizeof(vec3));
        /* determine ordering of principal moments of inertia */
        parts[i].o = get_order(c);
    }
    return;
}

static inline void update_parts(struct particle *parts, struct particle *pbest, int n_parts, int imin, float std_l_prior) {
    vec3 delta, nmid, lnext, ldiff;
    quat qmid, qdelta;
    float rbar, std_l, kappa, dsqr, d;
    /* update particles */
    for (int i = 0; i < n_parts; i++) {
        /* sample quaternion */

        /* projection of individual best onto global best */
        rbar = quat_dot(pbest[i].q0, pbest[imin].q0);
        if (rbar < 0.0f) {
            quat_negate_to(pbest[i].q0, pbest[i].q0);
            rbar = -rbar;
        }
        /* midpoint between individual best and global best */
        quat_midpoint(pbest[i].q0, pbest[imin].q0, qmid);

        /* for rbar = 1 kappa is undefined */
        if (rbar > 0.9999f) {
            dsqr = max(2.0f * (1.0f - min(rbar, 1.0f)), 0.0f);
            d = sqrtf(dsqr);
            randaxis(delta);
            quatv(qdelta, d, delta);
            quat_mul(qdelta, qmid, parts[i].q0);
        } else {
            kappa = (rbar * (4.0f - rbar * rbar)) / (1.0f - rbar * rbar);
            samplevmf4(parts[i].q0, qmid, kappa);
        }

        /* sample angular momentum vector */

        /* compute absolute deviation of individual best and global best */
        ldiff[0] = (pbest[i].l[0] * pbest[i].l[3] - pbest[imin].l[0] * pbest[imin].l[3]);
        ldiff[1] = (pbest[i].l[1] * pbest[i].l[3] - pbest[imin].l[1] * pbest[imin].l[3]);
        ldiff[2] = (pbest[i].l[2] * pbest[i].l[3] - pbest[imin].l[2] * pbest[imin].l[3]);
        std_l = vec3_norm(ldiff);

        if (std_l < std_l_prior) {
            /* sample from normal distribution */
            lnext[0] = (pbest[i].l[0] * pbest[i].l[3] + pbest[imin].l[0] * pbest[imin].l[3]) * 0.5f + randnrm() * std_l;
            lnext[1] = (pbest[i].l[1] * pbest[i].l[3] + pbest[imin].l[1] * pbest[imin].l[3]) * 0.5f + randnrm() * std_l;
            lnext[2] = (pbest[i].l[2] * pbest[i].l[3] + pbest[imin].l[2] * pbest[imin].l[3]) * 0.5f + randnrm() * std_l;

            parts[i].l[3] = vec3_norm(lnext);
            vec3_normalize_to(lnext, parts[i].l);
        } else {
            /* sample direction from vmf and magnitude from prior distribution */

            /* midpoint between individual best and global best */
            vec3_add(pbest[i].l, pbest[imin].l, nmid);
            vec3_normalize(nmid);

            /* projection of individual best onto global best */
            rbar = vec3_dot(pbest[i].l, pbest[imin].l);

            /* kappa = 1.0 handled by samplevmf3*/
            kappa = (rbar * (3.0f - rbar * rbar)) / (1.0f - rbar * rbar);
            samplevmf3(parts[i].l, nmid, kappa);
            parts[i].l[3] = (pbest[i].l[3] + pbest[imin].l[3]) * 0.5f + randnrm() * std_l_prior;
        }
    }
}

static inline float optimfunc(struct measurement *meas, int n_meas, struct particle *p, double dt) {
    int num_prop = 0; /* numerical propagation flag */
    vec3 m;           /* spin axis in body frame */
    vec3 batt;        /* body vector b rotated by current attitude qatt */
    vec3 nlb0;        /* intial direction of the angular momentum vector in the body frame */
    vec3 wb0;         /* initial angular velocity in the body frame */
    vec3 battxr;
    float phi, psi, phidot, psidot, cos_alpha, sin_alpha, residual;
    quat q0conj, qatt, qphi, qpsi, qspin;

    quat_conjugate(p->q0, q0conj);
    quat_rotatev(q0conj, p->l, nlb0);

    switch (p->o) {
        case SPHERE:  // I1 == I2 == I3
            wb0[0] = nlb0[0] * p->l[3];
            wb0[1] = nlb0[1] * p->l[3];
            wb0[2] = nlb0[2] * p->l[3];
            psidot = 0.0f;
            m[0] = 0.0f;
            m[1] = 0.0f;
            m[2] = 1.0f;
            break;
        case SYMZMAJ:  // I1 == I2 < I3
        case SYMZMIN:  // I1 == I2 > I3
            wb0[0] = nlb0[0] * p->l[3];
            wb0[1] = nlb0[1] * p->l[3];
            wb0[2] = nlb0[2] * p->l[3] / (1.0f - p->c[0]);
            psidot = p->c[0] * wb0[2];
            m[0] = 0.0f;
            m[1] = 0.0f;
            m[2] = 1.0f;
            break;

        case SYMXMAJ:  // I1 > I2 == I3
        case SYMXMIN:  // I1 < I2 == I3
            wb0[0] = nlb0[0] * p->l[3] / (1.0f - p->c[1]);
            wb0[1] = nlb0[1] * p->l[3];
            wb0[2] = nlb0[2] * p->l[3];
            psidot = p->c[1] * wb0[0];
            m[0] = 1.0f;
            m[1] = 0.0f;
            m[2] = 0.0f;
            break;

        case SYMYMAJ:  // I1 == I3 < I2
        case SYMYMIN:  // I1 == I3 > I2
            wb0[0] = nlb0[0] * p->l[3];
            wb0[1] = nlb0[1] * p->l[3] / (1.0f - p->c[2]);
            wb0[2] = nlb0[2] * p->l[3];
            psidot = p->c[2] * wb0[1];
            m[0] = 0.0f;
            m[1] = 1.0f;
            m[2] = 0.0f;
            break;

        case IXIZIY:  // I1 < I3 < I2
        case IYIZIX:  // i2 < I3 < I1
            wb0[0] = nlb0[0] * p->l[3] * (p->c[0] * p->c[1] + 1.0f) / (1.0f - p->c[1]);
            wb0[1] = nlb0[1] * p->l[3] * (p->c[0] * p->c[1] + 1.0f) / (1.0f + p->c[0]);
            wb0[2] = nlb0[2] * p->l[3];
            num_prop = 1;
            break;
        case IZIXIY:  // I3 < I1 < I2
        case IYIXIZ:  // I2 < I1 < I3
            wb0[0] = nlb0[0] * p->l[3];
            wb0[1] = nlb0[1] * p->l[3] * (p->c[1] * p->c[2] + 1.0f) / (1.0f - p->c[2]);
            wb0[2] = nlb0[2] * p->l[3] * (p->c[1] * p->c[2] + 1.0f) / (1.0f + p->c[1]);
            num_prop = 1;
            break;
        case IXIYIZ:  // I1 < I2 < I3
        case IZIYIX:  // I3 < I2 < I1
            wb0[0] = nlb0[0] * p->l[3] * (p->c[2] * p->c[0] + 1.0f) / (1.0f + p->c[2]);
            wb0[1] = nlb0[1] * p->l[3];
            wb0[2] = nlb0[2] * p->l[3] * (p->c[2] * p->c[0] + 1.0f) / (1.0f - p->c[0]);
            num_prop = 1;
            break;
        default:
            break;
    }

    float mres = 0;

    if (num_prop) {
        /* use double precision for integration */
        double qattdp[4];
        double wbdp[3];
        double cdp[3];

        qattdp[0] = (double)p->q0[0];
        qattdp[1] = (double)p->q0[1];
        qattdp[2] = (double)p->q0[2];
        qattdp[3] = (double)p->q0[3];

        wbdp[0] = (double)wb0[0];
        wbdp[1] = (double)wb0[1];
        wbdp[2] = (double)wb0[2];

        cdp[0] = (double)p->c[0];
        cdp[1] = (double)p->c[1];
        cdp[2] = (double)p->c[2];

        double trel = 0.0;

        double trelnext, dtsnap;

        int j = 0;

        while (j < n_meas) {
            trelnext = (double)meas[j].t;

            /* propagate to next time stamp */
            while (trel < trelnext) {
                dtsnap = ((trel + dt) > trelnext) ? (trelnext - trel) : dt;
                propagate_cg5(qattdp, wbdp, cdp, dtsnap);
                trel += dtsnap;
            }

            /* cast quaternion to float */
            qatt[0] = (float)qattdp[0];
            qatt[1] = (float)qattdp[1];
            qatt[2] = (float)qattdp[2];
            qatt[3] = (float)qattdp[3];

            /* compute residual */
            quat_rotatev(qatt, meas[j].b, batt);
            cos_alpha = sat(vec3_dot(batt, meas[j].r), -1.0f, 1.0f);
            vec3_cross(batt, meas[j].r, battxr);
            sin_alpha = vec3_norm(battxr);
            residual = atan2f(sin_alpha, cos_alpha);
            mres += residual;

            j++;
        }

    } else {
        /* use analytic solution for symmetric inertia tensor */
        phidot = p->l[3];
        for (int j = 0; j < n_meas; j++) {
            phi = phidot * meas[j].t;
            quatv(qphi, phi, p->l);
            psi = psidot * meas[j].t;
            quatv(qpsi, psi, m);
            quat_mul(p->q0, qpsi, qspin);
            quat_mul(qphi, qspin, qatt);

            quat_rotatev(qatt, meas[j].b, batt);
            cos_alpha = sat(vec3_dot(batt, meas[j].r), -1.0f, 1.0f);
            vec3_cross(batt, meas[j].r, battxr);
            sin_alpha = vec3_norm(battxr);
            residual = atan2f(sin_alpha, cos_alpha);
            mres += residual;
        }
    }

    mres /= (float)n_meas;

    return mres;
}

void bbpso(struct measurement *meas, int n_meas, struct particle *parts, int n_parts, int n_iter, float mag_l_prior, float std_l_prior, vec3 c, double dt_numint) {
    int i, k;
    int imin;
    float rad2deg = 180 / 3.14159265359f;

    struct particle *pbest = NULL;
    pbest = (struct particle *)malloc(n_parts * sizeof(struct particle));

    if (!pbest) {
        fprintf(stderr, "could not allocate %d particles\n", n_parts);
        return;
    }

    /* initialization */
    init_parts(parts, n_parts, mag_l_prior, std_l_prior, c);

    /* initialize particle values */
    for (i = 0; i < n_parts; i++) {
        parts[i].v = optimfunc(meas, n_meas, &parts[i], dt_numint);
    }
    /* initialize memory of best position */
    for (i = 0; i < n_parts; i++) {
        memcpy(&pbest[i], &parts[i], sizeof(struct particle));
    }
    /* initialize global best particle */
    imin = 0;
    for (i = 1; i < n_parts; i++) {
        if (pbest[i].v < pbest[imin].v) {
            imin = i;
        }
    }

    /* optimization */
    for (k = 0; k < n_iter; k++) {
        /* update particles */
        update_parts(parts, pbest, n_parts, imin, std_l_prior);

        /* compute particle values */
        for (i = 0; i < n_parts; i++) {
            parts[i].v = optimfunc(meas, n_meas, &parts[i], dt_numint);
        }

        /* update memory of best position */
        for (i = 0; i < n_parts; i++) {
            if (parts[i].v < pbest[i].v)
                memcpy(&pbest[i], &parts[i], sizeof(struct particle));
        }

        float v_avg = 0.0;

        /* update global best particle */
        for (i = 0; i < n_parts; i++) {
            if (pbest[i].v < pbest[imin].v) {
                imin = i;
            }
            v_avg += pbest[i].v;
        }
        v_avg /= (float)n_parts;

        printf("iteration %d/%d global best mean residual (degrees): %f\n", k + 1, n_iter, pbest[imin].v * rad2deg);
    }

    /* return best position of each particle */
    for (i = 0; i < n_parts; i++) {
        memcpy(&parts[i], &pbest[i], sizeof(struct particle));
    }

    free(pbest);
    return;
}

static int parsemeas(const char *line, struct measurement *meas) {
    int vals = sscanf(line, "%f %f %f %f %f %f %f\n",
                      &meas->t,
                      &meas->r[0],
                      &meas->r[1],
                      &meas->r[2],
                      &meas->b[0],
                      &meas->b[1],
                      &meas->b[2]);

    return (vals == 7);
}

#define LINE_BUF_SIZE 1024

static size_t procfilelbl(const char *filename, struct measurement **meas, size_t *n_meas_max) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "could not open file %s\n", filename);
        return -1;
    }
    char line[LINE_BUF_SIZE];
    size_t n_meas = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (n_meas >= *n_meas_max) {
            size_t new_n_meas_max = 2 * (*n_meas_max);
            struct measurement *tmp = NULL;
            tmp = (struct measurement *)realloc(*meas, new_n_meas_max * sizeof(struct measurement));
            if (tmp == NULL) {
                fprintf(stderr, "memory allocation failed\n");
                free(*meas);
                *meas = NULL;
                fclose(fp);
                return -1;
            }

            *meas = tmp;
            *n_meas_max = new_n_meas_max;
        }

        if (parsemeas(line, &(*meas)[n_meas])) {
            n_meas++;
        } else {
            break;
        }
    }
    fclose(fp);
    return n_meas;
}

int main(int argc, char *argv[]) {
    int n_parts = 100000;
    int n_iter = 100;
    double dt_numint = 0.01;
    float rad2deg = 180 / 3.14159265359f;

    if (argc < 7) {
        fprintf(stdout, "usage: %s <inputfile> <magnitude_l> <sigma_l> <c1> <c2> <c3> [-t <time step for numerical integration> -n <number of particles> -i <number of iterations> -o <outputfile>]\n", argv[0]);
        return 0;
    }

    const char *input_filename = argv[1];

    float mag_l_prior = atof(argv[2]);
    float std_l_prior = atof(argv[3]);
    vec3 c;
    c[0] = atof(argv[4]);
    c[1] = atof(argv[5]);
    c[2] = atof(argv[6]);

    char *output_filename = "particles.txt";

    for (int i = 7; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            n_parts = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            n_iter = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_filename = argv[++i];
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            dt_numint = atof(argv[++i]);
        }
    }

    size_t n_meas_max = 128;
    struct measurement *meas = NULL;
    meas = (struct measurement *)malloc(n_meas_max * sizeof(struct measurement));
    if (meas == NULL) {
        fprintf(stderr, "memory allocation failed\n");
        return 1;
    }

    struct particle *parts = NULL;
    parts = (struct particle *)malloc(n_parts * sizeof(struct particle));
    if (parts == NULL) {
        fprintf(stderr, "memory allocation failed\n");
        free(meas);
        return 1;
    }

    int n_meas = procfilelbl(input_filename, &meas, &n_meas_max);
    if (n_meas < 0) {
        free(meas);
        free(parts);
        return 1;
    }

    bbpso(meas, n_meas, parts, n_parts, n_iter, mag_l_prior, std_l_prior, c, dt_numint);

    FILE *fp = fopen(output_filename, "w");
    if (!fp) {
        fprintf(stderr, "could not open file %s\n", output_filename);
        free(meas);
        free(parts);
        return 1;
    }

    for (int i = 0; i < n_parts; i++) {
        fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", parts[i].v * rad2deg, parts[i].q0[3], parts[i].q0[0], parts[i].q0[1], parts[i].q0[2], parts[i].l[0], parts[i].l[1], parts[i].l[2], parts[i].l[3]);
    }

    fclose(fp);

    if (meas) free(meas);
    if (parts) free(parts);

    return 0;
}