#include "wmm_compensation.h"
#include "math.h"

static void rodrigues_formula(Vector v, Vector k, const double angle) {
    double sangle = sin(angle);
    double cangle = cos(angle);
    double vx = cangle*VEC_X(v) - sangle*(VEC_Y(k)*VEC_Z(v) - VEC_Z(k)*VEC_Y(v)) - VEC_X(k)*(cangle - 1)*(VEC_X(k)*VEC_X(v) + VEC_Y(k)*VEC_Y(v) + VEC_Z(k)*VEC_Z(v));
    double vy = sangle*(VEC_X(k)*VEC_Z(v) - VEC_Z(k)*VEC_X(v)) + cangle*VEC_Y(v) - VEC_Y(k)*(cangle - 1)*(VEC_X(k)*VEC_X(v) + VEC_Y(k)*VEC_Y(v) + VEC_Z(k)*VEC_Z(v));
    VEC_Z(v) = cangle*VEC_Z(v) - sangle*(VEC_X(k)*VEC_Y(v) - VEC_Y(k)*VEC_X(v)) - VEC_Z(k)*(cangle - 1)*(VEC_X(k)*VEC_X(v) + VEC_Y(k)*VEC_Y(v) + VEC_Z(k)*VEC_Z(v));
    VEC_X(v) = vx;
    VEC_Y(v) = vy;
}

static void wmm_compensate_dip(Vector v, Vector g, const double dip) {
    Vector perp = vec_cross_product(v, g);
    vec_normalize(perp);
    rodrigues_formula(v, perp, dip);
    vec_free(perp);
}

void wmm_compensate(Vector vecA, Vector gravity, const double dec, const double dip) {
    vec_normalize(vecA);
    vec_normalize(gravity);
    rodrigues_formula(vecA, gravity, dec);
    wmm_compensate_dip(vecA, gravity, dip);
}

void wmm_compensate_dec_only(Vector vecA, Vector gravity, const double dec) {
    vec_normalize(vecA);
    vec_normalize(gravity);
    rodrigues_formula(vecA, gravity, dec);
}

void wmm_compensate_dip_only(Vector vecA, Vector gravity, const double dip) {
    vec_normalize(vecA);
    vec_normalize(gravity);
    wmm_compensate_dip(vecA, gravity, dip);
}