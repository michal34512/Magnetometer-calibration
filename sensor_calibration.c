#include "sensor_calibration.h"
#include "math.h"
#include "ellipsoid_fit.h"
#include "eigen.h"

Matrix calib_ellipsoid_matrix(Vector coefA) {
    Matrix res = mat_new(3, 3);
    MAT_ELEM(res, 0, 0) = VEC_ELEM(coefA, 0);
    MAT_ELEM(res, 0, 1) = VEC_ELEM(coefA, 5);
    MAT_ELEM(res, 0, 2) = VEC_ELEM(coefA, 4);

    MAT_ELEM(res, 1, 0) = VEC_ELEM(coefA, 5);
    MAT_ELEM(res, 1, 1) = VEC_ELEM(coefA, 1);
    MAT_ELEM(res, 1, 2) = VEC_ELEM(coefA, 3);

    MAT_ELEM(res, 2, 0) = VEC_ELEM(coefA, 4);
    MAT_ELEM(res, 2, 1) = VEC_ELEM(coefA, 3);
    MAT_ELEM(res, 2, 2) = VEC_ELEM(coefA, 2);
    return res;
}

Matrix calib_calibrate_rotation(Matrix ellipMat, Vector prq, double d) {
    Eigen_t eig = eig_solve(ellipMat);
    Matrix rotation = mat_transpose(eig.eigenvectors);
    Matrix rotatedEllipMat = mat_multiply(rotation, ellipMat);
    Matrix diagonizedEllipMat = mat_multiply(rotatedEllipMat, eig.eigenvectors);
    mat_free(rotatedEllipMat);
    Vector rotPrq = vec_multiply_mat(prq, eig.eigenvectors);
    eig_free(eig);
    double a_ = MAT_ELEM(diagonizedEllipMat, 0, 0);
    double b_ = MAT_ELEM(diagonizedEllipMat, 1, 1);
    double c_ = MAT_ELEM(diagonizedEllipMat, 2, 2);
    mat_free(diagonizedEllipMat);
    double p_ = VEC_ELEM(rotPrq, 0);
    double q_ = VEC_ELEM(rotPrq, 1);
    double r_ = VEC_ELEM(rotPrq, 2);
    vec_free(rotPrq);
    Matrix gain = mat_new(3, 3);
    MAT_ELEM(gain , 0, 0) = 1/sqrt((p_*p_)/(a_*a_) + (q_*q_)/(a_*b_) + (r_*r_)/(a_*c_) - d/a_);
    MAT_ELEM(gain , 1, 1) = 1/sqrt((p_*p_)/(a_*b_) + (q_*q_)/(b_*b_) + (r_*r_)/(b_*c_) - d/b_);
    MAT_ELEM(gain , 2, 2) = 1/sqrt((p_*p_)/(a_*c_) + (q_*q_)/(b_*c_) + (r_*r_)/(c_*c_) - d/c_);
    Matrix res = mat_multiply(gain, rotation);
    mat_free(gain);
    mat_free(rotation);
    return res;
}

Vector calib_calibrate_offset(Matrix ellipMat, Vector prq) {
    mat_inv_3x3(ellipMat);
    Vector res = mat_multiply_vec(ellipMat, prq);
    vec_negate(res);
    return res;
}

Callibration_t calib_calibrate_sensor(Vector x, Vector y, Vector z) {
    Callibration_t res;
    res.offset = (Vector)0;
    res.rotation = (Matrix)0;
    if (x->size < 6) return res;
    Ellipsoid_t ellip = ellipsoid_fit(x, y, z);
    Matrix M = calib_ellipsoid_matrix(ellip.coefA);
    Vector prq = vec_copy_subvec(ellip.coefB, 0, 3);
    double d = ellip.coefB->data[3];
    ellipsoid_free(ellip);
    res.rotation = calib_calibrate_rotation(M, prq, d);
    res.offset = calib_calibrate_offset(M, prq);
    mat_free(M);
    vec_free(prq);
    return res;
}

void calib_free(Callibration_t calibA) {
    vec_free(calibA.offset);
    mat_free(calibA.rotation);
}
