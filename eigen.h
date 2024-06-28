#ifndef COMPONENTS_GEOMETRY_EIGEN_H_
#define COMPONENTS_GEOMETRY_EIGEN_H_

#include "Matrix.h"
#include "Vector.h"

typedef struct {
    Vector eigenvalues;
    Matrix eigenvectors;
} Eigen_t;

Eigen_t eig_solve(Matrix matA);
Vector eig_eigvec_of_largest_eigval(Eigen_t eigA);
void eig_free(Eigen_t eigA);

#endif  // COMPONENTS_GEOMETRY_EIGEN_H_
