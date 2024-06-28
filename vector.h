#ifndef COMPONENTS_GEOMETRY_VECTOR_H_
#define COMPONENTS_GEOMETRY_VECTOR_H_

#include "stdbool.h"

#define VEC_ELEM(VEC, ID) VEC->data[ID]
#define VEC_X(VEC) VEC->data[0]
#define VEC_Y(VEC) VEC->data[1]
#define VEC_Z(VEC) VEC->data[2]

typedef struct {
    unsigned int size;
    double *data;
} Vector_t, *Vector;

Vector vec_new(unsigned int size);
Vector vec_from_array(double *array, unsigned int size);
Vector vec_copy_subvec(Vector vecA, unsigned int felem, unsigned int elems);
void vec_replace(Vector vecA, Vector vecB);

void vec_fill(Vector vecA, double val);

double vec_dot_product(Vector vecA, Vector vecB);
void vec_multiply_scalar(Vector vecA, double val);
void vec_sub(Vector vecA, Vector vecB);
double vec_norm(Vector vecA);

float vec_angle_between(Vector vecA, Vector vecB);
float vec_angle_between_2D(Vector vecA, Vector vecB);

bool vec_normalize(Vector vecA);
void vec_negate(Vector vecA);
void vec_rotate_x(Vector vec, float rad);
void vec_rotate_y(Vector vec, float rad);
void vec_rotate_z(Vector vec, float rad);

bool vec_contains(Vector vecA, double val, double tol);
bool vec_equal(Vector vecA, Vector vecB, double tol);

void vec_print(Vector vecA);

void vec_free(Vector vecA);
void vec_from_array_free(Vector vecA);
#endif  // COMPONENTS_GEOMETRY_VECTOR_H_
