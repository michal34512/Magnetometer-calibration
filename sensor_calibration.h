#ifndef COMPONENTS_GEOMETRY_SENSOR_CALIBRATION_H_
#define COMPONENTS_GEOMETRY_SENSOR_CALIBRATION_H_

#include "matrix.h"
#include "vector.h"

typedef struct {
    Vector offset;
    Matrix rotation;
} Callibration_t;

Callibration_t calib_calibrate_sensor(Vector x, Vector y, Vector z);

void calib_free(Callibration_t calibA);

#endif  // COMPONENTS_GEOMETRY_SENSOR_CALIBRATION_H_
