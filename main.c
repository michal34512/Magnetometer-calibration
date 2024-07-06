#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "geometry.h"

// THIS A DEMO OF MAGNETOMETER CALIBRATION. DEMO CONSISTS OF:
// 1. Generating random mangetometer data (which you would normally read from the sensor)
// 2. Calculating calibration data (which can be later used to calibrate new data points)
// 3. Checking if calculating calibration data was successful
// 4. Calibrating data points previously used for calibration process

typedef struct {
    Vector vx;
    Vector vy;
    Vector vz;
} mag_generated_data_t;

#define MIN_GAIN 0.5f
#define MAX_GAIN 10
#define MIN_OFFSET -10
#define MAX_OFFSET 10
#define MAX_NOISE 0.1f
#define MIN_DATA_LEN 15
#define MAX_DATA_LEN 100

mag_generated_data_t generate_mag_data(int n) {
    mag_generated_data_t res;
    res.vx = vec_new(n);
    res.vy = vec_new(n);
    res.vz = vec_new(n);
    
    //Random rotation matrix
    float angleX = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;
    float angleY = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;
    float angleZ = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;
    Matrix rotationX = mat_rotation_x(angleX);
    Matrix rotationY = mat_rotation_x(angleY);
    Matrix rotationZ = mat_rotation_x(angleZ);
    Matrix rotationXY = mat_multiply(rotationX, rotationY);
    Matrix rotation = mat_multiply(rotationXY, rotationZ);
    mat_free(rotationX);
    mat_free(rotationY);
    mat_free(rotationZ);
    mat_free(rotationXY);
    
    //Random gain matrix
    Matrix gain = mat_new(3, 3);
    MAT_ELEM(gain, 0, 0) = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;
    MAT_ELEM(gain, 1, 1) = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;
    MAT_ELEM(gain, 2, 2) = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;

    //Random offset
    Vector offset = vec_new(3);
    VEC_ELEM(offset, 0) = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;
    VEC_ELEM(offset, 1) = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;
    VEC_ELEM(offset, 2) = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;

    // Generate data
    Matrix rotationGain = mat_multiply(rotation, gain);
    mat_free(rotation);
    mat_free(gain);
    for(int i = 0; i < n; i++) {
        // Data on sphere
        Vector point = vec_new(3);
        VEC_ELEM(point, 0) = ((float)rand() / (float)RAND_MAX) * 2 - 1; // -1 ... 1
        VEC_ELEM(point, 1) = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        VEC_ELEM(point, 2) = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        vec_normalize(point);
        
        // Add noise
        vec_multiply_scalar(point, (((float)rand() / (float)RAND_MAX) * 2 - 1) * MAX_NOISE + 1);

        // Data on ellipsoid
        Vector pointElip = mat_multiply_vec(rotationGain, point);
        vec_add(pointElip, offset);
        VEC_ELEM(res.vx, i) = VEC_ELEM(pointElip, 0);
        VEC_ELEM(res.vy, i) = VEC_ELEM(pointElip, 1);
        VEC_ELEM(res.vz, i) = VEC_ELEM(pointElip, 2);
        vec_free(point);
        vec_free(pointElip);
    }
    mat_free(rotationGain);
    vec_free(offset);
    return res;
}

int main() {
    srand(time(NULL));
    
    // Simulate random magnetometer data (length range: MAX_DATA_LEN - MIN_DATA_LEN)
    mag_generated_data_t data = generate_mag_data((int)(((float)rand() / (float)RAND_MAX) * (MAX_DATA_LEN - MIN_DATA_LEN)) + MIN_DATA_LEN);
    
    // Calulate calibration data (offset vector and transformation matrix)
    Callibration_t calib = calib_calibrate_sensor(data.vx, data.vy, data.vz);

    // Calibration might sometimes be unsuccessful (when there is not enough data points, or when it's impossible to fit an ellipsoid to these points)
    printf("Calculating calibration data of data points was %s\n", calib_calibration_success(calib) ? "SUCCESFUL" : "UNSUCCESFUL");

    printf("Variance of distances before calibration: %f\n", square_distance_variance(data.vx, data.vy, data.vz));

    // Calibrate generated points
    calib_calibrate_multiple_points(calib, data.vx, data.vy, data.vz);

    printf("Variance of distances after calibration: %f\n", square_distance_variance(data.vx, data.vy, data.vz));

    // Free allocated data
    calib_free(calib);
    vec_free(data.vx);
    vec_free(data.vy);
    vec_free(data.vz);
    return 0;
}
