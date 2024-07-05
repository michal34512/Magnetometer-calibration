# Magnetometer calibration
 magnetometer calibration algorithm with light weight linear algebra library
![image](https://github.com/michal34512/Magnetometer-calibration/assets/136522993/7fc41cf1-9d9d-41d6-8e2c-825d0499b35d)

# How to use
```c
// data arrays

double x[231] = {18.923634, 20.785405, 23.265396, ...};
double y[231] = {94.943907, 103.939205, 113.472091, ...};
double z[231] = {27.951207, 32.958822, 31.772349, ...};

// Create vector based on data arrays
Vector vx = vec_from_array(x, 231);
Vector vy = vec_from_array(y, 231);
Vector vz = vec_from_array(z, 231);

// Fit best ellipsoid to the data and save offest vector, translation matrix   
Callibration_t calib = calib_calibrate_sensor(vx, vy, vz);

// Print the result
vec_print(calib.offset);
mat_print(calib.transform);

// Compensate data based on calib
// new data = calib.transform*(old data - calib.offset)
calib_calibrate_point(calib, new_data_vector);
```
# Algorithm perfomance
Average variance of points length in relation to calibration data noise:
![image](https://github.com/michal34512/Magnetometer-calibration/assets/136522993/a787f111-d35e-47f4-a695-2f152841c7c6)
