# Magnetometer calibration
 magnetometer calibration algorithm with light weight linear algebra library
![image](https://github.com/michal34512/Magnetometer-calibration/assets/136522993/7fc41cf1-9d9d-41d6-8e2c-825d0499b35d)

# How to use
```c
// Create vector based on data arrays
Vector vx = vec_from_array(x, 231);
Vector vy = vec_from_array(y, 231);
Vector vz = vec_from_array(z, 231);

// Fit best ellipsoid to the data and save offest vector, translation matrix   
Callibration_t calib = calib_calibrate_sensor(vx, vy, vz);

// Print the result
vec_print(calib.offset);
mat_print(calib.rotation);

// Compensate data based on calib
// new data = calib.rotation*(old data - calib.offset)
Vector new_data_vector = vec_from_array(new_data_array, 3);
vec_sub(new_data_vector, offset);
Vector compensated_data = mat_multiply_vec(calib.rotation, new_data_vector);
```
