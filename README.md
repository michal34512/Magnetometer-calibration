# Magnetometer calibration
Magnetometer calibration algorithm using ellipsoid fit and least squares method with light weight linear algebra library. Unlike similar algorithms, it <b>does not</b> require external linear-algebra library thus it's perfect for:
<ul>
  <li>ESP32</li>
  <li>Arduino</li>
  <li>Lightweight applications</li>
</ul>

<p align="center">
Visualization (graphs not included in the code)
<img src="https://github.com/michal34512/Magnetometer-calibration/assets/136522993/7fc41cf1-9d9d-41d6-8e2c-825d0499b35d" border="10"/>
</p>



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

// Compensate new data point based on calib
// new data = calib.transform*(old data - calib.offset)
Vector new_data_vector = vec_new(3);
VEC_X(new_data_vector) = 18.923634;
VEC_Y(new_data_vector) = 94.943907;
VEC_Z(new_data_vector) = 27.951207;
calib_calibrate_point(calib, new_data_vector);

// Print calibrated vector
vec_print(new_data_vector);
```
# Algorithm perfomance
Average variance of points length in relation to calibration data noise (random mag. 100 data points):
![image](https://github.com/michal34512/Magnetometer-calibration/assets/136522993/df930068-2c8a-4ca8-b1ef-d4f1a6730eee)

Average chance of algorithm failure in relation to data point count (random mag. data with 10% noise):
![image](https://github.com/michal34512/Magnetometer-calibration/assets/136522993/df975675-cb75-4c7f-963d-f025fb2569e8)

Example of data calibration with 50% noise:
![image](https://github.com/michal34512/Magnetometer-calibration/assets/136522993/9752cadd-9f7c-496e-8559-d9e351b3e218)

