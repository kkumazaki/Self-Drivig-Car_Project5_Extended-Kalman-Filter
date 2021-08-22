#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
//1. Initialize variables and matrices.
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix - laser
  H_laser_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0;        

  //Jacobian matrix - radar (set 0 here)
  Hj_ << 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0,        
         0.0, 0.0, 0.0, 0.0;

  //state transition matrix (set delta_t = 0.5sec here)
  ekf_.F_ << 1.0, 0.0, 0.5, 0.0,
             0.0, 1.0, 0.0, 0.5, 
             0.0, 0.0, 1.0, 0.0, 
             0.0, 0.0, 0.0, 1.0;            

  //process noises (set large numbers for initialization)
  ekf_.P_ << 1000.0, 0.0, 0.0, 0.0,
             0.0, 1000.0, 0.0, 0.0, 
             0.0, 0.0, 1000.0, 0.0, 
             0.0, 0.0, 0.0, 1000.0;    

  //process covariance matrix (set 0 here)
  ekf_.Q_ << 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 
             0.0, 0.0, 0.0, 0.0, 
             0.0, 0.0, 0.0, 0.0;   

  //measurement noises
  // Any more??
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
//2. Initialize Kalman Filter position vector with 1st sensor measurements.
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Initialize the state ekf_.x_ with the first measurement.
      // Convert radar from polar to cartesian coordinates

      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rhodot = measurement_pack.raw_measurements_(2);

      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rhodot * cos(phi);
      float vy = rhodot * sin(phi);
    
      ekf_.x_ << px, py, vx, vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize the state ekf_.x_ with the first measurement.
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);

      ekf_.x_ << px, py, 0.0, 0.0;
 
    }

    // update timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

//3. Modify F and Q prior to prediction step based on the elapsed time delta_t.
  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // micro second --> second

  float noise_ax = 9.;
  float noise_ay = 9.;

  float delta_t2 = delta_t * delta_t;
  float delta_t3 = delta_t2 * delta_t;
  float delta_t4 = delta_t3 * delta_t;
  
  ekf_.F_ << 1.0, 0.0, delta_t, 0.0,
            0.0, 1.0, 0.0, delta_t, 
            0.0, 0.0, 1.0, 0.0, 
            0.0, 0.0, 0.0, 1.0;

  ekf_.Q_ << delta_t4/4.0*noise_ax, 0.0, delta_t3/2.0*noise_ax, 0.0,
        0.0, delta_t4/4.0*noise_ay, 0.0, delta_t3/2.0*noise_ay, 
        delta_t3/2.0*noise_ax, 0.0, delta_t2*noise_ax, 0.0, 
        0.0, delta_t3/2.0*noise_ay, 0.0, delta_t2*noise_ay;

  // update timestamp
  previous_timestamp_ = measurement_pack.timestamp_;

//4. Call update step for either Radar and Laser sensor measurements.
//   There are different functions for updating Radar and Laser.

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    float px = ekf_.x_[0];
    float py = ekf_.x_[1];

    float pxpy = px*px + py*py;

    // check the validity that px, py should not be zero 
    //assert(pxpy > 0.0001);
    if (pxpy <= 0.0001){
      cout << "check the validity that px, py should not be zero ";
      assert(pxpy > 0.0001);
    }

    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;  
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;  
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
