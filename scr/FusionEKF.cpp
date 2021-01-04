 #include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  // previous_timestamp_ = 0;   Changed uncommended

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  //observer covariance matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //observer Jakobian covariance matrix - rada
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;
  //state transition matrix
  ekf_.F_ =  MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  std::cout << "Enter FusionEKF::ProcessMeasurement" << std::endl;
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    std::cout << "Enter FusionEKF::ProcessMeasurement Initialization" << std::endl;
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      // Extract Data
      //std::cout << "Enter FusionEKF::ProcessMeasurement Initialization Radar" << std::endl;
      float roh     = measurement_pack.raw_measurements_[0];
      float phi     = measurement_pack.raw_measurements_[1];
      float roh_dot = measurement_pack.raw_measurements_[2];
      // Calculate point
      float x = roh * cos(phi);
      float y = roh * sin(phi);
      float vx = roh_dot * cos(phi);
      float vy = roh_dot * sin(phi);
      // Updating
      ekf_.x_ <<  x, y, vx, vy;
      //std::cout << "Exit FusionEKF::ProcessMeasurement Initialization Radar" << std::endl;
        }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      //std::cout << "Enter FusionEKF::ProcessMeasurement Initialization Laser" << std::endl;
      float x = measurement_pack.raw_measurements_[0];
      float y = measurement_pack.raw_measurements_[1];
      //std::cout << "writing ekf_.x_" << std::endl;
      ekf_.x_ << x, y, 1, 1;
      //std::cout << "Exit FusionEKF::ProcessMeasurement Initialization Laser" << std::endl;
    }
    // Setting the timestep
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;

    //std::cout << "initialization succesfull" << std::endl;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * 
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

   // Updating dt
  //std::cout << "Updating dt" << std::endl;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Updating F according to the new elapsed time.
  //std::cout << "Updating F" << std::endl;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // the process noise as a gaussian distribution with mean zero and covariance Q
  float dt2 = dt * dt;
  float dt3 = dt * dt * dt;
  float dt4 = dt * dt * dt * dt;

  noise_ax = 9;
  noise_ay = 9;

  float sig2ax = noise_ax;
  float sig2ay = noise_ay;
  //std::cout << "Creating Q Matrix" << std::endl;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4 / 4 * sig2ax, 0, dt3 / 2 * sig2ax, 0,
      0, dt4 / 4 * sig2ay, 0, dt3 / 2 * sig2ay,
      dt3 / 2 * sig2ax, 0, dt2 * sig2ax, 0,
      0, dt3 / 2 * sig2ay, 0, dt2 * sig2ay;

  // Measurement noise refers to uncertainty in sensor measurements R

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
    //std::cout << "Enter Update Radar" << std::endl;
    Tools tools;
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    //std::cout << "Exit  Update Radar" << std::endl;

  } else {
    // TODO: Laser updates
    //std::cout << "Enter Update Laser" << std::endl;
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    ekf_.Update( measurement_pack.raw_measurements_);
    //std::cout << "Exit Update Radar" << std::endl;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
