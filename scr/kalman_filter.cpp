#include "kalman_filter.h"
#include <math.h> /* atan2 */
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  //std::cout << "Enter KalmanFilter::Predict()" << std::endl;
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

  //std::cout << "prediction succesfull" << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  //std::cout << "Enter KalmanFilter::Update" << std::endl;

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  //std::cout << "update (laser) succesfull" << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  /*VectorXd z_cartesian(4);
  double roh = z[0];
  double phi = z[1];
  double roh_dot = z[2];
  // Calculate point
  double x = roh * cos(phi);
  double y = roh * sin(phi);
  double vx = roh * cos(phi);
  double vy = roh * sin(phi);
  // Updating
  z_cartesian << x, y, vx, vy;
  ekf_.Hj_ = tools.CalculateJacobian(&z_cartesian);
*/
  //std::cout << "Enter KalmanFilter::UpdateEKF" << std::endl;

  double px = x_[0];
  double py = x_[1];
  double px_dot = x_[2];
  double py_dot = x_[3];
  
  // Calculate z pred
  double roh = sqrt((px*px) + (py*py));
  double phi = atan2(py, px);
  double roh_dot = ((px * px_dot) + (py * py_dot)) / roh;

  VectorXd hx(3);
  hx << roh,
        phi,
        roh_dot;

 // Updating
  VectorXd y = z - hx;
 // namalization of the y, hence to the tranformation 
  
  if (y(1) < -M_PI){
    y(1) += 2 *M_PI;
    }
  else if (y(1) > M_PI){
    y(1) -= 2*M_PI;
    }
 
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  //std::cout<<"update (radar) succesfull"<<std::endl;
}
