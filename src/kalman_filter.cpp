#include "kalman_filter.h"
#include "tools.h"
#include <bits/stdc++.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::cout;

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Preparation
  MatrixXd I = MatrixXd::Identity(4, 4);

  // Measurement Update
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // New estimate
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Preparation
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c4 = atan2(py,px);
  float c5 = (px*vx + py*vy)/c2;

  //MatrixXd Hj = tools.CalculateJacobian(x_); //don't use here...
  MatrixXd I = MatrixXd::Identity(4, 4);

  // Measurement Update
  VectorXd z_pred(3); // Different from Lasar, nonlinear Vector.
  z_pred << c2, c4, c5; 
  VectorXd y = z - z_pred;
  
  if(y(1) < -M_PI){
    y(1) += 2*M_PI;
  }
  else if(y(1) > M_PI){
    y(1) -= 2*M_PI;
  }

  //MatrixXd Ht = Hj.transpose(); // Different from Lasar, use Jacobian. //don't use here...
  //MatrixXd S = Hj * P_ * Ht + R_; // Different from Lasar, use Jacobian. //don't use here...
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse(); 
  MatrixXd PHt = P_ * Ht; 
  MatrixXd K = PHt * Si;

  // New estimate
  x_ = x_ + K * y;
  //P_ = (I - K * Hj) * P_;  // Different from Lasar, use Jacobian.//don't use here...
  P_ = (I - K * H_) * P_; 
}
