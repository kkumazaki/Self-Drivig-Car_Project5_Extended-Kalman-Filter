#include "tools.h"
#include <iostream>
#include <cassert>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Initialize RMSE value
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  //assert(estimations.size() != 0);
  //assert(estimations.size() == ground_truth.size()); 

  if (estimations.size() == 0){
    cout << "the estimation vector size should not be zero";
    assert(estimations.size() != 0);
  }
  if (estimations.size() != ground_truth.size()){
    cout << "the estimation vector size should equal ground truth vector size";
    assert(estimations.size() == ground_truth.size()); 
  }

  // Accumulate squared residuals
  for (int i=0; i<estimations.size(); i++){
     VectorXd tmp = estimations[i] - ground_truth[i];
     tmp = tmp.array()*tmp.array();
     rmse += tmp;
  }

  // Calculate the mean
  rmse /= estimations.size();

  // Calculate the squared root
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  // Recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Pre-compute a set of terms
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = c1*c2;

  // Check division by zero // It's done in FunctionEKF.cpp
  //assert(fabs(c1) < 0.0001);
  //if (fabs(c1) <= 0.0001){
  //  cout << "jacobian should not be divided by zero";
  //  assert(fabs(c1) > 0.0001);
  //}

    // compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

   return Hj;
}
