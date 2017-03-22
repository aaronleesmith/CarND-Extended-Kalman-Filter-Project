#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // This code is basically right from the Udacity lectures.
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  double px, py, vx, vy;

  MatrixXd Hj(3,4);
  //recover state parameters
  px = x_state[0];
  py = x_state[1];
  vx = x_state[2];
  vy = x_state[3];

  double px_2 = pow(px, 2);
  double py_2 = pow(py, 2);
  double norm = sqrt(px_2 + py_2);

  //check division by zero
  if (fabs(norm) < 0.0001) {
    std::cout << "Error. Divide by zero.";
    throw 1;
  } else {

    Hj << px / norm,                                      py / norm,                                        0,          0,
          -1 * py / (px_2 + py_2),                        px / (px_2 + py_2),                               0,          0,
          py * (vx*py - vy*px) / pow(px_2 + py_2, 3/2.),  px * (vy*px - vx*py) / pow(px_2 + py_2, 3/2.),    px / norm,  py / norm;
  }


  return Hj;
}
