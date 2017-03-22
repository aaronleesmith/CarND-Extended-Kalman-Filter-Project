#include "kalman_filter.h"

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in,
                        MatrixXd &I_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  I_ = I_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

/**
 * We want this Kalman filter class to be independent of laser / radar.
 * So, the calculation of z must be done in the FusionEKF class.
 * @param z
 */
void KalmanFilter::Update(const VectorXd &measurement_error) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = (P_ * Ht) * S.inverse();
  x_ = x_ + (K * measurement_error);
  P_ = (I_ - (K * H_)) * P_;
}