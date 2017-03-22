#include "FusionEKF.h"
#include "tools.h"
#include "Eigen"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#pragma clean diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  ekf_.x_ = VectorXd(4);
  ekf_.I_ = MatrixXd::Identity(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

  // Value initialization
  ekf_.F_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;

  ekf_.Q_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;

  R_laser_ << 0.0225, 0,
              0,      0.0225;

  R_radar_ << 0.0225  ,0,     0,
              0,      0.0225, 0,
              0,      0,      0.0225;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  noise_ax_ = 9;
  noise_ay_ = 9;

  ekf_.P_ <<  1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Init(const MeasurementPackage &measurement_pack) {
  double px, py, vx, vy;

  /**
   * Initialize x based on initial measurement.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    double range = measurement_pack.raw_measurements_[0];
    double bearing = measurement_pack.raw_measurements_[1];
    double range_rate = measurement_pack.raw_measurements_[2];

    px = range * cos(bearing);
    py = range * sin(bearing);
    vx = range_rate * cos(bearing);
    vy = range_rate * sin(bearing);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
    Initialize state.
    */
    px = measurement_pack.raw_measurements_[0];
    py = measurement_pack.raw_measurements_[1];
    vx = 0;
    vy = 0;
  }

  // Initialize the state vector.
  ekf_.x_ << px, py, vx, vy;

  previous_timestamp_ = measurement_pack.timestamp_;
  is_initialized_ = true;
  return;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    Init(measurement_pack);
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  double dt = CalculateDt(measurement_pack.timestamp_);
  previous_timestamp_ = measurement_pack.timestamp_;

  UpdateStateTransitionMatrix(dt);
  UpdateProcessCovarianceMatrix(dt);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   * Measurement error is calculated based on either radar or laser, then passed into the single update function.
   */
  VectorXd measurement_error;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    try {
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      measurement_error = measurement_pack.raw_measurements_ - CalculateRadarMeasurement();
    } catch(int e) {
      return;
    }
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    measurement_error = measurement_pack.raw_measurements_ - CalculateLaserMeasurement();
  }

  ekf_.Update(measurement_error);

  // print the output
  cout << "==== New predicted state and covariance ====" << endl;
  cout << "---- x_ ----" << endl << ekf_.x_ << endl << endl;
  cout << "-----P_ ----" << endl << ekf_.P_ << endl << endl;
  cout << "=============================================" << endl;
}

double FusionEKF::CalculateDt(double new_timestamp) {
  double dt = (new_timestamp - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  if (dt < .001) { dt = 0.00001; }
  return dt;
}

void FusionEKF::UpdateProcessCovarianceMatrix(double dt) {
  double dt_4 = pow(dt, 4);
  double dt_3 = pow(dt, 3);
  double dt_2 = pow(dt, 1);

  ekf_.Q_ << dt_4 / 4 * noise_ax_, 0,                      dt_3 / 2 *noise_ax_,       0,
             0,                    dt_4 / 4 * noise_ay_,   0,                         dt_3 / 2 *noise_ay_,
             dt_3 /2 * noise_ax_,  0,                      pow(dt, 2) * noise_ax_,    0,
             0,                    dt_3 /2*noise_ay_,      0,                         dt_2 * noise_ay_;
}

void FusionEKF::UpdateStateTransitionMatrix(double dt) {
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
}

VectorXd FusionEKF::CalculateLaserMeasurement() {
  return ekf_.H_ * ekf_.x_;
}

VectorXd FusionEKF::CalculateRadarMeasurement() {
  double px, py, vx, vy, range, bearing, range_rate;

  px = ekf_.x_[0];
  py = ekf_.x_[1];
  vx = ekf_.x_[2];
  vy = ekf_.x_[3];

  range = sqrt(pow(px, 2) + pow(py, 2));
  bearing = NormalizedAtan(px, py);
  range_rate = (px * vx + py * vy) / range;

  VectorXd z = VectorXd(3);
  z << range, bearing, range_rate;
  return z;
}

double FusionEKF::NormalizedAtan(double px, double py) {
  double result = atan(py / px);
  while (result < -M_PI || result > M_PI) {
    result -= 2 * M_PI;
  }
  return result;
}