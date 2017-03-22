#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"
#include <cmath>


class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
   * Handle initial value and initialization. Only called once.
   * @param measurement_pack
   */
  void Init(const MeasurementPackage &measurement_pack);

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

private:
  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  MatrixXd R_laser_;
  MatrixXd R_radar_;
  MatrixXd H_laser_;
  MatrixXd Hj_;

  double noise_ax_;
  double noise_ay_;

  double CalculateDt(double new_timestamp);
  void UpdateProcessCovarianceMatrix(double dt);
  void UpdateStateTransitionMatrix(double dt);
  VectorXd CalculateLaserMeasurement();
  VectorXd CalculateRadarMeasurement();
  double NormalizedAtan(double px, double py);
};

#endif /* FusionEKF_H_ */
