#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"

class FusionEKF {
public:
  FusionEKF();
  virtual ~FusionEKF();
  const Eigen::VectorXd& GetStateVector();
  void ProcessMeasurement(const MeasurementPackage &m);

private:
  long long previous_timestamp_;
  Eigen::VectorXd x_;
  Eigen::MatrixXd P_;

  void Initialize(const MeasurementPackage &m);
  void Predict(const float dt);
  void LaserUpdate(const Eigen::VectorXd &z);
  void RadarUpdate(const Eigen::VectorXd &z);
  Eigen::VectorXd GetPolarStateVector();
  float NormalizeAngle(float);
  Eigen::MatrixXd CalculateJacobian();
};

#endif /* FusionEKF_H_ */
