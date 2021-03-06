#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>
#include <stdexcept>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

FusionEKF::FusionEKF() {
  x_ = VectorXd(4);
  P_ = MatrixXd(4,4);
  previous_timestamp_ = 0;
}

FusionEKF::~FusionEKF() {}

const VectorXd& FusionEKF::GetStateVector() {
  return x_;
}


void FusionEKF::Initialize(const MeasurementPackage &m) {
  //initialize state vector with first measurement
  if (m.sensor_type_ == MeasurementPackage::RADAR) {
    float ro = m.raw_measurements_[0];
    float theta = m.raw_measurements_[1];
    float ro_dot = m.raw_measurements_[2];
    x_ << ro*cos(theta), ro*sin(theta), 0, 0;
  }
  else if (m.sensor_type_ == MeasurementPackage::LASER) {
    float px = m.raw_measurements_[0];
    float py = m.raw_measurements_[1];
    x_ << px, py, 0, 0;
  }
  else {
    throw std::invalid_argument("Invalid sensor type");
  }

  //initialize state covariance matrix
	P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  //record timestamp of first measurement
  previous_timestamp_ = m.timestamp_;
}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &m) {
  //initialize state with first measurement
  if (previous_timestamp_ == 0) {
    Initialize(m);
    return;
  }

  //compute the time elapsed between the current and previous measurements
 	float dt = (m.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
 	previous_timestamp_ = m.timestamp_;

  //prediction
  Predict(dt);

  //measurement update
  if (m.sensor_type_ == MeasurementPackage::RADAR) {
    float ro = m.raw_measurements_[0];
    float theta = m.raw_measurements_[1];
    float ro_dot = m.raw_measurements_[2];
    VectorXd z = VectorXd(3);
    z << ro, theta, ro_dot;
    RadarUpdate(z);
  }
  else if (m.sensor_type_ == MeasurementPackage::LASER) {
    float px = m.raw_measurements_[0];
    float py = m.raw_measurements_[1];
    VectorXd z = VectorXd(2);
    z << px, py;
    LaserUpdate(z);
  }
  else {
    throw std::invalid_argument("Invalid sensor type");
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}


void FusionEKF::Predict(const float dt) {
  //calculate state transition matrix
  MatrixXd F = MatrixXd(4,4);
  F << 1, 0, dt, 0,
      0, 1, 0, dt,
      0, 0, 1, 0,
      0, 0, 0, 1;

  //calculate process noise covariance matrix
  MatrixXd Q = MatrixXd(4,4);
  MatrixXd G = MatrixXd(4,2);
  MatrixXd Qv = MatrixXd(2,2);
  float noise_ax = 9;
  float noise_ay = 9;

  G << dt*dt/2, 0,
      0, dt*dt/2,
      dt, 0,
      0, dt;
  Qv << noise_ax, 0,
      0, noise_ay;
  Q = G * Qv * G.transpose();

  //update state vector & covariance matrix
  x_ = F * x_;
  P_ = F * P_ * F.transpose() + Q;
}


void FusionEKF::LaserUpdate(const VectorXd &z) {
  MatrixXd H(2,4);
  H << 1, 0, 0, 0,
      0, 1, 0, 0;

  MatrixXd R(2,2);
  R << 0.0225, 0,
      0, 0.0225;

  MatrixXd Ht = H.transpose();
  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}


void FusionEKF::RadarUpdate(const VectorXd &z) {
  MatrixXd H = CalculateJacobian();

  MatrixXd R(3,3);
  R << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  MatrixXd Ht = H.transpose();
  VectorXd y = z - GetPolarStateVector();
  y(1) = NormalizeAngle(y(1));
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}


VectorXd FusionEKF::GetPolarStateVector() {
  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  //calculate h(x)
  VectorXd h(3);
  float ro = sqrt(px*px + py*py);
  h << ro, atan2(py, px), (px*vx + py*vy) / ro;

  return h;
}


float FusionEKF::NormalizeAngle(float rad) {
  const float PI = 3.14159265358979f;
  while (rad < -PI) rad += 2*PI;
  while (rad > PI) rad -= 2*PI;
  return rad;
}


MatrixXd FusionEKF::CalculateJacobian() {
  MatrixXd Hj(3,4);

  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero
  if (fabs(c1) < 0.0001) {
  	cout << "CalculateJacobian () - Error - Division by Zero" << endl;
  	return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
  	  -(py/c1), (px/c1), 0, 0,
  	  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
